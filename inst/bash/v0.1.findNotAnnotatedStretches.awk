# goal: evaluate prediction vs annotation; we want to give higher scores to predictions that get "a lot of introns n a row"
# by Hagen; I think to have seen scoring schemes like this somewhere ... not sure where exactly ... shoudl think about this for possible citations
# requires -v anno=gffFileName -v pred=gffFileName -v keyword
function max(a,b){
    if(a>b){return(a);}
    return(b);
}
function allOptionsThere(){
    if(keyword!="exon" && keyword!="intron"){
	print keyword" as a keyword is not allowed; exiting" > "/dev/stderr";
	exit;
    }
    if(!anno){
	print "-v anno=gffFileName must be provided; exiting" > "/dev/stderr";
	exit;
    }
    if(!pred){
	print "-v pred=gffFileName must be provided; exiting" > "/dev/stderr";
	exit;
    }
    if(!outfile){
	print "-v outfile=<fileNameTrunk> must be provided; exiting" > "/dev/stderr";
	exit;
    }
    if(!trIDcolumn){
	print "-v trIDcolumn=<int> must be provided; exiting" > "/dev/stderr";
	exit;
    }
    

}

BEGIN{

    print "## 0. preliminaries" > "/dev/stderr";
    # input/output variables
    FS="\t";
    OFS="\t";

    # need to have the required options
    allOptionsThere();
    
    print "## 1. annotation: assuming non-zipped, genomically ordered gff in the format Nicholas provided" > "/dev/stderr";
    print "## 1a. parsing" > "/dev/stderr";
    comm="zcat "anno;
    while(comm | getline){

	# only considering lines for elements of teh desired type (exon|intron)
	if($3!=keyword){continue;}

	# for summary stats (and sanity checks)
	annotatedChroms[$1]++;
	annotatedStrands[$7]++;
	
	transcriptID=$trIDcolumn;
	elementID=$1"_"$4"_"$5"_"$7;
	elementsAnno[elementID]=1;

	if(transcriptID in lastBorderAnno){
	    if($4 <= lastBorderAnno[transcriptID]){
		print "ERROR: input not ordered: transcriptID="transcriptID"; lastBorderAnno[transcriptID]=" lastBorderAnno[transcriptID]"; line="$0 > "/dev/stderr";
		exit;
	    }
	    lastBorderAnno[transcriptID]=$5;
	    elementNumberAnno[transcriptID]++;
	    transcript2ElementsAnno[transcriptID"\t"elementNumberAnno[transcriptID]]=elementID;	   
	}
	else{
	    lastBorderAnno[transcriptID]=$5;
	    elementNumberAnno[transcriptID]=1;
	    transcript2ElementsAnno[transcriptID"\t"elementNumberAnno[transcriptID]]=elementID;	   	    
	}	 	           
    }
    print "## 1b. getting all the stretches" > "/dev/stderr";
    # considering all transcripts
    for(tr in lastBorderAnno){

	# considering all stretch lengths
	for(stretchLen=1;stretchLen<=elementNumberAnno[tr];stretchLen++){
	    
	    # considering all start-indexes for stretches of this lengths	   
	    for(i=1;i<=elementNumberAnno[tr]-stretchLen+1;i++){
		stretch="";
		for(j=i;j<=i+stretchLen-1;j++){
		    stretch=stretch";%;"transcript2ElementsAnno[tr"\t"j];
		}
		if(stretch in annotatedStretch){
		    if(annotatedStretch[stretch]!~tr){
			annotatedStretch[stretch]=annotatedStretch[stretch]";%;"tr;
		    }
		}
		else{
		    annotatedStretch[stretch]=tr;
		}
	    }
	}
    }

    if(check){
	  print "## 1b.checks printing all annotated stretches" > "/dev/stderr";
	  file=check".annotated.stretches";
	  for(k in score){
	      print k,score[k] > file;
	  }    
    }

#parses the tmp/rnaseq.introns.gff file
#counts the number of chromosomes and strands
#transcript ID is in column 9
#elementID is the chr_start_stop_strand ID
#the elementsPred for all elementIDs is 1 by default
#lastBorderPred stores, for each transcript ID the final stop position of the transcriptID
#if this is the first time the transcriptID has been looked at...
#the elementNumberPred for the transcript ID is 1
#the transcript2ElementsPred gains an entry for the "transcriptID elementnumberPred[transcriptID](1 in this case) equal to the elementID (chr_start_stop_strand)
# Now, if the transcriptID is already in lastBorderPred...
# check that the start position is currently greater than the lastBorderPred of the transcriptID, otherwise error
# update the lastBorderPred to the stop column of the annotation
# add 1 to the elementNumberPred for the transcriptID, so it now becomes 2, 3, 4, etc...
# add an entry to transcript2ElementsPred for this transcript ID and elementNumberPred, equal to the new elementID (chr_start_stop_strand)
    print "## 2. prediction: assuming non-zipped, genomically ordered gff in the format Nicholas provided" > "/dev/stderr";
    print "## 2a. parsing"
    count=0;
    comm="zcat "pred;
    while(comm | getline){

	# only considering lines for elements of teh desired type (exon|intron)
	if($3!=keyword){continue;}
	#count++;
	#if(count>=100000){break;}
	# for summary stats (and sanity checks)
	predictedChroms[$1]++;
	predictedStrands[$7]++;
	
	transcriptID=$trIDcolumn;
	elementID=$1"_"$4"_"$5"_"$7;
	elementsPred[elementID]=1;

	if(transcriptID in lastBorderPred){
	    if($4 <= lastBorderPred[transcriptID]){
		print "ERROR: input not ordered: transcriptID="transcriptID"; lastBorderPred[transcriptID]=" lastBorderPred[transcriptID]"; line="$0 > "/dev/stderr";
		exit;
	    }
	    lastBorderPred[transcriptID]=$5;
	    elementNumberPred[transcriptID]++;
	    transcript2ElementsPred[transcriptID"\t"elementNumberPred[transcriptID]]=elementID;	   
	}
	else{
	    lastBorderPred[transcriptID]=$5;
	    elementNumberPred[transcriptID]=1;
	    transcript2ElementsPred[transcriptID"\t"elementNumberPred[transcriptID]]=elementID;	   	    
	}	 	           
    }
    
# for each transcript stored in lastBorderPred...
# the stretch starts as ""
# for each elementNumberPred (1, 2, 3....) for that transcript...
# the stretch becomes the elementID from transcript2ElementsPred of that transcript ID and elementNumberPred (ie. chr_start_stop_strand)
# n becomes the length of splitting "stretch" into its individual stretches
# deduct 1 from n
# if the stretch is found in annotatedStretch...
# print out that it is known
# otherwise print out that it is novel
    print "## 2b. are the RNAseq stretches annotated ?" > "/dev/stderr";
    # considering all transcripts  
    for(tr in lastBorderPred){	
	stretch="";
	for(j=1;j<=elementNumberPred[tr];j++){
	    stretch=stretch";%;"transcript2ElementsPred[tr"\t"j];
	}
	n=split(stretch,a,";%;");
	n--;
	if(stretch in annotatedStretch){
	    print n"\tknown\t"tr"\t"stretch > outfile;
	}
	else{
	    print n"\tnovel\t"tr"\t"stretch > outfile;
	}
    }
}

 
