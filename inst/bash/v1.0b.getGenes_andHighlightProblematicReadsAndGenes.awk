# by Hagen Tilgner for mapping of pacBio reads : Jan 8th, 2011
# expects to get 
# - an annotation file sorted by position piped in
# -v sortedAnno=
# -v intronGFF=
# -v feature=
# -v transIdColumn=
# -v task=



function d(x,y){
    split(x,X,"_"); 
    split(y,Y,"_");

    if(X[1]!=Y[1]){
	print "ERROR:"X[1]"!="Y[1] > "/dev/stderr";
	exit(0);
    }
    if(X[4]!=Y[4]){
	print "ERROR:"X[4]"!="Y[4] > "/dev/stderr";
	exit(0);
    }    

    if(X[3]<Y[2]){
	return(Y[2]-X[3]-1);
    }
    if(Y[3]<X[2]){
	return(X[2]-Y[3]-1);
    }
    return(-1);
      
} 

# treating the annotation
BEGIN{

    print "### executing pB.getGenes.awk" > "/dev/stderr";
    if(!sortedAnno){print "ERROR: no value for sortedAnno" > "/dev/stderr";exit(0);}
    if(!intronGFF){print "ERROR: no value for intronGFF" > "/dev/stderr";exit(0);}
    if(!feature){print "ERROR: no value for feature" > "/dev/stderr";exit(0);}
    if(!transIdColumn){print "ERROR: no value for transIdColumn" > "/dev/stderr";exit(0);}
    if(!geneIdColumn){print "ERROR: no value for geneIdColumn" > "/dev/stderr";exit(0);}
    

    print "## A. parsing file2="intronGFF > "/dev/stderr"
    comm="cat "intronGFF;
    while(comm | getline){	

#get read id and set read2Gene value for
#read id to ""
	split($9,a,"transcript_id_with_chr=");
	split(a[2],b,"@");
	readID=b[1];	
	read2Gene[readID]="";

#l1 is start position (of the intron!) minus 1
#key1 is chromosome, start pos intron, strand
	l1=$4-1;
	key1=$1"_"l1"_"$7;
	
#if this is the first time the key1 is encountered
#set endSite2Read[key1] to the read id
#otherwise append the new read id to endSite2Read[key1]
# so "endSite" refers to the end of exon?

	if(key1 in endSite2Read){endSite2Read[key1]=endSite2Read[key1]";"readID}
	else{endSite2Read[key1]=readID;}
       
#l2 is stop position of intron(!) plus 1
#key2 is chromosome, stop pos intron, strand
	l2=$5+1;
	key2=$1"_"l2"_"$7;
	
#as above (create list of read_ids matching key2
	if(key2 in startSite2Read){startSite2Read[key2]=startSite2Read[key2]";"readID}
	else{startSite2Read[key2]=readID;}

    }


    


    print "## B. parsing annotation: " > "/dev/stderr"
    comm="cat "sortedAnno;
    while(comm | getline){
#if column 3 does not equal "exon", skip this section
#if lastColumn4 exists and current column4 is less than previous column 4, error and exit
	if(feature && $3!=feature){continue;}
        if(lastColumn4 && $4<lastColumn4){
            print "ERROR:lastColumn4="lastColumn4"; and $0="$0 > "/dv/stderr";
            exit(0);
        }
#if column 7 does not equal + or -, error and exit
	if($7!="+" && $7!="-"){
	    print "ERROR: cannot deal with strand="$7 > "/dev/stderr";
	    exit(0);
	}
#if the transcript id in column 12 is in strand and strand[$12] does not equal column 7
# call an error because the strands do not match
	if($transIdColumn in strand && strand[$transIdColumn]!=$7){
	    print "ERROR: strands do no match:"strand[$transIdColumn]"!="$7 > "/dev/stderr";
	    exit(0);
	}
# if the transcript id in column 12 is in chr and chr[$12] does not equal column 1, error and exit
	if($transIdColumn in chr && chr[$transIdColumn]!=$1){
	    print "ERROR: chroms do no match:"chr[$transIdColumn]"!="$1 > "/dev/stderr";
	    exit(0);
	}
# split column 10 (gene id) into variable G, removing "
# gene[current transcript id] equals the gene name
# add 1 to the count for the current transcript id
# exon[transcript id + count for transcript id] = chr_start_stop_strand
# ^ in other words each exon is uniquely identified by the chr_start_stop_strand
# strand[transcript id] = strand column 7
# chr[transcript id] = chr column 1
# lastColumn4 becomes current column4 (start position)
	split($geneIdColumn,G,"\"");
	gene[$transIdColumn]=G[2];
	n[$transIdColumn]++;
	exon[$transIdColumn"\t"n[$transIdColumn]]=$1"_"$4"_"$5"_"$7;	
	strand[$transIdColumn]=$7;
	chr[$transIdColumn]=$1;
	lastColumn4=$4;
    }
    
    print "## C. going over all annotated transcripts: " > "/dev/stderr";
# for every transcript id listed in the strand variable
    for(tr in strand){
# iterate over every instance of a transcript
# going over all exons of this transcript
	for(i=1;i<=n[tr];i++){
# split each exon for a transcript into "a", separating _
# keyStart becomes chrom_start_strand
# keyStop becomes chrom_stop_strand
	    split(exon[tr"\t"i],a,"_");
	    keyStart=a[1]"_"a[2]"_"a[4];
	    keyEnd=a[1]"_"a[3]"_"a[4];
	    
# if keyEnd is in endSite2Read and this is not the last exon of the transcript...
# and the transcript's gene "\t" keyEnd "\t" "end" is not in geneSpliceSitePair
# m becomes length of the split of endSite2Read for this keyEnd, splitting based on ";", into b
# endSite2Read was generated in part A: It contains a ";" separated list of read_ids for each key.
# for every part of b (for each read id)
# read2Gene[read_id] = gene[read_id], appending any new gene id
	    if(keyEnd in endSite2Read && i<n[tr] && !(gene[tr]"\t"keyEnd"\t""end" in geneSpliceSitePair) ){				
		m=split(endSite2Read[keyEnd],b,";");
		for(j=1;j<=m;j++){
		    read2Gene[b[j]]=read2Gene[b[j]]";"gene[tr];
		}
		#sitesEnd[keyEnd"\t"gene[tr]]=1;
	    }
	    if(i<n[tr]){
		geneSpliceSitePair[gene[tr]"\t"keyEnd"\t""end"]=1;
	    }
# if keyStart is in startSite2Read and...
# if the number of exons is greater than 1 and...
# if the gene_keyStart_start is not in geneSpliceSitePair

	    if(keyStart in startSite2Read && i>1 && !(gene[tr]"\t"keyStart"\t""start" in geneSpliceSitePair) ){
		m=split(startSite2Read[keyStart],b,";");
		for(j=1;j<=m;j++){
		    read2Gene[b[j]]=read2Gene[b[j]]";"gene[tr];
		}
		#sitesStart[keyStart"\t"gene[tr]]=1;
	    }
# if this is the second or greater exon
# geneSpliceSitePair for this gene at this keyStart becomes 1
	    if(i>1){
		geneSpliceSitePair[gene[tr]"\t"keyStart"\t""start"]=1;		
	    }
       			    
	}	    
    }
    print "## D. counting the number of splice sites per gene:" > "/dev/stderr";
#for every geneSpliceSitePair entry
#split the gene id, keyEnd and "end" based on "\t"
#if the length of the split does not equal 3, error and exit
#add 1 to the count of spliceSiteNumber for each unique gene id
#if the keyEnd and strand are not in spliceSite2Gene...
#add an entry for the keyEnd and strand and set the value to the gene id
#if the keyEnd and strand are in spliceSite2Gene...
#add the gene id to problematicGene and output its name to the problematicGeneOutFile
    for(k in geneSpliceSitePair){
	x=split(k,a,"\t");
	if(x!=3){
	    print "ERROR: "k" is not a valid key in geneSpliceSitePair" > "/dev/stderr";
	    exit(0);
	}
	spliceSiteNumber[a[1]]++;

	if(a[2]"\t"a[3] in spliceSite2Gene){
	    problematicGene[a[1]]=""; #print a[1] > problematicGeneOutFile;
	    problematicGene[spliceSite2Gene[a[2]"\t"a[3]]]=""; #print spliceSite2Gene[a[2]"\t"a[3]] > problematicGeneOutFile;
	}
	else{
	    spliceSite2Gene[a[2]"\t"a[3]]=a[1];
	}
    }
    

    print "## E. checking whether a read has equal numbers of splice sites with multiple genes: " > "/dev/stderr"    
#for every read in read2Gene
# delete all entries in h
# split the value of read2Gene[read id] based on ";"
# so separate each gene_id associated with a read (there should only be 1?)
# add 1 to the h[gene] for the number of second, third etc genes
# g = "none"
# v = -1
# set problemRead and problemGene to "fine" as default
# for every gene in h...
# if the count of that genes is greater than -1
# problemRead stays as fineRead
# if k is in problematicGene, 
# problematicGene is set as "problematicGene"
# otherwise keep problematicGene as "fineGene"
# set g to the current gene
# set v to the count for the current gene
# move onto the next gene
# perhaps variables are -1 by default?
# because then we check whether the count for a gene id equals -1
#
#
# finally prints out read, gene, whether it is a problemRead and whether it is a problemGene
    for(r in read2Gene){
	#if(r!="m121212_085615_00126_c100418772550000001523036412191280_s1_p0/47403/ccs.path1"){
	#    continue;
	#}
	for(k in h){
	    #print "deleting "k" from h"; 
	    delete h[k];
	}
	m=split(read2Gene[r],b,";");
	#print "m="m;
	for(i=2;i<=m;i++){
	    #print "i="i"; b[i]="b[i]
	    h[b[i]]++;
	}	
	g="none";
	v=-1;
	problemRead="fineRead";
	problemGene="fineGene";
	for(k in h){
	    if(h[k]>v){
		problemRead="fineRead";
		if(k in problematicGene){
		    problemGene="problematicGene";
		}
		else{
		    problemGene="fineGene";
		}
		g=k;
		v=h[k];
		continue;
	    }
	    if(h[k]==v){
		problemRead="problematicRead";
		if(k in problematicGene){
		    problemGene="problematicGene";
		}
		g=g"@;@"k;
		v=h[k];
	    }
	}	
	print r"\t"g"\t"problemRead"\t"problemGene;
    }



}


