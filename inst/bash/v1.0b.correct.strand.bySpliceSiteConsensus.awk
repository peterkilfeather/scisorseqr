# by Hagen Tilgner for mapping of pacBio reads : Jan 1st, 2011
# expects -v intronTypeFile=??? (teh file that gives the splice site consensus) -v bestMatchFile=??? (the mapping file) -v outputCorrectFile=??? -v outputUnClassifiableFile=???
function reverseStrand(x){
    n=split(x,a,"\t");
    if(a[7]!="+" && a[7]!="-"){print "STRAND-ERROR:"x > "/dev/stderr";}
    if(a[7]=="+"){strand="-";}else{strand="+";}
    s=a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]"\t"strand;
    for(i=8;i<=NF;i++){s=s"\t"$i;}
    return(s);
}

function getID(x){
    split($9,a,"transcript_id_with_chr=");
    split(a[2],b,"@");
    return(b[1]);
}

function getType(x){
    split(x,a,"=");
    return(a[2]);
}

# reads in the intron consensus file from earlier
# for each read (trID=getID($9))
#checks whether it has already been classified (is found in cStrand, wStrand
# or notClassifiable)
#sets the value of c/w/not/Strand for that read to
#zero if not
#gets the splice junction type (the 4 letter consensus)
#checks whether it belongs to the "correct" (cStrand)
#or "wrong" (wStrand)
#and adds a count to the trID for that group
BEGIN{
    comm="zcat "intronTypeFile;
    while(comm | getline){
	trID=getID($9);
	if(!(trID in cStrand)){cStrand[trID]=0;}
	if(!(trID in wStrand)){wStrand[trID]=0;}
	if(!(trID in notClassifiable)){notClassifiable[trID]=0;}
	t=getType($10);
	if(t=="GTAG" || t=="GCAG" || t=="ATAC"){cStrand[trID]++; continue;}
	if(t=="CTAC" || t=="CTGC" || t=="GTAT"){wStrand[trID]++; continue;}
	notClassifiable[trID]++;
    }

#reads the best mapped filtered gff file
#split the last column based on "\""
#pull the read id from the split
#if the read id is found in wStrand (which it should be
#if it was processed in the previous while loop)
#if it has correct strand counts and no
#incorrect strand counts or notClassifiable counts
#the line is printed to output
#otherwise, if it has wrong strand counts, no correct counts
#and not notClassifiable counts, it is subject to the
#reverseStrand function and then output.
#If the notClassifiable count is positive, the line is
#output to the outputUnClassifiableFile
    while(getline<bestMatchFile>0){
	split($10,a,"\"");
	readID=a[2];
	if(readID in wStrand){
	    if(cStrand[readID]>0 && (wStrand[readID]==0 && notClassifiable[readID]==0)){print $0 > outputCorrectFile ;continue;}
	    if(wStrand[readID]>0 && (cStrand[readID]==0 && notClassifiable[readID]==0)){print reverseStrand($0) > outputCorrectFile ;continue;}
	    if(notClassifiable[readID]>0){print $0 > outputUnClassifiableFile ;continue;}
	}
    }
}
