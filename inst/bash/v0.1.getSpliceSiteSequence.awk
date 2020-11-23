#convert FASTA sequence to reverse complement
function comp(x){
    x=toupper(x) ;
    if(x=="A"){return("T");}
    if(x=="T"){return("A");} 
    if(x=="C"){return("G");} 
    if(x=="G"){return("C");} 
    if(x=="N"){return("N");}   
} 
function revcomp(x){
    y=""; 
    for(i=length(x);i>=1;i--){y=y""comp(substr(x,i,1));} 
    return(y);
}

#seq is everything after the first column
{seq=$2;}

#read formatted gff
#check that column 1 equals the variable "chr"
#if it does not, move to next line
#create a key based on the gff columns 1, 4, 5 and 7
#if the key is in "c", go to next line (consensus already found?)
#print the line with site_consensus equal to the c[key] entry and
#site_consensus_ext equal to the c2[key] entry
#define compare as the fasta sequence from gff column 4, two nucleotides
#to column 5 minus 1, two nucleotides
#define compare2 as column 4 minus 5, for 12 nucleotides
#and column 5 minus 1 minus 5, for 12 nucleotides (wider)
#with "NNNNN" inserted between
#if reverse stranded, create reverse complement of
#compare
#print whole line of gff with site_consensus as compare
#and site consensus ext as compare2
#set c[key] to compare
#set c2[key] to compare2
#go to next line
END{
    while(getline<file>0){
	if($1!=chr){continue;} 
	key=$1"_"$4"_"$5"_"$7; 
	if(key in c){ 
	    print $0"\tsite_consensus="c[key]"\tsite_consensus_ext="c2[key]; 
	    continue;
	} 
	compare=toupper(substr(seq,$4,2))""toupper(substr(seq,$5-1,2)); 
	compare2=toupper(substr(seq,$4-5,2+10))"NNNNN"toupper(substr(seq,$5-1-5,2+10));  
	if($7=="-"){
	    compare=revcomp(compare);
	} 
	print $0"\tsite_consensus="compare"\tsite_consensus_ext="compare2; 
	c[key]=compare; 
	c2[key]=compare2;   
    }
}
