# goal: remove polyA-tails and remove revcomp sequence to RNA if necessary
# by Hagen;

# exampel of ccs read data: 
#nuvol:~ htilgner$ gzcat /Users/htilgner/data/trios/input/pacBioIsoforms/v1/yorubian/Gm19238/round1/ccs/fa/m120413_003026_00126_c100277812550000001523007807041240_s2_p0.ccs.fasta.gz | head -2
#>m120413_003026_00126_c100277812550000001523007807041240_s2_p0/57/ccs
#TTTTTTTTGAAGGTTCTCAGGTCTTTATTTGCTCTCTCAACTTCCAGGAATTGACTTATTTAATTAATCC

# exampel of subread data
#nuvol:~ htilgner$ gzcat /Users/htilgner/data/trios/input/pacBioIsoforms/v1/yorubian/Gm19238/round1/subreads/fa/GM19238_poly-AcDNA_opt_smrtanalysis_common_jobs_016_016525_data_filtered_subreads.fasta.gz| head -2
#>m120527_035935_00126_c100333442550000001523018909161283_s1_p0/7/0_334
#TTTGCGACCTTTCGCCACGACGTCAGCGCGCTTTCCCGGGCAGACGCCTC


function max(a,b){
    if(a>b){return(a);}
    return(b);
}
function min(a,b){
    if(a<b){return(a);}
    return(b);
}
function allOptionsThere(){  
    if(!exonMapping){
	print "ERROR: -v exonMapping=<string> must be provided; exiting" > "/dev/stderr";
	exit;
    }
    if(!readsVsGenes){
	print "ERROR: -v readsVsGenes=<string> must be provided; exiting" > "/dev/stderr";
	exit;
    }
    if(!anno){
	print "ERROR: -v anno=<string> must be provided; exiting" > "/dev/stderr";
	exit;
    }
    # from when on do we actually use
    if(!(cageCutoff)){
	print "ERROR: -v cageCutoff must be provided;  exiting" > "/dev/stderr";
	exit;	
    }
    if(!(unzipCommand)){
	print "ERROR: -v unzipCommand must be provided;  exiting" > "/dev/stderr";
	exit;	
    }

}

function comp(x){
    x=toupper(x) ;
    if(x=="A"){return("T");}
    if(x=="T"){return("A");} 
    if(x=="C"){return("G");} 
    if(x=="G"){return("C");} 
    if(x=="N"){return("N");}   
    
    print "ERROR:"x" is not a vlaid symbol" > "/dev/stderr";
    exit(0);
} 

function revcomp(x){
    y=""; 
    for(j=length(x);j>=1;j--){
	y=y""comp(substr(x,j,1));
    } 
    return(y);
}

function abs(x){
    if(x<0){x=(-1)*x;} 
    return(x);
}

# Read the cage bed file
# if the peak score is greater than or equal to 50
# xP1 = peak start site + 1
# for i in xP1 to peak stop site...
# TSS[chrom_baseCoordinate_strand]=chrom_xP1_peakStopSite_strand
#
# next start reading read2Gene file
# chuck out any problem or "none" reads
# Concatenate any gene_ids that belong to the same read (using ;)

BEGIN{

    # need to have the required options
    allOptionsThere();
    OFS="\t";
    if(verbose>=1)
	print "# 1. reading annotation" > "/dev/stderr";

    comm=unzipCommand" "anno;
    while(comm|getline){
	
	if($5>=cageCutoff){
	    xP1=$2+1;
	    for(i=$2+1;i<=$3;i++)
		TSS[$1"_"i"_"$6]=$1"_"xP1"_"$3"_"$6;	
	}
    } 

    if(verbose>=1)
	print "# 2. reading readsVsGenes" > "/dev/stderr";
  
    comm=unzipCommand" "readsVsGenes; 
    while(comm | getline){
	if($2=="none" || $3~/problem/ || $4~/problem/){skippedReads[$1]=1; continue;}
	read2Gene[$1]="\""$2"\";";
    }
#
# next read exon mapping file (lists the exon coordinates of reads mapped)
# remove the " from read_id column
# r equals the read_id
# if r is not in read2Gene and not in skippedReads, print an error and exit
# if r is in read2ReadEnd...
# if on the positive strand...
# read2ReadEnd[r] equals the minimum of either the current read2ReadEnd[r] or the start site of exon mapping
# if on the negative strand...
# read2ReadEnd[r] equals the max of either the current read2ReadEnd[r] or the stop site of exon mapping
# ...if r is not in read2ReadEnd...
# if on the positive strand...
# read2ReadEnd[r] equals the start site of exon mapping
# if on the negative strand...
# read2ReadEnd[r] equals the stop site of exon mapping
# readStrand[r] equals the current strand
# readChrom[r] equals the current chrom

    if(verbose>=1)
	print "# 3. reading exonMapping" > "/dev/stderr";

    comm=unzipCommand" "exonMapping;
    while(comm|getline){
	split($10,a,"\""); 
	r=a[2]; 
	if(!(r in read2Gene)){if(!(r in skippedReads)){print "ERROR"r > "/dev/stderr"; exit;}else{continue;}}    
	if(r in read2ReadEnd){
	    if($7=="+"){read2ReadEnd[r]=min(read2ReadEnd[r],$4);}
	    if($7=="-"){read2ReadEnd[r]=max(read2ReadEnd[r],$5);}   
	}
	else{
	    if($7=="+"){read2ReadEnd[r]=$4;}
	    if($7=="-"){read2ReadEnd[r]=$5;}
	    readStrand[r]=$7;
	    readChrom[r]=$1;
	} 
    }

    if(verbose>=1)
	print "# 4. finding the closest TSS" > "/dev/stderr";
    
    # readCounter equals 0
    # for r in read2ReadEnd...
    # add 1 to the readCounter
    # if the read counter modulo 1000000 == 0, and verbose is on, print the readCounter value to stderr
    # if a read is not in read2Gene but IS in readStrand, create and error and exit
    # minDDownstream, minDDownstreamTSS, minDUpstream, minDUpstreamTSS are set up
    # TSS["none--"] equals "none--"
    # for read end position minus 50 to read end position plus 50...
    # s equals readChrom[r]_basePosition_readStrand[r]
    # if s in TSS...
    # if read is on the positive strand and basePosition is before read START and...
    # the absolute difference between current base position and read START is...
    # less than the current minDUpstream value...
    # minDUpstreamTSS becomes equal to s
    # the same approach is applied to negative strand reads, except comparing to read START upstream... 
    # (so in the 5' direction along the - strand)
    # the reverse is applied to find downstream TSS
    
    readCounter=0;
    print "#readID\tminDUpstream\tminDUpstreamTSS\tminDUpstream\tminDUpstreamTSS\tGene\treadChrom\treadStrand"; 
    for(r in read2ReadEnd){
	readCounter++;
	if(readCounter % 1000000 == 0 && verbose>=2){
	    print readCounter > "/dev/stderr";
	}
	if(!(r in read2Gene && r in readStrand)){print "ERROR2" > "/dev/stderr"; exit;}
	minDDownstream=1000000000;
	minDDownstreamTSS="none--";

	minDUpstream=1000000000;
	minDUpstreamTSS="none--";

	TSS["none--"]="none--";
	#print r"\t"read2ReadEnd[r];
	#for(s in TSS){
	#    print "TSS\t"s;
	#}
	
	for(i=read2ReadEnd[r]-50;i<=read2ReadEnd[r]+50;i++){
	    s=readChrom[r]"_"i"_"readStrand[r];
	    if(s in TSS){
		
		# closest upstream
		if(readStrand[r]=="+" && i<=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDUpstream){
		    minDUpstream=abs(read2ReadEnd[r]-i);
		    minDUpstreamTSS=s;			     
		}
		if(readStrand[r]=="-" && i>=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDUpstream){
		    minDUpstream=abs(read2ReadEnd[r]-i);
		    minDUpstreamTSS=s;			     
		}
		# closest downstream
		if(readStrand[r]=="+" && i>=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDDownstream){
		    minDDownstream=abs(read2ReadEnd[r]-i);
		    minDDownstreamTSS=s;	
		}		
		if(readStrand[r]=="-" && i<=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDDownstream){
		    minDDownstream=abs(read2ReadEnd[r]-i);
		    minDDownstreamTSS=s;	
		}

	    }	    
	}
	print r"\t"minDUpstream"\t"TSS[minDUpstreamTSS]"\t"minDDownstream"\t"TSS[minDDownstreamTSS]"\t"read2Gene[r]"\t"readChrom[r]"\t"readStrand[r];          
    }
}

 
