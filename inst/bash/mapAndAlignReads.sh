#!/bin/bash
# by Hagen Tilgner for mapping of pacBio reads : December 31st, 2011; updated July 2012
# Edited AJ 04/2020

{ echo "#################################";
  echo "# RUNNING [$0]";
  echo "# Current date: `date`";
  echo "# Current dir: `pwd`"; } 1>&2;


###############
# 0. the arguments
echo "+++++++++++++++++++++++++ 1. arguments";

fastqGuide=${1}
echo "fastqGuide="$fastqGuide

outdir=${2}
echo "outdir="$outdir

tmpDir=${3}
echo "tmpDir="$tmpDir

numThreads=${4}
echo "numThreads="$numThreads

seqDirectory=${5}
echo "seqDirectory="$seqDirectory

annoGZ=${6};
echo "annoGZ="$annoGZ

mappingBAMGuide=${7}
echo "mappingBAMGuide="$mappingBAMGuide

scriptDir=${8}
echo "scriptDir="$scriptDir

pyScriptDir=${9}

echo "++++++++++++++++++ 1b. deduced from arguments";
tmpdir1=$tmpDir"/tmp."$$
echo $tmpdir1
mkdir $tmpdir1
echo "tmpdir1="$tmpdir1


echo "++++++++++++++++++ 1c. input check";
if [ ! -f $tmpDir"/*/mapping.bestperRead.bam" ]
then
mappingBAM=$tmpdir1"/mapping.bam";
samtools merge $mappingBAM `cat $mappingBAMGuide`;

echo "+++++++++++++++++++++ 2b. sorting annotation and projection";

echo "+++++++++++++++++++++ 2b1 sorting annotation";
time zcat $annoGZ | awk '$3=="exon"' | sort -gk4,4 > $tmpdir1"/sortedAnno";

echo "++++++++++++++++++ 2c. how many were well mapped";
mappingWellMappedSAM=$tmpdir1"/mapping.bestperRead.sam";
mappingWellMappedBAM=$tmpdir1"/mapping.bestperRead.bam";

echo "+++++++++++++++ 2c1 finding best hits"
samtools view -H $mappingBAM > $mappingWellMappedSAM
#apply process described in https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.4259/MediaObjects/41587_2018_BFnbt4259_MOESM44_ESM.pdf
#awk to remove * aligning reads
#awk TO FILL OUT
for i in `cat $mappingBAMGuide`; do samtools view $i | awk '$3!="*"' | awk -v file=$i \
'BEGIN{u1=0;u2=0;d1=0;d2=0;}{n[$1]++; v[$1"\t"n[$1]]=$5; line[$1"\t"n[$1]]=$0;}END \
{for(k in n){if(n[k]==1){if(v[k"\t1"]>=20){print line[k"\t1"]; \
u1++; continue;}else{d1++; continue;}} \
max1=0; max2=0; argmax1=0; argmax2=0; for(i=1;i<=n[k];i++){if(v[k"\t"i]>=max1) \
{max2=max1; argmax2=argmax1; max1=v[k"\t"i];argmax1=i; continue;} \
if(v[k"\t"i]>=max2) {max2=v[k"\t"i]; argmax2=i; continue;}} \
if(max1>=20 && max1>=max2+20 && max2<=30){print line[k"\t"argmax1]; u2++;}else{d2++;}} \
x=(u1+u2)/(u1+u2+d1+d2);  print "wellMapped:\t"file"\t"u1"\t"u2"\t"d1"\t"d2"\t"x > "/dev/stderr";}';
done >> $mappingWellMappedSAM

echo "+++++++++++++++ 2c2 sam to bam conversion"
samtools view -b -S $mappingWellMappedSAM > $mappingWellMappedBAM

else
echo "---> Well mapped bam already exists, moving on to next step"
fi

echo "+++++++++++++++++++++ 2d. removing ribosomal RNA hits";
if [ ! -f $outdir"/mapping.bestperRead.noRiboRNA.bam" ]
then
echo "++++++++++++++++++ 2d1. finding ribosomal RNA hits";
#awk for rRNA annotations and prepare a ribosomal bed file with bedtools intersect
zcat $annoGZ | awk '$3=="exon" && ($14~/rRNA/ || $14~/rRNA_pseudogene/ || $14~/Mt_rRNA/)' > $tmpdir1"/tmp.AnnoRibo.gff";
bedtools intersect -abam $mappingWellMappedBAM -b $tmpdir1"/tmp.AnnoRibo.gff" -split -bed > $tmpdir1"/reads2Remove.bed";

echo "++++++++++++++++++ 2d2. removing ribosomal RNA hits";
#read bam file and reads2remove.bed and remove line if matches based on column 4/1
#export as bam
samtools view -h $mappingWellMappedBAM | awk -v hits=$tmpdir1"/reads2Remove.bed" \
'BEGIN{comm="cat "hits; while(comm|getline){rm[$4]=1;}}{if($1 in rm){next;} print;}' \
| samtools view -b -S - > $outdir"/mapping.bestperRead.noRiboRNA.bam";
else
echo "---> ribo-depleted bam already exists, moving on to next step"
fi

echo "++++++++++++++++++++++++ 3. further analysis";

echo "+++++++++++++++++++++ 3a. checking consensus";

if [ ! -f $outdir"/mapping.bestperRead.noRiboRNA.introns.gff.gz" ]
then
echo "++++++++++++++++++ 3a1. getting introns and exons in gff format";
#convert a bam into a bed file and if a read has a split alignment, represent
#each exon of alignment as a separate bed feature.
#start awk with tab output field sep
#if column 4 has already been seen, split column 11 into a variable called "blockLen"
#and split column 12 into a variable called "offset"
#n and m equal the length of blockLen and offset respectively
#if blockLen and offset are not equal (in length?) ERROR and exit
#if n < 2, then move on: There is no alignment split
#set the end of the exon to the start coordinate + the first blockLen value
#count up from i=2, but is less than or equal to the number of blockLen entries
#set the start of the exon to the start coordinate + offset + 1 because of bed format
#print our output of introns, by printing exEnd first (+1 for start of intron) and exStart second (-1 for intron end)
#set the exEnd to the end of the first exon for the next intron to be printed
time bedtools bamtobed -i $outdir"/mapping.bestperRead.noRiboRNA.bam" -bed12 | awk 'BEGIN{OFS="\t";}\
{seen[$4]++; n=split($11,blockLen,","); m=split($12,offset,","); \
if(m!=n){print "ERROR after bamtobed:m="m"!=n="n > "/dev/stderr";\
exit(0);} if(n<2){next;} exEnd=$2+blockLen[1];  for(i=2;i<=n;i++){ exStart=$2+offset[i]+1; \
print $1,"hg19","intron",exEnd+1,exStart-1,".",$6,".","transcript_id_with_chr="$4".path"seen[$4]"@"$1; \
exEnd=exStart+blockLen[i]-1;}}' > $outdir"/mapping.bestperRead.noRiboRNA.introns.gff"

#similar to above, except using normal exon reporting
time bedtools bamtobed -i $outdir"/mapping.bestperRead.noRiboRNA.bam" -bed12 | awk 'BEGIN{OFS="\t";}\
{seen[$4]++; n=split($11,blockLen,","); m=split($12,offset,","); \
if(m!=n){print "ERROR after bamtobed:m="m"!=n="n > "/dev/stderr"; exit(0);} \
for(i=1;i<=n;i++){ exStart=$2+offset[i]+1; exEnd=exStart+blockLen[i]-1; \
print $1,"hg19","cDNA_match",exStart,exEnd,".",$6,".","read_id \""$4".path"seen[$4]"\";"; }}' \
 > $outdir"/mapping.bestperRead.noRiboRNA.gff"

else
echo "---> intron.gff file already exists, moving on"
fi


echo "++++++++++++++++++ 3a2. getting dinucleotides at intron borders";
echo "++++++++++++++++++ 3a2.a. preparing commands for parallel execution";
if [ ! -f $outdir"/parallel.comp.anno.guide.3.siteSeq" ]
then

#for each chromosome in the introns gff file,
#print line if starts with chr, otherwise prepend chr and print line
#output seqfile is the chromosome fasta gzipped
#calls a f2tFunc.sh script to get fasta sequence names
#and a v0.1.getSpliceSiteSequence.awk script to get consensus splice sequences
#output to a site_sequence file specific to each chromosome
#echo the str command to a parallel document
#define collectCommand and rmCommand
#shuffle parallel commands

collectCommand="cat "
rmCommand="rm "
for i in `cat $outdir"/mapping.bestperRead.noRiboRNA.introns.gff" | awk '{if(!($1 in s))print $1;s[$1]=1;}' | sort -u | \
awk '{if($1~/chr|\./) {print } else {print "chr"$1}}'`; do
    seqfile=$seqDirectory"/"$i".fa.gz";
    f2tFunc=$scriptDir/"f2t.sh"
    echo $seqfile >&2;
    str="source "$f2tFunc"; zcat "$seqfile" | FastaToTbl | awk -v chr="$i" \
-v file="$outdir"/mapping.bestperRead.noRiboRNA.introns.gff -f \
"$scriptDir"/v0.1.getSpliceSiteSequence.awk > \
"$tmpdir1"/site_sequence."$i
    echo $str >> $outdir"/parallel.comp.anno.guide.3.siteSeq";
    collectCommand=$collectCommand" "$tmpdir1"/site_sequence."$i
    rmCommand=$rmCommand" "$tmpdir1"/site_sequence."$i
done
cat $outdir"/parallel.comp.anno.guide.3.siteSeq" | shuf > $outdir"/tmp";
mv $outdir"/tmp" $outdir"/parallel.comp.anno.guide.3.siteSeq"

echo "+++++++++++++++ 3a2b execution ";
#execute commands in parallel
#changed to GNU parallel: Couldn't get the python script to work
time cat $outdir"/parallel.comp.anno.guide.3.siteSeq" | parallel

echo "+++++++++++++++ 3a2c collecting results and removing temporary files";
#collate sequences
echo $collectCommand;
echo $rmCommand;
`$collectCommand > $outdir"/mapping.bestperRead.RNAdirection.introns.type.gff"`;
`$rmCommand`;

time gzip $outdir"/mapping.bestperRead.RNAdirection.introns.type.gff"

else
echo "---> parallel process already done, moving on"
fi

echo "+++++++++++++++++++++ 3b. strand correction and follow up analysis";

echo "++++++++++++++++++ 3b1. getting mapping.bestperRead.noRiboRNA.bam in gff format ";

if [ ! -f $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz" ]
then
#use exon gff,
#consensus intron site sequence gff
#set output consensus introns gff filename
#set unclassified output name
#run v1.0b.correct.strand.bySpliceSiteConsensus.awk script
#to generate consensus split mapped molecules
time awk -v bestMatchFile=$outdir"/mapping.bestperRead.noRiboRNA.gff" \
-v intronTypeFile=$outdir"/mapping.bestperRead.RNAdirection.introns.type.gff.gz" \
-v outputCorrectFile=$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff" \
-v outputUnClassifiableFile=$outdir"/mapping.bestperRead.RNAdirection.withWeirdIntron.gff" \
-f $scriptDir/v1.0b.correct.strand.bySpliceSiteConsensus.awk

time gzip $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff"

#remove path and \" from read_id column
#make a list of read_ids
# add all PG, HD, SQ lines to "use" (status = 1)
# then filter bam file for read_id being present in use
# export as bam
samtools view -h $outdir"/mapping.bestperRead.noRiboRNA.bam" | \
awk -v CSMM=$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz" \
'BEGIN{comm="zcat "CSMM; while(comm|getline){split($10,a,/(\.path|\")/); \
use[a[2]]=1;} use["@PG"]=1; use["@HD"]=1; use["@SQ"]=1;}{if(!($1 in use)){next;} print;}' \
| samtools view -b -S - > \
$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.bam";

echo "++++++++++++++++++ 3b3. getting introns";

# This replaces the strand sign with the strandsign of the consensusIntrons gff
# strip \" from read_id
# store the strand sign (+/-) of each read_id
# in the "introns" gff, extract the read_id
# if the read_id is in use,
# print the introns.gff with the consensusIntrons strand sign and any remaining fields
time zcat $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz" | \
awk -v introns=$outdir"/mapping.bestperRead.noRiboRNA.introns.gff" \
'{split($10,a,"\""); use[a[2]]=$7;}END{comm="cat "introns; \
while(comm|getline){split($9,a,/(transcript_id_with_chr=|\@)/); \
if(a[2] in use){s=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"use[a[2]]; \
for(i=8;i<=NF;i++){s=s"\t"$i;} print s;}}}' > \
$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.introns.gff";

#There is a bug here requiring each chromosome to be split (from mapping.intron.gff file and sortedAnno)
#And each split mapping file to be split into 100,000 rows at a time and processed
#the splitter_script.sh will be uploaded and called here instead of these lines below
echo "++++++++++++++++++ 3b4. getting genes";
time awk -v sortedAnno=$tmpdir1"/sortedAnno" -v \
intronGFF=$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.introns.gff" \
-v feature=exon -v transIdColumn=12 -v geneIdColumn=10 \
-f $scriptDir"/"v1.0b.getGenes_andHighlightProblematicReadsAndGenes.awk \
| gzip -c > $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz";

echo "+++++++++++++++++++++++++ 4. zipping what has not be zipped before";
for i in `ls $outdir"/"mapping* | perl -ne 'print unless /\.gz|bam/'`; do gzip $i; done

else
echo "---> transcriptWise.genes.gz already exists, skipping this step and moving on"
fi


echo "+++++++++++++++++++++++++ 5. new isoforms";
if [ ! -d $outdir"/newIsoforms_vs_Anno_ignoreAnno" ]
then
	mkdir $outdir"/newIsoforms_vs_Anno_ignoreAnno"
fi

if [ ! -f $outdir"/newIsoforms_vs_Anno_ignoreAnno/stretches.gz" ]
then
echo "++++++++++++++++++++++ 5a. annotation";
time zcat $annoGZ | awk -v feature=exon -v transcol=12 -f $scriptDir"/"exon.gtf.to.intron.gff3.awk | \
sort -k1,1 -gk4,4 > $tmpdir1"/anno.introns.gff"
gzip $tmpdir1"/anno.introns.gff"


echo "++++++++++++++++++++++ 5b. RNAseq";
time zcat $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.introns.gff.gz" | \
sort -k1,1 -gk4,4 | gzip -c > $tmpdir1"/rnaSeq.introns.gff.gz"

# prints the exons covered by each read?
time zcat $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz" | \
awk '{exon=$1"_"$4"_"$5"_"$7; if(!a[$10]) {a[$10]=$10"\t;%;"exon} else {a[$10]=a[$10]";%;"exon} }END \
{for(i in a) {print a[i]}}' | \
gzip -c > $outdir"/newIsoforms_vs_Anno_ignoreAnno/exonStretches.gz"

echo "++++++++++++++++++++++ 5c. comparison";
awk -v keyword=intron -v anno=$tmpdir1"/anno.introns.gff.gz" -v pred=$tmpdir1"/rnaSeq.introns.gff.gz" \
-v outfile=$outdir"/newIsoforms_vs_Anno_ignoreAnno/stretches" -v trIDcolumn=9 \
-f $scriptDir/v0.1.findNotAnnotatedStretches.awk

gzip $outdir"/newIsoforms_vs_Anno_ignoreAnno/stretches"

else
echo "---> stretches file already exists, finishing up."
fi



echo "+++++++++++++++++++++++++ 6. done";
mv $tmpdir1"/mapping.bestperRead.bam" $outdir
rm $tmpdir1"/"*;
rmdir $tmpdir1;
