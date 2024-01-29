
##### 3'REAP (3’ Reads Enrichment using Annotated PolyA sites) pipeline, using QuantSeq Pool read2 data as input. 

### by Luyang Wang, lwang@wistar.org


######### before running the code below, please create a data name list, DataNameList.txt, which contains only one column of each raw sequencing data name.
######### for example, if the raw sequencing data name is 12345-02-01-01_S168_L008_R1_001.fastq, the name should be put into the txt file should be "12345-02-01-01_S168_L008_R1_001".


#### use TrimGalore to trim adapter
mkdir TrimGalore
while IFS=$' \t\r\n' read -r sample; do
trim_galore ./rawfastq/${sample}.fastq -o ./TrimGalore/
done < DataNameList.txt
echo Time is `date`, TrimGalore is finished




### remove 5' most 15 nt 
mkdir clipped_15nt
while IFS=$' \t\r\n' read -r sample; do
cutadapt -u 15 ./TrimGalore/${sample}_trimmed.fq -o ./clipped_15nt/${sample}.clipped.fastq
done < DataNameList.txt
echo Time is `date`, cut_N is finished


### remove remaining 5'Ts 
python ./assistScripts/trim_5T.py --rawfastq_dir ./clipped_15nt --project_dir ./

echo Time is `date`, trim_5T.py is finished



#generate genome index, based on polyA_DB annotation, with window of -100nt, PAS, +25nt
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir ./star_index/hg19_polyA_DB_v4.1s_100.25/ --genomeFastaFiles ./assistScripts/hg19_polyA_DB_v4.1s_100.25.fasta --limitGenomeGenerateRAM 330000000000

echo Time is `date`, genome index is generated

### mapping
mkdir star_out

while IFS=$' \t\r\n' read -r sample; do
STAR --runThreadN 2 --genomeDir ./star_index/hg19_polyA_DB_v4.1s_100.25/ --readFilesIn ./clipped_15nt/${sample}.clipped.trimmed.fastq --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./star_out/${sample} --limitBAMsortRAM 2000000000
done < DataNameList.txt

echo Time is `date`, star_out is finished



# define LAP (last aligned position)
echo Time is `date`, running is started

# Working directory
WORK_DIR=.
# Input directory
INPUT_DIR=./star_out
# Output directory
OUTPUT_DIR=$WORK_DIR/LAP

## Create output directory
if [ ! -d $OUTPUT_DIR ]; then
  mkdir $OUTPUT_DIR
fi

cd $OUTPUT_DIR
while IFS=$' \t\r\n' read -r sample; do

bedtools bamtobed -cigar -i $INPUT_DIR/${sample}Aligned.sortedByCoord.out.bam > ${sample}.bed
sort -k 1,1 ${sample}.bed > ${sample}.bed.sorted
done < DataNameList.txt

wc -l *.bed.sorted

cd $WORK_DIR

echo Time is `date`, LAP is finished





mkdir result
mkdir result/csv

while IFS=$' \t\r\n' read -r sample; do
Rscript ./assistScripts/polyA_DB_match.R -bedLAP ./LAP/${sample}.bed.sorted -out ./result/csv/${sample}
done < DataNameList.txt
Rscript ./assistScripts/combine_all_sample_PAS_count_tables.R -csv ./result/csv -out ./result/cluster.all.reads.csv
echo Time is `date`, polyA_DB_match is finished



# generate UCSC genome browser bigwig files, based on LAP of each read.

# Working directory
WORK_DIR="."
# Input directory
INPUT_DIR=$WORK_DIR/result/csv
# Output directory
OUTPUT_DIR=$WORK_DIR/bigwig_LAP
# Chromosome sizes dir/file
chromsizes=hg19.chrom.sizes

## Create output directory
if [ ! -d $OUTPUT_DIR ]; then
  mkdir $OUTPUT_DIR
fi

cd $OUTPUT_DIR

while IFS=$' \t\r\n' read -r sample; do
  echo "Working on $sample..."
  ## Count total read number
  totalReadNum=`wc -l $INPUT_DIR/${sample}_PASS_bw.bed | sed s/[[:blank:]].*//`
  echo "step2, for file ${sample}_PASS_bw.bed, TotalReadNum=$totalReadNum"

  sort -k 1,1 $INPUT_DIR/${sample}_PASS_bw.bed > ./${sample}_PASS_bw.bed.sorted

  ## it is strand-specific 
  ## if reverse: + on genomeCoverageBed is "minus", and - is "plus"
  ## if forward: + on genomeCoverageBed is "plus", and - is "minus"
  ## Generate bedgraph file
  echo "Generate bedgraph files for + and - strands..."
  genomeCoverageBed -bg -split -i ./${sample}_PASS_bw.bed.sorted -strand '-' -g $chromsizes > $sample.plus.bedgraph
  genomeCoverageBed -bg -split -i ./${sample}_PASS_bw.bed.sorted -strand '+' -g $chromsizes > $sample.minus.bedgraph

  ## Normalize bedgraph counts
  echo "Normalize bedgraph counts..."
  norm_bedgraph.pl -t $totalReadNum -i "$sample.plus.bedgraph"
  norm_bedgraph.pl -t $totalReadNum -i "$sample.minus.bedgraph"

	## give minus strand negative value
  awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2, $3, -$4}' $sample.minus.bedgraph.normolized > $sample.minus.bedgraph.normolized1 

  ## Convert to bigwig file
  echo "Convert to bigwig file..."
  bedGraphToBigWig $sample.plus.bedgraph.normolized  $chromsizes $sample.plus.bw
  bedGraphToBigWig $sample.minus.bedgraph.normolized1  $chromsizes $sample.minus.bw

  echo Time is `date`
done < DataNameList.txt


cd $WORK_DIR

echo Time is `date`, bigwig_LAP is finished


