#! /usr/bin/env bash 

data="/Users/rachel/chip-seq-analysis/data"
TFBS="$data/encode.tfbs.chr22.bed.gz"
H3="$data/encode.h3k4me3.hela.chr22.bed.gz"
HG19="$data/hg19.genome"
TSS="$data/tss.hg19.chr22.bed.gz"
GENES="$data/genes.hg19.bed.gz"
FACT="$data/factorx.hela.chr22.fq.gz"
CHR22="$data/hg19.chr22.fa"
CTCF="$data/ctcf.hela.chr22.bg.gz"

#1 largest overlap between CTCF and H3K4me3

answer_1=$(bedtools intersect -u -a $TFBS -b $H3 \
| awk '($4=="CTCF")' \
| awk 'BEGIN {OFS="\t"} {print $3-$2}'\
| sort -nr \
| head -1)

echo "answer-1: $answer_1"


#2 GC content of nts 19,000,000 to 19,000,500 on chr22 of hg19, as frac.

answer_2=$(echo -e "chr22\t19000000\t19000500" > interval.bed \
| bedtools nuc -fi $CHR22 -bed interval.bed \
| cut -f 5 \
| tail -1)

echo "answer-2: $answer_2"


#3 length of CTCF ChIP-seq peak with the largest mean signal

answer_3=$(bedtools map -b $CTCF -a $TFBS -c 4 -o mean \
| sort -k5nr \
|awk '($4=="CTCF")' \
|head -1 \
|awk '{OFS="\t"} {print $3-$2}')

echo "answer-3: $answer_3"


#4 identify the gene promoter with the highest median signal in ctcf.hela.chr22.bg.gz, reporting the gene name

answer_4=$(bedtools slop -i $TSS -g $HG19 -l 1000 -r 0 -s \
| sort \
| bedtools map -a - -b $CTCF -c 4 -o median \
| sort -k7nr \
| head -1 \
| cut -f 4)

echo "answer-4: $answer_4"

#5 identify the longest interval on chr22 not covered by genes.hg19.bed.gz, reporting the gene:interval.

answer_5=$(bedtools sort -i $GENES \
| bedtools complement -i - -g $HG19 \
| awk 'BEGIN {OFS="\t"} ($1=="chr22") {print $1, $2, $3, $3-$2}' \
| sort -k4nr \
| head -1 \
| awk '{print $1":"$2"-"$3}')

echo "answer-5: $answer_5"

#6 bedtools cluster to find number of clusters in $GENES on chr22

answer_6=$(bedtools cluster -s -i $GENES \
| sort -k7n \
| awk '($1=="chr22")' \
| sort -u -k7 \
| wc -l)

echo "answer-6: $answer_6"
