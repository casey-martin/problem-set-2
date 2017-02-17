#! /usr/bin bash

mydir='/home/dalgarno/Documents/GAW/data-sets'

myanswer1=$(bedtools intersect -a $mydir/bed/encode.tfbs.chr22.bed.gz \
	-b $mydir/bed/encode.h3k4me3.hela.chr22.bed.gz \
	| awk '($4 == "CTCF") {print $3-$2}' \
	| sort -nr \
	| head -n 1) 


echo answer-1: $myanswer1

#generate temporary bed file with interval specified. current script leaves bed file in tmp
#for inspection in case of debugging. uncomment the "rm" command to delete tmp file.

echo -e "chr22\t19000000\t19000500\tfoo\n" > /tmp/answer2.bed

myanswer2=$(bedtools nuc -fi $mydir/fasta/hg19.chr22.fa -bed /tmp/answer2.bed \
	| cut -f 6 \
	| grep -v "^[0-9]_")


echo answer-2: $myanswer2
#rm /tmp/answer2.bed


# is it more efficient to filter out non-CTCF peaks prior to bedtools 
# operations? if so, is speed up worth it? need to find out why 
# bedtools indexing is better than nested "for loops"


myanswer3=$(bedtools map -wo -o mean -c 4 -a $mydir/bed/encode.tfbs.chr22.bed.gz \
	 -b $mydir/bedtools/ctcf.hela.chr22.bg.gz \
	| awk '($4 == "CTCF")' \
	| sort -k5nr \
	| head -n 1 \
	| awk 'BEGIN {OFS = "\t"} {print $3-$2}')


echo answer-3: $myanswer3

# create 1000 bp interval upstream with respect to strandedness. -l denotes left interval
# and -r right. -s takes +/- orientation. this file must be bedtools sort-ed for map ops.
# save to tmp/tss.bed. uncomment rm to remove tmp file.
bedtools slop -i $mydir/bed/tss.hg19.chr22.bed.gz -g $mydir/genome/hg19.genome \
	-l 1000 -r 0 -s \
	| bedtools sort > /tmp/tss.bed
#rm /tmp/tss.bed

# find median values of peaks in tss intervals. 
myanswer4=$(bedtools map -o median -c 4 -a /tmp/tss.bed \
	-b $mydir/bedtools/ctcf.hela.chr22.bg.gz \
	| sort -k7nr \
	| head -n 1 \
	| cut -f 4)

echo answer-4: $myanswer4



myanswer5=$(bedtools complement -i /tmp/genes.hg19.bed.gz -g $mydir/genome/hg19.genome \
	| awk '($1 == "chr22") {print $0, $3 - $2}' \
	| sort -k4nr \
	| head -n 1 \
	| awk '{print $1":" $2 "-" $3}')

echo answer-5: $myanswer5


echo -e "\tBonus Round: Where a cluster is defined as overlapped or bookended intervals on 
\tthe same DNA strand, return the number of genes in the largest gene cluster in hg19."


mybonus=$(bedtools cluster -i $mydir/bed/genes.hg19.bed.gz -s \
	| cut -f 7 \
	| sort \
	| uniq -c \
	| sort -k1nr \
	| awk '(NR==1) {print $1}')

echo Bonus Answer: $mybonus

