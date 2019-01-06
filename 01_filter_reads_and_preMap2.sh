
echo "-> Started handling reads - "$(date +"%D %T");
echo "";
## FILTER BY ARTIFACTS AND LEFT-TAIL QUALITY TRIM ##
echo "  -> Started filtering reads in " $1 " and "  $2 " and "  $3 " and "  $4 " - "$(date +"%D %T");

cat $1 $2 | fastx_artifacts_filter -Q33 -v | fastq_quality_trimmer -Q33 -t 28 -l 50 -v > side_1.fastq;
cat $3 $4 | fastx_artifacts_filter -Q33 -v | fastq_quality_trimmer -Q33 -t 28 -l 50 -v > side_2.fastq;

echo "  -> Finished filtering reads in " $1 " and "  $2 " and "  $3 " and "  $4 " - "$(date +"%D %T");
echo "";

## FILTER BY DUPLICATES, UNKNOWN NUCLEOTIDES AND LEFT-TAIL QUALITY TRIM ##
echo "  -> Started deduplicating reads with prinseq - "$(date +"%D %T");

prinseq-lite.pl -ns_max_n 0 -trim_qual_left 28 -derep 15 -min_len 50 -out_format 3 -out_good readsGood -out_bad null -no_qual_header -log filter_prinseq.log -fastq side_1.fastq -fastq2 side_2.fastq;

echo "  -> Finished deduplicating reads with prinseq - "$(date +"%D %T");
echo "";
echo "  -> Cleaning up and merging read files - "$(date +"%D %T");
rm side_1.fastq side_2.fastq;
## MERGE BOTH SINGLE-ENDS FILES INTO ONE ##

cat readsGood_1_singletons.fastq readsGood_2_singletons.fastq > readsGood_SE.fastq;
rm readsGood_1_singletons.fastq readsGood_2_singletons.fastq;

echo "  -> Finised cleaning up and merging read files - "$(date +"%D %T");
echo "";

## MAP WITH BOWTIE2 ##

mkdir BWT_INDEX;
cd BWT_INDEX;
bowtie2-build --quiet $5 backbone_pre;
cd ..;
#echo "  -> Started mapping SE reads with bowtie - "$(date +"%D %T");
#bowtie2 -q --phred33 --very-sensitive-local --no-unal -p $6 --reorder -x BWT_INDEX/backbone_pre -U readsGood_SE.fastq -S readsGood_SE_vs_backbone_pre.sam;
#echo "  -> Finished mapping SE reads with bowtie - "$(date +"%D %T");
#echo "";
echo "  -> Started mapping PE reads with bowtie - "$(date +"%D %T");
bowtie2 -q --phred33 --very-sensitive-local -I 100 -X 600 --fr --no-mixed --no-discordant --no-dovetail --no-unal -p $6 --reorder -x BWT_INDEX/backbone_pre -1 readsGood_1.fastq -2 readsGood_2.fastq -S readsGood_PE_vs_backbone_pre.sam;
echo "  -> Finished mapping PE reads with bowtie - "$(date +"%D %T");
echo "";
#echo "  -> Started mapping MIXED reads with bowtie - "$(date +"%D %T");
#bowtie2 -q --phred33 --very-sensitive-local --no-discordant --no-dovetail --no-unal -p $6 --reorder -x BWT_INDEX/backbone_pre -U 'readsGood_SE.fastq,readsGood_1.fastq,readsGood_2.fastq' -S readsGood_MIX_vs_backbone_pre.sam;
#echo "  -> Finished mapping MIXED reads with bowtie - "$(date +"%D %T");
#echo "";

## MAP SINGLE AND PAIRED-ENDS READS TO CONSENSUS ##
#perl ~/software/GapFiller/GapFiller.pl -l libFile.txt -s tev_gfp_gap.fasta -o 50 -d 50 -t 0 -g 1 -T 3 -b fillTest_1

echo "-> Finished handling reads - "$(date +"%D %T");
