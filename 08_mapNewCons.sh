#!/usr/bin/sh

#$1 full path to corrected backbone
#$2 full path to cleaned paired ends reads 1
#$3 full path to cleaned paired ends reads 2
#$4 number of threads

## MAP WITH BOWTIE2 ##
mkdir BWT_INDEX_post;
cd BWT_INDEX_post;
bowtie2-build --quiet $1 backbone_post;
cd ..;

echo "  -> Started mapping PE reads with bowtie - "$(date +"%D %T");
bowtie2 -q --phred33 --very-sensitive-local -I 100 -X 600 --fr --no-mixed --no-discordant --no-dovetail --no-unal -p $4 --reorder -x BWT_INDEX_post/backbone_post -1 $2 -2 $3 -S readsGood_PE_vs_backbone_post.sam --un-conc unmapped_post;
echo "  -> Finished mapping PE reads with bowtie - "$(date +"%D %T");
echo "";
echo "-> Finished handling reads - "$(date +"%D %T");


