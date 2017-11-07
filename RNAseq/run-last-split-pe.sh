echo "***** [$0] start " `date +'%Y/%m/%d %H:%M:%S'` " *****"
cat preprocessed_t3r1.maf | ./last-split-pe/RNAseq/last-split-pe -f 7.4 -s 1.6 --sam > result_t3r1.sam
echo "***** [$0] end " `date +'%Y/%m/%d %H:%M:%S'` " *****"
