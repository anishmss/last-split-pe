echo "***** [$0] start " `date +'%Y/%m/%d %H:%M:%S'` " *****"
lastal -r6 -q18 -a21 -b9 -d96 -e120 -Q0 -i1 db interleaved_t3r1.fa | last-split -g db -m 0.9 -n > preprocessed_t3r1.maf 
echo "***** [$0] end " `date +'%Y/%m/%d %H:%M:%S'` " *****"
