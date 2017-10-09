# cat Aligned.out.sam | awk -f extractSplit.awk | sort -V > split
BEGIN {
    isSplit = 0;
    readNameOld = "";
};
{
if (substr($1,1,1)!="@") {
    seqName = "";
    if(match($1, /seq\.[0-9]+/)) {
        seqName = substr($0, RSTART, RLENGTH);
    }
    if (seqName != readNameOld) {
       if(isSplit == 1) {
           SJ[chunk] = 1;
       }
       chunk = ""; 
       isSplit = 0;
    }
    readNameOld = seqName;
    if(chunk == "") {
        chunk = $0;
    } else {
        chunk = chunk "\n" $0;
    }
    if (match($6, /N/)) { 
        isSplit = 1;
    }
};
};

END {

for (ii in SJ) {
    print ii
};
};
