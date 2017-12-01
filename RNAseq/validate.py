import sys
import re

cigFile = open(sys.argv[1])
samFile = open(sys.argv[2])

while True:
    samLine = samFile.readline()
    if samLine[0] != '@':
        break
samLine = samLine.rstrip().split('\t')
samSeqName = samLine[0]
total = 0
one = 0
both = 0
noreport = 0
for cigLine in cigFile:
    cigLine = cigLine.rstrip().split('\t')
    cigSeqName = cigLine[0]
    cigRefName = cigLine[1]
    cigCigar = cigLine[4]
    if samSeqName != cigSeqName:
        isTooShort = False
        if re.match(r"^\d+M\d+N\d+M$", cigCigar):
            for match in re.findall(r"\d+[A-Z]", cigCigar):
                length = int(re.search(r"\d+", match).group())
                code = match[-1]
                if code == 'M' and length < 20:
                    isTooShort = True
            if not isTooShort:
                total += 1
                noreport += 1
                #print(cigSeqName, "is not reported.")
        continue
    assert(samSeqName == cigSeqName)
    # if true CIGAR is split, i.e. *M*N*M
    isTooShortMatch = False
    if re.match(r"^\d+M\d+N\d+M$", cigCigar):
        readLength = len(cigLine[-1])
        trueAlignInfo = [None] * readLength
        isForward = (cigLine[-2] == '+')
        consumeDir = 1 if isForward else -1
        refPos = int(cigLine[2])
        readPos = 0 if isForward else readLength-1
        for match in re.findall(r"\d+[A-Z]", cigCigar):
            length = int(re.search(r"\d+", match).group())
            code = match[-1]
            strand = cigLine[-2]
            if code == 'M':
                if length < 20:
                    isTooShortMatch = True
                for i in range(length):
                    trueAlignInfo[readPos] = (cigRefName, refPos, strand)
                    readPos += consumeDir
                    refPos += 1
            elif code == 'N':
                refPos += length

        predAlignInfo = [None] * len(cigLine[-1])
        while True:
            samCigar = samLine[5]
            refName = samLine[2]
            refPos = int(samLine[3])
            flag = int(samLine[1])
            isForward = ((flag & 0x10) == 0)
            consumeDir = 1 if isForward else -1
            refPos = int(samLine[3])
            readPos = 0 if isForward else readLength-1
            strand = '+' if isForward else '-'
            for match in re.findall(r"\d+[A-Z]", samCigar):
                length = int(re.search(r"\d+", match).group())
                code = match[-1]
                strand = cigLine[-2]
                if code == 'M' or code == '=' or code == 'X':
                    for i in range(length):
                        predAlignInfo[readPos] = (refName, refPos, strand)
                        refPos += 1
                        readPos += consumeDir
                elif code == 'H' or code == 'I' or code == 'S':
                    readPos += length * consumeDir
                elif code == 'N' or code == 'D':
                    refPos += length * consumeDir
            samLine = samFile.readline().rstrip().split('\t')
            samSeqName = samLine[0]
            if samSeqName != cigSeqName:
                break
        #for i in range(len(trueAlignInfo)):
        #    print(trueAlignInfo[i], predAlignInfo[i])
        idx = 0
        isCheckLeft = True
        isMatchLeft = False
        isMatchRight = False
        isForward = (cigLine[-2] == '+')
        for match in re.findall(r"\d+[A-Z]", cigCigar):
            length = int(re.search(r"\d+", match).group())
            code = match[-1]
            strand = cigLine[-2]
            if code == 'M':
                if isCheckLeft:
                    isCheckLeft = False
                    if isForward:
                        startRight = length
                        for i in range(length):
                            if predAlignInfo[i] == trueAlignInfo[i]:
                               isMatchLeft = True
                               break
                    else:
                        startRight = readLength-length-1
                        for i in range(readLength-1, readLength-length-1, -1):
                            if predAlignInfo[i] == trueAlignInfo[i]:
                               isMatchLeft = True
                               break
                else:
                    if isForward:
                        for i in range(startRight, startRight+length):
                            if predAlignInfo[i] == trueAlignInfo[i]:
                               isMatchRight = True
                               break
                    else:
                        for i in range(startRight, startRight-length, -1):
                            if predAlignInfo[i] == trueAlignInfo[i]:
                               isMatchRight = True
                               break
        if not isTooShortMatch:
            #print(cigSeqName, isMatchLeft, isMatchRight)
            total += 1
            if isMatchLeft and isMatchRight:
                both += 1
            else:
                if isMatchLeft or isMatchRight:
                    one += 1
#            print(cigSeqName, isMatchLeft, isMatchRight)
    else:
        while True:
            samLine = samFile.readline().rstrip().split('\t')
            samSeqName = samLine[0]
            if samSeqName != cigSeqName:
                break
print("total:{}, both correct: {}, one correct: {}, noreport: {}".format(total, both, one, noreport))
cigFile.close()
