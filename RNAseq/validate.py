import sys
import re

cigFile = open(sys.argv[1])
samFile = open(sys.argv[2])

while True:
    samLine = samFile.readline()
    if samLine[0] != '@':
        break
samLine = samLine.split('\t')
samSeqName = samLine[0]
total = 0
right = 0
left = 0
both = 0
for cigLine in cigFile:
    cigLine = cigLine.split('\t')
    cigSeqName = cigLine[0]
    cigRefName = cigLine[1]
    cigCigar = cigLine[4]
    if samSeqName != cigSeqName:
        continue
    assert(samSeqName == cigSeqName)
    # if true CIGAR is split, i.e. *M*N*M
    if re.match(r"^\d+M\d+N\d+M$", cigCigar):
        trueAlignInfo = []
        currentPos = int(cigLine[2])
        for match in re.findall(r"\d+[A-Z]", cigCigar):
            length = int(re.search(r"\d+", match).group())
            code = match[-1]
            strand = cigLine[-2]
            if code == 'M':
                for i in range(length):
                    trueAlignInfo.append((cigRefName, currentPos, strand))
                    currentPos += 1
            elif code == 'N':
                currentPos += length
        predAlignInfo = [None] * len(cigLine[-1])
        while True:
            samCigar = samLine[5]
            refName = samLine[2]
            refPos = int(samLine[3])
            seqPos = 0
            strand = '?'
            for match in re.findall(r"\d+[A-Z]", samCigar):
                length = int(re.search(r"\d+", match).group())
                code = match[-1]
                strand = cigLine[-2]
                if code == 'M':
                    for i in range(length):
                        predAlignInfo[seqPos] = (refName, refPos, strand)
                        refPos += 1
                        seqPos += 1
                elif code == 'N':
                    refPos += length
                elif code == 'H' or code == 'S':
                    seqPos += length
                else:
                    print('Unsupported CIGAR')
                    break
            samLine = samFile.readline().split('\t')
            samSeqName = samLine[0]
            if samSeqName != cigSeqName:
                break
        idx = 0
        isCheckLeft = True
        isMatchLeft = False
        isMatchRight = False
        for match in re.findall(r"\d+[A-Z]", cigCigar):
            length = int(re.search(r"\d+", match).group())
            code = match[-1]
            strand = cigLine[-2]
            if code == 'M':
                if isCheckLeft:
                    startRight = length
                    isCheckLeft = False
                    for i in range(length):
                        if predAlignInfo[i] is not None \
                           and trueAlignInfo[i][0] == predAlignInfo[i][0]\
                           and trueAlignInfo[i][1] == predAlignInfo[i][1]:
                               isMatchLeft = True
                else:
                    startRight -= 1
                    for i in range(length):
                        if predAlignInfo[startRight+i] is not None \
                           and trueAlignInfo[startRight+i][0] == predAlignInfo[startRight+i][0]\
                           and trueAlignInfo[startRight+i][1] == predAlignInfo[startRight+i][1]:
                               isMatchRight = True
        print(cigSeqName, isMatchLeft, isMatchRight)
        total += 1
        if isMatchLeft: left += 1
        if isMatchLeft: right += 1
        if isMatchLeft and isMatchRight: both += 1
    else:
        while True:
            samLine = samFile.readline().split('\t')
            samSeqName = samLine[0]
            if samSeqName != cigSeqName:
                break
print("total:{}, both correct: {}, left correct: {}, right correct: {}".format(total, both, left, right))
cigFile.close()
