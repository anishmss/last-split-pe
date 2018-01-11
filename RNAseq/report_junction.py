import sys
import re

samFile = open(sys.argv[1])
minIntronSize = int(sys.argv[2])
#threshold = -1
#if len(sys.argv) > 2:
#    threshold = int(sys.argv[2])
#>>>>>>> fd2fd3532ad21a59c8e1fa76a42900b126ab8a01

while True:
    samLine = samFile.readline()
    if samLine[0] != '@':
        break
samLine = samLine.rstrip().split('\t')
junctionDict = {} #???
while True:
    matchDict = {}
    while True: #loop for lines of same read
        #print samLine
        samSeqName = samLine[0] 
        samCigar = samLine[5]
        refName = samLine[2]
        refPos = int(samLine[3])
        flag = int(samLine[1])
        isForward = ((flag & 0x10) == 0)
        consumeDir = 1 if isForward else -1
        refPos = int(samLine[3])
        readLength = 0
        
        # why is this needed? just read in readlength?
        for match in re.findall(r"\d+[A-Z]", samCigar):
            #print match
            #print re.search(r"\d+", match).group()
            readLength += int(re.search(r"\d+", match).group()) # this is wrong if cigar contains D.
        #print "readLength:",readLength
        readPos = 0 if isForward else readLength-1
        
        for match in re.findall(r"\d+[A-Z]", samCigar):
            length = int(re.search(r"\d+", match).group())
            code = match[-1]
            if code == 'M' or code == '=' or code == 'X':
                for i in range(length):
                    if refName not in matchDict:
                        matchDict[refName] = {}
                    matchDict[refName][refPos] = 'M'
                    refPos += 1
                    readPos += consumeDir
            elif code == 'D':
                for i in range(length):
                    if refName not in matchDict:
                        matchDict[refName] = {}
                    matchDict[refName][refPos] = 'D'
                    refPos += consumeDir
            elif code == 'H' or code == 'I' or code == 'S':
                readPos += length * consumeDir
            elif code == 'N':
                refPos += length * consumeDir
        
    #while True:
        #samMapQ = int(samLine[4])
        #samSeqName = samLine[0]
        #if samMapQ > threshold:
            #samCigar = samLine[5]
            #refName = samLine[2]
            #refPos = int(samLine[3])
            #flag = int(samLine[1])
            #isForward = ((flag & 0x10) == 0)
            #consumeDir = 1 if isForward else -1
            #refPos = int(samLine[3])
            #readLength = 0
            #for match in re.findall(r"\d+[A-Z]", samCigar):
                #readLength += int(re.search(r"\d+", match).group())
            #readPos = 0 if isForward else readLength-1
            #for match in re.findall(r"\d+[A-Z=]", samCigar):
                #length = int(re.search(r"\d+", match).group())
                #code = match[-1]
                #if code == 'M' or code == '=' or code == 'X':
                    #for i in range(length):
                        #if refName not in matchDict:
                            #matchDict[refName] = {}
                        #matchDict[refName][refPos] = 'M'
                        #refPos += 1
                        #readPos += consumeDir
                #elif code == 'D':
                    #for i in range(length):
                        #if refName not in matchDict:
                            #matchDict[refName] = {}
                        #matchDict[refName][refPos] = 'D'
                        #refPos += 1
                #elif code == 'H' or code == 'I' or code == 'S':
                    #readPos += length * consumeDir
                #elif code == 'N':
                    #refPos += length
#>>>>>>> fd2fd3532ad21a59c8e1fa76a42900b126ab8a01
        samLine = samFile.readline().rstrip().split('\t')
        if samSeqName != samLine[0]:
            break
    for refName in matchDict.keys():
        matchPos = sorted(matchDict[refName].keys())
        prePos = matchPos[0]
        for pos in matchPos[1:]:
            if pos-prePos > 1 and matchDict[refName][prePos] == 'M' and matchDict[refName][pos]=='M':
                if pos-prePos-1 >= minIntronSize:
                    print("{}\t{}:{}-{}".format(samSeqName, refName, prePos, pos))
            prePos = pos
    if len(samLine) <= 1:
        break
