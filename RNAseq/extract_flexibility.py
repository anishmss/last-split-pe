from __future__ import print_function

import sys
import re

# grand truth
trueJunctions = {}
with open(sys.argv[2]) as f:
    for line in f:
        readName, refName, start, end = re.split('[\t:-]', line.rstrip())
        info = (refName, int(start), int(end))
        if readName in trueJunctions:
            trueJunctions[readName].append(info)
        else:
            trueJunctions[readName] = [info]
# ambiguity
ambiguity = {}
with open(sys.argv[3]) as f:
    for line in f:
        name, start, end, prefix, suffix = line.rstrip().split('\t')
        ambiguity[(name, int(start), int(end))] = (int(prefix), int(suffix))

correct = 0
incorrect = 0
crossed = 0
notcrossed = 0

with open(sys.argv[1]) as f:
    for line in f:
        readName, refName, start, end = re.split('[\t:-]', line.rstrip())
        start = int(start)
        end = int(end)
        isCorrect = False
        isCrossed = False
        r = -1
        l = -1
        if readName in trueJunctions:
            for info in trueJunctions[readName]:
                if refName != info[0]:
                    continue
                else:
                    # is edge correct?
                    trueIntronStart = info[1] + 1
                    trueIntronEnd   = info[2] - 1
                    intronInfo = (info[0], trueIntronStart, trueIntronEnd)
                    if intronInfo in ambiguity:
                        prefixAmbiguity, suffixAmbiguity = ambiguity[intronInfo]
                        PIS = start + 1 # predicted intron start-point
                        PIE = end - 1 # predicted intron end-point
                        if (PIS >= trueIntronStart - suffixAmbiguity - 1 and \
                            PIS <= trueIntronStart + prefixAmbiguity + 1) and \
                           (PIE >= trueIntronEnd - suffixAmbiguity -1 and \
                            PIE <= trueIntronEnd + prefixAmbiguity +1 ):
                               isCorrect = True
                               break
                    # is crossed?
                    for flexibility in [0,1,2,3,4,5]:
                        if (PIS >= trueIntronStart - suffixAmbiguity - flexibility and \
                            PIS <= trueIntronStart + prefixAmbiguity + flexibility) and \
                           (PIE >= trueIntronEnd - suffixAmbiguity - flexibility and \
                            PIE <= trueIntronEnd + prefixAmbiguity + flexibility):
                            isCrossed = True
        if isCorrect:
            correct += 1
            #print(line.rstrip())
        else:
            incorrect += 1
            #print(line.rstrip(), file=sys.stderr)
        if isCrossed:
            crossed += 1
        else:
            notcrossed += 1
        if isCrossed and not isCorrect:
            print(line.rstrip())
            #print("{},{},{}".format(readName, r, l))
print("precision: {}% correct|incorrect|total = {}|{}|{}".format(correct/float(correct+incorrect), correct, incorrect, (correct+incorrect)))
print("crossed precision: {}% crossed|not crossed|total = {}|{}|{}".format(crossed/float(crossed+notcrossed), crossed, notcrossed, (crossed+notcrossed)))
