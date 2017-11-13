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

with open(sys.argv[1]) as f:
    for line in f:
        readName, refName, start, end = re.split('[\t:-]', line.rstrip())
        start = int(start)
        end = int(end)
        isCorrect = False
        if readName in trueJunctions:
            for info in trueJunctions[readName]:
                if refName != info[0]:
                    continue
                else:
                    trueIntronStart = info[1] + 1
                    trueIntronEnd   = info[2] - 1
                    intronInfo = (info[0], trueIntronStart, trueIntronEnd)
                    if  intronInfo in ambiguity:
                        prefixAmbiguity, suffixAmbiguity = ambiguity[intronInfo]
                        PIS = start + 1 # predicted intron start-point
                        PIE = end - 1 # predicted intron end-point
                        if (PIS >= trueIntronStart - suffixAmbiguity and \
                            PIS <= trueIntronStart + prefixAmbiguity) and \
                           (PIE >= trueIntronEnd - suffixAmbiguity and \
                            PIE <= trueIntronEnd + prefixAmbiguity):
                               isCorrect = True
                               break
        if isCorrect:
            correct += 1
        else:
            incorrect += 1
        print(line, isCorrect)
print("correct|incorrect = {}|{}".format(correct, incorrect))
