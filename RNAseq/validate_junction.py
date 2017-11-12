import sys
import re

predFile = open(sys.argv[1])
trueFile = open(sys.argv[2])

none = 0
left = 0
right = 0
both = 0
for predLine in predFile:
    pReadName, pRefName, pStart, pEnd = re.split('[\t:-]', predLine.rstrip())
    pReadNum = int(re.search(r'(\d+)', pReadName).group(1))
    while True:
        trueLine = trueFile.readline().rstrip()
        tReadName, tRefName, tStart, tEnd = re.split('[\t:-]', trueLine)
        if tReadName == pReadName:
            if pStart == tStart and pEnd == tEnd:
                both += 1
            elif pStart == tStart and pEnd != tEnd:
                left += 1
            elif pStart != tStart and pEnd == tEnd:
                right += 1
            else:
                none += 1
            break
        else:
            tReadNum = int(re.search(r'(\d+)', tReadName).group(1))
            if tReadNum > pReadNum:
                none += 1
                break
print("Junctions Sides (none|left|right|both): {}|{}|{}|{}".format(none, left, right, both))
