import math
import sys
likeFile = sys.argv[1]
logFile = sys.argv[2]
likeFileO = open(likeFile, 'UC')
logFileO = open(logFile, 'U')
included = None
def sumLikes(likeFileO, included):
    for line in likeFileO: 
        print line,
        if line.startswith('Tree	-lnL	Site	-lnL'): 
            break
    sum = 0.0
    logSum = 0.0
    for line in likeFileO: 
        if line[0].isdigit(): 
            break
        x = line.split()
        ind int(x[0])
        lnL = -float(x[1])
        like = math.exp(lnL)
        print like
        if (included is None) or (ind in included):
            sum += like
        logSum += lnL
    print sum
    print 1-sum
    print logSum
    return sum
