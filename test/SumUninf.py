import math
import sys
def readLogFile(logFileO, charStatusList, charStateList):
    for line in logFileO: 
        #print line,
        if line.startswith('Character        Type         Status      Weight   States'): 
            break
    Uninformative = []
    logFileO.next()
    for line in logFileO: 
        if line.strip()=='': 
            break
        x = line.split()
        charStatus = x[2]
        ind = int(x[0])
        charState = x[4]
        if (charStatus in charStatusList) and (charState in charStateList):
            Uninformative.append(ind)
    return Uninformative
def sumLikes(likeFileO, included):
    for line in likeFileO: 
        #print line,
        if line.startswith('Tree	-lnL	Site	-lnL'): 
            break
    sum = 0.0
    logSum = 0.0
    for line in likeFileO: 
        if line[0].isdigit(): 
            break
        x = line.split()
        ind = int(x[0])
        lnL = -float(x[1])
        like = math.exp(lnL)
        #print like
        if (included is None) or (ind in included):
            sum += like
        logSum += lnL
    #print sum
    #print 1-sum
    #print logSum
    return sum

likeFile = sys.argv[1]
logFile = sys.argv[2]
likeFileO = open(likeFile, 'rU')
logFileO = open(logFile, 'rU')

nonconst = ['AC', 'AG', 'AT', 'ACG', 'ACT', 'AGT', 'ACGT', 'CGT', 'CT', 'CG', 'GT' ]
for states in nonconst:
    logFileO.seek(0)
    likeFileO.seek(0)
    u = readLogFile(logFileO, ['U'], [states]) #loop over all possible outcomes
    sl = sumLikes(likeFileO, u)
    print 'U {t:4s} {p:.15f}'.format(t=states, p=sl)

for states in 'ACGT':
    logFileO.seek(0)
    likeFileO.seek(0)
    u = readLogFile(logFileO, ['UC'], [states])
    sl = sumLikes(likeFileO, u)
    print 'c {t:4s} {p:.15f}'.format(t=states, p=sl) # c = constant

for states in nonconst:
    logFileO.seek(0)
    likeFileO.seek(0)
    u = readLogFile(logFileO, ['-'], [states])
    sl = sumLikes(likeFileO, u)
    print 'i {t:4s} {p:.15f}'.format(t=states, p=sl) #i = informative
