import sys
likeFile = sys.argv[1]
logFile = sys.argv[2]
likeFileO = open(likeFile, 'rU')
logFileO = open(logFile, 'rU')
def readLogFile(logFileO):
    for line in logFileO: 
        print line,
        if line.startswith('Character        Type         Status      Weight   States'): 
            break
    Uninformative = []
    logFileO.next()
    for line in logFileO: 
        if line.strip()=='': 
            break
        x = line.split()
        print x[2]
        ind = int(x[0])
    return Uninformative