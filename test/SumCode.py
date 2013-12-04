for column in logFileO: 
        if column.strip()=='': 
            break
        c = column.split()
        charStatus = x[5]
        ind = int(x[0])
        if charStatus in charStatusList:
            Uninformative.append(ind)
            return Uninformative
def sumProbs(probFileO, included):
    for column in likeFileO: 
        if column.startswith('ACGT'): 
            break
        if column.startswith('ACT'):
            break
        if column.startswith('AGT'):
            break
        if column.startswith('AT'):
            break
        if column.startswith('CGT'):
            break
        if column.startswith('CT'):
            break
        if column.startswith('GT'):
            break
    sum = 0.0
    logSum = 0.0
    for line in likeFileO: 
        if line[0].isdigit(): 
            break
        x = line.split()
        ind = int(x[0])
        lnL = -float(x[1])
        print like
        if (included is None) or (ind in included):
            sum += like
        logSum += lnL
    print sum
    print 1-sum
    print logSum
    return sum