import sys
our_prob = dict()
ourLikeFile = sys.argv[1]
paupLikeFile = sys.argv[2]

for line in open(ourLikeFile, "rU"):
    if line.startswith('Prob('):
        ls = line.split()
        states = ls[5]
        prob = float(ls[-1])
        
        if states in our_prob:
            our_prob[states] = our_prob[states]+prob
        else:
            our_prob[states] = prob
#        our_prob[states] = prob + our_prob.get(states,0.0)

paup_prob = dict()
for line in open(paupLikeFile, 'rU'):
    if line.strip() and not line.startswith('i'): #line.strip skips whitespace
        ls = line.split()
        states = ls[1] #3 
        prob = float(ls[2])
        paup_prob[states] = prob

pk = paup_prob.keys()
ok = our_prob.keys()
pk.sort()
ok.sort()
print pk
print ok
assert pk == ok

for k in pk:
    pp = paup_prob[k]
    op = our_prob[k]
    print k, pp, op
    assert abs(pp-op)<0.00001