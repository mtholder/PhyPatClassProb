def b():
    assert False
    
def c():
    b()

pk = paup_probs.keys()
ok = our_probs.keys()
pk.sort()
ok.sort()
assert pk == ok

for k in pk:
    pp = paup_prob[k]
    op = our_prob[k]
    print k, pp, op
    assert abs(pp-op)<0.00001