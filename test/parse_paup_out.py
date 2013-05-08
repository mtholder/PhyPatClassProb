#!/usr/bin/env python
import sys, re
from math import exp

site_like_pattern = re.compile(r'\t\t\d+\s+([0-9.]+)\s+')
fn = sys.argv[1]
inp = open(fn, 'rU')
like_list = []
for line in inp:
    m = site_like_pattern.match(line)
    if m:
        like_list.append(exp(-float(m.group(1))))

sum_obs_by_rep_dict = {}
states = 'ACGT'
for t1 in states:
    it1 =  states.index(t1)
    bt1 = (1 << it1)
    for t2 in states:
        it2 =  states.index(t2)
        bt2 = (1 << it2)
        for t3 in states:
            it3 =  states.index(t3)
            bt3 = (1 << it3)
            for t4 in states:
                it4 =  states.index(t4)
                bt4 = (1 << it4)
                o_list = [bt1, bt2, bt3, bt4]
                obs = sum(set([i for i in o_list]))
                rep = set()
                i_list = [it1, it2, it3, it4]
                for i in i_list:
                    if i_list.count(i) > 1:
                        rep.add(i)
                rep_list = list(rep)
                rep_list.sort()
                rep_t = tuple(rep_list)
                like = like_list.pop(0)
                if obs in sum_obs_by_rep_dict:
                    o_d = sum_obs_by_rep_dict[obs]
                else:
                    o_d = {}
                    sum_obs_by_rep_dict[obs] = o_d

                if rep_t in o_d:
                    o_d[rep_t] += like
                else:
                    o_d[rep_t] = like
                
                print '%s%s%s%s %f %d %s %s' % (t1, t2, t3, t4, like, obs, bin(obs), str(rep_t))

for o in range(1, 16):
    o_d = sum_obs_by_rep_dict[o]
    k = o_d.keys()
    k.sort()
    for nk in k:
        print 'obs=', o, '  rep = ', nk, ' sum of likelihoods =', o_d[nk]
print sum(like_list)