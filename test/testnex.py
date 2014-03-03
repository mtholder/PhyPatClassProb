#!/usr/bin/env python
import sys
def write_column(char_list, stream_list):
    for i,c in enumerate(char_list):
        stream_list[i].write(c)

def add_next(col_list,levels,stream_list):
    #print "    "*(2 - levels), "add_next col_list = ", str(col_list), " levels = ", levels
    for letter in 'ACGT':
        c = list(col_list)
        c.append(letter)
        if levels>0:
            add_next(c,levels-1,stream_list)
        else:
            write_column(c,stream_list)
from cStringIO import StringIO
stream_list=[]

shape = 'p'
if len(sys.argv) > 2:
    shape = sys.argv[2].lower()
    if shape not in 'bp':
        sys.exit('Shape must be "b" or "p"\n')
if shape == 'b':
    ex = int(sys.argv[1])
    ns = 2**ex
else:
    ns = int(sys.argv[1])
for i in range(ns):
    stream_list.append(StringIO())

add_next([], ns - 1, stream_list)
length = len(stream_list[0].getvalue())

print '''#NEXUS
begin data;
dimensions ntax = '''+str(ns)+''' nchar = '''+str(length)+''';
    format datatype = dna;
    matrix'''

for i in range(ns):
    print 's'+str(i+1), stream_list[i].getvalue()
treestr = ''
def build_balanced(exponent, start, blen):
    if exponent==1:
        return ('(s{x}:{blen},s{i}:{blen}): {blen}'.format(x=start, blen=blen, i=start+1), start+2)
    else:
        left_str, startn = build_balanced(exponent-1, start, blen)
        right_str, startn = build_balanced(exponent-1, startn, blen)
        n='('+left_str+','+right_str+'): '+blen
        return n, startn
if shape == 'p':
    for i in range(ns): 
        if i<ns-1:
            treestr += '('
        treestr += 's' + str(i+1) + ':0.02'
        if i < ns-1:
            treestr += ','
    treestr += '):0.02'*(ns-2)
    treestr += ')'
else:
    shape = 'b'
    if len(sys.argv) > 2:
        shape = sys.argv[2].lower()
        if shape not in 'bp':
            sys.exit('Shape must be "b" or "p"\n')
    if shape == 'b':
        treestr = build_balanced(ex, 1, '0.02') 
        treestr = treestr[0]

print ''';
end;
begin paup;
    set storebrlens;
    log start replace;
    CStatus full ;
end;
begin trees;
    tree mine = [&U] ''' + treestr + ''';
end;
begin paup;
    lscore / user nst = 1 basefreq = eq pinv = 0 rates = equal sitelike ScoreFile = paupsitelikes replace;
end;    

'''
