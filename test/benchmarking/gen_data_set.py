#!/usr/bin/env python
import sys
shape_arg = sys.argv[1].lower()
if shape_arg not in 'sp':
    sys.exit('Expecting "S" or "P" for symmetric or pectinate. Found "%s"\n' % sys.argv[1])
power = int(sys.argv[2])
assert power > 1
nt = 2**power
tree_string = 'TREE_HERE'
br_len = '0.05'

def gen_symm(tax_ind, depth):
	if depth > 1:
		left, tax_ind = gen_symm(tax_ind, depth - 1)
		right, tax_ind = gen_symm(tax_ind, depth - 1)
		ts = '(%s:%s,%s:%s)' % (left, br_len, right, br_len)
		return ts, tax_ind
	ts = '(t%d:%s,t%d:%s)' % (tax_ind, br_len, tax_ind + 1, br_len)
	return ts, tax_ind + 2
	
if shape_arg == 's':
	tree_string, ind = gen_symm(1, power)
elif shape_arg == 'p':
	tree_string = '(t1:%s,t2:%s)' % (br_len, br_len)
	for i in range(3, nt + 1):
		tree_string = '(%s:%s,t%d:%s)' % (tree_string, br_len, i, br_len)	

matrix = '\n'.join(['t%d A' % i for i in range(1, nt + 1)])
sys.stdout.write('''#NEXUS
BEGIN DATA;
    Dimensions NTax = %d NChar = 1;
    Format datatype = DNA ;
Matrix
%s
;
END;
BEGIN TREES; 
    Tree mod = [&U] %s ;
END;
''' % (nt, matrix, tree_string))
