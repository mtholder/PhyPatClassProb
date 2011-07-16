#! /usr/bin/env python
import sys, os, math

################################################################################
# Parse and validate command line
################################################################################
from optparse import OptionParser
import re
usage='''Creates a data set for testing PhyPatClassProb
Example:

python generate_test.py  -f.3,.4,.1 -r1.1,2,1.5,.8,3 -e0.2


'''
parser = OptionParser(usage=usage)
parser.add_option("-f", "--base-freq", dest="base_freq", default="0.25,0.25,0.25",
                  help="comma-separated list of the base frequencies: freqA,freqC,freqG (the frequency of T will be computed from the other three)")
parser.add_option("-r", "--rel-rates", dest="rel_rates", default="1.0,1.0,1.0,1.0,1.0",
                  help="comma-separated list of the exchangeability rates relative to G<->T: rAC,rAG,rAT,rCG,rCT")
parser.add_option("-t", "--tree", dest="tree", default=None,
                  help="Newick tree (with consecutively numbered leaves starting at 1)")
parser.add_option("-p", "--paup-file", dest="paup_file", default=None, type='str',
                  help="Optional file name for a NEXUS file that will check the probability calculations")

(options, args) = parser.parse_args()
if len(args) > 0:
    sys.exit('No arguments are expected. Use the -h option to see all of the options')

# Validate state frequencies
state_freq_str_list = options.base_freq.strip().split(',')
if len(state_freq_str_list) != 3:
    sys.exit('Expecting three base frequencies separated by commas (the frequency of T will be calculated from the others)')
try:
    state_freq = [float(i) for i in state_freq_str_list]
except:
    sys.exit('Error interpreting the base frequencies as numbers.')
if min(state_freq) <= 0.0:
    sys.exit('All base frequencies must be positive.')
s = sum(state_freq)
if s >= 1.0:
    sys.exit('The sum freqA + freqC +freqG must be less than 1.0')
state_freq.append(1.0-s)
freqA, freqC, freqG, freqT = state_freq

# Validate relrate matrix
rel_rates_str_list = options.rel_rates.strip().split(',')
if len(rel_rates_str_list) != 5:
    sys.exit('Expecting five relative rates to be separated by commas (the exchangeability of G<->T is set to 1, and the other rates are expressed relative to that rate)')
try:
    rel_rates = [float(i) for i in rel_rates_str_list]
except:
    sys.exit('Error interpreting the relative rates as numbers.')
if min(rel_rates) <= 0.0:
    sys.exit('All relative rates must be positive.')
rAC, rAG, rAT, rCG, rCT = rel_rates
rGT = 1.0

if not (options.tree):
    sys.exit('A tree must be specified')
# Validate the edge length
tax_ind_pat = re.compile(r'[,(](\d+):')
try:
    inds = [int(i) for i in tax_ind_pat.findall(options.tree)]
    assert(len(inds) > 2)
except:
    sys.exit('Error parsing taxon numbers out of the newick string (which should be consecutively numbered leaves starting at 1)')
num_tax = max(inds)
if num_tax != len(inds):
    sys.exit('max taxon number != the number of taxon numbers found (newick string should be consecutively numbered leaves starting at 1)')
if min(inds) != 1:
    sys.exit('Newick string should be consecutively numbered leaves starting at 1.')
if len(set(inds)) != num_tax:
    sys.exit('Newick string should be consecutively numbered leaves starting at 1 (repeated index found)')




sys.exit(0)



################################################################################
# call beagle in phython
################################################################################
from pytbeaglehon.disc_state_cont_time_model import RevDiscStateContTimeModel

m = RevDiscStateContTimeModel(state_freq=state_freq,
                              r_upper=[[rAC, rAG, rAT ],
                                       [     rCG, rCT ],
                                       [          rGT ]])
p_mat = m.prob_matrices(edge_length)
################################################################################

################################################################################
# Format output
################################################################################
p_mat_str = [[str(i) for i in row] for row in p_mat[0]]
max_str_len = max([max([len(i) for i in row]) for row in p_mat_str])
p_mat_str = [[i.ljust(max_str_len) for i in row] for row in p_mat_str]
for row in p_mat_str:
    print '\t'.join(row)

################################################################################
# Produce optional output
################################################################################
if options.paup_file:
    if os.path.exists(options.paup_file):
        verb = 'Overwriting'
    else:
        verb = 'Creating'
    sys.stderr.write(verb + ' ' + options.paup_file + '\n')
    
    expected_site_like_list = []
    site_n = 1
    for anc in range(4):
        state_f = state_freq[anc]
        for des in range(4):
            ti_prob = p_mat[0][anc][des]
            like = state_f*ti_prob
            lnL = math.log(like)
            s = 'Site %2d L = %f * %f = %f      lnL = %f' % (site_n, state_f, ti_prob, like, lnL)
            expected_site_like_list.append(s)
            site_n += 1
    expected_site_like = '\n'.join(expected_site_like_list)
            
            
    p = open(options.paup_file, 'w')
    p.write('''#NEXUS
Begin DATA;
    Dimensions NTax = 4 NChar = 16;
    Format Datatype = dna;
    Matrix 
        one     AAAA CCCC GGGG TTTT 
        two     AAAA CCCC GGGG TTTT 
        three   ACGT ACGT ACGT ACGT 
        four    ACGT ACGT ACGT ACGT 
    ;
End;
Begin PAUP;
    set storebrlens;
End;
Begin TREES;
    tree test = [&U] ( one : 0.0 , two : 0.0 , ( three : 0.0 , four : 0.0 ) : %(edge_length)f );
End;
Begin PAUP;
    Lscore / UserBrLens 
             NST=6
             BaseFreq=(%(freqA)f %(freqC)f %(freqG)f)
             rMatrix=(%(rAC)f %(rAG)f %(rAT)f %(rCG)f %(rCT)f)
             PInv=0
             Rates=Equal
             SiteLike;
End;    
[! 

We are expecting the following site likelihoods (produced via multiplying the 
    transition probability by the state frequency for the ancestral state:
%(expected_site_like)s
]
''' % {'edge_length' : edge_length, 
       'freqA' : freqA, 'freqC' : freqC, 'freqG' : freqG,
       'rAC': rAC, 'rAG': rAG, 'rAT': rAT, 'rCG': rCG, 'rCT': rCT, 
       'expected_site_like' : expected_site_like
       })
