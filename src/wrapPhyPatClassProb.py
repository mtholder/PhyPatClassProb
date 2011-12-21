#! /usr/bin/env python
import sys, os, math, re, copy, itertools

################################################################################
# Parse and validate command line
################################################################################
from optparse import OptionParser

run_invocation = False
usage='''Creates a data set for testing PhyPatClassProb
Example:

python wrapPhyPatClassProb.py  -f.3,.4,.1 -r1.1,2,1.5,.8,3 nexusfilename

to analyze the NEXUS file nexusfilename (which must contain a data matrix and 
    a tree).
'''
parser = OptionParser(usage=usage)
parser.add_option("-f", "--base-freq", dest="base_freq", default="0.25,0.25,0.25",
                  help="comma-separated list of the base frequencies: freqA,freqC,freqG (the frequency of T will be computed from the other three)")
parser.add_option("-r", "--rel-rates", dest="rel_rates", default="1.0,1.0,1.0,1.0,1.0",
                  help="comma-separated list of the exchangeability rates relative to G<->T: rAC,rAG,rAT,rCG,rCT")
parser.add_option("-i", "--prop-invar", dest="invar", default=None, type='float',
                  help="Proportion of invariant sites")
parser.add_option("-a", "--alpha", dest="alpha", default=None, type='float',
                  help="Shape parameter of the gamma distribution for gamma-distributed rates across sites (no rate gamma rate het is used if this is not specified)")
parser.add_option("-c", "--num-var-cat", dest="n_var_cat", default=4, type='int',
                  help="The number of discrete categories used to approximate the gamma-distributed rates.")
parser.add_option("-p", "--paup-file", dest="paup_file", default=None, type='str',
                  help="Optional file name for a NEXUS file that will check the probability calculations")
parser.add_option("-o", "--output-prefix", dest="ofprefix", default="test", type="str",
                  help="Prefix for output files for a prefix XYZ, the files XYZoutput.txt, XYZerror.txt, and XYZsummary.txt may be created/overwritten")
(options, args) = parser.parse_args()
if len(args) != 1:
    sys.exit('One arguments are expected. Use the -h option to see all of the options')
fn = args[0]


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

rates = []
rate_probs = []
var_rate_prob = 1.0
mean_var_rate = 1.0
if options.invar is not None:
    if options.invar < 0.0 or options.invar >= 1.0:
        sys.exit('proportion of invariant sites must be non-negative and less than 1.')
    var_rate_prob -= options.invar
    mean_var_rate = 1.0/var_rate_prob
    rates.append(0.0)
    rate_probs.append(options.invar)
if options.alpha is not None:
    if options.alpha <= 0.0:
        sys.exit('alpha most be a positive value')
    if options.n_var_cat < 2:
        sys.exit('num-var-cat must be > 1')
    try:
        from pytbeaglehon.asrv import GammaRateHetManager
    except:
        sys.exit('pytbeaglehon must be import-able (installed or on your PYTHONPATH variable)')
    grh = GammaRateHetManager(shape=options.alpha, num_categories=options.n_var_cat)
    raw_rates = grh.rates
    final_rates = [mean_var_rate*i for i in raw_rates]
    rates.extend(final_rates)
    rate_cat_probs = [var_rate_prob/options.n_var_cat] * options.n_var_cat
    rate_probs.extend(rate_cat_probs)
else:
    rate_probs.append(var_rate_prob)
    rates.append(mean_var_rate)
    
prefix = options.ofprefix
stdo, stde = prefix + 'output.txt', prefix + 'error.txt'
d = os.path.split(os.path.dirname(os.path.abspath(sys.argv[0])))[0]
opts = '-f%s -r%s' % (options.base_freq, options.rel_rates)
if len(rates) > 1:
    assert(len(rates) == len(rate_probs))
    opts = opts + ' -m%s -p%s' %(','.join([str(i) for i in rates]), ','.join([str(i) for i in rate_probs[:-1]]))
    

if run_invocation:
    invoc = '%s %s %s >%s 2>%s' % (os.path.join(d, 'src', 'PhyPatClassProb'), opts, fn, stdo, stde)
    sys.stderr.write('Creating/Overwriting %s and %s by running:\n   %s\n' % (stdo, stde, invoc))
    rc = os.system(invoc)
    if rc != 0:
        sys.stderr.write('Error! return code = %d\n' % rc)

summary_fn = prefix + 'summary.txt'
sys.stderr.write('Creating/Overwriting %s\n' % (summary_fn))
summary_f = open(summary_fn, 'w')
proc_out = open(stdo, 'rU')
import re
expPat = re.compile(r'Expected steps = (\d+) states = (\S+) prob = ([-e.0-9]+)')
pat_prob_table = []
for line in proc_out:
    m = expPat.match(line)
    if m:
        g = m.groups()
        ind, state, prob = int(g[0]), g[1], float(g[2])
        if ind >= len(pat_prob_table):
            pat_prob_table.append([])
        pat_prob_table[ind].append([state, prob])
    else:
        break
obsPat = re.compile(r'Observed steps = (\d+) states = (\S+) count = (\d+)')
pat_count_table = []
num_c = 0
for line in proc_out:
    m = obsPat.match(line)
    if m:
        g = m.groups()
        ind, state, count = int(g[0]), g[1], int(g[2])
        num_c += count
        if ind >= len(pat_count_table):
            pat_count_table.append([])
        second_ind = len(pat_count_table[ind])
        assert(pat_prob_table[ind][second_ind][0] == state)
        pat_count_table[ind].append([state, count])

expected_count_table = copy.deepcopy(pat_prob_table)
for row in expected_count_table:
    for ind in range(len(row)):
        p = row[ind][1]
        row[ind][1] = num_c*p

assert(len(expected_count_table) == len(pat_count_table))
summary_f.write('%10s\t%5s\t%10s\t%10s\t%10s\t%10s\n' % ('Num_steps', 'States', 'Expected', 'Observed', 'Difference', 'ChiSqStat'))
chi_sq = 0.0
for num_steps, rows in enumerate(itertools.izip(expected_count_table, pat_count_table)):
    e_row, o_row = rows
    assert(len(e_row) == len(o_row))
    for e_el, o_el in itertools.izip(e_row, o_row):
        state_set = e_el[0]
        assert(state_set == o_el[0]) # states should be in the same order, we store the state as element 0 so that we can double check
        e_c, o_c = e_el[1], o_el[1]
        d = o_c - e_c
        if e_c > 0.0:
            chi_sq_stat = d*d/e_c
            chi_sq += chi_sq_stat
            summary_f.write('%10d\t%5s\t%10.2f\t%10d\t%10.2f\t%10.2f\n' % (num_steps, state_set, e_c, o_c, d, chi_sq_stat))

summary_f.write('\n\nTotal Chi-Squared Stat = %10.2f\n' % chi_sq)



























