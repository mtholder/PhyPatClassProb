#!/bin/sh
inp="$1"
lf="$2"
nt="$3"
shape="$4"
python testnex.py $nt $shape > "$inp"
paup -n "$inp"
python SumUninf.py paupsitelikes "$lf" >paupSummary.txt
../src/UninformProb -f.25,0.25,0.25 -r1,1,1,1,1 "$inp" 2>/dev/null >ourOutput.txt || exit
python our_prob.py ourOutput.txt paupSummary.txt || exit
echo "Passed"