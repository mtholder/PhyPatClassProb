#!/bin/sh
invoc='../../src/PhyPatClassProb -f0.252818,0.284001,0.256973 -r1.86149,5.25704,0.72040,0.97137,6.83511  -m0,0.129385382,0.348650344,0.57517613,0.831849857,1.14122287,1.54443828,2.148571055,3.696423412 -p0.23193,0.09600875,0.09600875,0.09600875,0.09600875,0.09600875,0.09600875,0.09600875'
end_point=12
for ((i=3; i < ${end_point}; ++i))
do
    python  gen_data_set.py p ${i} > p${i}.nex
    echo pectinate${i}
    time $invoc ./p${i}.nex 2>p${i}.err.txt >>p${i}.out.txt
    if [ $? != 0 ]
    then
        cat p${i}.err.txt
        exit 1
    fi

    python  gen_data_set.py s ${i} > s${i}.nex
    echo symm${i}
    time $invoc ./s${i}.nex 2>s${i}.err.txt >>s${i}.out.txt
    if [ $? != 0 ]
    then
        cat s${i}.err.txt
        exit 1
    fi
done
