x = '''8    0m0.018 0m0.019
16  0m0.025 0m0.040
32  0m0.058 0m0.182
64  0m0.197 0m1.068
128 0m0.768 0m4.847
256 0m3.167 0m20.289
512 0m16.296        1m52.845
1024        1m14.894        12m11.998
2048        4m55.558        54m51.461
4096        18m59.904       212m51.392'''
for line in x.split('\n'):
    word_list = line.split()
    pw = word_list[1].split('m')
    sw = word_list[2].split('m')
    word_list[1] = str(int(pw[0])*60 + float(pw[1]))
    word_list[2] = str(int(sw[0])*60 + float(sw[1]))
    print '\t'.join(word_list)
