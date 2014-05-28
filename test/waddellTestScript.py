print '''
Begin characters;

	Dimensions ntax = 40 newtaxa '''

import math
import sys
repeated = sys.argv[1][0]
singleton = sys.argv[1][1:]
ns = len(sys.argv[1])
nt = 40

def calcNumPat(n, k): 
    return math.factorial(n)/(math.factorial(n-k))
nc = calcNumPat(nt, ns-1)
sys.stdout.write('nchar= {n};'.format(n=nc))

print '''	Format transpose datatype=nucleotide gap=-;

TaxLabels
	Ovis_canadenensis
	Okapia_johnstoni
	Lagenorhynchus_obscurus
	Balaenoptera_physalus
	Sus_scrofa
	Camelus_dromedarius
	Equus_grevyi
	Diceros_bicornis
	Tapirus_indicus
	Phoca_vitulina
	Ailurus_fulgens
	Panthera_uncia
	Manis_sp
	Manis_pentadactyla
	Myotis_vellifer
	Rhinolophus_creaghi
	Nyctimene_albiventer
	Erinaceus_europaeus
	Scalopus_aquaticus
	Talpa_sp
	Uropsilus_sp
	Cyncephalus_volans
	Homo_sapiens
	Tarsius_syrichta
	Mus_musculus
	Dolichotis_patagonum
	Hystrix_africaeaustralis
	Sylvilagus_sp
	Ochotona_sp
	Amblysomus_hottentotus
	Tenrec_sp
	Macroscelides_sp
	Orycteropus_afer
	Procavia_capensis
	Dugong_dugon
	Trichechus_manatus
	Elephas_maximus
	Cyclopes_didactylus
	Choloepus_hoffmani
	Dasypus_novemcinctus
	;
matrix
	'''
count = 1

if ns == 1:
    p = [repeated]*nt
    pattern = ''.join(p)
    print 'char_{j} {p}'.format(j=(1), p=pattern)
if ns == 2:
    for i in range(40):
        p = [repeated]*nt
        p[i] = singleton
        pattern = ''.join(p)
        print 'char_{j} {p}'.format(j=(1+i), p=pattern)
if ns == 3:
    for i in range(40):
        for j in range(40):
            if i != j: 
                p = [repeated]*nt
                p[i] = singleton[0]
                p[j] = singleton[1]
                pattern = ''.join(p)
                print 'char_{j:<4} {p}'.format(j=(count), p=pattern)
                count += 1            
if ns == 4:
    for i in range(40):
        for j in range(40): 
            if i == j:
                continue
            for k in range(40):
                if k == i or k == j:
                    continue
                p = [repeated]*nt
                p[i] = singleton[0]
                p[j] = singleton[1]
                p[k] = singleton[2]
                pattern = ''.join(p)
                print 'char_{j:<5} {p}'.format(j=(count), p=pattern)
                count += 1
print ''';
End;

begin paup;
    set warnreset=no autoclose=yes storebr = yes ;
end;
Begin trees;  [Tree copied by MTH from the fromAppendix1MLTree.tre referenced above]
    tree PAUP_1 = [&U] ((Ovis_canadenensis:0.032394,Okapia_johnstoni:0.017635):0.01,((Lagenorhynchus_obscurus:0.027034,Balaenoptera_physalus:0.006414):0.032389,(Sus_scrofa:0.066640,(Camelus_dromedarius:0.052440,(((((Equus_grevyi:0.034393,(Diceros_bicornis:0.029295,Tapirus_indicus:0.031382):0.010001):0.007005,((Myotis_vellifer:0.040913,Nyctimene_albiventer:0.035773):0.010546,Rhinolophus_creaghi:0.032744):0.010211):0.001519,(Manis_sp:0.042264,Manis_pentadactyla:0.014133):0.048143):0.001972,((Erinaceus_europaeus:0.100990,((Scalopus_aquaticus:0.032319,Talpa_sp:0.029934):0.014280,Uropsilus_sp:0.042655):0.038770):0.017019,(((Cyncephalus_volans:0.068755,(Homo_sapiens:0.029957,Tarsius_syrichta:0.118360):0.014135):0,((Mus_musculus:0.142037,(Dolichotis_patagonum:0.052944,Hystrix_africaeaustralis:0.039950):0.048807):0.010413,(Sylvilagus_sp:0.048551,Ochotona_sp:0.073138):0.034040):0):0.008779,(((((Amblysomus_hottentotus:0.060760,Tenrec_sp:0.165549):0.001972,Macroscelides_sp:0.131638):0.002285,Orycteropus_afer:0.058117):0.003145,((Procavia_capensis:0.058589,(Dugong_dugon:0.008694,Trichechus_manatus:0.002496):0.024717):0.003100,Elephas_maximus:0.027923):0.004373):0.021244,((Cyclopes_didactylus:0.036600,Choloepus_hoffmani:0.035142):0.009473,Dasypus_novemcinctus:0.048489):0.013952):0.006732):0.004735):0):0.002019,((Phoca_vitulina:0.025428,Ailurus_fulgens:0.023456):0.013763,Panthera_uncia:0.041133):0.021334):0.012236):0.002302):0.008760):0.044426);
End;
begin paup;
    lscore / sitelike UserBrLens nst = 1 basefreq = equal rates = equal pinv = 0 ScoreFile = likeScores.txt;
end;

'''

