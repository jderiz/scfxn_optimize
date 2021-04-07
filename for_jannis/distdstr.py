#!/usr/bin/python

import sys,glob,os,copy
from math import log,sqrt

CO_DIS2CUT = 64.0
MINSEP = 10
ROSETTADB = '/home/muellebk/rosetta/Rosetta_vanilla/Rosetta/main/database'
ATOMFILE = '%s/chemical/atom_type_sets/fa_standard/atom_properties.txt'%ROSETTADB
DBPATH = '%s/chemical/residue_type_sets/fa_standard/residue_types/l-caa'%ROSETTADB
REFFILE = 'distr.REF'

MINVAL = 1.5
MAXVAL = 5.0
BINSIZE = 0.05
NBINS = int((MAXVAL-MINVAL)/BINSIZE)
P_PSEUDO = 0.01

AAS = ['MET','LEU','ILE','VAL','ALA', #0~4
       'TRP','PHE','TYR', #5~7
       'CYS','PRO','GLY', #8~10
       'SER','THR','ASN','GLN','HIS', #11~15
       'ASP','GLU','ARG','LYS', #16~19
       ]

atmtypes = ['CNH2','COO','CH1','CH2','CH3','aroC','Ntrp','Nhis',
            'NH2O','Nlys','Narg','Npro','OH','ONH2','OOC','S',
            'Nbb','CAbb','CObb','OCbb',
            'Hpol','Haro','HNbb'
            ] #,'Hapo'

#################### misc functions ################################
def distance(mol1crd,mol2crd):
    displ=[]
    for i in range(3):
        displ.append(abs(float(mol1crd[i])-float(mol2crd[i])))
    dist=displ[0]**2+displ[1]**2+displ[2]**2
    dist=dist**0.5
    return dist

def d2(mol1crd,mol2crd):
    displ=[]
    for i in range(3):
        displ.append(abs(float(mol1crd[i])-float(mol2crd[i])))
    dist2=displ[0]**2+displ[1]**2+displ[2]**2
    return dist2

def pdb2crd(pdbfile,opt,res_in=[],ignore_insertion=False,chaindef=[]):
    pdbcont=file(pdbfile)

    if res_in == []:
        res_defined = False
    else:
        res_defined = True

    bbatomlist=[' C  ',' CA ', ' N  ',' O  ',' H  ']
    bblist = [' CA ', ' C  ', ' N  ']
    mclist = [' CA ', ' C  ', ' N  ',' O  ']
    crd={}
    for line in pdbcont:
        if line[:4]!='ATOM':
            continue
        resno = int(line[22:26])
        restype = line[16:20].strip()
        chain = line[21]
        if ignore_insertion and line[26] != ' ':
            continue
        if res_defined and resno not in res_in:
            continue
        if chaindef != [] and chain not in chaindef:
            continue

        atmtype = line[12:16].strip()
        if opt == 'CA':
            if line[12:16] == ' CA ':
                if resno in crd:
                    continue
                crd[resno] = [float(line[30+i*8:38+i*8]) for i in range(3)]
        elif opt == 'CB':
            if (restype == 'GLY' and line[12:16] == ' CA ') \
                   or line[12:16] == ' CB ':
                if resno in crd:
                    continue
                crd[resno] = [float(line[30+i*8:38+i*8]) for i in range(3)]
        else:
            if resno not in crd:
                crd[resno] = {}
            crd[resno][atmtype] = [float(line[30+i*8:38+i*8]) for i in range(3)]
    pdbcont.close()
    return crd

def pdb2res(pdbfile,bychain=False,chaindef=[]):
    pdbcont = file(pdbfile)
    restype = {}
    for line in pdbcont:
        if line[:4]!='ATOM':
            continue

        if line[12:16].strip() == 'CA':
            res = int(line[22:26])
            chain = line[21]
            if line[26] != ' ':
                continue
            if chaindef != [] and chain not in chaindef:
                continue

            char = line[17:20]
            if char in ['HIP','HID','HIE']:
                char = 'HIS'
            elif char == 'CSS':
                char = 'CYS'

            if bychain:
                if chain not in restype:
                    restype[chain] = {}
                restype[chain][res] = char
            else:
                restype[res] = char
    return restype

#############################################################################
class COtype:
    def __init__(self,minR,atype1,atype2):
        self.minR = minR
        self.dcut = minR+0.5
        self.ds = []
        self.atype1 = atype1
        self.atype2 = atype2

class resCO:
    def __init__(self,ind,res1,res2,aa1,aa2):
        self.ind = ind
        self.res1 = res1
        self.res2 = res2
        self.aa1 = aa1
        self.aa2 = aa2

class DistributionClass:
    def __init__(self,pairtype):
        self.distr = []
        self.pairtype = pairtype
        self.n = 0

    def dat2distribution(self,datlist,normalize=True):
        self.distr=[0.0 for k in range(NBINS)]
        for i_dat in datlist:
            scale_dat=int((i_dat-MINVAL)/BINSIZE)

            if scale_dat >= NBINS:
                self.distr[-1] += 1.0
            #elif scale_dat < 0:
            #    return
            else:
                self.distr[scale_dat]+=1.0

        if normalize:
            vsum = sum(self.distr)
            self.distr = [float(val)/vsum for val in self.distr]
        self.n = len(datlist)

    def read_dstr_ref(self,infile):
        for line in file(infile):
            words = line[:-1].split()
            pairtype = words[0].strip()
            if pairtype != self.pairtype:
                continue
            self.n = int(words[1])
            self.distr = [float(word) for word in words[2:]]

    def smoothen(self):
        #weighting = [0.1,0.2,0.4,0.2,0.1]
        weighting = [0.04,0.11,0.21,0.28,0.21,0.11,0.04]
        unweighted_dstr = copy.copy(self.distr)
        self.distr = []

        for i,val in enumerate(unweighted_dstr):
            wsum = 0.0
            valsum = 0.0
            for k,w in enumerate(weighting):
                if i-3+k >= 0 and i-3+k < len(unweighted_dstr):
                    wsum += w
                    valsum += w*unweighted_dstr[i-3+k]
            valsum /= wsum
            self.distr.append(valsum)

    def deltasum(self,distr2,report=False):
        dsum = 0.0
        for i,val1 in enumerate(self.distr):
            val2 = distr2[i]
            dsum += abs(val1-val2)
            if report:
                print i, val1, val2, dsum
        return dsum

    def get_kldivergence(self,distr2):
        kldiv = 0.0
        for i,val1 in enumerate(self.distr):
            val1p = val1 + P_PSEUDO
            val2 = distr2[i] + P_PSEUDO
            kldiv += val1p*log(val1p/val2)
        return kldiv

    def report(self,report_header=False):
        if report_header:
            header = '%4s %4s %5s'%('Atm1','Atm2','N')
            for k in range(NBINS):
                header += ' %6.3f'%(MINVAL+k*BINSIZE)

        print '%8s %5d'%(self.pairtype,self.n)+' %9.7f'*len(self.distr)%tuple(self.distr)

def get_atmpairs():
    atmpairs = {}
    typeorder = []
    k = 0
    for i,type1 in enumerate(atmtypes):
        for j,type2 in enumerate(atmtypes[i:]):
            atmpairs[type1+'_'+type2] = k
            typeorder.append(type1+'_'+type2)
            if j > 0:
                atmpairs[type2+'_'+type1] = k
            k += 1
    return atmpairs, typeorder

def read_atmtypemap(resfile):
    atmtypemap = {}
    for line in file(resfile):
        if line.startswith('ATOM'):
            words = line[:-1].split()
            atmname = words[1]
            atype = words[2]
            if atype not in atmtypes : continue
            atmtypemap[atmname] = atype
    return atmtypemap

def read_radii(infile):
    radii = {}
    for line in file(infile):
        words = line[:-1].split()
        atype = words[0]
        if atype in atmtypes:
            radius = float(words[2])
            radii[atype] = radius
    return radii

def get_fullCOs(typeorder):
    COs = []
    radii = read_radii(ATOMFILE)
    for apair in typeorder:
        atype1,atype2 = apair.split('_')
        COs.append(COtype(radii[atype1]+radii[atype2],
                        atype1,atype2))
    return COs

def get_contacts(pdb):
    crd = pdb2crd(pdb,'CB',chaindef=['A'])
    aa = pdb2res(pdb,chaindef=['A'])
    reslist = crd.keys()
    cos = []

    ind = pdb.split('/')[-1].split('.')[0].replace('_0001','')
    dist = []
    for i,res1 in enumerate(reslist[:-MINSEP]):
        for res2 in reslist[i+MINSEP:]:
            dis2 = d2(crd[res1],crd[res2])
            if dis2 < CO_DIS2CUT:
                cos.append(resCO(ind,res1,res2,aa[res1],aa[res2]))
    return cos

def get_dist(crd,rescos,COs,atype,atmpairs):
    for co in rescos:
        for atm1 in crd[co.res1]:
            if atm1 not in atype[co.aa1]: continue
            for atm2 in crd[co.res2]:
                if atm2 not in atype[co.aa2]: continue
                a1 = atype[co.aa1][atm1]
                a2 = atype[co.aa2][atm2]

                ipair = atmpairs[a1+'_'+a2]
                d = distance(crd[co.res1][atm1],crd[co.res2][atm2])
                CO = COs[ipair]
                if d < CO.dcut:
                    CO.ds.append(d)

def run_pdb(pdb,COs,atype,atmpairs):
    crd = pdb2crd(pdb,'',chaindef=['A'])
    resco = get_contacts(pdb)
    get_dist(crd,resco,COs,atype,atmpairs)

def main(mode='run_terse'):
    direc = sys.argv[1]
    pdbs = glob.glob('%s/*.pdb'%(direc))

    atmpairs,typeorder = get_atmpairs()
    COs     = get_fullCOs(typeorder)
    KLsum = 0
    counttot = 0

    resatms = {}
    for aa in AAS:
        resfile = '%s/%s.params'%(DBPATH,aa)
        resatms[aa] = read_atmtypemap(resfile)

    for pdb in pdbs:
        run_pdb(pdb,COs,resatms,atmpairs)

    if mode == 'run':
        print '%12s %6s %8s %8s'%('#PairType','Ndat','SumDelta','KLDiv')

    for i,CO in enumerate(COs):
        if len(CO.ds) == 0: continue

        dstr = DistributionClass(typeorder[i])
        dstr.dat2distribution(CO.ds)
        dstr.smoothen()

        if mode == 'report':
            dstr.report(report_header=(i==0))
        elif (mode == 'run') or (mode == 'run_terse'):
            dstr_ref = DistributionClass(typeorder[i])
            dstr_ref.read_dstr_ref(REFFILE)

            deltasum = dstr_ref.deltasum(dstr.distr)
            kldiv = dstr_ref.get_kldivergence(dstr.distr)

            if mode == 'run':
	            print '%-12s %6d %8.5f %8.5f'%(typeorder[i],dstr_ref.n,deltasum,kldiv)
            if dstr_ref.n > 200:
                KLsum += kldiv * kldiv
                counttot += 1
    if (mode == 'run') or (mode == 'run_terse'):
        print '%f %f'%(sqrt(KLsum/counttot),counttot)

if __name__ == "__main__":
    main()
