'''
Created on Jun 26, 2017
for ELink linker,data from realtime computing,1:1 and only have 1 linker for each complex
@author: WG
'''
from __future__ import division
import re,pickle,fcntl
import numpy as np
import math,datetime,time
from collections import Iterable
from pyteomics import mgf,auxiliary
from operator import mul
from Bio import SeqIO
from multiprocessing import Pool


'''
mass of all aa
'''
aaMasses = { #without water, with ICH2CONH2-CH2CONH2
                         'A' :  71.0371103,
                         'C' : 160.0306443,
                         'D' : 115.0269385, 
                         'E' : 129.0425877, 
                         'F' : 147.0684087, 
                         'G' :  57.0214611, 
                         'H' : 137.0589059, 
                         'I' : 113.0840579, 
                         'K' : 128.0949557, 
                         'L' : 113.0840579, 
                         'M' : 131.0404787, 
                         'N' : 114.0429222, 
                         'P' :  97.0527595, 
                         'Q' : 128.0585714, 
                         'R' : 156.1011021, 
                         'S' :  87.0320244, 
                         'T' : 101.0476736, 
                         'V' :  99.0684087, 
                         'W' : 186.0793065, 
                         'Y' : 163.0633228, 
                         'O' : 207.29944, 
                         'U' : 150.04244
}
rK = re.compile('Y')
h = 1.007825
h2o = 18.010564
co = 27.994914
oh = 17.002739
nh = 15.010898
nh3 = 17.026548
nh2 = 16.018723


def digestseq(seqs,enyname,maxl,minl,miscut,ko):
    """ add enzyme target lable object.seqs = sequese, enyname means the name of enzyme,ko=1 means only save the fragment with Y.
    return the datesit of aa fragment,maxl means the max length of peptide want to save,minl means the min length, miscut means the max number of site to be miss.
    """
    spliseq = []
    spliseqk= []
    #AA=['A','R','N','D','C','E','Q','G','H','J','L','K','M','F','P','S','T','W','Y','V']
    if str(enyname) == 'Trypsin/p' :
        seqs = re.sub(r'RA', "ROA", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'RR', "ROR", str(seqs))
        seqs = re.sub(r'RN', "RON", str(seqs))
        seqs = re.sub(r'RD', "ROD", str(seqs))
        seqs = re.sub(r'RC', "ROC", str(seqs))
        seqs = re.sub(r'RE', "ROE", str(seqs))
        seqs = re.sub(r'RQ', "ROQ", str(seqs))
        seqs = re.sub(r'RG', "ROG", str(seqs))
        seqs = re.sub(r'RH', "ROH", str(seqs))
        seqs = re.sub(r'RJ', "ROJ", str(seqs))
        seqs = re.sub(r'RL', "ROL", str(seqs))
        seqs = re.sub(r'RK', "ROK", str(seqs))
        seqs = re.sub(r'RM', "ROM", str(seqs))
        seqs = re.sub(r'RF', "ROF", str(seqs))
        seqs = re.sub(r'RP', "ROP", str(seqs))
        seqs = re.sub(r'RS', "ROS", str(seqs))
        seqs = re.sub(r'RT', "ROT", str(seqs))
        seqs = re.sub(r'RW', "ROW", str(seqs))
        seqs = re.sub(r'RY', "ROY", str(seqs))
        seqs = re.sub(r'RV', "ROV", str(seqs))
        seqs = re.sub(r'KA', "KOA", str(seqs))#the cut site2 of enzyme
        seqs = re.sub(r'KR', "KOR", str(seqs))
        seqs = re.sub(r'KN', "KON", str(seqs))
        seqs = re.sub(r'KD', "KOD", str(seqs))
        seqs = re.sub(r'KC', "KOC", str(seqs))
        seqs = re.sub(r'KE', "KOE", str(seqs))
        seqs = re.sub(r'KQ', "KOQ", str(seqs))
        seqs = re.sub(r'KG', "KOG", str(seqs))
        seqs = re.sub(r'KH', "KOH", str(seqs))
        seqs = re.sub(r'KJ', "KOJ", str(seqs))
        seqs = re.sub(r'KL', "KOL", str(seqs))
        seqs = re.sub(r'KK', "KOK", str(seqs))
        seqs = re.sub(r'KM', "KOM", str(seqs))
        seqs = re.sub(r'KF', "KOF", str(seqs))
        seqs = re.sub(r'KP', "KOP", str(seqs))
        seqs = re.sub(r'KS', "KOS", str(seqs))
        seqs = re.sub(r'KT', "KOT", str(seqs))
        seqs = re.sub(r'KW', "KOW", str(seqs))
        seqs = re.sub(r'KY', "KOY", str(seqs))
        seqs = re.sub(r'KV', "KOV", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
    elif str(enyname) == 'Trypsin':
        seqs = re.sub(r'RA', "ROA", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'RR', "ROR", str(seqs))
        seqs = re.sub(r'RN', "RON", str(seqs))
        seqs = re.sub(r'RD', "ROD", str(seqs))
        seqs = re.sub(r'RC', "ROC", str(seqs))
        seqs = re.sub(r'RE', "ROE", str(seqs))
        seqs = re.sub(r'RQ', "ROQ", str(seqs))
        seqs = re.sub(r'RG', "ROG", str(seqs))
        seqs = re.sub(r'RH', "ROH", str(seqs))
        seqs = re.sub(r'RJ', "ROJ", str(seqs))
        seqs = re.sub(r'RL', "ROL", str(seqs))
        seqs = re.sub(r'RK', "ROK", str(seqs))
        seqs = re.sub(r'RM', "ROM", str(seqs))
        seqs = re.sub(r'RF', "ROF", str(seqs))
        #seqs = re.sub(r'RP', "ROP", str(seqs))
        seqs = re.sub(r'RS', "ROS", str(seqs))
        seqs = re.sub(r'RT', "ROT", str(seqs))
        seqs = re.sub(r'RW', "ROW", str(seqs))
        seqs = re.sub(r'RY', "ROY", str(seqs))
        seqs = re.sub(r'RV', "ROV", str(seqs))
        seqs = re.sub(r'KA', "KOA", str(seqs))#the cut site2 of enzyme
        seqs = re.sub(r'KR', "KOR", str(seqs))
        seqs = re.sub(r'KN', "KON", str(seqs))
        seqs = re.sub(r'KD', "KOD", str(seqs))
        seqs = re.sub(r'KC', "KOC", str(seqs))
        seqs = re.sub(r'KE', "KOE", str(seqs))
        seqs = re.sub(r'KQ', "KOQ", str(seqs))
        seqs = re.sub(r'KG', "KOG", str(seqs))
        seqs = re.sub(r'KH', "KOH", str(seqs))
        seqs = re.sub(r'KJ', "KOJ", str(seqs))
        seqs = re.sub(r'KL', "KOL", str(seqs))
        seqs = re.sub(r'KK', "KOK", str(seqs))
        seqs = re.sub(r'KM', "KOM", str(seqs))
        seqs = re.sub(r'KF', "KOF", str(seqs))
        #seqs = re.sub(r'KP', "KOP", str(seqs))
        seqs = re.sub(r'KS', "KOS", str(seqs))
        seqs = re.sub(r'KT', "KOT", str(seqs))
        seqs = re.sub(r'KW', "KOW", str(seqs))
        seqs = re.sub(r'KY', "KOY", str(seqs))
        seqs = re.sub(r'KV', "KOV", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
    elif str(enyname) == 'ArgC':
        seqs = re.sub(r'RA', "ROA", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'RR', "ROR", str(seqs))
        seqs = re.sub(r'RN', "RON", str(seqs))
        seqs = re.sub(r'RD', "ROD", str(seqs))
        seqs = re.sub(r'RC', "ROC", str(seqs))
        seqs = re.sub(r'RE', "ROE", str(seqs))
        seqs = re.sub(r'RQ', "ROQ", str(seqs))
        seqs = re.sub(r'RG', "ROG", str(seqs))
        seqs = re.sub(r'RH', "ROH", str(seqs))
        seqs = re.sub(r'RJ', "ROJ", str(seqs))
        seqs = re.sub(r'RL', "ROL", str(seqs))
        seqs = re.sub(r'RK', "ROK", str(seqs))
        seqs = re.sub(r'RM', "ROM", str(seqs))
        seqs = re.sub(r'RF', "ROF", str(seqs))
        seqs = re.sub(r'RP', "ROP", str(seqs))
        seqs = re.sub(r'RS', "ROS", str(seqs))
        seqs = re.sub(r'RT', "ROT", str(seqs))
        seqs = re.sub(r'RW', "ROW", str(seqs))
        seqs = re.sub(r'RY', "ROY", str(seqs))
        seqs = re.sub(r'RV', "ROV", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
    elif str(enyname) == 'ArgC/AspN':
        seqs = re.sub(r'RA', "ROA", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'RR', "ROR", str(seqs))
        seqs = re.sub(r'RN', "RON", str(seqs))
        seqs = re.sub(r'RD', "ROD", str(seqs))
        seqs = re.sub(r'RC', "ROC", str(seqs))
        seqs = re.sub(r'RE', "ROE", str(seqs))
        seqs = re.sub(r'RQ', "ROQ", str(seqs))
        seqs = re.sub(r'RG', "ROG", str(seqs))
        seqs = re.sub(r'RH', "ROH", str(seqs))
        seqs = re.sub(r'RJ', "ROJ", str(seqs))
        seqs = re.sub(r'RL', "ROL", str(seqs))
        seqs = re.sub(r'RK', "ROK", str(seqs))
        seqs = re.sub(r'RM', "ROM", str(seqs))
        seqs = re.sub(r'RF', "ROF", str(seqs))
        seqs = re.sub(r'RP', "ROP", str(seqs))
        seqs = re.sub(r'RS', "ROS", str(seqs))
        seqs = re.sub(r'RT', "ROT", str(seqs))
        seqs = re.sub(r'RW', "ROW", str(seqs))
        seqs = re.sub(r'RY', "ROY", str(seqs))
        seqs = re.sub(r'RV', "ROV", str(seqs))
        seqs = re.sub(r'AD', "AOD", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'RD', "ROD", str(seqs))
        seqs = re.sub(r'ND', "NOD", str(seqs))
        seqs = re.sub(r'DD', "DOD", str(seqs))
        seqs = re.sub(r'CD', "COD", str(seqs))
        seqs = re.sub(r'ED', "EOD", str(seqs))
        seqs = re.sub(r'QD', "QOD", str(seqs))
        seqs = re.sub(r'GD', "GOD", str(seqs))
        seqs = re.sub(r'HD', "HOD", str(seqs))
        seqs = re.sub(r'JD', "JOD", str(seqs))
        seqs = re.sub(r'LD', "LOD", str(seqs))
        seqs = re.sub(r'KD', "KOD", str(seqs))
        seqs = re.sub(r'MD', "MOD", str(seqs))
        seqs = re.sub(r'FD', "FOD", str(seqs))
        seqs = re.sub(r'PD', "POD", str(seqs))
        seqs = re.sub(r'SD', "SOD", str(seqs))
        seqs = re.sub(r'TD', "TOD", str(seqs))
        seqs = re.sub(r'WD', "WOD", str(seqs))
        seqs = re.sub(r'YD', "YOD", str(seqs))
        seqs = re.sub(r'VD', "VOD", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
    elif str(enyname) == 'AspC':
        seqs = re.sub(r'DA', "DOA", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'DR', "DOR", str(seqs))
        seqs = re.sub(r'DN', "DON", str(seqs))
        seqs = re.sub(r'DD', "DOD", str(seqs))
        seqs = re.sub(r'DC', "DOC", str(seqs))
        seqs = re.sub(r'DE', "DOE", str(seqs))
        seqs = re.sub(r'DQ', "DOQ", str(seqs))
        seqs = re.sub(r'DG', "DOG", str(seqs))
        seqs = re.sub(r'DH', "DOH", str(seqs))
        seqs = re.sub(r'DJ', "DOJ", str(seqs))
        seqs = re.sub(r'DL', "DOL", str(seqs))
        seqs = re.sub(r'DK', "DOK", str(seqs))
        seqs = re.sub(r'DM', "DOM", str(seqs))
        seqs = re.sub(r'DF', "DOF", str(seqs))
        seqs = re.sub(r'DP', "DOP", str(seqs))
        seqs = re.sub(r'DS', "DOS", str(seqs))
        seqs = re.sub(r'DT', "DOT", str(seqs))
        seqs = re.sub(r'DW', "DOW", str(seqs))
        seqs = re.sub(r'DY', "DOY", str(seqs))
        seqs = re.sub(r'DV', "DOV", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
    elif str(enyname) == 'GluC':
        seqs = re.sub(r'EA', "EOA", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'ER', "EOR", str(seqs))
        seqs = re.sub(r'EN', "EON", str(seqs))
        seqs = re.sub(r'ED', "EOD", str(seqs))
        seqs = re.sub(r'EC', "EOC", str(seqs))
        seqs = re.sub(r'EE', "EOE", str(seqs))
        seqs = re.sub(r'EQ', "EOQ", str(seqs))
        seqs = re.sub(r'EG', "EOG", str(seqs))
        seqs = re.sub(r'EH', "EOH", str(seqs))
        seqs = re.sub(r'EJ', "EOJ", str(seqs))
        seqs = re.sub(r'EL', "EOL", str(seqs))
        seqs = re.sub(r'EK', "EOK", str(seqs))
        seqs = re.sub(r'EM', "EOM", str(seqs))
        seqs = re.sub(r'EF', "EOF", str(seqs))
        seqs = re.sub(r'EP', "EOP", str(seqs))
        seqs = re.sub(r'ES', "EOS", str(seqs))
        seqs = re.sub(r'ET', "EOT", str(seqs))
        seqs = re.sub(r'EW', "EOW", str(seqs))
        seqs = re.sub(r'EY', "EOY", str(seqs))
        seqs = re.sub(r'EV', "EOV", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
    elif str(enyname) == 'GluN':
        seqs = re.sub(r'AE', "AOE", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'RE', "ROE", str(seqs))
        seqs = re.sub(r'NE', "NOE", str(seqs))
        seqs = re.sub(r'DE', "DOE", str(seqs))
        seqs = re.sub(r'CE', "COE", str(seqs))
        seqs = re.sub(r'EE', "EOE", str(seqs))
        seqs = re.sub(r'QE', "QOE", str(seqs))
        seqs = re.sub(r'GE', "GOE", str(seqs))
        seqs = re.sub(r'HE', "HOE", str(seqs))
        seqs = re.sub(r'JE', "JOE", str(seqs))
        seqs = re.sub(r'LE', "LOE", str(seqs))
        seqs = re.sub(r'KE', "KOE", str(seqs))
        seqs = re.sub(r'ME', "MOE", str(seqs))
        seqs = re.sub(r'FE', "FOE", str(seqs))
        seqs = re.sub(r'PE', "POE", str(seqs))
        seqs = re.sub(r'SE', "SOE", str(seqs))
        seqs = re.sub(r'TE', "TOE", str(seqs))
        seqs = re.sub(r'WE', "WOE", str(seqs))
        seqs = re.sub(r'YE', "YOE", str(seqs))
        seqs = re.sub(r'VE', "VOE", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
    elif str(enyname) == 'AspN':
        '''
        Metalloprotease that hydrolyzes peptide bonds at the amino side of Asp and cysteic acid. If cysteine is reduced or alkylated, only -c-Asp-X is cleaved.
        '''
        seqs = re.sub(r'AD', "AOD", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'RD', "ROD", str(seqs))
        seqs = re.sub(r'ND', "NOD", str(seqs))
        seqs = re.sub(r'DD', "DOD", str(seqs))
        seqs = re.sub(r'CD', "COD", str(seqs))
        seqs = re.sub(r'ED', "EOD", str(seqs))
        seqs = re.sub(r'QD', "QOD", str(seqs))
        seqs = re.sub(r'GD', "GOD", str(seqs))
        seqs = re.sub(r'HD', "HOD", str(seqs))
        seqs = re.sub(r'JD', "JOD", str(seqs))
        seqs = re.sub(r'LD', "LOD", str(seqs))
        seqs = re.sub(r'KD', "KOD", str(seqs))
        seqs = re.sub(r'MD', "MOD", str(seqs))
        seqs = re.sub(r'FD', "FOD", str(seqs))
        seqs = re.sub(r'PD', "POD", str(seqs))
        seqs = re.sub(r'SD', "SOD", str(seqs))
        seqs = re.sub(r'TD', "TOD", str(seqs))
        seqs = re.sub(r'WD', "WOD", str(seqs))
        seqs = re.sub(r'YD', "YOD", str(seqs))
        seqs = re.sub(r'VD', "VOD", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
    elif str(enyname) == 'LysN':
        seqs = re.sub(r'AK', "AOK", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'RK', "ROK", str(seqs))
        seqs = re.sub(r'NK', "NOK", str(seqs))
        seqs = re.sub(r'DK', "DOK", str(seqs))
        seqs = re.sub(r'CK', "COK", str(seqs))
        seqs = re.sub(r'EK', "EOK", str(seqs))
        seqs = re.sub(r'QK', "QOK", str(seqs))
        seqs = re.sub(r'GK', "GOK", str(seqs))
        seqs = re.sub(r'HK', "HOK", str(seqs))
        seqs = re.sub(r'JK', "JOK", str(seqs))
        seqs = re.sub(r'LK', "LOK", str(seqs))
        seqs = re.sub(r'KK', "KOK", str(seqs))
        seqs = re.sub(r'MK', "MOK", str(seqs))
        seqs = re.sub(r'FK', "FOK", str(seqs))
        seqs = re.sub(r'PK', "POK", str(seqs))
        seqs = re.sub(r'SK', "SOK", str(seqs))
        seqs = re.sub(r'TK', "TOK", str(seqs))
        seqs = re.sub(r'WK', "WOK", str(seqs))
        seqs = re.sub(r'YK', "YOK", str(seqs))
        seqs = re.sub(r'VK', "VOK", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
    elif str(enyname) == 'Chymotrypsin+':
        seqs = re.sub(r'LA', "LOA", str(seqs))#the cut site of enzyme
        seqs = re.sub(r'LR', "LOR", str(seqs))
        seqs = re.sub(r'LN', "LON", str(seqs))
        seqs = re.sub(r'LD', "LOD", str(seqs))
        seqs = re.sub(r'LC', "LOC", str(seqs))
        seqs = re.sub(r'LE', "LOE", str(seqs))
        seqs = re.sub(r'LQ', "LOQ", str(seqs))
        seqs = re.sub(r'LG', "LOG", str(seqs))
        seqs = re.sub(r'LH', "LOH", str(seqs))
        seqs = re.sub(r'LJ', "LOJ", str(seqs))
        seqs = re.sub(r'LL', "LOL", str(seqs))
        seqs = re.sub(r'LK', "LOK", str(seqs))
        seqs = re.sub(r'LM', "LOM", str(seqs))
        seqs = re.sub(r'LF', "LOF", str(seqs))
        seqs = re.sub(r'LP', "LOP", str(seqs))
        seqs = re.sub(r'LS', "LOS", str(seqs))
        seqs = re.sub(r'LT', "LOT", str(seqs))
        seqs = re.sub(r'LW', "LOW", str(seqs))
        seqs = re.sub(r'LY', "LOY", str(seqs))
        seqs = re.sub(r'LV', "LOV", str(seqs))
        seqs = re.sub(r'MA', "MOA", str(seqs))#the cut site2 of enzyme
        seqs = re.sub(r'MR', "MOR", str(seqs))
        seqs = re.sub(r'MN', "MON", str(seqs))
        seqs = re.sub(r'MD', "MOD", str(seqs))
        seqs = re.sub(r'MC', "MOC", str(seqs))
        seqs = re.sub(r'ME', "MOE", str(seqs))
        seqs = re.sub(r'MQ', "MOQ", str(seqs))
        seqs = re.sub(r'MG', "MOG", str(seqs))
        seqs = re.sub(r'MH', "MOH", str(seqs))
        seqs = re.sub(r'MJ', "MOJ", str(seqs))
        seqs = re.sub(r'ML', "MOL", str(seqs))
        seqs = re.sub(r'MK', "MOK", str(seqs))
        seqs = re.sub(r'MM', "MOM", str(seqs))
        seqs = re.sub(r'MF', "MOF", str(seqs))
        seqs = re.sub(r'MP', "MOP", str(seqs))
        seqs = re.sub(r'MS', "MOS", str(seqs))
        seqs = re.sub(r'MT', "MOT", str(seqs))
        seqs = re.sub(r'MW', "MOW", str(seqs))
        seqs = re.sub(r'MY', "MOY", str(seqs))
        seqs = re.sub(r'MV', "MOV", str(seqs))
        seqs = re.sub(r'FA', "FOA", str(seqs))#the cut site3 of enzyme
        seqs = re.sub(r'FR', "FOR", str(seqs))
        seqs = re.sub(r'FN', "FON", str(seqs))
        seqs = re.sub(r'FD', "FOD", str(seqs))
        seqs = re.sub(r'FC', "FOC", str(seqs))
        seqs = re.sub(r'FE', "FOE", str(seqs))
        seqs = re.sub(r'FQ', "FOQ", str(seqs))
        seqs = re.sub(r'FG', "FOG", str(seqs))
        seqs = re.sub(r'FH', "FOH", str(seqs))
        seqs = re.sub(r'FJ', "FOJ", str(seqs))
        seqs = re.sub(r'FL', "FOL", str(seqs))
        seqs = re.sub(r'FK', "FOK", str(seqs))
        seqs = re.sub(r'FM', "FOM", str(seqs))
        seqs = re.sub(r'FF', "FOF", str(seqs))
        seqs = re.sub(r'FP', "FOP", str(seqs))
        seqs = re.sub(r'FS', "FOS", str(seqs))
        seqs = re.sub(r'FT', "FOT", str(seqs))
        seqs = re.sub(r'FW', "FOW", str(seqs))
        seqs = re.sub(r'FY', "FOY", str(seqs))
        seqs = re.sub(r'FV', "FOV", str(seqs))
        seqs = re.sub(r'WA', "WOA", str(seqs))#the cut site4 of enzyme
        seqs = re.sub(r'WR', "WOR", str(seqs))
        seqs = re.sub(r'WN', "WON", str(seqs))
        seqs = re.sub(r'WD', "WOD", str(seqs))
        seqs = re.sub(r'WC', "WOC", str(seqs))
        seqs = re.sub(r'WE', "WOE", str(seqs))
        seqs = re.sub(r'WQ', "WOQ", str(seqs))
        seqs = re.sub(r'WG', "WOG", str(seqs))
        seqs = re.sub(r'WH', "WOH", str(seqs))
        seqs = re.sub(r'WJ', "WOJ", str(seqs))
        seqs = re.sub(r'WL', "WOL", str(seqs))
        seqs = re.sub(r'WK', "WOK", str(seqs))
        seqs = re.sub(r'WM', "WOM", str(seqs))
        seqs = re.sub(r'WF', "WOF", str(seqs))
        seqs = re.sub(r'WP', "WOP", str(seqs))
        seqs = re.sub(r'WS', "WOS", str(seqs))
        seqs = re.sub(r'WT', "WOT", str(seqs))
        seqs = re.sub(r'WW', "WOW", str(seqs))
        seqs = re.sub(r'WY', "WOY", str(seqs))
        seqs = re.sub(r'WV', "WOV", str(seqs))
        seqs = re.sub(r'YA', "YOA", str(seqs))#the cut site5 of enzyme
        seqs = re.sub(r'YR', "YOR", str(seqs))
        seqs = re.sub(r'YN', "YON", str(seqs))
        seqs = re.sub(r'YD', "YOD", str(seqs))
        seqs = re.sub(r'YC', "YOC", str(seqs))
        seqs = re.sub(r'YE', "YOE", str(seqs))
        seqs = re.sub(r'YQ', "YOQ", str(seqs))
        seqs = re.sub(r'YG', "YOG", str(seqs))
        seqs = re.sub(r'YH', "YOH", str(seqs))
        seqs = re.sub(r'YJ', "YOJ", str(seqs))
        seqs = re.sub(r'YL', "YOL", str(seqs))
        seqs = re.sub(r'YK', "YOK", str(seqs))
        seqs = re.sub(r'YM', "YOM", str(seqs))
        seqs = re.sub(r'YF', "YOF", str(seqs))
        seqs = re.sub(r'YP', "YOP", str(seqs))
        seqs = re.sub(r'YS', "YOS", str(seqs))
        seqs = re.sub(r'YT', "YOT", str(seqs))
        seqs = re.sub(r'YW', "YOW", str(seqs))
        seqs = re.sub(r'YY', "YOY", str(seqs))
        seqs = re.sub(r'YV', "YOV", str(seqs))
        p = re.compile('O')
        se = str(seqs)
        spliseq = p.split(se)
        #start create the missed cleavages sequence
        misspl=spliseq
        msd=[]
        for i in range(miscut):
            wds=[]
            for z in range(len(misspl)):
                if (z+1)==len(misspl):
                    break
                if i < 1:
                    wds.append(misspl[z]+misspl[z+1])
                else:
                    zl=misspl[z]+spliseq[z+1+i]
                    wds.append(zl)
            for u in range(len(wds)):
                msd.append(wds[u])
            misspl=wds
        for i in range(len(msd)):
            spliseq.append(msd[i])
            
    reseql=[]
    if ko == 1:#get the sequence which have Y
        ps = re.compile('Y')
        for sseq in spliseq:
            match = ps.search(sseq)
            if match:
                spliseqk.append(sseq)
        spliseq=[]
        spliseq=spliseqk
    #remove the sequence which too long or too short,and remove unclear AA
    for yy in spliseq:
        
        if minl<len(yy)<maxl and 'B' not in yy and 'X' not in yy and 'J' not in yy and 'Z' not in yy:
            reseql.append(yy)

    #print(spliseq)
    return reseql
def SEQTOMASS(seq):
    '''
    calculate the mass of peptide
    '''
    M=0
    for aa in seq:
        M+=aaMasses[aa]
    return M+18.010564
def SavePdata(sarray,namef):
    '''
    save the all result to pickle, includ 12 parameters: 1. index; 2. file name; 3. linked type(1:1,1:2,2:1); 4. e charge; 5. modify; 6. number of linker; 7.sequence(ssZssXssZss);
    8. number of matched; 9. p value; 10. FDR; 11. protein or peptide name; 11. position of peptied in protein.  
    '''
    with open(namef,'ab') as f:
        fcntl.flock(f.fileno(),fcntl.LOCK_EX)
        pickle.dump(sarray,f)
        #print('save to scanconb file..'),protocol=pickle.HIGHEST_PROTOCOL
    return
def SaveWdata(sarray,namef):
    '''
    save the all result to pickle, includ 12 parameters: 1. index; 2. file name; 3. linked type(1:1,1:2,2:1); 4. e charge; 5. modify; 6. number of linker; 7.sequence(ssZssXssZss);
    8. number of matched; 9. p value; 10. FDR; 11. protein or peptide name; 11. position of peptied in protein.  
    '''
    with open(namef,'a+') as f:
        f.writelines(str(sarray))
        #print('save to scanconb file..')
    return

#def SaveJdata(sarray,namef):
    '''
    save the all result with joblib, includ 12 parameters: 1. index; 2. file name; 3. linked type(1:1,1:2,2:1); 4. e charge; 5. modify; 6. number of linker; 7.sequence(ssZssXssZss);
    8. number of matched; 9. p value; 10. FDR; 11. protein or peptide name; 11. position of peptied in protein.  
    '''
    #joblib.dump(sarray,filename=namef)

def InitA(spct):
    outpep= open("pepuseq.data", 'rb')
    outpro= open("prouseq.data", 'rb')
    outprom = open("prom.data", 'rb')
    outpepm = open("pepm.data", 'rb')
    outpepdic = open("pepdic.data", 'rb')
    link=[593.28225, 369.15494]#the mass of linker after linked, after MS-cut
    pepuseq=pickle.load(outpep)
    outpep.close()
    prouseq=pickle.load(outpro)#alnopep for complex; protein for interaction
    outpro.close()
    mpep=pickle.load(outpepm)
    outpepm.close()
    mpro=pickle.load(outprom)
    outprom.close()
    dict_pep=pickle.load(outpepdic)
    outpepdic.close()
    ppm=10
    maxms1charge=-1#max number of charges been considered in MS1, -1 means use the charge number in mgf file 
    msmt=0.1
    #the resized errormsmt=0.1
    #wl=len(mpep)
    if maxms1charge <=0 :
        permass=spct[0]*int(spct[3])-int(spct[3])*h# permass * e
        #get the mass for MS1
        bego=[]
        begocm=[]
        for linkv in link:
            #different linker mass
            #z=0
            for k,v in dict_pep.items():
                mmp=permass-v-linkv
                idd=[i for i,x in enumerate(mpro) if x<=mmp+msmt and x>=mmp-msmt]
                #iddpe=[i for i,x in enumerate(mpep[z:wl]) if x<=mmp+msmt and x>=mmp-msmt] the #marked only for conplex
                if idd != []:
                    for i in idd:
                        bego.append(prouseq[i]+'J'+k)
                        begocm.append(round((mpro[i]+v+linkv),7))
                #if iddpe != []:
                    #for i in iddpe:
                        #bego.append(pepuseq[i+z]+'J'+k)
                        #begocm.append(mpep[i+z]+v+linkv)
                #z+=1
        if len(bego) !=0:
            fil=[]
            fil.append(list(bego))
            fil.append(list(begocm))
            fil.append(list(spct[1]))#ms2list
            fil.append(list(spct[2]))#3 ms2intenslist
            fil.append(str(spct[4]))#rttime
            fil.append(int(spct[0]))#5 permass
            fil.append(int(spct[3]))#num_charge
            fil.append(str(spct[5]))#7 spectra_ref
            time.sleep(1)
            print('save to scanconb file..')
            SavePdata(fil,'scanconb')
            #del fil
        #else:
            #print('no match found...')
    else:
        #consider the diffent e charge for MS1
        bego=[]
        begocm=[]
        eseq=[]#the charge of hited MS1
        for ie in range(maxms1charge):
            permass=spct[0]*int(ie+1)-int(ie+1)*h# permass * e
            #get the mass for MS1
            for linkv in link:
                #z=0
                for k,v in dict_pep.items():
                    mmp=permass-v-linkv#remove the mass of linker and peptide, look for the paired protein
                    idd=[i for i,x in enumerate(mpro) if x<=mmp+msmt and x>=mmp-msmt]
                    #iddpe=[i for i,x in enumerate(mpep[z:wl]) if x<=mmp+msmt and x>=mmp-msmt] the #marked only for conplex
                    if idd != []:
                        for i in idd:
                            bego.append(prouseq[i]+'J'+k)
                            begocm.append(round((mpro[i]+v+linkv),7))
                            eseq.append(ie+1)
                    #if iddpe != []:
                        #for i in iddpe:
                            #bego.append(pepuseq[i+z]+'J'+k)
                            #begocm.append(mpep[i+z]+v+linkv)
                    #z+=1
        if len(bego) !=0:
            fil=[]
            fil.append(list(bego))
            fil.append(list(begocm))
            fil.append(list(spct[1]))#ms2list
            fil.append(list(spct[2]))#3 ms2intenslist
            fil.append(str(spct[4]))#rttime
            fil.append(int(spct[0]))#5 permass
            fil.append(list(eseq))#num_charge
            fil.append(str(spct[5]))#7 spectra_ref
            time.sleep(1)
            print('save to scanconb file..')
            SavePdata(fil,'scanconb')
            #del fil
        else:
            print('no match found...')
    return



if __name__=='__main__':

    mgff='2-allmod.mgf'#the 
    fhpro = open("uniprot-proteome-rat-UP000002494_fltd.fasta")#the protein fasta data
    fhpep = open("shanghai_sv_mature_200_mix_fltd.fasta")#the peptide fasta data
    # Name of the output file
    outpro = open("prouseq.data", 'wb')#the protein digested
    outpep = open("pepuseq.data", 'wb')#the peptides digested
    outprom = open("prom.data", 'wb')#the mass of digested proteins
    outpepm = open("pepm.data", 'wb')#the mass of digested peptides
    outprodic = open("prodic.data", 'wb')#the diction of digested protein
    outpepdic = open("pepdic.data", 'wb')#the diction of digested peptide
    
    link=[593.28225, 369.15494]#the mass of linker after linked, after MS-cut
    ppm2=10#ppm of ms2
    ppm=20#ppm of ms1
    maxms1charge=-1#max number of charges been considered in MS1, -1 means use the charge number in mgf file 
    proteinsq=[]
    peptidesq=[]
    #built the digest sequence fragments
    print('digest the seq data...')
    for pro in SeqIO.parse(fhpro,'fasta'):
        dseq=digestseq(pro.seq,'Trypsin',25,4,3,0)
        proteinsq.append(dseq)
        
    fhpro.close()
    for pep in SeqIO.parse(fhpep,'fasta'):
        peptidesq.append(digestseq(pep.seq,'Trypsin',25,4,3,1))
    fhpep.close()
    #make it flatten
    prosqfl=[]
    pepsqfl=[]
    for id1 in proteinsq:
        for id2 in id1:
            prosqfl.append(id2)
    for id1 in peptidesq:
        for id2 in id1:
            pepsqfl.append(id2)
    
    print(len(prosqfl),len(pepsqfl))
    print('remove the repeat value..')
    prouseq=list(set(prosqfl))
    pepuseq=list(set(pepsqfl))
    alnopep=list(set(prouseq)^set(pepuseq))#Symmetric difference set (unique to both) for complex,change the prouseq below to alnopep
    print(len(prouseq),len(pepuseq),len(alnopep))
    #bulit diction
    dict_pro={}
    mpro=[]
    dict_pep={}
    mpep=[]
    for zz in prouseq:
        dict_pro[zz]=SEQTOMASS(zz)
        mpro.append(SEQTOMASS(zz))
    for zz in pepuseq:
        dict_pep[zz]=SEQTOMASS(zz)
        mpep.append(SEQTOMASS(zz))
    
    
    pickle.dump(prouseq,outpro)#the nonK protein,when complex change it to alnopep
    pickle.dump(pepuseq,outpep)#the peptide with K
    pickle.dump(dict_pro,outprodic)#the nonK protein diction
    pickle.dump(dict_pep,outpepdic)#the peptide with K diction
    pickle.dump(mpro,outprom)#the mass of nonK protein
    pickle.dump(mpep,outpepm)#the mass of peptide with K
    outpro.close()
    outpep.close()
    outprom.close()
    outpepm.close()
    outprodic.close()
    outpepdic.close()
    '''
    mgff = the file of mgf peak list (eg. 'D:/TPP/data/OR20080317_S_SILAC-LH_1-1_01.mgf'); rfhpro = the proteins data with fasta format; rfhpep = the peptides data with fasta fromat
    link = the mass of linker after reaction; ppm2 = the errors of ms2; ppm = the errors of ms1; maxms1charge = the max number of charge to consider, wile -1 means use the mgf record number; maxms1charge = the max number of charge to consider when ms2 list
    '''
    header = mgf.read(mgff)
    print('scan the mgf file...')
    mblink=[226.14298,181.14532,191.12967,210.12425,209.14024]#the list of ELink fragments in MSMS
    zzc=0#the number of spectrum which has linker characteristic peak
    zc=0
    erm=0.5#the error range of round-off
    ctof=1#the min msms fragments needed for cut off  
    mops=[]#the data for pool
    mopsct=[]
    for spectrum in header:
        z=0
        expmass=spectrum['params']['pepmass'][0]#0
        ms2list=spectrum['m/z array']#1
        ms2intenslist=spectrum['intensity array']#2
        try:
            num_charge=spectrum['params']['charge'][0]#3
        except:
            num_charge=1
        rttime=spectrum['params']['rtinseconds']#4
        spectra_ref=spectrum['params']['title']#5
        bb=np.around(ms2list, decimals=4)
        tg=0
        zc+=1#the number of ms spectrum
        for a in mblink:
            for b in bb:
                if a<b+erm+1 and a>b-erm-1:
                    tg+=1
        if tg >= ctof:
            #list the MSdata for search
            zzc+=1
            mssct=[]
            mssct.append(expmass)#0
            mssct.append(ms2list)#1
            mssct.append(ms2intenslist)#2
            mssct.append(num_charge)#3
            mssct.append(rttime)#4
            mssct.append(spectra_ref)#5
            
            mopsct.append(mssct)
            
    print('listing the combination...')
    beg=time.time()
    p=Pool(40)
    p.map(InitA,mopsct)
    end=time.time()
    SavePdata(zzc, 'zzc.pkl')
    p.close()
    p.join()
    print(end-beg,'virtual combination finished!')
