#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:30:33 2022

@author: hilmi
"""

import pandas as pd
import statsmodels.api as sm
import numpy as np
import itertools
import time
import os,sys
import argparse
import scipy.stats
from decimal import Decimal

#sys.path.append('../localpackage')
import bio

#python3 olskmer.py test-in/Mus_musculus\|NA\|Unpublished\|Zfp24.txt test-out -k 3 -d 1 -g 1 -p 1

from sklearn import linear_model

nucleotides = ['A','C','G','T']

def adjustscr(score,shift=1000):  # Log transformation yaparak normalize etme  
    # log trans işlemi, yüksek değerdeki pbmdeki fluorescent sinyal miktarlarını aşağıya çekmek için kullanılır böylece linear model performansını artırılır
    minscr = score[score.idxmin()] # pbm datasetindeki min değeri 
    score_frame = pd.DataFrame(np.log(score - minscr + shift),columns=['score']) # log scaling   Log (SCORE - Minimum score + shift miktarı)
    return score_frame # shift değeri min değerini sıfırdan kurtarmak için eklenilen pay

def read_pbm(filename,kmer,nonrev_list,gappos=0,gapsize=0): # reading PBM input and apply transformation for scores and binary seq
    tbl = pd.read_csv(filename,names=['score','sequence'],delim_whitespace=True) #score,sequence ..
    score = adjustscr(tbl['score'],1000) # log transformation for fluorescent signals
    seqbin = [bio.seqtoi(x) for x in tbl['sequence']] #  PBM içindeki her bir sekansı binary gösterimine çevirir 
    oligfreq = bio.nonr_olig_freq(seqbin,kmer,nonrev_list,gappos=gappos,gapsize=gapsize) # feature vs sekans içeren count table oluşturur
    return pd.concat([score,oligfreq],axis=1)

def print_full(x):
    y = pd.DataFrame(x)
    #y = pd.DataFrame([bio.itoseq(z) for z in x])
    pd.set_option('display.max_rows', len(y))
    y.columns = [bio.itoseq(x) for x in y.columns]
    print(y)
    pd.reset_option('display.max_rows')

if __name__ == "__main__": # python olskmer.py $filein $output -k $kmer -d $chunk
    parser = argparse.ArgumentParser(description='Generate predicitons for all possible k-mer permutations.')
    parser.add_argument('input', type=str, nargs=2, help="Input for the script, format: <tf-path> <out-path>")
    parser.add_argument('-k','--kmer', type=int, help="The value of k that we want", default=6)
    parser.add_argument('-d','--div', type=int, help="Sometimes data can be too big, div is the number of chunks required", default=32)
    parser.add_argument('-g','--gapsize', type=int, help="The number of gaps", default=0)
    parser.add_argument('-p','--gappos', type=int, help="The position of gaps", default=0)
    args = parser.parse_args()

    tfpath = args.input[0]
    outpath = args.input[1]

    filename = os.path.splitext(os.path.basename(tfpath))[0]
    print("Got input file: " + filename)
    k = 6
    nonrev_list = bio.gen_nonreversed_kmer(k)
    start_time = time.time()
    df = read_pbm("/home/hilmi/qbic/PBM_inputs/Homo_sapiens_M00497_1.94d_Zoo_01_3027.txt",k,nonrev_list,0,0) # for k=3
    print(df)
    print("---Time to read pbm file: %s seconds ---" % (time.time() - start_time)) #5minutes
    lm = sm.OLS(df['score'],df.drop('score',axis=1)).fit()
    
    print(lm.summary()) # feature'ların ağırlıkları ve OLS summary statistics
    print(min(lm.params), max(lm.params)) # min ve max feature ağırlıkları
    bio.itoseq(6677) #TAATTA
    print("--- Finish training models: %s seconds ---" % (time.time() - start_time))
    #mutated context
    mutated = dict() # creating empty mutation dictionary,covering all possible mutations in all possible contexts
    #gapmer = args.kmer + args.gapsize
    for base in bio.nucleotides: # each is sequence of length 2k-1 with the mutated in the middle
        mutated[base] = [bio.insert_pos(seq,base,k-1) for seq in bio.seq_permutation(2*(k-1))] # for k =3 
    # bio.seq_permutation(4) --> range of seq represented in binary 
    # bio.seq_permutation(2*(6-1)) # all possible permutation of 2*(k-1) Nucleotides
    # # we have 1048576 combinations from 4**10 for k = 6
    # # bio.itoseq(1048576) --> AAAAA*AAAAA    
    # bio.itoseq((bio.insert_pos(1048576, "A", 5))) # mutasyonu ortada sabit tutup, yandakileri değiştiriyor
    # bio.itoseq((bio.insert_pos(1048577, "A", 5))) 
    # bio.itoseq(2097151) #--> TTTTT*TTTTT 
    # bio.itoseq((bio.insert_pos(2097151, "A", 5))) # This is last element of "A" nucleotide key in the dictionary
   
    # #For center A
    # bio.itoseq(mutated["A"][0]) # 11 len in bp
    # AAAAA*A*AAAAA 
    # bio.itoseq(mutated["A"][1])
    # AAAAA*A*AAAAC # middle one* is the mutated one that is inserted (total length = 10+1= 11 for k=6)
    # bio.itoseq(mutated["A"][2])
    # bio.itoseq(mutated["A"][4])
    # bio.itoseq(mutated["A"][-1])
    # #For center C
    # bio.itoseq(mutated["C"][0])
    # AAAAA*C*AAAAA  
    # bio.itoseq(mutated["C"][1])
    # bio.itoseq(mutated["C"][2])
    # bio.itoseq(mutated["C"][3])
    # bio.itoseq(mutated["C"][-1])
    # #For T
    # bio.itoseq(mutated["T"][0])
    # bio.itoseq(mutated["T"][-1])
    # #For G
    # bio.itoseq(mutated["G"][0])
    # bio.itoseq(mutated["G"][-1])
    # mutated dictionary includes --> 1048576 (4**10) possible mutataion context for each bases(A,C,T,G)
    print(len(mutated["A"])) # for k=3 , there are 4**4 possible context for each base
    div= 128 
    len(mutated["A"]) == len(bio.nucleotides)**(2*(6-1))
    chunk = len(bio.nucleotides)**(2*(k-1)) // div # 1048576 context --> 128 * 8192-chunk for k = 6
    # for k = 6 --> 1048576 context --> 128* 8192(chunk)
    print("Chunk size: " + str(chunk))
    output_all = pd.DataFrame(columns=['dna_seq','diff','t'])
    for i in range(0,div):
        i = 0
        print("Processing chunk-"+str(i))
        mutated_part = dict()
        count = dict()
        # count oligonucleotides in all k-mer combinations of sequence length k with mutated in the middle.
        # row = all sequences input, col = all possible permutations
        # Therefore we have 4**2k-1 rows (i.e. combinations)
        for base in bio.nucleotides:
            mutated_part[base] = mutated[base][i*chunk:(i+1)*chunk] #16384 # for i =0 mutated["A"][0:8192], mutated["A"][8192:8192+8192]
            count[base] = bio.nonr_olig_freq(mutated_part[base],k,nonrev_list,gappos=0,gapsize=0) #count the frequency of kmer in the mutation ##save this?
            #count --> contains #div rows x #features for each mutated base 
            # for k == 6: 1048576 mutated context per base x 2080 features 
        # dff = pd.DataFrame(count['A'])
        # dff['xyz'] = [bio.itoseq(ddd) for ddd in mutated_part["A"]] 
        # print_full(dff) # TTTATTT = 4095 (last element)
        diff_count = dict()  # c'(c.mutation-c.wild, sıklık farkı)      
        diff = dict() # c'B --> dot(c',weights)
        for b1 in nucleotides:  # A,C G T # her bir olası mutasyonel countu farkını simule etmiş oluyoruz
            single = []
            # AAAAA A AAAAA --> A wild type (reference)
            # AAAAA C AAAAA --> mutated to C (variant)
            # count["C"][:1] - count["A"][:1] --> COUNT DIFF
            for b2 in nucleotides: # FOR 1. A-C A-G- A-T >> 2. C-A, C-G,C-T >> 3. G-A, G-C, G-T >> 4. T-A,T-C,T-G
                if b1 != b2: # evaluate only when there is a base change
                    single.append(count[b2] - count[b1]) # middle base change 
                    # TTTTT*G*TTTTT - TTTTT*T*TTTTT  --> diff
            diff_count[b1] = pd.concat(single,axis=0,ignore_index=True) # merge these chunk*3 rows 
            diff[b1] = np.dot(diff_count[b1],lm.params) #dim: (49152, 2080) (2080,) .. diff_count == c' in the paper?
        # len(diff["A"]) --> chunk * 3 -- > 8192 * 3
        # sum(diff_count["A"].iloc[0]* lm.params) --> A to C change --> WEIGHTED DIFF
        # print("{} {}".format(diff_count[b1].shape,lm.params.shape))
        del count
        
        sd_diff = dict() 
        t = dict() #for t-test
        # null hypothesis = c'B = 0
        # alt hypo = c'B != 0
        for base in nucleotides: #math.sqrt((diff_count[base].transpose() * np.dot(lm.cov_params(),diff_count[base].transpose())).sum(axis=0))
            sd_diff[base] = (diff_count[base].transpose() * np.dot(lm.cov_params(),diff_count[base].transpose())).sum(axis=0).apply(np.sqrt)
            # standart deviation = T_diff_count * np.dot(covariance,T_diff_count).sum(allrows) and taking square root
            t[base] = diff[base] / sd_diff[base]
        # t["A"][0] --> t value of A to C change     
        # p value --> scipy.stats.norm.sf(abs(t["A"][0]))*2, failing to reject

        dna_seq = [] # to store 12-mer file
        diff_all = []
        t_all = []
        for b1 in nucleotides:
            diff_all += diff[b1].tolist() #conversion into list # all c'B values (3*chunk)*4
            t_all += t[b1].tolist() # all t values (including for each base) , same size with diff_all
            for b2 in nucleotides:
                if b1 != b2:  # | bio.nucleotides[b2] --> place mutated nucleotide to rightmost side
                    dna_seq += [ ((x << 2) | bio.nucleotides[b2]) for x in mutated_part[b1]] # 12-mer (including both mut/ref)
        #[bio.itoseq(ddd) for ddd in dna_seq]            
        newout = pd.DataFrame({'dna_seq':[bio.itoseq(x) for x in dna_seq],
                                        'diff':diff_all,
                                        't':t_all
                                },columns=['dna_seq','diff','t'])
        output_all = output_all.append(newout,ignore_index = True)
        # This will going on until all mutational context in mutated dict is compared 
        # 8192 * 128
    na_entries = pd.DataFrame(columns=['dna_seq','diff','t'])
    
    for base in nucleotides:
        na_entries = na_entries.append(pd.DataFrame({ # sağ tarafa mutated nucleotide i ekleyerek (12-mer format)
            'dna_seq': [bio.itoseq((bio.insert_pos(x,base,k-1) << 2) | bio.nucleotides[base]) for x in bio.seq_permutation(2*(k-1))],
            'diff':np.nan,
            't':np.nan
        },columns=['dna_seq','diff','t']))
    
    output_all = output_all.append(na_entries,ignore_index = True).sort_values(['dna_seq'],ascending=True) #replace(np.nan, 'NaN', regex=True)
    print(output_all[:12])
    #print([scipy.stats.norm.sf(abs(t))*2 for t in t_all]) # 'p':scipy.stats.norm.sf(abs(t_all))*2
    output_all.to_csv("{}/prediction{}mer.{}.txt".format(outpath,args.kmer,filename),columns=['diff','t'],sep=' ',index=None,float_format="%.5f")

    pvals = [scipy.stats.norm.sf(abs(x))*2 for x in output_all['t'].tolist()] # by using t values, calculate two sided p value
    if not os.path.exists("{}/pvals".format(outpath)):
        os.makedirs("{}/pvals".format(outpath))
    with open("{}/pvals/pval{}mer.{}.txt".format(outpath,args.kmer,filename),'w') as f:
        f.write("\n".join('%.4e' % Decimal(p) for p in pvals))

    print("--- Total time: %s seconds ---" % (time.time() - start_time))
#4372.6606secs