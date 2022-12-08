#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:39:41 2022

EXPLAINATION FOR "Bio" FILE

@author: hilmi
"""


import itertools
import pandas as pd
import numpy as np

# all permutations are already reverse-deleted
# all sequences are represented in binary

nucleotides = {'A':0,'C':1,'G':2,'T':3} # numerical conversion (seqtoi dönüşümünde kullanılıo)
#örnegin:
# [bin(0),bin(1),bin(2),bin(3)]
numtonuc = {0:'A',1:'C',2:'G',3:'T'} # string conversion (itoseq dönüşümünde kullanılıyor)

complement = {0:3,3:0,1:2,2:1} # complement for DNA strand # revcomp belirleme ve silmede kullanılıyor

def window(fseq, window_size): 
    for i in range(len(fseq) - window_size + 1):  # L - k + 1 (sum of the counts for any row i)
        yield fseq[i:i+window_size] # her biri belirtilen windows sizedaki kombinasyonları oluşturur ( a generator )

# sample input"":
# count=0
# for x in window("TAACCTCGACAAAATCGTCCTTATAAGACGACCAATGTCTGTGTTCCGTTGTCCGTGCTGT",11):
#     print(x)
#     count += 1

# return the first or the last number representation:
def seqpos(kmer,last): # TRUE = LAST ( AAA'dan TTT'ye kadar olan bütün k-mer permutasyonların sınırı)
    return 1 <<  (1 + 2 * kmer) if last else 1 << 2 * kmer; # left shift (2 right zero)
#sample input for k=6:
# first = seqpos(3,True)
# bin(8192-1)  --> "AAAAAAA" -1 
# itoseq(8191)
# last = seqpos(6,True)
# itoseq(last-1)
# itoseq(64 = seqpos(3,False)) --> "AAA" (1000000) 2*kmer
#itoseq(127 = seqpos(3,True)) --> "TTT" (1111111)

def seq_permutation(seqlen): # seqlen = k = 3 için AAA (64) = 1000000 'den TTT (1111111) range --> BÜTÜN K-MERli featureleri çıkartıyor(revcomp içinde!) 
    return (range(seqpos(seqlen,False),seqpos(seqlen,True))) #position range(4096,8192) for k=6
# seq_permutation(6) # gives a range (START TO END - ALL POSSIBLE elements)

def gen_nonreversed_kmer(k): # it removes comp reverse within the all possible values  --> it creates features according to k value
    nonrevk = list() # unique integer k-mer listemiz
    for i in range(seqpos(k,False),seqpos(k,True)): 
        if i <= revcomp(i): # revcomp değeri kücükse daha önce eklenmiş olduğundan dolayı listeye i değeri eklenmez
            nonrevk.append(i)
    return nonrevk
# gen_nonreversed_kmer(3)
# 75 <= revcomp(75) -->  "AGT" (75)
# 71<= revcomp(71) --> "ACT" (71) 
# k=4
# odd_formula = 4**k/2
# even_formula = (4**k - 2**k)/2 + 2**k --> 4**k/2 + 2**k/2
# len(gen_nonreversed_kmer(6)) == (4**6 - 2**6)/2 + 2**6
# len(gen_nonreversed_kmer(4))
# mylist2 = gen_nonreversed_kmer(3)
# for i in mylist2:
#     print(itoseq(i))
# for k=6 which is even
# revcompstr("ATATAT")  
#for k=5 which is odd
# revcompstr("ATATA")
# revcompstr("GCGCGC")
# revcompstr("GCGCG")
#sonuc olarak tek sayılı k-merlerde eşit sayıda revcomp elementleri varken (4**k/2)
#çift sayılı k-merlerde bazılarının revcomp elementi kendisine eşit olduğundan, +2**k/2 kadar daha çok featureları var
def itoseq(seqint): # seqtoi functionundan çevrilen integerı sekans bilgisine yeniden dönüştürüyor (decode ediyor)
    if type(seqint) is not int:
        return seqint
    seq = ""
    mask = 3
    copy = int(seqint) # prevent changing the original value
    while(copy) != 1: # binrep = 1 olana kadar (seqtoi tersi decoding etmeye devam eder)
        seq = numtonuc[copy&mask] + seq
        copy >>= 2 # right shift (sağa kaydırıp son 2'li tabanda son 2 haneyi silme) 
        if copy == 0:
            print("Could not find the append-left on the input sequence")
            return 0
    return seq # sekansı string olarak veriyor
# örnek:
# itoseq(89)
def seqtoi(seq,gappos=0,gapsize=0): # sekansları integer formlarına cevırıyor
    # due to various seqlengths, this project always needs append 1 to the left
    binrep = 1 
    gaps = range(gappos,gappos+gapsize) # gap kısmını hesaba katıyor  AAA---TAT --> AATTT
    for i in range(0,len(seq)): 
        if i in gaps: 
            continue
        binrep <<= 2 
        binrep |= nucleotides[seq[i]] # hangi base denk geliyorsa onun binarydeki karşılığını alıyor
    return binrep # binary representation 
# myinteger = seqtoi("ATGC")
# bin(myinteger) # 1(start point)-00-11-10

def revcomp(seqbin): # sekans binary representationına bakarak onun reverse complement binary rep çeviriyor
    rev = 1
    mask = 3
    copy = int(seqbin)

    while copy != 1:
        rev <<= 2
        rev |= complement[copy&mask]
        copy >>= 2
        if copy == 0:
            print("Could not find the append-left on the input sequence")
            return 0
    return rev
# revcomp(64) #"AAA" --> "TTT"
# itoseq(127)
# itoseq(revcomp(127))
# revcomp(127) # "TTT" --> "AAA"
itoseq(5876)

def revcompstr(seq): # taking reverse complementary --> string based revcomp dönüşüm 
    rev = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([rev[base] for base in reversed(seq)])

def insert_pos(seqint,base,pos): # pos is position from the right   --> insertion yaparak integer output veriyor
    return ((seqint << 2) & ~(2**(2*pos+2)-1)) | ((seqint & 2**(2*pos)-1) | (nucleotides[base] << pos*2))
    #return (seqint << 2) | (seqint & 2**pos-1) & ~(3 << (pos*2)) | (nucleotides[base] << pos*2)
# insert_pos(4096,"G",0) # "AAA" --> "ATAA" (second position)
# bin(insert_pos(4096,"G",0))
# itoseq(insert_pos(4096,"T",2))

# this function already counts without its reverse complement,
# i.e. oligfreq + reverse merge in the original R code
# Input: panda list and kmer length
# Output: oligonucleotide count with reverse removed
def nonr_olig_freq(seqtbl,kmer,nonrev_list,gappos=0,gapsize=0): # featureların(non_reverse list) her bir observationdaki COUNT TABLE i 
    # with the gapmodel, our model become gapsize + kmer
    gapmer = kmer+gapsize # 3 + 3 
    # separator, since this is binary, the number is counted from the right
    rightseparator = kmer-gappos  # 3-1 = 2  A---,AA right
    leftseparator = rightseparator+gapsize #2+3=5  A,---AA  left
    olig_df =  {k: [0] * len(seqtbl) for k in nonrev_list} # use dictionary first to avoid slow indexing from panda data frame
    for i in range(0,len(seqtbl)): #22s for 3000
        mask = (4**gapmer)-1
        cpy = int(seqtbl[i]) # copy
        while cpy > (4**gapmer)-1:
            # gap calculation here
            cur = cpy & mask
            right = cur & ((4**rightseparator)-1) # right 
            left = (cur >> 2*leftseparator) << 2*rightseparator # left
            gappedseqint = left | right  

            r = (1<<(2*kmer))|gappedseqint # append 1
            rc = revcomp(r)
            if r > rc:
                r = rc
            # 392secs with loc,434 secs with the regression. R time, 10secs for allocation, 3.97mins for linreg
            # with 'at', only 23secs! -- 254secs total for 6mer
            olig_df[r][i] += 1
            cpy >>= 2
    return pd.DataFrame(olig_df)

nonr_olig_freq([seqtoi("TAAAGCCTAATTAGCTCAAATTGAGGACCCCTGCG"),seqtoi("ACGTACTGCTGATGCTTAATTAACGCTTTCATCAC"),seqtoi("ACATGACAGCAGTATGTAATTAATTTCGGAATACA"),seqtoi("TAAAGCCTAATTAGCTCAAATTGAGGACCCCTGCG")], 3, gen_nonreversed_kmer(3))

nonr_olig_freq([seqtoi("TAGCAAATGCATT"),seqtoi("ACGTACTGCTGATGCTTAATTAACGCTTTCATCAC"),seqtoi("ACATGACAGCAGTATGTAATTAATTTCGGAATACA"),seqtoi("TAAAGCCTAATTAGCTCAAATTGAGGACCCCTGCG")], 3, gen_nonreversed_kmer(3), 1, 3)


# gapli kmer-->  pos : the position of the gap within the k-mer.
# ACATGACAGCAGTATGTAATTAATTTCGGAATACA  --> string
#              A---AA      gapmodel = kmer + gapsize = 3+3= 6
#                  A---AA
#                    T---TT
#                     T---TT                                                                                                            
# Totalde 3 gap koydugumuzda 4 counts çıkıyor (nongap versiyonu = 1)  
# itoseq(67) #2 sayım
# itoseq(revcomp(67)) # 2 sayım --> 2 +2 =3 For feature=67("AAT=ATT")
