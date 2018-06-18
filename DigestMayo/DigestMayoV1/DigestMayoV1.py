# -*- coding :Latin -1 -*
#INITIALIZATION
#TODO : reduce imports for size of standalone exe
import os
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
amb = IUPACAmbiguousDNA()
from bandwagon import (BandsPattern, BandsPatternsSet, LADDER_100_to_4k,
                       compute_digestion_bands)

#-----IMPORTATION OF SEQUENCES 1 AND 2-----#
dna_vec = open('vector.txt','r') #Vector importation
org_vec = dna_vec.read()
vec = Seq(org_vec,amb)

dna_constr = open('construction.txt','r') #Construction importation
org_constr = dna_constr.read()
constr = Seq(org_constr,amb)

print('Sequences importated!')

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

nbenz = file_len('enzymes.txt')
print('Nb of enzymes:',nbenz)

#Initialization RestrictionBatch
enzlist = open('enzymes.txt','r')
rb = RestrictionBatch() #RestrictionBatch = Enzymes to scan

for loop in range(nbenz):
    enz=enzlist.readline()
    #print(enz)
    rb.add(enz)
    #print('Contenu de RB:',rb)

print('Enzymes loading complete')

print(rb)

#-----COMPLETE ANALYSIS OF SEQUENCES 1 AND 2-----#
print('Vector analysis:')
Ana_vec = Analysis(rb, vec, linear=False) #Analysis of vector to find rb enzymes
Ana_vec.print_as('number')
Ana_vec.print_that()

print('Construction analysis:')
Ana_constr = Analysis(rb, constr, linear=False) #Analysis of my_seq to find rb enzymes
Ana_constr.print_as('number')
Ana_constr.print_that()


unisites_vec = Ana_vec.with_N_sites(1)
print('Enzymes that cut only once in Vector:',unisites_vec)

unisites_constr = Ana_constr.with_N_sites(1)
print('Enzymes that cut only once in Construction:',unisites_constr)

#-----ENZYMES COMPARISON-----#
#Finding enzymes cutting once in BOTH sequences
univec_set = set(unisites_vec.keys())
uniconstr_set = set(unisites_constr.keys())
comenz = set.intersection(univec_set,uniconstr_set)

print('Enzymes cutting once in each sequence:',comenz)

#Finding enzymes corresponding positions
print('Sites position of enzymes in VECTOR:')
for elt in comenz:
    print(elt, unisites_vec[elt])

print('Sites position of enzymes in CONSTRUCTION:')
for elt in comenz:
    print(elt, unisites_constr[elt])

print(' ')

#Finding enzymes cutting in VECTOR only
nocut_constr = Ana_constr.without_site()
nocutconstr_set = set(nocut_constr)
only_vec = set.intersection(univec_set, nocutconstr_set)
print('Enzymes cutting 1x in VECTOR only:', only_vec)
print('Sites position of enzymes in VECTOR:')
for elt in only_vec:
    print(elt, unisites_vec[elt])

#Find enzymes cutting in CONSTRUCTION only
nocut_vec = Ana_vec.without_site()
nocutvec_set = set(nocut_vec)
only_constr = set.intersection(uniconstr_set, nocutvec_set)
print('Enzymes cutting 1x in CONSTRUCTION only:', only_constr) #TODO : do a IF in case 0 enzymes are in this set (instead of printing 'set()'
print('Sites position of enzymes in CONSTRUCTION only:')
for elt in only_vec:
    print(elt, unisites_vec[elt])

print(' ')
ddigquestion = input('Would you like to do a double digestion? (Y/N)')
ddigquestion = ddigquestion.upper()
if ddigquestion == 'Y':
    # -----DOUBLE DIGESTION-----#
    print(' ')
    print('----------DOUBLE DIGESTION----------')
    print(' ')
    for elt_com in comenz:
        for elt_vec in only_vec:
            # Vector digestion
            poscom_vec = int(''.join(map(str, unisites_vec[elt_com])))
            posonly_vec = int(''.join(map(str, unisites_vec[elt_vec])))
            band1 = abs(posonly_vec - poscom_vec)  # Fragment size calculations
            band2 = len(vec) - band1

            if (abs(band1 - band2) < 700) or ((band1 + band2) > 2000 and abs(band1 - band2) < 1000) or (
                    band1 < 800 or band2 < 800):  # Filter for observable bands on gel
                continue
            # Construction digestion
            band3 = len(vec)
            print(' ')
            print('Digestion with ', elt_com, ' + ', elt_vec, ':')
            print('Vector only: ', band1, ' + ', band2)
            print('Construction: ', band3)

        for elt_constr in only_constr:
            # Construction digestion
            poscom_constr = int(''.join(map(str, unisites_constr[elt_com])))  # position of enzyme
            posonly_constr = int(''.join(map(str, unisites_constr[elt_constr])))
            band1 = abs(posonly_constr - poscom_constr)
            band2 = len(vec) - band1

            if (abs(band1 - band2) < 700) or ((band1 + band2) > 2000 and abs(band1 - band2) < 1000) or (
                    band1 < 800 or band2 < 800):  # Filter for observable bands on gel
                continue

            # Vector digestion: only 1 enzyme = 1 fragment
            band3 = len(vec)
            print(' ')
            print('Digestion with ', elt_com, ' + ', elt_constr, ':')
            print('Vector only: ', band3)
            print('Construction: ', band1, ' + ', band2)




print(' ')
print('PROGRAM COMPLETED!')


#postcom_vec_list = [4560]
#for eltx in postcom_vec_list:
#    print(poscom_vec_list)




#TODO: predict gels and propose best digestions

