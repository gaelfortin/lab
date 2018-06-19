from os import system
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
amb = IUPACAmbiguousDNA()

print('----------------------------------')
print('Welcome to Cloning 1.0!')
print(' ')
print('SEQUENCES IMPORTATION...')
print('Importing backbones...')

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

bb_dict = {}
#Backbones importation
with open('backbones.txt') as bb:
    nb_bb = int(file_len('backbones.txt')/2)

    for item in range(nb_bb):
       bb_dict[Seq(bb.readline(),amb)] = bb.readline() #Warning! keys and values are imported in the reverse way (unknown reason)

bb_dict = dict(zip(bb_dict.values(), bb_dict.keys()))#putting keys and values in the good order
print('The following backbones will be used for cloning: ',bb_dict.keys())
print(' ')
dna_insert = open('insert.txt','r')
insert = Seq(dna_insert.read(), amb)

#Initialization RestrictionBatch
enzlist = open('enzymes.txt','r')
rb = RestrictionBatch() #RestrictionBatch = Enzymes to scan
nbenz = file_len('ENZYMES.txt')
print('Nb of enzymes:',nbenz)

for loop in range(nbenz):
    enz=enzlist.readline()
    rb.add(enz)

print('Enzymes loading complete')

print('Loaded enzymes:', rb)

#----INPUT OF REGIONS TO CLONE-----#
print("""----------LET'S PREPARE THAT CLONING!----------
Please provide for the following sequences the regions (in bp) to use for cloning.
A help file is available in the application folder for asked positions.""")
#TODO : Do pdf help file to explain regions asked
#TODO: detect errors in given positions (impossible values, 2nd position < 1st position...)
print(' ')
print('---About the insert---')
insert_pos = {}
insert_pos['keep1'] = int(input("Give the starting position of the INSERT (keep1): "))
insert_pos['keep2'] = int(input("Give the ending position of the INSERT (keep2): "))
insert_pos['frame1'] = int(input("Give the beginning of the cloning region of the INSERT (frame1): "))
insert_pos['frame2'] = int(input("Give the end of the cloning region of the INSERT (frame2): "))

print(' ')
print('---About the backbone(s)---')
bb_pos_k1 = {}
bb_pos_k2 = {}
bb_pos_f1 = {}
bb_pos_f2 = {}

#TODO dans l'aide : si on a pas de région à éliminer du backbone, alors mettre mêmes valeurs pour keep et frame NON ne marche pas
#Todo : ask if region to remove, if not simplify code
for elt in bb_dict:
    print(' ')
    print("Let's talk about the backbone ",elt)
    bb_pos_k1[elt] = int(input("Give the starting position for the cloning of the BACKBONE (keep1):"))
    bb_pos_k2[elt] = int(input("Give the ending position for the cloning of the BACKBONE (keep2):"))
    bb_pos_f1[elt] = int(input("Give the beginning of the cloning region for the BACKBONE (frame1):"))
    bb_pos_f2[elt] = int(input("Give the end of the cloning region for the BACKBONE (frame2):"))


print(' ')
review = input('All positions have been saved! Would you like to review them? (Y/N) ')
review = review.upper()
if review == 'Y':
    print("-----Positions of the insert----")
    print('Following names are explained in the help file.')
    print(insert_pos)
    print(' ')
    print("-----Positions of the backbones-----")
    for elt in bb_dict:
        print("Positions of the backbone: ", elt)
        print('Keep1', bb_pos_k1[elt])
        print('Keep2', bb_pos_k2[elt])
        print('Frame1', bb_pos_f1[elt])
        print('Frame2', bb_pos_f2[elt])
        print(' ')
    print('If the above positions are not correct, please restart the software and enter correct positions.')
    system("pause")

print(' ')
print("""----------LET'S BEGIN!----------""")
#Defining regions usable for cloning:
#insert_r1 = insert[insert_pos['frame1']:insert_pos['keep1']]
#insert_r2 = insert[insert_pos['keep2']:insert_pos['frame2']]
#IS ALL THIS USEFUL ??? WE SCAN THE WHOLE PLASMID FOR THE ENZYMES BECAUSE THEY MUST NOT CUT ELSEWHERE!
#bb_r1 = {}
#bb_r2 = {}

#for elt in bb_dict:
    #bb_r1[elt] = bb_dict[elt][bb_pos_f1[elt]: bb_pos_k1[elt]]
    #bb_r2[elt] = bb_dict[elt][bb_pos_k2[elt] : bb_pos_f2[elt]]

#Enzymes analysis
print('Insert analysis:')
Ana_insert = Analysis(rb, insert, linear= False) #rb enzymes
Ana_insert.print_as('number')
Ana_insert.print_that()

print(' ')
print('Backbones analysis:')
Ana_bb = {}
for elt in bb_dict: #backbones analysis
    print('Analysis of ', elt)
    Ana_bb[elt] = Analysis(rb, bb_dict[elt], linear=False)
    Ana_bb[elt].print_as('number')
    Ana_bb[elt].print_that()

#Enzymes cutting once in INSERT and only in R1 or R2
insert_enz_once = Ana_insert.with_N_sites(1)
insert_enz_r1 = {}
insert_enz_r2 = {}
for elt in insert_enz_once:
    if insert_pos['frame1'] <= int(''.join(map(str, insert_enz_once[elt]))) <= insert_pos['keep1']:
        insert_enz_r1[elt] = int(''.join(map(str, insert_enz_once[elt])))
    elif insert_pos['keep2'] <= int(''.join(map(str, insert_enz_once[elt]))) <= insert_pos['frame2']:
        insert_enz_r2[elt] = int(''.join(map(str, insert_enz_once[elt])))

print('Enzymes cutting once in the plasmid and located in R1: ', insert_enz_r1) #TODO: organize the 2 dicts to give enzymes in growing positions
print('Enzymes cutting once in the plasmid and located in R2: ', insert_enz_r2)

#Enzymes cutting once in BACKBONES in R1 and R2
bb_enz_once = {} #dict that will stocks the Ana.with_N_sites(1) for each backbone
bb_enz_r1 = {}
bb_enz_r2 = {}


for elt in bb_dict:
    bb_enz_once[elt] = Ana_bb[elt].with_N_sites(1)

for bb_n in bb_dict: #loop that extract enzymes cutting once in the BACKBONE and in the correct region
    bb_enz_r1_temp = [] #temp LISTS storing enzymes for 1 backbone only
    bb_enz_r2_temp = []
    for elt in bb_enz_once[bb_n]:
        if bb_pos_f1[bb_n] <= int(''.join(map(str, bb_enz_once[bb_n][elt]))) <= bb_pos_k1[bb_n]: #WARNING! if no region to keep, enzyme will be in R1 AND R2
            bb_enz_r1_temp.append(elt)
        if bb_pos_k2[bb_n] <= int(''.join(map(str, bb_enz_once[bb_n][elt]))) <= bb_pos_f2[bb_n]:
            bb_enz_r2_temp.append(elt)
    bb_enz_r1[bb_n] = bb_enz_r1_temp #temp lists are stored in the final dict with item = [backbone : list of enzymes in the region]
    bb_enz_r2[bb_n] = bb_enz_r2_temp



#-----Enzymes comparison-----#
print(' ')
print('Comparing enzymes...')
print(' ')

insert_r1_set = set(insert_enz_r1) #creating sets used for comparing insert and backbones enzymes
insert_r2_set = set(insert_enz_r2)
insert_r1_list = list(insert_enz_r1)
insert_r2_list = list(insert_enz_r2)

for bb_n in bb_dict:
    print('------Cloning with the backbone ',bb_n, '-----')
    print(' ')
    print('----Enzymes for the region R1----')
    #exact enzymes for R1
    enz_comp_r1 = set.intersection(insert_r1_set, set(bb_enz_r1[bb_n]))
    print('Exact enzymes for R1 (before the insert): ', enz_comp_r1)
    print(' ')
    print('Enzymes in R1 compatible between the backbone and the insert:')
    print('Enzyme in the backbone | Enzyme in the insert')
    for bb_enz in bb_enz_r1[bb_n]:    #compatible enzymes for R1
        for ins_enz in insert_r1_list:
            if ins_enz%bb_enz == True and ins_enz != bb_enz:
                print(bb_enz, '   |   ', ins_enz)
    #---------------------------------------------------
    print("""
    ----Enzymes for the region R2----""")
    # exact enzymes for R2
    enz_comp_r2 = set.intersection(insert_r2_set, set(bb_enz_r2[bb_n]))
    print('Exact enzymes for R2 (after the insert): ', enz_comp_r2)
    print(' ')
    print('Enzymes in R2 compatible between the backbone and the insert:')
    print('Enzyme in the backbone | Enzyme in the insert')
    for bb_enz in bb_enz_r2[bb_n]:  # compatible enzymes for R2
        for ins_enz in insert_r2_list:
            if ins_enz % bb_enz == True and ins_enz != bb_enz:
                print(bb_enz, '   |   ', ins_enz)
print("""
--------------------------------------
Done! 
All combinations have been tested but don't forget that ORFs and fragments orders have not been verified.

By Gael Fortin. June 2018.""")


#When looking for enzymes, be sure that 1st enzyme position<2nd enzyme position in BACKBONE to avoid inverted clonings in case of no seq to keep)

#TODO: check for ORFs
#TODO: export DNA final sequence
#TODO: cloning with 3 enzymes
#TODO: 2-step cloning

dna_insert.close()

system('pause')