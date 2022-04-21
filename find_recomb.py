#!/usr/bin/env python3
import argparse
import re
#bga_test="/srv/nfs/ngs-stockage/NGS_Virologie/seqmet/ncov/220325_BH2TLJDRX2_seqmet_varcall_ncov/varcall/bga/MN908947/22Pl259-722001483601_S36.bed"

#./find_recomb.py -r /srv/nfs/ngs-stockage/NGS_Virologie/boyerth/recombinant/22Pl257-722001337801_S274_recomb.tsv -s /srv/nfs/ngs-stockage/NGS_Virologie/boyerth/recombinant/22Pl257-722001337801_S274_variant.tsv
#./find_recomb.py -r /srv/scratch/ONT/gridion20220318/test/BC02-022047050501-ONT_recomb.tsv -s /srv/scratch/ONT/gridion20220318/test/BC02-022047050501-ONT_variant.tsv -d /srv/scratch/ONT/gridion20220318/bam_files/BC02-022047050501-ONT_depth.tsv -o /srv/scratch/ONT/gridion20220318/test/

def rename(name_lineage):
    """
    Cette fonction permet de regrouper les differents noms des linéages
    """
    if  re.search("AY",name_lineage) or re.search("B_1_",name_lineage):
        return "D"
    elif re.search("BA",name_lineage):
        return "O"
    else:
        return "U" #U pour unknown

def read_tsv(tsv_file):
    """
    lis la sortie du script recombV0
    """
    parsed_tsv=[]   
    with open(tsv_file,"r") as file:
        for x in file:
            if x.strip("\n").split("\t") != "":
                parsed_tsv.append(x.strip("\n").split("\t"))
    return parsed_tsv

def dico_mutation(parsed_tsv):
    lineage=""
    dico_out={}
    for element in parsed_tsv :
    #si l'element est un lineage 
        if len(element) == 1: #Possible problème si la liste des mutations n'est que de 1
            lineage=element[0]
        else :
            element.remove('') #enleve le champ vide si il y en a 
            
            for mutation in element :
                pos_mut=mutation.split(',')
                if int(pos_mut[1]) not in dico_out : #si la position n'existe pas on l'ajoute
                    if lineage != "" : #cas avec lineage
                        #creation d'un dico clef = position; valeur = [variant_brut],[liste variants]
                        dico_out[int(pos_mut[1])]=[pos_mut[:],[lineage]] 
                    elif lineage == "" : #cas sans lineage
                        #creation d'un dico clef = position; valeur = [variant_brut][]
                        dico_out[int(pos_mut[1])]=pos_mut[:],[]

                elif int(pos_mut[1]) in dico_out : 
                    if lineage != "" : 
                        if lineage not in dico_out[int(pos_mut[1])][1]: #si le linéage n'est pas déjà dans le dico
                            dico_out[int(pos_mut[1])][1].append(lineage)
    return dico_out

def add_lineage(lst):
    final_call=""
    for lineage in lst :
        if re.search(rename(lineage),final_call) == None:
            final_call+=rename(lineage)
    return final_call

def merge_dico(dico1,dico2):
    merged_dico={}
    print(dico1.keys())
    for x in dico2:
        if x not in dico1.keys():
            merged_dico[x]=dico2[x]
        else :
            merged_dico[x]=dico1[x]
    return merged_dico

def list_bgapos(bga):
    """
    Takes a parsed bga file and output a list with the reference and the depth at every position 
    """
    bgapos = []
    for i in range(len(bga)):
        #verifie que la ref n'existe pas deja, sinon l'ajoute
        if bga[i][0] not in [ x[0] for x in bgapos ]:
            bgapos.append([bga[i][0], []])
            loop_index = [ x[0] for x in bgapos ].index(bga[i][0])
        else:
            loop_index = [ x[0] for x in bgapos ].index(bga[i][0])
        #pour cette ref ajoute la profondeur pour chaque intervalle du fichier bed
        bgapos[loop_index][1].append([int(bga[i][3]) for x in range(int(bga[i][1]),int(bga[i][2]))])
    return bgapos

def matrice_lineage(jesaispas):
    """
    Cette fonction permet de creer un matrice afin de savoir afin de savoir qu'elles linéages sont les plus présents et donc quel profil sont à comparer
    """
    return jesaispas

def amplicon(idontknow):
    return idontknow

parser = argparse.ArgumentParser(description='Create a file to summarise recombinant detection')
parser.add_argument('-r', '--recomb') #le croisement des vcfs au format vcf
parser.add_argument('-s', '--sample') #le tsv des variants de notre echantillon
parser.add_argument('-d', '--depth') #le fichier de profondeur pour chaque position
parser.add_argument('-p', '--phasing') #optionnel : actif si phasing
parser.add_argument('-a', '--amplicon') #optionnel : actif si phasing
parser.add_argument('-o', '--output' ) #le dossier de sortie

args = parser.parse_args()

recomb=read_tsv(args.recomb)
sample=read_tsv(args.sample)

dico_recomb=dico_mutation(recomb) #clef = mutation : val = lst caract mut et lst lineage associées
dico_sample=dico_mutation(sample)

#print(dico_recomb)
#print(dico_sample)

#but calculer le taux d'apparition des linéages, on comparera ceux avec les taux les plus élevé
i=0
lin=""
tuple_freq=[]
del recomb[0]
for x in recomb :
    if i%2 == 0: #élément linéage
        lin=x[0]
    else :
       tuple_freq.append((len(x),lin))
    i+=1

tuple_freq=sorted(tuple_freq, key=lambda tup: tup[0]) #sort le dico pour connaitre les linéages avec le plus de mutations
print(tuple_freq)


#Detection des variations 
#Ajout des attendues dans un dictionnaire 
dico_attendue={}
with open("/srv/nfs/ngs-stockage/NGS_Virologie/boyerth/recombinant/lst_formated.txt",'r') as attendu_file:
    for x in attendu_file :
        x=x.strip("\n").split("\\")
        mut=x[6].split(">")
        dico_attendue[int(x[2])-1]=[x[0],int(x[2])-1,mut[0],mut[1],x[7]],[] #clef = pos : valeur = Ref, pos, nc ref, nc mut, lineage

for x in dico_recomb:
    dico_recomb[x][0].append(add_lineage(dico_recomb[x][1]))
print(dico_recomb)

if args.depth != None:
    dico_depth={}
    with open(args.depth,"r") as depth_file:
        for y in depth_file:
            current=y.strip("\n").split("\t")
            dico_depth[int(current[1])]=int(current[2])

#Liste finale qui sera utilisé pour faire 
lst_final=[]#lineage    position    nc ref  nc mut  %   depth
for clef1 in dico_recomb.keys():
    if clef1 in dico_attendue.keys():              
        lst_final.append([dico_recomb[clef1][0][4],int(dico_recomb[clef1][0][1]),dico_recomb[clef1][0][2],dico_recomb[clef1][0][3],float(dico_sample[clef1][0][4]),dico_depth[clef1]])

for clef2 in dico_attendue.keys():
    if clef2 not in dico_recomb.keys():
        lst_final.append([dico_attendue[clef2][0][4],int(dico_attendue[clef2][0][1]),dico_attendue[clef2][0][2],dico_attendue[clef2][0][3],float(0),dico_depth[clef2]])        

for clef3 in dico_sample.keys():
    if (clef3 not in dico_recomb.keys()) and (clef3 not in dico_attendue.keys()):
        lst_final.append(["U",int(dico_sample[clef3][0][1]),dico_sample[clef3][0][2],dico_sample[clef3][0][3],float(dico_sample[clef3][0][4]),dico_depth[clef3]]) 

lst_final=sorted(lst_final, key=lambda pos:pos[1])

#Recherche d'une/des recombinaison(s)
previous=lst_final[0]
cur_lineage=lst_final[0][0]
toprint=[]
for profil in lst_final:
    #On ne prend pas en compte les elements inconnues pour laquelle la profondeur est inferieure à 20
    if profil[0] != "U" and profil[5]>20:
        #Si ce n'est pas le même lineage on analyse
        if (re.search(profil[0],previous[0]) == None) and (re.search(previous[0],profil[0])==None):
            #Si changement de current lineage pour un variant majoritaire
            if (profil[0] == "D" and profil[4] > 0.5) and cur_lineage == "O" :
                toprint.append(f"Possible recombinaison entre la position {previous[1]} et {profil[1]}\n{previous}\n{profil}\n")

            elif (profil[0] == "O" and profil[4] > 0.5) and cur_lineage == "D":
                toprint.append(f"Possible recombinaison entre la position {previous[1]} et {profil[1]}\n{previous}\n{profil}\n")

            #cas absence mutation
            elif previous[4]==0 and profil[4]==0:
                toprint.append(f"Possible recombinaison entre la position {previous[1]} et {profil[1]} par absence de deux mutations attendues\n{previous}\n{profil}\n")
        previous=profil
        if profil[4]>0.5:
            cur_lineage=profil[0]
    
with open(f"{args.output}_summary_recomb.txt","w") as output_file:
    for x in toprint:
        output_file.write(x)

with open(f"{args.output}_call_recomb.tsv","w") as output_file:
    output_file.write(f"lineage\tposition\tnucleotide ref\tnucleotide_mut\tpourcentage_variation\tprofondeur\n")
    for x in lst_final:
        output_file.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\t{x[5]}\n")