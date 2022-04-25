#!/usr/bin/env python3
import argparse
import re
#bga_test="/srv/nfs/ngs-stockage/NGS_Virologie/seqmet/ncov/220325_BH2TLJDRX2_seqmet_varcall_ncov/varcall/bga/MN908947/22Pl259-722001483601_S36.bed"

#./find_recomb.py -r /srv/nfs/ngs-stockage/NGS_Virologie/boyerth/recombinant/22Pl257-722001337801_S274_recomb.tsv -s /srv/nfs/ngs-stockage/NGS_Virologie/boyerth/recombinant/22Pl257-722001337801_S274_variant.tsv
#./find_recomb.py -r /srv/scratch/ONT/gridion20220318/test/BC02-022047050501-ONT_recomb.tsv -s /srv/scratch/ONT/gridion20220318/test/BC02-022047050501-ONT_variant.tsv -d /srv/scratch/ONT/gridion20220318/bam_files/BC02-022047050501-ONT_depth.tsv -o /srv/scratch/ONT/gridion20220318/test/
lst_delta=['AY_1','AY_10','AY_100','AY_101','AY_102','AY_102_1','AY_102_2','AY_103','AY_104','AY_105','AY_106','AY_107','AY_108','AY_109','AY_11','AY_110','AY_111','AY_112','AY_112_1','AY_113','AY_114','AY_116','AY_116_1','AY_117','AY_118','AY_119','AY_119_1','AY_119_2','AY_120','AY_120_1','AY_120_2','AY_120_2_1','AY_121','AY_121_1','AY_122','AY_122_1','AY_122_2','AY_122_3','AY_123','AY_123_1','AY_124','AY_124_1','AY_125','AY_126','AY_127','AY_127_1','AY_128','AY_129','AY_13','AY_131','AY_132','AY_133','AY_14','AY_15','AY_16','AY_16_1','AY_17','AY_18','AY_19','AY_2','AY_20','AY_20_1','AY_21','AY_22','AY_23','AY_23_1','AY_23_2','AY_24','AY_25','AY_25_1','AY_25_1_1','AY_26','AY_26_1','AY_27','AY_28','AY_29','AY_29_1','AY_3','AY_30','AY_31','AY_32','AY_33','AY_33_1','AY_34','AY_34_1','AY_34_1_1','AY_34_2','AY_35','AY_36','AY_37','AY_38','AY_39','AY_39_1','AY_39_1_1','AY_39_1_2','AY_39_1_3','AY_39_2','AY_3_1','AY_3_2','AY_3_3','AY_4','AY_40','AY_41','AY_42','AY_42_1','AY_43','AY_43_1','AY_43_2','AY_43_3','AY_43_4','AY_43_5','AY_43_6','AY_43_7','AY_44','AY_45','AY_46','AY_46_1','AY_46_2','AY_46_3','AY_46_4','AY_46_5','AY_46_6','AY_46_6_1','AY_47','AY_48','AY_49','AY_4_1','AY_4_10','AY_4_2','AY_4_2_1','AY_4_2_2','AY_4_2_3','AY_4_3','AY_4_4','AY_4_5','AY_4_6','AY_4_7','AY_4_8','AY_4_9','AY_5','AY_50','AY_51','AY_52','AY_53','AY_54','AY_55','AY_56','AY_57','AY_58','AY_59','AY_5_1','AY_5_2','AY_5_3','AY_5_4','AY_5_5','AY_6','AY_60','AY_61','AY_62','AY_63','AY_64','AY_65','AY_66','AY_67','AY_68','AY_69','AY_7','AY_70','AY_71','AY_72','AY_73','AY_74','AY_75','AY_75_2','AY_75_3','AY_76','AY_77','AY_78','AY_79','AY_7_1','AY_7_2','AY_8','AY_80','AY_81','AY_82','AY_83','AY_84','AY_85','AY_86','AY_87','AY_88','AY_9','AY_90','AY_91','AY_91_1','AY_92','AY_93','AY_94','AY_95','AY_96','AY_98','AY_98_1','AY_99','AY_99_1','AY_99_2','AY_9_2','AY_9_2_1','AY_9_2_2','B_1_617_2','B_1_640_1']
lst_omicron=['BA_1','BA_1_1','BA_2','BA_3']

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
    lis la sortie du script recombV0 et retourne une liste par ligne
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

def searchvcf(vcf_name):
    myfile=f"/srv/scratch/iai/seqmet/db/ref/ncov/MN908947/vcf/{vcf_name}.vcf"
    match=[]
    with open(myfile,"r") as vcf:
        lines=vcf.readlines()
        for element in lines :
            test=element.rstrip("\n").split("\t")
            if re.search("^(?!#)",test[0]): #si pas un hearder du vcf
                match.append(f"{test[0]},{int(test[1])-1},{test[3]},{test[4]}") #test[1]-1 car il semblerait que le script précédent mette -1 à toutes les positions
    return match

def dico_attendu(dico_sample,mutation_monitored,mutation_vcf,lst_lin):
    for mut in set(mutation_vcf): #vérification des mutations présentent dans les vcf et absentes des mutations de l'échantillons
        if mut not in set(mutation_monitored):
            mut=mut.split(",")
            if int(mut[1]) not in dico_sample.keys() :
                dico_sample[int(mut[1])]=[[mut[0],mut[1],mut[2],mut[3],0],[]]
                for var in lst_lin :
                    dico_sample[int(mut[1])][1].append(var)
            else :
                for var in lst_lin :
                    dico_sample[int(mut[1])][1].append(var)
        else :
            mut=mut.split(",")
            for var in lst_lin :
                dico_sample[int(mut[1])][1].append(var)
    return dico_sample

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

###but calculer le taux d'apparition des linéages, on comparera ceux avec les taux les plus élevé
i=0
lin=""
dico_lin={}
tuple_freq=[]
max_lin_D=[]
max_lin_O=[]
del recomb[0]
for x in recomb :
    if i%2 == 0: #élément linéage
        lin=x[0] #On garde le lineage en mémoire
    else :
        tuple_freq.append((len(x),lin,x))
        if lin in lst_delta: #séparation des éléments omicrons et delta
            max_lin_D.append(len(x))
            dico_lin[lin]=x
        elif lin in lst_omicron:
            max_lin_O.append(len(x))
            dico_lin[lin]=x
    i+=1

###

###recupération des vcf dont les mutations delta et omicron sont les plus présentes 
### Rq: pas de prise en charge de la comparaison Omicron-Omicron --> Z:\redmine\14181\run_220404 = liste echantillon recombinant pour tests
mutation_monitored=[] #mutations retrouvés dans échantillons
vcf_x=[]
vcf_y=[]
lst_vcf_1=[]
lst_vcf_2=[]


#Situation Delta - Omicron
for lineage in tuple_freq:
    if lineage[0] == max(max_lin_D) and lineage[1] in lst_delta: #cas des mutations delta --> récupération des linéages delta les plus représentatifs
        lst_vcf_1.append(lineage[1]) #ajout du linéage dans une liste --> à utiliser en output
        vcf_x+=searchvcf(lineage[1]) #ensemble des mutations présentes dans les différents linéages majoritaires
        for element in lineage[2]:
            mutation_monitored.append(element)
    elif lineage[0] == max(max_lin_O) and lineage[1] in lst_omicron: #cas des mutations omicrons --> récupération des linéages omicron les plus représentatifs
        lst_vcf_2.append(lineage[1]) #ajout du linéage dans une liste --> à utiliser en output
        vcf_y+=searchvcf(lineage[1]) #ensemble des mutations présentes dans les différents linéages majoritaires
        for element in lineage[2]:
            mutation_monitored.append(element)
###

###
### Prise en compte des recombs Omicron - Omicron --> Faire une matrice avec toutes les linéages Omicrons
### Comparer les mutations retrouvé dans l'échantillon avec chaques linéages, si un linéage a une mutation unique (genre BA.2 ou BA.3 alors ya moyen qu'il y ai plusieurs recombinants)
### 

### creation liste mutations présentes et absentes
### permet d'attribuer un pourcentage à 0

#ajout des linéages attendues et détectés
dico_x=dico_attendu(dico_sample,mutation_monitored,vcf_x,lst_vcf_1)
dico_attendu=dico_attendu(dico_x,mutation_monitored,vcf_y,lst_vcf_2)


# Detection des variations
# Ajout des attendues dans un dictionnaire
# dico_attendue2={}
# with open("/srv/nfs/ngs-stockage/NGS_Virologie/boyerth/recombinant/lst_formated.txt",'r') as attendu_file: #Creer un fichier avant avec les attendues défini par les vcf correspondants
#     for x in attendu_file :
#         x=x.strip("\n").split("\\")
#         mut=x[6].split(">")
#         dico_attendue2[int(x[2])-1]=[x[0],int(x[2])-1,mut[0],mut[1],x[7]],[] #clef = pos : valeur = Ref, pos, nc ref, nc mut, lineage

# for x in dico_recomb:
#     dico_recomb[x][0].append(add_lineage(dico_recomb[x][1]))


###ajout de la profondeur
if args.depth != None:
    dico_depth={}
    with open(args.depth,"r") as depth_file:
        for y in depth_file:
            current=y.strip("\n").split("\t")
            dico_depth[int(current[1])]=int(current[2])

for clef in dico_attendu.keys():
    dico_attendu[clef][0].append(dico_depth[clef-1])

def diminutif(lst_vcf_1,lst_vcf_2,current_lst):
    """
    cette fonction prend en entrée deux listes de linéage et une liste de linéage et retourne un diminutif en fonction de la présence ou non de ces linéages
    """
    diminutif=[]
    if len(current_lst)==0:
        diminutif.append("0")
    else :
        for lineage in current_lst :
            if lineage in lst_vcf_1 :
                diminutif.append("1")
            elif lineage in lst_vcf_2 :
                diminutif.append("2")
    diminutif=set(diminutif)
    return diminutif

#recherche de zones de cassures
cur_lineage=["0"]
for mutation in dico_attendu.values():
    profil=diminutif(lst_vcf_1,lst_vcf_2,mutation[1])
    if int(mutation[0][5])>=20 and "0" not in profil:
        if "1" in profil and "2" in profil: #cas {1,2} --> Changement de linéage impossible on continue avec une nouvelle boucle (évite aussi de démarrer avec un profil{1,2})
            continue

        elif "0" in cur_lineage and ("1" in profil or "2" in profil): #début du script, détection du premier linéage
            cur_lineage=profil
            previous=mutation[0][1]

        elif "1" in cur_lineage and ("2" in profil and float(mutation[0][4])>0.5): #cas {2} sur un profil {1} possible cassure
            print(f"ancien profil : {lst_vcf_1}\n{previous} --> {mutation[0][1]}\nnouveau profil : {lst_vcf_2}")
            previous=mutation[0][1]
            cur_lineage=profil
            

        elif "2" in cur_lineage and ("1" in profil and float(mutation[0][4])>0.5): #cas {1} sur un profil {2} possible cassure
            print(f"ancien profil : {lst_vcf_2}\n{previous} --> {mutation[0][1]}\nnouveau profil : {lst_vcf_1}")
            previous=mutation[0][1]
            cur_lineage=profil



#TROUVER LE DERNIER Delta ou Omicron pour complèter les zones
#Sauvegarder les zones pour rechercher des absences d'Omicron / Delta en même temps 

#Recherche d'une/des recombinaison(s)
# previous=""
# cur_lineage=""
# toprint=[]
# for profil in lst_final:
#     #On ne prend pas en compte les elements inconnues pour laquelle la profondeur est inferieure à 20
#     if profil[0] != "U" and profil[5]>20:
#         #Si ce n'est pas le même lineage on analyse
#         if (re.search(profil[0],previous[0]) == None) and (re.search(previous[0],profil[0])==None):
#             #Si changement de current lineage pour un variant majoritaire
#             if (profil[0] == "D" and profil[4] > 0.5) and cur_lineage == "O" :
#                 toprint.append(f"Possible recombinaison entre la position {previous[1]} et {profil[1]}\n{previous}\n{profil}\n")

#             elif (profil[0] == "O" and profil[4] > 0.5) and cur_lineage == "D":
#                 toprint.append(f"Possible recombinaison entre la position {previous[1]} et {profil[1]}\n{previous}\n{profil}\n")

#             #cas absence mutation
#             elif previous[4]==0 and profil[4]==0:
#                 toprint.append(f"Possible recombinaison entre la position {previous[1]} et {profil[1]} par absence de deux mutations attendues\n{previous}\n{profil}\n")
#         previous=profil
#         if profil[4]>0.5:
#             cur_lineage=profil[0]


# #Liste finale qui sera utilisé pour faire
# lst_final=[]#lineage    position    nc ref  nc mut  %   depth
# for clef1 in dico_recomb.keys():
#     if clef1 in dico_attendue.keys():
#         lst_final.append([dico_recomb[clef1][0][4],int(dico_recomb[clef1][0][1]),dico_recomb[clef1][0][2],dico_recomb[clef1][0][3],float(dico_sample[clef1][0][4]),dico_depth[clef1]])

# for clef2 in dico_attendue.keys():
#     if clef2 not in dico_recomb.keys():
#         lst_final.append([dico_attendue[clef2][0][4],int(dico_attendue[clef2][0][1]),dico_attendue[clef2][0][2],dico_attendue[clef2][0][3],float(0),dico_depth[clef2]])

# for clef3 in dico_sample.keys():
#     if (clef3 not in dico_recomb.keys()) and (clef3 not in dico_attendue.keys()):
#         lst_final.append(["U",int(dico_sample[clef3][0][1]),dico_sample[clef3][0][2],dico_sample[clef3][0][3],float(dico_sample[clef3][0][4]),dico_depth[clef3]])

# lst_final=sorted(lst_final, key=lambda pos:pos[1])

#Recherche d'une/des recombinaison(s)
# previous=""
# cur_lineage=""
# toprint=[]
# for profil in lst_final:
#     #On ne prend pas en compte les elements inconnues pour laquelle la profondeur est inferieure à 20
#     if profil[0] != "U" and profil[5]>20:
#         #Si ce n'est pas le même lineage on analyse
#         if (re.search(profil[0],previous[0]) == None) and (re.search(previous[0],profil[0])==None):
#             #Si changement de current lineage pour un variant majoritaire
#             if (profil[0] == "D" and profil[4] > 0.5) and cur_lineage == "O" :
#                 toprint.append(f"Possible recombinaison entre la position {previous[1]} et {profil[1]}\n{previous}\n{profil}\n")

#             elif (profil[0] == "O" and profil[4] > 0.5) and cur_lineage == "D":
#                 toprint.append(f"Possible recombinaison entre la position {previous[1]} et {profil[1]}\n{previous}\n{profil}\n")

#             #cas absence mutation
#             elif previous[4]==0 and profil[4]==0:
#                 toprint.append(f"Possible recombinaison entre la position {previous[1]} et {profil[1]} par absence de deux mutations attendues\n{previous}\n{profil}\n")
#         previous=profil
#         if profil[4]>0.5:
#             cur_lineage=profil[0]

# with open(f"{args.output}_summary_recomb.txt","w") as output_file:
#     for x in toprint:
#         output_file.write(x)

# with open(f"{args.output}_call_recomb.tsv","w") as output_file:
#     output_file.write(f"lineage\tposition\tnucleotide ref\tnucleotide_mut\tpourcentage_variation\tprofondeur\n")
#     for x in lst_final:
#         output_file.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\t{x[5]}\n")