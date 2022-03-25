#!/usr/bin/env python3
import argparse
import re
#./find_recomb.py -r /srv/scratch/ONT/gridion20220318/test/BC02-022047050501-ONT_recomb.tsv -v /srv/scratch/ONT/gridion20220318/test/BC02-022047050501-ONT_variant.tsv

parser = argparse.ArgumentParser(description='Create a file to summarise recombinant detection')
parser.add_argument('-r', '--recomb')
parser.add_argument('-v', '--variant')

def rename(name_mutation):
    if  re.search("AY",name_mutation) :
        return "Delta"
    elif re.search("B_1_",name_mutation):
        return "Delta"
    else :
        return name_mutation

args = parser.parse_args()
recomb_file=[]
with open(args.recomb,"r") as file:
    for x in file:
        recomb_file.append(x.strip("\n").split("\t"))

variant_file=[]
with open(args.variant,"r") as file:
    for x in file:
        #print(x)
        variant_file.append(x.strip("\n").split("\t"))

variant=""
dico_mutation={}
for element in recomb_file:
    #print(element)
    if len(element) == 1 and element[0] != "":
        variant=element[0]
    elif element[0]=="":
        continue
    else :
        for mutation in element :
            if mutation not in dico_mutation:
                dico_mutation[mutation]=[rename(variant)]
            else :
                if rename(variant) not in dico_mutation[mutation]:
                    dico_mutation[mutation].append(rename(variant))
                    
print(variant_file)
#print(dico_mutation)       

