#!/usr/bin/env python3
import argparse
from pathlib import Path, PurePath
from pickle import MARK
import re
##Count the number of minor variants in a target vcf reported as major variant in a reference vcf
#v0.0.6

def list_intersect(containedl, containingl, mode, threshold):
    """
    return what is missing or matching between two list of mutation
    """
    if mode == 'match':
        return [ x for x in containedl if x in containingl ]
    elif mode == 'missing':
        return [ x for x in containedl if x not in containingl ]

def count_commented(file):
    """
    count line with a "#"
    """
    lines = open(file, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')
    count = 0
    for line in lines:
        if line[0] == "#":
            count += 1
    return count

def list_flatten(toflat_list):
    """
    Recreate a flat list with a list of list 
    """
    toflat_list = [item for ilist in toflat_list for item in ilist]
    return toflat_list

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

def list_bedpos(bed):
    """
    Takes a parsed bed file and output a list with the reference and every positions for each amplicon
    """
    bedpos = []
    for i in range(len(bed)):
        #verifie que la ref n'existe pas deja, sinon l'ajoute
        if bed[i][0] not in [ x[0] for x in bedpos ]:
            bedpos.append([bed[i][0], []])
            loop_index = [ x[0] for x in bedpos ].index(bed[i][0])
        else:
            loop_index = [ x[0] for x in bedpos ].index(bed[i][0])
        #pour cette ref ajoute les positions entre chaque intervalle du fichier bed
        bedpos[loop_index][1].append([int(x) for x in range(int(bed[i][1]),int(bed[i][2]))])
    return bedpos

def vcf_parse(vcf, fields):
    """
    takes a parsed VCF file and output for every line in a list  chrom, pos, ref, alt  + specific(s) field(s) , dp
    """
    if len(vcf) == 0:
        return []
    fields_header = [ x.split('=')[0] for x in vcf[0].split('\t')[7].split(';') if '=' in x ]
    fields_index = [ fields_header.index(x) for x in fields.split(',') ]
    var = [ [x.split('\t')[0], int(x.split('\t')[1])-1, x.split('\t')[3], x.split('\t')[4]] for x in vcf]
    for i in range(len(var)):
        for j in range(len(fields_index)):
            var[i].append(vcf[i].split('\t')[7].split(';')[fields_index[j]].split('=')[1])
        if var[i][3][0] == '-':
            var_temp = var[i][2]
            var[i][2] = var[i][2] + var[i][3][1:]
            var[i][3] = var_temp
        if var[i][3][0] == '+':
            var[i][3] = var[i][2] + var[i][3][1:]
    return var

parser = argparse.ArgumentParser(description='Count the number of minor variants in a target vcf reported as major variant in a reference vcf')
debugmode = parser.add_mutually_exclusive_group()
debugmode.add_argument('-v', '--verbose', action='store_true')
debugmode.add_argument('-q', '--quiet', action='store_true')
parser.add_argument('--version', action='version', version='0.0.6')
parser.add_argument('-t', '--target', help='the base')
parser.add_argument('-r', '--reference', help='the base')
parser.add_argument('-e', '--exclusion', help='the base')
parser.add_argument('-m', '--mode', help='the base')
parser.add_argument('-d', '--depth', help='the base')
parser.add_argument('-b', '--bed', help='the base')
parser.add_argument('--min_depth', type=int, default=100, help='the base')
parser.add_argument('--min_freq', type=float, default=0.05, help='the base')
parser.add_argument('-o', '--output', help='output file', default='./')

if __name__ == '__main__':
    args = parser.parse_args()

#seq = [[x.split('\n')[0], ''.join(x.split('\n')[1:]).replace(' ','')] for x in open(args.reference, 'r').read().replace('\r\n','\n').rstrip('\n').split('>')[1:]]
tvcf = open(args.target, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.target):] #liste de chaque ligne
rvcf = open(args.reference, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.reference):]
bga = [x.split('\t') for x in open(args.depth, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')]
bed = [x.split('\t') for x in open(args.bed, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')]

depth = [ [x[0], list_flatten(x[1])] for x in list_bgapos(bga) ] #creation d'une liste [ ref, [ toute les profondeurs pour chaque pos]]
region = [ [x[0], list_flatten(x[1])] for x in list_bedpos(bed) ] #creation d'une liste [ref, [ pos amplicons ] ]

#chrom, pos, ref, alt, af, dp
tvar = vcf_parse(tvcf, 'AF') #recuperation des champs [ chrom, pos, ref, alt, af, dp ]

rvar = vcf_parse(rvcf, 'AF') #recuperation des champs [ chrom, pos, ref, alt, af, dp ]

depth_chrom = [y[0] for y in depth] #recuperation du "chromosome/ref"
region_chrom = [y[0] for y in region] #recuperation du "chromosome/ref"

#recupere les variants majoritaires avec une profondeur suffisante et dont la position correspond à un amplicon
#ref,pos,nc_orig,nc_mut,freq
tvar_major = [ x for x in tvar if float(x[4])>=0.5 and depth[depth_chrom.index(x[0])][1][x[1]]>args.min_depth and x[1] in region[region_chrom.index(x[0])][1] ]
#recupere les variants minoritaire avec une profondeur et une frequence suffisante et dont la position correspond à un amplicon
tvar_minor = [ x for x in tvar if float(x[4])<0.5 and float(x[4])>args.min_freq and depth[depth_chrom.index(x[0])][1][x[1]]>args.min_depth and x[1] in region[region_chrom.index(x[0])][1] ]
#recupere tous les variants majoritaires ayant une profondeur suffisante
rvar_major = [ x for x in rvar if float(x[4])>=0.5 and depth[depth_chrom.index(x[0])][1][x[1]]>args.min_depth ]
#rvar_minor = [ x for x in rvar if float(x[4])<0.5 and depth[depth_chrom.index(x[0])][1][x[1]]>args.min_depth ]


#même chose que tvar mais sans la frequence
tvar_major_profile = [ [ x[0], x[1], x[2], x[3] ] for x in tvar_major ] 
#permet d'obtenir le profil reverse pour les snps (pas les insersion/deletions) (inversion des champs 2 et 3)
tvar_major_profile_reversed = [ [ x[0], x[1], x[3], x[2] ] for x in tvar_major if len(x[2]) == 1 and len(x[3]) == 1]
tvar_minor_profile = [ [ x[0], x[1], x[2], x[3] ] for x in tvar_minor ]
rvar_major_profile = [ [ x[0], x[1], x[2], x[3] ] for x in rvar_major ]
rvar_major_profile_reversed = [ [ x[0], x[1], x[3], x[2] ] for x in rvar_major if len(x[2]) == 1 and len(x[3]) == 1]

expected = []
common = []

if args.mode == 'raw':
    expected += list_intersect( rvar_major_profile, tvar_major_profile, 'missing', 100 )
    expected += list_intersect( tvar_major_profile_reversed, rvar_major_profile_reversed , 'missing', 100 )
if args.mode == 'min':
    expected += list_intersect( rvar_major_profile, tvar_major_profile, 'missing', 100 )
    #expected += list_intersect( tvar_major_profile_reversed, rvar_major_profile_reversed , 'missing', 100 )
elif args.mode == 'maj':
    expected += rvar_major_profile
    expected += list_intersect( tvar_major_profile_reversed, rvar_major_profile_reversed , 'missing', 100 )

if args.exclusion:
    evcf = open(args.exclusion, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.exclusion):]
    evar = vcf_parse(evcf, 'AF')
    evar_major = [ x for x in evar if float(x[4])>=0.5 ]
    evar_major_profile = [ [ x[0], x[1], x[2], x[3] ] for x in evar_major ]
    evar_major_profile_reversed = [ [ x[0], x[1], x[3], x[2] ] for x in evar_major if len(x[2]) == 1 and len(x[3]) == 1]
    expected = list_intersect( expected, evar_major_profile , 'missing', 100 )
    expected = list_intersect( expected, evar_major_profile_reversed , 'missing', 100 )

expected = [list(x) for x in set(tuple(x) for x in expected)]
expected.sort()

if args.mode == 'min' or args.mode == 'raw':
    common += list_intersect( tvar_minor_profile , expected, 'match', 100 )
elif args.mode == 'maj':
    common += list_intersect( tvar_minor_profile , expected, 'match', 100 )
    common += list_intersect( tvar_major_profile , expected, 'match', 100 )


if args.mode == 'recomb':
    recomb=[]
    #ref,pos,nc_orig,nc_mut,freq
    #sort les variants communs entre la target et la ref
    filename=(Path(args.reference).stem)
    recomb.append(filename)
    recomb += list_intersect(tvar_major_profile,rvar_major_profile,'match', 100)
    #on compare le profil mineur avec le profil majeur des variant
    recomb += list_intersect(tvar_minor_profile,rvar_major_profile,'match', 100)

common = [list(x) for x in set(tuple(x) for x in common)]
common.sort()

common_length = len(common)
expected_length = len(expected)

#print(expected)
#print(common)
#print(expected_length)
#print(common_length)

refname = args.reference.split('/')[-1].split('.')[0]
w = open(args.output, 'a+')

if args.mode == 'raw':
    w.write(refname + "\t" + str(common_length) + "/" + str(expected_length) + "\n")
if args.mode == 'min' or args.mode == 'maj':
    if expected_length > 0:
        w.write(refname + "\t" + str(common) + "\t" + str(expected) + "\t" + str(common_length) + "\t" + str(expected_length) + "\t" + str(round(float(common_length)/float(expected_length), 3)) + "\n")
    else:
        w.write(refname + "\t" + str(common) + "\t" + str(expected) + "\t" + str(common_length) + "\t" + str(expected_length) + "\t0.000\n")

if args.mode == 'recomb': 
    for nb in range(len(recomb)):
        if re.search("MN908947|WG",recomb[nb][0]):
            w.write(f"{recomb[nb][0]},{recomb[nb][1]},{recomb[nb][2]},{recomb[nb][3]}\t")
        else :
            w.write(f"\n{recomb[nb]}\n")
    x=open(f"{PurePath(args.output).parent}/{Path(args.target).stem}_variant.tsv","w")
    for nb in range(len(tvar)):
        if re.search("MN908947|WG",tvar[nb][0]):
            x.write(f"{tvar[nb][0]},{tvar[nb][1]},{tvar[nb][2]},{tvar[nb][3]},{tvar[nb][3]}\t")
    x.close
w.close()

# REF="/srv/scratch/iai/seqmet/db/ref/ncov/MN908947"
# FOLDER="/srv/scratch/ONT/gridion20220318";
# for x in  $FOLDER/freebayes/BC*ONT.vcf;do
# for z in $REF/vcf/*;do 
# python3 /srv/scratch/ONT/script/recomb_V0.py -r ${z} -t "${x}" -d "$FOLDER/bam_files/${x}_bga_depth.bed" -b "$REF/midnight/MN908947_amplicon.bed"  --min_depth 20 --min_freq 0.05 --mode "recomb" --output "$FOLDER/test/${x}_recomb.tsv";
# done
# done