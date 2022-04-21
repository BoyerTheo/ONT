#!/usr/bin/env bash

#Lancement du pipeline : 
#1 : Creer un dossier pour le run
#2 : Transferer le dossier fastq_pass dans le dossier du run et le renommer en demultiplex_folder
#3 : Transformer la samplesheet en csv et la transferer dans le dossier du run
    #imperatif pour la sample samplesheet : col1:plate position col2:echantillons col3:index
#4 : changer les variables "RUN_NAME" avec le nom du dossier
#5 : changer la variable "SAMPLESHEET" par le nom de la samplesheet
#6 : lancer ce script
#7 : une fois termine, supprimer les fichiers dans le tmp avant de transferer sur ngs_viro
#8 : Si ca a plante, supprimer les fichiers genere a cette etape, et relancer (il passera les etapes deja realisees)
#9 : c'est encore complique pour les fichiers de sorties (il y en a un peu partout) et il manque deux etapes pour obtenir le fichier de validation 

CURRENT_DIR="/srv/scratch/ONT"
RUN_NAME="Grid20220415" #changer nom du run (doit correspondre au dossier)
DEMUX_FOLDER="/demultiplex_folder"
FOLDER="/srv/scratch/ONT/${RUN_NAME}" 
SCRIPT_SEQMET="/srv/scratch/iai/seqmet/script"
SCRIPT_ONT="/srv/scratch/ONT/script"
SINGU_SEQMET="/srv/scratch/iai/seqmet/singularity"
SINGU_ONT="/srv/scratch/ONT/singularity"
SAMPLESHEET="${FOLDER}/Grid20220415.csv" #mettre le csv dans le dossier du run et changer le nom du csv
NGS_VIRO="/srv/nfs/ngs-stockage/NGS_Virologie"
REF="/srv/scratch/iai/seqmet/db/ref/ncov/MN908947"
MODEL="r941_min_sup_variant_g507" #r104_e81_sup_g5015  r941_min_sup_variant_g507
KIT_ID="SQK-PBK004" #SQK-NBD112-24 SQK-RBK110-96  SQK-PBK004

#arborescence dossier :
# medaka_rerun
#   - guppyfastq
#   - bam_files
#   - freebayes
#   - medaka_output
#   - coinf
#   - conta
#   - tmp

start=$SECONDS
#Mise en place des fichiers communs
if [ ! -d $FOLDER/tmp ]; then
    mkdir -p $FOLDER/guppyfastq
    mkdir -p $FOLDER/bam_files
    mkdir -p $FOLDER/medaka_output
    mkdir -p $FOLDER/freebayes
    mkdir -p $FOLDER/tmp
    mkdir -p $FOLDER/coinf
    mkdir -p $FOLDER/fastq_filtered
    mkdir -p $FOLDER/conta
    mkdir -p $FOLDER/newval
    mkdir -p $FOLDER/recomb
fi

echo "Formatage Nom"
LST_GLIMS=()
#Formatage des noms d'echantillon pour l'ensemble du process qualité --> TypeManip+NumeroGLIMs+N°Barcode
#Récupération numéro GLIMs 
max=$(awk '{print $1}' $SAMPLESHEET | wc -l) #nombre de ligne dans la samplesheet
for num in `seq 2 $max`;
do  
    glims=$(sed -n ${num}p $SAMPLESHEET | awk -F ';' '{print $2}')  #Donne le code Glims pour chaque ligne de la samplesheet
    barcode=$(sed -n ${num}p $SAMPLESHEET | awk -F ';' '{print $3}') #Donne le Barcode pour chaque ligne de la samplesheet
    if [[ ${#barcode} -ge 4 ]]; #Verification que le champ barcode n'est pas compose de plusieurs elements (ex BC48 / RB48)
    then
        barcode=$(echo $barcode | awk -F ' / ' '{print $1}') #Recuperation du premier champ
    fi
    if [ ! -z ${glims} ] && [ ! -z ${barcode} ]; #verification que les chaines ne sont pas vides
    then
        barcode=${barcode//[!0-9]/} #formatage du barcode pour enlever les char associés ex: RB43 --> 43
        new_name="BC${barcode}-${glims}-ONT" 
        LST_GLIMS=(${LST_GLIMS[@]} $new_name)
    else 
        echo "ERREUR : Le numero glims et/ou le barcode sont vides"
    fi
done 


echo "remove barecode & primers"
for x in  ${LST_GLIMS[@]};
do
    #Création des fastas 
    if [ ! -e $FOLDER/guppyfastq/${x}.fastq ]; 
    then
        num=$(echo $x | awk -F 'BC|NC' '{print $2}' | awk -F '-' '{print $1}')
        $CURRENT_DIR/ont-guppy-cpu_5.1.12/bin/guppy_barcoder --barcode_kits $KIT_ID --input_path "$FOLDER/demultiplex_folder/barcode$num" --save_path "$FOLDER/fastq_filtered/" --recursive --records_per_fastq 0 --trim_barcodes --trim_primers
    fi
done

echo "guppyplex"
for x in  ${LST_GLIMS[@]};
do
    #Création des fastas 
    if [ ! -e $FOLDER/guppyfastq/${x}.fastq ]; 
    then
        num=$(echo $x | awk -F 'BC|NC' '{print $2}' | awk -F '-' '{print $1}')
        if [[ $KIT_ID == "SQK-RBK110-96" ]];
        then
            echo "rapide_barcoding : running with min_length=150, max_length=1400"
            singularity exec $SINGU_ONT/artic_qc.sif artic guppyplex --min-length 150 --max-length 1400 --directory $FOLDER/fastq_filtered/barcode${num} --output $FOLDER/guppyfastq/${x}.fastq
        else 
            echo "ligation : running with min_length=1000, max_length=1600"
            singularity exec $SINGU_ONT/artic_qc.sif artic guppyplex --min-length 1000 --max-length 1600 --directory $FOLDER/fastq_filtered/barcode${num} --output $FOLDER/guppyfastq/${x}.fastq
        fi
    fi
    #Mapping lectures sur génome de ref + virer les amorces d'amplicons
    if [ ! -e $FOLDER/bam_files/${x}_depth.tsv ]; 
    then
        singularity exec $SINGU_ONT/medaka-1.5.0.sif minimap2 -ax sr "$REF/MN908947.fna" "$FOLDER/guppyfastq/${x}.fastq" | samtools view --threads 1 -b -F 4 - |  tee "$FOLDER/bam_files/${x}_notrim.bam" | singularity exec $SINGU_SEQMET/samtools-1.13.sif samtools ampliconclip --hard-clip -u -o "$FOLDER/tmp/${x}_hardtrimed.bam" -b "$CURRENT_DIR/nCoV-2019.primer_2.bed" -         
        samtools sort "$FOLDER/tmp/${x}_hardtrimed.bam" -o "$FOLDER/bam_files/${x}_hardtrimed_sorted.bam"
        samtools index "$FOLDER/bam_files/${x}_hardtrimed_sorted.bam"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bedtools genomecov -ibam "$FOLDER/bam_files/${x}_hardtrimed_sorted.bam" -bga > "$FOLDER/bam_files/${x}_bga_depth.bed"
        samtools depth -aa "$FOLDER/bam_files/${x}_hardtrimed_sorted.bam" -o "$FOLDER/bam_files/${x}_depth.tsv"
    fi
done


echo "freebayes"
for x in ${LST_GLIMS[@]};
do
    if  [ ! -e $FOLDER/freebayes/${x}_cons.vcf ]; 
    then
        #freebayes pour variants majoritaires
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif freebayes --theta 0.001 --ploidy 1 --report-all-haplotype-alleles --pooled-continuous --use-best-n-alleles 2 --allele-balance-priors-off --haplotype-length 1 --use-duplicate-reads --genotyping-max-iterations 10 --genotyping-max-banddepth 4 --min-mapping-quality 1 --min-base-quality 0 -F 0.5 -C 2 --min-coverage 9 -f "$REF/MN908947.fna" -b "$FOLDER/bam_files/${x}_hardtrimed_sorted.bam" > "$FOLDER/freebayes/${x}_cons.vcf"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif vt decompose -s "$FOLDER/freebayes/${x}_cons.vcf" -o "$FOLDER/tmp/${x}_decomposed_snps.vcf"
        #Recalcul des valeurs 
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif python2.7 $SCRIPT_SEQMET/recalc_fbvcf-0.0.1.py -i "$FOLDER/tmp/${x}_decomposed_snps.vcf" -o "$FOLDER/tmp/${x}_recalc_snps.vcf" -n 0.0 -x 1.0 -t snp,del,ins,mnp,complex
    fi
done

echo "coinf"
for x in ${LST_GLIMS[@]};
do
    if [ ! -e "$FOLDER/coinf/${x}_normalized_var.vcf" ]; 
    then
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools filter -i "STRLEN(REF)==STRLEN(ALT) & AF >= 0.5 & ((SAF>=2 & SAR>=2) | AF >= 0.9) & INFO/QA>=80" "$FOLDER/tmp/${x}_recalc_snps.vcf" -o "$FOLDER/tmp/${x}_filtered_snps.vcf"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif vt normalize "$FOLDER/tmp/${x}_filtered_snps.vcf" -r "$REF/MN908947.fna" -o "$FOLDER/freebayes/${x}_normalized_snps.vcf"

        #Recréation ref avec les snps majoritaires
        echo "Nouveau consensus"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bgzip < "$FOLDER/freebayes/${x}_normalized_snps.vcf" > "$FOLDER/tmp/${x}_snps.vcf.gz"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools index "$FOLDER/tmp/${x}_snps.vcf.gz"
        cat "$REF/MN908947.fna" | singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools consensus "$FOLDER/tmp/${x}_snps.vcf.gz" > "$FOLDER/tmp/${x}_snps.fna"

        #vérifie que les snps différents de 0.5
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools filter -i 'STRLEN(REF)!=STRLEN(ALT) & AF >= 0.5 & ((SAF>=2 & SAR>=2) | AF >= 0.9) & INFO/QA>=80' "$FOLDER/tmp/${x}_recalc_snps.vcf" -o "$FOLDER/tmp/${x}_filtered_complex.vcf"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif vt normalize "$FOLDER/tmp/${x}_filtered_complex.vcf" -r "$REF/MN908947.fna" -o "$FOLDER/freebayes/${x}_normalized_complex.vcf"

        #recherche des variants de faibles fréquences
        echo "Faible frequence"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif freebayes --theta 0.001 --ploidy 1 --report-all-haplotype-alleles --pooled-continuous --use-best-n-alleles 2 --allele-balance-priors-off --haplotype-length 1 --use-duplicate-reads --genotyping-max-iterations 10 --genotyping-max-banddepth 4 --min-mapping-quality 1 --min-base-quality 0 -F 0.05 -C 2 --min-coverage 9 -f "$FOLDER/tmp/${x}_snps.fna" -b "$FOLDER/bam_files/${x}_hardtrimed_sorted.bam" > "$FOLDER/tmp/${x}_var.vcf"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif vt decompose -s "$FOLDER/tmp/${x}_var.vcf" -o "$FOLDER/tmp/${x}_decomposed_var.vcf"
        #Recalcul des valeurs 
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif python2.7 $SCRIPT_SEQMET/recalc_fbvcf-0.0.1.py -i "$FOLDER/tmp/${x}_decomposed_var.vcf" -o "$FOLDER/coinf/${x}_recalc_var.vcf" -n 0.0 -x 0.499 -t snp,del,ins,mnp,complex
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools filter -i 'AF < 0.5 & AF >= 0.05 & ( (SAF/(SAF+SAR))>=((SRF/(SRF+SRR))-0.49) & (SAF/(SAF+SAR))<=((SRF/(SRF+SRR))+0.49) & (SAF/(SAF+SAR))>=0.01 & (SAF/(SAF+SAR))<=0.99 ) & SAF>=2 & SAR>=2 & INFO/QA>=80' "$FOLDER/coinf/${x}_recalc_var.vcf" -o "$FOLDER/tmp/${x}_filtered_var.vcf"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif vt normalize "$FOLDER/tmp/${x}_filtered_var.vcf" -r "$FOLDER/tmp/${x}_snps.fna" -o "$FOLDER/coinf/${x}_normalized_var.vcf"
    fi
done

echo "preparation fichier final"
for x in ${LST_GLIMS[@]};
do 
    if [ ! -e "$FOLDER/freebayes/${x}.vcf" ];
    then
    cp "$FOLDER/freebayes/${x}_normalized_snps.vcf" "$FOLDER/tmp/${x}_unsorted_cons.vcf"
    #On récupère les variants majoritaires
    grep -v '#' "$FOLDER/freebayes/${x}_normalized_complex.vcf" >> "$FOLDER/tmp/${x}_unsorted_cons.vcf"
    singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools sort -Ov --temp-dir "$FOLDER/tmp/" -o "$FOLDER/tmp/${x}_sorted_cons.vcf" "$FOLDER/tmp/${x}_unsorted_cons.vcf"

    cp "$FOLDER/tmp/${x}_sorted_cons.vcf" "$FOLDER/tmp/${x}_sorted_concat.vcf"
    #On ajoute les variants de faibles proportions
    grep -v '#' "$FOLDER/coinf/${x}_normalized_var.vcf" >> "$FOLDER/tmp/${x}_sorted_concat.vcf"
    #On trie pour que ca soit dans le bonne ordre
    singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools sort -Ov --temp-dir "$FOLDER/tmp/" -o "$FOLDER/tmp/${x}_sorted.vcf" "$FOLDER/tmp/${x}_sorted_concat.vcf"
    #on vire d'éventuel doublon
    singularity exec $SINGU_SEQMET/varcall-0.0.5.sif vcfuniq "$FOLDER/tmp/${x}_sorted.vcf" > "$FOLDER/freebayes/${x}.vcf"
    fi
done 

echo "vcf2tsv"
for x in ${LST_GLIMS[@]};
do  
    if [ ! -e "$FOLDER/freebayes/${x}.tsv" ]; then 
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif vcf2tsv -n "." "$FOLDER/freebayes/${x}.vcf" > "$FOLDER/freebayes/${x}.tsv"
    fi
done 

echo "jolie tableau pour compter le nombre de variants"

for x in ${LST_GLIMS[@]};
do  
    sed '1d' "$FOLDER/freebayes/${x}.tsv" > "$FOLDER/tmp/${x}_varcount.tsv"
    echo -e "MN908947\t5-10%\t$(awk 'BEGIN {bp=0} $13/$15>=0.05 && $13/$15<0.1 {bp+=1} {print bp}' "$FOLDER/tmp/${x}_varcount.tsv" | tail -1)\t${x}" > "$FOLDER/freebayes/${x}_varcount.tsv"
    echo -e "MN908947\t10-20%\t$(awk 'BEGIN {bp=0} $13/$15>=0.1 && $13/$15<0.2 {bp+=1} {print bp}' "$FOLDER/tmp/${x}_varcount.tsv" | tail -1)\t${x}" >> "$FOLDER/freebayes/${x}_varcount.tsv"
    echo -e "MN908947\t20-50%\t$(awk 'BEGIN {bp=0} $13/$15>=0.2 && $13/$15<0.5 {bp+=1} {print bp}' "$FOLDER/tmp/${x}_varcount.tsv" | tail -1)\t${x}" >> "$FOLDER/freebayes/${x}_varcount.tsv"
    echo -e "MN908947\tvarcount\t$(awk 'BEGIN {bp=0} $13/$15>=0.05 && $13/$15<0.5 {bp+=1} {print bp}' "$FOLDER/tmp/${x}_varcount.tsv" | tail -1)\t${x}" >> "$FOLDER/freebayes/${x}_varcount.tsv"
done 

rm -f "$FOLDER/all_varcount.tsv"
for x in ${LST_GLIMS[@]};
do
    cat "$FOLDER/freebayes/${x}_varcount.tsv" >> "$FOLDER/all_varcount.tsv"
done 

echo "Variant table"
level=0
for x in ${LST_GLIMS[@]};
do  
    if [ ! -e $FOLDER/freebayes/${x}.txt ]; then
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif python2.7 $SCRIPT_SEQMET/vartsv_table-0.0.2.py -r "$REF/MN908947.fna" -g "$REF/MN908947.gff3" -t "$FOLDER/freebayes/${x}.tsv"  -o "$FOLDER/tmp/${x}_vartable.tsv" -d "freq"
        awk -F'[\t]' '!seen[$1]++' "$FOLDER/tmp/${x}_vartable.tsv" > "$FOLDER/freebayes/${x}_vartable.tsv"
        
        #merge table
        if [ ${level} != 0 ]; 
        then 
            touch $FOLDER/freebayes/${x}.txt
            for it in $FOLDER/freebayes/${x}_vartable.tsv; 
            do 
                sampleId="${it%%_vartable.tsv}" ; 
                cat "${it}" | while read line ; do echo -n ${line} | cut -f 1  | cut -d '|' -f 1-${level} | tr -d $'n' >> " $FOLDER/freebayes/${x}.txt"; echo -n -e '\t' >> " $FOLDER/freebayes/${x}.txt"; echo "${line}" | cut -f 2 >> " $FOLDER/freebayes/${x}.txt"; done; 
            done
            else for it in $FOLDER/freebayes/${x}_vartable.tsv; 
            do 
                sampleId="${it%%_vartable.tsv}"
                cat "${it}" >> "$FOLDER/freebayes/${x}.txt" 
            done
        fi
    fi
done

rm -f $FOLDER/all_vartable.tsv
if [ ! -e $FOLDER/all_vartable.tsv ]; then
    singularity exec $SINGU_SEQMET/varcall-0.0.5.sif python2.7 $SCRIPT_SEQMET/merge_tables-0.0.1.py $(ls $FOLDER/freebayes/*vartable.tsv) > "$FOLDER/all_vartable.tsv"
fi

echo "Recherche de conta"
echo "Creation db conta"
if [ ! -e "$FOLDER/conta/concat_full.vcf" ];
then
    for x in ${LST_GLIMS[@]};
    do
        if [[ ! -s "$FOLDER/tmp/concat_full.vcf" ]]; then
            #Création de la première ligne du fichier
            grep '#' "$FOLDER/freebayes/${x}.vcf" > "$FOLDER/tmp/concat_full.vcf"
        fi
        if [[ -s "$FOLDER/freebayes/${x}.vcf" ]]; 
        then
            #Ajout des lignes ayant une fréquence > 95%
            singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools filter -i 'AF > 0.95' "$FOLDER/freebayes/${x}.vcf" -o "$FOLDER/tmp/${x}temp.vcf"
            grep -v '#' "$FOLDER/tmp/${x}temp.vcf" >> "$FOLDER/tmp/concat_full.vcf"
        fi
    done
    if [[ ! -s "$FOLDER/tmp/concat_full.vcf" ]]; 
    then
        touch "$FOLDER/conta/${x}.vcf"
        exit 0
    fi
singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools sort -Ov -o "$FOLDER/tmp/concat_sorted.vcf" "$FOLDER/tmp/concat_full.vcf"
singularity exec $SINGU_SEQMET/varcall-0.0.5.sif vcfuniq "$FOLDER/tmp/concat_sorted.vcf" > "$FOLDER/conta/concat_full.vcf"
fi

echo "Recherche de contaminant"
if [ ! -e "$FOLDER/conta/${x}_conta_raw.tsv" ]
then
    for x in ${LST_GLIMS[@]};
    do
        for z in ${LST_GLIMS[@]};
        do
            if [[ -s "$FOLDER/freebayes/${x}.vcf" ]]; 
            then
                singularity exec $SINGU_SEQMET/varcall-0.0.5.sif python2.7 $SCRIPT_SEQMET/compare_vcf-0.0.6.py -t "$FOLDER/freebayes/${x}.vcf" -r "$FOLDER/freebayes/${z}.vcf" -d "$FOLDER/bam_files/${x}_bga_depth.bed" -b "$CURRENT_DIR/nCoV-2019.primer_2.bed" --min_depth 20 --min_freq 0.05 -m raw -o "$FOLDER/conta/${x}_conta_raw.tsv"
            fi
        done
        touch "$FOLDER/conta/${x}_conta_raw.tsv"
    done
fi 

rm -f "$FOLDER/all_conta.tsv"
for x in ${LST_GLIMS[@]};
do
    singularity exec $SINGU_SEQMET/varcall-0.0.5.sif python2.7 $SCRIPT_SEQMET/merge_tables-0.0.1.py $(ls $FOLDER/conta/*_conta_raw.tsv) > "$FOLDER/all_conta.tsv"
done 

echo "Recherche de coinfection"
for x in ${LST_GLIMS[@]};
do
    if [ ! -e "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv" ];
    then 
        for z in $REF/vcf/*.vcf; 
        do
            python3 $SCRIPT_SEQMET/compare_vcf-0.0.6.py -r ${z} -t "$FOLDER/freebayes/${x}.vcf" -d "$FOLDER/bam_files/${x}_bga_depth.bed" -b "$REF/midnight/MN908947_amplicon.bed" --min_depth 20 --min_freq 0.05 --mode "maj" --output "$FOLDER/coinf/${x}_amplicon_maj.tsv"
        done

        echo -e "MN908947\tcoinf_maj_match\t$(awk -F "\t" 'BEGIN {OFS="\t"}; $4>=6 {print $0} $4<6 {$6=0.0; print $0}' "$FOLDER/coinf/${x}_amplicon_maj.tsv" | sort -t$'\t' -k6,6 -k4,4 -nr | awk -F "\t" '{print $1}' | head -n 1 | tr -d '\n')\t${x}" >> "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv"
        echo -e "MN908947\tcoinf_maj_common\t$(awk -F "\t" 'BEGIN {OFS="\t"}; $4>=6 {print $0} $4<6 {$6=0.0; print $0}' "$FOLDER/coinf/${x}_amplicon_maj.tsv" | sort -t$'\t' -k6,6 -k4,4 -nr | awk -F "\t" '{print $4}' | head -n 1 | tr -d '\n')\t${x}" >> "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv"
        echo -e "MN908947\tcoinf_maj_ratio\t$(awk -F "\t" 'BEGIN {OFS="\t"}; $4>=6 {print $0} $4<6 {$6=0.0; print $0}' "$FOLDER/coinf/${x}_amplicon_maj.tsv" | sort -t$'\t' -k6,6 -k4,4 -nr | awk -F "\t" '{print $6}' | head -n 1 | tr -d '\n')\t${x}" >> "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv"

        MAJ_MATCH=$(cat "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv" | grep 'coinf_maj_match' | awk -F "\t" '{print $3}')

        for y in $REF/vcf/*.vcf; 
        do
                python3 $SCRIPT_SEQMET/compare_vcf-0.0.6.py -r ${y} -t "$FOLDER/freebayes/${x}.vcf" --exclusion "$REF/vcf/${MAJ_MATCH}.vcf" -d "$FOLDER/bam_files/${x}_bga_depth.bed" -b "$REF/midnight/MN908947_amplicon.bed" --min_depth 20 --min_freq 0.05 --mode "maj" --output "$FOLDER/coinf/${x}_amplicon_min.tsv"
        done

        echo -e "MN908947\tcoinf_min_match\t$(awk -F "\t" 'BEGIN {OFS="\t"}; $4>=6 {print $0} $4<6 {$6=0.0; print $0}' "$FOLDER/coinf/${x}_amplicon_min.tsv" | sort -t$'\t' -k6,6 -k4,4 -nr | awk -F "\t" '{print $1}' | head -n 1 | tr -d '\n')\t${x}" >> "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv"
        echo -e "MN908947\tcoinf_min_common\t$(awk -F "\t" 'BEGIN {OFS="\t"}; $4>=6 {print $0} $4<6 {$6=0.0; print $0}' "$FOLDER/coinf/${x}_amplicon_min.tsv" | sort -t$'\t' -k6,6 -k4,4 -nr | awk -F "\t" '{print $4}' | head -n 1 | tr -d '\n')\t${x}" >> "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv"
        echo -e "MN908947\tcoinf_min_ratio\t$(awk -F "\t" 'BEGIN {OFS="\t"}; $4>=6 {print $0} $4<6 {$6=0.0; print $0}' "$FOLDER/coinf/${x}_amplicon_min.tsv" | sort -t$'\t' -k6,6 -k4,4 -nr | awk -F "\t" '{print $6}' | head -n 1 | tr -d '\n')\t${x}" >> "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv"

        touch "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv" "$FOLDER/coinf/${x}_amplicon_min.tsv" "$FOLDER/coinf/${x}_amplicon_maj.tsv"
    fi
done 


rm -f "$FOLDER/all_coinf.tsv"
for x in ${LST_GLIMS[@]};
do
    cat "$FOLDER/coinf/${x}_MN908947_amplicon_coinf.tsv" >> "$FOLDER/all_coinf.tsv"
done 

echo "generation consensus"
for x in ${LST_GLIMS[@]};
do  
    if [ ! -e "$FOLDER/medaka_output/${x}/${x}_final_consensus.fasta" ]; then 
        #generation d'une séquence consensus brut 
        singularity exec $SINGU_ONT/medaka-1.5.0.sif medaka_consensus -i "$FOLDER/guppyfastq/${x}.fastq" -d "$REF/MN908947.fna" -o "$FOLDER/medaka_output/${x}" -t 16 -m $MODEL 
       
        #récupération des positions avec une couverture insuffisante
        awk '$4 < 20 {print $0}' "$FOLDER/bam_files/${x}_bga_depth.bed" > "$FOLDER/tmp/${x}_bga_depth.bed"
        
        #Convertion nculeotide du fichier fasta
        #"bedtools maskfasta"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bedtools maskfasta -soft -fi "$REF/MN908947.fna" -bed "$FOLDER/tmp/${x}_bga_depth.bed" -fo "$FOLDER/tmp/${x}_consensus_masked.fasta"
        #"bcftools filter : les variants majoritaires de notre vcf complet"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools filter -i 'AF >= 0.5' "$FOLDER/freebayes/${x}.vcf" -o "$FOLDER/tmp/${x}_unsorted_cons.vcf"
        #"sort on rearrange"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools sort -Ov -o "$FOLDER/tmp/${x}_sorted_cons.vcf" "$FOLDER/tmp/${x}_unsorted_cons.vcf"
        #"gzip"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bgzip < "$FOLDER/tmp/${x}_sorted_cons.vcf" > "$FOLDER/tmp/${x}_cons.vcf.gz"
        #"bcftools index"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools index "$FOLDER/tmp/${x}_cons.vcf.gz"
        #"bcftools consensus"
        cat "$FOLDER/tmp/${x}_consensus_masked.fasta" | singularity exec $SINGU_SEQMET/varcall-0.0.5.sif bcftools consensus "$FOLDER/tmp/${x}_cons.vcf.gz" > "$FOLDER/tmp/${x}_consensus.fasta"
        #"reheading" #A terme changer le header pour que ca soit le nom formaté
        sed "s/^>M[^ ]*/>${x}/" "$FOLDER/tmp/${x}_consensus.fasta" > "$FOLDER/tmp/${x}_rehead.fna"
        #"script python"
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif python2.7 $SCRIPT_SEQMET/fasta_masklowercase-0.0.1.py -i "$FOLDER/tmp/${x}_rehead.fna" -o "$FOLDER/medaka_output/${x}/${x}_final_consensus.fasta"
    fi
done

echo "Quality"
echo "Samtools coverage"
for x in ${LST_GLIMS[@]};
do  
    if [ ! -e $FOLDER/bam_files/${x}.cov ];then
        for line in $(samtools coverage $FOLDER/bam_files/${x}_hardtrimed_sorted.bam | tail -n +2 | tr '\t' ','); do SEG="${line%%,*}"; echo -e "MN908947\tmean_depth_${SEG}\t$(echo $line | cut -d ',' -f 7)\t${x}" >> "${FOLDER}/bam_files/${x}.cov"; done;
    fi
done

rm -f $FOLDER/all_cov.cov
cat $FOLDER/bam_files/*.cov >> $FOLDER/all_cov.cov

echo "Consensus"
rm -f $FOLDER/all_consensus.fasta 
cat $FOLDER/medaka_output/*/*_final_consensus.fasta >> $FOLDER/all_consensus.fasta                               

echo "Summary"
if [ ! -e $FOLDER/all_summary.tsv ];
then
    for x in ${LST_GLIMS[@]};
    do
        singularity exec $SINGU_SEQMET/varcall-0.0.5.sif /usr/bin/Rscript $SCRIPT_SEQMET/make_summary-0.0.6.R --output "$FOLDER/tmp/${x}_summary.tsv" --runname $RUN_NAME --coverage "$FOLDER/bam_files/${x}.cov" --cons "$FOLDER/medaka_output/${x}/${x}_final_consensus.fasta" --count-table "$FOLDER/freebayes/${x}_varcount.tsv" --conta "$FOLDER/conta/${x}_conta_raw.tsv" --coinf-table "$FOLDER/all_coinf.tsv"
    done
    echo "Merge summary"
    echo -e 'sample\treference\tmean_depth\tpercCOV\thasdp\thasposc\t5-10%\t10-20%\t20-50%\tvarcount\tcoinf_maj_match\tcoinf_maj_common\tcoinf_maj_ratio\tcoinf_min_match\tcoinf_min_common\tcoinf_min_ratio\tRUN'> $FOLDER/all_summary.tsv
    for x in ${LST_GLIMS[@]};
    do
        awk -v value="NA\tFAILED" -v row=2 -v col=5 'BEGIN{FS=OFS="\t"} NR==row {$col=value}1' $FOLDER/tmp/${x}_summary.tsv >> $FOLDER/tmp/${x}_new_summary.tsv
        cat $FOLDER/tmp/${x}_new_summary.tsv | awk FNR==2 >> $FOLDER/all_summary.tsv
    done
fi

echo "Nextclade"
if [ ! -e $FOLDER/nextclade_report.tsv ]; then
    singularity exec $SINGU_SEQMET/nextclade-1.11.0.sif nextclade run --input-dataset /srv/scratch/iai/seqmet/db/nextclade/sars-cov-2_220408/ --include-reference --in-order --input-fasta $FOLDER/all_consensus.fasta --output-tsv $FOLDER/nextclade_report.tsv -d $FOLDER/tmp/
fi

echo "Pangolin"
if [ ! -e $FOLDER/pangolin_report.csv ]; then 
    singularity exec $SINGU_SEQMET/pangolin-4.0.1.sif pangolin $FOLDER/all_consensus.fasta --outfile $FOLDER/'pangolin_report.csv' --use-assignment-cache --assignment-cache /srv/scratch/iai/seqmet/db/pangolin/usher_assignments-v1.2.133.cache.csv.gz --tempdir $FOLDER/tmp/
fi

echo "Validation"
singularity exec $SINGU_SEQMET/varcall-0.0.5.sif python3 $SCRIPT_SEQMET/gen_techval-0.1.0.py --nextcladeversion "(Nextclade v.1.11.0)" --rundate "22020415" --mode "ncov" --outdir "$FOLDER/newval/" --varcount_threshold 13 --dp_threshold 8 --cov_minok 90 --cov_maxneg 5 --summary $FOLDER/all_summary.tsv --nextclade $FOLDER/nextclade_report.tsv --pangolin $FOLDER/pangolin_report.csv --vartable $FOLDER/all_vartable.tsv --expectedmatrix $REF/matricemutlineage_ncov_expected_v1-220404.tsv --likelymatrix $REF/matricemutlineage_ncov_likely_v1-220404.tsv --gff $REF/MN908947.gff3 
#old=gen_techval-0.0.7.py

echo "Recherche recombinaison"
for x in ${LST_GLIMS[@]}; 
do
    for z in $REF/vcf/*.vcf; 
    do
        python3 $SCRIPT_ONT/recomb_V0.py -r ${z} -t "$FOLDER/freebayes/${x}.vcf" -d "$FOLDER/bam_files/${x}_bga_depth.bed" -b "$REF/midnight/MN908947_amplicon.bed" --min_depth 20 --min_freq 0.05 --mode recomb --output "$FOLDER/recomb/${x}_recomb.tsv"
    done
    python3 $SCRIPT_ONT/find_recomb.py -r "$FOLDER/recomb/${x}_recomb.tsv" -s "$FOLDER/recomb/${x}_variant.tsv" -d "$FOLDER/bam_files/${x}_depth.tsv" -o "$FOLDER/recomb/${x}";
done


#rm -rf $FOLDER/tmp/*

duration=$(( SECONDS - start ))
# Affichage du temps 
if (( $duration > 3600 )) ; then
    let "hours=SECONDS/3600"
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo "Completed in $hours hour(s), $minutes minute(s) and $seconds second(s)" 
elif (( $duration > 60 )) ; then
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo "Completed in $minutes minute(s) and $seconds second(s)"
else
    echo "Completed in $duration seconds"
fi