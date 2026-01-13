#!/bin/bash
set -e
#1. Data QC
#-s bash

#qc and remove polyg tails

array=("C_06h_f1" "C_06h_f2" "C_06h_f3" "C_06h_f4" \
        "C_24h_f1" "C_24h_f2" "C_24h_f3" "C_24h_f4" \
        "C_72h_f1" "C_72h_f2" "C_72h_f3" "C_72h_f4" \
        "E_06h_f1" "E_06h_f2" "E_06h_f3" "E_06h_f4" \
        "E_24h_f1" "E_24h_f2" "E_24h_f3" "E_24h_f4" \
        "E_72h_f1" "E_72h_f2" "E_72h_f3" "E_72h_f4")
for i in "${array[@]}";
do 
fastp -g -i \
    -i ../reads/"$i"_R1.fastq \
    -I ../reads/"$i"_R2.fastq \
    -o "$i"_c_R1.fastq.gz \
    -O "$i"_c_R2.fastq.gz \
    -h "$i".html -j "$i".json =w 9
done
conda deactivate

#create a kallisto index from a ffn file of both of my genomes
conda activate kallisto
kallisto index -i dual.index.idx peri.micro.ffn

array=("C_06h_f1" "C_06h_f2" "C_06h_f3" "C_06h_f4" \
        "C_24h_f1" "C_24h_f2" "C_24h_f3" "C_24h_f4" \
        "C_72h_f1" "C_72h_f2" "C_72h_f3" "C_72h_f4" \
        "E_06h_f1" "E_06h_f2" "E_06h_f3" "E_06h_f4" \
        "E_24h_f1" "E_24h_f2" "E_24h_f3" "E_24h_f4" \
        "E_72h_f1" "E_72h_f2" "E_72h_f3" "E_72h_f4")
for i in "${array[@]}";
do 
kallisto quant \
    -i dual.index.idx \
    -o "$i" \
    -t 14 -b 100 \
    ../clean/"$i"_c_R1.fastq.gz ../clean/"$i"_c_R2.fastq.gz ;
done
conda deactivate


mkdir ../abundances

#copy all abundance files to one folder and rename them per sample name
for f in ../output/*
do
	base=$(basename $f)
	cp "$f"/abundance.tsv ../abundances/$base.abund.tsv
done

#made an extra file that contains the reads for each set of replicates
paste C_06h_f1.abund.tsv C_06h_f2.abund.tsv C_06h_f3.abund.tsv C_06h_f4.abund.tsv > C_06h.abund.tsv
paste C_24h_f1.abund.tsv C_24h_f2.abund.tsv C_24h_f3.abund.tsv C_24h_f4.abund.tsv > C_24h.abund.tsv
paste C_72h_f1.abund.tsv C_72h_f2.abund.tsv C_72h_f3.abund.tsv C_72h_f4.abund.tsv > C_72h.abund.tsv
paste E_06h_f1.abund.tsv E_06h_f2.abund.tsv E_06h_f3.abund.tsv E_06h_f4.abund.tsv > E_06h.abund.tsv
paste E_24h_f1.abund.tsv E_24h_f2.abund.tsv E_24h_f3.abund.tsv E_24h_f4.abund.tsv > E_24h.abund.tsv
paste E_72h_f1.abund.tsv E_72h_f2.abund.tsv E_72h_f3.abund.tsv E_72h_f4.abund.tsv > E_72h.abund.tsv