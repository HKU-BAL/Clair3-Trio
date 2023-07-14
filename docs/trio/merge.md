# Merge trio VCF files

Clair3-Trio outputs variants in separate VCF files and requires additional steps to merge multiple outputs. Merging multiple VCF files is not a trivial task. The difficulty arises from two main issues:

1) determining whether a non-called site should be considered a reference call (0/0) or an undetermined call (./.) based on the sequencing data;
2) merging multiple alleles and overlapping variants between individuals and positions.

To address the challenge of merging VCF files, we tried different methods. This page will describe the methods we applied to merge VCF files for Clair3-Trio.

When merging VCF, We recommend **turning on `--print_ref_calls`** for Clair3-Trio to output all reference calls for better merging multiple VCF files.

## Contents

* [Option 1. Merge VCF with unknown to "0/0"](#option-1-merge-vcf-with-unknown-to-00)
* [Option 2. Merge VCF with unknown to "./."](#option-2-merge-vcf-with-unknown-to-)
* [Results and differents between Option 1 and Option 2](#results-and-differents-between-option-1-and-option-2)
* [Option 3. Merging with GVCF files](#option-3-merging-with-gvcf-files)
* [Quality of de novo variants](#quality-of-de-novo-variants)

## Option 1. Merge VCF with unknown to "0/0"

Merge the VCF files, and set all non-called sites to "0/0".

This is the method that we benchmarked in our paper, and it has a high recall to detect the de novo (0/0+0/0 -> 0/1) variants.

```
BCFTOOLS=bcftools
RTG=rtg
C_VCF=[CHILD_VCF]			#e.g. HG002.vcf.gz
P1_VCF=[PARENET_1_VCF]		#e.g. HG003.vcf.gz
P2_VCF=[PARENET_2_VCF]		#e.g. HG004.vcf.gz

_TRIO_PED=[TRIO PED FILE]
# e.g. cat $_TRIO_PED
##PED format pedigree
##fam-id ind-id pat-id mat-id sex phen
#1 HG002 HG003 HG004 1 0
#1 HG003 0 0 1 0
#1 HG004 0 0 2 0

M_VCF=TRIO.vcf.gz
M_VCF_annotated=TRIO_ann.vcf.gz

N_C_VCF=[CHILD NANE]_n.vcf.gz 
N_P1_VCF=[P1 NAME]_n.vcf.gz 
N_P2_VCF=[P2 NAME]_n.vcf.gz

# filter Lowqual call           
${BCFTOOLS} filter -e 'FILTER="LowQual" | FORMAT/GQ < 2' ${C_VCF} | ${BCFTOOLS} view -O z -o ${N_C_VCF}
${BCFTOOLS} index ${N_C_VCF}
${BCFTOOLS} filter -e 'FILTER="LowQual" | FORMAT/GQ < 2' ${P1_VCF} | ${BCFTOOLS} view -O z -o ${N_P1_VCF}
${BCFTOOLS} index ${N_P1_VCF}
${BCFTOOLS} filter -e 'FILTER="LowQual" | FORMAT/GQ < 2' ${P2_VCF} | ${BCFTOOLS} view -O z -o ${N_P2_VCF}
${BCFTOOLS} index ${N_P2_VCF}
         
# merge vcf
${BCFTOOLS} merge ${N_C_VCF} ${N_P1_VCF} ${N_P2_VCF} --threads 32 -m all -F x -0 | ${BCFTOOLS} view -O z -o ${M_VCF}
${BCFTOOLS} index ${M_VCF}

# add annotation
REF_SDF_FILE_PATH=[REFERENCE SDF FILE]
#${RTG} mendelian --all-records -i ${M_VCF} -o ${M_VCF_annotated} --pedigree ${_TRIO_PED} -t ${REF_SDF_FILE_PATH}
```


## Option 2. Merge VCF with unknown to "./."

Merge the VCF files, and set all non-called / unknown sites to "./.".

Keeping unknown sites as "./." will have a high precision rate to detect de novo (0/0+0/0 -> 0/1) variants. But may leave some sites uncomparable for incomplete records.

To keep unknown as "./." in the merged file, we can do this for Clair3-Trio:

```
BCFTOOLS=bcftools
RTG=rtg
C_VCF=[CHILD_VCF]			#e.g. HG002.vcf.gz
P1_VCF=[PARENET_1_VCF]		#e.g. HG003.vcf.gz
P2_VCF=[PARENET_2_VCF]		#e.g. HG004.vcf.gz

_TRIO_PED=[TRIO PED FILE]
# e.g. cat $_TRIO_PED
##PED format pedigree
##fam-id ind-id pat-id mat-id sex phen
#1 HG002 HG003 HG004 1 0
#1 HG003 0 0 1 0
#1 HG004 0 0 2 0

M_VCF=TRIO.vcf.gz
M_VCF_annotated=TRIO_ann.vcf.gz

N_C_VCF=[CHILD NANE]_n.vcf.gz 
N_P1_VCF=[P1 NAME]_n.vcf.gz 
N_P2_VCF=[P2 NAME]_n.vcf.gz

# mask Lowqual call        
${BCFTOOLS} filter -S . -e 'FILTER="LowQual" | FORMAT/GQ < 2' ${C_VCF} | ${BCFTOOLS} view -O z -o ${N_C_VCF}
${BCFTOOLS} index ${N_C_VCF}
${BCFTOOLS} filter -S . -e 'FILTER="LowQual" | FORMAT/GQ < 2' ${P1_VCF} | ${BCFTOOLS} view -O z -o ${N_P1_VCF}
${BCFTOOLS} index ${N_P1_VCF}
${BCFTOOLS} filter -S . -e 'FILTER="LowQual" | FORMAT/GQ < 2' ${P2_VCF} | ${BCFTOOLS} view -O z -o ${N_P2_VCF}
${BCFTOOLS} index ${N_P2_VCF}

# merge
${BCFTOOLS} merge ${N_C_VCF} ${N_P1_VCF} ${N_P2_VCF} --threads 32 -m all | ${BCFTOOLS} view -O z -o ${M_VCF}
${BCFTOOLS} index ${M_VCF}

# add annotation
REF_SDF_FILE_PATH=[REFERENCE SDF FILE]
#${RTG} mendelian --all-records -i ${M_VCF} -o ${M_VCF_annotated} --pedigree ${_TRIO_PED} -t ${REF_SDF_FILE_PATH}
```


## Results and differents between Option 1 and Option 2


Here is an example of a merged VCF: [chr20_trio_clair3_trio_unknown_to_0z](http://www.bio8.cs.hku.hk/clair3_trio/analysis_result/7_merge/merge_unknown_to_0/).

Here is an example of a merged VCF: [chr20_trio_clair3_trio_unknown_to_dot](http://www.bio8.cs.hku.hk/clair3_trio/analysis_result/7_merge/merge_unknown_to_dot/).

Below is the results of using different merge methods tested in HG002 trio, chr20.

|Tool|Sample|Contig|Coverage|Option 1: unknow to 0/0| | | |Option 2: unknow to ./.| | | |
|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|
| | | | |Number of MCV|Number of MCU|Number compared|de novo|Number of MCV|Number of MCU|Number compared|de novo|
|Clair3-Trio|HG002 trio|chr20|10x|7049|0|172899|1991|5064|8967|177844|1015|
|Clair3-Trio|HG002 trio|chr20|20x|6310|0|174932|1240|3233|6783|176541|412|
|Clair3-Trio|HG002 trio|chr20|30x|6041|0|174241|1102|3026|6452|175756|374|
|Clair3-Trio|HG002 trio|chr20|40x|6513|0|172255|935|3005|6489|173287|316|
|Clair3-Trio|HG002 trio|chr20|50x|6269|0|170096|915|2847|6334|171205|299|
|Clair3-Trio|HG002 trio|chr20|60x|6043|0|168922|971|2966|6652|170816|380|
|Clair3-Trio|HG002 trio|chr20|70x|5881|0|167697|976|2883|7035|170053|377|

Using Option 2 compared to Option 1, de novo variants were reduced by 40% to 60%.


## Option 3. Merging with GVCF files

GVCF keeps all positions in files and provides more information for the not called sites. However, when tested merging GVCF with available tools, we found many cases are problematic, which requires further finetuning for the results to make it useable. We are currently exploring available methods for merging GVCF.                                                                 
TBC

## Quality of de novo variants

All de novo variants in the merged files are labeled with the "MCV" INFO tag in the merged files.

We can use `${BCFTOOLS} view -i 'INFO/MCV="HG002:0/0+0/0->0/1"' ${M_VCF_annotated} | less` to check all de novo variants.

And the quality of de novo variants can be accessed by aggregating (by max/mean or other methods) each individual's genotype quality (FORMAT/GQ).


