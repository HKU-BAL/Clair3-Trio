## ONT Variant Calling Quick Demo
Here is a quick demo for the Oxford Nanopore (ONT) variant calling using GIAB HG002, HG003 and HG004 chromosome 20 data. 

```bash
Sample:     GIAB HG002, HG003, HG004
Coverage:   ~60x
Reference:  GRCh38_no_alt
Region:     chr20:100000-300000
Basecaller: Guppy 4.2.2
```

### **Download data**

```bash
# Parameters
_INPUT_DIR="${HOME}/clair3_trio_quick_demo"         # note that Absolute path is needed
_OUTPUT_DIR="${_INPUT_DIR}/output"
_TAR_URL="http://www.bio8.cs.hku.hk/clair3_trio/demo/"

mkdir -p ${_INPUT_DIR}
mkdir -p ${_OUTPUT_DIR}

# Download quick demo data
# GRCh38_no_alt reference
wget -P ${_INPUT_DIR} ${_TAR_URL}/GRCh38_no_alt_chr20.fa
wget -P ${_INPUT_DIR} ${_TAR_URL}/GRCh38_no_alt_chr20.fa.fai
# Child's BAM chr20:100000-300000
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG002.bam
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG002.bam.bai
# Parent1's BAM chr20:100000-300000
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG003.bam
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG003.bam.bai
# Parent2's BAM chr20:100000-300000
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG004.bam
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG004.bam.bai
# GIAB Truth VCF and BED for child
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG002_GRCh38_20_v4.2.1_benchmark.vcf.gz
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG002_GRCh38_20_v4.2.1_benchmark.vcf.gz.tbi
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG002_GRCh38_20_v4.2.1_benchmark_noinconsistent.bed
# GIAB Truth VCF and BED for parent1
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG003_GRCh38_20_v4.2.1_benchmark.vcf.gz
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG003_GRCh38_20_v4.2.1_benchmark.vcf.gz.tbi
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG003_GRCh38_20_v4.2.1_benchmark_noinconsistent.bed
# GIAB Truth VCF and BED for parent2
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG004_GRCh38_20_v4.2.1_benchmark.vcf.gz
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG004_GRCh38_20_v4.2.1_benchmark.vcf.gz.tbi
wget -P ${_INPUT_DIR} ${_TAR_URL}/HG004_GRCh38_20_v4.2.1_benchmark_noinconsistent.bed

_BAM_C=${_INPUT_DIR}/HG002.bam        
_BAM_P1=${_INPUT_DIR}/HG003.bam          
_BAM_P2=${_INPUT_DIR}/HG004.bam          
_SAMPLE_C="HG002"
_SAMPLE_P1="HG003"
_SAMPLE_P2="HG004"
_REF=${_INPUT_DIR}/GRCh38_no_alt_chr20.fa
_THREADS="36"                

_CONTIGS="chr20"
START_POS=100000
END_POS=300000
echo -e "${_CONTIGS}\t${START_POS}\t${END_POS}" > ${_INPUT_DIR}/quick_demo.bed
```

### **Run Clair3-Trio**

```bash
MODEL_C3='ont'
MODEL_C3T='c3t_hg002_g422'

docker run -it \
  -v ${_INPUT_DIR}:${_INPUT_DIR} \
  -v ${_OUTPUT_DIR}:${_OUTPUT_DIR} \
  hkubal/clair3-trio:latest \
  /opt/bin/run_clair3_trio.sh \
  --ref_fn=${_REF} \
  --bam_fn_c=${_BAM_C} \
  --bam_fn_p1=${_BAM_P1} \
  --bam_fn_p2=${_BAM_P2} \
  --sample_name_c=${_SAMPLE_C} \
  --sample_name_p1=${_SAMPLE_P1} \
  --sample_name_p2=${_SAMPLE_P2} \
  --threads=${_THREADS} \
  --model_path_clair3="/opt/models/clair3_models/${MODEL_C3}" \
  --model_path_clair3_trio="/opt/models/clair3_trio_models/${MODEL_C3T}" \
  --output=${_OUTPUT_DIR}

```

### **Run hap.py for benchmarking (optional)**

```bash
BASELINE_VCF_FILE_PATH_C="HG002_GRCh38_20_v4.2.1_benchmark.vcf.gz"
BASELINE_BED_FILE_PATH_C="HG002_GRCh38_20_v4.2.1_benchmark_noinconsistent.bed"
BASELINE_VCF_FILE_PATH_P1="HG003_GRCh38_20_v4.2.1_benchmark.vcf.gz"
BASELINE_BED_FILE_PATH_P1="HG003_GRCh38_20_v4.2.1_benchmark_noinconsistent.bed"
BASELINE_VCF_FILE_PATH_P2="HG004_GRCh38_20_v4.2.1_benchmark.vcf.gz"
BASELINE_BED_FILE_PATH_P2="HG004_GRCh38_20_v4.2.1_benchmark_noinconsistent.bed"
OUTPUT_VCF_FILE_C="HG002.vcf.gz"
OUTPUT_VCF_FILE_P1="HG003.vcf.gz"
OUTPUT_VCF_FILE_P2="HG004.vcf.gz"



# Benchmark variant calls against truth set with hap.py
mkdir -p ${_OUTPUT_DIR}/hap

BASELINE_VCF_FILE_PATH=${BASELINE_VCF_FILE_PATH_C}
BASELINE_BED_FILE_PATH=${BASELINE_BED_FILE_PATH_C}
OUTPUT_VCF_FILE_PATH=${OUTPUT_VCF_FILE_C}
HAPPY_PATH=hap/${_SAMPLE_C}_happy

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
${_OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
-f "${_INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
-r "${_REF}" \
-o "${_OUTPUT_DIR}/${HAPPY_PATH}" \
-l ${_CONTIGS}:${START_POS}-${END_POS} \
--engine=vcfeval \
--threads="8" \
--pass-only

BASELINE_VCF_FILE_PATH=${BASELINE_VCF_FILE_PATH_P1}
BASELINE_BED_FILE_PATH=${BASELINE_BED_FILE_PATH_P1}
OUTPUT_VCF_FILE_PATH=${OUTPUT_VCF_FILE_P1}
HAPPY_PATH=hap/${_SAMPLE_P1}_happy

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
${_OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
-f "${_INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
-r "${_REF}" \
-o "${_OUTPUT_DIR}/${HAPPY_PATH}" \
-l ${_CONTIGS}:${START_POS}-${END_POS} \
--engine=vcfeval \
--threads="8" \
--pass-only


BASELINE_VCF_FILE_PATH=${BASELINE_VCF_FILE_PATH_P2}
BASELINE_BED_FILE_PATH=${BASELINE_BED_FILE_PATH_P2}
OUTPUT_VCF_FILE_PATH=${OUTPUT_VCF_FILE_P2}
HAPPY_PATH=hap/${_SAMPLE_P2}_happy

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
${_OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
-f "${_INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
-r "${_REF}" \
-o "${_OUTPUT_DIR}/${HAPPY_PATH}" \
-l ${_CONTIGS}:${START_POS}-${END_POS} \
--engine=vcfeval \
--engine=vcfeval \
--threads="8" \
--pass-only

```

### **Calculate number of Mendelian violation and de novo variants**

```
mkdir -p ${_OUTPUT_DIR}/trio
M_VCF=${_OUTPUT_DIR}/trio/${_SAMPLE_C}_TRIO.vcf.gz
M_VCF_annotated=${_OUTPUT_DIR}/trio/${_SAMPLE_C}_TRIO_ann.vcf.gz

# merge predicted VCFs
docker run \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
staphb/bcftools:1.12 bcftools merge \
${_OUTPUT_DIR}/${_SAMPLE_C}.vcf.gz \
${_OUTPUT_DIR}/${_SAMPLE_P1}.vcf.gz \
${_OUTPUT_DIR}/${_SAMPLE_P2}.vcf.gz \
--threads 8 -f PASS -0 -m all -O z -o ${M_VCF}

docker run \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
staphb/bcftools:1.12 bcftools index ${M_VCF}

# data preparing
docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
realtimegenomics/rtg-tools:3.12.1 format \
-o ${_INPUT_DIR}/GRCh38_no_alt_analysis_set.sdf "${_REF}"

# generarte BED file
FILE="${_OUTPUT_DIR}/trio.ped"
cat <<EOM >$FILE
#PED format pedigree
#fam-id ind-id pat-id mat-id sex phen
1 HG002 HG003 HG004 1 0
1 HG003 0 0 1 0
1 HG004 0 0 2 0
EOM

# get Mendelian vilations
_TRIO_PED=${_OUTPUT_DIR}/trio.ped
REF_SDF_FILE_PATH=${_INPUT_DIR}/GRCh38_no_alt_analysis_set.sdf 
docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
realtimegenomics/rtg-tools:3.12.1 mendelian \
-i ${M_VCF} -o ${M_VCF_annotated} --pedigree ${_TRIO_PED} -t ${REF_SDF_FILE_PATH} |& tee ${_OUTPUT_DIR}/trio/MDL.log

# benchmark de novo variants
# note that checking de novo variants require a region specific in ${_CONTIGS}

# merge trio's bed file
BED2=${_INPUT_DIR}/${BASELINE_BED_FILE_PATH_C}
BED3=${_INPUT_DIR}/${BASELINE_BED_FILE_PATH_P1}
BED4=${_INPUT_DIR}/${BASELINE_BED_FILE_PATH_P2}
_TRIO_BED_PATH=${_INPUT_DIR}/020304.bed

docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools intersect -a ${BED2} -b ${BED3} > ${_INPUT_DIR}/tmp_out
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools sort -i ${_INPUT_DIR}/tmp_out > ${_INPUT_DIR}/0203.bed
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools intersect -a ${_INPUT_DIR}/0203.bed -b ${BED4} > ${_INPUT_DIR}/tmp_out
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools sort -i ${_INPUT_DIR}/tmp_out > ${_TRIO_BED_PATH}

# merge trio's true set
_TRIO_GIAB_MERGED=${_INPUT_DIR}/${_SAMPLE_C}_TRIO.vcf.gz

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
staphb/bcftools:1.12 bcftools merge \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH_C} \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH_P1} \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH_P2} \
--threads 8 -f PASS -0 -m all -O z -o ${_TRIO_GIAB_MERGED}

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
staphb/bcftools:1.12 bcftools index ${_TRIO_GIAB_MERGED}

# get de nove variants
docker run -it \
  -v ${_INPUT_DIR}:${_INPUT_DIR} \
  -v ${_OUTPUT_DIR}:${_OUTPUT_DIR} \
  hkubal/clair3-trio:latest \
  python /opt/bin/clair3.py \
  Check_de_novo --call_vcf ${M_VCF} --ctgName ${_CONTIGS} --bed_fn ${_TRIO_BED_PATH} --true_vcf ${_TRIO_GIAB_MERGED} 
```



### Hap.py Expected output:


```bash
Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio
INDEL,ALL,55,38,17,52,3,11,0,3,0.690909,0.926829,0.211538,0.791667,,,2.375,2.25
INDEL,PASS,55,38,17,52,3,11,0,3,0.690909,0.926829,0.211538,0.791667,,,2.375,2.25
SNP,ALL,389,389,0,420,1,29,0,1,1.0,0.997442,0.069048,0.99872,2.5045045045,2.38709677419,3.98717948718,4.12195121951
SNP,PASS,389,389,0,420,1,29,0,1,1.0,0.997442,0.069048,0.99872,2.5045045045,2.38709677419,3.98717948718,4.12195121951

Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio
INDEL,ALL,59,43,16,63,1,19,0,1,0.728814,0.977273,0.301587,0.834951,,,4.09090909091,3.76923076923
INDEL,PASS,59,43,16,63,1,19,0,1,0.728814,0.977273,0.301587,0.834951,,,4.09090909091,3.76923076923
SNP,ALL,402,402,0,466,1,62,0,0,1.0,0.997525,0.133047,0.998761,2.5350877193,2.31205673759,4.64788732394,4.47058823529
SNP,PASS,402,402,0,466,1,62,0,0,1.0,0.997525,0.133047,0.998761,2.5350877193,2.31205673759,4.64788732394,4.47058823529

Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio
INDEL,ALL,37,25,12,37,2,9,0,1,0.675676,0.928571,0.243243,0.78219,,,0.478260869565,0.384615384615
INDEL,PASS,37,25,12,37,2,9,0,1,0.675676,0.928571,0.243243,0.78219,,,0.478260869565,0.384615384615
SNP,ALL,213,212,1,224,0,12,0,0,0.995305,1.0,0.053571,0.997647,2.73684210526,2.61290322581,0.468965517241,0.493333333333
SNP,PASS,213,212,1,224,0,12,0,0,0.995305,1.0,0.053571,0.997647,2.73684210526,2.61290322581,0.468965517241,0.493333333333
```

### Check Mendelian and de novo variants Expected output:

```bash
Family: [HG003 + HG004] -> [HG002]
Concordance HG002: F:565/566 (99.82%)  M:558/566 (98.59%)  F+M:554/566 (97.88%)
0/566 (0.00%) records did not conform to expected call ploidy
12/566 (2.12%) records contained a violation of Mendelian constraints

in chr20, called de novo sites 1
in chr20, true de novo sites 39
FN, TP, FP: 39,0,1
```
