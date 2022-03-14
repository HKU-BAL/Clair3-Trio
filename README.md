# Clair3-Trio: high-performance Nanopore long-reads variant calling in family trios with Trio-to-Trio deep neural networks


[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)  

Contact: Ruibang Luo, Junhao Su

Email: rbluo@cs.hku.hk, jhsu@cs.hku.hk  

----

## Introduction

Clair3-Trio is a variants caller tailored for family trios from nanopore long-reads. Clair3-Trio employs a Trio-to-Trio deep neural network model that allows it to input all trio’s sequencing information and output all trio’s predicted variants within a single model, to perform far better variant calling. We also present MCVLoss, the first loss function that can improve variants calling in trios by leveraging the explicitly encoding of the priors of the Mendelian inheritance in trios. Clair3-Trio showed comprehensive improvement in experiments. It predicted much fewer Mendelian inheritance violation variations than current state-of-the-art methods.

Clair3-Trio is fully based on [Clair3](https://github.com/HKU-BAL/Clair3)

----

## Contents

* [Introduction](#introduction)
* [Latest Updates](#latest-updates)
* [What's New in Clair3](#whats-new-in-clair3)
* [Quick Demo](#quick-demo)
* [Pre-trained Models](#pre-trained-models)
  * [Guppy4 Model](#pre-trained-models)
* [Installation](#installation)
  + [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  + [Option 2. Singularity](#option-2-singularity)
  + [Option 3. Bioconda](#option-3--bioconda)
  + [Option 4. Build an anaconda virtual environment](#option-4-build-an-anaconda-virtual-environment)
  + [Option 5. Docker Dockerfile](#option-5-docker-dockerfile)
* [Usage](#usage)
* [Folder Structure and Submodule Descriptions](#folder-structure-and-submodule-descriptions)
* [Training Data](docs/training_data.md)
* [VCF/GVCF Output Formats](#vcfgvcf-output-formats)
* [Pileup Model Training](docs/pileup_training.md)
* [Full-Alignment Model Training](docs/full_alignment_training_r1.md)
* [Representation Unification](docs/representation_unification.md)
* [Visualization](docs)
  * [Model Input](docs/model_input_visualization.md)
  * [Representation Unification](docs/representation_unification_visualization.md)

----

## Latest Updates

*v0.1 (March 12, 2022)*: Initial release.

---


## What's New in Clair3-Trio

* **New Architecture.** Clair3-Trio employs a Trio-to-Trio deep neural network model that allows it to input all trio’s sequencing information and output all trio’s predicted variants within a single model.
* **Mendelian violations - Aware.**  Clair3-Trio uses MCVLoss to improve variants calling in trios by leveraging the explicitly encoding of the priors of the Mendelian inheritance in trios. 
* **Improved Performance.** Using HG002 trio 10-fold coverage ONT data for benchmarking, Clair-Trio achieved 97.30% SNP F1-score and 56.48% Indel F1-score. Compare to Clair3, Clair3-Trio reduced SNP errors by **~78%**,  and Indel errors by **~22%**. Clair3-Trio also reduced Mendelian violations from 48345 to 7072 in the family.

----

## Quick Demo

*   see [Trio Quick Demo](docs/trio/trio_quick_demo.md).

----

## Pre-trained Models

### HKU-provided Models

Download models from [here](http://www.bio8.cs.hku.hk/clair3_trio/clair3_trio_models/) or click on the links below.

In a docker installation, models are in `/opt/models/`. In a bioconda installation, models are in `{CONDA_PREFIX}/bin/models/`.

|           Model name           |  Platform   |                       Training samples                       | Included in the bioconda package | Included in the docker image |   Date   |  Basecaller  | File                                |                             Link                             |
| :----------------------------: | :---------: | :----------------------------------------------------------: | -------------------------------- | :--------------------------: | :------: | :----------: | ----------------------------------- | :----------------------------------------------------------: |
|    r941_prom_hac_g422     |     ONT     |                         HG001,2,3       | Yes                              |             Yes              | 20220312 | Guppy4 hac | r941_prom_hac_g422.tar.gz      | [Download](http://www.bio8.cs.hku.hk/clair3_trio/clair3_trio_models/r941_prom_hac_g422.tar.gz) |


----

## Installation

### Option 1. Build an anaconda virtual environment

**Anaconda install**:

Please install anaconda using the official [guide](https://docs.anaconda.com/anaconda/install) or using the commands below:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh
```

**Install Clair3 env using anaconda step by step:**


```bash
# create and activate an environment named clair3
conda create -n clair3 python=3.6.10 -y
source activate clair3

# install pypy and packages in the environemnt
conda install -c conda-forge pypy3.6 -y
pypy3 -m ensurepip
pypy3 -m pip install mpmath==1.2.1

# install python packages in environment
pip3 install tensorflow==2.2.0
pip3 install tensorflow-addons==0.11.2 tables==3.6.1
conda install -c anaconda pigz==2.4 -y
conda install -c conda-forge parallel=20191122 zstd=1.4.4 -y
conda install -c conda-forge -c bioconda samtools=1.10 -y
conda install -c conda-forge -c bioconda whatshap=1.0 -y

# clone Clair3-Trio
git clone https://github.com/HKU-BAL/Clair3-Trio.git
cd Clair3-Trio

# download Clair3's pre-trained models
mkdir clair3_models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
tar -zxvf clair3_models.tar.gz -C ./cliar3_models


# download Clair3-Trio's pre-trained models
mkdir clair3_trio_models
wget http://www.bio8.cs.hku.hk/clair3_trio/clair3_trio_models/clair3_trio_models.tar.gz 
tar -zxvf clair3_trio_models.tar.gz -C ./cliar3_trio_models


# run clair3-trio
_INPUT_DIR="[YOUR_INPUT_FOLDER]"            # e.g. ./input
_BAM_C=${_INPUT_DIR}/input_child.bam        # chnage your child's bam file name here
_BAM_P1=${_INPUT_DIR}/input_parent1.bam     # chnage your parenet1's bam file name here
_BAM_P2=${_INPUT_DIR}/input_parent2.bam     # chnage your parenet2's bam file name here
_SAMPLE_C="[Child sample ID]"               # child sample ID, e.g. HG002
_SAMPLE_P1="[Parent1 sample ID]"            # parent1 sample ID, e.g. HG003
_SAMPLE_P2="[Parent2 sample ID]"            # parent2 sample ID, e.g. HG004
_REF=${_INPUT_DIR}/ref.fa                   # change your reference file name here
_OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"          # e.g. ./output
_THREADS="[MAXIMUM_THREADS]"                # e.g. 8
_MODEL_DIR_C3="[Clair3 MODEL NAME]"         # e.g. ont
_MODEL_DIR_C3T="[Clair3-Trio MODEL NAME]"   # e.g. c3t_hg002_g422

./run_clair3_trio.sh \
  --bam_fn_c=${_BAM_C} \    
  --bam_fn_p1=${_BAM_P1} \
  --bam_fn_p2=${_BAM_P2} \
  --output=${_OUTPUT_DIR} \
  --ref_fn=${_REF} \
  --threads=${_THREADS} \
  --model_path_clair3="${_MODEL_DIR_C3}" \
  --model_path_clair3_trio="${_MODEL_DIR_C3T}" \
  --sample_name_c=${_SAMPLE_C} \
  --sample_name_p1=${_SAMPLE_P1} \
  --sample_name_p2=${_SAMPLE_P2}

```

## Usage

### General Usage

**Caution**:  Use `=value` for optional parameters, e.g. `--bed_fn=fn.bed` instead of `--bed_fn fn.bed`.

```bash
./run_clair3_trio.sh \
  --bam_fn_c=${_BAM_C} \
  --bam_fn_p1=${_BAM_P1} \
  --bam_fn_p2=${_BAM_P2} \
  --output=${_OUTPUT_DIR} \
  --ref_fn=${_REF} \
  --threads=${_THREADS} \
  --model_path_clair3="${_MODEL_DIR_C3}" \
  --model_path_clair3_trio="${_MODEL_DIR_C3T}" \
  --bed_fn=${_INPUT_DIR}/quick_demo.bed \
  --sample_name_c=${_SAMPLE_C} \
  --sample_name_p1=${_SAMPLE_P1} \
  --sample_name_p2=${_SAMPLE_P2}

```

### Options

**Required parameters:**

```bash
  -b, --bam_fn=FILE             BAM file input. The input file must be samtools indexed.
  -f, --ref_fn=FILE             FASTA reference file input. The input file must be samtools indexed.
  -m, --model_path=STR          The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00002, pileup.data-00001-of-00002 pileup.index, full_alignment.data-00000-of-00002, full_alignment.data-00001-of-00002  and full_alignment.index).
  -t, --threads=INT             Max threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The chunks being processed simultaneously is ceil($threads/4)*3. 3 is the overloading factor.
  -p, --platform=STR            Select the sequencing platform of the input. Possible options: {ont,hifi,ilmn}.
  -o, --output=PATH             VCF/GVCF output directory.
```

**Other parameters:**

 **Caution**:  Use `=value` for optional parameters, e.g., `--bed_fn=fn.bed` instead of `--bed_fn fn.bed`

```bash
      --bed_fn=FILE             Call variants only in the provided bed regions.
      --vcf_fn=FILE             Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.
      --ctg_name=STR            The name of the sequence to be processed.
      --sample_name=STR         Define the sample name to be shown in the VCF file.
      --qual=INT                If set, variants with >$qual will be marked PASS, or LowQual otherwise.
      --samtools=STR            Path of samtools, samtools version >= 1.10 is required.
      --python=STR              Path of python, python3 >= 3.6 is required.
      --pypy=STR                Path of pypy3, pypy3 >= 3.6 is required.
      --parallel=STR            Path of parallel, parallel >= 20191122 is required.
      --whatshap=STR            Path of whatshap, whatshap >= 1.0 is required.
      --chunk_size=INT          The size of each chuck for parallel processing, default: 5Mbp.
      --pileup_only             Use the pileup model only when calling, default: disable.
      --print_ref_calls         Show reference calls (0/0) in vcf file, default: disable.
      --include_all_ctgs        Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}, default: disable.
      --gvcf                    Enable GVCF output, default: disable.
      --enable_phasing          Output phased variants using whatshap, default: disable.
      --remove_intermediate_dir Remove intermediate directory, including intermediate phased BAM, pileup and full-alignment results. default: disable.
      --snp_min_af=FLOAT        Minimum SNP AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.08,hifi:0.08,ilmn:0.08.
      --indel_min_af=FLOAT      Minimum INDEL AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.15,hifi:0.08,ilmn:0.08.
      --var_pct_full=FLOAT      EXPERIMENTAL: Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.3.
      --ref_pct_full=FLOAT      EXPERIMENTAL: Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for ilmn and hifi, 0.1 for ont.
      --var_pct_phasing=FLOAT   EXPERIMENTAL: Specify an expected percentage of high quality 0/1 variants used in WhatsHap phasing, default: 0.8 for ont guppy5 and 0.7 for other platforms.
      --pileup_model_prefix=STR EXPERIMENTAL: Model prefix in pileup calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index. default: pileup.
      --fa_model_prefix=STR     EXPERIMENTAL: Model prefix in full-alignment calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index, default: full_alignment.
      --fast_mode               EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: disable.
      --haploid_precise         EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant, default: disable.
      --haploid_sensitive       EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant, default: disable.
      --no_phasing_for_fa       EXPERIMENTAL: Call variants without whatshap phasing in full alignment calling, default: disable.
      --call_snp_only           EXPERIMENTAL: Call candidates pass SNP minimum AF only, ignore Indel candidates, default: disable.
      --enable_long_indel       EXPERIMENTAL: Call long Indel variants(>50 bp), default: disable.
```


----

## Folder Structure and Submodule Descriptions

Submodules in __`clair3/`__ are for variant calling and model training. Submodules in __`preprocess`__ are for data preparation.

*For all the submodules listed below, you can use `-h` or `--help` for available options.*

`clair3/` | Note: submodules under this folder are pypy incompatible, please run using python
---: | ---
`CallVariants` | Call variants using a trained model and tensors of candidate variants.
`CallVarBam` | Call variants using a trained model and a BAM file.
`Train` | Training a model using the `RectifiedAdam` optimizer. We also use the `Lookahead` optimizer to adjust the `RectifiedAdam` parameters dynamically. The initial learning rate is `1e-3` with `0.1` learning rate warm-up. Input a binary containing tensors created by `Tensor2Bin`. 



`preprocess/` | Note: submodules under this folder is Pypy compatible unless specified.
---: | ---
`CheckEnvs`| Check the environment and  validity of the input variables, preprocess the BED input if necessary, `--chunk_size` sets the chuck size to be processed per parallel job. 
`CreateTensorPileup`| Generate variant candidate tensors in pileup format for training or calling. 
`CreateTensorFullAlignment`| Generate variant candidate tensors in phased full-alignment format for training or calling. 
`GetTruth`| Extract the variants from a truth VCF. Input: VCF; Reference FASTA if the VCF contains asterisks in ALT field.
`MergeVcf` | Merge pileup and full-alignment VCF/GVCF.
`RealignReads` | Reads local realignment for Illumina platform.
`SelectCandidates`| Select pileup candidates for full-alignment calling.
`SelectHetSnp` | Select heterozygous SNP candidates for whatshap phasing.
`SelectQual` | Select a quality cutoff using the pileup calling results. Variants below the cutoff are included in phasing and full-alignment calling. 
`SortVcf` | Sort VCF file. 
`SplitExtendBed` | Split BED file regions according to the contig names and extend bed region by 33bp by default for variant calling. 
`UnifyRepresentation` | Representation unification between candidate sites and true variants. 
`MergeBin` | Combine tensor binaries into a single file. 
`CreateTrainingTensor` | Create tensor binaries for pileup or full-alignment training. 
`Tensor2Bin` | Combine the variant and non-variant tensors and convert them to a binary, using `blosc:lz4hc` meta-compressor, the overall training memory is 10~15G (pypy incompatible). 

----

## Training Data

Clair3 trained both its pileup and full-alignment models using four GIAB samples (HG001, HG002, HG004 and HG005), excluded HG003. On ONT, we also trained a model using HG001, 2, 3, and 5, excluded HG004. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

|  Platform   |   Reference   |      Aligner      | Training samples |
| :---------: | :-----------: | :---------------: | :--------------: |
|     ONT     | GRCh38_no_alt |     minimap2      | HG001,2,(3\|4),5 |
| PacBio HiFi | GRCh38_no_alt |       pbmm2       |   HG001,2,4,5    |
|  Illumina   |    GRCh38     | BWA-MEM/NovoAlign |   HG001,2,4,5    |

Please find more details about the training data and links at [Training Data](docs/training_data.md).

----

