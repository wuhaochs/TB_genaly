# TB_genaly
A Multicenter Evaluated Whole Genome Sequencing Pipeline for Tuberculosis Drug Resistance Assessment and Transmission Analysis
## 1. Project Introduction
This repository contains a professional, fully automated whole-genome sequencing (WGS) analysis pipeline for **Mycobacterium tuberculosis (MTB)**. It consists of 3 core Python scripts, which implement a complete analytical workflow from raw sequencing data preprocessing to advanced genomic analysis, including variant calling, drug resistance annotation, species identification, lineage typing, genetic distance calculation, transmission cluster analysis, and phylogenetic tree construction.

This pipeline is developed for clinical and research scenarios of tuberculosis genomic surveillance, drug resistance detection, and transmission source tracking, and is fully compliant with the WHO guidelines for tuberculosis drug resistance detection and genomic epidemiology analysis.

| Pipeline Information | Details |
|----------------------|---------|
| Main Pipeline Version | v1.6.6 |
| Developer | wsh |
| Supported Platform | Linux (CentOS 7+/Ubuntu 18.04+, Linux-only) |
| Compatible Input Types | Illumina paired-end short-read data, ONT long-read data (raw fast5 signal included), assembled FASTA genome, aligned BAM file, variant VCF file |
| Core Features | One-stop automated analysis, breakpoint restart support, multi-input format compatibility, clinical drug resistance interpretation and medication guidance, visualization report generation |

---
## 2. Table of Contents
1. [Project Introduction](sslocal://flow/file_open?url=%231-project-introduction&flow_extra=eyJsaW5rX3R5cGUiOiJjb2RlX2ludGVycHJldGVyIn0=)
2. [Environment & Dependencies](sslocal://flow/file_open?url=%232-environment--dependencies&flow_extra=eyJsaW5rX3R5cGUiOiJjb2RlX2ludGVycHJldGVyIn0=)
3. [Core Scripts - Detailed Description & Usage](sslocal://flow/file_open?url=%233-core-scripts---detailed-description--usage&flow_extra=eyJsaW5rX3R5cGUiOiJjb2RlX2ludGVycHJldGVyIn0=)
4. [Standard End-to-End Analysis Workflow](sslocal://flow/file_open?url=%234-standard-end-to-end-analysis-workflow&flow_extra=eyJsaW5rX3R5cGUiOiJjb2RlX2ludGVycHJldGVyIn0=)
5. [Output File Interpretation](sslocal://flow/file_open?url=%235-output-file-interpretation&flow_extra=eyJsaW5rX3R5cGUiOiJjb2RlX2ludGVycHJldGVyIn0=)
6. [Important Notes & Precautions](sslocal://flow/file_open?url=%236-important-notes--precautions&flow_extra=eyJsaW5rX3R5cGUiOiJjb2RlX2ludGVycHJldGVyIn0=)
7. [Version Changelog](sslocal://flow/file_open?url=%237-version-changelog&flow_extra=eyJsaW5rX3R5cGUiOiJjb2RlX2ludGVycHJldGVyIn0=)
8. [Troubleshooting & FAQ](sslocal://flow/file_open?url=%238-troubleshooting--faq&flow_extra=eyJsaW5rX3R5cGUiOiJjb2RlX2ludGVycHJldGVyIn0=)

---
## 2. Environment & Dependencies
### 2.1 Basic Runtime Environment
- **Operating System**: Linux (CentOS 7+/Ubuntu 18.04+, the pipeline only supports Linux system)
- **Python Interpreter**: Python 3.7+ (Conda environment is highly recommended, the default environment for the pipeline is `TB_pip2`)
- **Hardware Minimum Requirements**:
  - Single sample short-read analysis: 16GB RAM, 10 CPU threads
  - Multi-sample phylogenetic and transmission analysis: 32GB RAM, 20 CPU threads
  - ONT fast5 basecalling: NVIDIA GPU (CPU mode is supported but with significantly slower speed)
  - Full microbial database species identification: 64GB RAM, sufficient disk space for database storage

### 2.2 Python Dependencies
Install the required Python packages with the following command:
```bash
pip install pandas==1.5.3 biopython pysam pymysql argparse regex
```

### 2.3 Required Bioinformatics Software
All the following software must be pre-installed and added to the system `PATH`, or updated to the absolute path of the software in the corresponding script.

| Functional Category | Required Software |
|---------------------|-------------------|
| Data Preprocessing & QC | fastp, fastqc, seqkit, unzip, gzip |
| Sequence Alignment | bwa, minimap2, samtools (v1.15+) |
| Variant Calling & Processing | freebayes, bcftools, vcftools, vt, snpEff, snippy |
| Coverage & Depth Calculation | mosdepth, qualimap, bedtools |
| Species Identification | kraken2, bracken, KronaTools |
| MTB Lineage Typing | TB-profiler |
| Phylogenetic & Genetic Distance Analysis | iqtree2, grapetree, snp-dists, fastANI |
| Long-read Data Processing | guppy_basecaller |
| Image Visualization | ImageMagick (convert tool) |
| Report Generation | R (with required packages for report rendering) |

### 2.4 Reference Database & Required Files
The pipeline uses absolute paths for all reference files by default. You must prepare the following files and verify the path consistency before running the pipeline.

| File Type | Default Absolute Path | Description |
|-----------|------------------------|-------------|
| MTB Reference Genome | `/data/deploy/TB_soft/ref/TB/ref.fa` | MTB H37Rv reference genome, with matching .fai index, bwa index, and gff3 annotation file |
| MTB Reference GFF | `/data/deploy/TB_soft/ref/TB/ref.gff` | GFF3 annotation file of the MTB reference genome |
| snpEff Config File | `/data/deploy/TB_soft/ref/TB/snpeff.config` | Configuration file for snpEff variant annotation |
| Drug Resistance Annotation Database | `/data/deploy/TB_soft/ref/TB/` | Includes WHO-standard drug resistance locus table, MIC value table, clinical medication guidance file, and MySQL drug resistance database `tb_dr_database` |
| Blacklist Region BED | `/data/deploy/TB_soft/ref/TB/black_list.bed` | BED file of hypervariable regions, recombinant regions, and PE/PPE gene regions excluded from variant calling and phylogenetic analysis |
| Virulence Gene Annotation | `/data/deploy/TB_soft/ref/TB/VB_anno.bed` | BED file of MTB virulence gene annotation |
| Kraken2 Taxonomic Database | `/data/deploy/TB_soft/Database/` | Two pre-built kraken2 databases: `TB_library` (MTB-specific library, fast analysis) and `GT_DB` (full microbial library, comprehensive species identification) |
| Taxonomic Annotation File | `/data/deploy/TB_soft/other_soft/taxa_info_0717.txt` | Taxonomic information file for species identification result annotation |

---
## 3. Core Scripts - Detailed Description & Usage
### 3.1 Main Pipeline Script: `Tcu_super_analysis.py`
This is the core script of the pipeline, which implements one-stop full workflow analysis of MTB WGS data, from raw data input to final result report generation.

#### 3.1.1 Core Functions
1. **Multi-format Input Support**: Compatible with Illumina paired-end fastq, ONT single-end fastq, raw fast5 directory, assembled FASTA genome, aligned BAM file, and variant VCF file
2. **Data Preprocessing**: Automated quality filtering with fastp, fastqc quality report generation and visualization, ONT fast5 basecalling, and FASTA to fastq conversion
3. **Reference Genome Preparation**: Automated index building for reference genome, snpEff database construction, and required file preparation
4. **Sequence Alignment & Variant Calling**: Support bwa for short-read alignment, minimap2 for long-read alignment; support freebayes for variant calling, and optional clockwork professional MTB variant calling workflow
5. **Variant Filtering & Annotation**: Multi-dimensional filtering of variants (depth, quality, MAF), snpEff functional annotation, and correction of multi-mutation in the same codon
6. **Drug Resistance Analysis**: WHO-standard drug resistance locus matching, MIC value annotation, drug resistance phenotype classification (Sensitive/RR-TB/HR-TB/MDR-TB/Pre-XDR-TB/XDR-TB), and clinical medication guidance generation
7. **Species & Lineage Analysis**: kraken2 + bracken species identification, Krona interactive visualization, and TB-profiler MTB lineage typing
8. **Result Summary & Report Generation**: Automated quality control summary, multi-sample result merging, HTML visualization report generation, and one-click packaging of core results
9. **Checkpoint Restart**: Support breakpoint restart, automatically skip completed steps when re-running the same command

#### 3.1.2 Full Command Line Parameters
| Parameter | Short Flag | Type | Default Value | Required | Detailed Description |
|-----------|------------|------|---------------|----------|----------------------|
| `--list` | `-l` | string | - | ✅ Yes | Tab-delimited input sample list file, see section 3.1.3 for format details |
| `--output` | `-o` | string | - | ✅ Yes | Output directory, will be created automatically if it does not exist |
| `--thread` | `-t` | integer | 10 | ❌ No | Number of CPU threads used for the entire analysis workflow |
| `--transcutoff` | `-u` | string | `5,12` | ❌ No | SNP number thresholds for transmission cluster analysis, comma-separated |
| `--vcf_dpt` | `-v` | integer | 5 | ❌ No | Minimum sequencing depth threshold for variant calling and filtering |
| `--db` | `-d` | string | `TB_library` | ❌ No | Kraken2 database selection: `TB_library` (MTB-specific library, fast speed) / `GT_DB` (full microbial library, slow speed) |
| `--mode` | `-m` | string | `hac` | ❌ No | ONT fast5 basecalling mode: `hac` (high accuracy mode) / `fast` (fast mode) |
| `--maf` | `-mf` | float | 0.25 | ❌ No | Minor allele frequency (MAF) filtering threshold for variant calling |
| `--clockwork` | `-ck` | string | `False` | ❌ No | Whether to enable the clockwork workflow for professional MTB variant calling, supported values: `True/False/1/0/F/f` |

#### 3.1.3 Input Sample List Format
The file specified by `--list` must be a **tab-delimited** plain text file, with one sample per line. The following formats are supported:
```
# 1. Illumina paired-end short-read data (3 columns: sample_name read1_path read2_path)
Sample_001  /data/raw_data/Sample_001_R1.fastq.gz  /data/raw_data/Sample_001_R2.fastq.gz
Sample_002  /data/raw_data/Sample_002_R1.fastq.gz  /data/raw_data/Sample_002_R2.fastq.gz

# 2. ONT long-read single-end data (2 columns: sample_name fastq_path)
ONT_Sample_001  /data/ont_data/ONT_Sample_001.fastq.gz

# 3. ONT raw fast5 directory (2 columns: sample_name fast5_directory_path)
ONT_Raw_Sample_001  /data/ont_fast5/ONT_Raw_Sample_001/

# 4. Assembled FASTA genome (2 columns: sample_name fasta_path)
Assembly_Sample_001  /data/assembly/Sample_001_genome.fasta

# 5. Aligned BAM file (2 columns: sample_name bam_path)
BAM_Sample_001  /data/aligned/Sample_001_sorted.bam

# 6. Variant VCF file (2 columns: sample_name vcf_path)
VCF_Sample_001  /data/variant/Sample_001_variants.vcf
```

#### 3.1.4 Usage Examples
```bash
# Basic usage: Illumina short-read data, MTB-specific database, 10 threads
python Tcu_super_analysis.py \
  --list sample_list.txt \
  --output ./MTB_analysis_result \
  --thread 10

# Enable clockwork workflow for high-precision variant calling, 20 threads
python Tcu_super_analysis.py \
  -l sample_list.txt \
  -o ./MTB_clockwork_analysis \
  -t 20 \
  -ck True

# ONT fast5 raw data analysis, fast basecalling mode, minimum depth 10
python Tcu_super_analysis.py \
  -l ont_sample_list.txt \
  -o ./ONT_MTB_analysis \
  -t 16 \
  -m fast \
  -v 10

# Full microbial database for species identification, custom MAF threshold 0.1
python Tcu_super_analysis.py \
  -l sample_list.txt \
  -o ./MTB_full_species_analysis \
  -t 20 \
  -d GT_DB \
  -mf 0.1
```

#### 3.1.5 Output Directory Structure
```
<output_directory>/
├── 0.QC/                          # Alignment quality control results
│   ├── *mapping_summary.tsv       # Mapping statistics for each sample
│   ├── *_quaimap/                 # Qualimap BAM QC report for each sample
│   └── *_ok                        # QC completion flag for each sample
├── 1.snp_calling/                 # Variant calling results
│   ├── <sample_name>/             # Separate directory for each sample
│   │   ├── snps.bam               # Sorted and deduplicated BAM file
│   │   ├── snps.bam.bai           # BAM index file
│   │   ├── snps.raw.vcf           # Raw unfiltered variant VCF file
│   │   ├── snps.filt.vcf          # Filtered variant VCF file
│   │   ├── snps.vcf.gz             # Compressed and indexed final VCF file
│   │   ├── snps.consensus.fa       # Consensus genome sequence
│   │   ├── Chin_snpdr.tsv          # Chinese version drug resistance annotation table
│   │   ├── snpdr.tsv               # Standard drug resistance annotation table
│   │   ├── DrugGuide.tsv           # Clinical medication guidance
│   │   ├── sam_dr.tsv               # Drug resistance phenotype result
│   │   ├── index_<sample>.html      # Interactive IGV visualization HTML
│   │   └── map_OK                   # Analysis completion flag
│   ├── all_sample.full.aln         # Full genome multi-sequence alignment of all samples
│   ├── all_sample.vcf               # Core genome variant VCF of all samples
│   └── all_sample_dr.tsv            # Merged drug resistance results of all samples
├── 2.Tree/                          # Phylogenetic tree construction directory
├── 3.Tb_profiler/                   # TB-profiler lineage typing results
│   ├── results/                     # Detailed typing results for each sample
│   ├── tbprofiler.txt               # Merged lineage typing results of all samples
│   └── tbpro.tsv                     # Simplified lineage result table
├── 4.kraken_taxonomic/              # Species identification results
│   ├── *.kraken2.txt                 # Kraken2 report for each sample
│   ├── *.bracken2.txt                # Bracken abundance estimation result
│   ├── *.krona.html                  # Interactive Krona visualization HTML
│   ├── *.list2.txt                   # Species abundance table for each sample
│   └── *_ok                          # Analysis completion flag
├── 5.virulence_genes/                # Virulence gene annotation results
├── 6.Snpit/                          # SNP-IT lineage typing results
├── 7.snpcluster/                     # Transmission cluster analysis results
├── fq_file/                          # Post-QC clean fastq data
│   ├── *_clean_1.gz / *_clean_2.gz  # Clean paired-end fastq
│   ├── *.fastp.json / *.fastp.log   # fastp QC report and log
│   ├── result_base.txt               # Merged QC statistics of all samples
│   └── *_OK                          # QC completion flag
├── ref/                               # Reference genome index and annotation files
├── main_results/                      # Summarized core results
│   ├── *.consensus.fasta             # Consensus genome of each sample
│   ├── *_snp.vcf                      # Final variant VCF of each sample
│   ├── *_Chin_snpdr.tsv               # Drug resistance annotation of each sample
│   ├── *_Spe.tsv                       # Species identification result of each sample
│   ├── lineage.tsv                     # Merged lineage typing result
│   └── *mapping_summary.tsv            # Mapping statistics
├── main_results.zip                    # One-click compressed package of core results
├── CombineAll.tsv                      # Comprehensive analysis summary table of all samples
├── sample_results.txt                  # Overall run statistics
├── *.html                               # Final visualization analysis report
├── trim_fqlist                          # List of clean fastq files
├── Samplelist.txt                       # List of all successfully analyzed samples
└── task.log / *.log                     # Run log files
```

---
### 3.2 Phylogenetic Tree Construction Script: `Treebuild.py`
This script implements automated core genome alignment and phylogenetic tree construction for MTB genome sequences, optimized for MTB transmission source tracking and evolutionary analysis.

#### 3.2.1 Core Functions
1. Automated core genome alignment with snippy, support for custom exclusion of blacklist regions
2. Support multiple tree construction methods: Maximum Likelihood (ML, iqtree2), Neighbor-Joining (NJ), and grapetree
3. Support bootstrap validation for phylogenetic tree reliability assessment
4. Support merging user-provided reference sequences and target samples for joint tree construction
5. Automatic metadata annotation for phylogenetic tree visualization

#### 3.2.2 Full Command Line Parameters
| Parameter | Short Flag | Type | Default Value | Required | Detailed Description |
|-----------|------------|------|---------------|----------|----------------------|
| `--inputlist` | `-i1` | string | False | ❌ No | Tab-delimited pre-defined sample list, format: `SampleName\tfasta_path` |
| `--inf2` | `-i2` | string | False | ❌ No | Fasta sequence file of a single target sample |
| `--inf3` | `-i3` | string | 0 | ❌ No | Path to a directory containing multiple fasta files for batch analysis |
| `--output` | `-o` | string | `Treebuild` | ❌ No | Output directory for tree construction results |
| `--boots` | `-b` | integer | 100 | ❌ No | Number of bootstrap replicates for phylogenetic tree validation |
| `--Samname` | `-s` | string | `focus_Sam` | ❌ No | Sample name for the single target sample (corresponding to --inf2) |
| `--Lin` | `-l` | integer | 0 | ❌ No | Lineage filter parameter |
| `--method` | `-mt` | string | `ml` | ❌ No | Tree construction method: `ml` (Maximum Likelihood, iqtree2) / `nj` (Neighbor-Joining) / `grapetree` |
| `--exbed` | `-eb` | string | `/data/deploy/TB_soft/ref/TB/black_list.bed` | ❌ No | BED file of regions excluded from core genome alignment and tree construction |

#### 3.2.3 Usage Examples
```bash
# Maximum Likelihood tree construction with 1000 bootstrap replicates
python Treebuild.py \
  --inputlist sample_fasta_list.txt \
  --output ./ML_tree_result \
  --method ml \
  --boots 1000

# Neighbor-Joining tree construction with custom excluded regions
python Treebuild.py \
  -i1 sample_fasta_list.txt \
  -o ./NJ_tree_result \
  -mt nj \
  -eb ./custom_excluded_regions.bed

# Joint tree construction for a single sample with reference database
python Treebuild.py \
  -i1 MTB_reference_list.txt \
  -i2 ./my_sample.fasta \
  -s my_sample \
  -o ./sample_reference_tree \
  -b 500

# Grapetree rapid tree construction for large sample size
python Treebuild.py \
  -i1 large_sample_fasta_list.txt \
  -o ./grapetree_result \
  -mt grapetree
```

#### 3.2.4 Core Output Files
| File Name | Description |
|-----------|-------------|
| `core.aln` | Core genome multi-sequence alignment result |
| `core.aln.contree` | Final consensus phylogenetic tree in Newick format (compatible with FigTree, iTOL, MEGA, etc.) |
| `core.aln.treefile` | Original ML tree file from iqtree2 |
| `runme.sh` | Snippy core genome alignment script |
| `task.log` | Detailed run log of the tree construction process |
| `tbprofiler.txt` | Metadata for sample annotation |

---
### 3.3 Pairwise SNP Distance & Transmission Cluster Analysis Script: `pairwise.py`
This script implements pairwise SNP distance calculation, genetic distance estimation, and transmission cluster analysis for MTB samples, which is the core tool for MTB genomic epidemiology and transmission chain tracking.

#### 3.3.1 Core Functions
1. Automated core genome alignment and pairwise SNP distance calculation between all samples
2. Support multiple genetic distance calculation methods: SNP distance, TN93 model, and ANI (Average Nucleotide Identity)
3. Transmission cluster analysis based on custom SNP threshold (WHO recommended 5/12 SNP for recent transmission)
4. Binned statistics of SNP distance distribution
5. Output of distance matrix and detailed clustering results for downstream visualization and analysis

#### 3.3.2 Full Command Line Parameters
| Parameter | Short Flag | Type | Default Value | Required | Detailed Description |
|-----------|------------|------|---------------|----------|----------------------|
| `--inputlist` | `-l` | string | - | ✅ Yes | Tab-delimited sample list, format: `SampleName\tfasta_path` |
| `--output` | `-o` | string | - | ✅ Yes | Output directory for analysis results |
| `--transcutoff` | `-u` | string | `5` | ❌ No | SNP threshold for transmission clustering (samples with SNP distance ≤ this value are grouped into the same cluster, WHO recommended 5/12) |
| `--snpdis` | `-s` | string | `10,20,50,100,200,500,1000` | ❌ No | Comma-separated intervals for SNP distance binned statistics |
| `--method` | `-m` | string | `TN93` | ❌ No | Genetic distance calculation method: `TN93` / `SNP` / `ANI` |
| `--ref` | `-r` | string | `/data/deploy/TB_soft/ref/TB/ref.fa` | ❌ No | MTB reference genome path for core genome alignment |

#### 3.3.3 Usage Examples
```bash
# Basic transmission cluster analysis with 5 SNP threshold
python pairwise.py \
  --inputlist sample_fasta_list.txt \
  --output ./snp_distance_result \
  --transcutoff 5

# 12 SNP threshold for transmission clustering, custom binned statistics
python pairwise.py \
  -l sample_fasta_list.txt \
  -o ./transmission_cluster_result \
  -u 12 \
  -s 5,12,50,100,200

# ANI method for genome similarity calculation
python pairwise.py \
  -l sample_fasta_list.txt \
  -o ./ANI_analysis_result \
  -m ANI

# SNP distance calculation with custom reference genome
python pairwise.py \
  -l sample_fasta_list.txt \
  -o ./custom_ref_result \
  -r ./custom_MTB_ref.fa
```

#### 3.3.4 Core Output Files
| File Name | Description |
|-----------|-------------|
| `dis.mat.txt` | Symmetric pairwise SNP distance matrix of all samples |
| `t_dis.tsv` | Detailed pairwise SNP difference count for each sample pair |
| `dis_bin.tsv` | Binned statistics of SNP distance distribution |
| `Cluster.tsv` | Transmission clustering result, including cluster ID, sample count, sample list, max/min/mean SNP distance per cluster |
| `Gdis.txt` / `Gdis_core.txt` | Genetic distance calculation result for the selected method |
| `core.aln` | Core genome multi-sequence alignment result |
| `snp.log` | Detailed run log of the analysis process |

---
## 4. Standard End-to-End Analysis Workflow
### Step 1: Prepare Input Data and Sample List
- Prepare your sequencing data (fastq/fast5/fasta/bam/vcf)
- Create a tab-delimited sample list file according to the format in section 3.1.3
- Verify all file paths are correct and accessible

### Step 2: Run the Main Analysis Pipeline
Use `Tcu_super_analysis.py` to process raw data, complete quality control, variant calling, drug resistance annotation, species and lineage identification. This step will generate consensus genome sequences for all samples.

### Step 3: Transmission Cluster Analysis
Extract the consensus fasta sequences from the main pipeline result, create a sample list, and use `pairwise.py` to calculate pairwise SNP distances and perform transmission cluster analysis to identify potential recent transmission chains.

### Step 4: Phylogenetic Tree Construction
Use `Treebuild.py` to construct a phylogenetic tree based on the core genome alignment of all samples, combined with clustering results, drug resistance phenotypes, and lineage information to visualize the evolutionary and transmission relationships.

### Step 5: Result Interpretation and Report Generation
The main pipeline will automatically generate a comprehensive HTML report and a summary table of all results. You can further visualize and analyze the results based on the output files.

---
## 5. Output File Interpretation
### Key Result Files
1. **`CombineAll.tsv`**: Comprehensive summary table of all samples, including sequencing platform, data volume, QC result, species, lineage, drug resistance phenotype, and detailed drug resistance mutation information for each sample.
2. **`main_results.zip`**: Compressed package of all core results, including consensus genome, variant VCF, drug resistance annotation, species identification, and lineage typing results for each sample.
3. **`1.snp_calling/<sample_name>/Chin_snpdr.tsv`**: Detailed drug resistance annotation table for each sample, including variant position, mutation frequency, gene, mutation effect, drug name, resistance correlation, MIC value, etc.
4. **`1.snp_calling/<sample_name>/DrugGuide.tsv`**: Clinical medication guidance for each sample, based on the drug resistance phenotype and WHO treatment guidelines.
5. **`4.kraken_taxonomic/<sample_name>.krona.html`**: Interactive Krona visualization of species identification results, which can be opened in any browser.
6. **`pairwise.py/Cluster.tsv`**: Transmission clustering result, which can be used to identify samples with close genetic distance and potential recent transmission.
7. **`Treebuild.py/core.aln.contree`**: Phylogenetic tree file in Newick format, which can be imported into FigTree, iTOL, MEGA and other tools for visualization and annotation.

---
## 6. Important Notes & Precautions
1. **Path Adaptation**: All software and database paths in the scripts use absolute paths by default. You must verify the actual deployment paths of all software and reference files before running the pipeline, and update the corresponding paths in the scripts if they are inconsistent.
2. **Input Specification**: The sample list file must be **tab-delimited** (spaces are not allowed as separators). Sample names must not contain special characters such as spaces, `/`, `\`, `@`, `#`, etc., to avoid workflow errors.
3. **Permission Requirements**: You must have read permission for all input files, write permission for the output directory, and execution permission for all dependent bioinformatics software.
4. **Clinical Disclaimer**: The drug resistance annotation, MIC value interpretation, and clinical medication guidance output by this pipeline are for research reference only, and shall NOT be used as the sole basis for clinical medication prescription. All clinical decisions must be made by professional clinicians combined with phenotypic drug susceptibility testing results.
5. **Checkpoint Restart**: The pipeline supports breakpoint restart. If the run is interrupted due to server failure or other reasons, you can re-execute the same command, and the pipeline will automatically skip the completed steps and continue the unfinished analysis.
6. **Database Update**: The drug resistance database and species identification database need to be updated regularly to ensure the accuracy of the analysis results.
7. **Resource Usage**: For large sample size analysis, please ensure sufficient CPU threads, memory, and disk space to avoid pipeline failure due to insufficient resources.

---
## 7. License
This pipeline is developed for academic and research use only. For commercial use, please contact the developer for authorization.
