# GATK + BWA SLURM Pipeline

End-to-end whole genome/exome variant calling pipeline using **BWA**, **GATK**, and **SLURM**.

## Features
- FastQC (pre & post trimming)
- Adapter trimming with Trim Galore
- BWA-MEM alignment
- Duplicate marking
- Base Quality Score Recalibration (BQSR)
- Variant calling with GATK HaplotypeCaller
- SLURM-ready

## Requirements
- BWA
- Samtools
- Picard
- GATK (>=4.x)
- FastQC
- Trim Galore
- MultiQC
- SLURM

## Usage
Edit paths and sample name in the script:
```bash
sample_name="YOUR_SAMPLE"
input_dir="/path/to/fastq"
result="/path/to/output"


Submit to SLURM:

sbatch gatk_bwa_slurm_pipeline.sh
