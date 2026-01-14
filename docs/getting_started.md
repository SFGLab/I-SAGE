# Getting Started with I-SAGE

This section explains how to set up and run the I-SAGE pipeline, from environment requirements to executing a first analysis.

I-SAGE is designed primarily for **HPC environments** and assumes familiarity with command-line tools and batch systems (e.g., SLURM).

## System Requirements

### Software

- **Nextflow** ≥ 22.x
- **Java** ≥ 11
- **Python** ≥ 3.9
- **Bash / core UNIX utilities**

Nextflow manages workflow execution and logging; Python scripts implement statistical and validation logic.

### Bioinformatics Tools
The following tools must be available in the execution environment (typically via modules or Conda):

- **samtools** — BAM processing, sorting, indexing
- **pysam** — Python bindings for BAM/CRAM access
- **bwa** (or equivalent aligner) — read alignment
- **bedGraph / UCSC utilities** — bedGraph and bigWig handling

### Python Dependencies

The following Python packages are required:

- `numpy`
- `pandas`
- `scipy`
- `matplotlib`

These are typically provided via:
- A Conda environment
- A system-wide Python installation
- An HPC module

All plotting is headless-safe (no GUI required).

> **Note:**  
> I-SAGE assumes these tools are provided by the execution environment (HPC modules, Conda, or container).  
> The pipeline does not install system-level dependencies automatically.

## Repository Structure

After cloning the repository, the main structure is:

```
I-SAGE/
├── workflows/
│   └── iblesse_month2/
│       └── main.nf
├── modules/
│   ├── stats/
│   ├── validation/
│   ├── viz/
│   └── break_calling/
├── configs/
│   └── iblesse.yaml
├── docs/
├── README.md
└── nextflow.config
```

You will typically only modify:
- `configs/iblesse.yaml`
- SLURM submission scripts
- (Optionally) Nextflow profiles

## Input Data Requirements

### FASTQ Files

- Paired-end or single-end iBLESS FASTQ files
- File names must match the pattern defined in `sample_pattern`
- Each sample must correspond to a unique biological condition/replicate

Example:
```yaml
sample_pattern: "*_R1.fastq.gz"
```

### Reference Genome

You must provide:

- A reference genome FASTA
- An index compatible with the aligner used in the pipeline

Example:

```yaml
genome_fasta: "/path/to/hg38_reference.fna"
genome_index: "/path/to/hg38_reference.fna"
```

If EBV contigs are included (e.g., `chrEBV`), EBV-specific analyses can be enabled.

## Configuration File (`iblesse.yaml`)

All pipeline behavior is controlled via a single YAML file.

Key sections include:

- `break_calling`
- `viz`
- `normalization`
- `stats`
- `validation`

A fully working example is provided in `configs/iblesse.yaml`. Detailed explanations are provided in the **Configuration Guide** section of the documentation.

## Running the Pipeline

### Basic Command

From the repository root:

```bash
nextflow run workflows/iblesse_month2/main.nf \
  -profile eden_local \
  -params-file configs/iblesse.yaml
```

This will:

- Read FASTQ files
- Execute all enabled pipeline stages
- Write results to the directory specified by `outdir`

### Running on an HPC Cluster (SLURM)

I-SAGE is commonly run inside a SLURM allocation.

Typical workflow:

1. Write a SLURM submission script
2. Allocate resources
3. Run Nextflow inside the job

Example (simplified):

```bash
sbatch run_isage.slurm
```

Within the SLURM script, Nextflow is executed normally. Nextflow parallelizes tasks within the allocated resources.

## Execution Outputs

During and after execution, I-SAGE produces:

- **Nextflow report** (`report_*.html`)
- **Execution trace** (`trace_*.txt`)
- **Timeline view** (`timeline_*.html`)
- Structured output directories:
  - `viz/`
  - `stats/`
  - `validation/`

These files are critical for debugging and reproducibility.

## First-Time Run Checklist

Before running a full analysis, verify:

- FASTQ paths are correct
- Reference genome paths exist
- Output and log directories are writable
- The selected Nextflow profile matches your environment
- Bin sizes and contrasts are intentional (not defaults by accident)

## Common First-Run Pitfalls

- Using `-resume` unintentionally (skips updated steps)
- Misnamed sample IDs in contrasts
- Missing EBV contigs while EBV analysis is enabled
- Insufficient disk space for intermediate files
- Running bin-size sweeps without accounting for increased runtime

## Next Steps

After completing a successful run:

- Review output structure
- Inspect differential statistics and plots
- Consult the **Configuration Guide** to fine-tune parameters
- Proceed to **Pipeline Modules** for a deeper understanding of each stage

---

**Next:** See [Configuration Guide](configuration.md) for detailed parameter explanations.
