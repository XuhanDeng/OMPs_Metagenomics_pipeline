# OMPs Metagenomics Pipeline

A comprehensive Snakemake workflow for metagenomics analysis, specifically designed for Organometallic Pollutants (OMPs) research. This pipeline processes raw sequencing data through quality control, assembly, gene prediction, clustering, and functional annotation.

## Overview

This pipeline implements a complete metagenomics workflow including:

- **Quality Control**: Read trimming and filtering with fastp
- **Assembly**: Metagenomic assembly with SPAdes
- **Scaffold Processing**: Length filtering and reformatting with anvi'o
- **Gene Prediction**: ORF identification with Pyrodigal
- **Gene Clustering**: Non-redundant gene catalog creation with MMSeqs2
- **Taxonomic Annotation**: Species assignment with GTDB and Kraken2/Bracken
- **Functional Annotation**: Protein function assignment with EggNOG
- **Abundance Quantification**: Gene abundance calculation with CoverM
- **Taxonomic Profiling**: Read-based taxonomic classification with Kraken2/Bracken/Taxpasta

## Pipeline Workflow

### Main Assembly-Based Workflow
```
Raw Reads (FASTQ)
    ↓
fastp (Quality Control)
    ↓
SPAdes (Metagenomic Assembly)
    ↓
anvi'o (Scaffold Reformatting)
    ↓
Pyrodigal (ORF Prediction)
    ↓
MMSeqs2 (Gene Clustering)
    ↓
Sequence Extraction
    ↓
┌─────────────────────┬─────────────────────┐
│  GTDB (Taxonomy)    │  EggNOG (Function)  │
└─────────────────────┴─────────────────────┘
    ↓
CoverM (Abundance Quantification)
```

### Read-Based Taxonomic Profiling (Snakefile_kraken.smk)
```
Quality-Controlled Reads (fastp)
    ↓
Kraken2 (Taxonomic Classification)
    ↓
Bracken (Abundance Estimation at Multiple Levels)
    ├── Phylum
    ├── Class
    ├── Order
    ├── Family
    ├── Genus
    └── Species
    ↓
┌─────────────────────────┬──────────────────────────┐
│  Merge Reports (Pandas) │  Taxpasta (with Lineage) │
│  97_braken_report/      │  96_taxpasta/            │
└─────────────────────────┴──────────────────────────┘
```

## Features

- 🔧 **SLURM Cluster Compatible**: Optimized resource allocation for HPC environments
- 🐍 **Conda Integration**: Isolated environments for reproducible analysis
- 📊 **Sample Grouping**: Process samples in user-defined groups (lab/osen)
- ⚙️ **Configurable**: All parameters centralized in `config.yaml`
- 📁 **Organized Output**: Structured directory organization
- 🔄 **Scalable**: Handles multiple samples efficiently

## Requirements

### Software Dependencies
- [Snakemake](https://snakemake.readthedocs.io/) (≥7.0)
- [Conda/Mamba](https://docs.conda.io/en/latest/)
- SLURM workload manager (for cluster execution)

### Hardware Requirements
- **Memory**: Minimum 256GB RAM for SPAdes assembly
- **Storage**: ~500GB free space per sample
- **CPUs**: Multi-core system recommended (8+ cores)

## Installation

### 1. Clone Repository
```bash
git clone https://github.com/XuhanDeng/OMPs_Metagenomics_pipeline.git
cd OMPs_Metagenomics_pipeline
```

### 2. Install Snakemake
```bash
# Create and activate conda environment
conda create -n snakemake snakemake mamba -c conda-forge -c bioconda
conda activate snakemake
```

### 3. Setup Databases
Download required databases and update paths in `config/config.yaml`:

```bash
# GTDB database
wget https://data.gtdb.ecogenomic.org/releases/latest/

# EggNOG database
download_eggnog_data.py -y

# Kraken2 database (for taxonomic profiling)
kraken2-build --standard --db /path/to/kraken2_db

# Build Bracken database
bracken-build -d /path/to/kraken2_db -t 32 -k 35 -l 150

# NCBI Taxonomy database (for Taxpasta)
cd /path/to/taxa_DB
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
md5sum --check taxdump.tar.gz.md5
mkdir taxdump
tar -C taxdump -xzf taxdump.tar.gz
```

Update paths in `config/config.yaml`:
```yaml
databases:
  gtdb_diamond: "/scratch/qou/database/GTDB_214.1_diamond/protein_faa_reps/GTDB214_dimond.dmnd"
  gtdb_database: "/scratch/qou/database/GTDB_214.1_diamond/GTDB_database"
  eggnog_data_dir: "/scratch/qou/database/eggnog"
  kraken2_db: "/scratch/qou/database/kraken2"
  taxa_db: "/scratch/qou/database/taxa_DB/taxdump"
```

## Configuration

### Sample Setup
Edit `config/config.yaml` to specify your samples:

```yaml
# Sample information
samples:
  - sample1
  - sample2
  - sample3

# Sample groups for different analysis
sample_groups:
  lab:
    - sample1
    - sample2
  osen:
    - sample3
```

### Directory Structure
Organize your input data as follows:
```
result/0_rawdata/
├── sample1.R1.fq.gz
├── sample1.R2.fq.gz
├── sample2.R1.fq.gz
├── sample2.R2.fq.gz
├── sample3.R1.fq.gz
└── sample3.R2.fq.gz
```

### Resource Configuration
Adjust computational resources in `config/config.yaml`:

```yaml
# Common resources
account: "research-ceg-wm"
regular_partition: "compute"
regular_memory: "3900"  # MB per CPU

# Tool-specific parameters
fastp:
  threads: 4
  runtime: 120  # minutes

spades:
  threads: 32
  memory: 256  # GB total
  runtime: 1200  # minutes

kraken2:
  threads: 40
  confidence: 0.1
  runtime: 1200  # minutes

bracken:
  threads: 30
  kmer_length: 35
  read_length: 150
  runtime: 600  # minutes
```

## Usage

### Main Workflow (Assembly-based)

#### Local Execution (Small datasets)
```bash
# Dry run to check workflow
snakemake --dry-run

# Execute with conda environments
snakemake --use-conda --cores 8
```

#### SLURM Cluster Execution
```bash
# Create SLURM profile (first time only)
cookiecutter https://github.com/Snakemake-Profiles/slurm.git

# Execute on cluster
snakemake --profile slurm --use-conda
```

### Kraken2/Bracken Taxonomic Profiling

```bash
# Dry run
snakemake -s Snakefile_kraken.smk --dry-run

# Local execution
snakemake -s Snakefile_kraken.smk --use-conda --cores 40

# SLURM cluster execution
snakemake -s Snakefile_kraken.smk --profile slurm --use-conda
```

### Target Specific Outputs
```bash
# Run only quality control
snakemake --use-conda --cores 8 \
  1_fastp/sample1/sample1_1P.fq.gz

# Run only assembly for specific sample  
snakemake --use-conda --cores 8 \
  2_spades_result/sample1/contigs.fasta

# Generate all final outputs
snakemake --use-conda --cores 8 all
```

## Output Structure

### Main Assembly-Based Pipeline
```
result/
├── 1_fastp/                    # Quality-controlled reads
│   └── {sample}/
│       ├── {sample}_1P.fq.gz
│       ├── {sample}_2P.fq.gz
│       ├── {sample}.fastp.html
│       └── {sample}.fastp.json
├── 2_spades_result/            # Assembly results
│   └── {sample}/
│       └── contigs.fasta
├── 3_reformated_scaffolds/     # Filtered scaffolds
│   └── {sample}/
│       ├── {sample}_reformated.fa
│       └── {sample}_renaming.tsv
├── 4_pyrodigal/                # Gene predictions
│   └── {sample}/
│       ├── gff/{sample}.gff
│       ├── amino/{sample}.faa
│       └── nucl/{sample}.nucl.ffn
├── 5_mmseqs2/                  # Gene clustering results
│   └── {group}/
│       ├── {group}_cluster_rep_seq.fasta
│       ├── unigene_id.txt
│       └── unigene_protein.fasta
├── 6_coverm/                   # Abundance quantification
│   └── {group}/
│       └── coverm_{sample}.TPM.tsv
├── 7_GTDB/                     # Taxonomic annotation
│   └── {group}/
│       └── annotation.txt
└── 8_functional_annotation/    # Functional annotation
    └── {group}_eggnogresult.emapper.annotations
```

### Kraken2/Bracken Taxonomic Profiling (Snakefile_kraken.smk)
```
result/
├── 99_kraken2/                 # Kraken2 and Bracken outputs
│   └── {sample}/
│       ├── {sample}_kraken2.report
│       ├── {sample}_kraken2.out
│       ├── {sample}_bracken_phylum.report
│       ├── {sample}_bracken_classes.report
│       ├── {sample}_bracken_orders.report
│       ├── {sample}_bracken_family.report
│       ├── {sample}_bracken_genus.report
│       └── {sample}_bracken_species.report
├── 98_bracken_kraken/          # Organized reports by taxonomy level
│   ├── krakenstykle/
│   │   ├── P/                  # Phylum level
│   │   ├── C/                  # Class level
│   │   ├── O/                  # Order level
│   │   ├── F/                  # Family level
│   │   ├── G/                  # Genus level
│   │   └── S/                  # Species level
│   └── brackenstykle/
│       └── [same structure]
├── 97_braken_report/           # Merged reports across samples
│   ├── phylum_rel_abund.csv
│   ├── phylum_read_numbers.csv
│   ├── classes_rel_abund.csv
│   ├── classes_read_numbers.csv
│   ├── orders_rel_abund.csv
│   ├── orders_read_numbers.csv
│   ├── family_rel_abund.csv
│   ├── family_read_numbers.csv
│   ├── genus_rel_abund.csv
│   ├── genus_read_numbers.csv
│   ├── species_rel_abund.csv
│   └── species_read_numbers.csv
└── 96_taxpasta/                # Taxpasta merged with lineage
    └── report_bracken_with_lineage.tsv

log/                            # Log files
├── kraken2/
├── bracken/
└── taxpasta/
```

## Tool Versions

The pipeline uses the following software versions (defined in `envs/` directory):

### Main Assembly-Based Pipeline
| Tool | Version | Purpose |
|------|---------|---------|
| fastp | 1.0.1 | Quality control |
| SPAdes | 4.2.0 | Metagenomic assembly |
| anvi'o | 7 | Scaffold processing |
| Pyrodigal | latest | Gene prediction |
| MMSeqs2 | 17.b804f | Gene clustering |
| seqkit | 2.10.1 | Sequence manipulation |
| Diamond | 2.1.13 | Sequence alignment |
| EggNOG-mapper | 2.1.13 | Functional annotation |
| CoverM | 0.7.0 | Abundance quantification |

### Taxonomic Profiling Pipeline
| Tool | Purpose | Environment |
|------|---------|-------------|
| Kraken2 | Taxonomic classification | kraken2.yaml |
| Bracken | Abundance re-estimation | kraken2.yaml |
| Pandas | Merge profiling reports | pandas.yaml |
| Taxpasta | Standardize & merge with lineage | taxpasta.yaml |

## Customization

### Adding New Samples
1. Add sample names to `config/config.yaml` under `samples:`
2. Assign samples to appropriate groups in `sample_groups:`
3. Ensure input files follow naming convention in `0_rawdata/`

### Modifying Parameters
Edit `config/config.yaml` to adjust:
- Thread counts for each tool
- Memory allocations
- Runtime limits
- Quality thresholds
- Database paths

### Adding New Rules
1. Create new conda environment in `envs/` if needed
2. Add rule to `Snakefile` with appropriate resources
3. Update `rule all:` inputs to include new outputs

## Troubleshooting

### Common Issues

**1. Memory Errors in SPAdes**
```bash
# Increase memory allocation in config.yaml
spades:
  memory: 512  # Increase from 256GB
```

**2. Database Path Errors**
```bash
# Verify database paths exist
ls -la /data/database/GTDB_214.1_diamond/
ls -la /data/database/eggnog/
```

**3. Conda Environment Issues**
```bash
# Clean conda environments
snakemake --use-conda --conda-cleanup-envs

# Force environment recreation
snakemake --use-conda --conda-create-envs-only
```

### Log Files
Check log files for detailed error messages:
- SPAdes logs: `log/2_spades/{sample}_spades.txt`
- Snakemake logs: `.snakemake/log/`

## Performance Tips

1. **Use SSD storage** for faster I/O operations
2. **Optimize thread allocation** based on available cores
3. **Monitor memory usage** during SPAdes assembly
4. **Use tmpdir** for intermediate files on high-speed storage
5. **Group similar samples** to maximize clustering efficiency

## Citation

If you use this pipeline in your research, please cite:

- **Snakemake**: Köster, J., & Rahmann, S. (2012). Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520-2522.
- **Individual tools**: Please cite the original publications for fastp, SPAdes, anvi'o, Prodigal, MMSeqs2, Diamond, EggNOG-mapper, and CoverM.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

For questions or issues:
- 📧 Contact: [Your Email]
- 🐛 Issues: [GitHub Issues](https://github.com/XuhanDeng/OMPs_Metagenomics_pipeline/issues)
- 📖 Documentation: This README

## Acknowledgments

This pipeline was developed for OMPs (Organometallic Pollutants) research and incorporates best practices from the metagenomics community.

---
*Generated with Claude Code*