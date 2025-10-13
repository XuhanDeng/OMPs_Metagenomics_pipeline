#!/bin/bash
#SBATCH --job-name=k2_test2
#SBATCH --output=k2_test_old.out
#SBATCH --error=k2_test_old.err
#SBATCH --partition=compute
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=3900   # 16 × 7500 MB ≈ 120 GB total
#SBATCH --time=12:00:00
#SBATCH --account=research-ceg-wm
#SBATCH --ntasks=1

# Activate your conda environment
source ~/.bashrc
conda activate kraken2

# Paths
DB=/scratch/qou/database/kraken2
R1=rawdata/test_1P.fq.gz
R2=rawdata/test_2P.fq.gz
OUTDIR=kraken2_old

cd /scratch/qou/test
mkdir -p $OUTDIR

# Run Kraken2 (via kraken2 binary)
kraken2 \
  --db $DB \
  --threads 30 \
  --confidence 0.1 \
  --paired $R1 $R2 \
  --report $OUTDIR/test.kraken2.report \
  --output $OUTDIR/test.kraken2.out \
  --classified-out $OUTDIR/test.classified_#.fq \
  --unclassified-out $OUTDIR/test.unclassified_#.fq


