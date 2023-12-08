#!/usr/bin/bash
#SBATCH --job-name=deeptools_snakemake
#SBATCH --account=rrg-jabado-ab
#SBATCH --cpus-per-task=20
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err

module load python/3.8.10

CONFIG="config.yaml"
SNAKEFILE="snakefile"

# script to launch snakemake pipeline and create rule graph at end
snakemake \
 -s ${SNAKEFILE} \
 --configfile ${CONFIG} \
 -c 4

#snakemake --forceall --dag | dot -Tpdf > dag.pdf
