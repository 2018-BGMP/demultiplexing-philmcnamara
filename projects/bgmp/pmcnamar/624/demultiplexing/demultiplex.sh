#! /bin/bash
#SBATCH --partition=long
#SBATCH --job-name=demultiplexing_philmcnamara
#SBATCH --output=demultiplex.output
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=phil.j.mcnamara@gmail.com

ml purge
ml slurm python3

./demultiplexing.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
