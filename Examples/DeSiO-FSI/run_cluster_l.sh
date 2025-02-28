#!/bin/bash -l
#SBATCH --job-name=desio_fsi
###SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=5G
#SBATCH --time=24:00:00
#SBATCH --partition=isd
#SBATCH --mail-user=c.hente@isd.uni-hannover.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output output.out
#SBATCH --error error.err

module load intel/2021a 

/bigwork/nhgehent/DeSiO/Compile/ifort/DeSiO
