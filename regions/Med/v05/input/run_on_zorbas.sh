#!/bin/bash

#SBATCH --job-name=MEDnp8           # Job name
#SBATCH --output=MEDnp8.out         # Standard output and error log
#SBATCH --error=MEDnp8.err          # Error log
#SBATCH --partition=batch           # Error log
#SBATCH --nodes=2                   # Number of nodes
#SBATCH --ntasks-per-node=4         # Number of MPI tasks per node (adjust as needed)
#SBATCH --partition=batch           # Partition name

# Run the MPI job
cd ~/darwin3/run
mpirun -np 8 ./mitgcmuv

# Submit using "sbatch run_on_zorbas.sh"
