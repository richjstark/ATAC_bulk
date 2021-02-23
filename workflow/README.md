# This is a bulk ATACseq analysis example
```git clone https://github.com/richjstark/ATAC_bulk
```

## Run locally, not using any batch system:
```snakemake --use-conda
```

## Run on a Cluster with a SLURM batch system:
### Create a conda environment for snakemake:
```conda create -name snakemake snakemake conda
```
### Run test job on cluster:
Handy reference for running on the cluster:
https://ulhpc-tutorials.readthedocs.io/en/latest/bio/snakemake/#cluster-configuration-for-snakemake
Run on Tinis with:
```source activate snakemake
conda activate snakemake
SLURM_ARGS="--partition {cluster.partition} --nodes {cluster.nodes} --ntasks {cluster.ntasks} --cpus-per-task {cluster.ncpus} --time {cluster.time} --job-name {cluster.job-name} --output {cluster.output} --error {cluster.error} --ntasks-per-node {cluster.ntaskspernode} --mem-per-cpu {cluster.mempercpu} --exclusive"
### Set -j argument for number of concurrent jobs to run:
snakemake -j 1 -pr --use-conda --keep-going --cluster-config workflow/cluster.yaml --cluster "sbatch $SLURM_ARGS"
```
## Report generation:
```snakemake --report
```
