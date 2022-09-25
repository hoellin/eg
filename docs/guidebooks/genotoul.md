# Genotoul cluster

Detailed information about Genotoul can be found on the [Genotoul website](http://www.genotoul.fr/en/).

## Basic usage

Connect through SSH to the Genotoul server:

```bash
ssh thoellinger@genologin.toulouse.inra.fr
```

Get account infos: `saccount_info thoellinger`  
List all available modules: `module av`  
Load one of them: `module load <name>`

Some typical modules:

```bash
module load bioinfo/samtools-1.9
module load system/Anaconda3-5.2.0
```

Install or update a Python package for user only (does not require sudo privilege):

```bash
pip install --user <package>
```

To update a package for user only: `pip install --user --upgrade <package>`

Remember to connect to a computation node before launching any resources-demanding program:

```bash
srun --pty bash
```

1 cpu core and 8 gb of ram by default with `srun` (2021). One may increase limits using, for instance:

```bash
srun --mem=64G --pty bash
```

## Submitting a script using `sbatch`

For more resources-demanding tasks or complex pipelines, it is better to write a dedicated script and to submit it to the queue using `sbatch`. Example:

```bash
echo "script content" | awk 'BEGIN{print "\#\!\/bin\/sh"} {print}' > script.sh
```

then launch the script on slurm:

```bash
sbatch --mem=4G --cpus-per-task=1 -J <job name> --mail-user=tristan.hoellinger@inserm.fr --mail-type=END,FAIL --workdir=$PWD --export=ALL -p workq script.sh
```

One can visualize what is in the queue with: `squeue -A thoellinger`

More sophisticated examples can be found on the notebooks over here, or on Genotoul's [wiki](https://wiki.genotoul.fr/wiki/Slurm).

[Resources available detailed here](http://bioinfo.genotoul.fr/index.php/faq/default_resources_faq/).
