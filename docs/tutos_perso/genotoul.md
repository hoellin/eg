# Genotoul

```bash
ssh thoellinger@genologin.toulouse.inra.fr
```

Account infos: `saccount_info thoellinger`  
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

Connect to a node before launching any resources-demanding program:

```bash
srun --pty bash
```

No more than 1 cpu core and 8 gb of ram by default with `srun` though. But one can increase limits using, for instance:

```bash
srun --mem=64G --pty bash
```

For more resources-demanding tasks, it might be better to write a dedicated script then to execute it using `sbatch`. Example:

```bash
echo "script content" | awk 'BEGIN{print "\#\!\/bin\/sh"} {print}' > script.sh
```

then launch the script on slurm:

```bash
sbatch --mem=4G --cpus-per-task=1 -J <job name> --mail-user=tristan.hoellinger@inserm.fr --mail-type=END,FAIL --workdir=$PWD --export=ALL -p workq script.sh
```

One can visualize what is in the queue with: `squeue -A thoellinger`

[Resources available detailed here](http://bioinfo.genotoul.fr/index.php/faq/default_resources_faq/).
