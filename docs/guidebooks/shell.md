# Shell basics

## File system

### Recursively find files

```bash
find -type f -name '*text*'
```

One can also delete them using:
```bash
find -type f -name '*text*' -delete
```

### Download files with `wget`
Recursively download files, without the ones with a .html extension:

```bash
wget -np -r -R html http://genoweb.toulouse.inra.fr/~sdjebali/material/
```

Note that even when there are no subfolders, the option -r is required to download the files. Also note that one can specify the depth limit with `--cut-dirs`.

### Parsing files

Sort a tabular-separated file wrt column 4:

```bash
sort -t$'\t' -nk4 <file>
```

Sort a ","-separated file wrt column 4:

```bash
sort -t, -nk4 <file>
```

or `sort -t ',' -nk4 <file>`.

Warning:

* use `-h` to sort human-readable numbers (1K, 2G, etc)
* **use `-n` to compare strings according to their numerical (base 10) value!**

Example:

```bash
awk '$4!="0"' example_chr22/ABC_output/Neighborhoods/EnhancerList.txt |sort -t$'\t' -nk4
```

### Copy files with `rsync`

When copying files from a distant machine to the local machine (resp. from distant to local) on a regular basis, each time with a few modifications, `rsync` is a lot faster than `scp`.

Example:

```bash
rsync -avz notes_perso/site/ thoellinger@genologin.toulouse.inra.fr:/home/thoellinger/public_html/notes_perso/
```

### Miscellaneous

Get size of file/repertory: 

```bash
du -sh <file>
```

The `h` in `-sh` is for human-readable size.

## SSH

### Basics

Connect to a Genotoul login node using:
```bash
ssh thoellinger@genologin.toulouse.inra.fr
```

### SCP - copy from local to distant

```bash
scp <file> <username>@<ip adress>:<destination directory>
```

For a directory:

```bash
scp -r <file> <username>@<ip adress>:<destination directory>
```

If using ipv6:

```bash
scp -6 <file> <username>@[ipv6 adress]:<destination directory>
```

```bash
scp -6r <file> <username>@[ipv6 adress]:<destination directory>
```

Example:

```bash
scp catalogue_GWAS.html thoellinger@genologin.toulouse.inra.fr:/home/thoellinger/public_html/autoupdate/GWAS_catalog.html
```

### SCP - copy from distant to local
```bash
scp <username>@<ip of distant machine>:<source directory> <destination directory on the local machine>
```

Example:

```bash
scp thoellinger@genologin.toulouse.inra.fr:/home/thoellinger/work/tutorial/DATA/*.zip .
```

## Jupyter notebook

### Increase I/O size limit

```bash
jupyter notebook --NotebookApp.iopub_data_rate_limit=1.0e10
```

## Useful aliases

I always include the following aliases in my `.bash_aliases`:

```bash
alias c='clear'
alias genologin='ssh -Y thoellinger@genologin.toulouse.inra.fr'
alias mgeno='sshfs thoellinger@genologin.toulouse.inra.fr:/home/thoellinger /home/hoellinger/mnt/ge\
notoul -o uid=1000,gid=1000,follow_symlinks,allow_other,workaround=rename'
```
