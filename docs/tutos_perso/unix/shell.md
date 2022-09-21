# Shell basics

## Download files with `wget`
Recursively download files, without the ones with a .html extension:

```bash
wget -np -r -R html http://genoweb.toulouse.inra.fr/~sdjebali/material/
```

Note that even when there are no subfolders, the option -r is required to download the files. One can specify the depth limit with `--cut-dirs`.

Recursively find and delete files that contain "text" in their name:

```bash
find -type f -name '*text*' -delete
```

## SSH, SCP
```bash
ssh thoellinger@genologin.toulouse.inra.fr
```

Or use the `genologin` alias.

### Copy from local to distant

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

### Copy from distant to local
```bash
scp <username>@<ip of distant machine>:<source directory> <destination directory on the local machine>
```

Example:

```bash
scp thoellinger@genologin.toulouse.inra.fr:/home/thoellinger/work/tutorial/DATA/*.zip .
```

## Copy files with `rsync`

Example:

```bash
rsync -avz notes_perso/site/ thoellinger@genologin.toulouse.inra.fr:/home/thoellinger/public_html/notes_perso/
```



## Parsing files

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
* **use `-n` to comparing strings according to their numerical (base 10) value!** 

Example:

```bash
awk '$4!="0"' example_chr22/ABC_output/Neighborhoods/EnhancerList.txt |sort -t$'\t' -nk4
```

## Misc

Get size of file/repertory: 

```bash
du -sh <file>
```

The `h` in `-sh` is for human-readable size.

### System and network

Get ip adresses: `ifconfig`

### Jupyter notebook

Increase I/O size limit:

```bash
jupyter notebook --NotebookApp.iopub_data_rate_limit=1.0e10
```
