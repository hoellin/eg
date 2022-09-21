# Conda basics

## Manage environments

Activate environment `env`:

`source env/bin/activate` (deprecated) or:

```bash
conda activate env
```

Deactivate current environment: `deactivate`.

## Install packages

Install <package> in <env> environment:

```bash
conda install -n <env> <package>
```

To specify in which channel to search the package:

```bash
conda install -n <env> --channel <channel> <package>
```

```bash
conda install -n <env> -c <channel> <package>
```

Example: one may use `conda-forge` or [`bioconda`](https://anaconda.org/bioconda).

[Bioconda documentation](https://anaconda.org/bioconda).

**Do not** use `pip`, whenever possible! Although, if there is no other solution, one can install packages inside a conda environment using the appropriate pip version (the one specific to that environment). This should be accomplishable using `python -m pip`. Sometimes `pip` directly does the job. If any doubt, check `python -m pip --version` vs `pip --version`.

## Get info about current conda environment

```bash
conda info
```

Example:

> ```bash
> $ conda info
> 
>      active environment : base
>     active env location : /usr/local/bioinfo/src/Anaconda/Anaconda3-5.2.0
>             shell level : 1
>        user config file : /home/thoellinger/.condarc
>  populated config files : 
>           conda version : 4.5.4
>     conda-build version : 3.10.5
>          python version : 3.6.5.final.0
>        base environment : /usr/local/bioinfo/src/Anaconda/Anaconda3-5.2.0  (read only)
>            channel URLs : https://repo.anaconda.com/pkgs/main/linux-64
>                           https://repo.anaconda.com/pkgs/main/noarch
>                           https://repo.anaconda.com/pkgs/free/linux-64
>                           https://repo.anaconda.com/pkgs/free/noarch
>                           https://repo.anaconda.com/pkgs/r/linux-64
>                           https://repo.anaconda.com/pkgs/r/noarch
>                           https://repo.anaconda.com/pkgs/pro/linux-64
>                           https://repo.anaconda.com/pkgs/pro/noarch
>           package cache : /usr/local/bioinfo/src/Anaconda/Anaconda3-5.2.0/pkgs
>                           /home/thoellinger/.conda/pkgs
>        envs directories : /home/thoellinger/.conda/envs
>                           /usr/local/bioinfo/src/Anaconda/Anaconda3-5.2.0/envs
>                platform : linux-64
>              user-agent : conda/4.5.4 requests/2.18.4 CPython/3.6.5 Linux/3.10.0-514.26.2.el7.x86_64 centos/7 glibc/2.17
>                 UID:GID : 15707:10265
>              netrc file : None
>            offline mode : False
> ```

