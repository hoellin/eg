# Conda basics

## Manage environments

Activate environment `env`:

```bash
conda activate env
```

Deactivate current environment: `deactivate`.

## Install packages

Install <package> in <env> environment:

```bash
conda install -n <env> <package>
```

To specify the channel where to find the package use `-c <channel>`

```bash
conda install -n <env> -c <channel> <package>
```

Example: one may use `conda-forge` or [`bioconda`](https://anaconda.org/bioconda).

**Do not** use `pip`, whenever possible! When no better solution is available, one can install packages inside a conda environment using the appropriate pip version using `python -m pip`. Sometimes (check `python -m pip --version` vs `pip --version`) `pip` directly does the job.

## Information about environements and packages

To get information about the current environment use:
```bash
conda info
```

To check the version of a Python package use:
```bash
pip freeze | grep <package>
```
