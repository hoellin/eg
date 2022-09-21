---
# metadata
hide:
#  - navigation # Hide navigation
#  - toc        # Hide table of contents
---

# Jupyter notebooks through ssh tunneling

One can launch a test server to check whether everything is working correctly, using the following command:

```bash
python -m http.server 8080
```

## Launch Jupyter

To start ssh tunneling, run the following on the local computer:

```bash
ssh -L 8080:localhost:8080 thoellinger@genologin.toulouse.inra.fr
```

Then launch Jupyter on the server (Genotoul here):

```bash
module load system/Anaconda3-5.2.0 && cd /work2/project/regenet/workspace/thoellinger/shared/notebooks && jupyter notebook --no-browser --port 8080
```

or try another port:

```bash
ssh -L 8081:localhost:8081 thoellinger@genologin.toulouse.inra.fr
```

```bash
module load system/Anaconda3-5.2.0 && cd /work2/project/regenet/workspace/thoellinger/shared/notebooks && jupyter notebook --no-browser --port 8081
```

Now you can open `localhost:8080` (resp. `localhost:8081`) in a browser and enter the required token to access the notebook.

Note that one cannot run the notebook from a node (`srun --pty bash` then `jupyter notebook` does not work, because of denied permission). Instead, one can launch the notebook from standard user, then specify in the notebooks the cells that shall be launched from a node or on slurm. For instance, this is what is done [here](../../notes_ABC/generic_notebooks/turnkey_notebook_to_run_ABC_with_example_over_GM12878/#use-abc-makecandidateregionspy-to-define-candidate-regions).

## Convert a notebook into html from shell

One can simply use:

```bash
jupyter nbconvert --to html test_notebooks/first_notebook.ipynb
```

Or use a dedicated script such as:

> ```bash
> #!/bin/sh                                                                          
> path_to_html="/home/thoellinger/public_html/autoupdate/"
> 
> module load system/Anaconda3-5.2.0
> jupyter nbconvert --to html /work2/project/regenet/workspace/thoellinger/notebooks\
> /first_notebook/first_notebook.ipynb --output "${path_to_html}first_notebook"
> ```
>

## Use a custom theme with `jt`

For instance I use:

```bash
jt -t monokai -f fira -fs 10 -nf ptsans -nfs 11 -N -kl -cursw 2 -cursc r -cellw 95% -T
```

One can always restore the default theme with:

```bash
jt -r
```

Note that one can execute these commands directly in a notebook / when jupyter is running, but the page has to be refreshed for the changes to apply.

More info available here: https://github.com/dunovank/jupyter-themes


