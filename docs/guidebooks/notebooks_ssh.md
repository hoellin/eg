---
# metadata
hide:
#  - navigation # Hide navigation
#  - toc        # Hide table of contents
---

# Jupyter notebooks through ssh tunneling

> One can launch a test server to check whether everything is working correctly, using the following command:
> 
> ```bash
> python -m http.server 8080
> ```

## Main steps

To start ssh tunneling, run the following on your (local) computer:

```bash
ssh -L 8080:localhost:8080 thoellinger@genologin.toulouse.inra.fr
```

Then launch Jupyter on the server (Genotoul here):

```bash
module load system/Anaconda3-5.2.0 && cd /work2/project/regenet/workspace/thoellinger/shared/notebooks && jupyter notebook --no-browser --port 8080
```

(One could also try with another port - for instance 8081).

Then, on your local computer, open a browser and go to `localhost:8080` (resp. `localhost:8081`). Enter the required token when prompted to access the notebook.

Note that one cannot run the notebook from a node (`srun --pty bash` then `jupyter notebook` does not work, because of denied permission). Instead, one can launch the notebook from standard user, then specify in the notebooks the cells that shall be launched from a node or on slurm. For instance, this is what is done [here](../../notes_ABC/generic_notebooks/turnkey_notebook_to_run_ABC_with_example_over_GM12878/#use-abc-makecandidateregionspy-to-define-candidate-regions).

## Convert a notebook into html from shell

One can simply use:

```bash
jupyter nbconvert --to html test_notebooks/first_notebook.ipynb
```

Or create a dedicated script (e.g. `convert_notebooks_to_html.sh`) with a content like:

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

I personally use thfollowing theme:

```bash
jt -t monokai -f fira -fs 10 -nf ptsans -nfs 11 -N -kl -cursw 2 -cursc r -cellw 95% -T
```

The default theme can be restored with:

```bash
jt -r
```

Note that these commands can be directly executed inside a notebook, using the `!` prefix, but the page needs to be reloaded to see the changes.

More info available here: https://github.com/dunovank/jupyter-themes


