# Git basics

[Any question? Take a look at this tutorial on OpenClassrooms](https://openclassrooms.com/fr/courses/1233741-gerez-vos-codes-source-avec-git).

## Create a new project

First create a directory, go within it, and initialize git using `git init`. If you want to be able to edit the project on this computer, do not use `git init --bare`, which is used to create a server in which there is only the `.git` repository.

Add files to git tracking with
```bash
git add <file name>
```

Then:
```bash
git commit -m "Initial comit"
```

## How I use Git to pull/push changes between my personal laptop, my professional laptot and the Genotoul cluster

The shared repository is named `shared`, and contains the `.git/`.

Get status: `git status`. Track new files: `git add <file name>`. To stop the tracking of a file: `git rm <file name>`. Remember to check whether new files are tracked before submitting a commit.

The following command could be used to commit all changes that have been done (except untracked files) and automatically untrack the files that have been deleted:

```bash
git commit -a -m "<short description>"
```

It is possible to commit only one or a few files using:

```bash
git commit -m "<short description>" <file1> <file2>
```

Undo last commit (`HEAD` designates the current commit, `HEAD^` the previous, `HEAD^^` or `HEAD~2` the n-2, etc):

```bash
git reset HEAD^
```

Verify commits logs: `git log` (or `git log --stat` for a summary).

Change description of the last commit (it opens a text editor):

```bash
git commit --amend
```

#### Branches

* `git branch` to see list of available branches
* `git checkout -b <branch_name>` to create a new branch
* `git checkout <branch_name>` to go on branch `branch_name`
* `git merge <branch_name>` to merge branch `branch_name` to current branch
* Delete a branch: `git branch -d <branch_to_delete>`

#### Untrack / retrack files

To start ignoring changes in a file:

```bash
git update-index --assume-unchanged notes_perso/site/
```

To start keeping track again:

```bash
 git update-index --no-assume-unchanged notes_perso/site/
```

### Git push (from local to Genotoul)

Only push once a day or so (when changes has been pushed, it's harder to reverse them than after a simple commit). **Push on the `swap` branch of Genotoul**, then go on the `master` branch on Genotoul and merge with `swap`. Note that it is not possible to push directly on master if the host is not configured as a server (`--bare` when initializing, which we shall not have used as we wanted to be able to edit files directly on Genotoul).

Speaking concretely, the following commands shall be used to push changes from my personal laptop to Genotoul:

1. `git st`

2. [...] (`git add <file or repository untracked yet>` if needed)

3. `git commit -a -m "<short description>"`

4. Then push:

   ```bash
   git push origin HEAD:swap
   ```

   `HEAD` designates the current commit (`HEAD^` the previous, `HEAD^^` or `HEAD~2` the n-2, etc).

Then on Genotoul:

5. `git checkout master` (if not currently on master - one can verify with `git branch`)

6. `git merge swap`

After step 6 is done, it is better to pull (Genotoul -> local), because since we did not push to `main` but to `swap`, it is as if we were still ahead by at least last commit of `main`.

So back on the local computer simply do:

7. `git pull`

### Git pull (on local computer to get changes done on Genotoul)

I should pull at 2 occasions: either when modification have been done on Genotoul, or when both Genotoul and local laptop are up-to-date in order not to induce lateness with respect to principal branch on Genotoul (when we push changes from local to Genotoul, it means we're several commits ahead of Genotoul as we do not push on main branch - see above).

```bash
git pull
```

## Connect to a remote project

To clone a project, first create a directory, go within it, and initialize git using `git init`.

Then connect the project to the remote repository using: 

```bash
git remote add <shortname> https://github.com/user/<reponame>`
```

 (`<shortname>` may be `upstream` or `origin` for instance).

Finally: `git pull <remote> <branch>` for instance `git pull upstream master`.

One can visualize all branches using `git branch`.

## "Push" the `mkdocs` site (outdated)

(Outdated guidelines - better use the custom push scripts now).<br />
(Outdated once again: there is no push script anymore, now we use the [`mkdocs gh-deploy` command](https://github.com/hoellin/hoellin.github.io)).

Better not to track the `site/`folder. Use `rsync` to push the `site/` folder to wherever you want it to be (e.g. on the Genotoul server).

```bash
rsync -avz notes_perso/site/ thoellinger@genologin.toulouse.inra.fr:/home/thoellinger/public_html/notes_perso/
```

In case `rsync` is not working, you can still delete all files on Genotoul and then use `scp`:

```bash
scp -r notes_perso/site/ thoellinger@genologin.toulouse.inra.fr:/home/thoellinger/public_html/notes_perso/
```

Note that scp is very slow + it does not overwrite files.