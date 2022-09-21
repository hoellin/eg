# Git basics

[Any question? Take a look at this tutorial on OpenClassrooms](https://openclassrooms.com/fr/courses/1233741-gerez-vos-codes-source-avec-git).

## Create a new project

First create a directory, go within it, and initialize git using `git init`. If you want to be able to edit the project on this computer, do not use `git init --bare`, which is used to create a server in which there is only the `.git` repository.

If necessary:

```bash
git config --global color.diff auto
git config --global color.status auto
git config --global color.branch auto
git config --global user.name "hoellinger"
git config --global user.email tristan.hoellinger@inserm.fr
git config --global alias.st "status"
```

To add files to git's tracking: `git add <file name>`.

Then: `git commit -m "Initial comit"`.

## Use git to pull/push changes between my 2 laptops and Genotoul

The shared repository shall be named `shared`, and contains a `.git/` repository (see "Create a new project").

Get status: `git status` or `git st` with the alias. To track new files: `git add <file name>`. To remove some: `git rm <file name>`. Remember to check whether new files are tracked before submitting a commit.

The following command shall be used to commit all changes that have been done (except untracked files) and automatically untrack the files that have been deleted:

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

### Git push (from local to genotoul)

Only push once a day or so (when changes has been pushed, it's a lot harder to reverse them than after a commit). **Push on the `swap` branch of genotoul**, then go on genotoul and merge. It is not possible to push directly on master if the host is not configured as a server (`--bare` when initializing, which we shall not have used as we wanted to be able to edit files directly on Genotoul).

On the local computer:

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

> ```bash
> thoellinger@genologin2 /work2/project/regenet/workspace/thoellinger/shared $ git checkout master 
> Switched to branch 'master'
> thoellinger@genologin2 /work2/project/regenet/workspace/thoellinger/shared $ git branch
> * master
>   swap
> thoellinger@genologin2 /work2/project/regenet/workspace/thoellinger/shared $ git merge swap
> Updating 9d23ca9..065f626
> Fast-forward
> notes_perso/git.md | 80 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---
> 1 file changed, 77 insertions(+), 3 deletions(-)
> ```
> 

After step 6 is done, it is better to pull (genotoul -> local), because as we did not push to `main` but to `swap`, git considers we are still ahead by at least last commit of `main`.

So back on the local computer:

7. `git pull`

### Git pull (on local computer to get changes done on Genotoul)

I should pull at 2 occasions: either when modification have been done on Genotoul, or when both Genotoul and local laptop are up-to-date in order not to induce lateness with respect to principal branch on Genotoul (when we push changes from local to Genotoul, it means we're several commits ahead of Genotoul as we do not push on main branch - see above).

```bash
git pull
```

> ```bash
> $ git pull
> thoellinger@genologin.toulouse.inra.fr's password: 
> remote: Counting objects: 7, done.
> remote: Compressing objects: 100% (4/4), done.
> remote: Total 4 (delta 3), reused 0 (delta 0)
> Unpacking objects: 100% (4/4), 589 bytes | 589.00 KiB/s, done.
> From genologin.toulouse.inra.fr:/work2/project/regenet/workspace/thoellinger/shared
>    e22e231..3791f57  master     -> origin/master
> Updating 065f626..3791f57
> Fast-forward
>  notes_perso/git.md | 13 +++++++++++++
>  1 file changed, 13 insertions(+)
> ```

## Connect to a remote project (found on Github for instance)

To clone a project, first create a directory, go within it, and initialize git using `git init`.

Then connect the project to the remote repository using: 

```bash
git remote add <shortname> https://github.com/user/<reponame>`
```

 (`<shortname>` may be `upstream` or `origin` for instance).

Finally: `git pull <remote> <branch>` for instance `git pull upstream master`.

One can visualize all branches using `git branch`.

## "Push" the `mkdocs` site

(Outdated guidelines - better use the custom push scripts now).

No need to track it at the moment (maybe we will track it when it starts being too heavy for ssh copy to be quick enough). Use `rsync` (the best) instead, or `scp` if ever there is a problem with `rsync`.

```bash
rsync -avz notes_perso/site/ thoellinger@genologin.toulouse.inra.fr:/home/thoellinger/public_html/notes_perso/
```

If there is a problem with `rsync`, first delete all files on genotoul, then:

```bash
scp -r notes_perso/site/ thoellinger@genologin.toulouse.inra.fr:/home/thoellinger/public_html/notes_perso/
```

Note that scp is very slow + it does not overwrite files.