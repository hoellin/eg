# Genome-wide identification of enhancer/gene interactions

## Introduction

This GitHub Pages summarize some of my investigations on recent methods for genome-wide cell type specific identification of enhancer/gene relationships.

## How to contribute?

First of all, clone this repository to your local computer. Then, all you have to do is to modify the `.md` files wherever you find it useful, to commit and to push. The deployement of the GitHub pages is then automatically handled by `mkdocs` and GitHub Action.

### What if I want to build the documentation without pushing it?

* `mkdocs serve` - Start a live-reloading docs server.
* `mkdocs build --no-directory-urls` - Build the documentation site (local).


### FAQ

#### How to install `mkdocs`?

Note that you need to install both `mkdocs` and `mkdocs-material` to build the documentation locally.

```bash
pip install mkdocs
pip install mkdocs-material
```

#### Syntax highlighting for `mkdocs`

To enable syntax highlighting, use [Pygments](https://squidfunk.github.io/mkdocs-material/reference/code-blocks/#installation) and extra css stylesheets (for instance choose one [here](https://highlightjs.org/static/demo/) and download it [here](https://github.com/highlightjs/highlight.js/tree/master/src/styles)).

```yml
extra_css:
 - stylesheets/extra.css
markdown_extensions:
  - pymdownx.highlight
  - pymdownx.superfences
```

## About

Powered by [mkdocs.org](https://www.mkdocs.org) and [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/).