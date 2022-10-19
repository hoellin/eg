# Mkdocs

For full documentation visit [mkdocs.org](https://www.mkdocs.org).

### Commands

* `mkdocs serve --no-directory-urls` - Start the live-reloading docs server (`--no-directory-urls` is needed to avoid inconsistent redirections)
* `mkdocs build --no-directory-urls` - Build the documentation site.

### Syntax highlighting

To enable syntax highlighting, use [Pygments](https://squidfunk.github.io/mkdocs-material/reference/code-blocks/#installation) and extra css stylesheets (for instance choose one [here](https://highlightjs.org/static/demo/) and download it [here](https://github.com/highlightjs/highlight.js/tree/master/src/styles)).

```yml
extra_css:
 - stylesheets/extra.css
markdown_extensions:
  - pymdownx.highlight
  - pymdownx.superfences
```

## Build options

Use the following:

```bash
mkdocs build --no-directory-urls
```

