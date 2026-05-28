# Sternheimer Electron–Defect T-matrix

Static report site for the **active/rest partitioning of the defect $T$-matrix** —
cRPA-style downfolding, a coupling-second-order effective potential, and a
$k$-decoupled Sternheimer solution. This is a *beyond-Born* extension of the
[EDI](https://github.com/rjguo1208/edi-dev) electron–defect workflow.

**Live site:** https://rjguo1208.github.io/claude-sternheimer/

At this stage the site presents the **theory & method** only; numerical tests are
forthcoming and are listed as *Planned* in the landing-page Test Catalog.

## Repository layout

| Path | Role |
|------|------|
| `docs/` | The published GitHub Pages site (this is what Pages serves). |
| `docs/index.html` | Landing page: **Test Catalog** + executive summary + method-at-a-glance. |
| `docs/pages/theory.html` | Full theory & method note, MathJax-rendered. |
| `docs/assets/style.css` | Shared stylesheet for all pages. |
| `research.md` | **Source** theory & method note (Markdown + LaTeX). |
| `tools/build_site.py` | Regenerates `docs/` from `research.md` (Python stdlib only). |
| `content/` | Reserved for additional per-test Markdown sources. |

## Build

The generator has **no third-party dependencies** (standard library only), so any
Python 3 works. On Anvil, use the project env:

```bash
ml anaconda
conda activate Tmat            # or: /home/x-rg47749/.conda/envs/2024.02-py311/Tmat/bin/python
python tools/build_site.py
```

This rewrites `docs/index.html` and `docs/pages/theory.html`. Open them locally in a
browser to verify the MathJax rendering before publishing.

## Math rendering

All equations are written in LaTeX and rendered client-side with **MathJax v3**
(`$…$` inline, `$$…$$` display). The build script *protects* every math span before
Markdown processing and restores it verbatim, so the source LaTeX is exactly what
MathJax sees. (This is the fix for the plain-ASCII-math problem in the template
this site is modeled on.)

## Publishing (GitHub Pages)

The site is served from the `main` branch `/docs` folder:

```bash
git add docs README.md CLAUDE.md .gitignore research.md tools content
git commit -m "Update site"
git push origin main
```

Then in **GitHub → Settings → Pages**: *Deploy from a branch* → **main** → **/docs**.
(The `gh` CLI is not installed on the build host, so Pages is enabled once, manually,
in the web UI.)

## Not published

Raw wavefunctions, QE `*.save/` dirs, cube/volumetric grids, `*.npy/*.npz/*.bin/*.h5/*.dat`,
scheduler logs, and any local/credential files are excluded by `.gitignore` and are only
ever *summarized* in tables on the site. See `CLAUDE.md` for the full convention.
