# CLAUDE.md — working notes & conventions for this repo

Guidance for Claude Code when maintaining the **claude-sternheimer** GitHub Pages site.
When the user says *"I added test X, update the website,"* follow this file to keep the
structure, formula style, and publish flow consistent.

Live site: https://rjguo1208.github.io/claude-sternheimer/ — served from `main` → `/docs`.

---

## 1. What this repo is

A static report site for the **active/rest defect $T$-matrix / Sternheimer** project
(a beyond-Born extension of EDI). Source of truth is Markdown (`research.md` and, later,
per-test files under `content/`); `tools/build_site.py` renders it to `docs/`.

Modeled on the structure & visual style of `github.com/rjguo1208/codex-sternheimer`,
**but with correct MathJax math** instead of plain-ASCII formulas.

## 2. Site structure (keep this shape)

- **`docs/index.html` — landing page. The first section is always the `Test Catalog`**:
  one row per piece of work = *Item · Type · Date · Status badge · one-line key result · link*.
  This is the most important element: a reader sees at a glance what was done and how it went.
  Then: Executive Summary (cards), Method at a Glance (key equation), Warnings/Not-published.
- **Larger items get their own subpage** in `docs/pages/<name>.html`, linked from the catalog.
- Every page: sticky top nav, a header (title + subtitle + meta pills), white "section" panels,
  consistent table styling, and a footer. Subpages also get an in-page **Contents** (TOC).

### Status badges
`badge ok` = Complete · `badge plan` = Planned (not yet run) · `badge prod` = headline production
result · `badge warn` = caveat. Planned catalog rows also get `class="planned"` (greyed).

## 3. Formula conventions (MUST render, not ASCII)

- Write **LaTeX**, never ASCII math. Inline `$…$`, display `$$…$$`.
- Every page loads MathJax v3 with this exact head snippet (the build script injects it):
  ```html
  <script>window.MathJax={tex:{inlineMath:[['$','$'],['\\(','\\)']],
    displayMath:[['$$','$$'],['\\[','\\]']]},
    options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']}};</script>
  <script id="MathJax-script" async
    src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
  ```
- Translation cheatsheet (ASCII → LaTeX): `omega0` → `$\omega_0$`;
  `||delta T^P_AA||_F` → `$\lVert\Delta T^{P}_{AA}\rVert_F$`; `U†U-I` → `$U^\dagger U-I$`;
  `g_AA` → `$g_{AA}$`; `T^P_AA` → `$T^{P}_{AA}$`; `DeltaV` → `$\Delta V$`; `Xi` → `$\Xi$`;
  `<m|V|n>` → `$\langle m|V|n\rangle$`.
- `\boxed`, `\lVert/\rVert`, `\dagger`, `\underbrace`, `\xrightarrow`, `\mathcal`, `\mathbb`
  are all supported by the default `tex-mml-chtml` (AMS) build.
- Put **code** (Fortran, shell, pseudocode) in fenced blocks ```` ```lang ````. MathJax skips
  `<pre>/<code>`, so `$`/`<`/`>` inside code are safe (the builder HTML-escapes code).
- **Self-check after every build:** open the page and confirm sub/superscripts, `\dagger`,
  and norm bars render; grep the HTML for leftover ASCII math in numbers/labels.

## 4. How the build works (`tools/build_site.py`)

Stdlib-only Markdown→HTML converter. Key invariant: it **protects** fenced code, `$$…$$`,
`$…$`, and `` `code` `` as placeholders *before* any Markdown processing, then restores math
**verbatim** (raw) and code **HTML-escaped**. So MathJax — not the converter — renders math.

- `convert_research()` turns `research.md` into the theory page (title/subtitle → header,
  each `## ` section → a `<section>` panel, `### ` → `<h3>`, pipe tables, blockquotes, lists,
  fenced code, display/inline math). Section ids are `sec-<N>`; an auto TOC is built from `## `.
- `CATALOG` (a Python list) drives the landing-page Test Catalog table.
- `build_index()` assembles the landing page from `CATALOG` + summary/glance/warnings blocks.
- Self-checks print section/eq/table counts and assert no placeholder leaked (`\x00`).

## 5. Recipe: add a new test (the common request)

1. Put the writeup in `content/<test>.md` (Markdown + LaTeX, same conventions as `research.md`).
   Summarize numbers in **tables**; copy any figures to `docs/assets/` and reference relatively.
2. In `tools/build_site.py`:
   - add a small `build_<test>()` (mirror `build_theory()`; or generalize it to take a source
     path + output path), writing `docs/pages/<test>.html`;
   - add a row to `CATALOG` with the right **Status** badge, date, one-line key result, and link.
3. Rebuild and verify:
   ```bash
   ml anaconda && conda activate Tmat
   python tools/build_site.py
   ```
4. Commit & push (see §6). Flip the catalog badge `plan → ok` (or `prod`) once the test passes.

## 6. Publish flow

```bash
git add docs README.md CLAUDE.md .gitignore research.md tools content
git commit -m "..."        # end with: Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
git push origin main
```
Pages source is `main` `/docs` (set once in Settings → Pages; `gh` CLI is absent on the host).
GitHub auth here is via SSH (`git@github.com`), already working.

## 7. Never publish (enforced by `.gitignore`)

Raw wavefunctions (`*.wfc`), QE `*.save/`, cubes (`*.cube`), arrays
(`*.npy/*.npz/*.bin/*.h5/*.dat`), logs (`*.out/*.err/*.log`), `.claude/`, and anything
key/credential-like. These are **summarized in tables** on the site, never committed.
The published payload is `docs/` plus the small Markdown/Python sources.

## 8. Environment notes

- Python: use `ml anaconda; conda activate Tmat` (env at
  `/home/x-rg47749/.conda/envs/2024.02-py311/Tmat`, Python 3.11). The build needs no
  third-party packages, so the env is only a convenience.
- Source project (theory note origin & EDI code): `research.md` here is derived from the
  working note in this HPC dir; the EDI plugin lives at
  `/anvil/projects/x-che190065/rjguo/qe-7.5/edi-dev`.
