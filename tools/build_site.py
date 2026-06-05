#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_site.py  —  Static-site generator for the Sternheimer electron-defect repo.

Converts the Markdown research note (research.md) into a MathJax-rendered theory
sub-page (docs/pages/theory.html) and regenerates the landing page
(docs/index.html) with the "Test Catalog" table.

Design goals
------------
* Math fidelity: $...$ and $$...$$ blocks are PROTECTED before any Markdown
  processing and restored verbatim, so MathJax (not the converter) renders them.
  This is the fix for the codex-sternheimer template's ASCII-math defect.
* Zero third-party dependencies: standard library only (works in any Python 3).
* Repeatable: to add a test page later, drop a Markdown file in content/ and
  register it in PAGES / CATALOG below, then re-run this script.

Usage
-----
    python tools/build_site.py
"""
import os, re, html, datetime

ROOT      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESEARCH  = os.path.join(ROOT, "research.md")
PLAN      = os.path.join(ROOT, "plan.md")
NOTE_KNORM = os.path.join(ROOT, "content", "note_kprime_norm.md")
NOTE_RESULTS = os.path.join(ROOT, "content", "note_tmatrix_results.md")
DOCS      = os.path.join(ROOT, "docs")
PAGES_DIR = os.path.join(DOCS, "pages")

GEN_DATE  = "2026-05-28"
SITE_TITLE = "Sternheimer Electron-Defect T-matrix"

# ======================================================================
#  Markdown -> HTML  (math-protected, stdlib only)
# ======================================================================
NUL = "\x00"

def _protect(md, store):
    """Replace fenced code, display math, inline math, inline code with
    null-delimited placeholders so Markdown processing cannot mangle them."""
    # 1) fenced code blocks ```lang \n ... \n```
    def fence(m):
        store["c"].append((m.group(1) or "", m.group(2)))
        return "%sC%d%s" % (NUL, len(store["c"]) - 1, NUL)
    md = re.sub(r"```[ \t]*([A-Za-z0-9_+-]*)[ \t]*\n(.*?)\n```", fence, md, flags=re.DOTALL)
    # 2) display math $$ ... $$  (protect BEFORE inline so inner $x$ survives)
    def disp(m):
        store["d"].append(m.group(1))
        return "%sD%d%s" % (NUL, len(store["d"]) - 1, NUL)
    md = re.sub(r"\$\$(.*?)\$\$", disp, md, flags=re.DOTALL)
    # 3) inline math $ ... $  (single line, no nested $)
    def inl(m):
        store["i"].append(m.group(1))
        return "%sI%d%s" % (NUL, len(store["i"]) - 1, NUL)
    md = re.sub(r"\$([^$\n]+?)\$", inl, md)
    # 4) inline code ` ... `
    def ic(m):
        store["k"].append(m.group(1))
        return "%sK%d%s" % (NUL, len(store["k"]) - 1, NUL)
    md = re.sub(r"`([^`]+?)`", ic, md)
    return md

def _restore(text, store):
    """Restore placeholders. Math is emitted RAW (for MathJax); code is escaped."""
    for n, (lang, code) in enumerate(store["c"]):
        cls = ' class="language-%s"' % lang if lang else ""
        repl = "<pre><code%s>%s</code></pre>" % (cls, html.escape(code))
        text = text.replace("%sC%d%s" % (NUL, n, NUL), repl)
    for n, tex in enumerate(store["d"]):
        repl = '<div class="math">\n$$%s$$\n</div>' % tex
        text = text.replace("%sD%d%s" % (NUL, n, NUL), repl)
    for n, tex in enumerate(store["i"]):
        text = text.replace("%sI%d%s" % (NUL, n, NUL), "$" + tex + "$")
    for n, code in enumerate(store["k"]):
        text = text.replace("%sK%d%s" % (NUL, n, NUL), "<code>" + html.escape(code) + "</code>")
    return text

def _inline(text):
    """Escape HTML, then apply images / links / bold / italic. Placeholders untouched."""
    text = html.escape(text, quote=False)
    # images ![alt](src)  (before links, since the syntax contains [..](..))
    text = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)",
                  r'<img src="\2" alt="\1" style="max-width:100%;height:auto;display:block;'
                  r'margin:1.2rem auto;border:1px solid #d0d7de;border-radius:6px">', text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<em>\1</em>", text)
    return text

def _slug(headtext):
    m = re.match(r"\s*(\d+)\.", headtext)
    if m:
        return "sec-" + m.group(1)
    s = re.sub(r"%s[A-Z]\d+%s" % (NUL, NUL), "", headtext)      # drop math placeholders
    s = re.sub(r"[^A-Za-z0-9]+", "-", s).strip("-").lower()
    return s or "sec"

def _render_table(header, body_rows):
    def cells(row):
        row = row.strip()
        if row.startswith("|"): row = row[1:]
        if row.endswith("|"):   row = row[:-1]
        return [c.strip() for c in row.split("|")]
    h = "".join("<th>%s</th>" % _inline(c) for c in cells(header))
    out = ['<div class="table-wrap"><table><thead><tr>%s</tr></thead><tbody>' % h]
    for r in body_rows:
        tds = "".join("<td>%s</td>" % _inline(c) for c in cells(r))
        out.append("<tr>%s</tr>" % tds)
    out.append("</tbody></table></div>")
    return "".join(out)

def _render_list(lines):
    ordered = bool(re.match(r"\s*\d+\.\s+", lines[0]))
    tag = "ol" if ordered else "ul"
    items, cur = [], None
    for ln in lines:
        m = re.match(r"\s*(?:[-*+]|\d+\.)\s+(.*)$", ln)
        if m:
            if cur is not None: items.append(cur)
            cur = m.group(1)
        else:                                   # continuation line
            cur = (cur + " " + ln.strip()) if cur is not None else ln.strip()
    if cur is not None: items.append(cur)
    # GitHub-style task list: "[ ] ..." / "[x] ..." -> disabled checkboxes
    has_task = any(re.match(r"\[[ xX]\]\s+", it) for it in items)
    lis = []
    for it in items:
        mt = re.match(r"\[([ xX])\]\s+(.*)$", it)
        if mt:
            chk = " checked" if mt.group(1) in ("x", "X") else ""
            lis.append('<li class="task"><input type="checkbox" disabled%s> %s</li>'
                       % (chk, _inline(mt.group(2))))
        else:
            lis.append("<li>%s</li>" % _inline(it))
    cls = ' class="tasklist"' if has_task else ""
    return "<%s%s>%s</%s>" % (tag, cls, "".join(lis), tag)

def _is_block_start(line):
    s = line.strip()
    if re.match(r"^#{1,6}\s", line): return True
    if re.match(r"^%s[DC]\d+%s$" % (NUL, NUL), s): return True
    if s.startswith(">"): return True
    if re.match(r"\s*(?:[-*+]|\d+\.)\s+", line): return True
    return False

def _parse_blocks(md):
    lines = md.split("\n")
    out, i, n = [], 0, len(md.split("\n"))
    while i < len(lines):
        line = lines[i]
        if line.strip() == "":
            i += 1; continue
        # headings (### / ####)
        m = re.match(r"^(#{1,6})\s+(.*)$", line)
        if m:
            lv = len(m.group(1)); out.append("<h%d>%s</h%d>" % (lv, _inline(m.group(2).strip()), lv))
            i += 1; continue
        # standalone display-math / fenced-code placeholder = its own block
        if re.match(r"^%s[DC]\d+%s$" % (NUL, NUL), line.strip()):
            out.append(line.strip()); i += 1; continue
        # blockquote
        if line.lstrip().startswith(">"):
            buf = []
            while i < len(lines) and lines[i].lstrip().startswith(">"):
                buf.append(re.sub(r"^\s*>\s?", "", lines[i])); i += 1
            out.append("<blockquote>%s</blockquote>" % _parse_blocks("\n".join(buf)))
            continue
        # table (header row followed by |---|--- separator)
        if "|" in line and i + 1 < len(lines) and re.match(r"^\s*\|?[\s:|-]*-[\s:|-]*\|?\s*$", lines[i+1]):
            header = line; i += 2; body = []
            while i < len(lines) and lines[i].strip() != "" and "|" in lines[i]:
                body.append(lines[i]); i += 1
            out.append(_render_table(header, body)); continue
        # list
        if re.match(r"^\s*(?:[-*+]|\d+\.)\s+", line):
            buf = []
            while i < len(lines) and lines[i].strip() != "" and \
                  (re.match(r"^\s*(?:[-*+]|\d+\.)\s+", lines[i]) or lines[i].startswith("  ")):
                buf.append(lines[i]); i += 1
            out.append(_render_list(buf)); continue
        # paragraph
        buf = [line]; i += 1
        while i < len(lines) and lines[i].strip() != "" and not _is_block_start(lines[i]):
            buf.append(lines[i]); i += 1
        out.append("<p>%s</p>" % _inline(" ".join(b.strip() for b in buf)))
    return "\n".join(out)

def convert_doc(md, want_subtitle=True):
    """Convert a Markdown doc (research.md / plan.md / per-test) to HTML pieces.
    Returns dict(title, subtitle, preamble, body, toc). 'preamble' is any content
    between the H1 (and optional H3 subtitle) and the first '## ' section."""
    lines = md.split("\n")
    title_md = lines[0][2:].strip() if lines[0].startswith("# ") else SITE_TITLE
    start, subtitle_md = 1, ""
    if want_subtitle:
        for j in range(1, min(8, len(lines))):
            if lines[j].startswith("### "):
                subtitle_md = lines[j][4:].strip(); start = j + 1; break
            if lines[j].startswith("## "):
                break
    rest_md = "\n".join(lines[start:])

    store = {"c": [], "d": [], "i": [], "k": []}
    rest_md = _protect(rest_md, store)

    # preamble (before first '## ') + H2 sections
    pre_lines, sections, cur = [], [], None
    for ln in rest_md.split("\n"):
        if ln.startswith("## "):
            if cur: sections.append(cur)
            cur = {"head": ln[3:].strip(), "body": []}
        elif re.match(r"^---+\s*$", ln.strip()):
            continue
        elif cur is None:
            pre_lines.append(ln)
        else:
            cur["body"].append(ln)
    if cur: sections.append(cur)

    preamble_html = _parse_blocks("\n".join(pre_lines)) if any(l.strip() for l in pre_lines) else ""

    toc, html_sections = [], []
    for s in sections:
        slug = _slug(s["head"]); head_html = _inline(s["head"])
        toc.append((slug, head_html))
        inner = _parse_blocks("\n".join(s["body"]))
        html_sections.append(
            '<section id="%s"><h2>%s</h2>\n%s\n</section>' % (slug, head_html, inner))
    body_html = "\n".join(html_sections)

    def inline_only(t):
        st = {"c": [], "d": [], "i": [], "k": []}
        return _restore(_inline(_protect(t, st)), st)

    return {
        "title": inline_only(title_md),
        "subtitle": inline_only(subtitle_md),
        "preamble": _restore(preamble_html, store),
        "body": _restore(body_html, store),
        "toc": [(sl, _restore(tx, store)) for sl, tx in toc],
    }

# ======================================================================
#  HTML templates
# ======================================================================
MATHJAX = (
    '<script>window.MathJax={tex:{inlineMath:[[\'$\',\'$\'],[\'\\\\(\',\'\\\\)\']],'
    'displayMath:[[\'$$\',\'$$\'],[\'\\\\[\',\'\\\\]\']]},'
    'options:{skipHtmlTags:[\'script\',\'noscript\',\'style\',\'textarea\',\'pre\',\'code\']}};</script>\n'
    '<script id="MathJax-script" async '
    'src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>'
)

def _topnav(active, prefix=""):
    def a(href, label, key):
        cls = ' style="color:#fff;text-decoration:underline"' if key == active else ""
        return '<a href="%s%s"%s>%s</a>' % (prefix, href, cls, label)
    return ('<nav class="topnav"><div class="inner">'
            '<span class="brand">Sternheimer&nbsp;EDI</span>'
            '%s%s%s%s%s</div></nav>' % (a("index.html", "Home", "home"),
                                      a("pages/theory.html", "Theory &amp; Method", "theory"),
                                      a("pages/plan.html", "Implementation Plan", "plan"),
                                      a("pages/results.html", "Results", "results"),
                                      a("pages/note-kprime-normalization.html", "Note: k&prime;-norm", "note")))

def page_shell(title, head_html, nav_html, body_html, css_href):
    return """<!doctype html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<link rel="stylesheet" href="{css}">
{mathjax}
</head><body>
{nav}
{head}
<main>
{body}
</main>
<footer>Generated {date} from <code>research.md</code> in the HPC working directory.
Static HTML; equations rendered client-side with MathJax v3. No raw wavefunctions, cubes,
binaries, or logs are published.</footer>
</body></html>""".format(title=html.escape(title), css=css_href, mathjax=MATHJAX,
                          nav=nav_html, head=head_html, body=body_html, date=GEN_DATE)

# ---------- theory page ----------
def build_theory():
    with open(RESEARCH, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=True)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">{s}</p>'
              '<div class="meta"><span class="pill">Theory &amp; Method note</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], s=r["subtitle"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Theory & Method",
                     header, _topnav("theory", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "theory.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_plan():
    with open(PLAN, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>EDT &mdash; Implementation Plan</h1>'
              '<p class="subtitle">How to build the downfolding + Sternheimer electron&ndash;defect '
              '$T$-matrix package on Quantum ESPRESSO, reusing the EDI code for $\\Delta V$, '
              'Wannier rotation/interpolation, and transport.</p>'
              '<div class="meta"><span class="pill">Implementation plan</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Implementation Plan",
                     header, _topnav("plan", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "plan.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_note():
    with open(NOTE_KNORM, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">EDT implementation note: fixing the rest-space '
              '$k\'$-sum normalization before quoting a physical $\\tilde V$ (P2 &rarr; P3).</p>'
              '<div class="meta"><span class="pill">Implementation note</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Note: k'-sum normalization",
                     header, _topnav("note", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "note-kprime-normalization.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_results():
    with open(NOTE_RESULTS, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">Numerical results on the MoS&#8322; S-vacancy: the full downfolded '
              'potential block $\\tilde V=M+\\Sigma$ (P3), the active-space $T$-matrix '
              '$T_{{PP}}=[1-\\tilde V G^A]^{{-1}}\\tilde V$ (P5-a), and the Wannier/Koster&ndash;Slater '
              'locality study (P5-b).</p>'
              '<div class="meta"><span class="pill">Numerical results</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Numerical results (P3–P5)",
                     header, _topnav("results", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "results.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

# ---------- landing page ----------
# Test Catalog rows: (item, type, date, badge_class, badge_label, summary, link_html)
CATALOG = [
    ("Active/Rest T-matrix &amp; Sternheimer theory", "Theory", GEN_DATE, "ok", "Complete",
     "Two-layer beyond-Born T-matrix: rest dressing via a $k$-decoupled Sternheimer ladder "
     "(exact to all orders in $V_{QQ}$) feeding an exact active-space inversion "
     "$T_{PP}=[1-\\tilde V G^A]^{-1}\\tilde V$.",
     '<a href="pages/theory.html">Open theory &amp; method &rarr;</a>'),
    ("EDT package implementation plan", "Plan", GEN_DATE, "ok", "Complete",
     "File-by-file outline + code snippets for the QE plug-in: reuses EDI for $\\Delta V$/Wannier/transport; "
     "new rest-space Sternheimer solve (QE <code>ccgsolve_all</code>), $V_{QQ}$ ladder, $\\tilde V$ assembly, "
     "and the small active inversion. Rest sum over the full BZ.",
     '<a href="pages/plan.html">Open implementation plan &rarr;</a>'),
    ("Note: $k'$-sum normalization (P2&rarr;P3)", "Note", GEN_DATE, "ok", "Resolved",
     "The rest-space $k'$-sum is a BZ-integral discretization ⇒ carries $1/N_k$, turning the naive "
     "$\\sum_{k'}\\!\\approx\\!-70$ Ry into the physical $\\Sigma_{nn}\\!\\approx\\!-0.5$ Ry. Closure "
     "sum rule + Born-limit mobility anchor to confirm (and expose any residual $N_{sc}$).",
     '<a href="pages/note-kprime-normalization.html">Open implementation note &rarr;</a>'),
    ("Downfolded potential $\\tilde V=M+\\Sigma$ (diagonal + full block)", "Result", "2026-06-01", "prod", "Complete",
     "Beyond-Born $\\tilde V$ on the MoS&#8322; active manifold from the per-$k'$ Sternheimer solve summed over the full BZ "
     "with $1/N_k$: diagonal (closure-validated) and the full $1584\\times1584$ block via a pool-parallel $k'$-sum "
     "(Hermitian to $9\\times10^{-12}$). Moderate rest dressing &mdash; VBM (band 13, $K$): $\\tilde V_{nn}=+0.112$ Ry "
     "($M=+0.246$, $\\Sigma=-0.134$); $\\lVert\\Sigma\\rVert/\\lVert M\\rVert=0.66$.",
     '<a href="pages/results.html">Open numerical results &rarr;</a>'),
    ("Active-space $T$-matrix $T_{PP}$ (P5-a)", "Result", "2026-06-01", "prod", "Complete",
     "Coarse-grid resummation $T_{PP}=[1-\\tilde V G^A]^{-1}\\tilde V$ ($G^A$ carries $1/N_k$). "
     "Born-limit validated ($T\\!\\to\\!\\tilde V$ to $10^{-4}$). On the VBM (band 13), rest-space <em>redistributes</em> "
     "the resummed scattering &mdash; screens the forward channel but enhances the norm "
     "($\\lVert T_{PP}\\rVert\\!=\\!177$ vs $\\lVert T_M\\rVert\\!=\\!70$ Ry); resonant at the band edge.",
     '<a href="pages/results.html#sec-2">Open results &rarr;</a>'),
    ("QE-Hamiltonian Sternheimer validation", "Test", "2026-05-31", "ok", "Complete",
     "Per-$k$ solve of $Q(\\omega_0-H_0)Q$ via projected PCG (QE <code>h_psi</code> matvec): "
     "$\\langle\\psi|H_0|\\psi\\rangle\\!=\\!\\varepsilon$ gate to $6\\times10^{-10}$ eV across all ranks; "
     "explicit rest-band sum converges to the all-band Sternheimer value (Born limit $T\\!\\to\\!V$ at $10^{-13}$ Ry).",
     '<a href="pages/plan.html">see P0&ndash;P3 in the plan &rarr;</a>'),
    ("Rest dressing ladder convergence", "Test", "&mdash;", "plan", "Planned",
     "Successive-ratio $\\lVert\\tilde V^{(m+1)}-\\tilde V^{(m)}\\rVert/\\lVert\\tilde V^{(m)}-"
     "\\tilde V^{(m-1)}\\rVert$ vs. rest $k$-grid and band cutoff (geometric rate "
     "$\\rho\\sim\\lVert V_{QQ}\\rVert/\\Delta$).", "&mdash;"),
    ("Wannier representation &amp; Koster&ndash;Slater (P5-b)", "Result", "2026-06-01", "ok", "Resolved",
     "Wannierizing $\\tilde V$ needs $U(k)$ in the same Bloch gauge as the evc that build $M$. The original "
     "<code>filukk</code> (a separate 17-band run) mismatched the 150-band NSCF evc, so $M^W(R_e;q{\\neq}0)$ "
     "came out flat. Fixed by re-Wannierizing on the 150-band NSCF (same Wannier space): $M^W$ now decays "
     "$\\sim\\!10^3\\times$, $\\tilde V^W$ is localized, and Koster&ndash;Slater converges by $R_{\\rm cut}{=}4$. "
     "The neutral defect is short-ranged &mdash; no supercell / range-separation issue.",
     '<a href="pages/results.html#sec-3">Open results &rarr;</a>'),
    ("Active-space dynamic resummation / transport", "Test", "&mdash;", "plan", "Planned",
     "Next (P6): frequency-dependent $T_{PP}(\\omega)$ on-shell, feeding $|T_{PP}|^2$ into the "
     "golden-rule rate in place of $|M|^2$ &mdash; beyond-Born vs Born mobility.", "&mdash;"),
]

def build_index():
    rows = ""
    for item, typ, date, bc, bl, summ, link in CATALOG:
        planned = ' class="planned"' if bc == "plan" else ""
        rows += ("<tr%s><td><strong>%s</strong></td><td>%s</td><td>%s</td>"
                 "<td><span class=\"badge %s\">%s</span></td><td>%s</td><td>%s</td></tr>\n"
                 % (planned, item, typ, date, bc, bl, summ, link))
    catalog = (
        '<section id="catalog"><h2>Test Catalog</h2>'
        '<p>One row per piece of work: what it is, when, status, the key result, and a link. '
        'Theory and plan are complete; the EDT package is implemented and the first '
        'numerical results are in &mdash; the downfolded potential $\\tilde V$ (P3) and the active-space '
        '$T$-matrix (P5) on the MoS&#8322; S-vacancy. See the '
        '<a href="pages/results.html">Numerical results</a> page.</p>'
        '<div class="table-wrap"><table><thead><tr>'
        '<th>Item</th><th>Type</th><th>Date</th><th>Status</th><th>Key result / summary</th><th>Link</th>'
        '</tr></thead><tbody>%s</tbody></table></div>'
        '<p class="small">Legend: '
        '<span class="badge ok">Complete</span> done &nbsp; '
        '<span class="badge plan">Planned</span> not yet run &nbsp; '
        '<span class="badge prod">Production</span> headline result.</p></section>'
        % rows)

    summary = (
        '<section id="summary"><h2>Executive Summary</h2>'
        '<div class="grid cards">'
        '<div class="card"><strong>Goal</strong><span>Electron&ndash;defect $T$-matrix beyond the '
        'first Born approximation, as a beyond-Born extension of the EDI workflow.</span></div>'
        '<div class="card"><strong>Layer 1 &mdash; rest</strong><span>Distant bands renormalize the '
        'bare defect potential $V$ into an effective $\\tilde V(\\omega_0)$ via a $k$-decoupled '
        'Sternheimer solve (Feshbach downfolding; exact to 2nd order in the active&ndash;rest '
        'coupling, iterable to all orders in $V_{QQ}$).</span></div>'
        '<div class="card"><strong>Layer 2 &mdash; active</strong><span>Full dynamical multiple '
        'scattering resummed by one small inversion '
        '$T_{PP}(\\omega)=[1-\\tilde V\\,G^A(\\omega)]^{-1}\\tilde V$.</span></div>'
        '<div class="card"><strong>Status</strong><span>Theory + plan complete; EDT implemented. '
        'First results: $\\tilde V$ block (P3) + active $T_{PP}$ (P5) on MoS&#8322;. '
        'Transport (P6) next. No raw data published.</span></div>'
        '</div>'
        '<p>The construction splits the host Green&rsquo;s function into an <strong>active</strong> '
        'block $A$ (bands near $E_F$, kept dynamical) and a <strong>rest</strong> block $R$ '
        '(distant bands, treated statically) &mdash; the constrained-RPA (cRPA) logic transported '
        'from screened interactions to multiple scattering. Neither layer ever forms or inverts a '
        'dense $(N_kN_b)\\times(N_kN_b)$ matrix.</p></section>')

    glance = (
        '<section id="glance"><h2>Method at a Glance</h2>'
        '<p>The two layers are cleanly separable and individually controlled:</p>'
        '<div class="math">\n$$\n'
        r'\underbrace{\;V\;\xrightarrow[\text{$k$-decoupled Sternheimer ladder}]'
        r'{\text{rest dressing, all orders in }V_{QQ}}\;\tilde V(\omega_0)\;}'
        r'_{\textbf{Layer 1: rest, statically renormalized}}'
        r'\;\xrightarrow[\text{exact active resummation}]{\quad\text{small inversion}\quad}\;'
        r'\underbrace{\;T_{PP}(\omega)=[1-\tilde V G^A(\omega)]^{-1}\tilde V\;}'
        r'_{\textbf{Layer 2: active dynamics}}'
        '\n$$\n</div>'
        '<p>See the <a href="pages/theory.html">Theory &amp; Method</a> page for the full '
        'derivation, the Sternheimer formulation, the $k$-decoupled ladder, convergence and '
        'failure modes, and the cRPA dictionary.</p></section>')

    warnings = (
        '<section id="not-published" class="warning"><h2>Warnings: Files Not Published</h2>'
        '<p>This GitHub Pages artifact contains only the static report under <code>docs/</code>. '
        'It intentionally excludes raw or large research data so the published site stays small '
        'and contains no credentials.</p>'
        '<div class="table-wrap"><table><thead><tr><th>Class</th><th>Policy</th></tr></thead><tbody>'
        '<tr><td>Wavefunctions, QE <code>*.save/</code></td><td>Excluded (large binary).</td></tr>'
        '<tr><td>Cube / volumetric grids <code>*.cube</code></td><td>Excluded (often &gt;100&nbsp;MB).</td></tr>'
        '<tr><td>Numerical arrays <code>*.npy/*.npz/*.bin/*.h5/*.dat</code></td><td>Excluded; summarized in tables only.</td></tr>'
        '<tr><td>Scheduler logs <code>*.out/*.err/*.log</code></td><td>Excluded.</td></tr>'
        '<tr><td>Local/private config <code>.claude/</code>, keys/tokens</td><td>Excluded.</td></tr>'
        '</tbody></table></div></section>')

    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">Active/rest partitioning, cRPA-style downfolding, and a '
              '$k$-decoupled Sternheimer solution for the electron&ndash;defect $T$-matrix '
              '&mdash; theory and (forthcoming) numerical tests.</p>'
              '<div class="meta">'
              '<span class="pill">Generated {d}</span>'
              '<span class="pill">GitHub Pages: docs/</span>'
              '<span class="pill">Branch: main</span>'
              '<span class="pill">MathJax v3</span></div></div></header>'
             ).format(t=SITE_TITLE, d=GEN_DATE)

    body = catalog + summary + glance + warnings
    out = page_shell(SITE_TITLE, header, _topnav("home"), body, "assets/style.css")
    with open(os.path.join(DOCS, "index.html"), "w", encoding="utf-8") as f:
        f.write(out)

# ======================================================================
#  main + self-checks
# ======================================================================
def main():
    os.makedirs(PAGES_DIR, exist_ok=True)
    build_theory()
    build_plan()
    build_note()
    build_results()
    build_index()

    th = open(os.path.join(PAGES_DIR, "theory.html"), encoding="utf-8").read()
    pl = open(os.path.join(PAGES_DIR, "plan.html"),   encoding="utf-8").read()
    nt = open(os.path.join(PAGES_DIR, "note-kprime-normalization.html"), encoding="utf-8").read()
    rs = open(os.path.join(PAGES_DIR, "results.html"), encoding="utf-8").read()
    ix = open(os.path.join(DOCS, "index.html"),       encoding="utf-8").read()

    def check(txt):
        problems = []
        if NUL in txt: problems.append("UNRESTORED placeholder (\\x00) present!")
        leftover = re.findall(r"%s[A-Z]\d+%s" % (NUL, NUL), txt)
        if leftover: problems.append("leftover tokens: %s" % leftover[:5])
        return problems

    def stats(name, txt):
        print("%-12s: %d bytes, %d sections, %d display-eq, %d tables, %d <pre>" % (
            name, len(txt), txt.count('<section id="sec-'),
            txt.count('class="math"'), txt.count("<table>"), txt.count("<pre>")))

    print("=== build_site.py ===")
    stats("theory.html", th)
    stats("plan.html",   pl)
    stats("note-knorm",  nt)
    stats("results.html", rs)
    print("results.html: %d <img>, %d tables" % (rs.count("<img "), rs.count("<table>")))
    print("index.html  : %d bytes, %d catalog rows" % (len(ix), ix.count("<tr")))
    for nm, txt in (("theory.html", th), ("plan.html", pl), ("note-knorm", nt),
                    ("results.html", rs), ("index.html", ix)):
        p = check(txt)
        print("  [%s] %s" % (nm, "OK" if not p else " ; ".join(p)))

if __name__ == "__main__":
    main()
