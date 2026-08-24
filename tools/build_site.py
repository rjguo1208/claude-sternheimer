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
NOTE_NPOL2 = os.path.join(ROOT, "content", "npol2cert.md")
NOTE_SOCZOOM = os.path.join(ROOT, "content", "soczoom.md")
NOTE_MODED   = os.path.join(ROOT, "content", "moded.md")
NOTE_PDCOO2  = os.path.join(ROOT, "content", "pdcoo2.md")
NOTE_PDVAC = os.path.join(ROOT, "content", "pdvac-coset.md")
NOTE_PTCOO2  = os.path.join(ROOT, "content", "ptcoo2.md")
NOTE_RESULTS = os.path.join(ROOT, "content", "note_tmatrix_results.md")
NOTE_LADDER = os.path.join(ROOT, "content", "note_sternheimer_ladder.md")
NOTE_KOSTER = os.path.join(ROOT, "content", "note_koster_slater.md")
NOTE_DEFLATED = os.path.join(ROOT, "content", "note_deflated_ladder.md")
NOTE_TLRES = os.path.join(ROOT, "content", "note_twolevel_results.md")
NOTE_METH = os.path.join(ROOT, "content", "note_methods.md")
NOTE_FESH = os.path.join(ROOT, "content", "note_feshbach_plan.md")
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
    # joined rather than "%s"*N so adding a page cannot desync the placeholder count
    return ('<nav class="topnav"><div class="inner">'
            '<span class="brand">Sternheimer&nbsp;EDI</span>'
            + "".join((a("index.html", "Home", "home"),
                                      a("pages/theory.html", "Theory &amp; Method", "theory"),
                                      a("pages/sternheimer-ladder.html", "Rest-space ladder", "ladder"),
                                      a("pages/deflated-ladder.html", "Deflated ladder", "deflated"),
                                      a("pages/twolevel-results.html", "Two-level results", "tlres"),
                                      a("pages/methods.html", "Method derivations", "methods"),
                                      a("pages/koster-slater.html", "Koster&ndash;Slater", "koster"),
                                      a("pages/feshbach-implementation.html", "Feshbach impl.", "fesh"),
                                      a("pages/plan.html", "Implementation Plan", "plan"),
                                      a("pages/results.html", "Results", "results"),
                                      a("pages/npol2-cert.html", "npol=2 cert.", "npol2"),
                                      a("pages/soc-zoom.html", "SOC zooms", "soczoom"),
                                      a("pages/mos2-mobility.html", "MoS&#8322; mobility", "mos2"),
                                      a("pages/mode-d-crosscheck.html", "MODE D check", "moded"),
                                      a("pages/pdcoo2-bands.html", "PdCoO2 bands", "pdcoo2"),
                                      a("pages/pdvac-coset.html", "Pd-vac 3D", "pdvac"),
                                      a("pages/ptcoo2.html", "PtCoO2 chain", "ptcoo2"),
                                      a("pages/note-kprime-normalization.html", "Note: k&prime;-norm", "note")))
            + '</div></nav>')

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

def build_npol2cert():
    with open(NOTE_NPOL2, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">How EDT&rsquo;s npol=2 (spinor) implementation was acquitted by three '
              'gauge-invariant criteria after the doubling test failed, how the EDI-direct noncolin matrix '
              'elements were convicted, and the EDT-produced born block and two-gate verdict that replaced '
              'them.</p>'
              '<div class="meta"><span class="pill">Certification</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — npol=2 certification",
                     header, _topnav("npol2", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "npol2-cert.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_soczoom():
    with open(NOTE_SOCZOOM, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">The first electron&ndash;defect spectral functions including SOC: '
              'K-valley CBM/VBM zooms, a 10-spinor-band active manifold, EDI-v8d born + an EDT r11 chain, '
              '$n_d=10^{{12}}$ cm$^{{-2}}$.</p>'
              '<div class="meta"><span class="pill">Production</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — SOC zoom spectral functions",
                     header, _topnav("soczoom", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "soc-zoom.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_moded():
    with open(NOTE_MODED, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=True)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">{s}</p>'
              '<div class="meta"><span class="pill">Cross-check</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], s=r.get("subtitle") or "", n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " \u2014 MODE D cross-check",
                     header, _topnav("moded", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "mode-d-crosscheck.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r


def build_pdvac():
    with open(NOTE_PDVAC, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=True)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">{s}</p>'
              '<div class="meta"><span class="pill">3D validation</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], s=r.get("subtitle") or "", n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " \u2014 PdCoO2 Pd-vacancy coset validation",
                     header, _topnav("pdvac", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "pdvac-coset.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r


def build_ptcoo2():
    with open(NOTE_PTCOO2, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=True)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">{s}</p>'
              '<div class="meta"><span class="pill">Production</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], s=r.get("subtitle") or "", n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — PtCoO2 electron-defect T-matrix chain",
                     header, _topnav("ptcoo2", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "ptcoo2.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r


def build_pdcoo2():
    with open(NOTE_PDCOO2, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=True)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">{s}</p>'
              '<div class="meta"><span class="pill">Standalone DFT</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], s=r.get("subtitle") or "", n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " \u2014 PdCoO2 band structure",
                     header, _topnav("pdcoo2", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "pdcoo2-bands.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r


def build_ladder():
    with open(NOTE_LADDER, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">Why the 2nd-order (Born) rest-space self-energy over-screens the '
              'MoS&#8322; S-vacancy $e$ level by $\\sim$1 eV, and the order-by-order Sternheimer ladder '
              '($Q(\\omega-H_0)Q$ reused each order; one solve buys two orders, the odd order free) '
              'that systematically corrects it.</p>'
              '<div class="meta"><span class="pill">Method / derivation</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Rest-space Sternheimer ladder",
                     header, _topnav("ladder", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "sternheimer-ladder.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_deflated():
    with open(NOTE_DEFLATED, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">The production formulation of the dressed rest: Schur-eliminate the '
              'explicit $\\le$150-band block into an exact all-order inverse $D_1$, leaving only '
              'tail-passing vertices $W_{{22}}$ for a short Neumann/GMRES ladder — plain iteration '
              'diverges ($\\rho=1.215$), the deflated one contracts at $\\rho\\approx0.15$.</p>'
              '<div class="meta"><span class="pill">Method / derivation</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Deflated ladder (explicit + tail)",
                     header, _topnav("deflated", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "deflated-ladder.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_tlres():
    with open(NOTE_TLRES, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">Level-by-level truth gates for the two-level (deflated) dressing on the '
              '6$\\times$6 grid: the Sdisp probe lands the CBM doublet at 2 meV, and the relaxed V$_S$ gap '
              '$e$-doublet — historically $+217$ meV (truncated) / $+41$ meV (best $\\chi$) — comes out at '
              '$+10/+39$ meV with no augmentation and no hyperparameters.</p>'
              '<div class="meta"><span class="pill">Results</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Two-level dressing results",
                     header, _topnav("tlres", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "twolevel-results.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_methods():
    with open(NOTE_METH, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">One page, all the formulas: the exact downfolding identity and each '
              'method as one choice of the rest self-energy $\\Sigma$ &mdash; M-only, bare second order, '
              'the two-level static ladder, the $\\omega$-resolved block-Lanczos continued fraction, and '
              '$\\chi$ augmentation.</p>'
              '<div class="meta"><span class="pill">Method</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Method derivations",
                     header, _topnav("methods", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "methods.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_koster():
    with open(NOTE_KOSTER, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">Defect levels as roots of $\\det_D[1-g(E)\\Delta V_D]=0$ in the small '
              'defect-localized block &mdash; the exact, all-bands host-Green&rsquo;s-function route '
              '(Koster&ndash;Slater), sidestepping the $P$/$Q$ rest dressing entirely; $\\sim$seconds vs '
              'explicit&rsquo;s 46 min / Feshbach&rsquo;s days.</p>'
              '<div class="meta"><span class="pill">Method / derivation</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Koster–Slater defect Green's function",
                     header, _topnav("koster", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "koster-slater.html"), "w", encoding="utf-8") as f:
        f.write(out)
    return r

def build_feshplan():
    with open(NOTE_FESH, encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">Static-$\\omega_0$ full-order rest dressing as ONE all-$k$ Sternheimer '
              'solve $[Q(H_0{{+}}\\Delta V{{-}}\\omega_0)Q{{+}}\\alpha P]\\,X_b=Q\\Delta V|b\\rangle$. '
              'v1 post-mortem (the isolated check caught a $1/N_k^2$ convention bug), the corrected '
              'cube-anchored $\\Delta V$ applier, an element-wise validation ladder anchored on the trusted '
              '$M$ block, a zero-communication source-parallel layout, and measured-anchor cost.</p>'
              '<div class="meta"><span class="pill">Implementation plan (v2)</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(SITE_TITLE + " — Full-order Feshbach implementation plan",
                     header, _topnav("fesh", prefix="../"), body, "../assets/style.css")
    with open(os.path.join(PAGES_DIR, "feshbach-implementation.html"), "w", encoding="utf-8") as f:
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
    ('MoS$_2$ defect-limited mobility: a four-level scattering-theory comparison + K-valley zoom spectral functions', 'Production',
     '2026-08-23', 'prod', 'Headline',
     # hand-inserted into docs/index.html before this was registered here,
     # which is why a plain rebuild of the site silently dropped it once
     'A 12-cell matrix of {full $T$, Born} &times; {$\\Sigma$ kept, dropped} on one and the same 144b22 chain '
     '(three defects &times; e/h, 300 K $\\eta\\to0$ IBTE, $n=10^{10}$ and $n_d=10^{12}$ cm$^{-2}$). '
     '<b>Born is wrong by 1&ndash;2 orders of magnitude for deep wells</b> (V$_\\mathrm{S}$ holes 53&times;, '
     'strict second order 190&times;); <b>$\\Sigma$&rsquo;s role flips with the defect</b>: it screens for '
     'V$_\\mathrm{S}$ (19.9&rarr;13.3 eV), builds a CBM resonance for O$_\\mathrm{S}$ ($\\Sigma$=42 eV, 4&times; '
     'the bare vertex, and dropping $\\Sigma$ inflates $\\mu_e$ by a spurious 72&times;), and nearly cancels for '
     'Se$_\\mathrm{S}$ (0.3&ndash;0.6 eV residual); only full $T$+$\\Sigma$ is consistent across the board. '
     'Also: K-valley zoom $A(k,\\omega)$ (O$_\\mathrm{S}$ CBM spectral weight pulled into the gap), '
     '<b>Wannier-decay validation of the edmat and downfolded matrices</b> (all 15 objects LOCALIZED, '
     'peak/floor 10&#8308;&ndash;10&#8311;; chain VAA &equiv; the edmat active block, bit-for-bit 0.0), and a '
     '2.2&times; vertex correction for Se$_\\mathrm{S}$ holes.',
     '<a href="pages/mos2-mobility.html">Open the mobility page &rarr;</a>'),
    ("PtCoO$_2$ full chain ($S0\\to S11$): is that factor of 2 a method systematic, or material-dependent?", "Production",
     "2026-08-23", "prod", "Headline",
     "The PdCoO$_2$ chain <b>ported verbatim to PtCoO$_2$</b> (same code, same production "
     "parameters, structures from the same refinement) to settle one discriminating question. "
     "<b>Answer: not consistent. Full-order $T$/IBTE gives a Frenkel slope of 13021 nΩ·cm/% "
     "= 1.235× the measured 10546, against 2.05× for PdCoO$_2$</b> \u2014 and the two "
     "<b>measured</b> slopes differ by only 1% (10546 vs 10654), so the whole difference sits "
     "in the theory: the factor of 2 is <b>material-dependent</b>, and PtCoO$_2$ actually "
     "agrees better. <b>The mechanism is identified too</b>: at identical production settings "
     "Pt/Pd = <b>0.593</b>, and <b>the two arms give that ratio independently</b> "
     "(vacancy 0.571 / interstitial 0.605) \u2014 a material-level effect, not an accident of "
     "defect configuration. With $\\Sigma$ <b>dropped, Pt/Pd rises to 0.83\u20130.87</b>: the "
     "materials differ by only ~15% in bare scattering but 41% once the rest-space self-energy "
     "is kept, so <b>$\\Sigma$ carries the difference</b> \u2014 Pt\u2019s 5d/6s leave the rest "
     "space 27% closer to $E_F$ (2.425 vs 3.311 eV), a prediction written down the moment S1 "
     "measured band 27 and later confirmed four independent ways. Structure fixed by a "
     "<b>Fermi-surface topology test</b> (both alternative parameter sets give two crossing "
     "bands, while PtCoO$_2$ is experimentally a single-band metal); S1 hard gate = manifold "
     "isolation 0.2295 eV (4.1× narrower than PdCoO$_2$ but positive, so no disentanglement). "
     "The <b>ordering reversal</b> the source document asks for does appear (interstitial 2.0× "
     "stronger at full order, vacancy 6.5× stronger at Born). <b>Convergence: 9 axes / 37 "
     "points, 18/18 byte-identical to production</b> from a different job script; a conservative "
     "±9.2% band gives 1.12–1.35×, <b>disjoint</b> from PdCoO$_2$\u2019s 2.05–2.06×. "
     "Counter-intuitive find: <b>n_lancz=24 is required for the interstitial arm</b> \u2014 the "
     "vacancy is within 0.42% at 12 steps, the interstitial is off by <b>97%</b> (dirty kernel, "
     "asym 0.41 = 5× the vacancy). Reported honestly: <b>the $G_0$ axis does not converge</b> "
     "(±20%, driven by the optical-theorem factor $s$ swinging 0.76\u20131.28), RCUT has only "
     "2 points, and the Fermi-surface topology is near-critical (band 23 top just 35 meV below "
     "$E_F$).",
     '<a href="pages/ptcoo2.html">Open the PtCoO$_2$ chain &rarr;</a>'),
    ("PdCoO$_2$ Pd vacancy: end-to-end validation of the 3D $T$-matrix pipeline (2&times;2&times;2 coset vs supercell ground truth)", "3D validation",
     "2026-08-14", "prod", "Verified",
     "The <b>first three-dimensional system the EDT pipeline has run</b>. Supercell &Delta;V (extract_pot, five "
     "orders of magnitude of decay within 4 &Aring;) &rarr; core alignment (two independent methods agreeing to "
     "8 meV) &rarr; all-band born (arbitrated by an independent plane-wave integrator) &rarr; MODE C chain "
     "(herm ~6e-15) &rarr; coset DOS reconciled against the vacancy supercell&rsquo;s &Gamma; ground truth: "
     "<b>integrated deviation 0.37% (gate-lane noise floor 0.13%), median QP level error 0.1 meV, and the 5 "
     "unpaired levels exactly matching the removed Pd&rsquo;s 4d count</b>. Along the way, a q-block fingerprint "
     "localized and fixed a hard 2D-only bug in fold_col (a q3 slot collision plus a missing g3 winding phase; "
     "qe-edt v2.2): the asymmetric pattern of exactly 16 q3=+1/2 blocks coming out precisely +1.000 pointed "
     "uniquely at three lines of code. 2D runs are bit-for-bit unaffected. The <b>6&sup3; production chain now "
     "runs</b> (SVD rank 288/3456=8.3%, herm ~3e-15), with <b>Wannier decay validation</b> complete: V/T decay "
     "two orders of magnitude within the first shell (the cluster method&rsquo;s licence to operate), the "
     "box-corner systematic is quantified, and a <b>&times;5 local resonance in T at E$_F$</b> was discovered. "
     "<b>Spectral functions landed</b>: A(k,&omega;) along the RHL1 path (gate-0 identity = N_k to 2.9e-13), "
     "&Gamma;(E$_F$) = 24 meV at 0.5% vacancies &rarr; &tau; = 27 fs and &ell; ~ 20 nm &mdash; a resonantly "
     "enhanced strong scatterer. <b>Transport comparison</b>: the SERTA residual-resistivity slope of "
     "4.4&ndash;4.7&times;10&sup3; n&Omega;&middot;cm/% is <b>about half the unitary limit (9.2&times;10&sup3;, "
     "with our own n) and half the irradiation experiment (9&ndash;10&times;10&sup3;)</b>, reproducing with zero "
     "free parameters the conclusion that the experiment sits at the unitary limit.",
     '<a href="pages/pdvac-coset.html">Open the 3D validation &rarr;</a>'),
    ("PdCoO$_2$ delafossite bands: PBE vs PBE+U, cross-validated in two cells", "Standalone DFT",
     "2026-08-14", "ok", "Verified",
     "A standalone calculation starting from one POSCAR (independent of this site&rsquo;s defect T-matrix main "
     "line). Two paths: the rhombohedral primitive cell (4 atoms, RHL1) and a programmatically constructed "
     "hexagonal conventional cell (12 atoms, &Gamma;-M-K-&Gamma;-A-L-H-A); <b>only one Pd band crosses "
     "E$_F$</b> (folding into 3 in the hexagonal cell, self-consistently), and it is almost perfectly flat along "
     "&Gamma;&rarr;A &mdash; the band-structure explanation for the extreme conductivity and the "
     "quasi-two-dimensionality. <b>+U (Co-3d, 4 eV) only moves the Co-3d occupied bands at &minus;1.5..&minus;4 eV, "
     "shifting E$_F$ by a mere 3 meV</b>: the transport physics is insensitive to U. Total energies per f.u. agree "
     "between the two cells to 10 meV (0.37 mRy with +U). <b>DFPT linear response gives U(Co-3d) = 7.95 eV</b> "
     "(screening ratio &chi;/&chi;&#8320; = 0.135, so Co-3d really is localized), but a spin-polarized check finds "
     "the moment strictly 0 (collapsing back to the non-magnetic solution, a 0.1 &mu;eV difference) &mdash; a "
     "closed shell means +U only gives a rigid shift, <b>so transport uses PBE and only spectroscopy uses "
     "U=7.95</b>. <b>16 MLWFs</b> reproduce the entire isolated p-d manifold (no disentanglement): two different "
     "projection sets converge to the same global minimum (&Omega;_I identical to 9 digits), with a point-by-point "
     "RMS against the DFT bands of <b>3.6 meV</b> and 2.8 meV near the Fermi level.",
     '<a href="pages/pdcoo2-bands.html">Open the PdCoO$_2$ bands &rarr;</a>'),
    ("MODE D (folded free-electron tail) cross-check: DOS reconciliation against MODE C + free-electron energy-space convergence", "Cross-check",
     "2026-08-13", "ok", "Verified",
     "An independent reproduction and reconciliation of the OPW / folded free-electron tail contributed by our "
     "collaborator cz (EDI `opw-tail`), on <b>three systems</b> in 6&times;6 non-SOC (ideal V$_S$, relaxed V$_S$, "
     "relaxed O$_S$). The M2 numbers reproduce one by one (drop 11.10 / ffree 1.32 / fexact 2.79 meV); <b>all 18 "
     "DOS peaks across the three systems shift by &le;5 meV relative to MODE C</b> (the mesh resolution), whereas "
     "drop-block, which discards the same tail, pushes defect levels by 10&ndash;90 meV; the in-gap $L_1$ error "
     "falls to 1/4.2&ndash;1/6.2. A newly swept free-electron energy-space convergence shows production &chi; needs "
     "only ~2500 plane waves (a 300 eV cutoff, costing 0.007&ndash;0.024 meV). <b>The production form now works "
     "end to end</b> (EDT gained a chi export, with the tail being the true complement up to the plane-wave "
     "cutoff): both MODE C and MODE D pass the phantom gate (O$_S$&rsquo;s gap is empty, matching ground truth), "
     "but MODE D puts V$_S$&rsquo;s deep doublet <b>100 meV</b> above MODE C &mdash; so the M2 self-consistency "
     "test (which reports 8 meV at the same $n_x$) systematically overstates the method&rsquo;s accuracy.",
     '<a href="pages/mode-d-crosscheck.html">Open the cross-check &rarr;</a>'),
    ("22-band 144k dilute-limit DOS: a three-way comparison over three defects (truth / array / dilute)", "Production", "2026-08-12", "prod", "Headline",
     "Dilute-limit DOS from a 12&times;12&times;22-spinor-band chain (start-block SVD `svd_tol=1e-4`, rank "
     "~1024&ndash;1080 = 32&ndash;34%, 4.6 h per chain on 4 Banff nodes). V$_S$&rsquo;s deep in-gap doublet: "
     "truth +1.085/+1.135 &rarr; array +1.115/+1.165 &rarr; dilute +1.120/+1.170 &mdash; <b>a dilution shift of "
     "&le;5 meV, so the deep states have already reached the isolated limit at 1/36</b>; O$_S$&rsquo;s CB "
     "resonance at +1.545 coincides across all three; and Se$_S$ gets its first downfolded DOS (a clean gap). "
     "Three SVD scales bracket the true rank in (1008, 1024], so the physical rank is set by the defect "
     "potential&rsquo;s spatial support and is independent of the defect chemistry.",
     '<a href="pages/soc-zoom.html#sec-4">Open the 3-way DOS &rarr;</a>'),
    ("First batch of SOC defect spectral functions: K-valley CBM/VBM zooms (O$_S$ + V$_S$)", "Production", "2026-08-10", "prod", "Headline",
     "The first $A(k,\\omega)$ including spin&ndash;orbit coupling: 12&times;12, a 10-spinor-band active manifold, "
     "$n_d=10^{12}$ cm$^{-2}$, $\\eta=5$ meV. The CBM window shows the conduction-band SOC doublet directly "
     "(3 meV splitting at K, opening along $K\\to M$); the VBM window shows the 149 meV K-point valence splitting "
     "together with the defect flat-band multiplet (a new +0.181 eV feature on V$_S$&rsquo;s gap side). "
     "born = EDI-v8d (cross-code 1e-14), chain = EDT r11 (unit test 2.5&ndash;2.7e-11), and the 10-band cluster "
     "brings the T cache down from hours to 5 minutes. Se$_S$ is running.",
     '<a href="pages/soc-zoom.html">Open the zooms &rarr;</a>'),
    ("npol=2 certification: the doubling-test verdict, EDI noncolin convicted, EDT&rsquo;s born block made self-sufficient", "Verification", "2026-08-09", "ok", "Certified",
     "A forensic adjudication triggered by a FAILing doubling test: the local arm matches an independent numpy "
     "referee at exactly 1.000 on 52/52 elements; degenerate-pair 2-vector parallelism reaches $|\\cos|=1.0000$ on "
     "25/25; and the scalar regression is $8\\times10^{-15}$ &mdash; EDT is innocent. The EDI-direct noncolin matrix "
     "elements were convicted (spin components lost under mixed spinors, with residuals reaching 0.50). After "
     "porting the born block into EDT: closure unit test $4.6\\times10^{-11}$, an SV audit over 1260 blocks with "
     "p90 $7.8\\times10^{-3}$, and a final verdict of 0.0000 meV spinor-pair degeneracy. Lesson: the median SV can "
     "be fooled by the $M_\\mathrm{loc}\\otimes I_2+M_\\mathrm{KB}\\otimes e_{11}$ structure &mdash; always look at "
     "the whole distribution.",
     '<a href="pages/npol2-cert.html">Open the adjudication &rarr;</a>'),
    ("Method derivations: $\\Sigma$, the $T$-matrix, the e&ndash;d self-energy", "Method", "2026-08-07", "ok", "Reference",
     "The exact downfolding identity and each method as one choice of $\\Sigma$: M-only / bare 2nd order / "
     "two-level static ladder (with its $\\rho(D_2W_{22})<1$ caveat) / $\\omega$-resolved block-Lanczos "
     "continued fraction / $\\chi$ augmentation; plus the fold-grid requirement and the route from the "
     "downfolded vertex to the $T$-matrix, the spectral function and the electron&ndash;defect self-energy.",
     '<a href="pages/methods.html">Open derivations &rarr;</a>'),
    ("Full reorthogonalization is not needed: 531 GB &rarr; 62 GB", "Method", "2026-08-09", "ok", "Measured",
     "Storing every Krylov block is what caps the chain length and blocks SOC, so the necessity of the sweep "
     "was measured rather than assumed. Orthogonalizing against only the last 2 blocks instead of all: "
     "orthogonality IDENTICAL to four figures ($6.038\\times10^{-9}$), bit-identical levels, chain matrices "
     "differing $3\\times10^{-14}$ on a scale of 31, zero Ritz ghosts in the gap. The reason is Paige&#39;s: "
     "loss of orthogonality needs a CONVERGED Ritz value, and 24 blocks of width 396 in $10^6$ dimensions "
     "converge none. Speed only 1.65x, but $Q_s$ drops from $N_S+1$ blocks to two &mdash; and an SOC "
     "12x12 run from eight nodes to one.",
     '<a href="pages/methods.html#sec-9">Open the test &rarr;</a>'),
    ("Methods manuscript (REVTeX / APS): the whole route, publication form", "Writeup", "2026-08-08", "ok", "Drafted",
     "The complete method chain written as a PRX-format Methods section: Feshbach downfolding, the fold and its "
     "commensurability condition, the free $k$-mesh, block-Lanczos $\\Sigma^R(\\omega)$, the pair-Wannier "
     "transform and its origin convention, the real-space Koster&ndash;Slater closure, and a validation table "
     "gating every stage. 6 pages, compiles clean under REVTeX 4.2.",
     '<a href="paper/methods.pdf">methods.pdf &rarr;</a> &nbsp; <a href="https://github.com/rjguo1208/claude-sternheimer/blob/main/paper/methods.tex">source</a>'),
    ("Spectral function A(k,$\\omega$) end to end, welded to the coarse-grid levels", "Result", "2026-08-08", "prod", "Headline",
     "The whole route run end to end: $\\omega$-resolved vertex &rarr; pair-Wannier &rarr; real-space "
     "Koster&ndash;Slater cluster &rarr; $\\Sigma^{ed}=n_d T_{kk}$ &rarr; $A(k,\\omega)$ on "
     "$\\Gamma$&ndash;M&ndash;K&ndash;$\\Gamma$, with NO fit in $\\omega$. Acceptance: the deep in-gap "
     "bound state, a root of a determinant with no $k$ in it, comes out at $+1.2018$ eV against $+1.2023$ "
     "from direct diagonalization &mdash; 0.5 meV. That welds Wannier gauge, pair kernel, cluster truncation "
     "and fine-grid $G^A$ onto the validated coarse-grid result. 40.7 min, 90% of it the continued fraction.",
     '<a href="pages/twolevel-results.html#sec-10">Open the spectral function &rarr;</a>'),
    ("$\\omega$-dependent vertex + active-space size: three acceptance tests", "Method", "2026-08-08", "ok", "Passed",
     "$\\omega_0$ sits BELOW the gap, so freezing $\\Sigma^R$ there costs 18&ndash;24 meV on the in-gap levels "
     "&mdash; more than every other error in the chain. Unnecessary: $\\Sigma^R$ is analytic across the gap "
     "(nearest pole &gt;3 eV), so 7 Chebyshev nodes + a 4th-order fit give 0.026 meV on a held-out node and "
     "0.201 meV on the quasiparticle fixed points. Leave-one-out on $\\tilde V(\\omega)$ passes at every "
     "$\\omega$ (0.61&ndash;0.88 meV). And the transform is LINEAR, so the $\\omega$ dependence rides on 5 "
     "scalar coefficients &mdash; one set of transforms, not one per $\\omega$. Separately: a 5-band Wannier "
     "space is 4.4x worse at interpolation ($\\Omega_I$ per WF 4.36 vs 1.61 &Aring;$^2$, gauge-invariant), so "
     "keep 11 bands on this route.",
     '<a href="pages/methods.html#sec-9">Open the tests &rarr;</a>'),
    ("Wannier interpolation of the defect vertex: leave-one-out gate", "Method", "2026-08-08", "ok", "Passed",
     "The $N_k^3$ wall is escapable only if the pair kernel $\\mathcal M(R_e,R_p)$ is local. Measured: it falls "
     "from 1.46 on site to $2\\times10^{-4}$ at $6a$, and a kernel built from the $\\Gamma$ coset alone (6.2% "
     "of pairs) predicts the other 93.8% to a MEDIAN 0.74 meV &mdash; below both the 7.6 meV downfold residual "
     "and the 4&ndash;11 meV image error of the reference. Both $R$ indices are ELECTRON Wannier positions "
     "(EDI, ft_convention.md and this project&#39;s earlier tmatrix_p6_wannier.py agree, up to a shear); the "
     "defect is only the origin they decay about, and getting that origin wrong on a finer lattice reads "
     "300 meV instead of 1.",
     '<a href="pages/methods.html#sec-9">Open the test &rarr;</a>'),
    ("Downfolding RESULTS: Sdisp, V$_S$, O$_S$, Se$_S$ truth gates (rev 4)", "Result", "2026-08-08", "prod", "Headline",
     "REVISION 4 adds a second, independent gate: splitting a 12x12 $k$-grid into its four cosets gives "
     "four shifted 6x6 runs at the SAME defect concentration, so the downfold can be tested at the "
     "supercell $M$ point against a dedicated 107-atom NSCF it has never seen. Raw offset $+25.0\\pm0.8$ "
     "meV at $M$ vs $+25.5\\pm1.7$ at $\\Gamma$ &mdash; the reference constant is a property of "
     "$\\Delta V$, not of $k$. Run directly instead, the same grid is a four-times more dilute defect "
     "with full $q$ resolution (the supercell is a box holding one isolated defect, not a periodicity, "
     "so the $k$-grid is free): that run puts the 6x6 cell\'s OWN periodic-image error at 4&ndash;11 meV, "
     "i.e. the downfold residual is no longer the largest error in the comparison. "
     "NET RESULT: $\\omega$-resolved MODE C matches the supercell truth on ALL three defects to $\\le$8 meV "
     "absolute and 0.7&ndash;3.8 meV in shape, with no fitting and no basis engineering.",
     '<a href="pages/twolevel-results.html">Open results &rarr;</a>'),
    ("Deflated ladder: explicit all-order rest ($\\le$150) + Sternheimer tail", "Method", "2026-08-05", "ok", "Derived &amp; certified",
     "Schur elimination locks the coherent low rest into an exact inverse $D_1$; only tail vertices "
     "$W_{22}=V_{22}+V_{21}D_1V_{12}$ iterate. Plain ladder diverges ($\\rho=1.215$, measured); deflated "
     "contracts $\\rho=0.574\\to0.323\\to0.153$ at splits $N_1=90/120/140$ &rArr; physical tail "
     "$\\rho\\approx0.1$&ndash;$0.2$, 4&ndash;5 rungs. Stage-B production formula.",
     '<a href="pages/deflated-ladder.html">Open derivation &rarr;</a>'),
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
    ("Defect levels recovered: explicit 21-band $T$-matrix", "Result", "2026-06-08", "prod", "Complete",
     "Benchmarked against the S-vacancy supercell (&Gamma;-point KS levels): the active-space + Born-$\\Sigma$ "
     "pipeline keeps the $a_1$ but <em>over-screens</em> the $e$ (drags it to the VBM). Keeping the conduction "
     "manifold <em>explicitly</em> (raw Bloch $M=\\langle nk|\\Delta V|mk'\\rangle$, $12\\times12$, no Wannier / no "
     "downfolding) recovers the full $C_{3v}$ pattern, and the $e$ is now <strong>band-converged onto DFT</strong>: "
     "$+1.495$ (11 bands) $\\to$ $+1.348$ (21) $\\to$ $+1.205$ (61, $N_A{=}8758$) vs DFT $+1.19$ &mdash; within "
     "$15$ meV. The 2nd-order rest dressing, not the active space, was the culprit.",
     '<a href="pages/results.html#sec-3">Open results &rarr;</a>'),
    ("Rest-space Sternheimer ladder: beyond 2nd order", "Theory", GEN_DATE, "ok", "Derived",
     "Why the 2nd-order (Born) rest dressing over-screens &mdash; Born series $\\rho\\sim\\mathcal O(1)$ for the deep "
     "vacancy, so it over-binds the $e$ by $\\sim$1 eV &mdash; and the order-by-order Sternheimer ladder that fixes it: "
     "$Q(\\omega\\!-\\!H_0)Q$ reused each order, $\\Sigma^{(2+p)}_{ab}=\\langle\\chi^i_a|\\Delta V_{QQ}|\\chi^j_b\\rangle$ "
     "with $i+j=p-1$ (one solve buys two orders; $\\Sigma^{(3)}$ free). Full Feshbach is the resummed endpoint. "
     "Derivation + nonlocal implementation plan; not yet coded.",
     '<a href="pages/sternheimer-ladder.html">Open derivation &rarr;</a>'),
    ("Koster&ndash;Slater defect Green&rsquo;s function: efficient defect levels", "Theory", GEN_DATE, "ok", "Derived",
     "Defect levels as roots of $\\det_D[1-g(E)\\Delta V_D]=0$ in the small (~tens-orbital) defect-localized block. "
     "The host GF $g(E)$ carries the bands (cheap eigenvalue $k$-sum, fine grid); $\\Delta V_D$ is short-ranged "
     "(Fig 17&ndash;20). Exact &mdash; no rest dressing, so no over-screening and no ladder divergence &mdash; and "
     "$\\sim$seconds vs explicit&rsquo;s 46 min / Feshbach&rsquo;s days. $C_{3v}$ blocks give rigorous $a_1$/$e$ "
     "labels; Krein&ndash;Friedel gives the $\\Delta$DOS. Groundwork (gauge-fixed $M^W$, $R_{\\rm cut}{=}4$) in place; not yet run.",
     '<a href="pages/koster-slater.html">Open derivation &rarr;</a>'),
    ("Full-order Feshbach Sternheimer: implementation plan (v2)", "Plan", GEN_DATE, "ok", "Drafted",
     "Static-$\\omega_0$ full-order rest dressing as ONE all-$k$ Sternheimer solve "
     "$[Q(H_0{+}\\Delta V{-}\\omega_0)Q{+}\\alpha P]X_b=Q\\Delta V|b\\rangle$. Piece-B v1's isolated check caught a "
     "$1/N_k^2$ convention bug &rarr; v2 fixes the three root causes (cube-anchored unfold, 36 cells + intra-cell "
     "Bloch phases; $\\Delta V$ field straight from <code>V_colin</code>), validates element-wise against the "
     "trusted $M$ block (V0), and runs source-parallel with a zero-communication matvec. "
     "Measured-anchor cost $\\sim$8.6 h/node, $\\div N$ nodes. Implementation pending (P-I&ndash;P-IV).",
     '<a href="pages/feshbach-implementation.html">Open plan &rarr;</a>'),
    ("QE-Hamiltonian Sternheimer validation", "Test", "2026-05-31", "ok", "Complete",
     "Per-$k$ solve of $Q(\\omega_0-H_0)Q$ via projected PCG (QE <code>h_psi</code> matvec): "
     "$\\langle\\psi|H_0|\\psi\\rangle\\!=\\!\\varepsilon$ gate to $6\\times10^{-10}$ eV across all ranks; "
     "explicit rest-band sum converges to the all-band Sternheimer value (Born limit $T\\!\\to\\!V$ at $10^{-13}$ Ry).",
     '<a href="pages/plan.html">see P0&ndash;P3 in the plan &rarr;</a>'),
    ("Rest dressing ladder convergence", "Test", "2026-06-11", "ok", "Measured",
     "Successive ratios MEASURED across rest sizes: $r_3{=}0.37$ ($Q$=18&ndash;28), "
     "$r_3{=}0.74,\\ r_4{=}0.91,\\ r_5{=}1.21$ ($Q$=18&ndash;70, crossing 1), $r_3{=}1.54$ (full rest, in-code) "
     "&mdash; the Born series for the rest dressing is <em>divergent</em>, and the partial sums oscillate "
     "($e$: $+1.50\\to+0.73\\to+1.60$). Only the full resummation lands correctly.",
     '<a href="pages/sternheimer-ladder.html">See the ladder page &rarr;</a>'),
    ("Full-order rest dressing: direct resolvent from explicit-60", "Result", "2026-06-11", "prod", "Complete",
     "The full-order static-$\\omega_0$ Feshbach dressing of the 11-band active space, computed by DIRECT "
     "resolvent inversion on the explicit-60 data ($Q$=bands 18&ndash;70, one $H_{QQ}$ eigendecomposition "
     "&mdash; no Sternheimer iteration). <strong>Self-consistent full-order $e$ = $+1.205$ = the explicit-60 "
     "all-band value, to the displayed digits</strong> (DFT $+1.19$): the downfolding identity closes exactly, "
     "and the over-screening of the 2nd order ($+0.73$) is cured. Dressed $\\tilde V$ written in the standard "
     "block format; full-order multiband spectral function on the results page (Fig 17 vs Fig 10).",
     '<a href="pages/results.html#sec-3">Open results &rarr;</a>'),
    ("Wannier representation &amp; Koster&ndash;Slater (P5-b)", "Result", "2026-06-01", "ok", "Resolved",
     "Wannierizing $\\tilde V$ needs $U(k)$ in the same Bloch gauge as the evc that build $M$. The original "
     "<code>filukk</code> (a separate 17-band run) mismatched the 150-band NSCF evc, so $M^W(R_e;q{\\neq}0)$ "
     "came out flat. Fixed by re-Wannierizing on the 150-band NSCF (same Wannier space): $M^W$ now decays "
     "$\\sim\\!10^3\\times$, $\\tilde V^W$ is localized, and Koster&ndash;Slater converges by $R_{\\rm cut}{=}4$. "
     "The neutral defect is short-ranged &mdash; no supercell / range-separation issue.",
     '<a href="pages/results.html#sec-3">Open results &rarr;</a>'),
    ("Near-edge scattering rates: 1st-Born $\\lvert M\\rvert^2$ vs full $T$ vs full-order $\\tilde V$", "Result", "2026-06-12", "prod", "Complete",
     "On-shell optical-theorem rates $\\hbar/\\tau_{nk}$ for all 677 states within 0.3 eV of the band edges "
     "(48&times;48 interpolation, $n_d{=}1\\%$), three treatments. <strong>The first-Born $\\lvert M\\rvert^2$ "
     "golden rule overestimates near-edge scattering by 1&ndash;1.5 orders of magnitude</strong> ($\\tau$ "
     "$\\approx$ 3 fs at the VBM); the full $T$-matrix in the same potential cuts it 29&ndash;45&times; "
     "($\\tau\\sim$100 fs) &mdash; the S-vacancy is a strong scatterer, so the Born series doesn't converge "
     "and resummation is essential. The full-order rest dressing then lowers it a further 2&ndash;3&times;. "
     "At the CBM the bare-$M$ $e$ resonance (0.16 eV below the edge) floods the rate; $\\tilde V$ moves it "
     "mid-gap and flattens it &mdash; the transport-level counterpart of Fig 15 panels (2) vs (4).",
     '<a href="pages/results.html#sec-5">Open results &rarr;</a>'),
    ("Spectral function at the supercell concentration ($n_d{=}1/36$)", "Result", "2026-06-14", "ok", "Complete",
     "S-vacancy $A(k,\\omega)$ at the density actually defined by the cell &mdash; one vacancy per "
     "$6\\times6$ = $n_d{=}1/36\\approx2.78\\%$ (Fig 24, full-order $\\tilde V$ vs bare $M$). Since "
     "$\\Sigma{=}n_d T$ is linear and $T(\\omega)$ is concentration-independent, $T$ is computed once and "
     "cached &mdash; after which <strong>any $n_d$ is a 2-second Dyson inversion</strong> (vs the ~7-min "
     "one-off $T$ cost). The in-gap $e$ resonance broadens linearly with $n_d$ at fixed position.",
     '<a href="pages/results.html#sec-5">Open results &rarr;</a>'),
    ("Substitutional defects O$_S$ / Se$_S$ (explicit-60 + spectral)", "Result", "2026-06-15", "warn", "Caveat",
     "The explicit-60 pipeline extended to two isovalent chalcogen substitutions with NO code change: the "
     "defect-species nonlocal $d_{\\mathrm{van}}$ is loaded by declaring O/Se as a <strong>zero-atom extra "
     "species in the host NSCF</strong> (&psi; unchanged). Full local+nonlocal $M$ from the relaxed cells, "
     "full-order rest dressing, spectral $A(k,\\omega)$ at $n_d{=}1/36$ (Figs 25&ndash;26). <strong>The two "
     "isovalent subs scatter oppositely</strong>: O$_S$ keeps the bands sharp (weak/localized, one mid-gap "
     "doublet $+0.73$); Se$_S$ broadens the whole valence edge across the BZ (strong/broadband; $\\sim$15 "
     "resonant near-VBM levels + a deep $+1.19$). A relaxation/alignment sensitivity study (unrelaxed + "
     "<code>pot_align=none</code>, Figs 27&ndash;28) shows the broadening &amp; dense manifolds are "
     "<strong>generated by the ionic relaxation</strong>, not the bare substitution &mdash; frozen-lattice "
     "gives sharp bands + isolated doublets ($9$&ndash;$13\\times$ taller in-gap peaks). "
     "<strong>CAVEAT (diagnosed):</strong> direct DFT of O$_S$ shows NO in-gap state; the substitution "
     "in-gap features (Figs 25&ndash;29) are a <strong>basis under-convergence artifact</strong> &mdash; "
     "the matrix element is correct, but the explicit-60 band basis is converged for the vacancy and NOT "
     "for the weak isovalent substitutions, whose apparent level descends steadily toward the VBM with "
     "band count (Fig 30) while the vacancy $e$ saturates at the DFT $+1.19$. Only the S-vacancy results "
     "are DFT-validated; the substitution in-gap levels are not physical.",
     '<a href="pages/results.html#sec-6">Open results &rarr;</a>'),
    ("Active-space dynamic resummation / transport", "Test", "&mdash;", "plan", "Planned",
     "Next (P6): the remaining mobility integral &mdash; velocity factors + BZ sum over the on-shell rates "
     "of &sect;5 (beyond-Born vs Born mobility), and the frequency-dependent $\\Sigma_{\\rm rest}(\\omega)$.", "&mdash;"),
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
def leak_scan():
    """placeholder-leak guard over EVERY generated page (a bare `$` in prose --
    e.g. a shell variable like $A -- can pair with a later one and swallow a
    protected span).  Cheap, and it catches the failure the per-builder checks
    do not see."""
    bad = []
    for root, _dirs, files in os.walk(os.path.join(ROOT, "docs")):
        for fn in files:
            if fn.endswith(".html"):
                fp = os.path.join(root, fn)
                with open(fp, encoding="utf-8") as f:
                    n = f.read().count("\x00")
                if n:
                    bad.append((os.path.relpath(fp, ROOT), n))
    if bad:
        for fp, n in bad:
            print("  [LEAK] %s: %d placeholder(s)" % (fp, n))
        raise SystemExit("placeholder leak -- fix the source markdown")
    print("  [leak scan] %s: clean" % "all pages")


def main():
    os.makedirs(PAGES_DIR, exist_ok=True)
    build_theory()
    build_plan()
    build_note()
    build_results()
    build_ladder()
    build_npol2cert()
    build_soczoom()
    build_moded()
    build_pdcoo2()
    build_pdvac()
    build_ptcoo2()
    build_deflated()
    build_tlres()
    build_methods()
    build_koster()
    build_feshplan()
    build_index()
    leak_scan()

    th = open(os.path.join(PAGES_DIR, "theory.html"), encoding="utf-8").read()
    pl = open(os.path.join(PAGES_DIR, "plan.html"),   encoding="utf-8").read()
    nt = open(os.path.join(PAGES_DIR, "note-kprime-normalization.html"), encoding="utf-8").read()
    rs = open(os.path.join(PAGES_DIR, "results.html"), encoding="utf-8").read()
    ld = open(os.path.join(PAGES_DIR, "sternheimer-ladder.html"), encoding="utf-8").read()
    ks = open(os.path.join(PAGES_DIR, "koster-slater.html"), encoding="utf-8").read()
    fp = open(os.path.join(PAGES_DIR, "feshbach-implementation.html"), encoding="utf-8").read()
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
    stats("ladder.html", ld)
    stats("koster.html", ks)
    stats("feshplan", fp)
    print("index.html  : %d bytes, %d catalog rows" % (len(ix), ix.count("<tr")))
    for nm, txt in (("theory.html", th), ("plan.html", pl), ("note-knorm", nt),
                    ("results.html", rs), ("ladder.html", ld), ("koster.html", ks),
                    ("feshplan", fp), ("index.html", ix)):
        p = check(txt)
        print("  [%s] %s" % (nm, "OK" if not p else " ; ".join(p)))

if __name__ == "__main__":
    main()
