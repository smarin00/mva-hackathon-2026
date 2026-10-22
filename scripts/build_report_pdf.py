#!/usr/bin/env python3
"""
Render the Track 2 report from Markdown to a styled PDF (A4).

  python3 scripts/build_report_pdf.py docs/smarin00_track2_report.md \\
      docs/smarin00_track2_report.pdf

Needs python-markdown and Playwright with Chromium (not part of the analysis
environment: `pip install markdown playwright && playwright install chromium`).
Figures referenced as SVG are inlined so their text is set in the page fonts;
PNG figures are embedded as data URIs. A figure is an image paragraph followed
by a paragraph starting with "**Figure N.**", which becomes its caption.
"""
import base64
import html
import os
import re
import sys

import markdown

CSS = r"""
@page { size: A4; margin: 17mm 17mm 16mm 17mm; }
:root {
  --ink: #0b0b0b; --ink2: #52514e; --muted: #898781; --rule: #e1e0d9;
  --wash: #f6f5f1; --accent: #2a78d6;
}
html { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
body { font-family: "Source Serif 4", "Source Serif Pro", Lora, Georgia, serif; font-size: 9.6pt;
       line-height: 1.43; color: var(--ink); margin: 0; background: #fff;
       hyphens: auto; -webkit-hyphens: auto; text-rendering: optimizeLegibility; }
h1, h2, h3, h4, th, figcaption, .meta, table, blockquote, code, svg { font-family: "Inter", "Inter Sys", sans-serif; }
h1 { font-size: 19pt; line-height: 1.18; font-weight: 700; letter-spacing: -0.012em;
     margin: 0 0 7pt 0; }
.meta { font-size: 8.2pt; color: var(--ink2); line-height: 1.5; margin: 0 0 12pt 0;
        padding-bottom: 9pt; border-bottom: 1.2pt solid var(--ink); }
.meta strong { color: var(--ink); font-weight: 600; }
.meta a { color: var(--accent); text-decoration: none; }
h2 { font-size: 12.2pt; font-weight: 700; margin: 17pt 0 5pt 0; padding-top: 5pt;
     border-top: 0.6pt solid var(--rule); break-after: avoid; page-break-after: avoid; }
h3 { font-size: 10.2pt; font-weight: 650; margin: 12pt 0 3pt 0; break-after: avoid;
     page-break-after: avoid; }
h2 + p, h3 + p { margin-top: 0; }
p { margin: 0 0 6.5pt 0; text-align: justify; orphans: 3; widows: 3; }
a { color: var(--accent); text-decoration: none; }
strong { font-weight: 700; }
code { font-size: 7.6pt; background: var(--wash); padding: 0.3pt 1.5pt; border-radius: 2pt;
       color: #2e2d2a; overflow-wrap: anywhere; }
ul, ol { margin: 0 0 7pt 0; padding-left: 15pt; }
li { margin: 0 0 2.5pt 0; }
blockquote { margin: 0 0 12pt 0; padding: 10pt 13pt 4pt 13pt; background: var(--wash);
             border-left: 2.4pt solid var(--accent); border-radius: 0 3pt 3pt 0;
             font-size: 8.7pt; line-height: 1.5; color: #1d1c1a; }
blockquote p { text-align: left; margin-bottom: 6pt; }
blockquote strong { font-weight: 650; color: var(--ink); }
table { border-collapse: collapse; width: 100%; margin: 4pt 0 10pt 0; font-size: 7.7pt;
        line-height: 1.38; break-inside: avoid; page-break-inside: avoid; }
thead th { text-align: left; font-weight: 600; color: var(--ink2); font-size: 7.2pt;
           text-transform: none; border-bottom: 0.9pt solid var(--ink); padding: 3pt 6pt 3pt 0; }
tbody td { border-bottom: 0.5pt solid var(--rule); padding: 3.2pt 6pt 3.2pt 0;
           vertical-align: top; }
tbody tr:last-child td { border-bottom: 0.9pt solid var(--ink); }
td strong { font-weight: 650; }
table.kv thead { display: none; }
table.kv td:first-child { width: 24%; font-weight: 600; color: var(--ink); }
table.kv tbody tr:first-child td { border-top: 0.9pt solid var(--ink); }
figure { margin: 8pt 0 12pt 0; break-inside: avoid; page-break-inside: avoid; }
figure svg, figure img { width: 100%; height: auto; display: block; }
figcaption { font-size: 7.6pt; line-height: 1.45; color: var(--ink2); margin-top: 5pt;
             text-align: left; }
figcaption strong { color: var(--ink); font-weight: 650; }
figcaption code { font-size: 7pt; }
ol.refs { font-size: 7.4pt; line-height: 1.38; padding-left: 16pt; color: #2e2d2a;
           columns: 2; column-gap: 16pt; }
ol.refs li { break-inside: avoid; }
section.keep { break-inside: avoid; page-break-inside: avoid; }
ol.refs li { margin-bottom: 2pt; text-align: left; }
ol.refs em { font-style: italic; }
.pb { break-before: page; page-break-before: always; }
"""


RANGES = {
    "latin": "U+0000-00FF, U+0131, U+0152-0153, U+02BB-02BC, U+02C6, U+02DA, U+02DC, "
             "U+0304, U+0308, U+0329, U+2000-206F, U+2074, U+20AC, U+2122, U+2191, U+2193, "
             "U+2212, U+2215, U+FEFF, U+FFFD",
    "greek": "U+0370-03FF",
}


def font_faces():
    """Optional: embed Source Serif 4 and Inter (woff2 from fontsource:
    `npm pack @fontsource/source-serif-4 @fontsource/inter`) from the folder in
    REPORT_FONTS; glyphs outside those subsets fall back to the system fonts."""
    d = os.environ.get("REPORT_FONTS")
    css = ['@font-face { font-family: "Inter Sys"; src: local("Inter"); }']
    if not d or not os.path.isdir(d):
        return "\n".join(css)
    for fam, prefix in (("Source Serif 4", "source-serif-4"), ("Inter", "inter")):
        for sub, rng in RANGES.items():
            for w in (400, 500, 600, 700):
                for st in ("normal", "italic"):
                    f = os.path.join(d, f"{prefix}-{sub}-{w}-{st}.woff2")
                    if not os.path.exists(f):
                        continue
                    b64 = base64.b64encode(open(f, "rb").read()).decode()
                    css.append(f'@font-face {{ font-family: "{fam}"; font-weight: {w}; '
                               f'font-style: {st}; unicode-range: {rng}; '
                               f'src: url(data:font/woff2;base64,{b64}) format("woff2"); }}')
    return "\n".join(css)


def inline_svg(path):
    svg = open(path).read()
    svg = re.sub(r"<\?xml[^>]*\?>", "", svg)
    svg = re.sub(r"<!DOCTYPE[^>]*>", "", svg)
    svg = re.sub(r'(<svg[^>]*?)\swidth="[^"]*"', r"\1", svg, count=1)
    svg = re.sub(r'(<svg[^>]*?)\sheight="[^"]*"', r"\1", svg, count=1)
    return svg.replace("'Inter'", "'Inter', 'Inter Sys', sans-serif")


def data_uri(path):
    ext = os.path.splitext(path)[1].lstrip(".").lower()
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg"}[ext]
    return f"data:{mime};base64," + base64.b64encode(open(path, "rb").read()).decode()


def sectionize(body):
    """Wrap each h2 section and h3 subsection in <section> so that short ones are
    kept on one page (break-inside: avoid; longer ones still break)."""
    head, *parts = re.split(r"(?=<h2>)", body)
    out = [head]
    for part in parts:
        sub_head, *subs = re.split(r"(?=<h3>)", part)
        inner = sub_head + "".join(f'<section class="s3">{x}</section>' for x in subs)
        out.append(f'<section class="s2{"" if subs else " leaf"}">{inner}</section>')
    return "".join(out)


def build_html(md_path):
    base = os.path.dirname(os.path.abspath(md_path))
    text = open(md_path).read()
    body = markdown.markdown(text, extensions=["tables", "sane_lists"])

    # title block: h1 + the paragraph after it
    body = re.sub(r"(<h1>.*?</h1>)\s*<p>(.*?)</p>",
                  lambda m: m.group(1) + f'<div class="meta">{m.group(2)}</div>', body,
                  count=1, flags=re.S)
    body = re.sub(r'(<div class="meta">)(.*?)(</div>)',
                  lambda m: m.group(1) + m.group(2).replace("\n", "<br>") + m.group(3), body,
                  count=1, flags=re.S)
    body = re.sub(r"(https://github\.com/[\w./-]+)(?![^<]*</a>)",
                  r'<a href="\1">\1</a>', body)

    # figures: image paragraph + caption paragraph
    def fig(m):
        src, cap = m.group(1), m.group(2)
        path = os.path.normpath(os.path.join(base, html.unescape(src)))
        if path.endswith(".svg"):
            media = inline_svg(path)
        else:
            media = f'<img src="{data_uri(path)}" alt="">'
        return f"<figure>{media}<figcaption>{cap}</figcaption></figure>"
    body = re.sub(r'<p><img alt="[^"]*" src="([^"]+)" ?/?></p>\s*<p>(<strong>Figure \d+\..*?)</p>',
                  fig, body, flags=re.S)

    # line-break opportunities inside code paths, so justified lines do not stretch
    body = re.sub(r"<code>(.*?)</code>",
                  lambda m: "<code>" + re.sub(r"([/_.])", r"\1<wbr>", m.group(1)) + "</code>",
                  body)

    # two-column key/value tables (empty header row)
    body = re.sub(r"<table>\s*<thead>\s*<tr>\s*<th></th>\s*<th></th>\s*</tr>\s*</thead>",
                  '<table class="kv"><thead><tr><th></th><th></th></tr></thead>', body)

    # references block
    body = re.sub(r"(<h2>References</h2>\s*)<ol>", r'\1<ol class="refs">', body)
    body = sectionize(body)
    title = re.search(r"<h1>(.*?)</h1>", body, re.S).group(1)
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            f"<title>{re.sub('<[^>]+>', '', title)}</title><style>{font_faces()}{CSS}</style></head>"
            f"<body>{body}</body></html>")


def main():
    md_path, pdf_path = sys.argv[1], sys.argv[2]
    page = build_html(md_path)
    import tempfile
    html_path = os.path.join(tempfile.mkdtemp(), "report.html")
    open(html_path, "w").write(page)
    from playwright.sync_api import sync_playwright
    foot = ("<div style='font-family:Inter,sans-serif;font-size:6.5pt;color:#898781;"
            "width:100%;padding:0 17mm;display:flex;justify-content:space-between'>"
            "<span>BUB1B-related MVA: from genotype to drug hypotheses · Track 2 · smarin00</span>"
            "<span><span class='pageNumber'></span> / <span class='totalPages'></span></span></div>")
    with sync_playwright() as p:
        b = p.chromium.launch()
        pg = b.new_page()
        pg.goto("file://" + os.path.abspath(html_path))
        pg.wait_for_load_state("networkidle")
        # keep a (sub)section on one page only if it is short: < 40% of the page
        # body height (264 mm); longer ones flow, so no page is left half empty
        pg.emulate_media(media="print")
        pg.evaluate("""() => {
            const page = 264 / 25.4 * 96;
            for (const s of document.querySelectorAll('section.s3, section.leaf'))
                if (s.getBoundingClientRect().height < 0.40 * page) s.classList.add('keep');
        }""")
        pg.pdf(path=pdf_path, format="A4", print_background=True,
               display_header_footer=True, header_template="<span></span>",
               footer_template=foot,
               margin={"top": "17mm", "bottom": "16mm", "left": "17mm", "right": "17mm"})
        b.close()
    print(f"wrote {pdf_path}")


if __name__ == "__main__":
    main()
