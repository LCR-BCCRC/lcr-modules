"""
Post-processing for igv-reports --tabulator HTML output.

igv.createBrowser() can reject its promise when it fails to load the
session for the first variant (e.g. bad reads, unusual CIGAR strings, or
annotation track issues at that specific locus). The default template has
no .catch() handler, so a rejection silently prevents initTable() from
ever running — the pileup partially renders but the variant table stays
empty.

This script patches the generated HTML to:
1. Add a .catch() on igv.createBrowser() so initTable() runs even on failure
2. Guard igvBrowser in the row-click handler so clicking rows is safe when
   igv.js failed to initialize
"""

import sys

html_file = sys.argv[1]

with open(html_file) as f:
    content = f.read()

# Patch 1: add .catch() after the .then() block so initTable() always runs
old_then = """.then(function (b) {
                igvBrowser = b
                initTable()
            })"""

new_then = """.then(function (b) {
                igvBrowser = b
                initTable()
            })
            .catch(function(err) {
                console.error("igv-reports: igv.createBrowser failed:", err)
                initTable()
            })"""

# Patch 2: guard igvBrowser in row-click handler so it doesn't throw when
# igv.js failed and igvBrowser was never set
old_click = """igvBrowser.loadSession({
                    url: session
                })"""

new_click = """if (igvBrowser) {
                    igvBrowser.loadSession({url: session})
                }"""

patched = content.replace(old_then, new_then).replace(old_click, new_click)

if patched == content:
    print(f"WARNING: patch_igv_reports_html.py found no matching patterns in {html_file}", file=sys.stderr)

with open(html_file, "w") as f:
    f.write(patched)
