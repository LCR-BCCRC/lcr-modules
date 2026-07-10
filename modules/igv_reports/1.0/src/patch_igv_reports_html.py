"""
Post-processing for igv-reports --tabulator HTML output.

igv.createBrowser() can reject its promise when it fails to load the
session for the first variant (e.g. bad reads, unusual CIGAR strings, or
annotation track issues at that locus). The default template bundles
browser creation and session loading into one promise, so a rejection
silently prevents initTable() from ever running.

This script patches the generated HTML to:
1. Decouple browser creation from initial session loading so igvBrowser is
   always set and row-click navigation always works.
2. Load the first variant's session separately after the browser is ready,
   with a .catch() so a bad first-variant session doesn't break the table.
"""

import sys

html_file = sys.argv[1]

with open(html_file) as f:
    content = f.read()

# Remove sessionURL from createBrowser options so browser creation is
# decoupled from (potentially-failing) session loading.
old_session_url_line = "            sessionURL: sessionDictionary[\"0\"],\n"
new_session_url_line = ""

# After igvBrowser is assigned and initTable() is called, load the first
# session separately with its own error handler.
old_then_body = """                igvBrowser = b
                initTable()
            })"""

new_then_body = """                igvBrowser = b
                initTable()
                b.loadSession({url: sessionDictionary["0"]}).catch(function(err) {
                    console.error("igv-reports: failed to load initial session:", err)
                })
            })"""

patched = (content
           .replace(old_session_url_line, new_session_url_line)
           .replace(old_then_body, new_then_body))

if patched == content:
    print(f"WARNING: patch_igv_reports_html.py found no matching patterns in {html_file}", file=sys.stderr)

with open(html_file, "w") as f:
    f.write(patched)
