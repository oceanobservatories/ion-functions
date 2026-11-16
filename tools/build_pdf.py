"""
Render the MkDocs print-site page to a single PDF.

Run after ``mkdocs build`` (the print-site plugin writes
``site/print_page/index.html``). Requires the Playwright Python bindings
(conda-forge ``playwright-python``) with Chromium installed
(``playwright install chromium``). MathJax is loaded from
cdnjs, so the machine needs network access while the page renders.

Usage::

    mkdocs build --strict
    python tools/build_pdf.py
"""
import argparse
from pathlib import Path

from playwright.sync_api import sync_playwright


def build_pdf(site: Path, out: Path, scale: float = 0.85) -> None:
    """
    Render ``<site>/print_page/index.html`` to a PDF.

    Parameters
    ----------
    site : Path
        MkDocs build directory.
    out : Path
        Output PDF path.
    scale : float, optional
        Rendering scale passed to Chromium (0.1 to 2). The default of
        0.85 brings the theme's 16 px body text to roughly 10 pt.
    """
    page_path = (site / "print_page" / "index.html").resolve()
    if not page_path.exists():
        raise SystemExit(f"{page_path} not found; run mkdocs build with the print-site plugin first")

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.goto(page_path.as_uri(), wait_until="networkidle")
        # wait for MathJax 3 to finish typesetting every equation
        page.wait_for_function("window.MathJax && MathJax.startup && MathJax.startup.promise")
        page.evaluate("MathJax.startup.promise")
        page.pdf(path=str(out), format="Letter", scale=scale, print_background=True, outline=True, tagged=True,
                 margin={"top": "0.6in", "bottom": "0.6in", "left": "0.6in", "right": "0.6in"})
        browser.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Render the MkDocs print-site page to a single PDF.")
    parser.add_argument("--site", type=Path, default=Path("site"), help="MkDocs build directory (default: site)")
    parser.add_argument("--out", type=Path, default=Path("site/ion-functions.pdf"),
                        help="output PDF path (default: site/ion-functions.pdf)")
    parser.add_argument("--scale", type=float, default=0.85,
                        help="rendering scale, 0.1 to 2 (default: 0.85)")
    args = parser.parse_args()
    if not 0.1 <= args.scale <= 2:
        parser.error("--scale must be between 0.1 and 2")
    build_pdf(args.site, args.out, args.scale)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
