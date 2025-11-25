from pathlib import Path
import re

REPO_ROOT = Path(__file__).resolve().parents[1]
DOC_PATHS = [REPO_ROOT / "README.md", *sorted((REPO_ROOT / "src" / "spectrumpy_flight" / "docs").glob("*.md"))]

LATEX_PATTERN = re.compile(r"(\$\$|\\\\\(|\\\\\)|\\\\\[|\\\\\]|\\begin\{)")


def test_docs_do_not_contain_raw_latex():
    offending = []
    for path in DOC_PATHS:
        text = path.read_text(encoding="utf-8")
        match = LATEX_PATTERN.search(text)
        if match:
            offending.append(f"{path.relative_to(REPO_ROOT)} -> {match.group(0)}")
    assert not offending, "Remove raw LaTeX delimiters from documentation: " + ", ".join(offending)
