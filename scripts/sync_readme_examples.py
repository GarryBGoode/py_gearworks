# Copyright 2026 Gergely Bencsik
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Keep README.md code blocks in sync with the example files they came from.

README.md marks each managed code block with a pair of HTML comments naming
the source file:

    <!-- BEGIN EXAMPLE: examples/readme_example.py -->
    ```python
    ...
    ```
    <!-- END EXAMPLE -->

The referenced file marks the portion to embed with matching markers:

    # --8<-- [start:snippet]
    ...
    # --8<-- [end:snippet]

Run with no arguments to rewrite README.md in place. Run with --check to
verify README.md is already up to date (used in CI); the process exits with
a non-zero status if it is not.
"""

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
README_PATH = REPO_ROOT / "README.md"

BEGIN_RE = re.compile(r"<!-- BEGIN EXAMPLE: (?P<path>.+?) -->")
END_MARKER = "<!-- END EXAMPLE -->"
SNIPPET_START = "--8<-- [start:snippet]"
SNIPPET_END = "--8<-- [end:snippet]"


def extract_snippet(source_path: Path) -> str:
    lines = source_path.read_text(encoding="utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if SNIPPET_START in line)
    end = next(i for i, line in enumerate(lines) if SNIPPET_END in line)
    return "\n".join(lines[start + 1 : end])


def sync(readme_text: str) -> str:
    out = []
    i = 0
    lines = readme_text.splitlines(keepends=True)
    while i < len(lines):
        line = lines[i]
        match = BEGIN_RE.search(line)
        if not match:
            out.append(line)
            i += 1
            continue

        source_path = REPO_ROOT / match.group("path")
        snippet = extract_snippet(source_path)
        out.append(line)
        out.append("```python\n")
        out.append(snippet + "\n")
        out.append("```\n")

        # skip over the old fenced block and any content up to END marker
        i += 1
        while i < len(lines) and END_MARKER not in lines[i]:
            i += 1
        out.append(lines[i])  # the END marker line itself
        i += 1

    return "".join(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify README.md is up to date instead of rewriting it",
    )
    args = parser.parse_args()

    original = README_PATH.read_text(encoding="utf-8")
    updated = sync(original)

    if args.check:
        if original != updated:
            print(
                "README.md is out of date with its example files. "
                "Run `python scripts/sync_readme_examples.py` to fix.",
                file=sys.stderr,
            )
            return 1
        print("README.md is up to date.")
        return 0

    if original != updated:
        README_PATH.write_text(updated, encoding="utf-8")
        print("README.md updated.")
    else:
        print("README.md already up to date.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
