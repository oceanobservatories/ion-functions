# Contributing to ion-functions

Thank you for your interest in contributing to ion-functions. This document
covers how to report bugs, propose changes, and submit pull requests.

## Table of contents

- [Code of conduct](#code-of-conduct)
- [Getting started](#getting-started)
- [Reporting bugs](#reporting-bugs)
- [Proposing changes](#proposing-changes)
- [Development workflow](#development-workflow)
- [Code style](#code-style)
- [Docstrings](#docstrings)
- [Tests](#tests)
- [Documentation](#documentation)
- [Commit messages](#commit-messages)
- [Pull requests](#pull-requests)

---

## Code of conduct

This project follows the
[Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/)
code of conduct. By participating, you agree to uphold a welcoming and
respectful environment for everyone.

---

## Getting started

1. Fork the repository on GitHub.
2. Clone your fork locally:

   ```bash
   git clone https://github.com/<your-username>/ion-functions.git
   cd ion-functions
   ```

3. Create the development environment:

   ```bash
   conda env create -f conda_env.yml
   conda activate ion
   pip install -e .
   ```

4. Verify the setup:

   ```bash
   pytest
   ```

---

## Reporting bugs

Before filing a bug report, search the
[issue tracker](https://github.com/oceanobservatories/ion-functions/issues)
to see if the problem has already been reported.

When opening a new issue, include:

- A clear, descriptive title.
- The Python version and operating system you are using.
- A minimal, self-contained code example that reproduces the problem.
- The actual output or error you observed and what you expected instead.

---

## Proposing changes

For small, focused fixes (typos, documentation corrections, minor bug fixes)
you can open a pull request directly. For larger changes — new instrument
families, algorithm updates, or structural refactors — open an issue first to
discuss the approach before investing time in implementation.

---

## Development workflow

Work on a feature branch off `master`:

```bash
git checkout master
git pull upstream master
git checkout -b my-feature-branch
```

Keep your branch focused on a single logical change. Open separate pull
requests for unrelated fixes rather than bundling them together.

---

## Code style

- **Python version:** 3.13+.
- **Line length:** wrap code and comments at 120 characters.
- **Formatting:** follow [PEP 8](https://peps.python.org/pep-0008/). Use
  `black` or `ruff` if you prefer an automated formatter, but do not reformat
  code outside the scope of your change.
- **Types:** type annotations are welcome but not required.
- **Imports:** standard library first, then third-party, then local — each
  group separated by a blank line.
- **No magic numbers:** use named constants or document the source of a
  numerical value in a comment.

### Function conventions

- All functions in `ion_functions/data/` accept and return NumPy arrays.
- Functions must be stateless (no class-level or module-level mutable state).
- Missing or fill values use `NaN` (float) or `-9999999` (integer).
- Do not alter the signature or logic of existing functions without discussion;
  OOI data processing pipelines depend on stable interfaces.

---

## Docstrings

All public functions use **NumPy-style docstrings**. Wrap docstrings at
79 characters.

Minimum required sections for a data-processing function:

```python
def my_function(arg1, arg2):
    """
    One-line summary of what this function computes.

    Parameters
    ----------
    arg1 : ndarray
        Description, including OOI variable name and units [units].
    arg2 : float
        Description.

    Returns
    -------
    result : ndarray
        Description, including OOI data product code and level [units].

    Notes
    -----
    Brief notes (2-5 lines). Extended content belongs in the docs site.
    """
```

Use plain ASCII in docstrings — no Unicode, Greek letters, or special
symbols. Write units as plain text (e.g., `degC`, `m s^-1`, `umol/kg`).

Do not add a `References` section to individual docstrings; citations belong
in the corresponding documentation page under `docs/api/`.

---

## Tests

Every function in `ion_functions/data/` should have a corresponding test in
`ion_functions/data/test/`. Tests use `unittest` and follow the existing
conventions in that directory.

- Test against known-good values derived from the OOI Data Product
  Specification (DPS) worked examples or from independent reference
  implementations where available.
- Include at least one test for a representative input array and one for
  scalar or edge-case inputs.
- Do not introduce test dependencies beyond those already in `conda_env.yml`.

Run the full test suite before opening a pull request:

```bash
pytest
```

---

## Documentation

The documentation site is built with
[MkDocs](https://www.mkdocs.org/) and
[mkdocstrings](https://mkdocstrings.github.io/). Sources live under `docs/`.

To build and preview locally:

```bash
mkdocs serve
```

Then open `http://127.0.0.1:8000` in your browser. The site rebuilds
automatically as you edit files.

If you add a new instrument family module, add a corresponding page under
`docs/api/` and add it to the `nav` section of `mkdocs.yml`.

---

## Commit messages

Use short, imperative-mood subject lines (50 characters or fewer):

```
Fix sign error in co2_functions t_coeff correction
Add NumPy docstrings for nit_functions
Update conda_env.yml to require gsw >= 3.6
```

For non-trivial changes, add a body paragraph separated from the subject
by a blank line that explains *why* the change was made, not just what.

---

## Pull requests

Before submitting:

- [ ] All existing tests pass (`pytest`).
- [ ] New or modified functions have NumPy docstrings.
- [ ] New functions have tests.
- [ ] If you changed the public API, the relevant `docs/api/` page is updated.
- [ ] Your branch is up to date with `master`.

In the pull request description:

- Summarize what changed and why.
- Reference any related issues (e.g., `Closes #123`).
- Note any algorithm decisions that reviewers should be aware of, especially
  if the implementation differs from the DPS or a cited source.

A maintainer will review your pull request and may request changes before
merging. Please respond to review comments promptly and keep the discussion
focused on the code.
