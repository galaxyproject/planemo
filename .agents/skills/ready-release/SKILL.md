---
name: ready-release
description: Walk a maintainer through cutting a planemo release (changelog, checks, make release, GitHub release).
---

Walk me through readying a planemo release. The full process is documented in `docs/developing.rst` — read that first and defer to it for details. This skill drives the interactive parts.

## Steps

1. **Confirm target version** - read `__version__` in `planemo/__init__.py`, show me the current `.devN` and the commits since the last tag, and ask me to confirm the target release version. `make release` mints the next patch, so a minor release needs `__version__` and the HISTORY `.devN` header edited first.

2. **Setup venv** - if `.venv` is missing, run `make setup-venv`.

3. **Populate changelog** - run `make add-history`, show me the new `HISTORY.rst` entries (flag any that duplicate hand-written entries), pause so I can edit, then commit the result (`make check-release` requires a clean tree).

4. **Pre-release checks** - run `make check-release`. If it fails, stop and report. Do not work around failures; fix the underlying issue.

5. **Release** - after I confirm, run `make release`. This pushes `master` and the tag to `UPSTREAM`, and the tag publishes to PyPI.

6. **Verify publish** - watch the tag's `Deploy` workflow run and confirm the version on PyPI.

7. **GitHub release** - after I confirm, create the release as in `docs/developing.rst`.

Stop after each step and report status before moving to the next.
