Walk me through readying a planemo release. Run each step, check for problems, and pause if anything looks wrong.

## Steps

1. **Check git status** - ensure working tree is clean, no missing files.

2. **Verify version** - read `planemo/__init__.py` and confirm `__version__` is a `.devN` variant of the intended release. Show me the current version and ask me to confirm the target release version.

3. **Setup venv** - check `.venv` exists. If not, run `make setup-venv`. Confirm dev-requirements are installed.

4. **Check UPSTREAM remote** - the Makefile defaults `UPSTREAM` to `galaxyproject`. Check if `$UPSTREAM` is set in the environment; if not, check if a git remote named `galaxyproject` exists. If it doesn't, check for `origin` or `upstream` remotes pointing to `galaxyproject/planemo` and offer to create a `galaxyproject` alias via `git remote add galaxyproject <url>`. This must be resolved before `make release` can push.

5. **Add history** - run `make add-history` to pull contributions into HISTORY.rst under the .dev0 entry. Show me the new HISTORY.rst additions for review.

6. **Lint** - run `make clean && make lint`. Report any failures.

7. **Review** - show me a summary of outstanding uncommitted changes (if any) and ask if I want to commit them before proceeding.

8. **Release** - after I confirm, run `make release` which does: commit-version, new-version, check-dist, push-release. This tags, bumps to next dev version, and pushes upstream.

Stop after each step and report status before moving to the next.
