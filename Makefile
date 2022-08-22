# Default tests run with make test and make quick-tests
TESTS?=tests planemo
# Default environment for make tox
ENV?=py27
# Extra arguments supplied to tox command
ARGS?=
# Location of virtualenv used for development.
VENV?=.venv
VENV3?=.venv3
# Open resource on Mac OS X or Linux
OPEN_RESOURCE=bash -c 'open $$0 || xdg-open $$0'
# Source virtualenv to execute command (flake8, sphinx, twine, etc...)
IN_VENV=if [ -f $(VENV)/bin/activate ]; then . $(VENV)/bin/activate; fi;
IN_VENV3=if [ -f $(VENV3)/bin/activate ]; then . $(VENV3)/bin/activate; fi;
# TODO: add this upstream as a remote if it doesn't already exist.
UPSTREAM?=galaxyproject
SOURCE_DIR?=planemo
SOURCE_DOC_EXCLUDE=$(SOURCE_DIR)/cwl/cwl2script
BUILD_SCRIPTS_DIR=scripts
VERSION?=$(shell python $(BUILD_SCRIPTS_DIR)/print_version_for_release.py $(SOURCE_DIR))
DOC_URL?=https://planemo.readthedocs.org
PROJECT_URL?=https://github.com/galaxyproject/planemo
PROJECT_NAME?=planemo
TEST_DIR?=tests
DOCS_DIR?=docs
BUILD_SLIDESHOW?=$(IN_VENV) python $(BUILD_SCRIPTS_DIR)/build_slideshow.py
SLIDESHOW_TO_PDF?=bash -c 'docker run --rm -v `pwd`:/cwd astefanutti/decktape /cwd/$$0 /cwd/`dirname $$0`/`basename -s .html $$0`.pdf'

.PHONY: clean-pyc clean-build docs clean

help:
	@egrep '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

submodule:
	git submodule init
	git submodule update

install: submodule ## install into Python envirnoment
	python setup.py install && cd cwl-runner && python setup.py install

setup-venv: ## setup a development virtualenv in current directory
	if [ ! -d $(VENV) ]; then virtualenv $(VENV); exit; fi;
	$(IN_VENV) pip install --upgrade pip && pip install -r dev-requirements.txt -r requirements.txt

setup-venv3: ## setup a development virtualenv in current directory
	if [ ! -d $(VENV3) ]; then virtualenv -p python3 $(VENV3); exit; fi;
	$(IN_VENV3) pip install --upgrade pip && pip install -r dev-requirements.txt -r requirements.txt

setup-git-hook-lint: ## setup precommit hook for linting project
	cp $(BUILD_SCRIPTS_DIR)/pre-commit-lint .git/hooks/pre-commit

setup-git-hook-lint-and-test: ## setup precommit hook for linting and testing project
	cp $(BUILD_SCRIPTS_DIR)/pre-commit-lint-and-test .git/hooks/pre-commit

flake8: ## check style using flake8 for current Python (faster than lint)
	$(IN_VENV) flake8 $(SOURCE_DIR) $(TEST_DIR)

lint: ## check style using tox and flake8 for Python 2 and Python 3
	$(IN_VENV) tox -e py37-lint

test: ## run tests with the default Python (faster than tox)
	$(IN_VENV) pytest $(TESTS)

quick-test: ## run quickest tests with the default Python
	$(IN_VENV) PLANEMO_SKIP_SLOW_TESTS=1 PLANEMO_SKIP_GALAXY_TESTS=1 pytest $(TESTS)

gxwf-test-test:  ## run test of workflow testing script
	bash $(BUILD_SCRIPTS_DIR)/test_workflow_tests.sh

tox: ## run tests with tox in the specified ENV, defaults to py27
	$(IN_VENV) tox -e $(ENV) -- $(ARGS)

_coverage-report: ## build coverage report with the default Python
	coverage run --source $(SOURCE_DIR) setup.py $(TEST_DIR)
	coverage report -m
	coverage html

_open-coverage: ## open coverage report using open
	open htmlcov/index.html || xdg-open htmlcov/index.html

coverage: _coverage-report open-coverage ## check code coverage quickly with the default Python

open-history:  # view HISTORY.rst as HTML.
	rst2html HISTORY.rst > /tmp/planemo_history.html
	$(OPEN_RESOURCE) /tmp/planemo_history.html

ready-docs:  ## rebuild docs folder ahead of running docs or lint-docs
	rm -f $(DOCS_DIR)/$(SOURCE_DIR).rst
	rm -f $(DOCS_DIR)/modules.rst
	$(BUILD_SLIDESHOW) 'Galaxy Tool Framework Changes' $(DOCS_DIR)/galaxy_changelog.md
	$(BUILD_SLIDESHOW) 'Planemo: A Scientific Workflow SDK' $(DOCS_DIR)/presentations/2016_workflows.md
	$(IN_VENV) sphinx-apidoc -f -o $(DOCS_DIR)/ $(SOURCE_DIR) $(SOURCE_DOC_EXCLUDE)
	$(IN_VENV) python scripts/commands_to_rst.py

docs: ready-docs ## generate Sphinx HTML documentation, including API docs
	$(IN_VENV) $(MAKE) -C $(DOCS_DIR) clean
	$(IN_VENV) $(MAKE) -C $(DOCS_DIR) html

lint-docs: ready-docs
	$(IN_VENV) $(MAKE) -C $(DOCS_DIR) clean
	$(IN_VENV) $(MAKE) -C $(DOCS_DIR) html 2>&1 | python $(BUILD_SCRIPTS_DIR)/lint_sphinx_output.py

_open-docs:
	$(OPEN_RESOURCE) $(DOCS_DIR)/_build/html/index.html

open-slides-galaxy-changelog: ready-docs
	$(OPEN_RESOURCE) $(DOCS_DIR)/galaxy_changelog.html

open-slides-workflows: ready-docs
	$(OPEN_RESOURCE) $(DOCS_DIR)/presentations/2016_workflows.html

export-slides-workflows: ready-docs
	$(SLIDESHOW_TO_PDF) $(DOCS_DIR)/presentations/2016_workflows.html

open-docs: docs _open-docs ## generate Sphinx HTML documentation and open in browser

open-rtd: docs ## open docs on readthedocs.org
	$(OPEN_RESOURCE) $(PROJECT_URL)

open-project: ## open project on github
	$(OPEN_RESOURCE) $(PROJECT_URL)

dist: clean submodule ## create and check packages
	$(IN_VENV) python setup.py sdist bdist_wheel
	$(IN_VENV) twine check dist/*
	ls -l dist

release-test-artifacts: dist
	$(IN_VENV) twine upload -r test dist/*
	$(OPEN_RESOURCE) https://testpypi.python.org/pypi/$(PROJECT_NAME)

release-artifacts: release-test-artifacts ## Package and Upload to PyPi
	@while [ -z "$$CONTINUE" ]; do \
		read -r -p "Have you executed release-test and reviewed results? [y/N]: " CONTINUE; \
	done ; \
	[ $$CONTINUE = "y" ] || [ $$CONTINUE = "Y" ] || (echo "Exiting."; exit 1;)
	@echo "Releasing"
	$(IN_VENV) twine upload dist/*

commit-version: ## Update version and history, commit and add tag
	$(IN_VENV) python $(BUILD_SCRIPTS_DIR)/commit_version.py $(SOURCE_DIR) $(VERSION)

new-version: ## Mint a new version
	$(IN_VENV) python $(BUILD_SCRIPTS_DIR)/new_version.py $(SOURCE_DIR) $(VERSION)

release-local: commit-version new-version

release-brew: ## Mint a new homebrew release
	bash $(BUILD_SCRIPTS_DIR)/update_planemo_recipe.bash $(VERSION)

push-release: ## Push a tagged release to github
	git push $(UPSTREAM) master
	git push --tags $(UPSTREAM)

release: release-local push-release ## package, review, and upload a release

add-history: ## Reformat HISTORY.rst with data from Github's API
	$(IN_VENV) python $(BUILD_SCRIPTS_DIR)/bootstrap_history.py $(ITEM)

update-extern: ## update external artifacts copied locally
	sh scripts/update_extern.sh
