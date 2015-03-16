.PHONY: clean-pyc clean-build docs clean

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "setup-venv - setup a development virutalenv in current directory."
	@echo "lint - check style with flake8"
	@echo "lint-readme - check README formatting for PyPI"
	@echo "test - run tests quickly with the default Python"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	@echo "release - package and upload a release"
	@echo "dist - package"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

setup-venv:
	if [ -f .venv ]; then virtualenv .venv; fi;
	. .venv/bin/activate && pip install -r requirements.txt && pip install -r dev-requirements.txt

lint:
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; flake8 --max-complexity 11 planemo tests

lint-readme:
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; python setup.py check -r -s

test:
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; nosetests tests planemo

coverage:
	coverage run --source planemo setup.py test
	coverage report -m
	coverage html
	open htmlcov/index.html || xdg-open htmlcov/index.htm

docs:
	rm -f docs/planemo.rst
	rm -f docs/planemo_ext.rst
	rm -f docs/modules.rst
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; sphinx-apidoc -f -o docs/ planemo_ext planemo_ext/galaxy/eggs
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; sphinx-apidoc -f -o docs/ planemo
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; python scripts/commands_to_rst.py
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; $(MAKE) -C docs clean
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; $(MAKE) -C docs html

open-docs: docs
	open docs/_build/html/index.html || xdg-open docs/_build/html/index.html

open-rtd: docs
	open https://planemo.readthedocs.org || xdg-open https://planemo.readthedocs.org

open-project:
	open https://github.com/galaxyproject/planemo || xdg-open https://github.com/galaxyproject/planemo

dist: clean
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; python setup.py sdist bdist_egg bdist_wheel
	ls -l dist

release-test: dist
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; twine upload -r test dist/*
	open https://testpypi.python.org/pypi/planemo || xdg-open https://testpypi.python.org/pypi/planemo

release:
	@while [ -z "$$CONTINUE" ]; do \
		read -r -p "Have you executed release-test and reviewed results? [y/N]: " CONTINUE; \
	done ; \
	[ $$CONTINUE = "y" ] || [ $$CONTINUE = "Y" ] || (echo "Exiting."; exit 1;)
	@echo "Releasing"
	if [ -f .venv/bin/activate ]; then . .venv/bin/activate; fi; twine upload dist/*

update-extern:
	sh scripts/update_extern.sh
