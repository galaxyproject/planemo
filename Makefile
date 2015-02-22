.PHONY: clean-pyc clean-build docs clean

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests quickly with the default Python"
	@echo "test-all - run tests on every Python version with tox"
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

lint:
	flake8 --max-complexity 11 planemo tests

lint-readme:
	python setup.py check -r -s

test:
	python setup.py test

test-all:
	tox

coverage:
	coverage run --source planemo setup.py test
	coverage report -m
	coverage html
	open htmlcov/index.html

docs:
	rm -f docs/planemo.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ planemo
	sphinx-apidoc -o docs/ planemo_ext
	python scripts/commands_to_rst.py
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	open docs/_build/html/index.html

dist: clean
	python setup.py sdist bdist_egg bdist_wheel
	ls -l dist

release-test: dist
	twine upload -r test dist/*
	echo "Review https://testpypi.python.org/pypi/planemo"

release:
	@while [ -z "$$CONTINUE" ]; do \
		read -r -p "Have you executed release-test and reviewed results? [y/N]: " CONTINUE; \
	done ; \
	[ $$CONTINUE = "y" ] || [ $$CONTINUE = "Y" ] || (echo "Exiting."; exit 1;)
	@echo "Releasing"
	twine upload dist/*
