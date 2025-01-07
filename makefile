sources = geomxaggregation

.PHONY: test format lint unittest coverage pre-commit clean
test: format lint unittest

format:
	poetry run isort $(sources) tests
	poetry run  ruff format $(sources) tests

lint:
	poetry run ruff check --fix $(sources) tests
	poetry run  mypy $(sources)

unittest:
	poetry run  pytest

coverage:
	poetry run  pytest --cov=$(sources) --cov-branch --cov-report=term-missing tests

pre-commit:
	pre-commit run --all-files

clean:
	rm -rf .mypy_cache .pytest_cache
	rm -rf *.egg-info
	rm -rf .tox dist site
	rm -rf coverage.xml .coverage
