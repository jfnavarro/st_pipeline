# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at https://github.com/jfnavarro/st_pipeline/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

### Write Documentation

ST Pipeline could always use more documentation, whether as part of the
official ST Pipeline docs, in docstrings, or even on the web in blog posts,
articles, and such.

### Submit Feedback

The best way to send feedback is to file an issue at https://github.com/jfnavarro/st_pipeline/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.

## Get Started!

Ready to contribute? Here's how to set up `ST Pipeline` for local development.

1. Fork the `ST Pipeline` repo on GitHub.
2. Clone your fork locally

``` console
git clone git@github.com:jfnavarro/st_pipeline.git
```

3. Ensure [poetry](https://python-poetry.org/docs/) is installed.
4. Install dependencies and start your virtualenv:

``` console
poetry install -E test -E doc -E dev
```

Note that you can use your own Python environment (e.g Anaconda) by
changing the default behaviour in poetry with this command:

``` console
poetry config virtualenvs.create false
```

5. Create a branch for local development:

``` console
git checkout -b name-of-your-bugfix-or-feature
```

Now you can make your changes locally.

6. When you're done making changes, check that your changes pass the
   tests, including testing other Python versions, with pytest:

``` console
poetry run pytest
```

7. Commit your changes and push your branch to GitHub:

``` console
git add .
git commit -m "Your detailed description of your changes."
git push origin name-of-your-bugfix-or-feature
```

8. Submit a pull request through the GitHub website.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.md.
3. The pull request should work for Python 3.10, 3.11 and 3.12. Check
   https://github.com/jfnavarro/st_pipeline/actions
   and make sure that the tests pass for all supported Python versions.

## Testing

You can run the tests with pytest:

``` console
poetry run pytest
```

Replace test_your_module.py with the actual name of your test file.

## Makefile

A `makefile`is included in the repo with the following actions:

To run formatting tools

```bash
make format
```

To run linting tools

```bash
make lint
```

To run the tests

```bash
make unittet
```

To run the tests with coverage

```bash
make coverage
```

To clean the temporary files and cache

```bash
make clean
```

## Deploying

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in CHANGELOG.md).
Then run:

``` console
poetry run bump2version patch # possible: major / minor / patch
git push
git push --tags
```

GitHub Actions will then create a release and publish documentation if tests pass.

You can also create the documentation manually by running:

``` console
poetry run mkdocs build
```
