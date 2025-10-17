# Contributing to Panacea

Thank you for your interest in contributing to Panacea! 
This document explains how to set up a development environment, 
coding standards, testing, and our pull request process.

## Getting started

- Fork the repository and create a feature branch from `main`.
- Ensure you are using Python 3.9+.
- Install development dependencies:

```
python -m pip install -U pip
python -m pip install -e .[dev]
```

Note on pyhetdex: If you need `pyhetdex`, install it from the private index:
```
pip install --extra-index-url https://gate.mpe.mpg.de/pypi/simple/ pyhetdex
```

## Running tests

We use pytest. Once tests are added, run:
```
pytest -q
```
To see coverage (if configured):
```
pytest --cov=panacea --cov-report=term-missing
```

## Coding style and docstrings

- Follow PEP8
- Public functions/classes/methods should include Google‑style docstrings. Example:

```
def example(a, b):
    """Add two integers.

    Args:
        a: First number.
        b: Second number.

    Returns:
        The sum of a and b.

    Raises:
        ValueError: If inputs are invalid.
    """
    return a + b
```

## Commit messages

- Use descriptive language.
- Reference issues if applicable.

## Pull requests

1. Ensure tests pass locally: `pytest`.
2. Update documentation.
3. Open a PR against `main` with a clear description of the change and motivation.

## Release process (summary)

- Use semantic versioning. Tag releases as `vMAJOR.MINOR.PATCH`.
- After creating a release, update CITATION.cff and README badges.

## Code of Conduct

By participating, you agree to abide by our Code of Conduct (CODE_OF_CONDUCT.md).
