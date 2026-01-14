# Convenience targets to mirror CI locally

.PHONY: help lint test coverage ci clean

help:
	@echo "Common targets:"
	@echo "  make lint      - Run Ruff lint (same as CI)"
	@echo "  make test      - Run pytest (quiet)"
	@echo "  make coverage  - Run tests under coverage and produce reports"
	@echo "  make ci        - Run the same steps as GitHub Actions (lint + coverage)"
	@echo "  make clean     - Remove caches and coverage artifacts"

lint:
	ruff check .

test:
	pytest -q

coverage:
	coverage run -m pytest -q
	coverage report -m
	coverage xml -o coverage.xml

ci: lint coverage

clean:
	rm -rf .pytest_cache .coverage coverage.xml
