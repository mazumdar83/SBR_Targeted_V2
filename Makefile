.PHONY: install test lint smoke clean

install:
	pip install -e ".[dev]"

test:
	pytest -q tests/

lint:
	ruff check src/ tests/

smoke:
	python -m pytest tests/test_smoke_live.py -v -m live

clean:
	rm -rf data/cache/* build/ dist/ *.egg-info
