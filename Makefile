.PHONY: test experiments all clean

test:
	python -m pytest

experiments:
	python experiments/run_all.py

all: test experiments

clean:
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
	rm -rf .pytest_cache src/*.egg-info
