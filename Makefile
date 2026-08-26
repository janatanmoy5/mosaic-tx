.PHONY: data test demo check
data:
	python3 scripts/build_public_data.py --dataset GSE50084
test:
	python3 -m unittest discover -s tests -v
demo:
	python3 -m mosaic_tx.cli rank --disease examples/synthetic/disease_signature.csv --perturbations examples/synthetic/perturbations.csv --safety examples/synthetic/safety.csv --output results/demo_rankings.csv --audit results/demo_audit.json
check: test demo
	@test -f data/processed/GSE50084/manifest.json || (echo "Run make data first"; exit 1)
	python3 -c 'from pathlib import Path; from mosaic_tx.geo import validate_compiled; validate_compiled(Path("data/processed/GSE50084"))'
