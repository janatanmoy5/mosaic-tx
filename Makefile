.PHONY: all data annotate bulk scores pilot manuscript-figure validate test check

all: data annotate bulk scores pilot manuscript-figure check validate

data:
	python3 scripts/build_public_data.py --dataset GSE50084
	python3 scripts/build_bulk_cohorts.py --dataset GSE36059
	python3 scripts/build_bulk_cohorts.py --dataset GSE275824

annotate:
	Rscript scripts/02b_annotate_bioc.R data/processed/GSE50084/expression_probe_matrix.csv.gz data/processed/GSE50084/samples.csv GPL6244 data/processed/GSE50084
	Rscript scripts/02b_annotate_bioc.R data/processed/GSE36059/expression_probe_matrix.csv.gz data/processed/GSE36059/samples.csv GPL570 data/processed/GSE36059
	Rscript scripts/02b_annotate_bioc.R data/processed/GSE275824/expression_probe_matrix.csv.gz data/processed/GSE275824/samples.csv GPL13667 data/processed/GSE275824

bulk:
	Rscript scripts/03_bulk_de_gsea.R data/processed/GSE50084/gene_expression.csv.gz data/processed/GSE50084/samples.csv results/bulk/GSE50084
	Rscript scripts/04_external_validation.R data/processed/GSE36059/gene_expression.csv.gz data/processed/GSE36059/samples.csv results/bulk/GSE36059
	Rscript scripts/05_treatment_response.R data/processed/GSE275824/gene_expression.csv.gz data/processed/GSE275824/samples.csv results/bulk/GSE275824

scores:
	Rscript scripts/06_published_score_trajectories.R data/processed/GSE275824/samples.csv results/bulk/GSE275824

pilot:
	Rscript scripts/08_build_robustness_atlas.R . results/pilot/atlas
	Rscript scripts/09_make_pilot_figures.R results/pilot/atlas results/pilot/figures

manuscript-figure:
	Rscript scripts/09b_make_manuscript_figure.R . results/pilot/atlas results/pilot/figures

validate:
	python3 scripts/14_validate_pipeline_outputs.py --root .

test:
	python3 -m unittest discover -s tests -v
	python3 -m py_compile scripts/07_score_mechanism_modules.py scripts/10_build_scanpy_object.py scripts/11_sc_qc_cluster.py scripts/12_patient_pseudobulk.py scripts/13_grouped_nested_ml.py
	@for f in scripts/*.R; do Rscript -e "parse(file='$$f')" >/dev/null || exit 1; done

check: test
	python3 -c 'from pathlib import Path; from mosaic_tx.geo import validate_compiled; validate_compiled(Path("data/processed/GSE50084"))'
	@test -f results/bulk/GSE275824/published_score_manifest.json
