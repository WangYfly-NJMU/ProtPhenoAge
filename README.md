# ProtPhenoAge

[中文](README_CN.md) | **English**

> A pipeline for computing protein-based phenotypic age (ProtPhenoAge) from plasma proteomics
> vbeta.20261116 (2026-11-16)

ProtPhenoAge predicts phenotypic age (PhenoAge) from Olink plasma protein measurements. The corresponding paper was published in *Journal of Advanced Research* (2026). This Python pipeline provides three workflows:

1. **UKB reproduction**: use the trained UK Biobank protein phenotypic age model, or reproduce its training on UKB data
2. **Independent-cohort training**: train a model on your own cohort (BorutaShap feature selection + Optuna tuning + XGBoost)
3. **Independent validation**: load a trained model and predict ProtPhenoAge for any new protein matrix

The UKB model under `model/` can be loaded for prediction directly, with no training data required.

## Processing of Olink Protein Matrix before input

1. **Long-to-wide**: spread Olink `result` by `protein_id`
2. **Deduplication and repeated measurements**: deduplicate by `eid`; for participants with repeated measurements, take the first non-missing value per protein
3. **Protein QC**: remove proteins failing quality control (`GLIPR1` in UKB); remove proteins with > 20% missing values (keep proteins with > 80% non-missing)
4. **Sample QC**: retain White ethnicity only (optional); exclude baseline cancer (optional)
5. **Missing-value imputation**: KNN (k = 10, scale = TRUE, weighted average)
6. **Normalization**: rank-based inverse normal transformation (RINT) per protein
7. **PhenoAge calculation**: when training on your own cohort, compute PhenoAge in the training set first; see [BioAge](https://github.com/dayoonkwon/BioAge/)

> Olink panels measure relative protein abundance, so internal retraining within each dataset is recommended. When that is not feasible, the trained UKB model (or its parameters) can be used directly.

## Repository layout

| Path | Description |
|---|---|
| `ProtPhenoAge4github.qmd` | Pipeline specification (Quarto) |
| `0.0 qc.r` | Reference UKB preprocessing implementation (R; paths are historical) |
| `1.0 train_model.py` | UKB reproduction training: internal-validation early stopping → refit on full data → save model/metadata/test predictions |
| `script/ProtPhenoAge4UKB.py` | Inference: load the UKB model and predict on any protein matrix |
| `script/ProtPhenoAgeTrain.py` | Independent-cohort training: feature selection → tuning → training → model saving |
| `model/` | Trained UKB model (`ProtPhenoAge_UKB.json`) and metadata |
| `data/` | Training/test data (UKB raw data are subject to data-use agreements and are not redistributed here) |

## Requirements

```bash
pip install pandas numpy scipy scikit-learn xgboost
# Optional, for the default independent-training workflow:
pip install optuna BorutaShap
```

## Quick start

### 1. Predict with the pretrained UKB model

```bash
python script/ProtPhenoAge4UKB.py \
  --input data/test_prot_1000.csv \
  --model model/ProtPhenoAge_UKB.json \
  --output output/test_prot_1000_protphenoage.csv
```

- `--input` is required; `--model` defaults to `model/ProtPhenoAge_UKB.json`; `--output` defaults to `<input_stem>_protphenoage.csv` (next to the input)
- Preserves an ID column when present; if a `phenoage` column exists, `protPA_accel` (PhenoAge − ProtPhenoAge) is also output
- `--keep-all`: output the full input matrix with an appended `protphenoage` column
- The model contains 185 protein features aligned by column name; missing model proteins are treated as missing values, with a warning

You can also load the model file directly in Python:

```python
import xgboost as xgb

booster = xgb.Booster(model_file='model/ProtPhenoAge_UKB.json')
features = booster.feature_names
pred = booster.predict(xgb.DMatrix(df.reindex(columns=features),
                                   feature_names=features))
```

> Input protein values must be on the same preprocessing scale as the training data (missing-value imputation + RINT within each cohort). Because Olink is relatively quantified, direct cross-cohort transfer may carry systematic bias; internal retraining within your own cohort is recommended when feasible.

### 2. Train on an independent cohort

- Coming soon

## Input data format

All input files are **CSV**: comma-separated (`,`), UTF-8 encoded, with the column names in the first row. Matrices follow the convention "one row per sample, one column per variable".

### Protein matrix (prediction `--input` / training `--x`)

- **Rows**: one sample (participant) per row; row order defines sample order
- **Row names (first column)**: a row-index column, header usually empty (`""`, read by pandas as `Unnamed: 0`). It may also be a sample ID column (e.g., `eid`, `f.eid`, `sample_id`), which is detected automatically and excluded from features. The column is optional (e.g., `data/test_prot_1000.csv`)
- **Columns**: one protein per remaining column, named by Olink protein symbol (e.g., `A1BG`, `GDF15`, `IL6`). Proteins are aligned by **column name**, so column order does not matter; extra columns are ignored and missing proteins are treated as missing values
- **Cells**: protein abundance values (Olink NPX), which must be on the same scale as the training data, i.e., after QC, missing-value imputation and RINT; missing values may be left blank or written as `NA`
- **Delimiter**: comma (`,`)

Example (from `data/x_train.csv`; 2,920 protein columns in total):

```csv
,A1BG,AAMDC,AARSD1,ABCA2,...,ZPR1
31180,0.0678,-0.95205,-0.91035,-0.1899,...
```

### Label file (training/evaluation `--y`)

- **Rows**: aligned one-to-one with the protein matrix. Without a shared ID, rows are aligned by order (both files must have the same number of rows); with a shared ID, samples are matched by ID
- **Columns**:
  - First column: row index (empty header; may be omitted)
  - `f.eid`: sample/participant identifier
  - `phenoage`: phenotypic age label (numeric, years), computed in advance with the [BioAge](https://github.com/dayoonkwon/BioAge/) method
- **Delimiter**: comma (`,`)

Example (`data/y_train.csv`):

```csv
,f.eid,phenoage
31180,5184712,57.7314167660788
```

### Prediction output

The output is `<input_stem>_protphenoage.csv` (CSV, comma-separated, with header):

- The first column is the input ID column when present, otherwise `row` (1-based row number)
- `protphenoage`: predicted protein phenotypic age (years)
- When the input has a `phenoage` column, `protPA_accel = phenoage − protphenoage` is also output

Example:

```csv
row,protphenoage
1,41.804337
2,56.852463
```

## Citation

If this pipeline is useful for your research, please cite:

> Wang Y, Sun Y, Yang F, et al. **ProtPhenoAge: Integrating plasma proteomics to predict Aging-Related disease Risks**. *Journal of Advanced Research*, 2026. doi:10.1016/j.jare.2026.05.022

[Paper link (ScienceDirect)](https://www.sciencedirect.com/science/article/pii/S2090123226004170)

```bibtex
@article{wang2026protophenoage,
  title   = {ProtPhenoAge: Integrating plasma proteomics to predict Aging-Related disease Risks},
  author  = {Wang, Yifan and Sun, Yue and Yang, Fan and others},
  journal = {Journal of Advanced Research},
  year    = {2026},
  doi     = {10.1016/j.jare.2026.05.022}
}
```

## Acknowledgements

- PhenoAge calculation: [BioAge](https://github.com/dayoonkwon/BioAge/)
- Feature selection: [BorutaShap](https://github.com/Ekeany/Boruta-Shap)
- Tuning and model: Optuna, XGBoost

## Contact

Wang Yifan — wangyfly123@163.com
