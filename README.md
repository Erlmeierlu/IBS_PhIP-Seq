# PhIP-Seq Antibody Repertoire Analysis in Irritable Bowel Syndrome (IBS)

This directory contains the analysis scripts used to process, analyze, and visualize phage-display immunoprecipitation sequencing (PhIP-Seq) antibody repertoire data in individuals with Irritable Bowel Syndrome (IBS) and population controls from the **Lifelines-DEEP** cohort.

---

## 1. Data Availability & Controlled Access

The PhIP-Seq profiles used in this study were originally described in:

1. **Andreu-Sánchez, S., Bourgonje, A. R., et al.** *Phage display sequencing reveals that genetic, environmental, and intrinsic factors influence variation of human antibody epitope repertoire.* Immunity 56, 1376–1392 (2023).

Data are available via managed access in the European Genome-Phenome Archive (EGA) under study accession number EGAS00001006999 ([https://ega-archive.org/studies/EGAS00001006999]). 

Phenotype and clinical data belong to the **Lifelines Biobank** and cannot be distributed directly in a public repository due to participant privacy regulations. Researchers can request access to the underlying Lifelines-DEEP datasets through the Lifelines Research Office (https://www.lifelines.nl/researcher(https://www.lifelines.nl/researcher)).


---

## 2. Pipeline Execution Order

Run the scripts from this directory, in numerical order:

1. **`0_utils.R`** *(sourced automatically by scripts 1–6)*
   - Defines the publication `ggplot2` theme (`theme_my()`) and the shared seroprevalence function (`calc_prev()`).

2. **`1_QC.R`**
   - Reads the baseline and follow-up PhIP-Seq tables, labels follow-up IBS trajectories, binarizes the matrix via `Rcpp`, filters low-richness samples (`< 200` peptides), filters peptides outside the `5%–95%` seroprevalence window, and deduplicates identical amino acid sequences.
   - **Outputs**: `1_phenos_ibs.rds`, `1_phenos_fup.rds`, `pi_ibs_markers.csv`, `fig1_qc_panels.rds`, `extended_figure_1.pdf`.

3. **`2_Exploratory.R`**
   - Performs confounder screening across study groups, builds the baseline characteristics table (`baseline_characteristics.xlsx`), runs mode-imputed PCA and k-means clustering (`k = 2`), tests PC1/PC2 associations with clinical and technical variables, and assembles **Figure 1**.
   - **Outputs**: `2_phenos_subtypes.rds`, `food_categories.csv`, `baseline_characteristics.xlsx`, `figure_1.pdf`.

4. **`3_Association.R`**
   - Runs per-peptide Firth's bias-reduced logistic regression (`logistf`) adjusting for Age, Sex, and sequencing plate across IBS subtypes, general IBS, baseline controls developing IBS at follow-up, IBS excluding IBS-M, Age quartiles, and Sex.
   - **Outputs**: `Association_*.csv`, `extended_table_9.xlsx`.

5. **`4_Prevalence_and_Enrichment.R`**
   - Conducts taxon-level and source-database seroprevalence ratio testing and fast Gene Set Enrichment Analysis (`fgsea`) across cross-sectional and longitudinal cohorts, pairwise Phi-coefficient repertoire stability permutation testing, PI-IBS and food allergen enrichment, and assembles **Figure 2** and **Figure 3**.
   - **Outputs**: `figure_2.pdf`, `figure_3.pdf`, `extended_figure_3.pdf`, `extended_figure_5.pdf`, `extended_figure_6.pdf`, `extended_figure_7.pdf`, and `extended_table_2.xlsx` through `extended_table_5.xlsx`.

6. **`5_Run_Classification.R`**
   - Trains, tunes (Latin hypercube grid + Bayesian optimization minimizing log-loss), and evaluates gradient-boosted decision trees (`xgboost` via `tidymodels`) with out-of-fold F1 threshold tuning and `shapviz` SHAP attribution for:
     1. Binary IBS Subtypes (`IBSC`, `IBSD`, `IBSM` vs `no IBS`)
     2. Binary General IBS (`IBS` vs `no IBS`)
     3. Multiclass IBS Subtypes (`noIBS`, `IBSC`, `IBSD`, `IBSM`)
   - **Outputs**: `best_models_logloss.rds`, `metrics_all_recipes_logloss.rds`, `best_IBS_models_logloss.rds`, `metrics_IBS_all_recipes_logloss.rds`, `models_multiclass_logloss.rds`, and `metrics_multiclass_logloss.rds`.

7. **`6_Visualize_ML.R`**
   - Computes stratified bootstrap confidence intervals for ROC-AUC and PR-AUC, generates **Figure 4** (`figure_4a.pdf`, `figure_4b.pdf`, `figure_4c.pdf`) and **Extended Data Figure 8** (`extended_figure_8.pdf`), and exports model hyperparameter, performance, and SHAP summary tables (`extended_table_1.xlsx`, `extended_table_6.xlsx`, `extended_table_7.xlsx`).
