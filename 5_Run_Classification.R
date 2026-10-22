# ==============================================================================
# 5_Run_Classification.R
# XGBoost classification pipelines (Binary Subtypes, Binary Any-IBS, Multiclass)
# ==============================================================================
#
# Purpose:
#   Combines the three gradient-boosted decision tree (XGBoost) classification
#   analyses into a unified, reproducible pipeline:
#     1. Binary IBS Subtype Classifiers (`IBSC`, `IBSD`, `IBSM` vs `no IBS`):
#        - Stratified 80/20 train/test split and repeated stratified k-fold CV
#          (k = 5/4/5, repeats = 7/9/7 to maintain >= 5 positive cases per fold).
#        - Evaluates 5 preprocessing recipes across covariate-only (`covars`),
#          peptide-only (`features`), and combined (`both`) feature sets.
#        - Optimizes hyperparameters via Latin hypercube space-filling grid
#          (30 candidates) followed by Bayesian optimization (22 iterations)
#          minimizing log-loss (`mn_log_loss`).
#        - Tunes two-class probability thresholds on out-of-fold (OOF) training
#          predictions maximizing F1-score (`f_meas`), avoiding test-set leakage.
#        - Evaluates finalized models on the held-out test split and computes
#          SHAP feature attributions (`shapviz`).
#     2. Binary General IBS Classifier (`IBS` vs `no IBS`):
#        - Same binary log-loss and OOF threshold-tuning pipeline with 5x7
#          repeated stratified CV.
#     3. Multiclass IBS Subtype Classifier (`noIBS`, `IBSC`, `IBSD`, `IBSM`):
#        - Uses inverse class-frequency `importance_weights` to account for
#          class imbalance across the 4 groups.
#        - Trains multiclass XGBoost (`objective = 'multi:softprob'`,
#          `eval_metric = 'mlogloss'`) with 5x10 repeated stratified CV.
#
# Expected Inputs:
#   - `annotation.csv`:
#     Peptide annotation table (`peptide`, `taxa`, `source_database`, `description`).
#   - `1_phenos_ibs.rds`: QC'd baseline PhIP-Seq and phenotype data.table.
#
# Generated Outputs:
#   - `best_models_logloss.rds`
#   - `metrics_all_recipes_logloss.rds`
#   - `best_IBS_models_logloss.rds`
#   - `metrics_IBS_all_recipes_logloss.rds`
#   - `metrics_multiclass_logloss.rds`
#   - `models_multiclass_logloss.rds`
# ==============================================================================

library(data.table)
library(tidymodels)
library(shapviz)
library(readr)
library(forcats)
library(themis)
library(purrr)
library(butcher)
library(colino)
library(future)
library(probably)
library(fs)

source("0_utils.R")

# 1. Reproducibility & Parallel Setup ------------------------------------------
# RNGkind("L'Ecuyer-CMRG") provides independent, seedable random substreams across
# parallel workers. Must be set prior to spawning workers via future::plan().
RNGkind("L'Ecuyer-CMRG")

workers <- 40
plan(multicore, workers = workers)

seed_main      <- 1337L
seed_split     <- 2024L
seed_folds     <- 4242L
seed_grid_init <- 5150L
seed_bayes     <- 6160L
seed_lastfit   <- 7070L
set.seed(seed_main)


# 2. Load Input Data -----------------------------------------------------------

anno_path   <- "annotation.csv"
phenos_path <- "1_phenos_ibs.rds"

annotation <- fread(anno_path)
phenos     <- read_rds(phenos_path)


# 3. Shared Helper Functions ---------------------------------------------------

# Disk-backed cache helper per classification task
cache_or <- function(name, expr, cache_dir, refresh = FALSE) {
  fs::dir_create(cache_dir)
  p <- file.path(cache_dir, paste0(name, ".rds"))
  if (!refresh && file.exists(p)) {
    message("[cache] reading ", p)
    return(read_rds(p))
  }
  message("[cache] computing ", name)
  val <- force(expr)
  write_rds(val, p)
  val
}

# Filter and format dataset for tidymodels classification
prepare_data_ml <- function(
  data,
  covariates,
  filter_by,
  filter_for = NULL,
  subset_col = NULL,
  subset = NULL,
  metadata = annotation
) {
  data <- copy(data)
  if (!is.null(subset)) {
    metadata <- metadata[get(subset_col) %in% subset]
  }
  peptide_cols <- metadata[peptide %in% names(data), unique(peptide)]

  data <- data[, c(covariates, filter_by, peptide_cols), with = FALSE]
  setnames(data, filter_by, "classification")

  data[, .SDcols = is.character, names(.SD) := lapply(.SD, as.factor)]

  if (!is.null(filter_for)) {
    data <- data[classification %in% filter_for] |>
      droplevels()
  }
  data[, classification := fct_relabel(classification, ~ gsub("\\-| ", "", .x))]
  data[, classification := fct_rev(classification)]
  data[]
}

# Check whether selected hyperparameters hit the boundary of the tuning range
warn_if_at_boundary <- function(ps, best, tol = 1e-6, rel = 0.001) {
  rngs <- map2_dfr(ps$object, ps$id, function(p, i) {
    r <- dials::range_get(p)
    tibble(id = i, lower = r$lower, upper = r$upper)
  })

  rngs <- rngs |> filter(id %in% names(best))

  vals <- best |>
    select(any_of(rngs$id)) |>
    dplyr::slice(1) |>
    tidyr::pivot_longer(everything(), names_to = "id", values_to = "value")

  out <- rngs |>
    left_join(vals, by = "id") |>
    mutate(
      width       = pmax(upper - lower, .Machine$double.eps),
      near_low    = (abs(value - lower) <= pmax(tol, rel * width)),
      near_high   = (abs(value - upper) <= pmax(tol, rel * width)),
      at_boundary = near_low | near_high
    )

  to_test <- out |> filter(id != "mtry")
  if (any(to_test$at_boundary)) {
    which <- out |>
      filter(at_boundary) |>
      mutate(side = ifelse(near_low, "min", "max")) |>
      transmute(
        msg = paste0(
          id,
          " at ",
          side,
          " (value=",
          signif(value, 3),
          ", range=[",
          signif(lower, 3),
          ", ",
          signif(upper, 3),
          "])"
        )
      ) |>
      pull(msg)
    warning(
      paste(c("Best config hit parameter boundary:", which), collapse = "\n - ")
    )
  }
  out
}

## 3a. Recipe Builders ---------------------------------------------------------
make_base <- function(data, covariates) {
  recipe(classification ~ ., data = data) |>
    update_role(
      -c(classification, has_role("case_weights")),
      new_role = "predictor"
    ) |>
    update_role(classification, new_role = "outcome") |>
    add_role(all_of(covariates), new_role = "covariates")
}

step_mod_impute <- function(rec) {
  rec |>
    step_impute_knn(c(!has_role("covariates") & has_role("predictor"))) |>
    step_impute_mean(c(has_role("covariates") & has_type("numeric")))
}

step_mod_norm <- function(rec) {
  rec |>
    step_normalize(has_role("covariates") & has_type("numeric"))
}

step_mod_only_base <- function(rec) {
  rec |>
    step_rm(
      !c(
        has_role("covariates") |
          has_role("outcome") |
          has_role("case_weights")
      )
    )
}

step_mod_only_features <- function(rec) {
  rec |>
    step_rm(has_role("covariates"))
}

step_mod_corr <- function(rec) {
  rec |>
    step_corr(
      all_numeric_predictors() & !has_role("covariates"),
      threshold = 0.9
    )
}

make_rec <- function(
  use,
  impute,
  balance_method,
  fisher,
  feat_roc,
  norm,
  corr,
  rec_name,
  data,
  covariates
) {
  rec <- make_base(data, covariates)
  if (use == "both") {
    rec <- rec |>
      step_dummy(all_nominal_predictors(), one_hot = FALSE)
  }
  if (use == "covars") {
    rec <- rec |>
      step_mod_only_base() |>
      step_dummy(all_nominal_predictors(), one_hot = FALSE)
  }
  if (use == "features") {
    rec <- rec |>
      step_mod_only_features()
  }
  rec <- rec |> step_zv(all_predictors())
  if (impute) {
    rec <- step_mod_impute(rec)
  }
  if (norm) {
    rec <- step_mod_norm(rec)
  }
  if (corr) {
    rec <- step_mod_corr(rec)
  }
  rec
}

# Candidate recipe grid shared across all classification tasks
opts <- crossing(
  use = c("both", "covars", "features"),
  impute = c(TRUE, FALSE),
  balance_method = c("none"),
  fisher = c(FALSE),
  feat_roc = c(FALSE),
  norm = TRUE,
  corr = c(TRUE, FALSE)
) |>
  filter(
    !(use == "covars" & (impute | fisher | feat_roc | balance_method != "none" | !norm)),
    !(use == "both" & !norm),
    !(use == "covars" & corr),
    !(use != "covars" & !corr)
  ) |>
  mutate(
    rec_name = paste0(
      "Use", use,
      "_Imp", impute,
      "_Bal", balance_method,
      "_Fisher", fisher,
      "_Fs", feat_roc,
      "_Norm", norm,
      "_Corr", corr
    )
  )

## 3b. Hyperparameter Search Space & Tuning Controls ---------------------------
add_param_ranges <- function(wflow, mtry_cap) {
  ps <- hardhat::extract_parameter_set_dials(wflow)
  ps <- update(
    ps,
    trees = trees() |> range_set(c(300L, 3500L)),
    tree_depth = tree_depth(range = c(2L, 8L)),
    learn_rate = learn_rate(range = c(-3, -1)),
    min_n = min_n(range = c(2L, 80L)),
    loss_reduction = loss_reduction(range = c(-10, 1)),
    mtry = mtry(range = c(1L, mtry_cap))
  )
  if ("threshold" %in% ps$id) {
    ps <- update(
      ps,
      threshold = threshold() |> range_set(c(0.55, 0.85))
    )
  }
  ps
}

ctrl_grid_obj <- control_grid(
  save_pred     = TRUE,
  save_workflow = FALSE,
  parallel_over = "everything",
  verbose       = FALSE
)

ctrl_bayes_obj <- control_bayes(
  save_pred     = TRUE,
  save_workflow = FALSE,
  parallel_over = "everything",
  no_improve    = 12L,
  verbose_iter  = TRUE
)

tune_one <- function(
  wf,
  param_set,
  folds,
  metrics,
  seed_init = seed_grid_init,
  seed_b = seed_bayes
) {
  set.seed(seed_init)
  init_grid <- grid_space_filling(param_set, size = 30, type = "latin_hypercube")

  set.seed(seed_init)
  initial_rs <- tune_grid(
    wf,
    resamples = folds,
    grid      = init_grid,
    metrics   = metrics,
    control   = ctrl_grid_obj
  )

  set.seed(seed_b)
  tune_bayes(
    wf,
    resamples  = folds,
    param_info = param_set,
    initial    = initial_rs,
    iter       = 22,
    metrics    = metrics,
    control    = ctrl_bayes_obj
  )
}

# Tune probability decision threshold on out-of-fold training predictions
choose_threshold <- function(
  tune_obj,
  best,
  st,
  thresholds = seq(0.05, 0.95, by = 0.005)
) {
  oof <- collect_predictions(tune_obj, parameters = best, summarize = FALSE)
  est_col <- paste0(".pred_", st)
  probably::threshold_perf(
    oof,
    truth      = classification,
    estimate   = !!rlang::sym(est_col),
    metrics    = metric_set(f_meas),
    thresholds = thresholds
  ) |>
    dplyr::slice_max(.estimate, n = 1, with_ties = FALSE) |>
    dplyr::pull(.threshold)
}


# 4. Binary Classification Runner (Subtypes & General IBS) ---------------------

run_binary_classification <- function(
  sub_data,
  k_by_sub,
  rep_by_sub,
  cache_dir,
  metrics_out_path,
  models_out_path
) {
  # 4a. Stratified train/test splits and repeated CV folds
  splits_folds <- cache_or("splits_folds", {
    set.seed(seed_split)
    split_tbl <- map(
      sub_data,
      ~ initial_split(.x, strata = "classification", prop = .8)
    )
    train_tbl <- map(split_tbl, training)
    test_tbl  <- map(split_tbl, testing)

    set.seed(seed_folds)
    folds_tbl <- imap(
      train_tbl,
      ~ vfold_cv(
        .x,
        v = k_by_sub[.y],
        repeats = rep_by_sub[.y],
        strata = classification
      )
    )

    list(
      split_tbl = split_tbl,
      train_tbl = train_tbl,
      test_tbl  = test_tbl,
      folds_tbl = folds_tbl
    )
  }, cache_dir = cache_dir)

  split_tbl <- splits_folds$split_tbl
  train_tbl <- splits_folds$train_tbl
  folds_tbl <- splits_folds$folds_tbl

  # 4b. Binary XGBoost specification & workflow set
  mod_bin <- boost_tree(
    trees = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    min_n = tune(),
    mtry = tune(),
    loss_reduction = tune()
  ) |>
    set_engine(
      "xgboost",
      nthread = 1,
      objective = "binary:logistic",
      max_delta_step = 5,
      eval_metric = "logloss"
    ) |>
    set_mode("classification")

  wf_sets <- imap(train_tbl, function(train, st, cov = c("Age", "Sex")) {
    rec_list <- pmap(opts, .f = make_rec, data = train, covariates = cov) |>
      set_names(opts$rec_name)
    models <- map(opts$balance_method, \(bm) mod_bin)

    workflow_set(rec_list, models, cross = FALSE) |>
      mutate(
        mtry_cap = if_else(opts$use == "covars", 2, ncol(train) - 1),
        use = opts$use
      )
  }) |>
    bind_rows(.id = "subtype")

  metrics_bin <- metric_set(mn_log_loss, pr_auc, roc_auc, brier_class)
  metric_to_optimize <- "mn_log_loss"

  wf_sets <- wf_sets |>
    mutate(wf = map(info, "workflow")) |>
    mutate(param_set = map2(wf, mtry_cap, ~ add_param_ranges(.x[[1]], .y)))

  # 4c. Hyperparameter tuning
  tuned_sets <- cache_or("tuned_sets", {
    wf_sets |>
      mutate(
        tuning = pmap(
          list(wf, param_set, subtype),
          \(wflow, ps, st) tune_one(wflow[[1]], ps, folds_tbl[[st]], metrics_bin)
        )
      )
  }, cache_dir = cache_dir)

  selected <- tuned_sets |>
    mutate(
      best = map(tuning, select_best, metric = metric_to_optimize),
      best_ll = map2_dbl(
        tuning,
        best,
        ~ collect_metrics(.x) |>
          dplyr::filter(.config == .y$.config, .metric == "mn_log_loss") |>
          dplyr::pull(mean)
      ),
      pr_at_best_ll = map2_dbl(
        tuning,
        best,
        ~ collect_metrics(.x) |>
          dplyr::filter(.config == .y$.config, .metric == "pr_auc") |>
          dplyr::pull(mean)
      ),
      roc_at_best_ll = map2_dbl(
        tuning,
        best,
        ~ collect_metrics(.x) |>
          dplyr::filter(.config == .y$.config, .metric == "roc_auc") |>
          dplyr::pull(mean)
      ),
      boundary_check = map2(param_set, best, warn_if_at_boundary)
    )

  # 4d. Finalize workflows and evaluate on held-out test set
  evaluated <- cache_or("evaluated", {
    selected |>
      mutate(
        wf_final = map2(wf, best, ~ finalize_workflow(.x[[1]], .y)),
        fin_fit = map2(wf_final, subtype, function(.x, .y) {
          set.seed(seed_lastfit)
          last_fit(.x, split_tbl[[.y]], metrics = metrics_bin)
        })
      )
  }, cache_dir = cache_dir)

  # 4e. OOF threshold tuning and test set scoring
  scored <- evaluated |>
    mutate(
      raw_preds = map(fin_fit, collect_predictions),
      best_thresh = pmap_dbl(
        list(tuning, best, subtype),
        ~ choose_threshold(..1, ..2, ..3)
      ),
      preds = pmap(
        list(raw_preds, best_thresh, subtype),
        function(.x, .y, .z) {
          est_col <- paste0(".pred_", .z)
          tibble(
            "{est_col}" := .x[[est_col]],
            .pred_class = probably::make_two_class_pred(
              estimate  = .x[[est_col]],
              levels    = levels(.x$classification),
              threshold = .y
            ) |>
              as.factor(),
            classification = .x$classification
          )
        }
      ),
      conf_mat = map(
        preds,
        ~ conf_mat(.x, truth = classification, estimate = .pred_class)
      ),
      mets = map(fin_fit, collect_metrics),
      test_ll  = map_dbl(mets, ~ .x |> filter(.metric == "mn_log_loss") |> pull(.estimate)),
      test_roc = map_dbl(mets, ~ .x |> filter(.metric == "roc_auc")     |> pull(.estimate)),
      test_pr  = map_dbl(mets, ~ .x |> filter(.metric == "pr_auc")      |> pull(.estimate))
    )

  # 4f. Save all-recipe summary metrics
  wf_out <- scored |>
    select(
      subtype, wflow_id, use,
      best_ll, pr_at_best_ll, roc_at_best_ll, best_thresh,
      test_ll, test_roc, test_pr,
      mets, conf_mat, boundary_check
    )

  fs::dir_create(dirname(metrics_out_path))
  write_rds(wf_out, metrics_out_path)

  # 4g. Select best recipe per (subtype, use), compute SHAP values, and export
  best_sets <- scored |>
    slice_min(order_by = best_ll, by = c(subtype, use), n = 1, with_ties = FALSE) |>
    mutate(
      shap_vals = map2(
        fin_fit,
        subtype,
        \(x, y) {
          train_df <- train_tbl[[y]]
          rec <- extract_recipe(x)
          X_pred <- bake(
            rec,
            has_role("predictor"),
            new_data = train_df,
            composition = "matrix"
          )
          train_df <- train_df |> rename(Sex_M = "Sex")

          workflows::extract_fit_parsnip(x)$fit |>
            shapviz(X_pred = X_pred, X = train_df)
        }
      ),
      model = map(fin_fit, ~ extract_workflow(.x) |> butcher())
    ) |>
    select(
      subtype, wflow_id, use,
      best_ll, pr_at_best_ll, roc_at_best_ll, best_thresh,
      test_ll, test_roc, test_pr,
      mets, preds, conf_mat, shap_vals, model, boundary_check
    )

  fs::dir_create(dirname(models_out_path))
  write_rds(best_sets, models_out_path)
  invisible(best_sets)
}


# 5. Multiclass Classification Runner (`noIBS`, `IBSC`, `IBSD`, `IBSM`) --------

run_multiclass_classification <- function(
  phenos_dt,
  cache_dir,
  metrics_out_path,
  models_out_path
) {
  # 5a. Prepare 4-class dataset with inverse class-frequency importance weights
  full_data <- prepare_data_ml(
    phenos_dt,
    covariates = c("Age", "Sex"),
    filter_by = "IBS_Subtype"
  )

  full_data <- full_data |>
    mutate(classification = relevel(full_data$classification, "noIBS"))

  class_counts <- table(full_data$classification)
  total_count  <- sum(class_counts)
  weights      <- total_count / (length(class_counts) * class_counts)

  full_data <- full_data |>
    mutate(imp_weight = as.double(weights[as.character(classification)])) |>
    mutate(imp_weight = importance_weights(imp_weight))

  # 5b. Stratified train/test split and 5x10 repeated CV folds
  splits_folds <- cache_or("splits_folds", {
    set.seed(seed_split)
    split_tbl <- initial_split(full_data, strata = "classification", prop = .8)
    train_tbl <- training(split_tbl)
    test_tbl  <- testing(split_tbl)

    set.seed(seed_folds)
    folds_tbl <- vfold_cv(
      train_tbl,
      v = 5,
      repeats = 10,
      strata = classification
    )

    list(
      split_tbl = split_tbl,
      train_tbl = train_tbl,
      test_tbl  = test_tbl,
      folds_tbl = folds_tbl
    )
  }, cache_dir = cache_dir)

  split_tbl <- splits_folds$split_tbl
  train_tbl <- splits_folds$train_tbl
  folds_tbl <- splits_folds$folds_tbl

  # 5c. Multiclass XGBoost specification & workflow set with case_weights
  mod_multi <- boost_tree(
    trees = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    min_n = tune(),
    mtry = tune(),
    loss_reduction = tune()
  ) |>
    set_engine(
      "xgboost",
      nthread = 1,
      objective = "multi:softprob",
      max_delta_step = 5,
      eval_metric = "mlogloss"
    ) |>
    set_mode("classification")

  rec_list <- pmap(
    opts,
    .f = make_rec,
    data = train_tbl,
    covariates = c("Age", "Sex")
  ) |>
    set_names(opts$rec_name)
  models <- map(opts$balance_method, \(bm) mod_multi)

  wf_sets <- workflow_set(
    rec_list,
    models,
    cross = FALSE,
    case_weights = imp_weight
  ) |>
    mutate(
      mtry_cap = if_else(opts$use == "covars", 2, ncol(train_tbl) - 1),
      use = opts$use
    )

  metrics_multi <- metric_set(mn_log_loss, roc_auc, pr_auc, accuracy)
  metric_to_optimize <- "mn_log_loss"

  wf_sets <- wf_sets |>
    mutate(wf = map(info, "workflow")) |>
    mutate(param_set = map2(wf, mtry_cap, ~ add_param_ranges(.x[[1]], .y)))

  # 5d. Hyperparameter tuning
  tuned_sets <- cache_or("tuned_sets", {
    wf_sets |>
      mutate(
        tuning = pmap(
          list(wf, param_set),
          \(wflow, ps) tune_one(wflow[[1]], ps, folds_tbl, metrics_multi)
        )
      )
  }, cache_dir = cache_dir)

  selected <- tuned_sets |>
    mutate(
      best = map(tuning, select_best, metric = metric_to_optimize),
      best_ll = map2_dbl(
        tuning,
        best,
        ~ collect_metrics(.x) |>
          dplyr::filter(.config == .y$.config, .metric == "mn_log_loss") |>
          dplyr::pull(mean)
      ),
      pr_at_best_ll = map2_dbl(
        tuning,
        best,
        ~ collect_metrics(.x) |>
          dplyr::filter(.config == .y$.config, .metric == "pr_auc") |>
          dplyr::pull(mean)
      ),
      roc_at_best_ll = map2_dbl(
        tuning,
        best,
        ~ collect_metrics(.x) |>
          dplyr::filter(.config == .y$.config, .metric == "roc_auc") |>
          dplyr::pull(mean)
      ),
      boundary_check = map2(param_set, best, warn_if_at_boundary)
    )

  # 5e. Finalize workflows and evaluate on held-out test set
  evaluated <- cache_or("evaluated", {
    selected |>
      mutate(
        wf_final = map2(wf, best, ~ finalize_workflow(.x[[1]], .y)),
        fin_fit = map(wf_final, function(.x) {
          set.seed(seed_lastfit)
          last_fit(.x, split_tbl, metrics = metrics_multi)
        })
      )
  }, cache_dir = cache_dir)

  scored <- evaluated |>
    mutate(
      preds = map(fin_fit, collect_predictions),
      conf_mat = map(
        preds,
        ~ conf_mat(.x, truth = classification, estimate = .pred_class)
      ),
      mets = map(fin_fit, collect_metrics),
      test_ll  = map_dbl(mets, ~ .x |> filter(.metric == "mn_log_loss") |> pull(.estimate)),
      test_roc = map_dbl(mets, ~ .x |> filter(.metric == "roc_auc")     |> pull(.estimate)),
      test_pr  = map_dbl(mets, ~ .x |> filter(.metric == "pr_auc")      |> pull(.estimate))
    )

  # 5f. Save summary metrics and best multiclass models
  wf_out <- scored |>
    select(
      wflow_id, use,
      best_ll, pr_at_best_ll, roc_at_best_ll,
      test_ll, test_roc, test_pr,
      mets, conf_mat, boundary_check
    )

  fs::dir_create(dirname(metrics_out_path))
  write_rds(wf_out, metrics_out_path)

  best_sets <- scored |>
    slice_min(order_by = best_ll, by = use, n = 1, with_ties = FALSE) |>
    mutate(
      shap_vals = map(
        fin_fit,
        \(x) {
          train_df <- train_tbl
          rec <- extract_recipe(x)
          X_pred <- bake(
            rec,
            has_role("predictor"),
            new_data = train_df,
            composition = "matrix"
          )
          train_df <- train_df |> rename(Sex_M = "Sex")

          workflows::extract_fit_parsnip(x)$fit |>
            shapviz(X_pred = X_pred, X = train_df)
        }
      ),
      model = map(fin_fit, ~ extract_workflow(.x) |> butcher())
    ) |>
    select(
      wflow_id, use,
      best_ll, pr_at_best_ll, roc_at_best_ll,
      test_ll, test_roc, test_pr,
      mets, preds, conf_mat, shap_vals, model, boundary_check
    )

  fs::dir_create(dirname(models_out_path))
  write_rds(best_sets, models_out_path)
  invisible(best_sets)
}


# 6. Execute All Three Classification Tasks ------------------------------------

## 6a. Task 1: Binary Subtype Classifiers (IBSC, IBSD, IBSM vs no IBS) ---------
sub_data_subtypes <- list(
  "IBSC" = prepare_data_ml(
    phenos,
    covariates = c("Age", "Sex"),
    filter_by  = "IBS_Subtype",
    filter_for = c("no IBS", "IBS-C")
  ),
  "IBSD" = prepare_data_ml(
    phenos,
    covariates = c("Age", "Sex"),
    filter_by  = "IBS_Subtype",
    filter_for = c("no IBS", "IBS-D")
  ),
  "IBSM" = prepare_data_ml(
    phenos,
    covariates = c("Age", "Sex"),
    filter_by  = "IBS_Subtype",
    filter_for = c("no IBS", "IBS-M")
  )
)

run_binary_classification(
  sub_data         = sub_data_subtypes,
  k_by_sub         = c(IBSC = 5, IBSD = 4, IBSM = 5),
  rep_by_sub       = c(IBSC = 7, IBSD = 9, IBSM = 7),
  cache_dir        = "cache_logloss",
  metrics_out_path = "metrics_all_recipes_logloss.rds",
  models_out_path  = "best_models_logloss.rds"
)

## 6b. Task 2: Binary General IBS Classifier (IBS vs no IBS) -------------------
sub_data_ibs <- list(
  "IBS" = prepare_data_ml(
    phenos,
    covariates = c("Age", "Sex"),
    filter_by  = "IBS_Status",
    filter_for = c("no IBS", "IBS")
  )
)

run_binary_classification(
  sub_data         = sub_data_ibs,
  k_by_sub         = c(IBS = 5),
  rep_by_sub       = c(IBS = 7),
  cache_dir        = "cache_IBS_logloss",
  metrics_out_path = "metrics_IBS_all_recipes_logloss.rds",
  models_out_path  = "best_IBS_models_logloss.rds"
)

## 6c. Task 3: Multiclass Subtype Classifier (noIBS, IBSC, IBSD, IBSM) ---------
run_multiclass_classification(
  phenos_dt        = phenos,
  cache_dir        = "cache_multiclass_logloss",
  metrics_out_path = "metrics_multiclass_logloss.rds",
  models_out_path  = "models_multiclass_logloss.rds"
)
