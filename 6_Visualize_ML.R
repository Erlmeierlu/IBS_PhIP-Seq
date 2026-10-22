# ==============================================================================
# 6_Visualize_ML.R
# Machine learning evaluation, bootstrap confidence intervals, and SHAP plots
# ==============================================================================
#
# Purpose:
#   1. Loads the trained binary subtype (`IBSC`, `IBSD`, `IBSM`), binary general
#      IBS (`IBS`), and 4-class (`multiclass`) XGBoost models and metrics.
#   2. Computes test-set ROC curves, sensitivity, specificity, F1-score, and
#      stratified bootstrap (100 replicates) ROC-AUC and PR-AUC confidence
#      intervals.
#   3. Generates Manuscript Figure 4 panels:
#      - Figure 4A (`figure_4a.pdf`): Test-set ROC curves with
#        bootstrap paths comparing Full, Covariate-only, and Peptide-only models.
#      - Figure 4B (`figure_4b.pdf`): Test-set confusion matrices.
#      - Figure 4C (`figure_4c.pdf`): IBS-C SHAP beeswarm plot
#        of the top 13 predictors.
#   4. Exports Extended Data Table 7 (`extended_table_7.xlsx`):
#      Top IBS-C SHAP peptide annotations and case/control seroprevalence.
#   5. Exports Extended Data Table 1 (`extended_table_1.xlsx`):
#      XGBoost hyperparameter search space, tuned hyperparameters, and selected
#      preprocessing steps across all models.
#   6. Generates Extended Data Figure 8 (`extended_figure_8.pdf`):
#      Multiclass one-vs-all ROC curves and 4x4 confusion matrix.
#   7. Exports Extended Data Table 6 (`extended_table_6.xlsx`):
#      Complete test-set performance summary (ROC-AUC, PR-AUC, Sensitivity,
#      Specificity, F-Measure) across binary and multiclass models.
#
# Expected Inputs:
#   - `annotation.csv`: Peptide annotation table.
#   - `2_phenos_subtypes.rds`: Processed baseline dataset (from 2_Exploratory.R).
#   - `best_models_logloss.rds` & `metrics_all_recipes_logloss.rds`
#   - `best_IBS_models_logloss.rds` & `metrics_IBS_all_recipes_logloss.rds`
#   - `models_multiclass_logloss.rds` & `metrics_multiclass_logloss.rds`
#
# Generated Outputs:
#   - `figure_4a.pdf`
#   - `figure_4b.pdf`
#   - `figure_4c.pdf`
#   - `extended_figure_8.pdf`
#   - `extended_table_1.xlsx`
#   - `extended_table_6.xlsx`
#   - `extended_table_7.xlsx`
# ==============================================================================

library(data.table)
library(shapviz)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(rsample)
library(yardstick)
library(ggplot2)
library(paletteer)
library(patchwork)
library(tidytext)
library(fgsea)
library(Matrix)
library(matrixStats)
library(forcats)
library(writexl)
library(rlang)

source("0_utils.R")

set.seed(2345)
tight <- theme(plot.margin = margin(2, 2, 2, 2), panel.spacing = unit(3, "pt"))

# 1. Load Data & Saved Model Objects -------------------------------------------

meta <- fread("annotation.csv")
phenos <- read_rds("2_phenos_subtypes.rds")
anno <- fread("annotation.csv")

best_mods <- read_rds("best_models_logloss.rds")
all_mets <- read_rds("metrics_all_recipes_logloss.rds")
best_mods_ibs <- read_rds("best_IBS_models_logloss.rds")
all_mets_ibs <- read_rds("metrics_IBS_all_recipes_logloss.rds")
best_mods_multi <- read_rds("models_multiclass_logloss.rds")
all_mets_multi <- read_rds("metrics_multiclass_logloss.rds")

best_mods_multi$subtype <- "multiclass"
all_mets_multi$subtype <- "multiclass"


# 2. Shared Evaluation & Bootstrap Helper Functions ----------------------------

# Computes test-set classification metrics and stratified bootstrap ROC/PR-AUC
evaluate_and_bootstrap_models <- function(mods_df, boot_times = 100) {
  dt <- mods_df |>
    select(subtype, use, preds, auc = test_roc, pr_auc = test_pr) |>
    as.data.table()

  dt[, preds := lapply(preds, as.data.table)]

  if (!"multiclass" %in% dt$subtype) {
    dt[,
      roc := map2(
        preds,
        subtype,
        ~ roc_curve(.x, classification, paste0(".pred_", .y))
      )
    ]
    dt[,
      sens := map2(
        preds,
        subtype,
        ~ sens(.x, classification, ".pred_class") |> pull(.estimate)
      ) |>
        unlist()
    ]
    dt[,
      spec := map2(
        preds,
        subtype,
        ~ spec(.x, classification, ".pred_class") |> pull(.estimate)
      ) |>
        unlist()
    ]
    dt[,
      fmeas := map2(
        preds,
        subtype,
        ~ f_meas(.x, classification, ".pred_class") |> pull(.estimate)
      ) |>
        unlist()
    ]
  } else {
    multi_pred_cols <- function(x) {
      setdiff(grep("^\\.pred_", names(x), value = TRUE), ".pred_class")
    }
    dt[,
      roc := map(preds, \(x) {
        roc_curve(x, classification, all_of(multi_pred_cols(x)))
      })
    ]
    dt[,
      sens := map2(
        preds,
        subtype,
        ~ sens(.x, classification, ".pred_class") |> pull(.estimate)
      ) |>
        unlist()
    ]
    dt[,
      spec := map2(
        preds,
        subtype,
        ~ spec(.x, classification, ".pred_class") |> pull(.estimate)
      ) |>
        unlist()
    ]
    dt[,
      fmeas := map2(
        preds,
        subtype,
        ~ f_meas(.x, classification, ".pred_class") |> pull(.estimate)
      ) |>
        unlist()
    ]
  }

  dt[,
    boots := lapply(preds, \(x) {
      as.data.table(bootstraps(x, times = boot_times, strata = classification))
    })
  ]

  dt[,
    boots := map2(boots, subtype, \(x, subtype) {
      x[, let(
        roc = lapply(x$splits, function(y) {
          rs <- analysis(y)
          tryCatch(
            if (subtype != "multiclass") {
              as.data.table(roc_curve(
                rs,
                classification,
                paste0(".pred_", subtype)
              ))
            } else {
              cols <- setdiff(
                grep("^\\.pred_", names(rs), value = TRUE),
                ".pred_class"
              )
              as.data.table(roc_curve(rs, classification, all_of(cols)))
            },
            error = function(e) NULL
          )
        }),
        pr = map2_dbl(x$splits, subtype, function(y, subtype) {
          rs <- analysis(y)
          tryCatch(
            if (subtype != "multiclass") {
              pr_auc(rs, classification, paste0(".pred_", subtype))[[
                ".estimate"
              ]]
            } else {
              pr_auc(rs, classification, contains("IBS"))[[".estimate"]]
            },
            error = function(e) NA_real_
          )
        }),
        auc = map2_dbl(x$splits, subtype, function(y, subtype) {
          rs <- analysis(y)
          tryCatch(
            if (subtype != "multiclass") {
              roc_auc(rs, classification, paste0(".pred_", subtype))[[
                ".estimate"
              ]]
            } else {
              roc_auc(rs, classification, contains("IBS"))[[".estimate"]]
            },
            error = function(e) NA_real_
          )
        }),
        splits = NULL
      )]
    })
  ]

  dt[,
    boots := lapply(boots, \(x) {
      if (is.null(x)) {
        return(NULL)
      }
      x <- x[vapply(x$roc, is.data.table, logical(1))]
      if (nrow(x) == 0) {
        return(NULL)
      }
      x[, by = .(id, auc, pr), rbindlist(roc)]
    })
  ]

  dt[,
    boots_mean_auc := vapply(
      boots,
      \(x) mean(x$auc, na.rm = TRUE),
      FUN.VALUE = double(1)
    )
  ]
  dt[,
    boots_sd_auc := vapply(
      boots,
      \(x) sd(x$auc, na.rm = TRUE),
      FUN.VALUE = double(1)
    )
  ]
  dt[,
    boots_quant_auc := lapply(
      boots,
      \(x) quantile(x$auc, c(0.025, 0.975), na.rm = TRUE)
    )
  ]
  dt[, boots_lower := boots_quant_auc |> lapply(as.vector) |> lapply(`[`, 1)]
  dt[, boots_higher := boots_quant_auc |> lapply(as.vector) |> lapply(`[`, 2)]

  dt[,
    boots_mean_pr := vapply(
      boots,
      \(x) mean(x$pr, na.rm = TRUE),
      FUN.VALUE = double(1)
    )
  ]
  dt[,
    boots_sd_pr := vapply(
      boots,
      \(x) sd(x$pr, na.rm = TRUE),
      FUN.VALUE = double(1)
    )
  ]
  dt[,
    boots_quant_pr := lapply(
      boots,
      \(x) quantile(x$pr, c(0.025, 0.975), na.rm = TRUE)
    )
  ]
  dt[, boots_pr_lower := boots_quant_pr |> lapply(as.vector) |> lapply(`[`, 1)]
  dt[, boots_pr_higher := boots_quant_pr |> lapply(as.vector) |> lapply(`[`, 2)]

  dt
}


# 3. Binary Model ROC Curves & Confusion Matrices (Figure 4A & 4B) -------------

mods <- bind_rows(best_mods, best_mods_ibs) |>
  evaluate_and_bootstrap_models(boot_times = 100)

## 3a. Figure 4A: Binary ROC Curves --------------------------------------------
preds_dt <- mods[, keyby = .(subtype, use, auc), rbindlist(roc)]
boots_dt <- mods[
  use == "both",
  by = .(subtype, use, boots_mean_auc, boots_sd_auc),
  rbindlist(boots)
]

klein_pal <- paletteer::paletteer_d("MoMAColors::Klein")
subtype_order <- c("IBS", "IBSC", "IBSD", "IBSM")
subtype_cols <- setNames(klein_pal[c(5, 1, 2, 3)], subtype_order)

boots_dt[,
  col := fcase(
    subtype == "IBSC" , "#ca8290" ,
    subtype == "IBSD" , "#6e8a8d" ,
    subtype == "IBSM" , "#a27850" ,
    subtype == "IBS"  , "#798958" ,
    default = "grey50"
  )
]

subtype_order <- subtype_order[
  subtype_order %in% unique(c(preds_dt$subtype, boots_dt$subtype))
]
preds_dt[, subtype := factor(subtype, levels = subtype_order)]
boots_dt[, subtype := factor(subtype, levels = subtype_order)]

auc_dt <- preds_dt[use == "both", by = subtype, .(auc = unique(auc))]

roc_plot <- ggplot(mapping = aes(x = 1 - specificity, y = sensitivity)) +
  facet_wrap(~subtype, ncol = 2) +
  geom_path(
    data = boots_dt,
    colour = boots_dt$col,
    alpha = .1,
    lwd = .5
  ) +
  scale_color_manual(values = subtype_cols, guide = "none") +
  geom_path(
    data = preds_dt[use != "both"],
    aes(linetype = use),
    colour = "grey30",
    lwd = 1,
    alpha = .85
  ) +
  geom_path(
    data = preds_dt[use == "both"],
    aes(linetype = use, colour = subtype),
    lwd = 1,
    alpha = .85
  ) +
  scale_linetype_manual(
    name = "Model type",
    values = c(
      both = "solid",
      covars = "dotted",
      features = "dashed"
    ),
    labels = c("Full", "Covariates", "Peptides")
  ) +
  geom_abline(intercept = 0, slope = 1, col = "grey30") +
  scale_x_continuous(
    expand = expansion(add = 0),
    name = expression(1 - Specificity)
  ) +
  scale_y_continuous(
    expand = expansion(add = 0),
    name = expression(Sensitivity)
  ) +
  coord_equal(clip = "off") +
  geom_text(
    data = auc_dt,
    aes(
      x = .75,
      y = 0,
      label = paste0("AUC = ", round(auc, 2))
    ),
    colour = "grey30",
    size = 3,
    vjust = -.5
  ) +
  theme_my(
    panel.grid.major.y = element_blank(),
    legend.title.position = "top",
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.box = "vertical",
    legend.location = "plot"
  ) +
  tight

## 3b. Figure 4B: Binary Confusion Matrices ------------------------------------
plot_conf_panel <- function(dt, st, col, show_x = TRUE, show_y = TRUE) {
  dt <- copy(dt)
  dt[,
    cell_type := fcase(
      Truth == Prediction , paste0("True ", Truth) ,
      default = paste0("False ", Prediction)
    )
  ]
  dt[, let(
    pct = round(100 * N / sum(N), 1),
    text_col = fifelse(N >= 0.55 * max(N), "light", "dark")
  )]
  dt[, label := paste0(cell_type, "\nn=", N, "\n", pct, "%")]
  dt[, Truth := factor(Truth, levels = c("IBS", "noIBS"))]
  dt[, Prediction := factor(Prediction, levels = c("noIBS", "IBS"))]

  ggplot(dt, aes(x = Prediction, y = Truth, fill = N)) +
    geom_tile(color = NA) +
    geom_text(
      aes(label = label, color = text_col),
      fontface = "bold",
      lineheight = .85,
      size = 3.2
    ) +
    scale_fill_gradient(
      name = "Number of cases\nin test set",
      low = "#f7f7f7",
      high = unname(col),
      limits = c(0, NA)
    ) +
    scale_color_manual(
      values = c(light = "white", dark = "grey15"),
      guide = "none"
    ) +
    scale_x_discrete(
      labels = c(noIBS = "Predicted noIBS", IBS = "Predicted IBS")
    ) +
    scale_y_discrete(labels = c(noIBS = "Actual noIBS", IBS = "Actual IBS")) +
    labs(
      title = paste0("Confusion Matrix: ", st, " vs. noIBS (test set)"),
      x = NULL,
      y = NULL
    ) +
    theme_my(
      legend.position = "right",
      legend.title = element_text(angle = 90),
      panel.grid.major.y = element_blank()
    ) +
    theme(
      axis.text.x = if (show_x) {
        element_text(angle = 45, hjust = 1, vjust = 1)
      } else {
        element_blank()
      },
      axis.text.y = if (show_y) {
        element_text(angle = 45, hjust = 1)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
    ) +
    coord_fixed() +
    tight
}

conf_dt <- best_mods |>
  bind_rows(best_mods_ibs) |>
  filter(use == "both") |>
  select(subtype, conf_mat)
setDT(conf_dt)

conf_dt[, conf_mat := lapply(conf_mat, \(x) as.data.table(x$table))]
conf_dt <- conf_dt[, by = subtype, rbindlist(conf_mat)]
conf_dt[, let(
  Prediction = sub("C$|D$|M$", "", Prediction),
  Truth = sub("C$|D$|M$", "", Truth)
)]
conf_dt[, names(.SD) := lapply(.SD, as.factor), .SDcols = 2:3]

conf_subtypes <- subtype_order[subtype_order %in% unique(conf_dt$subtype)]
conf_ncol <- 2L
conf_nrow <- ceiling(length(conf_subtypes) / conf_ncol)

conf_panels <- lapply(seq_along(conf_subtypes), \(i) {
  st <- conf_subtypes[i]
  col <- subtype_cols[[st]]
  if (is.null(col) || is.na(col)) {
    col <- "grey50"
  }
  row_i <- ceiling(i / conf_ncol)
  col_i <- (i - 1L) %% conf_ncol + 1L
  plot_conf_panel(
    conf_dt[subtype == st],
    st,
    col,
    show_x = row_i == conf_nrow,
    show_y = col_i == 1L
  )
})

conf_mat <- wrap_plots(conf_panels, ncol = conf_ncol)


# 4. SHAP Feature Attribution Analysis (Figure 4C & Extended Table 7) ----------

shap_transform <- function(object) {
  imp <- shapviz:::.get_imp(get_shap_values(object), sort_features = TRUE)

  if (ncol(object) > 13L) {
    imp <- imp[seq_len(13L)]
  }
  ord <- names(imp)
  object <- object[, ord]

  S <- get_shap_values(object)
  X <- shapviz:::.scale_X(get_feature_values(object))
  transform(
    as.data.frame.table(S, responseName = "value"),
    feature = factor(Var2, levels = rev(ord)),
    color = as.data.frame.table(X)$Freq
  )
}

shap_raw <- best_mods |> filter(use == "both") |> select(subtype, shap_vals)
setDT(shap_raw)

shap <- shap_raw |> copy()
shap[, shap_vals := lapply(shap_vals, shap_transform)]
shap <- shap[, by = subtype, rbindlist(shap_vals)]

bee_width <- 0.4
bee_adjust <- 0.5
color_bar_title <- "Feature value"
viridis_args <- list(option = "E", begin = 0, end = .85)

abs_sum <- function(x) {
  sum(abs(x))
}

lut_shap <- shap[
  !feature %like% "^Sex|Age|Stool" & subtype == "IBSC",
  .(
    unique(feature),
    short = c(
      "LysM (Streptococcus)",
      "Polyprotein (Rhinovirus)",
      "EBNA1 (EBV) 1",
      "VP1 (Rhinovirus)",
      "EBNA1 (EBV) 2",
      "IgA protease (N. meningitidis)",
      "LukL2 (Staphylococcus)",
      "HAP (H. influenzae)",
      "FnIII (Bacteroidales)",
      "Polyprotein (Enterovirus B)",
      "Unknown (Bacteroidales)",
      "Glycoprotein C (HHV-3)"
    )
  )
]

shap <- shap[subtype == "IBSC"]
shap[lut_shap, on = .(feature = V1), short := short]
shap[is.na(short), short := feature]
shap[, short := reorder_within(short, -value, subtype, fun = abs_sum)]
shap[, y_within := reorder_within(feature, -value, subtype, fun = abs_sum)]

rng <- as.data.table(shap)[,
  .(
    xmin = min(value, na.rm = TRUE),
    xmax = max(value, na.rm = TRUE),
    y_top = uniqueN(short) + 0.8
  ),
  by = subtype
][, let(x_left = pmin(0, xmin) * 0.03, x_right = pmax(0, -xmin) * 0.03)][]

shap_anno <- droplevels(unique(shap[
  feature %like% "^agi|^twi",
  .(subtype, feature)
]))
shap_anno[, by = subtype, rank := seq_len(.N)]
shap_anno <- shap_anno[anno, on = .(feature = peptide), nomatch = NULL]
setorder(shap_anno, subtype, rank)

shap_extended <- shap_anno |> copy()
auc_pass <- c("IBSC")

s <- shap[subtype %in% auc_pass] |>
  ggplot2::ggplot(ggplot2::aes(x = value, y = short)) +
  scale_y_reordered() +
  facet_wrap(~subtype, scales = "free", ncol = 1) +
  ggplot2::geom_vline(xintercept = 0, color = "darkgray") +
  ggplot2::geom_point(
    ggplot2::aes(color = color),
    stroke = 0,
    alpha = .7,
    size = .7,
    position = shapviz:::position_bee(width = bee_width, adjust = bee_adjust)
  ) +
  shapviz:::.get_color_scale(
    viridis_args = viridis_args,
    bar = !is.null(color_bar_title),
    ncol = length(unique(shap$color))
  ) +
  ggplot2::labs(
    x = "SHAP value",
    y = ggplot2::element_blank(),
    color = color_bar_title
  ) +
  ggplot2::theme(legend.box.spacing = grid::unit(0, "pt")) +
  theme_my(panel.grid.major.y = element_blank()) +
  coord_cartesian(clip = "off") +
  theme(
    legend.title.position = "left",
    legend.title = element_text(angle = 90)
  ) +
  guides(
    color = guide_colorbar(
      direction = "vertical",
      barwidth = unit(2, "mm"),
      barheight = unit(4, "cm")
    )
  ) +
  tight +
  geom_text(
    data = rng[subtype %in% auc_pass],
    inherit.aes = FALSE,
    aes(x = x_left, y = y_top, label = "\u2190 Towards noIBS"),
    hjust = 1,
    vjust = -0.3,
    size = 2.8,
    fontface = "bold",
    colour = "grey30"
  ) +
  geom_text(
    data = rng[subtype %in% auc_pass],
    inherit.aes = FALSE,
    aes(
      x = x_right,
      y = y_top,
      label = paste0("Towards ", subtype, " \u2192")
    ),
    hjust = 0,
    vjust = -0.3,
    size = 2.8,
    fontface = "bold",
    colour = "grey30"
  ) +
  theme(strip.text = element_blank())

# Save Figure 4 panels
ggsave(
  "figure_4a.pdf",
  plot = roc_plot,
  width = 5.416667,
  height = 7.5
)
ggsave(
  "figure_4b.pdf",
  plot = conf_mat,
  width = 8,
  height = 7.5
)
ggsave(
  "figure_4c.pdf",
  plot = s,
  width = 6,
  height = 3.5
)

## 4a. Export Extended Data Table 7 (Top IBS-C SHAP Features) ------------------
prev_dt <- calc_prev(phenos, group_var = "IBS_Subtype", wide = FALSE)
prev_dt[, IBS_Subtype := gsub("\\-| ", "", IBS_Subtype)]

ctrl <- prev_dt[, dcast(
  .SD,
  subset = .(IBS_Subtype == "noIBS"),
  peptide ~ IBS_Subtype,
  value.var = "V1"
)]
prev_dt <- prev_dt[ctrl, on = .(peptide)][IBS_Subtype != "noIBS"]
shap_extended[
  prev_dt,
  on = .(subtype = IBS_Subtype, feature = peptide),
  let(prev = V1, prev_ctrl = noIBS)
]
shap_extended[,
  Protein := c(
    "LysM domain-containing protein from Streptococcus",
    "Genome polyprotein from Rhinovirus",
    "Nuclear antigen 1 from EBV",
    "VP1 capsid protein from Rhinovirus",
    "Nuclear antigen 1 from EBV",
    "IgA endopeptidase from N. meningitidis",
    "Leukocidin-like protein 2 from Staphylococcus",
    "Adhesion and penetration protein from H. influenzae",
    "Fibronectin type III domain protein from Bacteroidales",
    "Genome polyprotein from Enterovirus B",
    "Uncharacterized protein from Bacteroidales",
    "Envelope glycoprotein C from HHV-3"
  )
]

extended_tab_5 <- shap_extended[, .(
  `IBS Type` = subtype,
  peptide = feature,
  Protein,
  `Subtype Frequency (%)` = prev,
  `Control Frequency (%)` = prev_ctrl
)]
extended_tab_5[, .SDcols = is.numeric, names(.SD) := round(.SD * 100, 1)]
extended_tab_5 <- extended_tab_5[`IBS Type` %in% auc_pass]

write_xlsx(
  extended_tab_5,
  path = "extended_table_7.xlsx"
)


# 5. Export Extended Data Table 1 (Tuned Hyperparameters & Preprocessing) ------

tab_1 <- data.table(
  Model = "Search Space (Range)",
  trees = "300-3500",
  tree_depth = "2-8",
  learn_rate = "0.001-0.1",
  min_n = "2-80",
  loss_reduction = "0 – 10",
  mtry = "1 - p",
  Preprocessing = NA
)

setDT(best_mods)
setDT(best_mods_ibs)
setDT(best_mods_multi)

extract_tuned_summary <- function(mods_dt, model_prefix = NULL) {
  tuned <- mods_dt[,
    by = .(subtype, use),
    as.data.table(lapply(.SD$model[[1]]$fit$actions$model$spec$args, eval_tidy))
  ]

  prefix <- if (is.null(model_prefix)) tuned$subtype else model_prefix
  tuned[,
    Model := paste(
      prefix,
      fcase(
        use == "both"     , "(Combined)"   ,
        use == "covars"   , "(Covariates)" ,
        use == "features" , "(Peptides)"
      )
    )
  ]

  prep <- mods_dt[,
    by = .(subtype, use),
    .(
      Preprocessing = .SD$model[[1]] |>
        workflows::extract_recipe() |>
        tidy() |>
        pull(type)
    )
  ]
  prep[,
    Preprocessing := fcase(
      Preprocessing == "zv"          , "ZV Filter"                ,
      Preprocessing == "impute_knn"  , "KNN Impute (Peptides)"    ,
      Preprocessing == "impute_mean" , "Mean Impute (Covariates)" ,
      Preprocessing == "corr"        , "Correlation Filter"       ,
      Preprocessing == "normalize"   , "Normalize"                ,
      Preprocessing == "dummy"       , "Dummy Vars"               ,
      Preprocessing == "center"      , "Center"                   ,
      Preprocessing == "scale"       , "Scale"
    )
  ]
  prep <- prep |> na.omit()
  prep <- prep[,
    by = .(subtype, use),
    .(Preprocessing = do.call(paste, list(Preprocessing, collapse = " + ")))
  ]
  prep$subtype <- prep$use <- NULL

  cbind(tuned, prep)
}

tuned       <- extract_tuned_summary(best_mods)
tuned_ibs   <- extract_tuned_summary(best_mods_ibs, model_prefix = "Any IBS")
tuned_multi <- extract_tuned_summary(best_mods_multi, model_prefix = "Multiclass")

tuned_all <- rbind(tuned, tuned_multi, tuned_ibs)
tuned_all$subtype <- tuned_all$use <- NULL
tuned_all[, let(
  learn_rate = formatC(learn_rate, digits = 2, format = "e"),
  loss_reduction = formatC(loss_reduction, format = "e", digits = 2)
)]
setcolorder(tuned_all, "Model", before = 1)

tab1 <- rbind(tab_1, tuned_all)
setnames(
  tab1,
  c(2, 3, 5, 7),
  c(
    "Number of Trees (trees)",
    "Max Tree Depth (tree_depth)",
    "Minimum Node Size (min_n)",
    "Predictors Sampled (mtry)"
  )
)

write_xlsx(
  tab1,
  path = "extended_table_1.xlsx"
)


# 6. Multiclass Evaluation (Extended Data Figure 8 & Extended Data Table 6) ----

mods_ibs <- evaluate_and_bootstrap_models(best_mods_ibs, boot_times = 100)

multi_metrics_boot <- function(x, times = 1000, seed = 1337) {
  set.seed(seed)
  classes <- levels(x$classification)

  calc_all_metrics <- function(df) {
    map_dfr(classes, function(cls) {
      binary_truth <- factor(
        ifelse(df$classification == cls, "Event", "Other"),
        levels = c("Event", "Other")
      )
      binary_pred <- factor(
        ifelse(df$.pred_class == cls, "Event", "Other"),
        levels = c("Event", "Other")
      )
      prob_vals <- df[[paste0(".pred_", cls)]]

      tibble(
        Subtype = cls,
        roc_auc = roc_auc_vec(truth = binary_truth, estimate = prob_vals),
        pr_auc  = pr_auc_vec(truth = binary_truth, estimate = prob_vals),
        sens    = sens_vec(truth = binary_truth, estimate = binary_pred),
        spec    = spec_vec(truth = binary_truth, estimate = binary_pred),
        f_meas  = f_meas_vec(truth = binary_truth, estimate = binary_pred)
      )
    })
  }

  calc_aucs_only <- function(df) {
    map_dfr(classes, function(cls) {
      binary_truth <- factor(
        ifelse(df$classification == cls, "Event", "Other"),
        levels = c("Event", "Other")
      )
      prob_vals <- df[[paste0(".pred_", cls)]]

      tibble(
        Subtype = cls,
        roc_auc = roc_auc_vec(truth = binary_truth, estimate = prob_vals),
        pr_auc  = pr_auc_vec(truth = binary_truth, estimate = prob_vals)
      )
    })
  }

  point_estimates <- calc_all_metrics(x)
  boots <- bootstraps(x, times = times)
  boot_results <- map_dfr(boots$splits, ~ calc_aucs_only(analysis(.x)))

  ci_results <- boot_results |>
    group_by(Subtype) |>
    summarize(
      roc_auc_lower = quantile(roc_auc, probs = 0.025, na.rm = TRUE),
      roc_auc_upper = quantile(roc_auc, probs = 0.975, na.rm = TRUE),
      pr_auc_lower  = quantile(pr_auc, probs = 0.025, na.rm = TRUE),
      pr_auc_upper  = quantile(pr_auc, probs = 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  point_estimates |>
    left_join(ci_results, by = "Subtype") |>
    select(
      Subtype,
      roc_auc, roc_auc_lower, roc_auc_upper,
      pr_auc, pr_auc_lower, pr_auc_upper,
      sens, spec, f_meas
    )
}

mods_mutli_individual <- map(best_mods_multi$preds, multi_metrics_boot) |>
  setNames(c("both", "covars", "features")) |>
  rbindlist(idcol = "use")

mods_multi <- evaluate_and_bootstrap_models(best_mods_multi, boot_times = 100)

## 6a. Extended Data Figure 8: Multiclass One-vs-All ROC & Confusion Matrix ----
multi_class_order <- c("IBSC", "IBSD", "IBSM", "noIBS")
multi_class_lab <- c(
  "IBSC"  = "IBS-C",
  "IBSD"  = "IBS-D",
  "IBSM"  = "IBS-M",
  "noIBS" = "noIBS"
)
multi_klein_cols <- setNames(klein_pal[c(1, 2, 3, 4)], multi_class_order)
multi_boot_cols <- c(
  IBSC  = "#ca8290",
  IBSD  = "#6e8a8d",
  IBSM  = "#a27850",
  noIBS = "#AE984BFF"
)

preds_multi_dt <- mods_multi[use == "both", keyby = use, rbindlist(roc)]
preds_multi_dt[, class := factor(.level, levels = multi_class_order)]

preds_multi_alt_dt <- mods_multi[use != "both", keyby = use, rbindlist(roc)]
preds_multi_alt_dt[, class := factor(.level, levels = multi_class_order)]

boots_multi_dt <- mods_multi[use == "both", rbindlist(boots)]
boots_multi_dt[, class := factor(.level, levels = multi_class_order)]
boots_multi_dt[, col := multi_boot_cols[as.character(class)]]
boots_multi_dt[is.na(col), col := "grey50"]

auc_multi_dt <- mods_multi[use == "both", {
  x <- preds[[1]]
  rbindlist(lapply(multi_class_order, \(cls) {
    binary_truth <- factor(
      fifelse(x$classification == cls, "Event", "Other"),
      levels = c("Event", "Other")
    )
    data.table(
      class = cls,
      auc = roc_auc_vec(
        truth = binary_truth,
        estimate = x[[paste0(".pred_", cls)]]
      )
    )
  }))
}]
auc_multi_dt[, class := factor(class, levels = multi_class_order)]

roc_multi_plot <- ggplot(mapping = aes(x = 1 - specificity, y = sensitivity)) +
  facet_wrap(~class, ncol = 2, labeller = as_labeller(multi_class_lab)) +
  geom_path(
    data = boots_multi_dt,
    aes(colour = col, group = interaction(id, class)),
    alpha = .1,
    linewidth = .5
  ) +
  geom_path(
    data = preds_multi_alt_dt,
    aes(linetype = use),
    colour = "grey30",
    linewidth = 1,
    alpha = .85
  ) +
  geom_path(
    data = preds_multi_dt,
    aes(colour = class, linetype = use),
    linewidth = 1,
    alpha = .85
  ) +
  scale_linetype_manual(
    name = "Model type",
    values = c(
      both = "solid",
      covars = "dotted",
      features = "dashed"
    ),
    labels = c("Full", "Covariates", "Peptides")
  ) +
  geom_abline(intercept = 0, slope = 1, colour = "grey30") +
  geom_text(
    data = auc_multi_dt,
    aes(
      x = .75,
      y = 0,
      label = paste0("AUC = ", round(auc, 2))
    ),
    inherit.aes = FALSE,
    colour = "grey30",
    size = 3,
    vjust = -.5
  ) +
  scale_color_manual(values = multi_klein_cols, guide = "none") +
  scale_x_continuous(
    expand = expansion(add = 0),
    name = expression(1 - Specificity)
  ) +
  scale_y_continuous(
    expand = expansion(add = 0),
    name = expression(Sensitivity)
  ) +
  coord_equal(clip = "off") +
  theme_my(
    panel.grid.major.y = element_blank(),
    legend.title.position = "top",
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.box = "vertical",
    legend.location = "plot"
  ) +
  tight

conf_multi_dt <- best_mods_multi[use == "both", conf_mat][[1]]$table |>
  as.data.table()
if (!"N" %in% names(conf_multi_dt)) {
  setnames(
    conf_multi_dt,
    intersect(names(conf_multi_dt), c("Freq", "n"))[1],
    "N"
  )
}
conf_multi_dt[, `:=`(
  Truth = factor(Truth, levels = multi_class_order),
  Prediction = factor(Prediction, levels = multi_class_order),
  correct = Prediction == Truth
)]

conf_multi_plot <- ggplot(
  conf_multi_dt,
  aes(x = Truth, y = Prediction, fill = correct)
) +
  geom_tile(colour = NA) +
  geom_text(aes(label = N), fontface = "bold", size = 5) +
  scale_x_discrete(labels = multi_class_lab) +
  scale_y_discrete(labels = multi_class_lab) +
  scale_fill_manual(
    values = c(`TRUE` = "#7CB8A6", `FALSE` = "#E8A4B4"),
    guide = "none"
  ) +
  labs(title = "multiclass", x = "Truth", y = "Prediction") +
  theme_my(panel.grid.major.y = element_blank(), legend.position = "none") +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  coord_fixed() +
  tight

multi_figure <- roc_multi_plot +
  conf_multi_plot +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "#fbf9f4", colour = NA)
  )

ggsave(
  "extended_figure_8.pdf",
  plot = multi_figure,
  width = 10.5
)

## 6b. Export Extended Data Table 6 (Model Performance Summary) ----------------
format_perf_cols <- function(dt, model_prefix = NULL) {
  prefix <- if (is.null(model_prefix)) dt$subtype else model_prefix
  dt[,
    Model := paste(
      prefix,
      fcase(
        use == "both"     , "(Combined)"   ,
        use == "covars"   , "(Covariates)" ,
        use == "features" , "(Peptides)"
      )
    )
  ]
  dt[, .(
    Model,
    auc,
    boots_lower,
    boots_higher,
    pr_auc,
    boots_pr_lower,
    boots_pr_higher,
    sens,
    spec,
    fmeas
  )]
}

mods_summary       <- format_perf_cols(mods)
mods_ibs_summary   <- format_perf_cols(mods_ibs, model_prefix = "Any IBS")
mods_multi_summary <- format_perf_cols(mods_multi, model_prefix = "Multiclass Macro")

mods_mutli_individual[Subtype == "noIBS", Subtype := "Ctrl"]
mods_mutli_individual[,
  Model := paste(
    "Mutliclass",
    Subtype,
    "specific",
    fcase(
      use == "both"     , "(Combined)"   ,
      use == "covars"   , "(Covariates)" ,
      use == "features" , "(Peptides)"
    )
  )
]
mods_mutli_individual <- mods_mutli_individual[, .(
  Model,
  auc = roc_auc,
  boots_lower = roc_auc_lower,
  boots_higher = roc_auc_upper,
  pr_auc,
  boots_pr_lower = pr_auc_lower,
  boots_pr_higher = pr_auc_upper,
  sens,
  spec,
  fmeas = f_meas
)]

extended_tab_4 <- rbind(
  mods_summary,
  mods_multi_summary,
  mods_mutli_individual,
  mods_ibs_summary
)
extended_tab_4[,
  names(.SD) := lapply(.SD, unlist),
  .SDcols = c(
    "boots_lower",
    "boots_higher",
    "boots_pr_lower",
    "boots_pr_higher"
  )
]
extended_tab_4[,
  names(.SD) := lapply(.SD, round, digits = 2),
  .SDcols = is.numeric
]
extended_tab_4[,
  `ROC-AUC (95% CI)` := paste0(auc, " (", boots_lower, "-", boots_higher, ")")
]
extended_tab_4[,
  `PR-AUC (95% CI)` := paste0(
    pr_auc,
    " (",
    boots_pr_lower,
    "-",
    boots_pr_higher,
    ")"
  )
]
extended_tab_4 <- extended_tab_4[, .(
  Model,
  `ROC-AUC (95% CI)`,
  `PR-AUC (95% CI)`,
  Sensitivity = sens,
  Specificity = spec,
  `F-Measure` = fmeas
)]

write_xlsx(
  extended_tab_4,
  path = "extended_table_6.xlsx"
)
