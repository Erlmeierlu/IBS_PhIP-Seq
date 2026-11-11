# ==============================================================================
# 2_Exploratory.R
# Baseline characteristics, confounder screening, PCA, and clustering
# ==============================================================================
#
# Purpose:
#   1. Reads the food-antigen table (`food_antigen_annotation.csv`) and exports
#      `food_categories.csv` for peptides present in the baseline data.
#   2. Screens demographic, technical, and clinical covariates (Age, Sex,
#      plate_id, CMV serostatus, Stool_Freq) for confounding across IBS status
#      and subtypes at baseline and follow-up.
#   3. Generates the baseline cohort characteristics summary table
#      (`baseline_characteristics.xlsx`).
#   4. Filters peptides with >5% missingness, performs mode imputation for
#      ordination, runs Principal Component Analysis (PCA), and assembles
#      Manuscript Figure 1 (`figure_1.pdf`).
#   5. Evaluates PC associations with IBS subtypes, performs k-means clustering
#      (k = 2) on PCA coordinates, and tests cluster drivers.
#   6. Tests PC1 and PC2 associations across IBS subtypes, CMV serostatus,
#      antibody binding load, and sequencing plate.
#
# Expected Inputs:
#   - `1_phenos_ibs.rds`: QC'd baseline PhIP-Seq dataset (from 1_QC.R),
#     including `stool_type`.
#   - `1_phenos_fup.rds`: QC'd follow-up PhIP-Seq dataset (from 1_QC.R),
#     including `followup_time`.
#   - `fig1_qc_panels.rds`: Figure 1 panels saved by 1_QC.R.
#   - `food_antigen_annotation.csv`: Food antigens used in the analysis
#     (`antigen`, `food`, `food subcategory`, `food item`, `description`).
#
# Generated Outputs:
#   - `food_categories.csv`: Food allergen category lookup table.
#   - `baseline_characteristics.xlsx`: Formatted baseline cohort table.
#   - `2_phenos_subtypes.rds`: Baseline dataset with cluster assignments
#     and missingness-filtered peptides (<5% NA).
#   - `figure_1.pdf`: Manuscript Figure 1.
#   - Confounder, Age/Sex distribution, PCA p-value histogram, and clustering
#     driver plots.
# ==============================================================================

library(data.table)
library(readr)
library(ggplot2)
library(broom)
library(purrr)
library(future.apply)
library(factoextra)
library(gt)
library(dunn.test)
library(paletteer)
library(patchwork)
library(gtsummary)

source("0_utils.R")

# 1. Load Data & Prepare Food Allergen Annotations -----------------------------

phenos <- read_rds("1_phenos_ibs.rds")
phenos_fup <- read_rds("1_phenos_fup.rds")

phenos_fup[, median(followup_time)]

# Food allergen category lookup
food_anno <- fread("food_antigen_annotation.csv")

peps <- names(phenos[, .SD, .SDcols = patterns("^agi|^twi")])

food <- food_anno[
  antigen %in% peps & `food subcategory` != "",
  .(
    peptide = antigen,
    category = `food subcategory`,
    item = `food item`
  )
]
fwrite(food, "food_categories.csv")


# 2. Confounder Screening & Baseline Cohort Table ------------------------------

## 2a. General IBS vs no IBS Confounder Screening ------------------------------
phenos[, shapiro.test(Age)]
phenos[, wilcox.test(Age ~ IBS_Status)]
phenos[, chisq.test(Sex, IBS_Status)]
phenos[, fisher.test(plate_id, IBS_Status, simulate.p.value = TRUE)]
phenos[, chisq.test(CMV_prediction, IBS_Status)]
phenos[, shapiro.test(Stool_Freq)]
phenos[, wilcox.test(Stool_Freq ~ IBS_Status)]

phenos[, table(Sex, IBS_Status)]
phenos[, chisq.test(Sex, IBS_Status)]$expected

## 2b. Baseline Cohort Characteristics Summary Table ---------------------------
char <- phenos[, .SDcols = !patterns("^agi|^twi"), .SD]
char[, fup := fifelse(ID %in% phenos_fup$ID, 1, 0)]
char[phenos_fup, on = .(ID), time_interval := followup_time]
char[, Age := as.integer(Age)]

tbl1_gs <- char |>
  tbl_summary(
    by = IBS_Subtype,
    include = c(
      Age,
      Sex,
      Stool_Freq,
      stool_type,
      n_distinct_peptides,
      fup,
      time_interval
    ),
    value = list(Sex ~ "F"),
    digits = all_continuous2() ~ 2,
    missing = "no",
    label = list(
      Stool_Freq = "Bowel movements per day",
      stool_type = "Stool consistency (BSFS score)",
      fup = "Follow-up samples",
      n_distinct_peptides = "Unique antigens bound",
      time_interval = "Time interval in years"
    )
  ) |>
  add_p(include = -c(fup)) |>
  modify_header(
    label ~ "**Baseline Characteristic**",
    stat_1 ~ "**HC**<br>N = {n}"
  ) |>
  add_stat_label(label = list(Sex ~ "Female (%)")) |>
  modify_spanning_header(c("stat_2", "stat_3", "stat_4") ~ "**IBS**") |>
  bold_p()

tbl1_n <- tbl1_gs |>
  modify_table_body(
    ~ .x %>%
      mutate(
        label = ifelse(row_type == "label", paste0("**", label, "**"), label)
      )
  ) |>
  as_gt() |>
  gt::fmt_markdown(columns = label)

tbl2_n <- tbl1_n |>
  tab_options(
    table.background.color = "#fbf9f4",
    table.border.top.color = "#669999",
    table.border.top.width = px(2),
    table_body.border.top.color = "#669999",
    table_body.border.top.width = px(2),
    table_body.border.bottom.color = "#669999",
    table_body.border.bottom.width = px(2),
    column_labels.border.top.color = "#669999",
    column_labels.border.top.width = px(2),
    column_labels.border.bottom.color = "#669999",
    column_labels.border.bottom.width = px(2)
  )

tbl2_xlsx <- tbl1_gs |>
  bold_labels() |>
  as_flex_table() |>
  flextable::bg(bg = "#fbf9f4", part = "all") |>
  flextable::border_outer(
    border = officer::fp_border(color = "#669999", width = 2)
  ) |>
  flextable::border_inner_h(
    border = officer::fp_border(color = "#669999", width = 1)
  ) |>
  flextable::border_inner_v(
    border = officer::fp_border(color = "#669999", width = 1)
  ) |>
  flextable::bold(part = "header", bold = TRUE) |>
  flextable::bold(
    i = which(tbl1_gs[["table_body"]][["row_type"]] == "label"),
    j = 1,
    bold = TRUE,
    part = "body"
  )

openxlsx2::wb_workbook() |>
  openxlsx2::wb_add_worksheet("Baseline characteristics") |>
  flexlsx::wb_add_flextable("Baseline characteristics", tbl2_xlsx) |>
  openxlsx2::wb_save("baseline_characteristics.xlsx")

## 2c. IBS Subtype Confounder Screening ----------------------------------------
phenos[, table(Sex, IBS_Subtype)]
phenos[, kruskal.test(Age ~ IBS_Subtype)]
phenos[, dunn.test(Age, IBS_Subtype, method = "bh")]
phenos[, chisq.test(Sex, IBS_Subtype)]
phenos[, fisher.test(plate_id, IBS_Subtype, simulate.p.value = TRUE)]
phenos[, fisher.test(CMV_prediction, IBS_Subtype, simulate.p.value = TRUE)]
phenos[, kruskal.test(Stool_Freq ~ IBS_Subtype)]
phenos[, dunn.test(Stool_Freq, IBS_Subtype, method = "bh")]
phenos[, pairwise.wilcox.test(Stool_Freq, IBS_Subtype)]

pvals <- c(
  phenos[, kruskal.test(Age ~ IBS_Subtype)$p.value],
  phenos[, chisq.test(Sex, IBS_Subtype)$p.value],
  phenos[, fisher.test(plate_id, IBS_Subtype, simulate.p.value = TRUE)$p.value],
  phenos[,
    fisher.test(CMV_prediction, IBS_Subtype, simulate.p.value = TRUE)$p.value
  ],
  phenos[, kruskal.test(Stool_Freq ~ IBS_Subtype)]$p.value
)

confounders <- data.table(
  var = c("Age", "Sex", "plate_id", "CMV_prediction", "Stool_Freq"),
  p_value = pvals
)
confounders |>
  ggplot(aes(reorder(var, -p_value), -log10(p_value), fill = p_value < 0.05)) +
  geom_col(col = "grey30") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "deepskyblue2"
  ) +
  coord_flip() +
  labs(
    title = "Confounder Screening: IBS Subtype",
    x = "Variable",
    y = expression(-log[10](pval))
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_my()

ggsave("confounders.pdf")

## 2d. Follow-Up Cohort Confounder Screening -----------------------------------
phenos_fup[, by = IBS_Subtype, shapiro.test(Age)$p.value]
phenos_fup[, kruskal.test(Age ~ IBS_Subtype)]
phenos_fup[, dunn.test(Age, IBS_Subtype, method = "bh")]
phenos_fup[, table(Sex, IBS_Subtype)]
phenos_fup[, fisher.test(Sex, IBS_Subtype)]
phenos_fup[, table(plate_id, IBS_Subtype)]
phenos_fup[, fisher.test(plate_id, IBS_Subtype, simulate.p.value = TRUE)]
phenos_fup[, table(CMV_prediction, IBS_Subtype)]
phenos_fup[, fisher.test(CMV_prediction, IBS_Subtype)]

## 2e. Age and Sex Distribution Plots ------------------------------------------
phenos |>
  ggplot() +
  geom_histogram(aes(Age, fill = Sex), col = "grey30", alpha = 0.9) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(expand = expand_scale(add = c(0, 3))) +
  scale_y_continuous(expand = expand_scale(add = c(0, 3))) +
  theme_my() +
  facet_wrap(~IBS_Subtype, scales = "free_y") +
  xlab("Age (years)") +
  ylab("")

ggsave("Age_Sex_distribution_Subtypes.pdf")

phenos |>
  ggplot() +
  geom_histogram(aes(Age, fill = Sex), col = "grey30", alpha = 0.9) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(expand = expand_scale(add = c(0, 3))) +
  scale_y_continuous(expand = expand_scale(add = c(0, 3))) +
  theme_my() +
  facet_wrap(~IBS_Status) +
  xlab("Age (years)") +
  ylab("")

ggsave("Age_Sex_distribution.pdf")

phenos_fup |>
  ggplot() +
  geom_histogram(aes(Age, fill = Sex), col = "grey30", alpha = 0.9) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(expand = expand_scale(add = c(0, 3))) +
  scale_y_continuous(expand = expand_scale(add = c(0, 3))) +
  theme_my() +
  facet_wrap(~IBS_Subtype, scales = "free_y") +
  xlab("Age (years)") +
  ylab("")

ggsave("Age_Sex_distribution_follow_up.pdf")

phenos_fup |>
  ggplot() +
  geom_histogram(aes(Age, fill = Sex), col = "grey30", alpha = 0.9) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(expand = expand_scale(add = c(0, 3))) +
  scale_y_continuous(expand = expand_scale(add = c(0, 3))) +
  theme_my() +
  facet_wrap(~IBS_Status_fup) +
  xlab("Age (years)") +
  ylab("")

ggsave("Age_Sex_distribution_ibs_follow_up.pdf")


# 3. Principal Component Analysis (PCA) & Figure 1 Assembly --------------------

plot_pca <- function(
  pca_object,
  data,
  colour_by = NULL,
  plot_ellipse = TRUE,
  pcs = c("PC1", "PC2")
) {
  colour_by <- match.arg(colour_by, colnames(data))
  pca_object <- summary(pca_object)
  pca_data <- pca_object$x |> as.data.table(keep.rownames = "sample_id")
  pca_data <- data[pca_data, on = .(sample_id)]
  pca_data <- pca_data[, .SD, .SDcols = !patterns("^agi|^twi")]

  x_data <- sort(pcs)[1]
  y_data <- sort(pcs)[2]

  x_label <- paste0(
    x_data,
    " (",
    round(pca_object$importance[2, x_data] * 100, 2),
    "%)"
  )
  y_label <- paste0(
    y_data,
    " (",
    round(pca_object$importance[2, y_data] * 100, 2),
    "%)"
  )

  fill_function <- if (is.null(colour_by)) {
    NULL
  } else if (scale_type(data[[colour_by]]) == "continuous") {
    scale_fill_fermenter(palette = "Set3")
  } else {
    scale_fill_paletteer_d(palette = "MoMAColors::Klein")
  }

  geom_function <- if (is.null(colour_by)) {
    geom_point(alpha = 0.5)
  } else {
    geom_point(aes(fill = !!ensym(colour_by)), shape = 21, alpha = 0.5)
  }

  ellipse_function <- if (isFALSE(plot_ellipse)) {
    NULL
  } else if (is.null(colour_by)) {
    NULL
  } else {
    stat_ellipse(
      aes(col = !!ensym(colour_by)),
      level = 0.95,
      linetype = "dashed",
      alpha = 1
    )
  }

  col_function <- if (is.null(colour_by)) {
    NULL
  } else if (scale_type(data[[colour_by]]) == "continuous") {
    scale_colour_fermenter(palette = "Set3")
  } else {
    scale_colour_paletteer_d(palette = "MoMAColors::Klein")
  }

  ggplot(pca_data, aes(x = !!ensym(x_data), y = !!ensym(y_data))) +
    geom_function +
    ellipse_function +
    fill_function +
    col_function +
    theme_my(title = element_blank()) +
    labs(x = x_label, y = y_label)
}

# Exclude peptides with > 5% missing values and mode-impute remaining NAs
peptides <- phenos[, .SD, .SDcols = patterns("^agi|^twi")]
(colSums(is.na(peptides)) / nrow(peptides)) |> round(digits = 2) |> table()

peptides <- peptides[,
  .SD,
  .SDcols = round(colSums(is.na(peptides)) / nrow(peptides), digits = 2) > 0.05
] |>
  names()

peptides <- setdiff(names(phenos), peptides)
phenos <- phenos[, ..peptides]

IBS_samples <- copy(phenos)
to_fill <- names(IBS_samples[, .SD, .SDcols = patterns("^agi|^twi")])

Mode <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}

plan(multisession)

IBS_samples[,
  names(.SD) := future_lapply(.SD, function(x) {
    x[is.na(x)] <- Mode(x)
    x
  }),
  .SDcols = to_fill
]

IBS_samples_peptides <- IBS_samples[,
  .SD,
  .SDcols = patterns("^agi|^twi|sample_id")
]
IBS_samples_peptides <- as.data.frame(IBS_samples_peptides)
row.names(IBS_samples_peptides) <- IBS_samples_peptides$sample_id
IBS_samples_peptides$sample_id <- NULL

# Run PCA (unscaled for binary presence/absence matrix)
pca_result <- prcomp(IBS_samples_peptides, scale. = FALSE)

# Number of PCs explaining > 90% cumulative variance
dimensionality <- summary(pca_result)$importance[,
  which(summary(pca_result)$importance[3, ] > 0.9)[1],
  drop = FALSE
]
print(dimensionality)

IBS_samples[,
  IBS_Subtype := factor(
    IBS_Subtype,
    levels = c("IBS-C", "IBS-D", "IBS-M", "no IBS")
  )
]
pca1 <- plot_pca(pca_result, data = IBS_samples, colour_by = "IBS_Subtype") +
  theme(panel.grid.major.y = element_blank()) +
  coord_cartesian(clip = "off") +
  guides(
    colour = guide_legend(title = "Study group"),
    fill = guide_legend(title = "Study group")
  )

IBS_samples[,
  CMV_prediction := fifelse(CMV_prediction == 1, "Positive", "Negative")
]
pca2 <- plot_pca(pca_result, data = IBS_samples, colour_by = "CMV") +
  theme(panel.grid.major.y = element_blank()) +
  guides(
    colour = guide_legend(title = "Serological\nCMV status"),
    fill = guide_legend(title = "Serological\nCMV status")
  ) +
  coord_cartesian(clip = "off")

# Load Figure 1 panels from 1_QC.R if not already in the environment
if (!all(c("oligo_info", "raincloud", "prev_filter_plot") %in% ls())) {
  fig1_panels <- read_rds("fig1_qc_panels.rds")
  oligo_info <- fig1_panels$oligo_info
  raincloud <- fig1_panels$raincloud
  prev_filter_plot <- fig1_panels$prev_filter_plot
}

free(
  oligo_info +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank()
    )
) /
  (raincloud +
    prev_filter_plot +
    ((pca1 / pca2) + plot_layout(axis_titles = "collect_x")) +
    plot_layout(widths = c(.6, .7, .7))) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "#fbf9f4", colour = NA)
  )

ggsave("figure_1.pdf", width = 15, height = 13)

# Test association of each PC with IBS subtype via ANOVA
Associations_PC <- as.data.table(pca_result$x) |>
  map_vec(
    \(x, y) {
      m1 <- lm(x ~ y)
      m0 <- lm(x ~ 1)
      P <- anova(m0, m1)
      P$`Pr(>F)`[2]
    },
    phenos$IBS_Subtype
  ) |>
  rbind() |>
  t() |>
  as.data.table(keep.rownames = "pc")

setnames(Associations_PC, "V1", "p")

Associations_PC |>
  ggplot() +
  geom_histogram(aes(p)) +
  theme_my()

ggsave("pca_association_pval_hist.pdf")

Associations_PC[, fdr := p.adjust(p, "fdr")]
Associations_PC[fdr < 0.05]
Associations_PC[order(p)]


# 4. K-Means Clustering on PCA Coordinates & Driver Analysis -------------------

fviz_nbclust(pca_result$x, kmeans, method = "wss")
fviz_nbclust(pca_result$x, kmeans, method = "silhouette")

k <- 2
clusters <- kmeans(pca_result$x, centers = k, nstart = 25)
phenos$cluster <- as.factor(clusters$cluster)

# Test whether clusters are driven by IBS status or covariates
phenos[,
  lapply(.SD, function(x) chisq.test(cluster, x)$p.value),
  .SDcols = c("CMV_prediction", "Sex", "IBS_Subtype", "IBS_Status", "plate_id")
]
phenos[, .(
  Age = kruskal.test(Age ~ cluster)$p.value,
  Stool_Freq = kruskal.test(Stool_Freq ~ cluster)$p.value
)]

drivers <- cbind(
  phenos[,
    lapply(.SD, function(x) chisq.test(cluster, x)$p.value),
    .SDcols = c(
      "CMV_prediction",
      "Sex",
      "IBS_Subtype",
      "IBS_Status",
      "plate_id"
    )
  ],
  phenos[, .(
    Age = kruskal.test(Age ~ cluster)$p.value,
    Stool_Freq = kruskal.test(Stool_Freq ~ cluster)$p.value
  )]
)[, melt(.SD)]

drivers |>
  ggplot(aes(reorder(variable, -value), -log10(value), fill = value < 0.05)) +
  geom_col(col = "grey30") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "deepskyblue2"
  ) +
  coord_flip() +
  labs(
    title = "Clustering Drivers",
    x = "Variable",
    y = expression(-log[10](pval))
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_my()

ggsave("clustering_drivers.pdf")


# 5. PC1 / PC2 Univariate & Multivariate Tests ---------------------------------

nnn <- phenos |> copy()
nnn[,
  `Enriched antibodies` := fcase(
    n_distinct_peptides <= 1000 , "1000" ,
    n_distinct_peptides <= 2000 , "2000" ,
    n_distinct_peptides <= 4000 , "3000"
  )
]

pca_12 <- pca_result$x[, 1:10]
pca_12 <- (pca_12 |> as.data.table(keep.rownames = "sample_id"))[
  nnn,
  on = .NATURAL
]

# IBS Subtype
pca_12[, kruskal.test(PC1, IBS_Subtype, method = "bh")]
pca_12[, kruskal.test(PC2, IBS_Subtype, method = "bh")]

# CMV Serostatus
pca_12[, wilcox.test(PC1 ~ CMV_prediction, method = "bh")]
pca_12[, wilcox.test(PC2 ~ CMV_prediction, method = "bh")]
manova_res <- manova(cbind(PC1, PC2) ~ CMV_prediction, data = pca_12)
summary(manova_res)

# Enriched Antibody Load & Sequencing Plate
pca_12[, kruskal.test(PC1, `Enriched antibodies`, method = "bh")]
pca_12[, kruskal.test(PC2, `Enriched antibodies`, method = "bh")]
pca_12[, kruskal.test(PC1, plate_id, method = "bh")]
pca_12[, kruskal.test(PC2, plate_id, method = "bh")]

# Summary statistics of bound peptide counts by IBS status and top 20 PC variance
nnn[, by = IBS_Status, median(as.integer(n_distinct_peptides))]
nnn[, by = IBS_Status, quantile(as.integer(n_distinct_peptides))[c(2, 4)]]
round(sum(summary(pca_result)$importance[2, 1:20]) * 100, 2)


# 6. Save Processed Dataset ----------------------------------------------------

write_rds(phenos, "2_phenos_subtypes.rds")
