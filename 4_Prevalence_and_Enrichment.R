# ==============================================================================
# 4_Prevalence_and_Enrichment.R
# Seroprevalence ratio testing, taxonomic/functional GSEA, and longitudinal stability
# ==============================================================================
#
# Purpose:
#   1. Exports Extended Data Tables 2, 3, 4, and 5.
#   2. Evaluates peptide seroprevalence ratios (Wilcoxon rank-sum test of each
#      taxon or source database vs all others) and taxonomic fast Gene Set
#      Enrichment Analysis (`fgsea`) for IBS subtypes vs controls and general
#      IBS vs controls, assembling Manuscript Figure 2 (`figure_2.pdf`).
#   3. Performs enrichment analyses for post-infectious IBS (PI-IBS) markers
#      (Extended Data Figure 5) and food allergen categories (Extended Data
#      Figure 6).
#   4. Quantifies longitudinal antibody repertoire stability using pairwise
#      Phi-coefficients between baseline and follow-up samples, runs a node-level
#      permutation test across subtypes, compares within- vs between-individual
#      similarity and bound peptide counts (Extended Data Figure 7), and builds
#      Manuscript Figure 3 (`figure_3.pdf`) combining subtype
#      trajectories (alluvial plot), the Phi heatmap, and prevalence/database
#      shifts in Persistent/switched and Improved cohorts.
#   5. Runs taxonomic enrichment (`fgsea`) for Age and Sex associations
#      (Extended Data Figure 3).
#   6. Runs sensitivity analyses excluding IBS-M and comparing baseline profiles
#      of controls who developed IBS at follow-up vs consistent controls.
#
# Expected Inputs:
#   - `annotation.csv`: Peptide taxonomic and source database annotations.
#   - `pi_ibs_markers.csv`: Post-infectious IBS marker annotations.
#   - `food_categories.csv`: Food allergen category annotations.
#   - `1_phenos_ibs.rds`: QC'd baseline PhIP-Seq dataset.
#   - `1_phenos_fup.rds`: QC'd follow-up PhIP-Seq dataset.
#   - `Association_IBS_subtypes.csv`: Subtype association results.
#   - `Association_ibs.csv`: General IBS association results.
#   - `Association_age.csv`: Age association results.
#   - `Association_sex.csv`: Sex association results.
#   - `Association_baseline_fup.csv`: Baseline developed-IBS results.
#   - `Association_IBS_no_d.csv`: IBS excluding IBS-M sensitivity results.
#
# Generated Outputs:
#   - `extended_table_2.xlsx` (IBS vs Control associations)
#   - `extended_table_3.xlsx` (IBS Subtypes vs Control associations)
#   - `extended_table_4.xlsx` (PI-IBS marker table)
#   - `extended_table_5.xlsx` (Food allergen category table)
#   - `prevalence_IBS_noIBS.csv`
#   - `prevalence_developed_ctrl_fup.csv`
#   - `figure_2.pdf`
#   - `figure_3.pdf`
#   - `extended_figure_3.pdf`
#   - `extended_figure_5.pdf`
#   - `extended_figure_6.pdf`
#   - `extended_figure_7.pdf`
#   - `fup_plot_all.pdf`
#   - `bacteroidales_no_ibsm_sensitivity.pdf`
#   - `ctrls_develop_baseline.pdf`
# ==============================================================================

library(data.table)
library(readr)
library(ggplot2)
library(fgsea)
library(broom)
library(patchwork)
library(ggrepel)
library(forcats)
library(ggsignif)
library(paletteer)
library(psych)
library(ggpubr)
library(coin)
library(writexl)
library(vegan)
library(ggnewscale)
library(ggalluvial)
library(stringr)

source("0_utils.R")

set.seed(69420)

# 1. Load Data & Export Extended Data Tables 2-5 -------------------------------

annotation <- fread("annotation.csv")
pi_markers <- fread("pi_ibs_markers.csv")
food <- fread("food_categories.csv")
phenos <- read_rds("1_phenos_ibs.rds")
phenos_fup <- read_rds("1_phenos_fup.rds")

res_sub <- fread("Association_IBS_subtypes.csv")
res_ibs <- fread("Association_ibs.csv")
res_age <- fread("Association_age.csv")
res_sex <- fread("Association_sex.csv")
res_base <- fread("Association_baseline_fup.csv")
res_nod <- fread("Association_IBS_no_d.csv")

phenos_fup_all <- phenos_fup |> copy()

write_xlsx(
  res_sub[baseline == "noIBS"],
  path = "extended_table_3.xlsx"
)
write_xlsx(
  res_ibs,
  path = "extended_table_2.xlsx"
)
write_xlsx(
  pi_markers,
  path = "extended_table_4.xlsx"
)
food <- food[category != ""]
write_xlsx(
  food,
  path = "extended_table_5.xlsx"
)


# 2. Prevalence, Ratio Testing, & Enrichment Helper Functions ------------------

calc_ratios <- function(
  prev_dat,
  by,
  base,
  compare = prevalence,
  min_size = 4
) {
  by <- deparse(substitute(by))
  base <- deparse(substitute(base))
  compare <- deparse(substitute(compare))

  prev_dat <- prev_dat[, c("peptide", "case", ..by, ..compare, ..base)] |>
    unique()
  prev_dat[, by = c("case", by), pass := .N > min_size]

  prev_dat[
    pass == TRUE,
    .(case, by, peptide, ratio = get(compare) / get(base)),
    env = list(by = by)
  ][is.finite(ratio)]
}

calc_statistics <- function(ratio_table) {
  by <- names(ratio_table)[2]

  ratio_table <- ratio_table[, by = c("case", by), {
    current_case <- .BY$case
    current_taxon <- .BY[[by]]

    rest <- ratio_table[
      case == current_case & ..by != current_taxon,
      .(ratio, group = factor("rest"))
    ]
    table <- .SD[, .(ratio, group = factor("test"))]
    table <- rbind(table, rest)

    stats <- wilcox_test(ratio ~ group, data = table)
    .(p = pvalue(stats), direction = sign(statistic(stats)))
  }]

  ratio_table[, by = case, FDR := p.adjust(p, "fdr")][order(FDR)]
}

filter_hits <- function(stats, p_cut = p_cutoff, top_n = 5) {
  stats[
    stats[
      FDR < p_cut,
      .I[frankv(FDR, ties.method = "first") <= top_n],
      by = case
    ]$V1
  ]
}

plot_prev <- function(
  prev_table,
  top_hits,
  base,
  by,
  compare = prevalence,
  hex = TRUE
) {
  base <- deparse(substitute(base))
  by <- deparse(substitute(by))
  compare <- deparse(substitute(compare))

  dat <- prev_table[, !"FDR"][top_hits, on = .NATURAL] |> droplevels()

  up <- top_hits[direction == 1]
  up <- up[
    order(FDR),
    keyby = .(case),
    .(get(by), FDR, direction = "\u2191", seq = seq(get(by)) / 14)
  ]

  down <- top_hits[direction == -1]
  down <- down[
    order(-FDR),
    keyby = .(case),
    .(get(by), FDR, direction = "\u2193", seq = seq(get(by)) / 14)
  ]

  bg <- if (hex) {
    geom_hex(data = prev_table[, !..by], fill = "gray", show.legend = FALSE)
  } else {
    geom_point(data = prev_table[, !..by], col = "gray", show.legend = FALSE)
  }

  dat |>
    ggplot(aes(x = (!!sym(base) * 100), y = 100 * !!sym(compare))) +
    bg +
    geom_point(
      aes(col = !!sym(by)),
      alpha = 0.8
    ) +
    geom_text(
      data = up,
      aes(
        x = 100 * .2,
        y = 100 * (1.1 - seq),
        label = paste0(
          direction,
          " FDR = ",
          formatC(FDR, format = "e", digits = 2)
        ),
        colour = V1
      ),
      size = 3,
      show.legend = FALSE
    ) +
    geom_text(
      data = down,
      aes(
        x = 100 * .75,
        y = 100 * (seq - .1),
        label = paste0(
          direction,
          " FDR = ",
          formatC(FDR, format = "e", digits = 2)
        ),
        colour = V1
      ),
      size = 3,
      show.legend = FALSE
    ) +
    facet_wrap(~case, ncol = 2) +
    scale_color_paletteer_d("MoMAColors::Klein") +
    theme_my() +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.225)
    ) +
    coord_cartesian(clip = "off")
}

run_fgsea <- function(results, peptide_list) {
  if ("statistic" %notin% names(results)) {
    results <- results[, .(
      statistic = sign(estimate) * -log10(p.value),
      peptide
    )] |>
      unique()
  } else {
    results <- results[, .(statistic, peptide)] |> unique()
  }

  set.seed(1)

  if (any(duplicated(results$statistic))) {
    dups <- duplicated(results$statistic) |
      duplicated(results$statistic, fromLast = TRUE)
    results$statistic[dups] <- results$statistic[dups] +
      runif(sum(dups), -1e-6, 1e-6)
  }

  results <- results[, .(stats = statistic, peptide)] |> unique()
  ranks <- setNames(results$stats, results$peptide)

  fgseaRes <- fgsea::fgsea(
    peptide_list,
    ranks,
    minSize = 4,
    nPermSimple = 100000
  )
  fgseaRes[order(padj)]
}

plot_fgsea <- function(fgsea_res, fgsea_top) {
  dat <- fgsea_res[fgsea_top, on = .(pathway), allow.cartesian = TRUE]
  dat |>
    ggplot() +
    scale_fill_distiller(palette = "BrBG") +
    geom_point(
      aes(case, pathway, fill = NES, size = pmin(5, -log10(FDR))),
      shape = 21,
      color = "grey30"
    ) +
    theme_my(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x = element_blank()
    ) +
    guides(
      size = guide_legend(
        title = expression(bold(-log[10])(bolditalic(FDR))),
        override.aes = list(shape = 16)
      )
    ) +
    coord_cartesian(clip = "off")
}

plot_db <- function(ratios, statistics, reorder = TRUE) {
  dat <- ratios[statistics, on = .NATURAL]
  dat[, ratio := log2(ratio)]
  if (reorder) {
    dat[,
      source_database := fct_reorder(
        source_database,
        ratio,
        .fun = function(x) mean(x[is.finite(x)])
      )
    ]
  }
  dat[, case := factor(case)]

  stats <- dat[
    FDR < .05,
    by = .(case, source_database, FDR),
    .(ratio = max(ratio))
  ] |>
    unique()
  stats[, let(
    start = as.numeric(source_database) - .2,
    end = as.numeric(source_database) + .2,
    label = fcase(
      FDR < .001 , "***" ,
      FDR < .01  , "**"  ,
      FDR < .05  , "*"
    ),
    y = ratio + .6
  )]

  dat[is.finite(ratio)] |>
    ggplot(aes(source_database, ratio)) +
    geom_boxplot(aes(fill = source_database), outliers = FALSE) +
    geom_jitter(height = 0, width = 0.3, size = .2, alpha = .3) +
    geom_signif(
      data = stats,
      aes(
        xmin = start,
        xmax = end,
        annotations = label,
        y_position = y,
        group = source_database,
        col = fifelse(label == "*", "grey30", "red")
      ),
      manual = TRUE,
      textsize = 4,
      vjust = .3,
      tip_length = NA
    ) +
    facet_wrap(~case, ncol = 1) +
    scale_fill_paletteer_d("MoMAColors::Klein") +
    theme_my(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        hjust = 0,
        vjust = -.3,
        margin = margin(b = 4)
      ),
      panel.border = element_blank(),
      axis.title.x = element_blank()
    ) +
    ylim(NA, 3.5) +
    scale_color_identity()
}

# Helper to maintain consistent taxonomic colors across figures
extend_taxa_palette <- function(
  existing_cols,
  target_taxa,
  fallback_cols = NULL
) {
  overlap <- existing_cols[names(existing_cols) %in% target_taxa]
  no_overlap <- sort(unique(target_taxa[target_taxa %notin% names(existing_cols)]))

  pal <- as.character(paletteer::paletteer_d("MoMAColors::Klein"))
  avail_col <- length(pal) - length(existing_cols) - length(no_overlap)
  if (avail_col < 0 && !is.null(fallback_cols)) {
    pal <- c(pal, fallback_cols)
  }

  remain_pal <- pal[pal %notin% existing_cols]
  new_cols <- setNames(remain_pal[seq_along(no_overlap)], no_overlap)

  list(
    fin_col = c(overlap, new_cols),
    new_cols = new_cols,
    all_cols = c(existing_cols, new_cols)
  )
}


# 3. Cross-Sectional IBS Subtypes vs Controls ----------------------------------

## 3a. Calculate Per-Subtype Seroprevalence ------------------------------------
prev_sub <- calc_prev(phenos, group_var = "IBS_Subtype")
prev_sub <- prev_sub[, melt(
  .SD,
  measure = patterns("^IBS"),
  var = "case",
  val = "prevalence"
)]
prev_sub <- prev_sub[annotation, on = .(peptide), nomatch = 0]

## 3b. Taxon-Level Prevalence Ratios -------------------------------------------
p_cutoff <- .01
ratios <- calc_ratios(prev_sub, taxa, `no IBS`) |> unique()
ratio_stats <- calc_statistics(ratios)
top_hits <- filter_hits(ratio_stats, top_n = 5)

prev_plot <- plot_prev(
  prev_sub,
  top_hits,
  base = `no IBS`,
  by = taxa,
  hex = FALSE
) +
  guides(
    color = guide_legend(
      title = "Taxonomic group",
      keywidth = .5,
      keyheight = .8,
      default.unit = "lines"
    )
  ) +
  xlim(c(0, 100)) +
  xlab(expression("Control individuals with significant peptide binding [%]")) +
  ylab(expression("IBS individuals with significant peptide binding [%]")) +
  ggtitle(expression(bold("IBS Subtypes vs. Controls")))

# Extract palette mapping so downstream plots color matching taxa identically
bld <- ggplot_build(prev_plot)
sc <- bld$plot$scales$get_scales("colour")
levs <- sc$get_limits()
colors <- setNames(sc$map(levs), levs)

## 3c. Taxonomic Gene Set Enrichment Analysis (fGSEA) --------------------------
pathway_list <- split(annotation, by = "taxa") |>
  lapply(function(x) unique(x[["peptide"]]))

fgsea_sub <- res_sub[
  baseline == "noIBS",
  by = term,
  run_fgsea(.SD, peptide_list = pathway_list)
]
setnames(fgsea_sub, c("term", "padj"), c("case", "FDR"))

top_fgsea <- filter_hits(fgsea_sub)
fgsea_plot <- plot_fgsea(fgsea_sub, top_fgsea) +
  ylab(expression("Top enriched taxonomic groups"))

## 3d. Source Database Prevalence Log-Ratios -----------------------------------
prev_sub[
  source_database == "Antibody-coated bacteria",
  source_database := "AB-coated bacteria"
]
ratios_db <- calc_ratios(prev_sub, source_database, `no IBS`)
ratio_stats_db <- calc_statistics(ratios_db)

db_plot <- plot_db(ratios_db, ratio_stats_db) +
  ylab(expression(atop(
    "Log-ratio of significantly bound antigens",
    log[2] ~ ("%-Subtype / %-Control")
  ))) +
  theme(axis.title.x = element_blank())

xorder <- levels(db_plot$data$source_database)


# 4. Post-Infectious IBS Markers & Food Allergens (Extended Figures 5 & 6) -----

## 4a. PI-IBS Marker Enrichment (Extended Data Figure 5) -----------------------
pi_list <- split(pi_markers, by = "source_database") |>
  lapply(function(x) unique(x[["peptide"]]))
pi_list_full <- append(setNames(list(pi_markers$peptide), "PI-IBS marker"), pi_list)

fgsea_pi <- res_sub[
  baseline == "noIBS",
  by = term,
  run_fgsea(.SD, peptide_list = pi_list_full)
]
setnames(fgsea_pi, c("term", "padj"), c("case", "FDR"))

fgsea_pi_ibs <- res_ibs[
  baseline == "noIBS",
  by = term,
  run_fgsea(.SD, peptide_list = pi_list_full)
]
setnames(fgsea_pi_ibs, c("term", "padj"), c("case", "FDR"))

fgsea_pi <- rbind(fgsea_pi, fgsea_pi_ibs)
top_fgsea_pi <- filter_hits(fgsea_pi, p = 1)
fgsea_plot_pi <- plot_fgsea(fgsea_pi, top_fgsea_pi)

ggsave(
  plot = fgsea_plot_pi,
  "extended_figure_5.pdf",
  width = 5
)

fgsea_pi[, leadingEdge := NULL]
setorder(fgsea_pi, "FDR")

## 4b. Food Allergen Category Enrichment (Extended Data Figure 6) --------------
food_cat <- split(food, by = "category") |>
  lapply(function(x) unique(x[["peptide"]]))

fgsea_food <- res_sub[
  baseline == "noIBS",
  by = term,
  run_fgsea(.SD, peptide_list = food_cat)
]
setnames(fgsea_food, c("term", "padj"), c("case", "FDR"))

top_fgsea_food <- filter_hits(fgsea_food, p = 1)
fgsea_plot_food <- plot_fgsea(fgsea_food, top_fgsea_food)

ggsave(
  plot = fgsea_plot_food,
  "extended_figure_6.pdf",
  width = 5
)

fgsea_food[, leadingEdge := NULL]
setorder(fgsea_food, "FDR")


# 5. Longitudinal Repertoire Stability & Trajectory Analysis -------------------

## 5a. Pairwise Phi-Coefficient Matrix Between Baseline and Follow-Up ----------
cols <- intersect(names(phenos_fup), names(phenos))

combined <- rbindlist(
  list(
    fup = phenos_fup[, ..cols],
    base = phenos[ID %in% phenos_fup$ID, ..cols]
  ),
  fill = TRUE,
  use.names = TRUE,
  idcol = "tp"
)

full <- combined[, .SD, .SDcols = patterns("^ID|^agi|^twi|^tp")]
full <- full[, melt(
  .SD,
  id = c("ID", "tp"),
  var = "peptide",
  variable.factor = FALSE
)][!is.na(value)]
full[, value := factor(value)]

d_base <- full[tp == "base", .(ID_base = ID, peptide, value_base = value)]
d_fup <- full[tp == "fup", .(ID_fup = ID, peptide, value_fup = value)]
setkey(d_base, peptide)
setkey(d_fup, peptide)

full <- d_base[d_fup, on = "peptide", allow.cartesian = TRUE, nomatch = 0L]
d_base <- d_fup <- NULL

setcolorder(full, c("ID_base", "ID_fup", "peptide", "value_base", "value_fup"))
setkey(full, ID_base, ID_fup)

mat <- full[,
  keyby = .(ID_base, ID_fup),
  .(phi = phi(table(value_fup, value_base), digits = 3))
]

# Order samples by subtype so subtype blocks are contiguous along both axes
samp <- unique(phenos_fup[, .(ID, IBS_Subtype)])
samp[,
  IBS_Subtype := factor(
    IBS_Subtype,
    levels = c("IBS-C", "IBS-D", "IBS-M", "no IBS")
  )
]
setorder(samp, IBS_Subtype, ID)
samp[, `:=`(x = .I, y = .I)]

hm_dat <- copy(mat)
hm_dat[samp, on = .(ID_base = ID), x := i.x]
hm_dat[samp, on = .(ID_fup = ID), y := i.y]

n_id <- nrow(samp)
blk_x <- samp[, .(lo = min(x) - .5, hi = max(x) + .5), by = IBS_Subtype]
blk_y <- samp[, .(lo = min(y) - .5, hi = max(y) + .5), by = IBS_Subtype]

bar <- max(1.5, n_id * 0.02)
edge <- 0

fup_hm <- ggplot() +
  geom_tile(data = hm_dat, aes(x, y, fill = phi)) +
  scale_fill_distiller(
    palette = "BrBG",
    name = expression(bold("Phi-coef") ~ (Phi)),
    guide = guide_colourbar(
      position = "top",
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  ggnewscale::new_scale_fill() +
  geom_rect(
    data = blk_x,
    aes(
      xmin = lo,
      xmax = hi,
      ymin = edge - bar,
      ymax = edge,
      fill = IBS_Subtype
    ),
    colour = "grey30",
    linewidth = 0
  ) +
  geom_rect(
    data = blk_y,
    aes(
      ymin = lo,
      ymax = hi,
      xmin = edge - bar,
      xmax = edge,
      fill = IBS_Subtype
    ),
    colour = "grey30",
    linewidth = 0
  ) +
  scale_fill_paletteer_d(
    palette = "MoMAColors::Klein",
    name = "Follow-up IBS Status",
    guide = guide_legend(position = "bottom", nrow = 1, title.position = "left")
  ) +
  theme_my() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.frame = element_rect(color = "grey30"),
    legend.ticks = element_line(color = "grey30"),
    axis.line = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(t = 10, r = 5, b = 10, l = 5)
  ) +
  coord_equal(clip = "off") +
  xlab("Baseline antibody repertoire") +
  ylab("Follow-up antibody repertoire") +
  ggtitle("Antibody repertoire stability")

## 5b. Node-Level Permutation Test & Stability Comparisons (Extended Fig 7) ----
mat[,
  same := fifelse(
    ID_base == ID_fup,
    "Same individuals",
    "Different individuals"
  )
]
mat[
  phenos_fup,
  on = .(ID_base = ID),
  let(
    IBS_Subtype = factor(
      IBS_Subtype,
      levels = c("IBS-C", "IBS-D", "IBS-M", "no IBS")
    ),
    follow_up_grp = follow_up_grp
  )
]

mat[
  same == "Different individuals",
  pvalue(independence_test(
    phi ~ IBS_Subtype,
    distribution = approximate(nresample = 10000)
  ))
]

paired <- mat[same == "Different individuals"]
paired[, distance := 1 - phi]

obs_medians <- paired[, .(med = median(phi)), by = IBS_Subtype]$med
obs_stat <- var(obs_medians)

individuals <- unique(paired[, .(ID_base, IBS_Subtype)])

set.seed(123)
n_perms <- 1000
perm_stats <- numeric(n_perms)

for (i in 1:n_perms) {
  shuffled_individuals <- copy(individuals)
  shuffled_individuals[, perm_IBS := sample(IBS_Subtype)]

  temp_pairs <- merge(
    paired[, .(ID_base, ID_fup, phi)],
    shuffled_individuals[, .(ID_base, perm_IBS)],
    by = "ID_base"
  )

  fake_medians <- temp_pairs[, .(med = median(phi)), by = perm_IBS]$med
  perm_stats[i] <- var(fake_medians)
}

p_value <- sum(perm_stats >= obs_stat) / n_perms
cat("Observed Variance of Medians:", obs_stat, "\n")
cat("Global Permuted P-value:", p_value, "\n")

mat[, keyby = .(same, IBS_Subtype), .(mean(phi), sd(phi))]

ex_long_1 <- mat |>
  ggplot(aes(y = phi, fill = same, x = IBS_Subtype)) +
  ggdist::stat_halfeye(
    aes(slab_colour = after_scale(fill)),
    adjust = .6,
    slab_alpha = .4,
    .width = 0,
    height = .8,
    point_colour = NA,
    justification = -.1,
    show.legend = FALSE
  ) +
  theme_my(panel.grid.major.y = element_blank()) +
  geom_boxplot(
    color = "grey30",
    width = 0.12,
    outlier.color = NA,
    alpha = 0.4
  ) +
  stat_compare_means(
    aes(x = IBS_Subtype, y = phi, group = same),
    paired = FALSE,
    method = "wilcox.test",
    bracket.size = NA,
    vjust = 0,
    hjust = -.75,
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    ),
    hide.ns = FALSE
  ) +
  scale_fill_paletteer_d(palette = "MoMAColors::Klein") +
  scale_color_paletteer_d(palette = "MoMAColors::Klein") +
  guides(
    fill = guide_legend(
      title = "Comparison",
      nrow = 1,
      byrow = TRUE,
      override.aes = list(colour = NA),
      position = "bottom"
    ),
    colour = guide_legend(title = "Study group")
  ) +
  scale_x_discrete(expand = expansion(add = c(0.25, 0))) +
  scale_y_continuous(expand = expansion(mult = c(.01, 0.1))) +
  coord_flip(clip = "off") +
  labs(
    y = expression(
      "Similarity between follow-up and baseline samples" ~ ("phi-coef" ~ Phi)
    ),
    x = element_blank()
  )

base_fup <- phenos_fup |> copy()
base_fup[phenos, on = .(ID), n_baseline := i.n_distinct_peptides]
base_fup <- base_fup[, .(ID, IBS_Subtype, n_distinct_peptides, n_baseline)]
base_fup <- base_fup[, melt(.SD, measure = patterns("^n_"))]
base_fup[,
  IBS_Subtype := factor(
    IBS_Subtype,
    levels = c("IBS-C", "IBS-D", "IBS-M", "no IBS")
  )
]

ex_long_3 <- base_fup |>
  ggplot(aes(y = value, fill = variable, x = IBS_Subtype)) +
  geom_boxplot(
    width = .3,
    position = position_dodge(width = .4),
    outliers = FALSE
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = .25,
      jitter.height = 0,
      dodge.width = .4
    ),
    size = .2,
    alpha = .3
  ) +
  theme_my(panel.grid.major.y = element_blank()) +
  scale_fill_paletteer_d(
    palette = "MoMAColors::Klein",
    labels = c("Follow-Up", "Baseline")
  ) +
  stat_compare_means(
    aes(y = value, x = IBS_Subtype, group = variable),
    paired = TRUE,
    method = "wilcox.test",
    bracket.size = NA,
    vjust = 0,
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    ),
    hide.ns = FALSE
  ) +
  coord_flip() +
  guides(fill = guide_legend(title = "Timepoint", reverse = TRUE)) +
  xlab(element_blank()) +
  ylab("Distinct Antibodies Bound")

base_fup[,
  by = IBS_Subtype,
  wilcox.test(
    .SD[variable == "n_baseline", value],
    .SD[variable != "n_baseline", value],
    paired = TRUE
  )$p.value
]

ex_long_1 /
  ex_long_3 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave(
  "extended_figure_7.pdf",
  width = 7,
  height = 8
)

## 5c. IBS Subtype Trajectories Alluvial Plot (Baseline -> Follow-Up) ----------
traj_counts <- phenos_fup_all[, by = .(IBS_Subtype_base, IBS_Subtype), .N]

# Exclude stable no IBS -> no IBS so symptomatic trajectories are clearly visible
traj_allu <- traj_counts[
  !(IBS_Subtype_base == "no IBS" & IBS_Subtype == "no IBS")
]

ibs_levels <- c("IBS-C", "IBS-D", "IBS-M", "no IBS")
.extra_levels <- setdiff(
  union(
    as.character(traj_allu$IBS_Subtype_base),
    as.character(traj_allu$IBS_Subtype)
  ),
  ibs_levels
)
sub_levels <- c(ibs_levels, .extra_levels)

traj_allu[, IBS_Subtype_base := factor(IBS_Subtype_base, levels = sub_levels)]
traj_allu[, IBS_Subtype := factor(IBS_Subtype, levels = sub_levels)]

klein_pal <- paletteer::paletteer_d("MoMAColors::Klein")
fill_cols <- setNames(
  c(klein_pal[seq_along(ibs_levels)], klein_pal[5:(4 + length(.extra_levels))]),
  sub_levels
)

fup_grp_cols <- c(
  "Baseline" = "grey95",
  "Improved" = "#9CCB86",
  "Persistent / switched" = "#CF597E"
)

traj_alluvial_plot <- ggplot(
  traj_allu,
  aes(axis1 = IBS_Subtype_base, axis2 = IBS_Subtype, y = N)
) +
  geom_alluvium(aes(fill = IBS_Subtype_base), width = 1 / 8, alpha = 0.8) +
  scale_fill_manual(values = fill_cols, name = "Baseline subtype") +
  ggnewscale::new_scale_fill() +
  geom_stratum(
    aes(
      fill = after_stat(ifelse(
        x == 1,
        "Baseline",
        ifelse(
          as.character(stratum) == "no IBS",
          "Improved",
          "Persistent / switched"
        )
      ))
    ),
    width = 1 / 8,
    colour = "grey30"
  ) +
  scale_fill_manual(
    values = fup_grp_cols,
    breaks = c("Improved", "Persistent / switched"),
    name = "Follow-up group"
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.2) +
  scale_x_discrete(
    limits = c("Baseline", "Follow-up"),
    expand = expansion(add = c(0.12, 0.12))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme_my(panel.grid.major.y = element_blank(), legend.position = "bottom") +
  labs(x = element_blank(), y = "Individuals (n)") +
  coord_cartesian(clip = "off")

## 5d. Follow-Up Cohort 1: Persistent / Switched IBS ---------------------------
phenos_fup <- phenos_fup_all[follow_up_grp != "improvement"] |> copy()
cols <- intersect(names(phenos_fup), names(phenos))

phenos_base_any <- phenos |> copy()
phenos_base_any[, IBS_Subtype := "any"]

combined <- rbindlist(
  list(
    fup = phenos_fup[, ..cols],
    base = phenos_base_any[ID %in% phenos_fup$ID, ..cols]
  ),
  fill = TRUE,
  use.names = TRUE,
  idcol = "tp"
)

prev_fup <- calc_prev(combined, group_var = c("IBS_Subtype", "tp"))
prev_fup <- prev_fup[, melt(
  .SD,
  measure = patterns("IBS"),
  var = "case",
  val = "prevalence"
)]
prev_fup <- prev_fup[annotation, on = .(peptide), nomatch = 0] |> droplevels()

setnames(prev_fup, "any_base", "bl")
ratios_fup <- calc_ratios(prev_fup, taxa, base = bl) |> unique()
ratio_stats_fup <- calc_statistics(ratios_fup)
top_hits_fup <- filter_hits(ratio_stats_fup, p_cut = .01)

pal_fup1 <- extend_taxa_palette(
  existing_cols = colors,
  target_taxa = top_hits_fup$taxa,
  fallback_cols = c("deepskyblue2", "firebrick")
)
fin_col <- pal_fup1$fin_col
new_cols <- pal_fup1$new_cols

prev_plot_fup <- plot_prev(
  prev_fup[case %notin% "no IBS_fup"],
  top_hits_fup[case %notin% "no IBS_fup"],
  base = bl,
  by = taxa,
  hex = FALSE
) +
  guides(
    color = guide_legend(
      title = "Taxonomic group",
      keywidth = .5,
      keyheight = .8,
      default.unit = "lines"
    )
  ) +
  scale_color_manual(values = fin_col) +
  xlab(expression(
    "Baseline individuals with significant peptide binding [%]"
  )) +
  ylab(expression(
    "Follow-up individuals with significant peptide binding [%]"
  )) +
  ggtitle(expression(bold("Baseline vs. Follow-up"))) +
  coord_fixed(clip = "off")

prev_fup[
  source_database == "Antibody-coated bacteria",
  source_database := "AB-coated bacteria"
]
prev_fup[, source_database := factor(source_database, levels = xorder)]
ratios_db_fup <- calc_ratios(prev_fup, source_database, base = bl)
ratio_stats_db_fup <- calc_statistics(ratios_db_fup)

db_plot_fup <- plot_db(
  ratios_db_fup[case %notin% "no IBS_fup"],
  ratio_stats_db_fup[case %notin% "no IBS_fup"],
  reorder = FALSE
) +
  ylab(expression(atop(
    "Log-ratio of significantly bound antigens",
    log[2] ~ ("%-Follow-up / %-Baseline")
  )))

fup_plot_all <- prev_plot_fup + db_plot_fup
ggsave(
  plot = fup_plot_all,
  filename = "fup_plot_all.pdf",
  width = 9
)

## 5e. Follow-Up Cohort 2: Improved (Resolved IBS at Follow-Up) ----------------
phenos_fup <- phenos_fup_all[follow_up_grp == "improvement"]
cols <- intersect(names(phenos_fup), names(phenos))

combined <- rbindlist(
  list(
    fup = phenos_fup[, ..cols],
    base = phenos[ID %in% phenos_fup$ID, ..cols]
  ),
  fill = TRUE,
  use.names = TRUE,
  idcol = "tp"
)

prev_fup <- calc_prev(combined, group_var = c("IBS_Subtype", "tp"))
prev_fup <- prev_fup[, melt(
  .SD,
  measure = patterns("^IBS"),
  var = "case",
  val = "prevalence"
)]
prev_fup <- prev_fup[annotation, on = .(peptide), nomatch = 0] |> droplevels()

setnames(prev_fup, "no IBS_fup", "bl")
ratios_fup <- calc_ratios(prev_fup, taxa, base = prevalence, compare = bl) |>
  unique()
ratio_stats_fup <- calc_statistics(ratios_fup)
top_hits_fup <- filter_hits(ratio_stats_fup, p_cut = .05)

colors <- c(colors, new_cols)
pal_fup2 <- extend_taxa_palette(
  existing_cols = colors,
  target_taxa = top_hits_fup$taxa,
  fallback_cols = c("deepskyblue2", "firebrick", "springgreen2")
)
fin_col <- pal_fup2$fin_col
new_cols <- pal_fup2$new_cols

prev_plot_improve <- plot_prev(
  prev_fup,
  top_hits_fup,
  base = prevalence,
  compare = bl,
  by = taxa,
  hex = FALSE
) +
  guides(
    color = guide_legend(
      title = "Taxonomic group",
      keywidth = .5,
      keyheight = .8,
      default.unit = "lines"
    )
  ) +
  scale_color_manual(values = fin_col) +
  xlab(expression(
    "Baseline individuals with significant peptide binding [%]"
  )) +
  ylab(expression(
    "Follow-up individuals with significant peptide binding [%]"
  )) +
  ggtitle(expression(bold("Baseline vs. Follow-up"))) +
  coord_fixed(clip = "off")

prev_fup[
  source_database == "Antibody-coated bacteria",
  source_database := "AB-coated bacteria"
]
prev_fup[, source_database := factor(source_database, levels = xorder)]
ratios_db_fup <- calc_ratios(
  prev_fup,
  source_database,
  base = prevalence,
  compare = bl
)
ratio_stats_db_fup <- calc_statistics(ratios_db_fup)

db_plot_improve <- plot_db(ratios_db_fup, ratio_stats_db_fup, reorder = FALSE) +
  ylab(expression(atop(
    "Log-ratio of significantly bound antigens",
    log[2] ~ ("%-Follow-up / %-Baseline")
  )))


# 6. Age and Sex Taxonomic Enrichment (Extended Data Figure 3) -----------------

pathway_list <- split(annotation, by = "taxa") |>
  lapply("[[", "peptide")
res_age_sex <- rbindlist(list(age = res_age, sexM = res_sex), idcol = "case")
res_age_sex[p.value <= 0, p.value := .Machine$double.xmin * 10^258]

fgsea_age_sex <- res_age_sex[,
  by = case,
  run_fgsea(.SD, peptide_list = pathway_list)
]
setnames(fgsea_age_sex, c("padj"), c("FDR"))

p_cutoff <- .05
top_fgsea_age_sex <- filter_hits(fgsea_age_sex, top_n = 5)
fgsea_plot_age_sex <- plot_fgsea(fgsea_age_sex, top_fgsea_age_sex)

ggsave(
  plot = fgsea_plot_age_sex,
  "extended_figure_3.pdf"
)


# 7. Cross-Sectional General IBS vs Controls & Figure 2 / Figure 3 Assembly ----

## 7a. Per-Peptide Seroprevalence (IBS vs no IBS) ------------------------------
prev_ibs <- calc_prev(phenos, group_var = "IBS_Status")
prev_ibs <- prev_ibs[, melt(
  .SD,
  measure = patterns("^IBS"),
  var = "case",
  val = "prevalence"
)]
prev_ibs <- prev_ibs[annotation, on = .(peptide), nomatch = 0]

fwrite(prev_ibs, "prevalence_IBS_noIBS.csv")

## 7b. Taxon-Level Prevalence Ratios (IBS vs no IBS) ---------------------------
p_cutoff <- 0.01
ratios_ibs <- calc_ratios(prev_ibs, taxa, `no IBS`, min_size = 4)
ratio_stats_ibs <- calc_statistics(ratios_ibs)
top_hits_ibs <- filter_hits(ratio_stats_ibs)

new_cols_full <- c(colors, new_cols)

prev_plot_ibs <- plot_prev(
  prev_ibs,
  top_hits_ibs,
  base = `no IBS`,
  by = taxa,
  hex = FALSE
) +
  guides(
    color = guide_legend(
      title = "Taxonomic group",
      keywidth = .5,
      keyheight = .8,
      default.unit = "lines",
      position = "bottom"
    )
  ) +
  scale_color_manual(values = new_cols_full) +
  xlab(expression("Control individuals with significant peptide binding [%]")) +
  ylab(expression("IBS individuals with significant peptide binding [%]")) +
  ggtitle(expression(bold("IBS vs. Controls")))

## 7c. Taxonomic fGSEA & Bacteroidales Running Enrichment Plot -----------------
fgsea_ibs <- res_ibs[
  baseline == "noIBS",
  by = term,
  run_fgsea(.SD, peptide_list = pathway_list)
]
setnames(fgsea_ibs, c("term", "padj"), c("case", "FDR"))

top_fgsea_ibs <- filter_hits(fgsea_ibs, .01)
gradient_range <- ggplot_build(fgsea_plot)$plot$scales$get_scales(
  "fill"
)$range$range

fgsea_plot_ibs <- plot_fgsea(fgsea_ibs, top_fgsea_ibs) +
  ylab("Top enriched taxonomic groups") +
  scale_fill_distiller(palette = "BrBG", limits = gradient_range)

stats_ibs <- res_ibs[, .(
  stats = sign(estimate) * -log10(p.value),
  peptide
)] |>
  unique()
ranks <- setNames(stats_ibs$stats, stats_ibs$peptide)

enrich_plot_ibs <- plotEnrichment(
  pathway_list$Bacteroidales,
  ranks,
  ticksSize = .1
) +
  theme_my(panel.grid.major.y = element_blank()) +
  labs(
    title = "Bacteroidales Enrichment Score",
    x = "Rank in Ordered Dataset",
    y = "Enrichment Score (ES)"
  ) +
  scale_x_reverse() +
  coord_flip()

## 7d. Source Database Log-Ratios (IBS vs no IBS) ------------------------------
prev_ibs[
  source_database == "Antibody-coated bacteria",
  source_database := "AB-coated bacteria"
]
prev_ibs[, source_database := factor(source_database, levels = xorder)]
ratios_db_ibs <- calc_ratios(prev_ibs, source_database, `no IBS`)
ratio_stats_db_ibs <- calc_statistics(ratios_db_ibs)

db_plot_ibs <- plot_db(ratios_db_ibs, ratio_stats_db_ibs, reorder = FALSE) +
  ylab(expression(atop(
    "Log-ratio of significantly bound antigens",
    log[2] ~ ("%-Subtype / %-Control")
  ))) +
  theme(strip.text = element_blank())

## 7e. Assemble Manuscript Figure 2 (Cross-Sectional) --------------------------
fig2_row1 <- free(enrich_plot_ibs) +
  free(prev_plot_ibs) +
  db_plot_ibs +
  plot_layout(widths = c(1.2, 1, 1.1))

fig2_row2 <- fgsea_plot +
  free(prev_plot, type = "label") +
  free(db_plot) +
  plot_layout(widths = c(.25, 1.5, 1.1))

free(fig2_row1) /
  free(fig2_row2) +
  plot_layout(heights = c(.8, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "#fbf9f4", colour = NA)
  )

ggsave("figure_2.pdf", width = 15, height = 10)

## 7f. Assemble Manuscript Figure 3 (Longitudinal Follow-Up) -------------------
fup_gutter_label <- function(text, col) {
  wrap_elements(
    grid::textGrob(
      text,
      x = 0.5,
      y = 0.5,
      rot = 90,
      gp = grid::gpar(col = unname(col), fontface = "bold", cex = 0.95)
    )
  )
}

fig_fup_design <- "
#AB
GCD
HEF
"

fig_fup <-
  free(traj_alluvial_plot) +
  free(fup_hm) +
  free(prev_plot_fup) +
  db_plot_fup +
  free(prev_plot_improve) +
  db_plot_improve +
  fup_gutter_label(
    "Persistent / switched",
    fup_grp_cols["Persistent / switched"]
  ) +
  fup_gutter_label("Improved", fup_grp_cols["Improved"]) +
  plot_layout(
    design = fig_fup_design,
    widths = c(0.05, 1, 1),
    heights = c(1, 1, 1)
  ) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "", ""))) &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "#fbf9f4", colour = NA)
  )

ggsave(
  "figure_3.pdf",
  plot = fig_fup,
  width = 12.5,
  height = 18
)


# 8. Sensitivity Analyses ------------------------------------------------------

## 8a. IBS Excluding IBS-M Sensitivity Analysis --------------------------------
nod <- phenos[IBS_Subtype != "IBS-M"] |> droplevels()

prev_nod <- calc_prev(nod, group_var = "IBS_Status")
prev_nod <- prev_nod[, melt(
  .SD,
  measure = patterns("^IBS"),
  var = "case",
  val = "prevalence"
)]
prev_nod <- prev_nod[annotation, on = .(peptide), nomatch = 0]

p_cutoff <- 0.05
ratios_nod <- calc_ratios(prev_nod, taxa, `no IBS`, min_size = 4)
ratio_stats_nod <- calc_statistics(ratios_nod)
top_hits_nod <- filter_hits(ratio_stats_nod, top_n = 6)
top_hits_nod <- top_hits_nod[taxa == "Bacteroidales"]

prev_plot_nod <- plot_prev(
  prev_nod,
  top_hits_nod,
  base = `no IBS`,
  by = taxa,
  hex = FALSE
) +
  guides(
    color = guide_legend(
      title = "Taxonomic group",
      keywidth = .5,
      keyheight = .8,
      default.unit = "lines",
      position = "bottom"
    )
  ) +
  xlab(expression("Control individuals with significant peptide binding [%]")) +
  ylab(expression("IBS individuals with significant peptide binding [%]")) +
  ggtitle(expression(bold("IBS vs. Controls")))

fgsea_nod <- res_nod[
  baseline == "noIBS",
  by = term,
  run_fgsea(.SD, peptide_list = pathway_list)
]
setnames(fgsea_nod, c("term", "padj"), c("case", "FDR"))

top_fgsea_nod <- filter_hits(fgsea_nod, .05)
fgsea_plot_nod <- plot_fgsea(fgsea_nod, top_fgsea_nod) +
  ylab("Top enriched taxonomic groups") +
  scale_fill_distiller(palette = "BrBG")

stats_nod <- res_nod[, .(
  stats = sign(estimate) * -log10(p.value),
  peptide
)] |>
  unique()
ranks <- setNames(stats_nod$stats, stats_nod$peptide)

enrich_plot_nod <- plotEnrichment(
  pathway_list$Bacteroidales,
  ranks,
  ticksSize = .1
) +
  theme_my(panel.grid.major.y = element_blank()) +
  labs(
    title = "Bacteroidales Enrichment Score",
    x = "Rank in Ordered Dataset",
    y = "Enrichment Score (ES)"
  ) +
  scale_x_reverse() +
  coord_flip()

prev_plot_nod + fgsea_plot_nod + plot_layout(widths = c(2.2, .4))

ggsave(
  "bacteroidales_no_ibsm_sensitivity.pdf",
  width = 7
)

## 8b. Baseline Profiles of Controls Developing IBS at Follow-Up ---------------
ind <- phenos_fup[IBS_Subtype_base == "no IBS"] |> droplevels()
ind[, consistent := fifelse(IBS_Subtype == "no IBS", "consistent", "developed")]
ind[, consistent := factor(consistent)]

fup_base <- phenos[ind, on = .(ID)]

prev_base <- calc_prev(fup_base, group_var = "consistent")
prev_base <- prev_base[, melt(
  .SD,
  measure = patterns("^dev"),
  var = "case",
  val = "prevalence"
)]
prev_base <- prev_base[annotation, on = .(peptide), nomatch = 0]

fwrite(prev_base, "prevalence_developed_ctrl_fup.csv")

p_cutoff <- 0.05
ratios_base <- calc_ratios(prev_base, taxa, consistent, min_size = 4)
ratio_stats_base <- calc_statistics(ratios_base)
top_hits_base <- filter_hits(ratio_stats_base)

prev_plot_base <- plot_prev(
  prev_base,
  top_hits_base,
  base = consistent,
  by = taxa,
  hex = FALSE
) +
  guides(
    color = guide_legend(
      title = "Taxonomic group",
      keywidth = .5,
      keyheight = .8,
      default.unit = "lines",
      position = "bottom"
    )
  ) +
  xlab(expression(
    "Peptide binding at baseline in controls developing IBS [%]"
  )) +
  ylab(expression("Peptide binding at baseline in consistent conrols [%]")) +
  ggtitle(expression(bold("IBS vs. Controls")))

fgsea_base <- res_base[,
  by = term,
  run_fgsea(.SD, peptide_list = pathway_list)
]
setnames(fgsea_base, c("term", "padj"), c("case", "FDR"))

top_fgsea_base <- filter_hits(fgsea_base, .05)
fgsea_plot_base <- plot_fgsea(fgsea_base, top_fgsea_base) +
  ylab("Top enriched taxonomic groups") +
  scale_fill_distiller(palette = "BrBG")

prev_base[
  source_database == "Antibody-coated bacteria",
  source_database := "AB-coated bacteria"
]
prev_base[, source_database := factor(source_database)]
ratios_db_base <- calc_ratios(prev_base, source_database, consistent)
ratio_stats_db_base <- calc_statistics(ratios_db_base)

db_plot_base <- plot_db(ratios_db_base, ratio_stats_db_base) +
  ylab(expression(atop(
    "Log-ratio of significantly bound antigens",
    log[2] ~ ("%-Subtype / %-Control")
  ))) +
  theme(strip.text = element_blank())

(fgsea_plot_base + prev_plot_base) / db_plot_base
ggsave(
  "ctrls_develop_baseline.pdf",
  width = 10,
  height = 8
)
