# ==============================================================================
# 1_QC.R
# Quality control, binarization, and peptide filtering
# ==============================================================================
#
# Purpose:
#   1. Reads participant-level baseline and follow-up PhIP-Seq tables that
#      already include the covariates used downstream, and labels follow-up
#      IBS trajectories from the baseline and follow-up subtype.
#   2. Summarizes the PhIP-Seq antigen library composition (treemap).
#   3. Converts continuous PhIP-Seq enrichment scores to binary presence/absence
#      calls (>0 -> 1, -1 contamination -> NA, NA -> 0) via Rcpp.
#   4. Applies sample-level QC (>= 200 bound peptides per individual),
#      peptide prevalence filtering (5% <= prevalence < 95%), and amino acid
#      sequence deduplication (retaining the higher-prevalence peptide).
#   5. Compares peptide richness across study groups and generates Extended
#      Data Figure 1 (`extended_figure_1.pdf`).
#   6. Formats PI-IBS marker annotations and saves QC'd datasets for
#      downstream scripts. Peptide taxonomy in
#      `annotation.csv` is supplied with the scripts and
#      is not written here.
#
# Primary data (this directory):
#   Sample-to-participant crosswalks are not used here. Each PhIP-Seq table is
#   one row per participant. Covariate columns sit beside the peptide columns.
#   - `baseline_phipseq.csv`: `ID`, `sample_id`, `plate_id`, `Age`, `Sex`
#     (`F`/`M`), `Stool_Freq`, `stool_type`, `IBS_Subtype`
#     (`no IBS`, `IBS-C`, `IBS-D`, `IBS-M`), `CMV_prediction`, and peptide
#     columns (`agilent_*`, `twist_*`). Enrichment values are `> 0` (bound),
#     `-1` (contamination), or `0`/`NA` (unbound).
#   - `followup_phipseq.csv`: the same peptide columns plus `ID`,
#     `sample_id`, `plate_id`, `Age`, `Sex`, `IBS_Subtype`, `CMV_prediction`,
#     and `followup_time` (years).
#   - `peptide_library.csv`: peptide sequences used to drop identical
#     amino-acid sequences (`V1` = peptide ID, `aa_seq`).
#   - `pi_ibs_annotation.csv`: post-infectious IBS markers (`peptide`,
#     `taxa`, `Is_PI_IBS_Marker`, `PI_IBS_Category`).
#
# Generated Outputs:
#   - `pi_ibs_markers.csv`: Filtered PI-IBS marker lookup table.
#   - `1_phenos_ibs.rds`: QC'd baseline PhIP-Seq and phenotype data.table.
#   - `1_phenos_fup.rds`: QC'd follow-up PhIP-Seq and phenotype data.table.
#   - `fig1_qc_panels.rds`: Figure 1 panels (library treemap, baseline
#     richness raincloud, and prevalence histogram) used in 2_Exploratory.R.
#   - `extended_figure_1.pdf`: Follow-up QC summary figure.
# ==============================================================================

library(data.table)
library(readr)
library(stringr)
library(ggplot2)
library(ggdist)
library(paletteer)
library(scales)
library(patchwork)
library(treemapify)
library(Rcpp)

source("0_utils.R")

# 1. Load Primary Tables -------------------------------------------------------

raw_dat <- fread("baseline_phipseq.csv")
raw_fup <- fread("followup_phipseq.csv")
anno <- fread("peptide_library.csv")
setnames(anno, "V1", "peptide")

ibs_levels <- c("no IBS", "IBS-C", "IBS-D", "IBS-M")
pep_cols <- names(raw_dat)[grepl("^agi|^twi", names(raw_dat))]

raw_dat[, let(
  IBS_Subtype = factor(IBS_Subtype, levels = ibs_levels),
  Sex = factor(Sex, levels = c("F", "M")),
  CMV_prediction = factor(CMV_prediction)
)]
raw_dat <- raw_dat[!is.na(IBS_Subtype)]
raw_dat[,
  IBS_Status := factor(
    fifelse(IBS_Subtype == "no IBS", "no IBS", "IBS"),
    levels = c("no IBS", "IBS")
  )
]

raw_fup[, let(
  IBS_Subtype = factor(IBS_Subtype, levels = ibs_levels),
  Sex = factor(Sex, levels = c("F", "M")),
  CMV_prediction = factor(CMV_prediction)
)]
raw_fup <- raw_fup[!is.na(IBS_Subtype)]

# Label follow-up trajectories from the subtype already present at each timepoint
base_pheno <- raw_dat[, .(
  ID,
  IBS_Subtype_base = IBS_Subtype,
  IBS_Status_base = IBS_Status
)]
raw_fup <- base_pheno[raw_fup, on = .(ID), nomatch = 0]
raw_fup[,
  IBS_Status := factor(
    fifelse(IBS_Subtype == "no IBS", "no IBS", "IBS"),
    levels = c("no IBS", "IBS")
  )
]
raw_fup[,
  follow_up_grp := fcase(
    IBS_Subtype == IBS_Subtype_base                                              , "stable"       ,
    IBS_Status_base == "no IBS" & IBS_Subtype %in% c("IBS-C", "IBS-M", "IBS-D")  , "new_onset"    ,
    IBS_Subtype_base %in% c("IBS-C", "IBS-M", "IBS-D") & IBS_Subtype == "no IBS" , "improvement"  ,
    default = "class_switch"
  )
]
raw_fup[,
  follow_up_grp_ibs := fcase(
    IBS_Status == IBS_Status_base                     , "stable"      ,
    IBS_Status_base == "no IBS" & IBS_Status == "IBS" , "new_onset"   ,
    IBS_Status_base == "IBS" & IBS_Status == "no IBS" , "improvement"
  )
]


# 2. Plot Antigen Library Description ------------------------------------------

oligo_dat <- data.table(
  group = c(
    "Bacterial genes",
    "Bacterial strains",
    "Phages",
    "Controls",
    "Probiotic bacteria",
    "Pathogenic bacteria",
    "Immune epitopes",
    "Allergens",
    "Virulence factors",
    "Antibody-coated bacteria"
  ),
  count = c(
    122551,
    24510,
    40000,
    11525,
    14700,
    24500,
    30000,
    30000,
    24164,
    22050
  )
)

oligo_dat[, perc := round(count / sum(count) * 100)]
oligo_dat[,
  label := sprintf(
    "<b>%s</b><br><span style='font-size:8pt'><i>%s&nbsp;peptides&nbsp;(~%s%%)</i></span>",
    group,
    comma(count),
    perc
  )
]

layout <- treemapify(oligo_dat, area = "count", layout = "squarified")
layout[, let(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2)]

oligo_info <- layout |>
  ggplot() +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group),
    show.legend = FALSE,
    color = "grey30",
    alpha = .75
  ) +
  ggtext::geom_richtext(
    aes(x, y, label = label),
    colour = "black",
    label.color = NA,
    size = 3.5,
    fill = NA,
    vjust = .5,
    hjust = .5
  ) +
  theme_void() +
  scale_fill_paletteer_d(palette = "MoMAColors::Klein") +
  theme(plot.background = element_rect(fill = "#fbf9f4", colour = NA))


# 3. Binarize PhIP-Seq Enrichment Matrix ---------------------------------------
# Values > 0 are converted to 1 (bound); -1 indicates contamination and is set
# to NA; missing values are set to 0 prior to conversion.

cppFunction(
  "
void optimize_matrix(NumericMatrix x) {
  int n = x.size();
  for(int i = 0; i < n; ++i) {
    if (x[i] == -1) {
      x[i] = NA_REAL;
    } else if (x[i] > 0) {
      x[i] = 1;
    }
  }
}
"
)

convert_to_binary <- function(dt, cols) {
  setnafill(dt, cols = cols, fill = 0L)

  mat <- dt[, ..cols] |> as.matrix()
  optimize_matrix(mat)

  mat <- as.data.frame(mat)
  setDT(mat)

  dt <- dt[, !..pep_cols]
  cbind(dt, mat)
}

raw_fup <- convert_to_binary(raw_fup, pep_cols)
gc()
raw_dat <- convert_to_binary(raw_dat, pep_cols)
gc()

id_order <- raw_dat$ID
raw_fup <- raw_fup[match(id_order, raw_fup$ID)][!is.na(ID)]


# 4. Sample-Level & Peptide-Level Quality Control ------------------------------

# Filter out low-quality samples with < 200 uniquely bound peptides
sample_filter <- function(dt, cols = pep_cols) {
  rs <- rowSums(as.matrix(dt[, ..cols]), na.rm = TRUE)
  dt[, n_distinct_peptides := rs]
  dt[n_distinct_peptides >= 200]
}

raw_dat <- sample_filter(raw_dat)
raw_fup <- sample_filter(raw_fup)

# Compare number of bound peptides across IBS subtypes (Kruskal-Wallis & Dunn's test)
raw_dat[, kruskal.test(n_distinct_peptides ~ IBS_Subtype)]
raw_dat[, dunn.test::dunn.test(n_distinct_peptides, IBS_Subtype, method = "BH")]
raw_fup[, kruskal.test(n_distinct_peptides ~ IBS_Subtype)]
raw_fup[, dunn.test::dunn.test(n_distinct_peptides, IBS_Subtype, method = "BH")]

# Medians and IQRs per subtype and globally
raw_fup[,
  by = IBS_Subtype,
  print(setNames(.(round(quantile(n_distinct_peptides))), .BY[[1]]))
]
raw_fup[, round(quantile(n_distinct_peptides))]

raw_dat[,
  by = IBS_Subtype,
  print(setNames(.(round(quantile(n_distinct_peptides))), .BY[[1]]))
]
raw_dat[, round(quantile(n_distinct_peptides))]

raw_dat[, median(n_distinct_peptides)]
raw_dat[, quantile(as.integer(n_distinct_peptides))[c(2, 4)]]
raw_fup[, median(n_distinct_peptides)]
raw_fup[, quantile(as.integer(n_distinct_peptides))[c(2, 4)]]

raw_fup[, IBS_Subtype_base := relevel(IBS_Subtype_base, "no IBS")]

# Raincloud plots of bound antigen richness per subtype (baseline and follow-up)
raincloud <- raw_dat |>
  ggplot(aes(n_distinct_peptides, IBS_Subtype, fill = IBS_Subtype)) +
  stat_halfeye(
    aes(slab_colour = after_scale(fill)),
    adjust = .6,
    slab_alpha = .4,
    .width = 0,
    height = .8,
    point_colour = NA,
    justification = -.15,
    show.legend = FALSE
  ) +
  geom_boxplot(
    color = "grey30",
    width = 0.12,
    outlier.color = NA,
    alpha = 0.4
  ) +
  geom_vline(
    xintercept = round(median(raw_dat$n_distinct_peptides)),
    color = "grey30",
    linetype = "dashed",
    alpha = .6
  ) +
  annotate(
    "text",
    x = round(median(raw_dat$n_distinct_peptides)),
    y = max(as.numeric(raw_dat$IBS_Subtype)) + .85,
    label = expression(italic("Global Median")),
    size = 3,
    hjust = -.1,
    color = "grey30"
  ) +
  scale_fill_paletteer_d(palette = "MoMAColors::Klein") +
  scale_color_paletteer_d(palette = "MoMAColors::Klein") +
  stat_dots(
    aes(y = as.numeric(IBS_Subtype) - .125, col = IBS_Subtype),
    binwidth = 5
  ) +
  theme_my(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(margin = margin(r = 1)),
    axis.title.y = element_blank()
  ) +
  guides(
    fill = guide_legend(
      title = "Study group",
      nrow = 1,
      byrow = TRUE,
      override.aes = list(colour = NA)
    ),
    colour = guide_legend(title = "Study group")
  ) +
  scale_y_discrete(expand = expansion(add = c(0.25, 0))) +
  scale_x_continuous(expand = expansion(mult = c(.01, 0.01))) +
  coord_cartesian(clip = "off") +
  labs(x = "Distinct significantly-bound antigens (per individual)")

raw_fup[,
  IBS_Subtype := factor(
    IBS_Subtype,
    levels = c("IBS-C", "IBS-D", "IBS-M", "no IBS")
  )
]
raincloud_fup <- raw_fup |>
  ggplot(aes(n_distinct_peptides, IBS_Subtype, fill = IBS_Subtype)) +
  stat_halfeye(
    aes(slab_colour = after_scale(fill)),
    adjust = .6,
    slab_alpha = .4,
    .width = 0,
    height = .8,
    point_colour = NA,
    justification = -.15,
    show.legend = FALSE
  ) +
  geom_boxplot(
    color = "grey30",
    width = 0.12,
    outlier.color = NA,
    alpha = 0.4
  ) +
  geom_vline(
    xintercept = round(median(raw_fup$n_distinct_peptides)),
    color = "grey30",
    linetype = "dashed",
    alpha = .6
  ) +
  annotate(
    "text",
    x = round(median(raw_fup$n_distinct_peptides)),
    y = max(as.numeric(raw_dat$IBS_Subtype)) + .85,
    label = expression(italic("Global Median")),
    size = 3,
    hjust = -.1,
    color = "grey30"
  ) +
  scale_fill_paletteer_d(palette = "MoMAColors::Klein") +
  scale_color_paletteer_d(palette = "MoMAColors::Klein") +
  stat_dots(
    aes(y = as.numeric(IBS_Subtype) - .125, col = IBS_Subtype),
    binwidth = 5
  ) +
  theme_my(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(margin = margin(r = 1)),
    axis.title.y = element_blank()
  ) +
  guides(
    fill = guide_legend(
      title = "Study group",
      nrow = 1,
      byrow = TRUE,
      override.aes = list(colour = NA)
    ),
    colour = guide_legend(title = "Study group")
  ) +
  scale_y_discrete(expand = expansion(add = c(0.25, 0))) +
  scale_x_continuous(expand = expansion(mult = c(.01, 0.01))) +
  coord_cartesian(clip = "off") +
  labs(x = "Distinct significantly-bound antigens (per individual)")


# 5. Peptide Prevalence & Sequence Similarity Filtering ------------------------

cs <- colSums(as.matrix(raw_dat[, ..pep_cols]), na.rm = TRUE)
cs <- as.data.table(cs, keep.rownames = TRUE)

cs_fup <- colSums(as.matrix(raw_fup[, ..pep_cols]), na.rm = TRUE)
cs_fup <- as.data.table(cs_fup, keep.rownames = TRUE)

# Calculate overall seroprevalence per peptide
prevalence <- calc_prev(raw_dat, cols = pep_cols)
prev_fup <- calc_prev(raw_fup, cols = pep_cols)

cs[, prev := prevalence]
cs[, pass := between(prev, .05, .95)]

cs_fup[, prev := prev_fup]
cs_fup[, pass := between(prev, .05, .95)]

# Retain peptides bound in >= 5% and < 95% of individuals
prev_filter <- function(dt, prev) {
  rm_cols <- names(prev[prev < 0.05 | prev >= 0.95])
  dt[, !..rm_cols]
}
raw_dat <- prev_filter(raw_dat, prevalence)
raw_fup <- prev_filter(raw_fup, prev_fup)
gc()

# Deduplicate peptides sharing identical amino acid sequences, keeping the
# peptide with the highest seroprevalence
anno <- anno[peptide %in% union(names(raw_dat), names(raw_fup))]

simil_filter <- function(dt, prev, annotaton = anno) {
  prev <- as.data.table(prev, keep.rownames = TRUE)
  setnames(prev, c("peptide", "preva"))
  annotaton <- annotaton[, .(peptide, aa_seq = str_extract(aa_seq, "^\\w+"))]
  prev <- prev[annotaton, on = .(peptide), nomatch = 0]
  prev <- prev[, .SD[which.max(preva)], by = aa_seq]

  seq_omit <- names(dt)[grepl("^agi|^twi", names(dt))]
  seq_omit <- setdiff(seq_omit, prev$peptide)
  dt[, !..seq_omit]
}

raw_dat <- simil_filter(raw_dat, prevalence)
raw_fup <- simil_filter(raw_fup, prev_fup)

raw_dat[, IBS_Subtype := relevel(IBS_Subtype, "no IBS")]
raw_fup[, IBS_Subtype := relevel(IBS_Subtype, "no IBS")]


# 6. Plot Seroprevalence Filter Distributions & Extended Figure 1 --------------

fin_pass <- raw_dat[, length(.SD), .SDcols = patterns("^agi|^twi")]
fin_pass_fup <- raw_fup[, length(.SD), .SDcols = patterns("^agi|^twi")]

prev_filter_plot <-
  as.data.table(cs) |>
  ggplot() +
  geom_histogram(
    aes(cs, fill = pass),
    binwidth = 9,
    alpha = .9,
    position = "identity",
    bins = 70
  ) +
  annotate(
    "segment",
    x = cs[pass == TRUE, min(cs) + 20],
    xend = cs[pass == TRUE, max(cs) - 20],
    y = 800,
    yend = 800,
    colour = "grey30",
    linewidth = .6,
    linetype = "solid",
    arrow = arrow(
      ends = "both",
      type = "closed",
      length = unit(4, "pt")
    )
  ) +
  annotate(
    "text",
    x = cs[pass == TRUE, (min(cs) + max(cs)) / 2],
    label = paste0(
      "italic('",
      fin_pass,
      " antibody-bound peptides in 5-95% of individuals')"
    ),
    parse = TRUE,
    y = 800,
    vjust = -.3,
    size = 3.5
  ) +
  scale_fill_paletteer_d(palette = "MoMAColors::Klein") +
  theme_my(legend.position = "bottom", panel.grid.major.y = element_blank()) +
  scale_y_log10(
    expand = expansion(mult = c(0, 0.01)),
    breaks = breaks_log(base = 10, n = 6),
    labels = label_log()
  ) +
  scale_x_continuous(expand = expansion(mult = 0)) +
  guides(fill = guide_legend(title = "Prevalence filter pass")) +
  labs(
    x = "Number of individuals",
    y = "Number of seroprevalent antibody-bound peptides"
  ) +
  coord_cartesian(clip = "off")

prev_filter_plot_fup <-
  as.data.table(cs_fup) |>
  ggplot() +
  geom_histogram(
    aes(cs_fup, fill = pass),
    binwidth = 9,
    alpha = .9,
    position = "identity",
    bins = 70
  ) +
  annotate(
    "segment",
    x = cs_fup[pass == TRUE, min(cs_fup) + 20],
    xend = cs_fup[pass == TRUE, max(cs_fup) - 20],
    y = 800,
    yend = 800,
    colour = "grey30",
    linewidth = .6,
    linetype = "solid",
    arrow = arrow(
      ends = "both",
      type = "closed",
      length = unit(4, "pt")
    )
  ) +
  annotate(
    "text",
    x = cs_fup[pass == TRUE, (min(cs_fup) + max(cs_fup)) / 2],
    label = paste0(
      "italic('",
      fin_pass_fup,
      " antibody-bound peptides in 5-95% of individuals')"
    ),
    parse = TRUE,
    y = 800,
    vjust = -.3,
    size = 3.5
  ) +
  scale_fill_paletteer_d(palette = "MoMAColors::Klein") +
  theme_my(legend.position = "bottom", panel.grid.major.y = element_blank()) +
  scale_y_log10(
    expand = expansion(mult = c(0, 0.01)),
    breaks = breaks_log(base = 10, n = 6),
    labels = label_log()
  ) +
  scale_x_continuous(expand = expansion(mult = 0)) +
  guides(fill = guide_legend(title = "Prevalence filter pass")) +
  labs(
    x = "Number of individuals",
    y = "Number of seroprevalent antibody-bound peptides"
  ) +
  coord_cartesian(clip = "off")

raincloud_fup +
  prev_filter_plot_fup +
  plot_layout(widths = c(0.6, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "#fbf9f4")
  )

ggsave("extended_figure_1.pdf", width = 8)


# 7. Format PI-IBS Markers & Export QC Outputs ---------------------------------

pi_ibs_markers <- fread("pi_ibs_annotation.csv")[
  Is_PI_IBS_Marker == "Yes",
  .(peptide, taxa, PI_IBS_Category)
]
pi_ibs_markers[, source_database := paste(PI_IBS_Category, "PI-IBS marker")]
pi_ibs_markers[, PI_IBS_Category := NULL]
pi_ibs_markers <- unique(pi_ibs_markers)

fwrite(pi_ibs_markers, "pi_ibs_markers.csv")

write_rds(raw_dat, "1_phenos_ibs.rds")
write_rds(raw_fup, "1_phenos_fup.rds")

write_rds(
  list(
    oligo_info = oligo_info,
    raincloud = raincloud,
    prev_filter_plot = prev_filter_plot
  ),
  "fig1_qc_panels.rds"
)
