# ==============================================================================
# 0_utils.R
# Shared utilities, plotting theme, and helper functions
# ==============================================================================
#
# Purpose:
#   Provides the custom ggplot2 theme (`theme_my`) and the shared
#   seroprevalence calculation (`calc_prev`) used across scripts 1 through 6.
#
# Usage:
#   source("0_utils.R") from the same directory as the analysis scripts.
# ==============================================================================

library(data.table)
library(ggplot2)
library(forcats)

# Custom ggplot2 Theme ---------------------------------------------------------
theme_my <- function(...) {
  theme_bw() +
    theme(
      axis.ticks = element_line(colour = "grey50"),
      axis.line = element_line(colour = "grey50"),
      axis.line.x.top = element_line(colour = "black"),
      panel.grid = element_line(color = "#b4aea9"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linetype = "dashed"),
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      panel.border = element_rect(color = "grey50", linewidth = .5),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12),
      legend.background = element_rect(fill = "#fbf9f4"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        hjust = 0,
        vjust = -.3,
        margin = margin(b = 4)
      ),
      legend.frame = element_rect(colour = "grey30"),
      legend.ticks = element_line(colour = "grey30"),
      legend.title = element_text(face = "bold")
    ) +
    theme(...)
}

# Shared Prevalence Helper -----------------------------------------------------
# Computes per-peptide seroprevalence (mean of binary 0/1 binding calls).
# - If `group_var = NULL`: returns a named numeric vector of overall prevalence.
# - If `group_var` is provided and `wide = TRUE`: returns a wide data.table
#   with one column per group level (used in 4_Prevalence_and_Enrichment.R).
# - If `group_var` is provided and `wide = FALSE`: returns a long data.table
#   with columns `<group_var>`, `V1`, and `peptide` (used in 6_Visualize_ML.R).
calc_prev <- function(dt, group_var = NULL, cols = NULL, wide = TRUE) {
  if (is.null(cols)) {
    cols <- grep("^agi|^twi", names(dt), value = TRUE)
  }

  if (is.null(group_var)) {
    return(colMeans(as.matrix(dt[, ..cols]), na.rm = TRUE))
  }

  res <- dt[,
    keyby = group_var,
    .(colMeans(.SD, na.rm = TRUE), peptide = fct_inorder(cols)),
    .SDcols = cols
  ]

  if (wide) {
    return(res[, dcast(.SD, peptide ~ ..., value.var = "V1")])
  }

  res
}
