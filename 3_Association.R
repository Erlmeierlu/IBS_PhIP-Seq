# ==============================================================================
# 3_Association.R
# Per-peptide Firth logistic regression association analyses
# ==============================================================================
#
# Purpose:
#   Tests per-peptide antibody binding associations using Firth's bias-reduced
#   penalized-likelihood logistic regression (`logistf`), adjusting for Age,
#   Sex, and sequencing plate (`plate_id`), across:
#     1. IBS Subtypes (`IBS-C`, `IBS-D`, `IBS-M` vs `no IBS` and pairwise
#        between-subtype contrasts via `emmeans`).
#     2. Baseline controls who developed IBS at follow-up vs consistent controls
#        (Extended Data Table 9).
#     3. Sensitivity analysis of IBS vs no IBS excluding IBS-M (`res_nod`).
#     4. Age group extremes (oldest vs youngest quartile, `res_age`).
#     5. Sex (Male vs Female, `res_sex`).
#     6. General IBS vs no IBS (`res_ibs`).
#
# Expected Inputs:
#   - `annotation.csv`: Peptide taxonomic and source database annotations.
#   - `1_phenos_ibs.rds`: QC'd baseline PhIP-Seq dataset (from 1_QC.R).
#   - `1_phenos_fup.rds`: QC'd follow-up PhIP-Seq dataset (from 1_QC.R).
#
# Generated Outputs:
#   - `Association_IBS_subtypes.csv`
#   - `Association_ibs.csv`
#   - `Association_age.csv`
#   - `Association_sex.csv`
#   - `Association_baseline_fup.csv`
#   - `Association_IBS_no_d.csv`
#   - `extended_table_9.xlsx`
#   - `pval_histo_subtypes_ctrl.pdf`
#   - `pval_histo_between_subtypes.pdf`
# ==============================================================================

library(data.table)
library(readr)
library(future.apply)
library(rlang)
library(broom)
library(emmeans)
library(brglm2)
library(logistf)
library(ggplot2)
library(writexl)

source("0_utils.R")

plan(multisession)

# 1. Load Data -----------------------------------------------------------------

annotation <- fread("annotation.csv")
phenos <- read_rds("1_phenos_ibs.rds")
fup <- read_rds("1_phenos_fup.rds")


# 2. Association Regression Helper ---------------------------------------------

run_association <- function(
  DT,
  peptide,
  method = c("firth", "glm", "lm"),
  covariates = covariates,
  filter_by = c("IBS_Status", "IBS_Subtype")
) {
  DT <- DT[!is.na(get(filter_by))] |> droplevels()

  filter_by <- match.arg(filter_by, names(DT))
  method <- match.arg(method, c("firth", "glm", "lm"))

  comp <- DT[, unique(get(filter_by))]

  if (method == "glm") {
    .fun <- function(formula, data) {
      glm(formula, data = data, family = binomial(), method = "brglmFit")
    }
  } else if (method == "lm") {
    .fun <- lm
  } else if (method == "firth") {
    .fun <- function(formula, data) logistf::logistf(formula, data = data)
  }

  form <- formula(paste0(
    peptide,
    " ~ ",
    paste0(c(filter_by, covariates), collapse = "+")
  ))

  res <- DT[, .fun(formula = form, .SD)]

  ncomp <- length(comp)
  if (!ncomp > 2) {
    if (method == "firth") {
      res <- data.table(
        term = res$term,
        SE = sqrt(diag(vcov(res))),
        df = Inf,
        p.value = res$prob,
        estimate = res$coef,
        covariates = paste(covariates, collapse = "+"),
        baseline = gsub(" ", "", levels(comp)[1]),
        peptide = peptide
      )
      res <- res[grepl(filter_by, term)]
      res[, term := gsub(filter_by, "", term)]

      return(res)
    } else {
      res <- broom::tidy(res)
      setDT(res)
      res <- res[grepl(filter_by, term)]
      res[, let(
        peptide = peptide,
        covariates = paste(covariates, collapse = "+"),
        baseline = gsub(" ", "", levels(comp)[1]),
        term = gsub(filter_by, "", term)
      )]
      return(res)
    }
  }

  # Multi-group comparison: compute all pairwise contrasts via emmeans
  EMM <- emmeans(res, filter_by, adjust = "none")
  res <- contrast(EMM, method = "revpairwise", adjust = "none") |>
    as.data.frame()
  setDT(res)
  res[, let(
    peptide = peptide,
    covariates = paste(covariates, collapse = "+"),
    contrast = gsub("\\(|\\)", "", contrast)
  )]
  res[, c("term", "baseline") := tstrsplit(contrast, " - ")]
  res[, baseline := gsub(" ", "", baseline)]
  setnames(
    res,
    c("z.ratio", "t.ratio"),
    c("statistic", "statistic"),
    skip_absent = TRUE
  )

  contrasts <- ncomp - 1
  times <- rep(1 / contrasts, contrasts)

  res_2 <- contrast(
    EMM,
    list("noIBS_vs_IBSany" = c(-1, times)),
    adjust = "none"
  ) |>
    as.data.frame()
  setDT(res_2)

  res_2[, let(
    peptide = peptide,
    covariates = paste(covariates, collapse = "+")
  )]
  res_2[, c("baseline", "term") := tstrsplit(contrast, "_vs_")]
  setnames(
    res_2,
    c("z.ratio", "t.ratio"),
    c("statistic", "statistic"),
    skip_absent = TRUE
  )

  res <- rbind(res, res_2)
  res
}


# 3. IBS Subtype Associations --------------------------------------------------

peptides <- grep(c("^agi|^twi"), names(phenos), value = TRUE)
covariates_pool <- c("Age", "Sex", "plate_id", "Stool_Freq")
covariates <- c("Age", "Sex", "plate_id")

res_sub <- future_lapply(
  peptides,
  run_association,
  DT = phenos,
  covariates = covariates,
  method = "firth",
  filter_by = "IBS_Subtype",
  future.seed = 833432,
  future.packages = "logistf"
) |>
  rbindlist()

res_any <- res_sub[term == "IBSany"] |> droplevels()
res_any[, FDR := p.adjust(p.value, method = "fdr")]
res_any <- annotation[res_any, on = .(peptide), nomatch = 0][order(p.value)]

res_sub <- res_sub[term != "IBSany"][order(p.value)] |> droplevels()
res_sub[, let(
  comparison = paste(term, baseline, sep = "_vs_"),
  FDR = p.adjust(p.value, method = "fdr")
)]
res_sub <- annotation[res_sub, on = .(peptide), nomatch = 0]

## 3a. P-Value Histograms ------------------------------------------------------
res_sub[baseline == "noIBS"] |>
  ggplot() +
  geom_histogram(aes(p.value, fill = comparison), col = "grey30", alpha = 0.8) +
  theme_my(plot.title.position = "plot") +
  scale_fill_brewer(palette = "Set2")
ggsave("pval_histo_subtypes_ctrl.pdf")

res_sub[baseline != "noIBS"] |>
  ggplot() +
  geom_histogram(aes(p.value, fill = comparison), col = "grey30", alpha = 0.8) +
  theme_my(plot.title.position = "plot") +
  scale_fill_brewer(palette = "Set2")
ggsave("pval_histo_between_subtypes.pdf")


# 4. Baseline Controls Developing IBS at Follow-Up vs Consistent Controls ------

ind <- fup[IBS_Subtype_base == "no IBS"] |> droplevels()
ind[, consistent := fifelse(IBS_Subtype == "no IBS", "consistent", "developed")]
ind[, consistent := factor(consistent)]

fup_base <- phenos[ind, on = .(ID)]

res_base <- future_lapply(
  peptides,
  run_association,
  DT = fup_base,
  covariates = covariates,
  method = "firth",
  filter_by = "consistent",
  future.seed = 833432,
  future.packages = "logistf"
) |>
  rbindlist()

res_base <- res_base[order(p.value)] |> droplevels()
res_base[, let(
  comparison = paste(term, baseline, sep = "_vs_"),
  FDR = p.adjust(p.value, method = "fdr")
)]
res_base <- annotation[res_base, on = .(peptide), nomatch = 0]

write_xlsx(
  res_base,
  path = "extended_table_9.xlsx"
)


# 5. Sensitivity Analysis: IBS vs Controls Excluding IBS-M ---------------------

ibs <- phenos[IBS_Subtype != "IBS-M"] |> droplevels()

res_nod <- future_lapply(
  peptides,
  run_association,
  DT = ibs,
  covariates = covariates,
  method = "firth",
  filter_by = "IBS_Status",
  future.seed = 833432,
  future.packages = "logistf"
) |>
  rbindlist()

res_nod <- res_nod[order(p.value)] |> droplevels()
res_nod[, let(
  comparison = paste(term, baseline, sep = "_vs_"),
  FDR = p.adjust(p.value, method = "fdr")
)]
res_nod <- annotation[res_nod, on = .(peptide), nomatch = 0]


# 6. Age Associations (Youngest vs Oldest Quartile) ----------------------------

age_dat <- phenos |> copy()
qs <- age_dat[, quantile(Age)]

age_dat[,
  age_grp := fcase(
    Age < qs[2] , "youngest" ,
    Age > qs[4] , "oldest"
  )
]
age_dat <- age_dat[!is.na(age_grp)]
age_dat[, age_grp := factor(age_grp, levels = c("youngest", "oldest"))]

res_age <- future_lapply(
  peptides,
  \(x) {
    form <- formula(paste0("age_grp ~ ", x))
    res <- tryCatch(
      logistf::logistf(form, data = age_dat),
      error = function(e) return(NULL)
    )
    if (is.null(res)) {
      return(NULL)
    }
    data.table(
      term = res$terms,
      SE = sqrt(diag(vcov(res))),
      df = Inf,
      p.value = res$prob,
      ci.low = res$ci.lower,
      ci.up = res$ci.upper,
      estimate = res$coef,
      comp = paste(rev(levels(age_dat$age_grp)), collapse = " vs ")
    )
  },
  future.seed = 833432,
  future.packages = "logistf"
) |>
  rbindlist()

res_age <- res_age[!term %like% "Inter"]
res_age[, odds_ratio := exp(estimate)]
res_age[, direction := fifelse(odds_ratio > 1, "old", "young")]
res_age[, FDR := p.adjust(p.value, method = "BH")]
setorder(res_age, FDR)
sig_age <- res_age[FDR < .05]
res_age <- annotation[res_age, on = .(peptide = term), nomatch = 0]


# 7. Sex Associations (Male vs Female) -----------------------------------------

sex_dat <- phenos |> copy()
res_sex <- future_lapply(
  peptides,
  \(x) {
    form <- formula(paste0("Sex ~ ", x))
    res <- logistf::logistf(form, data = sex_dat)
    data.table(
      term = res$terms,
      SE = sqrt(diag(vcov(res))),
      df = Inf,
      p.value = res$prob,
      ci.low = res$ci.lower,
      ci.up = res$ci.upper,
      estimate = res$coef,
      comp = paste(rev(levels(sex_dat$Sex)), collapse = " vs ")
    )
  },
  future.seed = 833432,
  future.packages = "logistf"
) |>
  rbindlist()

res_sex <- res_sex[!term %like% "Inter"]
res_sex[, odds_ratio := exp(estimate)]
res_sex[, direction := fifelse(odds_ratio > 1, "M", "F")]
res_sex[, FDR := p.adjust(p.value, method = "BH")]
setorder(res_sex, FDR)
sig_sex <- res_sex[FDR < .05]
sig_sex[, by = direction, .N]
res_sex <- annotation[res_sex, on = .(peptide = term), nomatch = 0]


# 8. General IBS vs no IBS Associations ----------------------------------------

res_ibs <- future_lapply(
  peptides,
  run_association,
  DT = phenos,
  method = "firth",
  covariates = covariates,
  filter_by = "IBS_Status",
  future.seed = 833432,
  future.packages = "logistf"
) |>
  rbindlist()

res_ibs <- res_ibs[order(p.value)]
res_ibs[, let(
  comparison = paste(term, baseline, sep = "_vs_"),
  FDR = p.adjust(p.value, method = "fdr")
)]
res_ibs <- annotation[res_ibs, on = .(peptide)]


# 9. Save Summary Tables -------------------------------------------------------

fwrite(res_sub, "Association_IBS_subtypes.csv")
fwrite(res_ibs, "Association_ibs.csv")
fwrite(res_age, "Association_age.csv")
fwrite(res_sex, "Association_sex.csv")
fwrite(res_base, "Association_baseline_fup.csv")
fwrite(res_nod, "Association_IBS_no_d.csv")
