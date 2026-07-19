# Combined cell-type-specific two-sample MR analysis
# Exposures: cell-type-specific eQTLs
# Outcomes: depression and amygdala phenotypes

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(dplyr)
  library(stringr)
  library(ieugwasr)
  library(plinkbinr)
})

# -----------------------------------------------------------------------------
# 1. File paths and analysis parameters
# -----------------------------------------------------------------------------

gene_list_file <- "test/MR/cellQTL/gene.txt"
cell_type_file <- "test/MR/cell_types.txt"
amygdala_outcome_list_file <- "test/MR/amygdala_outcomes.txt"

mr_base_dir <- "test/MR"
amygdala_outcome_dir <- "test/MR/AMY"

depression_outcome_file <- paste0(
  "test/cFDR/traitfolder/",
  "depression_summary_statistics.ma"
)

plink_ref_bfile <- "reference/EUR/EUR"

exposure_p_threshold <- 1e-5
clump_kb <- 10000
clump_r2 <- 0.001
exposure_sample_size <- 192

# -----------------------------------------------------------------------------
# 2. Helper functions
# One IV: Wald ratio; two or more IVs: IVW only
# -----------------------------------------------------------------------------

calculate_instrument_strength <- function(dat, sample_size) {

  maf <- pmin(dat$eaf.exposure, 1 - dat$eaf.exposure)

  r2_numerator <- 2 * maf * (1 - maf) * dat$beta.exposure^2

  r2_denominator <- 2 * maf * (1 - maf) *
    (
      dat$beta.exposure^2 +
        sample_size * dat$se.exposure^2
    )

  r2 <- ifelse(
    r2_denominator == 0,
    0,
    r2_numerator / r2_denominator
  )

  f_statistic <- ifelse(
    r2 >= 1,
    Inf,
    r2 * (sample_size - 2) / (1 - r2)
  )

  data.frame(
    SNP = dat$SNP,
    F_value = f_statistic,
    stringsAsFactors = FALSE
  )
}


run_mr_for_outcome <- function(
    exposure_all,
    outcome_dat,
    outcome_name,
    cell_type,
    exposure_sample_size) {

  mr_results <- data.frame()
  f_results <- data.frame()
  heterogeneity_results <- data.frame()

  genes <- unique(exposure_all$exposure)

  for (gene in genes) {

    exposure_dat <- exposure_all %>%
      filter(exposure == gene) %>%
      as.data.frame()

    harmonised_dat <- tryCatch(
      {
        TwoSampleMR::harmonise_data(
          exposure_dat = exposure_dat,
          outcome_dat = outcome_dat,
          action = 2
        )
      },
      error = function(e) NULL
    )

    if (is.null(harmonised_dat) || nrow(harmonised_dat) == 0) {
      next
    }

    if ("mr_keep" %in% names(harmonised_dat)) {
      harmonised_dat <- harmonised_dat %>%
        filter(mr_keep)
    }

    number_of_snps <- nrow(harmonised_dat)

    if (number_of_snps == 0) {
      next
    }

    f_dat <- calculate_instrument_strength(
      harmonised_dat,
      exposure_sample_size
    )

    f_dat$cell_type <- cell_type
    f_dat$gene <- gene
    f_dat$outcome <- outcome_name

    f_results <- bind_rows(f_results, f_dat)

    mr_methods <- if (number_of_snps == 1) {
      "mr_wald_ratio"
    } else {
      "mr_ivw"
    }

    mr_result <- tryCatch(
      {
        TwoSampleMR::mr(
          harmonised_dat,
          method_list = mr_methods
        )
      },
      error = function(e) NULL
    )

    if (!is.null(mr_result) && nrow(mr_result) > 0) {
      mr_result$cell_type <- cell_type
      mr_result$gene <- gene
      mr_result$snp_count <- number_of_snps
      mr_result$outcome_name <- outcome_name

      mr_results <- bind_rows(
        mr_results,
        mr_result
      )
    }

    if (number_of_snps >= 2) {

      heterogeneity_result <- tryCatch(
        {
          TwoSampleMR::mr_heterogeneity(
            harmonised_dat,
            method_list = "mr_ivw"
          )
        },
        error = function(e) NULL
      )

      if (!is.null(heterogeneity_result) &&
          nrow(heterogeneity_result) > 0) {

        heterogeneity_result$cell_type <- cell_type
        heterogeneity_result$gene <- gene
        heterogeneity_result$outcome_name <- outcome_name

        heterogeneity_results <- bind_rows(
          heterogeneity_results,
          heterogeneity_result
        )
      }
    }
  }

  list(
    mr = mr_results,
    f = f_results,
    heterogeneity = heterogeneity_results
  )
}


write_result_if_available <- function(dat, output_file) {
  if (!is.null(dat) && nrow(dat) > 0) {
    write.csv(
      dat,
      output_file,
      row.names = FALSE
    )
  }
}

# -----------------------------------------------------------------------------
# 3. Read shared input files
# -----------------------------------------------------------------------------

gene_list <- fread(
  gene_list_file,
  header = FALSE,
  encoding = "UTF-8"
)[[1]]

if (length(gene_list) == 0) {
  stop("The gene list is empty.")
}

cell_types <- fread(
  cell_type_file,
  header = FALSE,
  encoding = "UTF-8"
)[[1]]

cell_types <- unique(
  cell_types[!is.na(cell_types) & cell_types != ""]
)

if (length(cell_types) == 0) {
  stop("The cell-type list is empty.")
}

amygdala_outcomes <- fread(
  amygdala_outcome_list_file,
  header = FALSE,
  encoding = "UTF-8"
)[[1]]

amygdala_outcomes <- unique(
  amygdala_outcomes[
    !is.na(amygdala_outcomes) &
      amygdala_outcomes != ""
  ]
)

if (length(amygdala_outcomes) == 0) {
  stop("The amygdala outcome list is empty.")
}

# -----------------------------------------------------------------------------
# 4. Format the depression outcome
# -----------------------------------------------------------------------------

raw_depression_outcome <- as.data.frame(
  fread(depression_outcome_file)
)

raw_depression_outcome$trait_id <- "depression"

depression_outcome_dat <- TwoSampleMR::format_data(
  dat = raw_depression_outcome,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FRQ",
  pval_col = "P",
  samplesize_col = "N",
  chr_col = "CHR",
  pos_col = "BP",
  phenotype_col = "trait_id"
)

if (nrow(depression_outcome_dat) == 0) {
  stop("The formatted depression outcome dataset is empty.")
}

# -----------------------------------------------------------------------------
# 5. Load amygdala outcome datasets
# -----------------------------------------------------------------------------

amygdala_outcome_data <- list()

for (outcome_id in amygdala_outcomes) {

  outcome_file <- file.path(
    amygdala_outcome_dir,
    paste0(outcome_id, ".rds")
  )

  if (!file.exists(outcome_file)) {
    warning("Amygdala outcome file not found: ", outcome_file)
    next
  }

  outcome_dat <- tryCatch(
    readRDS(outcome_file),
    error = function(e) NULL
  )

  if (!is.null(outcome_dat) && nrow(outcome_dat) > 0) {
    amygdala_outcome_data[[as.character(outcome_id)]] <- outcome_dat
  }
}

if (length(amygdala_outcome_data) == 0) {
  stop("No valid amygdala outcome datasets were loaded.")
}

# -----------------------------------------------------------------------------
# 6. Run exposure preparation and MR analyses for each cell type
# -----------------------------------------------------------------------------

for (cell_type in cell_types) {

  message("\nProcessing cell type: ", cell_type)

  cell_eqtl_dir <- file.path(
    mr_base_dir,
    cell_type
  )

  cell_eqtl_file <- file.path(
    cell_eqtl_dir,
    paste0(cell_type, ".with_eaf.txt")
  )

  filtered_eqtl_file <- file.path(
    cell_eqtl_dir,
    paste0(cell_type, ".filteredgene.txt")
  )

  exposure_file <- file.path(
    cell_eqtl_dir,
    paste0(cell_type, ".exposure_all_genes.csv")
  )

  depression_output_dir <- file.path(
    cell_eqtl_dir,
    "depression_out"
  )

  amygdala_output_dir <- file.path(
    cell_eqtl_dir,
    "amygdala_out"
  )

  dir.create(
    depression_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    amygdala_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!file.exists(cell_eqtl_file)) {
    warning("Input eQTL file not found: ", cell_eqtl_file)
    next
  }

  # ---------------------------------------------------------------------------
  # 6.1 Extract genes of interest
  # ---------------------------------------------------------------------------

  cell_eqtl <- fread(
    cell_eqtl_file,
    sep = "\t",
    stringsAsFactors = FALSE,
    encoding = "UTF-8"
  )

  if (!"Gene_id" %in% names(cell_eqtl)) {
    warning("Gene_id column not found in: ", cell_eqtl_file)
    next
  }

  cell_eqtl[, Ensembl_ID := str_extract(Gene_id, "ENSG[0-9]+")]

  filtered_eqtl <- cell_eqtl[
    Ensembl_ID %in% gene_list
  ]

  if (nrow(filtered_eqtl) == 0) {
    warning("No target genes were found for: ", cell_type)
    next
  }

  filtered_eqtl[, gene_order := match(Ensembl_ID, gene_list)]
  setorder(filtered_eqtl, gene_order)
  filtered_eqtl[, gene_order := NULL]

  fwrite(
    filtered_eqtl,
    filtered_eqtl_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "",
    encoding = "UTF-8"
  )

  # ---------------------------------------------------------------------------
  # 6.2 Select and LD-clump eQTL instruments
  # ---------------------------------------------------------------------------

  cell_eqtl_df <- as.data.frame(filtered_eqtl)
  gene_ids <- unique(cell_eqtl_df$Gene_id)
  exposure_list <- list()

  for (gene_id in gene_ids) {

    dat <- subset(
      cell_eqtl_df,
      Gene_id == gene_id &
        P < exposure_p_threshold
    )

    if (nrow(dat) == 0) {
      next
    }

    dat$trait_id <- gene_id

    clumped <- tryCatch(
      {
        ieugwasr::ld_clump(
          dplyr::tibble(
            rsid = dat$SNP,
            pval = dat$P,
            id = dat$trait_id
          ),
          plink_bin = plinkbinr::get_plink_exe(),
          bfile = plink_ref_bfile,
          clump_kb = clump_kb,
          clump_r2 = clump_r2
        )
      },
      error = function(e) NULL
    )

    if (is.null(clumped) || nrow(clumped) == 0) {
      next
    }

    dat <- dat[
      dat$SNP %in% clumped$rsid,
      ,
      drop = FALSE
    ]

    if (nrow(dat) == 0) {
      next
    }

    p_for_se <- pmax(
      dat$P,
      .Machine$double.xmin
    )

    z_score <- qnorm(
      p_for_se / 2,
      lower.tail = FALSE
    )

    dat$SE <- abs(dat$Beta) / z_score

    dat <- dat[
      is.finite(dat$SE) &
        !is.na(dat$SE) &
        dat$SE > 0,
      ,
      drop = FALSE
    ]

    if (nrow(dat) == 0) {
      next
    }

    exposure_dat <- tryCatch(
      {
        TwoSampleMR::format_data(
          dat = dat,
          type = "exposure",
          snp_col = "SNP",
          beta_col = "Beta",
          se_col = "SE",
          eaf_col = "eaf",
          effect_allele_col = "effect_allele",
          other_allele_col = "other_allele",
          pval_col = "P",
          chr_col = "CHR",
          pos_col = "BP",
          phenotype_col = "trait_id"
        )
      },
      error = function(e) NULL
    )

    if (!is.null(exposure_dat) &&
        nrow(exposure_dat) > 0) {
      exposure_list[[gene_id]] <- exposure_dat
    }
  }

  if (length(exposure_list) == 0) {
    warning("No valid exposure instruments for: ", cell_type)
    next
  }

  exposure_all <- bind_rows(exposure_list)

  write.csv(
    exposure_all,
    exposure_file,
    row.names = FALSE
  )

  # ---------------------------------------------------------------------------
  # 6.3 Depression as the outcome
  # ---------------------------------------------------------------------------

  depression_results <- run_mr_for_outcome(
    exposure_all = exposure_all,
    outcome_dat = depression_outcome_dat,
    outcome_name = "depression",
    cell_type = cell_type,
    exposure_sample_size = exposure_sample_size
  )

  write_result_if_available(
    depression_results$mr,
    file.path(
      depression_output_dir,
      paste0(cell_type, "_depression_MR_results.csv")
    )
  )

  write_result_if_available(
    depression_results$f,
    file.path(
      depression_output_dir,
      paste0(cell_type, "_depression_F_statistics.csv")
    )
  )

  write_result_if_available(
    depression_results$heterogeneity,
    file.path(
      depression_output_dir,
      paste0(cell_type, "_depression_heterogeneity.csv")
    )
  )

  # ---------------------------------------------------------------------------
  # 6.4 Amygdala phenotypes as outcomes
  # ---------------------------------------------------------------------------

  all_amygdala_mr <- data.frame()
  all_amygdala_f <- data.frame()
  all_amygdala_heterogeneity <- data.frame()

  for (outcome_name in names(amygdala_outcome_data)) {

    amygdala_results <- run_mr_for_outcome(
      exposure_all = exposure_all,
      outcome_dat = amygdala_outcome_data[[outcome_name]],
      outcome_name = outcome_name,
      cell_type = cell_type,
      exposure_sample_size = exposure_sample_size
    )

    all_amygdala_mr <- bind_rows(
      all_amygdala_mr,
      amygdala_results$mr
    )

    all_amygdala_f <- bind_rows(
      all_amygdala_f,
      amygdala_results$f
    )

    all_amygdala_heterogeneity <- bind_rows(
      all_amygdala_heterogeneity,
      amygdala_results$heterogeneity
    )

  }

  write_result_if_available(
    all_amygdala_mr,
    file.path(
      amygdala_output_dir,
      paste0(cell_type, "_amygdala_MR_results.csv")
    )
  )

  write_result_if_available(
    all_amygdala_f,
    file.path(
      amygdala_output_dir,
      paste0(cell_type, "_amygdala_F_statistics.csv")
    )
  )

  write_result_if_available(
    all_amygdala_heterogeneity,
    file.path(
      amygdala_output_dir,
      paste0(cell_type, "_amygdala_heterogeneity.csv")
    )
  )


  message("Completed cell type: ", cell_type)
}

message("\nAll depression and amygdala MR analyses are complete.")
