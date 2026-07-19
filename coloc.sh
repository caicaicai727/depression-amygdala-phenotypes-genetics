# Cell-type-specific colocalization analysis
#
# Steps:
# 1. Extract the minimum-P eQTL SNP and its +/-100-kb locus for each cell-gene pair.
# 2. Merge each eQTL locus with the corresponding depression GWAS locus.
# 3. Harmonize alleles and remove unmatched or palindromic SNPs.
# 4. Run coloc.abf and combine the posterior-probability results.

  library(data.table)
  library(coloc)


# =============================================================================
# 1. File paths and analysis parameters
# =============================================================================

cell_type_file <- "test/coloc/cell_types.txt"
cell_gene_file <- "test/coloc/cell_gene_list.txt"
snp_position_file <- "test/coloc/reference/snp_pos.txt"

eqtl_base_dir <- "test/MR"
eqtl_locus_dir <- "test/coloc/eqtl_loci"
depression_gwas_locus_dir <- "test/coloc/depression_GWAS_loci"
harmonised_dir <- "test/coloc/harmonised"
result_dir <- "test/coloc/results"

eqtl_input_suffix <- ".filteredgene.txt"
eqtl_sample_size <- 192
window_size <- 100000
minimum_overlap_snps <- 50

# Explicit coloc priors
prior_p1 <- 1e-4
prior_p2 <- 1e-4
prior_p12 <- 1e-5

dir.create(eqtl_locus_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(harmonised_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. Read cell types, selected genes, and hg19 SNP positions
# =============================================================================

cell_types <- fread(
  cell_type_file,
  header = FALSE,
  col.names = "cell"
)

cell_types <- unique(
  cell_types[
    !is.na(cell) &
      cell != "" &
      !grepl("^#", cell)
  ]$cell
)

cell_gene_list <- fread(
  cell_gene_file,
  header = FALSE,
  col.names = c("cell", "gene")
)

cell_gene_list <- unique(
  cell_gene_list[
    !is.na(cell) &
      !is.na(gene) &
      cell != "" &
      gene != "" &
      !grepl("^#", cell)
  ]
)

snp_position <- fread(snp_position_file)

required_position_columns <- c(
  "SNP",
  "chr",
  "pos_hg19"
)

missing_columns <- setdiff(
  required_position_columns,
  names(snp_position)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns in SNP-position file: ",
    paste(missing_columns, collapse = ", ")
  )
}

snp_position <- unique(
  snp_position[, ..required_position_columns],
  by = "SNP"
)

# =============================================================================
# 3. Extract the minimum-P eQTL locus for each selected gene
# =============================================================================

extract_lead_locus <- function(eqtl, gene_name) {

  gene_dat <- eqtl[
    GeneSymbol == gene_name &
      !is.na(P) &
      !is.na(BP)
  ]

  if (nrow(gene_dat) == 0) {
    return(NULL)
  }

  minimum_p <- min(
    gene_dat$P,
    na.rm = TRUE
  )

  lead_snps <- gene_dat[
    P == minimum_p
  ]

  locus <- rbindlist(
    lapply(seq_len(nrow(lead_snps)), function(i) {

      lead_bp <- lead_snps$BP[i]

      if ("CHR" %in% names(gene_dat) &&
          !is.na(lead_snps$CHR[i])) {

        gene_dat[
          CHR == lead_snps$CHR[i] &
            abs(BP - lead_bp) <= window_size
        ]

      } else {

        gene_dat[
          abs(BP - lead_bp) <= window_size
        ]
      }
    }),
    use.names = TRUE,
    fill = TRUE
  )

  locus <- unique(
    locus,
    by = "SNP"
  )

  old_position_columns <- intersect(
    c("chr", "pos_hg19"),
    names(locus)
  )

  if (length(old_position_columns) > 0) {
    locus[, (old_position_columns) := NULL]
  }

  merge(
    locus,
    snp_position,
    by = "SNP",
    all.x = TRUE,
    sort = FALSE
  )
}

locus_summary <- list()

for (cell_id in cell_types) {

  input_file <- file.path(
    eqtl_base_dir,
    cell_id,
    paste0(cell_id, eqtl_input_suffix)
  )

  if (!file.exists(input_file)) {
    warning("Input eQTL file not found: ", input_file)
    next
  }

  eqtl <- fread(input_file)

  required_eqtl_columns <- c(
    "Gene_id",
    "SNP",
    "BP",
    "P",
    "Beta",
    "effect_allele",
    "other_allele"
  )

  missing_columns <- setdiff(
    required_eqtl_columns,
    names(eqtl)
  )

  if (length(missing_columns) > 0) {
    warning(
      "Missing columns in ",
      input_file,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
    next
  }

  eqtl[
    ,
    GeneSymbol := sub(
      "_ENSG.*",
      "",
      Gene_id
    )
  ]

  selected_genes <- unique(
    cell_gene_list$gene[
      cell_gene_list$cell == cell_id
    ]
  )

  cell_output_dir <- file.path(
    eqtl_locus_dir,
    cell_id
  )

  dir.create(
    cell_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (gene_id in selected_genes) {

    locus <- extract_lead_locus(
      eqtl,
      gene_id
    )

    if (is.null(locus) || nrow(locus) == 0) {

      locus_summary[[length(locus_summary) + 1]] <- data.table(
        Cell_type = cell_id,
        Gene = gene_id,
        Status = "NOT_FOUND",
        Lead_P = NA_real_,
        Locus_SNPs = 0L
      )

      next
    }

    keep_columns <- intersect(
      c(
        "SNP",
        "Beta",
        "P",
        "effect_allele",
        "other_allele",
        "MAF",
        "eaf",
        "chr",
        "pos_hg19"
      ),
      names(locus)
    )

    locus_output <- locus[
      ,
      ..keep_columns
    ]

    fwrite(
      locus_output,
      file.path(
        cell_output_dir,
        paste0(gene_id, "_lead100kb.txt")
      ),
      sep = "\t",
      quote = FALSE,
      na = ""
    )

    locus_summary[[length(locus_summary) + 1]] <- data.table(
      Cell_type = cell_id,
      Gene = gene_id,
      Status = "EXTRACTED",
      Lead_P = min(locus$P, na.rm = TRUE),
      Locus_SNPs = nrow(locus),
      Missing_hg19 = sum(is.na(locus$pos_hg19))
    )
  }
}

if (length(locus_summary) > 0) {
  fwrite(
    rbindlist(locus_summary, use.names = TRUE, fill = TRUE),
    file.path(result_dir, "eqtl_locus_summary.txt"),
    sep = "\t",
    quote = FALSE,
    na = ""
  )
}

# =============================================================================
# 4. Harmonize eQTL and depression GWAS loci
# =============================================================================

harmonisation_summary <- list()

for (i in seq_len(nrow(cell_gene_list))) {

  cell_id <- cell_gene_list$cell[i]
  gene_id <- cell_gene_list$gene[i]

  eqtl_file <- file.path(
    eqtl_locus_dir,
    cell_id,
    paste0(gene_id, "_lead100kb.txt")
  )

  depression_file <- file.path(
    depression_gwas_locus_dir,
    cell_id,
    paste0(gene_id, "_GWAS.txt")
  )

  if (!file.exists(eqtl_file) ||
      !file.exists(depression_file)) {

    harmonisation_summary[[length(harmonisation_summary) + 1]] <- data.table(
      Cell_type = cell_id,
      Gene = gene_id,
      Status = "INPUT_MISSING",
      Overlap_SNPs = NA_integer_,
      Harmonised_SNPs = NA_integer_
    )

    next
  }

  eqtl <- fread(eqtl_file)
  depression <- fread(depression_file)

  required_eqtl_columns <- c(
    "SNP",
    "Beta",
    "P",
    "effect_allele",
    "other_allele"
  )

  required_depression_columns <- c(
    "SNP",
    "BETA",
    "SE",
    "A1",
    "A2",
    "FRQ",
    "NCAS",
    "NCON"
  )

  if (length(setdiff(required_eqtl_columns, names(eqtl))) > 0 ||
      length(setdiff(required_depression_columns, names(depression))) > 0) {

    harmonisation_summary[[length(harmonisation_summary) + 1]] <- data.table(
      Cell_type = cell_id,
      Gene = gene_id,
      Status = "COLUMN_MISSING",
      Overlap_SNPs = NA_integer_,
      Harmonised_SNPs = NA_integer_
    )

    next
  }

  eqtl <- eqtl[
    P > 0 &
      P < 1 &
      grepl("^rs", SNP)
  ]

  if ("MAF" %in% names(eqtl)) {
    eqtl[, MAF_eqtl := MAF]
  } else if ("eaf" %in% names(eqtl)) {
    eqtl[, MAF_eqtl := pmin(eaf, 1 - eaf)]
  } else {
    harmonisation_summary[[length(harmonisation_summary) + 1]] <- data.table(
      Cell_type = cell_id,
      Gene = gene_id,
      Status = "EAF_MISSING",
      Overlap_SNPs = NA_integer_,
      Harmonised_SNPs = NA_integer_
    )
    next
  }

  eqtl[
    ,
    SE_eqtl := abs(
      Beta /
        qnorm(P / 2, lower.tail = FALSE)
    )
  ]

  eqtl <- eqtl[
    is.finite(SE_eqtl) &
      SE_eqtl > 0 &
      !is.na(MAF_eqtl) &
      MAF_eqtl > 0 &
      MAF_eqtl < 1
  ]

  eqtl <- eqtl[
    ,
    .(
      SNP,
      Beta_eqtl = Beta,
      SE_eqtl,
      MAF_eqtl,
      EA_eqtl = toupper(effect_allele),
      OA_eqtl = toupper(other_allele)
    )
  ]

  depression <- depression[
    grepl("^rs", SNP) &
      !is.na(BETA) &
      !is.na(SE) &
      !is.na(FRQ) &
      FRQ > 0 &
      FRQ < 1
  ]

  depression <- depression[
    ,
    .(
      SNP,
      Beta_depression = BETA,
      SE_depression = SE,
      FRQ_depression = FRQ,
      EA_depression = toupper(A1),
      OA_depression = toupper(A2),
      NCAS,
      NCON
    )
  ]

  dat <- merge(
    eqtl,
    depression,
    by = "SNP"
  )

  overlap_snps <- nrow(dat)

  if (overlap_snps < minimum_overlap_snps) {

    harmonisation_summary[[length(harmonisation_summary) + 1]] <- data.table(
      Cell_type = cell_id,
      Gene = gene_id,
      Status = "TOO_FEW_OVERLAP_SNPS",
      Overlap_SNPs = overlap_snps,
      Harmonised_SNPs = NA_integer_
    )

    next
  }

  dat[
    ,
    allele_status := fifelse(
      EA_eqtl == EA_depression &
        OA_eqtl == OA_depression,
      "same",
      fifelse(
        EA_eqtl == OA_depression &
          OA_eqtl == EA_depression,
        "flip",
        "remove"
      )
    )
  ]

  dat <- dat[
    allele_status != "remove"
  ]

  dat[
    allele_status == "flip",
    Beta_depression := -Beta_depression
  ]

  # Remove palindromic SNPs
  dat <- dat[
    !(
      (EA_eqtl == "A" & OA_eqtl == "T") |
        (EA_eqtl == "T" & OA_eqtl == "A") |
        (EA_eqtl == "C" & OA_eqtl == "G") |
        (EA_eqtl == "G" & OA_eqtl == "C")
    )
  ]

  dat <- unique(
    dat,
    by = "SNP"
  )

  harmonised_snps <- nrow(dat)

  if (harmonised_snps < minimum_overlap_snps) {

    harmonisation_summary[[length(harmonisation_summary) + 1]] <- data.table(
      Cell_type = cell_id,
      Gene = gene_id,
      Status = "TOO_FEW_HARMONISED_SNPS",
      Overlap_SNPs = overlap_snps,
      Harmonised_SNPs = harmonised_snps
    )

    next
  }

  cell_output_dir <- file.path(
    harmonised_dir,
    cell_id
  )

  dir.create(
    cell_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  fwrite(
    dat,
    file.path(
      cell_output_dir,
      paste0(gene_id, "_harmonised.txt")
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
  )

  harmonisation_summary[[length(harmonisation_summary) + 1]] <- data.table(
    Cell_type = cell_id,
    Gene = gene_id,
    Status = "HARMONISED",
    Overlap_SNPs = overlap_snps,
    Harmonised_SNPs = harmonised_snps
  )
}

if (length(harmonisation_summary) > 0) {
  fwrite(
    rbindlist(
      harmonisation_summary,
      use.names = TRUE,
      fill = TRUE
    ),
    file.path(
      result_dir,
      "harmonisation_summary.txt"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
  )
}

# =============================================================================
# 5. Run coloc.abf
# =============================================================================

coloc_results <- list()

for (i in seq_len(nrow(cell_gene_list))) {

  cell_id <- cell_gene_list$cell[i]
  gene_id <- cell_gene_list$gene[i]

  harmonised_file <- file.path(
    harmonised_dir,
    cell_id,
    paste0(gene_id, "_harmonised.txt")
  )

  if (!file.exists(harmonised_file)) {
    next
  }

  dat <- fread(harmonised_file)

  dat <- dat[
    !is.na(MAF_eqtl) &
      MAF_eqtl > 0 &
      MAF_eqtl < 1 &
      !is.na(FRQ_depression) &
      FRQ_depression > 0 &
      FRQ_depression < 1
  ]

  if (nrow(dat) < minimum_overlap_snps) {
    next
  }

  dat[
    ,
    MAF_depression := pmin(
      FRQ_depression,
      1 - FRQ_depression
    )
  ]

  depression_n <- round(
    mean(
      dat$NCAS + dat$NCON,
      na.rm = TRUE
    )
  )

  depression_case_fraction <- mean(
    dat$NCAS,
    na.rm = TRUE
  ) / mean(
    dat$NCAS + dat$NCON,
    na.rm = TRUE
  )

  dataset_eqtl <- list(
    beta = dat$Beta_eqtl,
    varbeta = dat$SE_eqtl^2,
    snp = dat$SNP,
    MAF = dat$MAF_eqtl,
    N = eqtl_sample_size,
    type = "quant"
  )

  dataset_depression <- list(
    beta = dat$Beta_depression,
    varbeta = dat$SE_depression^2,
    snp = dat$SNP,
    MAF = dat$MAF_depression,
    N = depression_n,
    s = depression_case_fraction,
    type = "cc"
  )

  coloc_result <- tryCatch(
    coloc.abf(
      dataset1 = dataset_eqtl,
      dataset2 = dataset_depression,
      p1 = prior_p1,
      p2 = prior_p2,
      p12 = prior_p12
    ),
    error = function(e) NULL
  )

  if (is.null(coloc_result)) {
    next
  }

  result_row <- as.data.table(
    as.list(
      coloc_result$summary
    )
  )

  result_row[
    ,
    `:=`(
      Gene = gene_id,
      Cell_type = cell_id,
      Prior_p1 = prior_p1,
      Prior_p2 = prior_p2,
      Prior_p12 = prior_p12
    )
  ]

  setcolorder(
    result_row,
    c(
      "Gene",
      "Cell_type",
      "Prior_p1",
      "Prior_p2",
      "Prior_p12",
      setdiff(
        names(result_row),
        c(
          "Gene",
          "Cell_type",
          "Prior_p1",
          "Prior_p2",
          "Prior_p12"
        )
      )
    )
  )

  coloc_results[[length(coloc_results) + 1]] <- result_row
}

if (length(coloc_results) > 0) {

  final_results <- rbindlist(
    coloc_results,
    use.names = TRUE,
    fill = TRUE
  )

  if ("PP.H4.abf" %in% names(final_results)) {
    setorder(
      final_results,
      -PP.H4.abf
    )
  }

  fwrite(
    final_results,
    file.path(
      result_dir,
      "depression_coloc_all_results.txt"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
  )

} else {

  warning("No valid colocalization results were generated.")
}
