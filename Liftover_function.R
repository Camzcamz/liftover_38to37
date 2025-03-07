library(plotrix)
library(glue)
library(data.table)
library(dplyr)
library(rtracklayer)

setwd("")

# Function to perform liftover
perform_liftover <- function(input_file, output_file, chain_file,
                             chr_col, pos_col, rsid_col, effect_allele_col,
                             other_allele_col, af_col = NULL, beta_col, se_col, p_col, n_col) {
  
  # Read input summary statistics
  sumstats <- fread(input_file)
  
  # Ensure columns exist
  required_cols <- c(chr_col, pos_col, rsid_col, effect_allele_col, other_allele_col, beta_col, se_col, p_col, n_col)
  if (!all(required_cols %in% colnames(sumstats))) {
    stop("One or more required columns are missing in the input file.")
  }
  
  # Filter SNPs
  sumstats <- sumstats %>%
    filter(!!sym(effect_allele_col) != "D", !!sym(other_allele_col) != "I",
           startsWith(!!sym(rsid_col), "rs"),
           !is.na(!!sym(p_col)) & !!sym(p_col) != "",
           !!sym(chr_col) != "X")
  
  # Ensure chromosome column is formatted correctly
  sumstats[[chr_col]] <- as.character(sumstats[[chr_col]])
  sumstats[[chr_col]] <- gsub("^chr", "", sumstats[[chr_col]])  # Remove "chr" prefix if present
  sumstats[[chr_col]] <- as.numeric(sumstats[[chr_col]])        # Convert to numeric
  
  # Liftover preparation
  chainObject <- import.chain(chain_file)
  sumstats$old_BP <- sumstats[[pos_col]]
  sumstats$startField <- sumstats[[pos_col]]
  sumstats$endField <- sumstats[[pos_col]]
  sumstats$newChromosomeColumn <- paste0("chr", sumstats[[chr_col]])
  
  # Create genomic ranges
  grGWAS_SNPs <- makeGRangesFromDataFrame(
    sumstats,
    seqnames.field = "newChromosomeColumn",
    start.field = "startField",
    end.field = "endField", 
    keep.extra.columns = TRUE
  )
  
  # Perform liftover
  lifted_results <- liftOver(grGWAS_SNPs, chainObject)
  results <- unlist(lifted_results)  # Unlist before converting to a data frame
  results <- as.data.frame(results)
  
  # Ensure POS column exists
  if (!"start" %in% colnames(results)) {
    stop("Liftover failed: No positions were mapped.")
  }
  
  results$POS <- results$start 
  
  # Select relevant columns
  selected_columns <- c(
    "CHR" = chr_col, "POS", "RSID" = rsid_col,
    "EFFECT_ALLELE" = effect_allele_col, "OTHER_ALLELE" = other_allele_col,
    "BETA" = beta_col, "SE" = se_col, "P" = p_col, "N" = n_col
  )
  
  if (!is.null(af_col) && af_col %in% colnames(sumstats)) {
    selected_columns["AF_1000G"] <- af_col
  }
  
  # Ensure columns exist before selection
  selected_cols <- unname(selected_columns[selected_columns %in% colnames(results)])
  results <- dplyr::select(results, all_of(selected_cols))
  
  # Save output
  fwrite(results, output_file, sep = "\t")
  
  # Debugging output
  cat("Before liftover:", nrow(sumstats), " SNPs\n")
  cat("After liftover:", nrow(results), " SNPs\n")
}

# Apply function to input summary statistics
perform_liftover(
  input_file = "GWAS.txt",
  output_file = "lifted_GWAS.txt",
  chain_file = "hg38ToHg19.over.chain",
  chr_col = "CHR",
  pos_col = "POS",
  rsid_col = "RSID",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  af_col = "AF_1000G",
  beta_col = "BETA",
  se_col = "SE",
  p_col = "P",
  n_col = "N"
)
