# Illustrative input preparation: replace filenames and column names.
# These are examples, not files included with the project.
library(epicurious)

# Preferred: preprocessed M-value matrix, CpGs x samples, saved as RDS.
methylation <- readRDS("data/methylation.rds")
phenotypes <- readRDS("data/phenotypes.rds")
genotypes <- readRDS("data/genotypes.rds")
cpg_snp_map <- readRDS("data/cpg_snp_map.rds")

# If methylation is stored as a wide sample-by-CpG CSV:
# Read IDs as character so leading zeros survive. check.names = FALSE keeps IDs.
# wide <- read.csv("data/methylation.csv", check.names = FALSE,
#                  colClasses = c(sample_id = "character"), na.strings = c("", "NA"))
# methylation <- methylation_from_table(wide, sample_id = "sample_id",
#                                       cpg_cols = c("cg00000001", "cg00000002"))

# Recode categorical phenotype variables intentionally, including reference levels.
phenotypes$sex <- factor(phenotypes$sex, levels = c("female", "male"))

# Map any original study column names to the explicit standard schema:
# cpg_snp_map <- data.frame(cpg_id = original_map$CpG,
#                           snp_id = original_map$indepSNP)

# Optional input route matching the original bigsnpr-based script:
# big <- bigsnpr::snp_attach("data/genotypes.rds")
# snps <- unique(cpg_snp_map$snp_id[cpg_snp_map$cpg_id %in% rownames(methylation)])
# genotypes <- genotypes_from_bigsnp(big, snps,
#                                    sample_ids = phenotypes$sample_id)
# This adapter loads only requested SNPs into RAM; it does not run genotype QC.

fit <- ewas_mqtl(methylation, phenotypes, ~ exposure + age + sex,
                 exposure = "exposure", genotypes = genotypes,
                 cpg_snp_map = cpg_snp_map, methylation_scale = "M")
write_ewas(fit, "results/my_ewas")
