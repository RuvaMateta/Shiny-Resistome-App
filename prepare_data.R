# prepare_data.R
# Reads the final result tables from data/, checks sample names against the
# metadata and tree tip labels, joins the tables, and saves data/app_data.rds
# for the Shiny app to load once.

library(readr)
library(dplyr)
library(ape)

# ---- 1. Read the tables ----
meta  <- read_csv("data/metadata_final.csv", show_col_types = FALSE)
amr   <- read_tsv("data/amr_all.tsv", show_col_types = FALSE)
icont <- read_tsv("data/integron_contigs_all.tsv", show_col_types = FALSE)
iele  <- read_tsv("data/integron_elements_all.tsv", show_col_types = FALSE)
pcl   <- read_tsv("data/plasmid_cluster_amr_integron.tsv", show_col_types = FALSE)
tiers <- read_tsv("data/cross_species_tiers.tsv", show_col_types = FALSE)

# ---- 2. Check sample names against metadata (character(0) = all match) ----
print(setdiff(amr$sample, meta$sample))
print(setdiff(icont$sample, meta$sample))
print(setdiff(iele$sample, meta$sample))
print(setdiff(pcl$sample, meta$sample))
print(setdiff(c(tiers$qsample, tiers$ssample), meta$sample))

# ---- 3. Check tree tip labels ----
trees <- list(
  E_coli       = read.tree("data/E_coli.treefile"),
  Salmonella   = read.tree("data/Salmonella.treefile"),
  K_pneumoniae = read.tree("data/K_pneumoniae_no119.treefile")
)
for (n in names(trees)) {
  tips <- trees[[n]]$tip.label
  cat("\n==", n, "==", length(tips), "tips\n")
  cat("Tips not in meta:", tips[!tips %in% meta$sample], "\n")
  cat("Isolates missing from tree:",
      setdiff(meta$sample[meta$species == n], tips), "\n")
}

# ---- 4. Join the tables ----
amr_meta <- amr %>%
  left_join(meta, by = "sample")

integron_by_sample <- icont %>%
  group_by(sample) %>%
  summarise(n_integrons = n(), .groups = "drop")

meta_full <- meta %>%
  left_join(integron_by_sample, by = "sample") %>%
  mutate(n_integrons = ifelse(is.na(n_integrons), 0L, n_integrons))

# ---- 5. Save for the app ----
saveRDS(list(meta = meta_full, amr = amr_meta, icont = icont, iele = iele,
             pcl = pcl, tiers = tiers),
        "data/app_data.rds")

# ---- 6. Quick check of the saved file ----
d <- readRDS("data/app_data.rds")
print(table(d$meta$species, d$meta$n_integrons > 0))
print(sum(is.na(d$amr$species)))
