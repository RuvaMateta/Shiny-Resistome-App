# shiny-resistome-app

An interactive R Shiny app for exploring antimicrobial resistance genes, integrons, plasmids and phylogeny in 94 bacterial isolates (32 *E. coli*, 30 *Salmonella* and 32 *K. pneumoniae*) from Malawi, Zambia and South Africa.

## What the app shows

The app has four tabs.

| Tab | What it shows |
|---|---|
| **Overview** | Filter isolates by gene type (AMR, STRESS, VIRULENCE), species, country, source and year. A bar chart shows the top 20 genes, and a searchable table lists every gene hit. |
| **Resistance and integrons** | Pick a gene group (all AMR, beta-lactam, ESBL, AmpC, carbapenemase, quinolone, aminoglycoside and others). Plots show the percentage of isolates carrying the group by species and country, and the top genes in the group. A second plot shows integron frequency per species (28/32 *K. pneumoniae*, 19/30 *Salmonella*, 17/32 *E. coli*). |
| **Phylogeny** | One tree per species. Tips can be coloured by country, year or source, and isolates that carry a chosen gene or plasmid cluster are circled. |
| **Plasmid sharing** | An interactive network of isolates that share plasmids. Nodes are coloured by country, shaped by species, and lines are coloured by sharing tier. A table lists every shared pair with identity, coverage and sizes. |

Gene groups are picked by gene name. ESBL means *bla*CTX-M genes, AmpC means *bla*CMY, *bla*ACT and *bla*EC, and carbapenemase means KPC, NDM, VIM, IMP and OXA-48-like genes. Other beta-lactamases such as *bla*OXA-1 and *bla*OXA-10 are not counted as ESBLs.

The *K. pneumoniae* tree has 30 of the 32 isolates (ERR13744119 and ERR13744124 are not in it). Each tree also has a tip called `Reference`, the reference genome used to build the tree. It has no metadata, so it is drawn in its own colour.

## Packages needed

R 4.x and these packages:

```r
install.packages(c("shiny", "dplyr", "readr", "ggplot2", "DT", "plotly",
                   "ape", "visNetwork", "BiocManager"))
BiocManager::install(c("ggtree", "ggtreeExtra"))
```

## How to run it

1. Download the repository (green **Code** button, then **Download ZIP**) and unzip it.
2. Open `Shiny-Resistome-App.Rproj` in RStudio, or set the working directory to the unzipped folder.
3. Run:

```r
shiny::runApp()
```

The app reads `data/app_data.rds`. If that file is missing, run `prepare_data.R` once to create it from the result tables in `data/`.

## Files

```
app.R                  the Shiny app (all four tabs)
prepare_data.R         reads the result tables, checks sample names against
                       the metadata and tree tips, joins them and saves
                       data/app_data.rds
data/
  metadata_final.csv   isolate metadata (species, country, source, year)
  amr_all.tsv          AMRFinderPlus results for all isolates
  integron_contigs_all.tsv, integron_elements_all.tsv
                       integron results
  plasmid_cluster_amr_integron.tsv
                       plasmid clusters with their AMR genes and integrons
  cross_species_tiers.tsv
                       plasmids shared between isolates, with sharing tier
  E_coli.treefile, Salmonella.treefile, K_pneumoniae_no119.treefile
                       phylogenetic trees (Newick)
  app_data.rds         joined tables the app loads at start-up
```

Raw reads and assemblies are not included.
