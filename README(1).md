# Shiny-Resistome-App

An interactive R Shiny app for exploring antimicrobial resistance genes, integrons, plasmids and phylogeny in 94 bacterial isolates (32 *E. coli*, 30 *Salmonella* and 32 *K. pneumoniae*) from Malawi, Zambia and South Africa.

## What the app shows

The app has four tabs.

| Tab | What it shows |
|---|---|
| **Overview** | Filter isolates by gene type (AMR, STRESS, VIRULENCE), species, country, source and year. A bar chart shows the top 20 genes, and a searchable table lists every gene hit. |
| **Resistance and integrons** | Pick a gene group (all AMR, beta-lactam, ESBL, AmpC, carbapenemase, quinolone, aminoglycoside and others). Plots show the percentage of isolates carrying the group by species and country, and the top genes in the group. A second plot and table show, for each species, how many isolates carry any integron-related element, a complete integron, a CALIN or an In0 element. An isolate can carry more than one category, so the categories overlap and should not be added together. |
| **Phylogeny** | One tree per species. Tips can be coloured by country, year or source, and isolates that carry a chosen gene or plasmid cluster are circled. |
| **Plasmid sharing** | An interactive network of isolates that share plasmids. Nodes are coloured by country, shaped by species, and lines are coloured by sharing tier. A table lists every shared pair with identity, coverage and sizes. |

Gene groups are picked by gene name. ESBL means *bla*CTX-M genes, AmpC means *bla*CMY, *bla*ACT and *bla*EC, and carbapenemase means KPC, NDM, VIM, IMP and OXA-48-like genes. Other beta-lactamases such as *bla*OXA-1 and *bla*OXA-10 are not counted as ESBLs.

## Integron categories

Counts are unique isolates in the reviewed table, not independently validated biological totals. The categories overlap within an isolate, so the counts must not be added to calculate prevalence.

| Species | Any detected category | Complete integron | CALIN | In0 |
|---|---|---|---|---|
| *E. coli* | 17/32 (53.1%) | 16/32 (50.0%) | 8/32 (25.0%) | 0/32 (0%) |
| *K. pneumoniae* | 28/32 (87.5%) | 19/32 (59.4%) | 19/32 (59.4%) | 8/32 (25.0%) |
| *Salmonella* | 19/30 (63.3%) | 15/30 (50.0%) | 4/30 (13.3%) | 0/30 (0%) |

- **Complete integron:** an integrase and an attC site.
- **CALIN:** attC sites without a detected integron integrase.
- **In0:** an integrase without detected attC sites.

Isolates with no integron rows are counted as having no detections only because their IntegronFinder runs were confirmed as completed on Galaxy.

## Phylogeny notes

The *K. pneumoniae* tree has 30 of the 32 isolates. ERR13744119 was left out (the tree file is named `K_pneumoniae_no119`) and ERR13744124 was used as the reference genome. Each tree also has a tip called `Reference`, the reference genome used to build the tree. It has no metadata, so it is drawn in its own colour.

## Packages needed

R 4.x and these packages:

```r
install.packages(c("shiny", "dplyr", "readr", "tidyr", "ggplot2", "DT", "plotly",
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
                       the metadata and tree tips (it stops on a duplicated
                       or unmatched sample name), counts integron elements by
                       category, joins the tables and saves data/app_data.rds
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

The `data/` folder also holds other result tables (`contig_report_all.tsv`, `mobtyper_all.tsv`, `mge_report_all.tsv`, `cross_species_pairs.tsv`, `plasmid_pairs.tsv`, `shared_region_genes.tsv`) and the IQ-TREE reports (`.iqtree`). The app does not use them.

Raw reads and assemblies are not included.

## Live app

https://ruvarashe.shinyapps.io/shiny-resistome-app/
