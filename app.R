# app.R
# Resistome Shiny app. Run with the "Run App" button (do not use Source).
# Reads data/app_data.rds, which is made by prepare_data.R.

library(shiny)
library(dplyr)
library(DT)
library(ggplot2)
library(ape)
library(ggtree)
library(visNetwork)

# ---- Load the joined data once, when the app starts ----
d    <- readRDS("data/app_data.rds")
meta <- d$meta   # one row per isolate, with n_integrons
amr  <- d$amr    # one row per gene hit, with isolate metadata attached

# ---- Gene groups for the Resistance and integrons tab ----
# Each group is a rule that picks rows of the AMRFinderPlus table.
# AMRFinderPlus has no ESBL label, so ESBL, AmpC and carbapenemase groups
# are picked by gene name. ESBL here means blaCTX-M genes. Other beta-lactamases
# (blaOXA-1, blaOXA-10, blaSHV variants) are not counted as ESBL.
is_amr <- function(x) x$type %in% "AMR"
gene_groups <- list(
  "All antibiotic resistance genes"          = function(x) is_amr(x),
  "Beta-lactam (all)"                        = function(x) is_amr(x) & grepl("BETA-LACTAM", x$class),
  "ESBL (blaCTX-M genes)"                    = function(x) is_amr(x) & grepl("^blaCTX-M", x$gene),
  "AmpC (blaCMY, blaACT, blaEC)"             = function(x) is_amr(x) & grepl("^bla(CMY|ACT|EC)", x$gene),
  "Carbapenemase (KPC, NDM, VIM, IMP, OXA-48-like)" = function(x) is_amr(x) & grepl("^bla(KPC|NDM|VIM|IMP|OXA-(48|162|181|204|232|244))", x$gene),
  "Quinolone"                                = function(x) is_amr(x) & grepl("QUINOLONE", x$class),
  "Aminoglycoside"                           = function(x) is_amr(x) & grepl("AMINOGLYCOSIDE", x$class),
  "Tetracycline"                             = function(x) is_amr(x) & grepl("TETRACYCLINE", x$class),
  "Sulfonamide"                              = function(x) is_amr(x) & grepl("SULFONAMIDE", x$class),
  "Trimethoprim"                             = function(x) is_amr(x) & grepl("TRIMETHOPRIM", x$class),
  "Phenicol"                                 = function(x) is_amr(x) & grepl("PHENICOL", x$class),
  "Macrolide"                                = function(x) is_amr(x) & grepl("MACROLIDE", x$class),
  "Fosfomycin"                               = function(x) is_amr(x) & grepl("FOSFOMYCIN", x$class),
  "Metal and biocide resistance (STRESS)"    = function(x) x$type %in% "STRESS",
  "Virulence genes"                          = function(x) x$type %in% "VIRULENCE"
)

# ---- Integron frequency per species (from the metadata table) ----
integron_df <- meta %>%
  group_by(species) %>%
  summarise(isolates = n(),
            with_integron = sum(n_integrons > 0),
            .groups = "drop") %>%
  mutate(pct = 100 * with_integron / isolates,
         label = paste0(with_integron, "/", isolates))

# ---- Phylogenetic trees (one per species) ----
trees <- list(
  E_coli       = read.tree("data/E_coli.treefile"),
  Salmonella   = read.tree("data/Salmonella.treefile"),
  K_pneumoniae = read.tree("data/K_pneumoniae_no119.treefile")
)
pcl <- d$pcl   # plasmid clusters per isolate

# ---- Plasmid sharing (cross_species_tiers.tsv) ----
tiers <- d$tiers
tier_col <- grep("^tier", names(tiers), value = TRUE)[1]
tiers$tier_use <- if (is.na(tier_col)) "all" else as.character(tiers[[tier_col]])
tier_levels <- sort(unique(tiers$tier_use))
tier_cols <- setNames(hcl.colors(max(length(tier_levels), 3), "Set 2")[seq_along(tier_levels)],
                      tier_levels)
countries <- sort(unique(meta$country))
country_cols <- setNames(hcl.colors(max(length(countries), 3), "Dark 3")[seq_along(countries)],
                         countries)
species_shapes <- c(E_coli = "dot", Salmonella = "square", K_pneumoniae = "triangle")

# ---- User interface ----
ui <- fluidPage(
  titlePanel("Resistome App"),
  p(paste0("Antimicrobial resistance genes, integrons, plasmids and phylogeny of ",
           nrow(meta), " bacterial isolates (",
           paste(gsub("_", " ", sort(unique(meta$species))), collapse = ", "),
           ") from ", paste(sort(unique(meta$country)), collapse = ", "), ".")),
  tabsetPanel(

    # Tab 1: overview with filters
    tabPanel("Overview",
      sidebarLayout(
        sidebarPanel(
          selectInput("gene_type", "Gene type",
                      choices = c("All", "AMR", "STRESS", "VIRULENCE"),
                      selected = "AMR"),
          selectInput("species", "Species",
                      choices = c("All", sort(unique(meta$species)))),
          selectInput("country", "Country",
                      choices = c("All", sort(unique(meta$country)))),
          selectInput("source", "Source",
                      choices = c("All", sort(unique(meta$source_category)))),
          sliderInput("year", "Year",
                      min = min(meta$year_final, na.rm = TRUE),
                      max = max(meta$year_final, na.rm = TRUE),
                      value = range(meta$year_final, na.rm = TRUE),
                      step = 1, sep = ""),
          actionButton("reset", "Reset filters")
        ),
        mainPanel(
          h3("Genes found in the isolates"),
          p("Filter the isolates on the left. The chart shows how many isolates carry each gene, and the table lists every gene hit. AMR = antibiotic resistance, STRESS = metal and biocide resistance."),
          textOutput("n_isolates"),
          plotOutput("gene_plot", height = "450px"),
          DTOutput("amr_table")
        )
      )
    ),

    # Tab 2: resistance groups and integrons
    tabPanel("Resistance and integrons",
      sidebarLayout(
        sidebarPanel(
          selectInput("gene_group", "Gene group",
                      choices = names(gene_groups)),
          p("ESBL = blaCTX-M genes. AmpC = blaCMY, blaACT and blaEC. Carbapenemase = KPC, NDM, VIM, IMP and OXA-48-like genes. Groups are picked by gene name.")
        ),
        mainPanel(
          h3("Resistance groups and integrons"),
          p("Choose a gene group on the left to see which species and countries carry it, then see how common integrons are in each species."),
          h4("Isolates carrying a gene from the group, by species and country"),
          p("Each bar is the percentage of isolates in that species and country that carry at least one gene in the group. Labels show carriers / isolates."),
          plotOutput("group_plot", height = "400px"),
          h4("Genes in the group, by species"),
          plotOutput("gene_species_plot", height = "450px"),
          h4("Integron frequency per species"),
          p("Number of isolates with at least one integron, out of all isolates of that species."),
          plotOutput("integron_plot", height = "300px"),
          DTOutput("integron_table")
        )
      )
    ),

    # Tab 3: phylogeny
    tabPanel("Phylogeny",
      sidebarLayout(
        sidebarPanel(
          selectInput("tree_species", "Species",
                      choices = names(trees)),
          selectInput("tree_color", "Colour tips by",
                      choices = c("Country", "Year", "Source")),
          radioButtons("hl_type", "Circle isolates that carry a",
                       choices = c("None", "Gene", "Plasmid cluster"),
                       selected = "None"),
          uiOutput("hl_choice"),
          p("The tip called Reference is the reference genome used to build the tree. It has no metadata, so it is drawn in its own colour.")
        ),
        mainPanel(
          h3("Phylogeny"),
          p("One tree per species. Tips are coloured by country, year or source, and a black circle marks isolates that carry the chosen gene or plasmid cluster."),
          textOutput("tree_note"),
          textOutput("hl_note"),
          plotOutput("tree_plot", height = "750px")
        )
      )
    ),

    # Tab 4: plasmid sharing
    tabPanel("Plasmid sharing",
      sidebarLayout(
        sidebarPanel(
          checkboxGroupInput("tier_pick", "Sharing tier",
                             choices = tier_levels, selected = tier_levels),
          p("Each node is an isolate and each line is a plasmid shared by two isolates. Node colour is the country; node shape is the species (circle = E. coli, square = Salmonella, triangle = K. pneumoniae). Line colour is the sharing tier. Hover over a node or line for details, and drag nodes to rearrange them."),
          uiOutput("share_legend")
        ),
        mainPanel(
          h3("Plasmid sharing between isolates"),
          p("Isolates that share a plasmid are joined by a line. Untick a tier on the left to hide it."),
          visNetworkOutput("plasmid_net", height = "650px"),
          h4("Shared plasmid pairs"),
          DTOutput("tiers_table")
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # Reset button on the Overview tab
  observeEvent(input$reset, {
    updateSelectInput(session, "gene_type", selected = "AMR")
    updateSelectInput(session, "species", selected = "All")
    updateSelectInput(session, "country", selected = "All")
    updateSelectInput(session, "source", selected = "All")
    updateSliderInput(session, "year", value = range(meta$year_final, na.rm = TRUE))
  })

  # ----- Overview tab -----
  filtered <- reactive({
    df <- amr
    if (input$gene_type != "All") df <- df %>% filter(type == input$gene_type)
    if (input$species != "All") df <- df %>% filter(species == input$species)
    if (input$country != "All") df <- df %>% filter(country == input$country)
    if (input$source  != "All") df <- df %>% filter(source_category == input$source)
    df %>% filter(is.na(year_final) |
                  (year_final >= input$year[1] & year_final <= input$year[2]))
  })

  output$n_isolates <- renderText({
    df <- filtered()
    if (nrow(df) == 0) {
      return("No isolates match these filters. Try choosing All or widening the year range.")
    }
    paste(n_distinct(df$sample), "isolates with gene hits,",
          nrow(df), "gene hits")
  })

  output$gene_plot <- renderPlot({
    df <- filtered()
    validate(need(nrow(df) > 0, "No isolates match these filters."))
    top <- df %>%
      group_by(gene) %>%
      summarise(isolates = n_distinct(sample), .groups = "drop") %>%
      slice_max(isolates, n = 20, with_ties = FALSE)
    ggplot(top, aes(x = reorder(gene, isolates), y = isolates)) +
      geom_col(fill = "steelblue") +
      coord_flip() +
      labs(x = "Gene", y = "Number of isolates", title = "Top 20 genes") +
      theme_minimal(base_size = 14)
  })

  output$amr_table <- renderDT({
    df <- filtered()
    validate(need(nrow(df) > 0, "No isolates match these filters."))
    datatable(df %>% select(sample, species, country, gene, type, class,
                            subclass, pct_identity, pct_coverage),
              filter = "top", options = list(pageLength = 10))
  })

  # ----- Resistance and integrons tab -----
  group_hits <- reactive({
    rule <- gene_groups[[input$gene_group]]
    amr[rule(amr), ]
  })

  output$group_plot <- renderPlot({
    hits <- group_hits()
    carriers <- hits %>%
      distinct(species, country, sample) %>%
      count(species, country, name = "carriers")
    totals <- meta %>% count(species, country, name = "total")
    df <- totals %>%
      left_join(carriers, by = c("species", "country")) %>%
      mutate(carriers = coalesce(carriers, 0L),
             pct = 100 * carriers / total,
             label = paste0(carriers, "/", total))
    ggplot(df, aes(x = country, y = pct, fill = country)) +
      geom_col(show.legend = FALSE) +
      geom_text(aes(label = label), vjust = -0.4, size = 4) +
      facet_wrap(~ species) +
      scale_y_continuous(limits = c(0, 110)) +
      labs(x = NULL, y = "% of isolates", title = input$gene_group) +
      theme_minimal(base_size = 14)
  })

  output$gene_species_plot <- renderPlot({
    hits <- group_hits()
    validate(need(nrow(hits) > 0, "No genes from this group were found."))
    g <- hits %>%
      distinct(sample, species, gene) %>%
      count(gene, species, name = "isolates")
    top <- g %>%
      group_by(gene) %>%
      summarise(total = sum(isolates), .groups = "drop") %>%
      slice_max(total, n = 15, with_ties = FALSE)
    g <- g %>% filter(gene %in% top$gene)
    ggplot(g, aes(x = reorder(gene, isolates, FUN = sum), y = isolates, fill = species)) +
      geom_col() +
      coord_flip() +
      labs(x = "Gene", y = "Number of isolates", fill = "Species",
           title = "Top 15 genes in the group") +
      theme_minimal(base_size = 14)
  })

  output$integron_plot <- renderPlot({
    ggplot(integron_df, aes(x = species, y = pct, fill = species)) +
      geom_col(show.legend = FALSE) +
      geom_text(aes(label = label), vjust = -0.4, size = 5) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = NULL, y = "% of isolates with an integron") +
      theme_minimal(base_size = 14)
  })

  output$integron_table <- renderDT({
    datatable(integron_df %>%
                transmute(Species = species,
                          Isolates = isolates,
                          `With integron` = with_integron,
                          `Percent` = round(pct, 1)),
              rownames = FALSE, options = list(dom = "t"))
  })

  # ----- Phylogeny tab -----
  tree_data <- reactive({
    req(input$tree_species)
    tr  <- trees[[input$tree_species]]
    col <- c("Country" = "country", "Year" = "year_final",
             "Source" = "source_category")[[input$tree_color]]
    tipdf <- data.frame(label = tr$tip.label, stringsAsFactors = FALSE) %>%
      left_join(meta %>% select(sample, country, year_final, source_category),
                by = c("label" = "sample"))
    tipdf$colvar <- as.character(tipdf[[col]])
    tipdf$colvar[is.na(tipdf$colvar)] <- "No metadata (Reference)"
    list(tree = tr, tips = tipdf)
  })

  # Gene or plasmid cluster choices depend on the species
  output$hl_choice <- renderUI({
    sp <- input$tree_species
    if (input$hl_type == "Gene") {
      genes <- sort(unique(amr$gene[amr$species == sp & amr$type %in% "AMR"]))
      selectInput("hl_gene", "Gene", choices = genes)
    } else if (input$hl_type == "Plasmid cluster") {
      cl <- sort(unique(pcl$cluster[pcl$sample %in% meta$sample[meta$species == sp]]))
      selectInput("hl_cluster", "Plasmid cluster", choices = cl)
    }
  })

  # Isolates to circle on the tree
  carriers <- reactive({
    sp <- input$tree_species
    if (input$hl_type == "Gene" && !is.null(input$hl_gene)) {
      unique(amr$sample[amr$species == sp & amr$gene == input$hl_gene])
    } else if (input$hl_type == "Plasmid cluster" && !is.null(input$hl_cluster)) {
      unique(pcl$sample[pcl$cluster == input$hl_cluster])
    } else {
      character(0)
    }
  })

  output$tree_note <- renderText({
    tips <- trees[[input$tree_species]]$tip.label
    all_iso <- meta$sample[meta$species == input$tree_species]
    missing <- setdiff(all_iso, tips)
    paste0(length(intersect(tips, all_iso)), " of ", length(all_iso),
           " isolates are in this tree.",
           if (length(missing) > 0)
             paste0(" Not in the tree: ", paste(missing, collapse = ", "), "."))
  })

  output$hl_note <- renderText({
    if (input$hl_type == "None") return("")
    tips <- trees[[input$tree_species]]$tip.label
    paste0(length(intersect(carriers(), tips)),
           " isolates in this tree carry the selected ",
           tolower(input$hl_type), " (circled).")
  })

  output$tree_plot <- renderPlot({
    td <- tree_data()
    hl <- carriers()
    p <- ggtree(td$tree) %<+% td$tips +
      geom_tippoint(aes(colour = colvar), size = 3) +
      geom_tiplab(size = 3.2) +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.4))) +
      labs(colour = input$tree_color) +
      theme(legend.position = "right")
    if (length(hl) > 0) {
      p <- p + geom_point(
        data = function(x) x[x$isTip & x$label %in% hl, ],
        aes(x = x, y = y),
        shape = 21, size = 7, fill = NA, colour = "black", stroke = 1.3,
        inherit.aes = FALSE)
    }
    p
  })

  # ----- Plasmid sharing tab -----
  tier_filtered <- reactive({
    tiers %>% filter(tier_use %in% input$tier_pick)
  })

  output$plasmid_net <- renderVisNetwork({
    df <- tier_filtered()
    validate(need(nrow(df) > 0, "No plasmid pairs match these filters."))

    ids <- unique(c(df$qsample, df$ssample))
    nodes <- data.frame(id = ids, stringsAsFactors = FALSE) %>%
      left_join(meta %>% select(sample, species, country), by = c("id" = "sample"))
    nodes$label <- sub("^(E_coli|Salmonella|K_pneumoniae)_", "", nodes$id)
    nodes$color <- unname(country_cols[nodes$country])
    nodes$shape <- unname(species_shapes[nodes$species])
    nodes$title <- paste0("<b>", nodes$id, "</b><br>", nodes$species,
                          "<br>", nodes$country)

    edges <- data.frame(from = df$qsample, to = df$ssample,
                        stringsAsFactors = FALSE)
    edges$color <- unname(tier_cols[df$tier_use])
    edges$width <- 3
    edges$title <- paste0("Tier: ", df$tier_use,
                          "<br>", df$qcl, " and ", df$scl,
                          "<br>Identity: ", df$ident,
                          "<br>Coverage of smaller plasmid: ", df$cov_smaller,
                          "<br>Sizes: ", df$qsize, " and ", df$ssize)

    nodes$size <- 25
    nodes$font.size <- 16

    visNetwork(nodes, edges) %>%
      visOptions(highlightNearest = TRUE) %>%
      visEdges(smooth = FALSE)
  })

  # Legend drawn in the sidebar so it never overlaps the network
  output$share_legend <- renderUI({
    dot <- function(col, sym) {
      sprintf('<span style="color:%s;font-size:20px;">%s</span>', col, sym)
    }
    line <- function(col) {
      sprintf('<span style="display:inline-block;width:30px;height:5px;background:%s;vertical-align:middle;"></span>', col)
    }
    country_html <- paste0(dot(unname(country_cols), "&#9679;"), " ",
                           names(country_cols), collapse = "<br>")
    species_html <- paste0(
      dot("#555555", "&#9679;"), " E. coli<br>",
      dot("#555555", "&#9632;"), " Salmonella<br>",
      dot("#555555", "&#9650;"), " K. pneumoniae")
    tier_html <- paste0(line(unname(tier_cols)), " ", names(tier_cols),
                        collapse = "<br>")
    HTML(paste0(
      "<hr><b>Country (node colour)</b><br>", country_html,
      "<br><br><b>Species (node shape)</b><br>", species_html,
      "<br><br><b>Sharing tier (line colour)</b><br>", tier_html))
  })

  output$tiers_table <- renderDT({
    df <- tier_filtered()
    validate(need(nrow(df) > 0, "No plasmid pairs match these filters."))
    datatable(df %>%
                select(any_of(c("qsample", "qspecies", "q_country", "qcl", "qsize",
                                "ssample", "sspecies", "s_country", "scl", "ssize",
                                "ident", "cov_smaller", "tier_use"))) %>%
                rename(tier = tier_use),
              filter = "top", rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE))
  })
}

shinyApp(ui, server)
