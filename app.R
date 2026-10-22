# app.R
# Resistome Shiny app. Run with the "Run App" button (do not use Source).
# Reads data/app_data.rds, which is made by prepare_data.R.

library(shiny)
library(dplyr)
library(DT)
library(ggplot2)

# ---- Load the joined data once, when the app starts ----
d    <- readRDS("data/app_data.rds")
meta <- d$meta   # one row per isolate, with n_integrons
amr  <- d$amr    # one row per gene hit, with isolate metadata attached

# ---- Gene groups for the Resistance and integrons tab ----
# Each group is a rule that picks rows of the AMRFinderPlus table.
# AMRFinderPlus has no ESBL label, so the cephalosporin and carbapenem
# groups are the closest candidates; they are not confirmed ESBLs or
# carbapenemases.
is_amr <- function(x) x$type %in% "AMR"
gene_groups <- list(
  "All antibiotic resistance genes"          = function(x) is_amr(x),
  "Beta-lactam (all)"                        = function(x) is_amr(x) & grepl("BETA-LACTAM", x$class),
  "Cephalosporin (ESBL / AmpC candidates)"   = function(x) is_amr(x) & grepl("CEPHALOSPORIN", x$subclass),
  "Carbapenem (carbapenemase candidates)"    = function(x) is_amr(x) & grepl("CARBAPENEM", x$subclass),
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

# ---- User interface ----
ui <- fluidPage(
  titlePanel("Resistome App"),
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
                      step = 1, sep = "")
        ),
        mainPanel(
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
          p("Cephalosporin and carbapenem groups are candidates picked from AMRFinderPlus drug classes. They are not confirmed ESBLs or carbapenemases.")
        ),
        mainPanel(
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
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

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
    paste(n_distinct(filtered()$sample), "isolates with gene hits,",
          nrow(filtered()), "gene hits")
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
}

shinyApp(ui, server)
