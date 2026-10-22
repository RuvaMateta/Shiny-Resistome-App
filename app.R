# app.R
# Resistome Shiny app. Run with the "Run App" button (do not use Source).
# Reads data/app_data.rds, which is made by prepare_data.R.

library(shiny)
library(dplyr)
library(DT)
library(ggplot2)

# Load the joined data once, when the app starts
d    <- readRDS("data/app_data.rds")
meta <- d$meta
amr  <- d$amr   # one row per gene hit, with isolate metadata attached

ui <- fluidPage(
  titlePanel("Resistome App"),
  sidebarLayout(
    sidebarPanel(
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
      tabsetPanel(
        tabPanel("Overview",
                 p("Filter the isolates on the left. The chart shows how many isolates carry each resistance gene, and the table lists every gene hit."),
                 textOutput("n_isolates"),
                 plotOutput("gene_plot", height = "450px"),
                 DTOutput("amr_table"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Filtered AMR rows; updates whenever an input changes
  filtered <- reactive({
    df <- amr
    if (input$species != "All") df <- df %>% filter(species == input$species)
    if (input$country != "All") df <- df %>% filter(country == input$country)
    if (input$source  != "All") df <- df %>% filter(source_category == input$source)
    df %>% filter(is.na(year_final) |
                    (year_final >= input$year[1] & year_final <= input$year[2]))
  })
  
  output$n_isolates <- renderText({
    paste(n_distinct(filtered()$sample), "isolates with resistance gene hits,",
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
      labs(x = "Gene", y = "Number of isolates",
           title = "Top 20 resistance genes") +
      theme_minimal(base_size = 14)
  })
  
  output$amr_table <- renderDT({
    df <- filtered()
    validate(need(nrow(df) > 0, "No isolates match these filters."))
    datatable(df %>% select(sample, species, country, gene, class, subclass,
                            pct_identity, pct_coverage),
              filter = "top", options = list(pageLength = 10))
  })
}

shinyApp(ui, server)