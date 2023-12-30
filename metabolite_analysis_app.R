library(shiny)
library(readr)
library(DescTools)
library(DHARMa)
library(ggplot2)
library(dplyr)
library(multcomp)

# Function to perform ANOVA analysis
perform_anova <- function(data) {
  features <- data[3:ncol(data)]
  groups <- data[2, ]

  anova_results <- lapply(features, function(feature) {
    model <- aov(get(feature) ~ factor(groups))
    result <- summary(model)
    return(result)
  })

  return(anova_results)
}

# Function to perform post-hoc analysis
perform_posthoc <- function(data, selected_feature, threshold, method) {
  groups <- data[2, ]
  feature_data <- data[selected_feature, ]

  model <- aov(feature_data ~ factor(groups))
  posthoc_results <- glht(model, linfct = mcp(groups = method))
  confint_results <- confint(posthoc_results, level = 1 - threshold)

  return(confint_results)
}

# Function to display boxplot with compact letter display
display_boxplot <- function(data, selected_feature, posthoc_results) {
  boxplot_data <- data.frame(
    Groups = as.factor(data[2, ]),
    Value = as.numeric(data[selected_feature, ])
  )

  ggplot(boxplot_data, aes(x = Groups, y = Value, fill = Groups)) +
    geom_boxplot() +
    geom_text(data = as.data.frame(posthoc_results), aes(x = Groups1, y = Inf, label = Group1), vjust = -0.5) +
    labs(title = paste("Boxplot for", selected_feature),
         x = "Groups",
         y = selected_feature) +
    theme_minimal()
}

# Function to display individual metabolite boxplot
display_metabolite_boxplot <- function(data, selected_feature) {
  boxplot_data <- data.frame(
    Groups = as.factor(data[2, ]),
    Value = as.numeric(data[selected_feature, ])
  )

  ggplot(boxplot_data, aes(x = Groups, y = Value, fill = Groups)) +
    geom_boxplot() +
    labs(title = paste("Boxplot for", selected_feature),
         x = "Groups",
         y = selected_feature) +
    theme_minimal()
}

# Define UI
ui <- fluidPage(
  titlePanel("Metabolite Analysis App"),

  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = c(".csv")),
      radioButtons("analysis_tab", "Select Analysis", choices = c("File Information", "ANOVA Results", "Post-hoc Analysis", "Boxplot with Letters", "Individual Metabolite Boxplot"), selected = "File Information"),
      conditionalPanel(
        condition = "input.analysis_tab == 'Post-hoc Analysis'",
        selectInput("selected_feature_posthoc", "Select Feature for Post-hoc Analysis", choices = NULL),
        sliderInput("threshold", "Significance Threshold", min = 0.01, max = 0.5, step = 0.01, value = 0.05),
        selectInput("method", "Select Post-hoc Method", choices = c("bonferroni", "sidak", "holm", "fdr"))
      ),
      conditionalPanel(
        condition = "input.analysis_tab %in% c('Boxplot with Letters', 'Individual Metabolite Boxplot')",
        selectInput("selected_feature", "Select Feature", choices = NULL)
      )
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("File Information", tableOutput("file_info")),
        tabPanel("ANOVA Results", verbatimTextOutput("anova_results")),
        tabPanel("Post-hoc Analysis", tableOutput("posthoc_results")),
        tabPanel("Boxplot with Letters", plotOutput("boxplot_letters")),
        tabPanel("Individual Metabolite Boxplot", plotOutput("individual_boxplot"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  data <- reactive({
    req(input$file)
    read.csv(input$file$datapath, header = TRUE)
  })

  # Display file information
  output$file_info <- renderTable({
    data_info <- data()
    data.frame(
      "Number of samples" = ncol(data_info) - 2,
      "Number of groups" = length(table(data_info[2, ])),
      "Groups and sample counts" = table(data_info[2, ])
    )
  })

  # Perform ANOVA analysis
  output$anova_results <- renderPrint({
    anova_results <- perform_anova(data())
    lapply(names(anova_results), function(feature) {
      cat("Feature:", feature, "\n")
      print(anova_results[[feature]])
    })
  })

  # Update choices for selected feature inputs
  observe({
    updateSelectInput(session, "selected_feature_posthoc", choices = names(data())[3:ncol(data())])
    updateSelectInput(session, "selected_feature", choices = names(data())[3:ncol(data())])
  })

  # Perform post-hoc analysis
  output$posthoc_results <- renderTable({
    selected_feature_posthoc <- input$selected_feature_posthoc
    threshold <- input$threshold
    method <- input$method

    posthoc_results <- perform_posthoc(data(), selected_feature_posthoc, threshold, method)
    posthoc_results <- cbind.data.frame(posthoc_results$confint)
    posthoc_results
  })

  # Display boxplot with compact letter display
  output$boxplot_letters <- renderPlot({
    selected_feature_boxplot <- input$selected_feature
    posthoc_results_boxplot <- perform_posthoc(data(), selected_feature_boxplot, 0.05, "bonferroni")
    display_boxplot(data(), selected_feature_boxplot, posthoc_results_boxplot)
  })

  # Display individual metabolite boxplot
  output$individual_boxplot <- renderPlot({
    selected_metabolite <- input$selected_feature
    display_metabolite_boxplot(data(), selected_metabolite)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
