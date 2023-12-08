library(dplyr)
library(ggplot2)
library(bslib)
library(shiny)
library(colourpicker)

# Define the UI
ui <- fluidPage(
  titlePanel("Differential Expression Analysis"),
  tabsetPanel(
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(
                 fileInput("Browse", "Sample file", accept = ".csv")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            tableOutput("sample_info_summary_table")
                   ),
                   tabPanel("Table",
                            dataTableOutput("sample_info_datatable")
                   ),
                   tabPanel("Plots",
                            uiOutput("sample_info_plot")
                   )
                 )
               )
             )
    ),
    tabPanel("Counts",
             sidebarLayout(
               sidebarPanel(
                 fileInput("Browse", "Normalized Counts Data", accept = ".csv")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary for Tab2",
                            tableOutput("summary_table_tab2")
                   ),
                   tabPanel("Table for Tab2",
                            dataTableOutput("datatable_tab2")
                   ),
                   tabPanel("Plots for Tab2"
                            #, plotOutput("plot_tab2")
                   )
                 )
               )
             )
    ),
    tabPanel("DE",
             sidebarLayout(
               sidebarPanel(
                 fileInput("BrowseTab3", "Sample file for Tab3", accept = ".csv")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary for Tab3",
                            tableOutput("summary_table_tab3")
                   ),
                   tabPanel("Table for Tab3",
                            dataTableOutput("datatable_tab3")
                   ),
                   tabPanel("Plots for Tab3"
                            #, plotOutput("plot_tab3")
                   )
                 )
               )
             )
    ),
    tabPanel("Network",
             sidebarLayout(
               sidebarPanel(
                 fileInput("BrowseTab4", "Sample file for Tab4", accept = ".csv")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary for Tab4",
                            tableOutput("summary_table_tab4")
                   ),
                   tabPanel("Table for Tab4",
                            dataTableOutput("datatable_tab4")
                   ),
                   tabPanel("Plots for Tab4"
                            #, plotOutput("plot_tab4")
                   )
                 )
               )
             )
    )
  )
)

# Define the server function
server <- function(input, output, session) {
  
  # Read the sample info file
  sample_info <- reactive({
    req(input$Browse)
    data <- read.csv(input$Browse$datapath)
    return(data)
  })
  
  summary_table <- function(sample_data) {
    summary_df <- data.frame(
      'Column Name' = names(sample_data),
      'Type' = sapply(sample_data, function(x) class(x)[1]),
      'Mean(sd) or Distinct Values' = sapply(sample_data, function(x) {
        if (is.numeric(x)) {
          sprintf("%.2f (+/- %.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
        } else {
          x_vector <- as.vector(x)
          toString(unique(x_vector))
        }
      })
    )
    return(summary_df)
  }
  
  # Create a table of the sample info
  output$sample_info_summary_table <- renderTable({
    summary_table(sample_info())
  })
  
  # Create a DataTable of the sample info
  output$sample_info_datatable <- renderDataTable({
    sample_info()
  })
  
  # Render histograms dynamically
  output$sample_info_plot <- renderUI({
    req(names(sample_info()))  # Ensure data is available
    
    numeric_vars <- sapply(sample_info(), is.numeric)
    numeric_var_names <- names(sample_info())[numeric_vars]
    
    plots <- lapply(numeric_var_names, function(var_name) {
      hist_output_id <- paste0("hist_", var_name)
      plotOutput(hist_output_id)
    })
    
    tagList(plots)
  })
  
  # Plot histograms
  observe({
    req(names(sample_info()))  # Ensure data is available
    
    numeric_vars <- sapply(sample_info(), is.numeric)
    numeric_var_names <- names(sample_info())[numeric_vars]
    
    plots <- lapply(numeric_var_names, function(var_name) {
      hist_output_id <- paste0("hist_", var_name)
      output[[hist_output_id]] <- renderPlot({
        hist(sample_info()[[var_name]], main = paste("Histogram for", var_name),
             xlab = var_name, col = "lightblue", border = "black")
      })
      plotOutput(hist_output_id)
    })
    
    output$sample_info_plot <- renderUI({
      tagList(plots)
    })
  })
}

# Run the app
shinyApp(ui, server)
