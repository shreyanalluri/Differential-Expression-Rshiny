library(dplyr)
library(tidyr)
library(matrixStats)
library(ggplot2)
library(bslib)
library(shiny)
library(colourpicker)

# Define the UI
ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "litera"),
  titlePanel("Differential Expression Analysis"),
  tabsetPanel(
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(
                 fileInput("Browse", "Sample file", accept = c(".csv",".tsv"))
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
                 fileInput("Browsecounts", "Normalized Counts Data", accept = ".csv"),
                 # Input: Simple integer interval ----
                 sliderInput("nonzero_slider", "Select the threshold for number of non-zero samples/gene:",
                             min = 0, max = 1000,
                             value = 500),
                 
                 # Input: Decimal interval with step value ----
                 sliderInput("variance_slider", "Variance percentile:",
                             min = 0, max = 100,
                             value = 50, step = 1)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            tableOutput("counts_filter_summary_table")
                   ),
                   tabPanel("Scatterplot",
                            plotOutput("filtered_genes_plot")
                   ),
                   tabPanel("Heatmap"
                            #, plotOutput("plot_tab3")
                   ),
                   tabPanel("PCA"
                            #, plotOutput("plot_tab4")
                   )
                 )
               )
             )
    ),
    tabPanel("DE",
             sidebarLayout(
               sidebarPanel(
                 fileInput("BrowseDE", "Differential expression data", accept = c(".csv",".tsv")),
                 sliderInput("pval_slider", "p-value Threshold",
                             min = 0, max = 1,
                             value = 0.5,step = 0.001),
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Differential Expression Data",
                            dataTableOutput("DE_datatable")
                   ),
                   tabPanel("p-value Plot"
                            , plotOutput("DE_pval_plot")
                   ),
                   tabPanel("Log2 FC Plot"
                            , plotOutput("DE_log2fc_plot")
                   ),
                   tabPanel("Volcano Plot"
                            , plotOutput("DE_volcano_plot")
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

###################
### Summary Tab ###
###################
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
        hist(sample_info()[[var_name]], main = paste(var_name, "Histogram"),
             xlab = var_name, col = "lightblue", border = "black")
      })
      plotOutput(hist_output_id)
    })
    
    output$sample_info_plot <- renderUI({
      tagList(plots)
    })
  })
  
  
###################
### Counts Tab ###
###################
  
  # Read the normalized counts data
  normalized_counts <- reactive({
    req(input$Browsecounts)
    data <- read.csv(input$Browsecounts$datapath)
    return(data)
  })
  
  filtered_genes <- reactive({
    # Filter genes based on variance and non-zero samples
    variance_filtered <- apply(normalized_counts(), 1, function(row) var(row) >= input$variance_slider)
    nonzero_samples_filtered <- apply(normalized_counts() != 0, 1, function(row) sum(row) >= input$nonzero_slider)
    
    selected_genes <- rownames(normalized_counts())[variance_filtered & nonzero_samples_filtered]
    
    # Create a table for display
    filtered_genes_table <- data.frame(
      Gene = selected_genes,
      Variance = apply(normalized_counts()[selected_genes, ], 1, var),
      NonZeroSamples = apply(normalized_counts()[selected_genes, ] != 0, 1, sum)
    )
    
    filtered_genes_table
  })
  
  output$filtered_genes_table <- renderDataTable({
    filtered_genes()
  })
  
  observeEvent(c(input$variance_slider, input$nonzero_slider), {
    output$counts_filter_summary_table <- renderTable({
      data <- normalized_counts()
      total_genes <- nrow(data)
      total_samples <- ncol(data)
      filtered_genes_count <- nrow(filtered_genes())
      not_filtered_genes_count <- total_genes - filtered_genes_count
      
      cbind(
        "Total Samples" = total_samples,
        "Total Genes" = total_genes,
        "Genes Passing Filter" = filtered_genes_count,
        "Genes Not Passing Filter" = not_filtered_genes_count,
        "% Passing Filter" = sprintf("%.2f%%", 100 * filtered_genes_count / total_genes),
        "% Not Passing Filter" = sprintf("%.2f%%", 100 * not_filtered_genes_count / total_genes)
      )
    })
    
    
    ####For this table, I need the dataframe to be gene:median count, variance, nzeroes 
    # output$filtered_genes_plot <- renderPlot({
    #   variance_filtered <- apply(normalized_counts(), 1, function(row) var(row) >= input$variance_slider)
    #   nonzero_samples_filtered <- apply(normalized_counts() != 0, 1, function(row) sum(row) >= input$nonzero_slider)
    # 
    #   filtered_plot <- normalized_counts()[variance_filtered & nonzero_samples_filtered, ]
    # 
    #   gene_names <- colnames(filtered_plot)
    # 
    #   plot_tibble <- tibble(gene_name = gene_names, !!filtered_plot) %>%
    #     pivot_longer(cols = -gene_name, names_to = "sample_name", values_to = "count") %>%
    #     group_by(gene_name) %>%
    #     summarize(
    #       median_count = median(count),
    #       num_zeroes = sum(count == 0),
    #       variance = colVars(as.matrix(count))[1]
    #     )
    # 
    # })
    
    # output$filtered_genes_plot <- renderPlot({
    #   variance_filtered <- apply(normalized_counts(), 1, function(row) var(row) >= input$variance_slider)
    #   nonzero_samples_filtered <- apply(normalized_counts() != 0, 1, function(row) sum(row) >= input$nonzero_slider)
    #   
    #   filtered_plot <- normalized_counts()[variance_filtered & nonzero_samples_filtered, ]
    #   
    #   gene_names <- rownames(filtered_plot)
    #   
    #   # Convert all columns to numeric
    #   filtered_plot <- apply(filtered_plot, 2, as.numeric)
    #   
    #   plot_tibble <- tibble(gene_name = gene_names, !!filtered_plot) %>%
    #     pivot_longer(cols = -gene_name, names_to = "sample_name", values_to = "count") %>%
    #     group_by(gene_name) %>%
    #     summarize(
    #       median_count = median(count, na.rm = TRUE),
    #       num_zeroes = sum(count == 0, na.rm = TRUE),
    #       variance = colVars(as.matrix(count))[1, na.rm = TRUE]
    #     )
    #   
    #   # Scatterplot of median count vs variance
    #   ggplot(plot_tibble, aes(x = median_count, y = variance)) +
    #     geom_point(color = "darkblue", alpha = 0.7) +
    #     labs(title = "Scatterplot of Median Count vs Variance",
    #          x = "Median Count",
    #          y = "Variance") +
    #     theme_minimal()
    # })
    
    
    
  })
  
  
  ###################
  ### DE Tab ###
  ###################
  
  diff_exp_data <- reactive({
    req(input$BrowseDE)
    data <- read.csv(input$BrowseDE$datapath)
    return(data)
  })
  
  output$DE_datatable <- renderDataTable({
    diff_exp_data()
  })
  

  label_res <- function(deseq2_res, padj_threshold) {
    # Create a tibble from DESeq2 results
    deseq2_tibble <- as_tibble(deseq2_res)
    
    # Add volc_plot_status column based on your criteria using dplyr
    deseq2_tibble <- deseq2_tibble %>%
      mutate(
        volc_plot_status = ifelse(
          padj < padj_threshold & log2FoldChange > 0, "UP",
          ifelse(
            padj < padj_threshold & log2FoldChange < 0, "DOWN", "NS"
          )
        )
      )
    
    return(deseq2_tibble)
  }
  
  plot_pvals <- function(labeled_results) {
    pvalue_histogram <- ggplot(labeled_results, aes(x = pvalue)) +
      geom_histogram(binwidth = 0.02, fill = "pink", color = "magenta") +  # Adjust binwidth and colors as needed
      labs(x = "pvalue",
           y = "count") +
      theme_minimal()
    
    # Print the histogram
    return(pvalue_histogram)
  }
  
  plot_log2fc <- function(labeled_results, padj_threshold) {
    plot_data <- labeled_results %>%
      dplyr::filter(padj < padj_threshold)
    
    log2fc_histogram <- ggplot(plot_data, aes(x = log2FoldChange)) +
      geom_histogram(binwidth = 0.2, fill = "pink",color = "magenta") +  # Adjust binwidth and colors as needed
      labs(title = "Histogram of Log2FoldChanges for DE Genes",
           x = "log2FoldChange",
           y = "count") +
      theme_minimal()
    
    return(log2fc_histogram)
  }
  
  plot_volcano <- function(labeled_results) {
    plot_data <- labeled_results %>%
      drop_na(volc_plot_status)
    volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = volc_plot_status)) +
      labs(x = "log2FoldChange", y = "-log10(padj)") +
      ggtitle("Volcano plot of differential expression results") +
      geom_hline(yintercept= 0, linetype = 'dashed', col="black") +
      theme_minimal()
    
    return(volcano_plot)
  }
  
  observeEvent(input$pval_slider, {
    
    output$DE_pval_plot <- renderPlot({
      plot_pvals(label_res(diff_exp_data(), input$pval_slider))
    })
    
    output$DE_log2fc_plot <- renderPlot({
      plot_log2fc(label_res(diff_exp_data(), input$pval_slider), input$pval_slider)
    })
    
    output$DE_volcano_plot <- renderPlot({
      plot_volcano(label_res(diff_exp_data(), input$pval_slider))
      
    })
    
  })
}

# Run the app
shinyApp(ui, server)
