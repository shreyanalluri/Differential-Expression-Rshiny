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
                 fileInput("Browsecounts", "Normalized Counts Data", accept = c(".csv",".tsv")),
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
                   tabPanel("Table",
                            tableOutput("summary_filtered_genes_table")
                   ),
                   tabPanel("Scatterplot",
                            plotOutput("scatter_variance"),
                            plotOutput("scatter_zeros")
                            
                   ),
                   tabPanel("Heatmap"
                            , plotOutput("counts_heatmap")
                   ),
                   tabPanel("PCA", 
                            plotOutput("counts_pca")
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
  
 
  ###Create a filtered table 
  filter_counts <- function(norm_counts, var_threshold, zero_threshold) {
    colnames(norm_counts)[1] <- "gene"
    
    filtered_counts_table <- norm_counts %>%
      rowwise() %>%
      mutate(
        variance = var(c_across(-gene), na.rm = TRUE),
        zero_count = sum(c_across(-gene) == 0, na.rm = TRUE)
      ) %>%
      filter(
        variance > var_threshold,
        zero_count <= zero_threshold
      ) %>%
      select(-variance, -zero_count)
    return(filtered_counts_table)
  }
  ### create summary of filtered table 
  summary_filter_counts <- function(norm_counts, filtered_counts) {
    total_genes <- nrow(norm_counts)
    total_samples <- ncol(norm_counts) - 1
    passing_filter <- nrow(filtered_counts)
    
    summary_table <- data.frame(
      "Total Samples" = total_samples,
      "Total Genes" = total_genes,
      "Passing Filter Genes" = passing_filter,
      "Percent Passing Filter" = passing_filter/total_genes * 100,
      "Not Passing Filter Genes" = total_genes - passing_filter,
      "Percent Not Passing Filter" = (total_genes - passing_filter)/total_genes * 100
    )
    return(summary_table)
  }

  
  output$summary_filtered_genes_table <- renderTable({
    summary_filter_counts(normalized_counts(),filter_counts(normalized_counts(),input$variance_slider,input$nonzero_slider))
    })
  
  output$counts_med_v_var_plot <- renderPlot({
    plot_median_vs_variance(filter_counts(normalized_counts(),input$variance_slider,input$nonzero_slider))
  })
  
  # Scatter plot for median count vs variance
  output$scatter_variance <- renderPlot({
    norm_counts <- normalized_counts()
    colnames(norm_counts)[1] <- "gene"
    filtered_counts_table <- filter_counts(norm_counts, input$variance_slider, input$nonzero_slider)
    
    # Create a new data frame with calculated values for median and variance
    plot_data <- norm_counts %>%
      mutate(
        median_count = (median(c_across(-gene), na.rm = FALSE)),
        variance = (var(c_across(-gene), na.rm = FALSE)),
        status = ifelse(gene %in% filtered_counts_table$gene, "Pass", "Fail")
      )
    
    ggplot(plot_data, aes(x = median_count, y = variance, color = status)) +
      geom_point() +
      labs(title = "Median Count vs Variance", x = "Median Count", y = "Variance") +
      scale_color_manual(values = c("Pass" = "darkblue", "Fail" = "lightblue"))
  })
  
  # Scatter plot for median count vs number of zeros
  output$scatter_zeros <- renderPlot({
    norm_counts <- normalized_counts()
    colnames(norm_counts)[1] <- "gene"
    filtered_counts_table <- filter_counts(norm_counts, input$variance_slider, input$nonzero_slider)
    
    # Create a new data frame with calculated values for median and number of zeros
    plot_data <- norm_counts %>%
      mutate(
        median_count = (median(c_across(-gene), na.rm = FALSE)),
        zeros = sum(c_across(-gene) == 0, na.rm = FALSE),
        status = ifelse(gene %in% filtered_counts_table$gene, "Pass", "Fail")
      )
    
    ggplot(plot_data, aes(x = median_count, y = zeros, color = status)) +
      geom_point() +
      labs(title = "Median Count vs Number of Zeros", x = "Median Count", y = "Number of Zeros") +
      scale_color_manual(values = c("Pass" = "darkblue", "Fail" = "lightblue"))
  })
  
  output$counts_heatmap <- renderPlot({
    heatmap_data <- filter_counts(normalized_counts(),input$variance_slider,input$nonzero_slider)
    heatmap_data <- as.matrix(heatmap_data[-1, -1])
    heatmap(heatmap_data)
  })
  
  output$counts_pca <- renderPlot({
    pca_data <- filter_counts(normalized_counts(),input$variance_slider,input$nonzero_slider)
    pca <- prcomp(t(pca_data[-1]))
    
    pc1_variance <- round(100 * pca$sdev[1]^2 / sum(pca$sdev^2), 0)
    pc2_variance <- round(100 * pca$sdev[2]^2 / sum(pca$sdev^2), 0)
    
    # Create axis labels that include the variance information
    x_label <- paste("PC1: ", pc1_variance, "% variance", sep = "")
    y_label <- paste("PC2: ", pc2_variance, "% variance", sep = "")
    
    pca_scores <- as.data.frame(pca$x)
    pca_graph <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
      geom_point() +
      labs(x = x_label, y = y_label, title = 'PCA')
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
