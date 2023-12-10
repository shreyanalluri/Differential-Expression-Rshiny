library(bslib)
library(colourpicker)
library(tidyverse)
library(dplyr)
library(tidyr)
library(matrixStats)
library(ggplot2)
library(plotly)
library(ggbeeswarm)
library(shiny)
library(biomaRt)
library('fgsea')
library(shinyWidgets)

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
                            #sidebarLayout(
                              #sidebarPanel(
                                selectInput("sample_plot_type", "Plot type:",
                                            choices = c("Histogram", "Density Plot", "Violin Plot"),
                                            selected = "Histogram"),
                                selectInput("sample_plot_var", "Variable to Plot",
                                            choices = NULL, selected = NULL),
                                selectInput("sample_group_var", "Variable to Group by",
                                            choices = NULL, selected = NULL),
                              #),
                            uiOutput("sample_info_plot")
                     #)
                   )
                 )
               )
             )
    ),
    tabPanel("Counts",
             sidebarLayout(
               sidebarPanel(
                 fileInput("Browsecounts", "Normalized Counts Data", accept = c(".csv",".tsv")),
                 sliderInput("nonzero_slider", "Select the threshold for number of non-zero samples/gene:",
                             min = 0, max = 1000,
                             value = 500),
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
                   tabPanel("Heatmap",
                              switchInput("log_transform_switch", "Log Transform", value = FALSE),
                              plotOutput("counts_heatmap")
                   ),
                   tabPanel("PCA", 
                            selectInput("plot_type", "Select Plot Type", choices = c("Scatter Plot", "Beeswarm Plot"), selected = "Scatter Plot"),
                            conditionalPanel(
                              condition = "input.plot_type == 'Scatter Plot'",
                              selectInput("pc_selector", "Select Principal Components", choices = NULL, multiple = TRUE)
                            ),
                            conditionalPanel(
                              condition = "input.plot_type == 'Beeswarm Plot'",
                              numericInput("top_n", "Top N Principal Components", value = 2, min = 1)
                            ),
                            plotOutput("pca_plot")
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
                 fileInput("Browsefgsea", "fgsea file", accept = c(".csv",".tsv"))
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("fgsea results",
                            dataTableOutput("fgsea_datatable")
                   ),
                   tabPanel("NES plot",
                            plotlyOutput("fgsea_NES_plot")
                   ),
                   tabPanel("Scatterplot", 
                            plotOutput("fgsea_scatter_plot")
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
  
  # Update choices for selectInput based on the reactive sample_info
  observe({
    numeric_vars <- sapply(sample_info(), is.numeric)
    vars <- names(sample_info())[numeric_vars]
    all_vars <- names(sample_info())
    updateSelectInput(session, "sample_plot_var", choices = vars, selected = vars[1])
    updateSelectInput(session, "sample_group_var", choices = all_vars, selected = all_vars[1])
  })
  
  # Render selected plot dynamically
  output$sample_info_plot <- renderUI({
    req(sample_info())  # Ensure data is available
    
    plot_output_id <- paste0("plot_", input$sample_plot_var, "_", input$sample_group_var)
    plotOutput(plot_output_id)
  })
  
  # Plot selected plot
  observe({
    req(sample_info())  # Ensure data is available
    
    output[[paste0("plot_", input$sample_plot_var, "_", input$sample_group_var)]] <- renderPlot({
      ggplot(sample_info(), aes_string(x = input$sample_plot_var, fill = input$sample_group_var)) +
        switch(input$sample_plot_type,
               "Histogram" = geom_histogram(position = "identity", bins = 30, alpha = 0.7),
               "Violin Plot" = geom_violin(alpha = 0.7),
               "Density Plot" = geom_density(alpha = 0.7)
        ) +
        labs(title = paste(input$sample_plot_type, "of", input$sample_plot_var), x = input$sample_plot_var, y = ifelse(input$sample_plot_type == "Density Plot", "Density", "Count"))
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
      dplyr::rowwise() %>%
      mutate(
        variance = var(c_across(-gene), na.rm = TRUE),
        zero_count = sum(c_across(-gene) == 0, na.rm = TRUE)
      ) %>%
      filter(
        variance > var_threshold,
        zero_count <= zero_threshold
      ) %>%
      dplyr::select(-variance, -zero_count)
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
    heatmap_data <- filter_counts(normalized_counts(), input$variance_slider, input$nonzero_slider)
    heatmap_data <- as.matrix(heatmap_data[-1, -1])
    
    # Check if log transform switch is ON
    if (input$log_transform_switch) {
      # Log10 transform the data
      heatmap_data <- log10(heatmap_data + 1)  # Adding 1 to avoid log(0)
    }
    
    heatmap(heatmap_data)
  })
  
  observe({
    pca_data <- filter_counts(normalized_counts(), input$variance_slider, input$nonzero_slider)
    pca <- prcomp(t(pca_data[-1]))
    updateSelectInput(session, "pc_selector", choices = colnames(pca$x), selected = colnames(pca$x)[1:2])
  })
  
  output$pca_plot <- renderPlot({
    pca_data <- filter_counts(normalized_counts(), input$variance_slider, input$nonzero_slider)
    pca <- prcomp(t(pca_data[-1]))
    
    if (input$plot_type == "Scatter Plot") {
      selected_pcs <- input$pc_selector[1:2]
      
      # Scatter plot for selected PCs
      pca_scores <- as.data.frame(pca$x[, selected_pcs])
      ggplot(pca_scores, aes_string(x = selected_pcs[1], y = selected_pcs[2])) +
        geom_point() +
        labs(x = paste(selected_pcs[1], "(", round(100 * pca$sdev[which(colnames(pca$x) == selected_pcs[1])]^2 / sum(pca$sdev^2), 0), "% variance)"),
             y = paste(selected_pcs[2], "(", round(100 * pca$sdev[which(colnames(pca$x) == selected_pcs[2])]^2 / sum(pca$sdev^2), 0), "% variance)"),
             title = 'PCA')
    } else {
      # Beeswarm plot for top N PCs
      # Beeswarm plot for top N PCs
      top_n_pcs <- paste0("PC", seq_len(min(input$top_n, ncol(pca$x))))
      pca_scores <- as.data.frame(pca$x[, top_n_pcs])
      
      ggplot(pca_scores, aes(x = factor(1), y = pca_scores[[1]])) +
        geom_beeswarm() +
        labs(x = "Top Principal Components", y = "", title = 'Top N PCA')
    }
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
  
  ###################
  ###FGSEA Tab ###
  ###################
  
  fgsea_table <- function(diff_exp_df) {
    # Map Ensembl gene ids to HGNC symbols
    diff_exp_gsea <- as.data.frame(diff_exp_df)  # Convert to data frame if needed
    #ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    # gene_map <- as.data.frame(
    #   getBM(
    #     attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
    #     mart = ensembl
    #   )
    # )
    
    #diff_exp_gsea <- left_join(diff_exp_gsea, gene_map, by = c("row" = "ensembl_gene_id"))
    
    # Simplify table for analysis, rank, and turn into a list for fgsea
    diff_exp_gsea_trim <- diff_exp_gsea %>% 
      dplyr::select(symbol, log2FoldChange) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(symbol) %>% 
      summarize(log2FoldChange = mean(log2FoldChange)) %>%
      arrange(desc(log2FoldChange)) %>%
      deframe()
    
    # Perform fgsea
    pathways_hallmark <- fgsea::gmtPathways('h.all.v2023.2.Hs.symbols.gmt')
    fgsea_results <- fgsea::fgsea(pathways_hallmark, diff_exp_gsea_trim, minSize = 15, maxSize = 500)
    
    # Convert results to a data frame if needed
    fgsea_results <- as.data.frame(fgsea_results)
    
    return(fgsea_results)
  }
  
  ###plot top pathways as barplot 
  top_pathways <- function(fgsea_results, num_paths) {
    
    # barplot() function is used to
    # plot the bar and horiz field is
    # used to plot bar horizontally
    
    plt_top <- fgsea_results %>%
      arrange(NES)%>%
      head(num_paths)
    plt_bottom <- fgsea_results %>%
      arrange(NES) %>%
      tail(num_paths)
    
    plt_data <- rbind(plt_top,plt_bottom)
    
    plt_data$color <- ifelse(plt_data$NES < 0, "Negative", "Positive") 
    plt_data <- plt_data %>%
      arrange(NES)
    
    plt <- plot_ly(data=plt_data, x=~NES, y=~reorder(pathway, +NES), type="bar", marker=list(color=~color)) %>%
      layout(title="fgsea results for Hallmark MSigDB gene set",
             xaxis=list(title="Normalized Enrichment Score (NES)"),
             yaxis=list(title=""))
    
    # plt <-ggplot(data=plt_data, aes(x=NES, y=reorder(pathway,+NES),fill = color)) +
    #   geom_bar(stat="identity")+
    #   theme(axis.text = element_text(size = 4),  # Adjust the font size of axis labels
    #         axis.title = element_text(size = 6)) + 
    #   guides(fill = 'none')+
    #   labs(x = "Normalized Enrichment Score (NES)", y = '',title = "fgsea results for Hallmark MSigDB gene set")
    
    
    return(plt)
  }
  
  ###Scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis, with gene sets below threshold in grey color
  fgsea_scatter <- function(fgsea_results, padj_threshold) {
    significant_gene_sets <- fgsea_results$pathway[fgsea_results$padj < padj_threshold]
    
    # Scatter plot
    fgsea_scatter_plot <- ggplot(fgsea_results, 
                                 aes(x = NES, y = -log10(padj), 
                                     color = pathway %in% significant_gene_sets)) +
      geom_point() +
      scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "red")) +  # Adjust color for significance
      labs(x = "Normalized Enrichment Score (NES)", y = "-log10(Adjusted P-value)",
           title = "Scatter plot of NES vs -log10 adjusted p-value") +
      theme_minimal()
    return(fgsea_scatter_plot)
  }
  
  
  output$fgsea_datatable <- renderDataTable({
    fgsea_table(diff_exp_data())
  })
  
  output$fgsea_NES_plot <- renderPlotly({
    top_pathways(fgsea_table(diff_exp_data()),6)
  })
  
  output$fgsea_scatter_plot <- renderPlot({
    fgsea_scatter(fgsea_table(diff_exp_data()),0.05)
  })
  
}

# Run the app
shinyApp(ui, server)
