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
library(RColorBrewer)
library(ComplexHeatmap)
library(data.table)

# Define the UI
ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "litera"),
  titlePanel("Differential Expression Analysis"),
  tabsetPanel(
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(
                 fileInput("Browse", "Upload Sample Metadata File", accept = c(".csv",".tsv"))
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary Table ",
                            tableOutput("sample_info_summary_table")
                   ),
                   tabPanel("Metadata Table",
                            dataTableOutput("sample_info_datatable")
                   ),
                   tabPanel("Sample Distribution Plot",
                            #sidebarLayout(
                              #sidebarPanel(
                                # selectInput("sample_plot_type", "Plot type:",
                                #             choices = c("Histogram", "Density Plot", "Violin Plot"),
                                #             selected = "Histogram"),
                                selectInput("sample_plot_var", "Select Column to Plot",
                                            choices = NULL, selected = NULL),
                                selectInput("sample_group_var", "Select Column to Group by",
                                            choices = NULL, selected = NULL),
                              #),
                            plotOutput("sample_info_plot")
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
                 sliderInput("nonzero_slider", "Select the threshold for number of zero counts/gene:",
                             min = 0, max = 100,
                             value = 50, step = 1),
                 sliderInput("variance_slider", "Variance percentile:",
                             min = 0, max = 100,
                             value = 50, step = 1)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Filtering Results",
                            tableOutput("summary_filtered_genes_table")
                   ),
                   tabPanel("Scatterplot",
                            plotOutput("scatter_variance"),
                            plotOutput("scatter_zeros"),
                            switchInput("scatter_log_switch", "Log Transform", value = FALSE)
                            
                   ),
                   tabPanel("Heatmap",
                              plotOutput("counts_heatmap"),
                              switchInput("log_transform_switch", "Log Transform", value = FALSE),
                   ),
                   tabPanel("PCA", 
                              selectInput("pc_selector", "Select Principal Components", choices = NULL, multiple = TRUE),
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
    tabPanel("GSEA",
               mainPanel(
                 tabsetPanel(
                   tabPanel("GSEA Results",
                            #sidebarPanel(
                              fileInput("Browsefgsea", "Upload GSEA result file", accept = c(".csv",".tsv")),
                              sliderInput("fgsea_padj_slider", "Select Adjusted p-value threshold",
                                          min = 0, max = 1,
                                          value = 0.5, step = 0.001),
                              radioButtons("nesPathways", "Select NES Pathways:",
                                           choices = c("All", "Positive", "Negative"),
                                           selected = "All"),
                              downloadButton("fgsea_download_csv", "Download GSEA results"),
                              
                            #),
                            dataTableOutput("fgsea_datatable")
                   ),
                   tabPanel("Pathway Barplot",
                            sliderInput("NES_padj_slider", "Select Adjusted p-value threshold",
                                        min = 0, max = 1,
                                        value = 0.5, step = 0.001),
                            sliderInput("num_paths_slider", "Number of Pathways", min = 0, max = 50, value = 10, step = 2),
                            plotlyOutput("fgsea_NES_plot")
                   ),
                   tabPanel("Pathway Scatterplot", 
                            sliderInput("fgseascatter_padj_slider", "Select Adjusted p-value threshold",
                                        min = 0, max = 1,
                                        value = 0.5, step = 0.001),
                            plotOutput("fgsea_scatter_plot")
                   )
                 )
               )
             
    )
  )
)

# Define the server function
server <- function(input, output, session) {

  options(shiny.maxRequestSize = 100*1024^2)

###################
### Summary Tab ###
###################
  # Read the sample info file
  sample_info <- reactive({
    req(input$Browse)
    data <- read.csv(input$Browse$datapath)
    return(data)
  })
  
  #Create summary table for metadata 
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
      }),
      check.names = FALSE
    )
    return(summary_df)
  }
  
  # output summary table 
  output$sample_info_summary_table <- renderTable({
    summary_table(sample_info())
  })
  
  # output a DataTable of the sample info
  output$sample_info_datatable <- renderDataTable({
    sample_info()
  })
  
  # Update choices for plot choices based on the reactive sample_info
  observe({
    numeric_vars <- sapply(sample_info(), is.numeric)
    vars <- names(sample_info())[numeric_vars]
    all_vars <- names(sample_info())
    updateSelectInput(session, "sample_plot_var", choices = vars, selected = vars[1])
    updateSelectInput(session, "sample_group_var", choices = all_vars, selected = all_vars[1])
  })
  
  # Render selected plot dynamically
  output$sample_info_plot <- renderPlot({
    req(sample_info())  
    
    ggplot(sample_info(), aes_string(x = input$sample_plot_var, fill = input$sample_group_var)) +
      geom_histogram(position = "identity", bins = 30, alpha = 0.4) + 
      labs(title = paste("Histogram of", input$sample_plot_var), x = input$sample_plot_var, y = 'Frequency')
  })
  
  # Plot selected plot
  # observe({
  #   req(sample_info())  # Ensure data is available
  #   
  #   output[[paste0("plot_", input$sample_plot_var, "_", input$sample_group_var)]] <- renderPlot({
  #     ggplot(sample_info(), aes_string(x = input$sample_plot_var, fill = input$sample_group_var)) +
  #       geom_histogram(position = "identity", bins = 30, alpha = 0.7)
  #       # switch(input$sample_plot_type,
  #              # "Histogram" = geom_histogram(position = "identity", bins = 30, alpha = 0.7),
  #              # "Violin Plot" = geom_violin(alpha = 0.7),
  #              # "Density Plot" = geom_density(alpha = 0.7)
  #       #) +
  #       labs(title = paste("Histogram of", input$sample_plot_var), x = input$sample_plot_var)
  #   })
  # })
###################
### Counts Tab ###
###################
  
  # Read the normalized counts data
  normalized_counts <- reactive({
    req(input$Browsecounts)
    data <- read.csv(input$Browsecounts$datapath)
    return(data)
  })
  
 
  ###Create a table filtered by variance and zero threshold 
 # filter_counts <- function(norm_counts, var_threshold, zero_threshold) {
 #  colnames(norm_counts)[1] <- "gene"
 #  
 #  variance_and_zero_count <- apply(norm_counts[, -1], 1, function(row) {
 #    variance <- var(row, na.rm = TRUE)
 #    zero_count <- sum(row == 0, na.rm = TRUE)
 #    return(c(variance = variance, zero_count = zero_count))
 #  })
 #  
 #  variance_and_zero_count <- t(variance_and_zero_count)
 #  
 #  filtered_counts_table <- norm_counts[
 #    variance_and_zero_count[, "variance"] >= var_threshold &
 #    variance_and_zero_count[, "zero_count"] <= zero_threshold,
 #  ]
 #  
 #  return(filtered_counts_table)
 # }
 
 filter_counts <- function(norm_counts, var_percentile, zero_threshold) {
   colnames(norm_counts)[1] <- "gene"
   
   variance_and_zero_count <- apply(norm_counts[, -1], 1, function(row) {
     variance <- var(row, na.rm = TRUE)
     zero_count <- sum(row == 0, na.rm = TRUE)
     return(c(variance = variance, zero_count = zero_count))
   })
   
   variance_and_zero_count <- t(variance_and_zero_count)
   
   var_threshold <- quantile(variance_and_zero_count[, "variance"], (var_percentile/100))
   
   filtered_counts_table <- norm_counts[
     variance_and_zero_count[, "variance"] >= var_threshold &
       variance_and_zero_count[, "zero_count"] <= zero_threshold,
   ]
   
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
      "Percent Not Passing Filter" = (total_genes - passing_filter)/total_genes * 100,
      check.names = FALSE
    )
  }

  #output summary of fitlered counts data 
  output$summary_filtered_genes_table <- renderTable({
    summary_filter_counts(normalized_counts(),filter_counts(normalized_counts(),input$variance_slider,input$nonzero_slider))
    })
  
  # 
  # # Scatter plot for median count vs variance
  # output$scatter_variance <- renderPlot({
  #   norm_counts <- normalized_counts()
  #   colnames(norm_counts)[1] <- "gene"
  #   filtered_counts_table <- filter_counts(norm_counts, input$variance_slider, input$nonzero_slider)
  #   
  #   # Create a new data frame with calculated values for median and variance
  #   plot_data <- norm_counts %>%
  #     rowwise()%>%
  #     mutate(
  #       median_count = (median(c_across(-gene), na.rm = TRUE)),
  #       variance = (var(c_across(-gene), na.rm = TRUE)),
  #       status = ifelse(gene %in% filtered_counts_table$gene, "Pass", "Fail")
  #     )
  #   
  #   ggplot(plot_data, aes(x = median_count, y = variance, color = status)) +
  #     geom_point() +
  #     labs(title = "Median Count vs Variance", x = "Median Count", y = "Variance") +
  #     scale_color_manual(values = c("Pass" = "darkblue", "Fail" = "lightblue"))
  # })
  # 
  # # Scatter plot for median count vs number of zeros
  # output$scatter_zeros <- renderPlot({
  #   norm_counts <- normalized_counts()
  #   colnames(norm_counts)[1] <- "gene"
  #   filtered_counts_table <- filter_counts(norm_counts, input$variance_slider, input$nonzero_slider)
  #   
  #   # Create a new data frame with calculated values for median and number of zeros
  #   plot_data <- norm_counts %>%
  #     rowwise()%>%
  #     mutate(
  #       median_count = (median(c_across(-gene), na.rm = FALSE)),
  #       zeros = sum(c_across(-gene) == 0, na.rm = FALSE),
  #       status = ifelse(gene %in% filtered_counts_table$gene, "Pass", "Fail")
  #     )
  #   ggplot(plot_data, aes(x = median_count, y = zeros, color = status)) +
  #     geom_point() +
  #     labs(title = "Median Count vs Number of Zeros", x = "Median Count", y = "Number of Zeros") +
  #     scale_color_manual(values = c("Pass" = "darkblue", "Fail" = "lightblue"))
  # })
  
  # Scatter plot for median count vs variance
  output$scatter_variance <- renderPlot({
    norm_counts <- normalized_counts()
    colnames(norm_counts)[1] <- "gene"
    filtered_counts_table <- filter_counts(norm_counts, input$variance_slider, input$nonzero_slider)
    
    # Create a new data frame with calculated values for median and variance
    plot_data <- norm_counts %>%
      mutate(
        median_count = apply(.[-1], 1, function(x) median(x, na.rm = TRUE)),
        variance = apply(.[-1], 1, function(x) var(x, na.rm = TRUE)),
        status = ifelse(gene %in% filtered_counts_table$gene, "Pass", "Fail")
      )
    
    if (input$scatter_log_switch) {
      plot_data$variance <- log10(plot_data$variance)
      y_title <- "Log10(variance)"
    }
    else {
      y_title = 'Variance'
    }
    
    ggplot(plot_data, aes(x = median_count, y = variance, color = status)) +
      geom_point() +
      labs(title = "Median Count vs Variance", x = "Median Count", y = y_title) +
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
        median_count = apply(.[-1], 1, function(x) median(x, na.rm = FALSE)),
        zeros = apply(.[-1], 1, function(x) sum(x == 0, na.rm = FALSE)),
        status = ifelse(gene %in% filtered_counts_table$gene, "Pass", "Fail")
      )
    
    if (input$scatter_log_switch) {
      plot_data$zeros <- log10(plot_data$zeros)
      y_title <- "Log10(number of zeros)"
    }
    else {
      y_title = 'Number of Zeros'
    }
    
    ggplot(plot_data, aes(x = median_count, y = zeros, color = status)) +
      geom_point() +
      labs(title = "Median Count vs Number of Zeros", x = "Median Count", y = y_title) +
      scale_color_manual(values = c("Pass" = "darkblue", "Fail" = "lightblue"))
  })
  
  
  
  
  
  
  
  #output heatmap 
  output$counts_heatmap <- renderPlot({
    heatmap_data <- filter_counts(normalized_counts(), input$variance_slider, input$nonzero_slider)
    heatmap_data_plot <- as.matrix(heatmap_data[, -1])
    
    # Check if log transform switch is ON
    if (input$log_transform_switch) {
      # Log10 transform the data
      heatmap_data_plot <- log10(heatmap_data_plot + 1)  # Adding 1 to avoid log(0)
      legend_title = 'Log10(counts)'
    }
    else {
      legend_title = 'counts'
    }
    Heatmap(heatmap_data_plot,
            name = "Heatmap",
            col = brewer.pal(11, 'RdBu'),
            show_row_names = FALSE,
            show_column_names = TRUE,
            cluster_columns = TRUE,
            cluster_rows = TRUE,
            heatmap_legend_param = list(title = legend_title))
  })
  
  #Allow user to select the principal components that they want to plot 
  observe({
    pca_data <- filter_counts(normalized_counts(), input$variance_slider, input$nonzero_slider)
    pca <- prcomp(t(pca_data[-1]))
    updateSelectInput(session, "pc_selector", choices = colnames(pca$x), selected = colnames(pca$x)[1:2])
  })
  
  #output pca plot of selected PCs 
  output$pca_plot <- renderPlot({
    pca_data <- filter_counts(normalized_counts(), input$variance_slider, input$nonzero_slider)
    pca <- prcomp(t(pca_data[-1]))
    
    #if (input$plot_type == "Scatter Plot") {
      selected_pcs <- input$pc_selector[1:2]
      
      # Scatter plot for selected PCs
      pca_scores <- as.data.frame(pca$x[, selected_pcs])
      ggplot(pca_scores, aes_string(x = selected_pcs[1], y = selected_pcs[2])) +
        geom_point() +
        labs(x = paste(selected_pcs[1], "(", round(100 * pca$sdev[which(colnames(pca$x) == selected_pcs[1])]^2 / sum(pca$sdev^2), 0), "% variance)"),
             y = paste(selected_pcs[2], "(", round(100 * pca$sdev[which(colnames(pca$x) == selected_pcs[2])]^2 / sum(pca$sdev^2), 0), "% variance)"),
             title = 'PCA')
    # } else {
    #   # Beeswarm plot for top N PCs
    #   top_n_pcs <- paste0("PC", seq((ncol(pca$x) - input$top_n + 1), ncol(pca$x)))
    #   pca_scores <- as.data.frame(pca$x[, top_n_pcs])
    # 
    #   # Create a column for x-axis variable
    #   pca_scores$x_axis <- factor(top_n_pcs)
    # 
    #   ggplot(pca_scores, aes(x = x_axis, y = pca_scores[[1]])) +
    #     geom_beeswarm() +
    #     labs(x = "Top Principal Components", y = "", title = 'Top N PCA')
    # }
  })
  
  ###################
  ### DE Tab ###
  ###################
  
  #Read differential expression data 
  diff_exp_data <- reactive({
    req(input$BrowseDE)
    data <- read.csv(input$BrowseDE$datapath)
    return(data)
  })
  #output differential expression data 
  output$DE_datatable <- renderDataTable({
    diff_exp_data()
  })
  
  #create dataframe for volcano plot 
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
  
  #plot p value historgram 
  plot_pvals <- function(labeled_results) {
    pvalue_histogram <- ggplot(labeled_results, aes(x = pvalue)) +
      geom_histogram(binwidth = 0.02, fill = "lightblue", color = "darkblue") +  # Adjust binwidth and colors as needed
      labs(x = "p-value",
           y = "Frequency") +
      theme_minimal()

    return(pvalue_histogram)
  }
  
  #plot for log2fc histogram 
  plot_log2fc <- function(labeled_results, padj_threshold) {
    plot_data <- labeled_results %>%
      dplyr::filter(padj < padj_threshold)
    
    log2fc_histogram <- ggplot(plot_data, aes(x = log2FoldChange)) +
      geom_histogram(binwidth = 0.2, fill = "lightblue",color = "darkblue") +  # Adjust binwidth and colors as needed
      labs(title = "Histogram of Log2FoldChanges for DE Genes",
           x = "Log2FoldChange",
           y = "Frequency") +
      theme_minimal()
    
    return(log2fc_histogram)
  }
  
  #volcano plot 
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
    diff_exp_gsea <- as.data.frame(diff_exp_df) 
    
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
    pathways_hallmark <- fgsea::gmtPathways('c2.cp.v2023.2.Hs.symbols.gmt')
    fgsea_results <- fgsea::fgsea(pathways_hallmark, diff_exp_gsea_trim, minSize = 15, maxSize = 500)
    
    # Convert results to a data frame if needed
    fgsea_results_df <- as.data.frame(fgsea_results)
    
    return(fgsea_results_df)
  }
  
  
  
  #filter fgsea table for plots and table 
  filtered_fgsea <- function(fgsea, padj_threshold, pathway_selection) {
    filtered_table <- fgsea%>%
      filter(padj < padj_threshold)
    
    if (pathway_selection == "All") {
      return(filtered_table)
    } else if (pathway_selection == "Positive") {
      return(filtered_table %>% filter(NES > 0))
    } else if (pathway_selection == "Negative") {
      return(filtered_table %>% filter(NES < 0))
    }
    
    return(filtered_table)
    
  }
  
  ###plot top pathways as barplot with interactive bars 
  top_pathways <- function(fgsea_results, num_paths) {
    
    plt_top <- fgsea_results %>%
      arrange(NES)%>%
      head(num_paths/2)
    plt_bottom <- fgsea_results %>%
      arrange(NES) %>%
      tail(num_paths/2)
    
    plt_data <- rbind(plt_top,plt_bottom)
    
    plt_data$color <- ifelse(plt_data$NES < 0, "Negative", "Positive") 
    # plt_data <- plt_data %>%
    #   arrange(NES)
    
    plt <- plot_ly(data=plt_data, x=~NES, y=~reorder(pathway, +NES), type="bar", marker=list(color=ifelse(plt_data$color == "Negative", "red", "blue"))) %>%
      layout(title="fgsea results for Canonical MSigDB gene set",
             xaxis=list(title="NES"),
             yaxis=list(title=""))
  
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
      scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "darkblue")) +  # Adjust color for significance
      labs(x = "NES", y = "-log10(Adjusted P-value)",
           title = "Scatter plot of NES vs -log10 adjusted p-value") +
      labs(color = "Significant Pathway")+
      theme_minimal()
    return(fgsea_scatter_plot)
  }
  
  
  output$fgsea_datatable <- renderDataTable({
    filtered_fgsea(fgsea_table(diff_exp_data()),input$fgsea_padj_slider,input$nesPathways)
  })
  
  output$fgsea_download_csv <- downloadHandler(
    
    filename = function() {
      paste("fgsea_results_",".csv", sep = "")
    },
    content = function(file) {
      fgsea_table <- filtered_fgsea(fgsea_table(diff_exp_data()),input$fgsea_padj_slider,input$nesPathways)
      #fgsea_table <- as.data.frame(fgsea_table)
      fwrite(res, file=filename, sep=",", sep2=c("", " ", ""))
    }
  )
  
  
 output$fgsea_NES_plot <- renderPlotly({
    num_paths <- input$num_paths_slider
    top_pathways(filtered_fgsea(fgsea_table(diff_exp_data()),input$NES_padj_slider,"All"), num_paths)
  })
  
  output$fgsea_scatter_plot <- renderPlot({
    fgsea_scatter(fgsea_table(diff_exp_data()),input$fgseascatter_padj_slider)
  })
  
}

# Run the app
shinyApp(ui, server)
