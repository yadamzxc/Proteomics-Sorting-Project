
> # ============================================================================
> # CRISPR Screen Cross-Reference Shiny App
> # ============================================================================
> 
> library(shiny)
Warning message:
package ‘shiny’ was built under R version 4.5.2 
> library(readr)
Warning message:
package ‘readr’ was built under R version 4.5.2 
> library(readxl)
Warning message:
package ‘readxl’ was built under R version 4.5.2 
> library(dplyr)

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Warning message:
package ‘dplyr’ was built under R version 4.5.2 
> library(stringr)
Warning message:
package ‘stringr’ was built under R version 4.5.2 
> library(DT)

Attaching package: ‘DT’

The following objects are masked from ‘package:shiny’:

    dataTableOutput, renderDataTable

Warning message:
package ‘DT’ was built under R version 4.5.2 
> library(ggplot2)
Warning message:
package ‘ggplot2’ was built under R version 4.5.2 
> library(tidyr)
Warning message:
package ‘tidyr’ was built under R version 4.5.2 
> library(ComplexUpset)
Warning message:
package ‘ComplexUpset’ was built under R version 4.5.2 
> 
> # ============================================================================
> # HELPER FUNCTIONS
> # ============================================================================
> 
> # Clean gene IDs: uppercase, trim, remove NA/blank, deduplicate
> clean_gene_ids <- function(genes) {
+   genes_char <- as.character(genes)
+   genes_clean <- str_trim(genes_char)
+   genes_upper <- str_to_upper(genes_clean)
+   genes_nonempty <- genes_upper[!is.na(genes_upper) & genes_upper != ""]
+   unique(genes_nonempty)
+ }
> 
> # Load and clean a gene list file
> load_gene_file <- function(file_path) {
+   ext <- tools::file_ext(file_path)
+   
+   if (ext %in% c("csv", "txt")) {
+     df <- read_csv(file_path, show_col_types = FALSE)
+   } else if (ext %in% c("xlsx", "xls")) {
+     df <- read_excel(file_path)
+   } else {
+     stop("Unsupported file format. Use CSV or XLSX.")
+   }
+   
+   # Check for gene_id column (case-insensitive)
+   col_names <- names(df)
+   gene_col_idx <- which(str_to_lower(col_names) == "gene_id")
+   
+   # If no gene_id column found, use the FIRST column
+   if (length(gene_col_idx) == 0) {
+     gene_col_idx <- 1
+   }
+   
+   # Extract and clean gene column
+   genes_raw <- df[[gene_col_idx[1]]]
+   genes_clean <- clean_gene_ids(genes_raw)
+   
+   # Create new dataframe
+   df_clean <- data.frame(gene_id = genes_clean, stringsAsFactors = FALSE)
+   
+   # Add back other columns if they exist
+   if (ncol(df) > 1) {
+     other_cols <- df[!duplicated(df[[gene_col_idx[1]]]), -gene_col_idx[1], drop = FALSE]
+     if (nrow(other_cols) == nrow(df_clean)) {
+       df_clean <- cbind(df_clean, other_cols)
+     }
+   }
+   
+   return(df_clean)
+ }
> 
> # Load ranked screen (detects positive/negative columns)
> load_ranked_screen <- function(file_path) {
+   df <- load_gene_file(file_path)
+   col_names_lower <- str_to_lower(names(df))
+   
+   # Look for rank columns with pos/neg indicators
+   pos_cols <- grep("pos|positive", col_names_lower)
+   neg_cols <- grep("neg|negative", col_names_lower)
+   
+   # Filter to only rank columns
+   rank_pattern <- "rank|score"
+   pos_cols <- pos_cols[grep(rank_pattern, col_names_lower[pos_cols])]
+   neg_cols <- neg_cols[grep(rank_pattern, col_names_lower[neg_cols])]
+   
+   result <- list(positive = NULL, negative = NULL)
+   
+   if (length(pos_cols) > 0) {
+     pos_df <- df[, c(1, pos_cols[1])]
+     names(pos_df) <- c("gene_id", "rank")
+     pos_df$rank <- as.numeric(pos_df$rank)
+     pos_df <- pos_df[!is.na(pos_df$rank), ]
+     pos_df <- pos_df[order(pos_df$rank), ]
+     result$positive <- pos_df
+   }
+   
+   if (length(neg_cols) > 0) {
+     neg_df <- df[, c(1, neg_cols[1])]
+     names(neg_df) <- c("gene_id", "rank")
+     neg_df$rank <- as.numeric(neg_df$rank)
+     neg_df <- neg_df[!is.na(neg_df$rank), ]
+     neg_df <- neg_df[order(neg_df$rank), ]
+     result$negative <- neg_df
+   }
+   
+   if (is.null(result$positive) && is.null(result$negative)) {
+     stop("Could not find positive or negative rank columns.")
+   }
+   
+   return(result)
+ }
> 
> # Load two-sheet Excel file
> load_two_sheet_excel <- function(file_path) {
+   sheets <- excel_sheets(file_path)
+   
+   if (length(sheets) < 2) {
+     stop("Excel file must contain at least 2 sheets.")
+   }
+   
+   # Read both sheets
+   sheet1 <- read_excel(file_path, sheet = 1)
+   sheet2 <- read_excel(file_path, sheet = 2)
+   
+   # Process sheet 1
+   gene_col_idx <- which(str_to_lower(names(sheet1)) == "gene_id")
+   if (length(gene_col_idx) == 0) gene_col_idx <- 1
+   s1_genes <- clean_gene_ids(sheet1[[gene_col_idx[1]]])
+   
+   # Process sheet 2
+   gene_col_idx <- which(str_to_lower(names(sheet2)) == "gene_id")
+   if (length(gene_col_idx) == 0) gene_col_idx <- 1
+   s2_genes <- clean_gene_ids(sheet2[[gene_col_idx[1]]])
+   
+   return(list(
+     s1 = data.frame(gene_id = s1_genes, stringsAsFactors = FALSE),
+     s2 = data.frame(gene_id = s2_genes, stringsAsFactors = FALSE),
+     names = sheets[1:2]
+   ))
+ }
> 
> # Load MAGeCK output
> load_mageck <- function(file_path) {
+   df <- load_gene_file(file_path)
+   col_names_lower <- str_to_lower(names(df))
+   
+   # Look for FDR column
+   fdr_cols <- grep("fdr", col_names_lower)
+   
+   if (length(fdr_cols) == 0) {
+     stop("MAGeCK file must contain an FDR column.")
+   }
+   
+   names(df)[fdr_cols[1]] <- "fdr_mageck"
+   df$fdr_mageck <- as.numeric(df$fdr_mageck)
+   df <- df[!is.na(df$fdr_mageck), ]
+   
+   return(df)
+ }
> 
> # ============================================================================
> # UI
> # ============================================================================
> 
> ui <- fluidPage(
+   titlePanel("CRISPR Screen Cross-Reference"),
+   
+   sidebarLayout(
+     sidebarPanel(
+       width = 3,
+       
+       fileInput("orig", "Original Gene List", accept = c(".csv", ".xlsx", ".xls")),
+       hr(),
+       
+       fileInput("s1", "Screen 1 (Ranked)", accept = c(".csv", ".xlsx", ".xls")),
+       conditionalPanel(
+         condition = "output.s1_loaded",
+         sliderInput("s1n", "Top N genes:", min = 50, max = 1000, value = 200, step = 50)
+       ),
+       hr(),
+       
+       fileInput("s2", "Screen 2 (Two-sheet Excel)", accept = c(".xlsx", ".xls")),
+       hr(),
+       
+       fileInput("s3", "Screen 3 (MAGeCK)", accept = c(".csv", ".xlsx", ".xls")),
+       conditionalPanel(
+         condition = "output.s3_loaded",
+         sliderInput("fdr", "FDR threshold:", min = 0.001, max = 1, value = 0.1, step = 0.01)
+       ),
+       hr(),
+       
+       fileInput("s4", "Screen 4 (Ranked)", accept = c(".csv", ".xlsx", ".xls")),
+       conditionalPanel(
+         condition = "output.s4_loaded",
+         sliderInput("s4n", "Top N genes:", min = 50, max = 1000, value = 200, step = 50)
+       ),
+       hr(),
+       
+       fileInput("s5", "Screen 5 (Ranked)", accept = c(".csv", ".xlsx", ".xls")),
+       conditionalPanel(
+         condition = "output.s5_loaded",
+         sliderInput("s5n", "Top N genes:", min = 50, max = 1000, value = 200, step = 50)
+       ),
+       hr(),
+       
+       fileInput("s6", "Screen 6 (Fixed List)", accept = c(".csv", ".xlsx", ".xls")),
+       hr(),
+       
+       checkboxInput("intersect", "Show only genes in Original list", FALSE)
+     ),
+     
+     mainPanel(
+       width = 9,
+       tabsetPanel(
+         tabPanel("Overview", 
+                  br(),
+                  DTOutput("overview"), 
+                  br(),
+                  plotOutput("bar", height = "400px")
+         ),
+         tabPanel("Combined Results", 
+                  br(),
+                  sliderInput("min_support", "Minimum screens supporting:", 
+                             min = 1, max = 10, value = 1, step = 1),
+                  DTOutput("combined"),
+                  br(),
+                  downloadButton("download_combined", "Download Combined Results")
+         ),
+         tabPanel("Visualizations", 
+                  br(),
+                  plotOutput("upset", height = "600px")
+         )
+       )
+     )
+   )
+ )
> 
> # ============================================================================
> # SERVER
> # ============================================================================
> 
> server <- function(input, output, session) {
+   
+   # Reactive values for loaded data
+   orig_data <- reactiveVal(NULL)
+   s1_data <- reactiveVal(NULL)
+   s2_data <- reactiveVal(NULL)
+   s3_data <- reactiveVal(NULL)
+   s4_data <- reactiveVal(NULL)
+   s5_data <- reactiveVal(NULL)
+   s6_data <- reactiveVal(NULL)
+   
+   # Load original list
+   observeEvent(input$orig, {
+     req(input$orig)
+     tryCatch({
+       df <- load_gene_file(input$orig$datapath)
+       orig_data(df$gene_id)
+       showNotification("Original list loaded!", type = "message")
+     }, error = function(e) {
+       showNotification(paste("Error:", e$message), type = "error")
+     })
+   })
+   
+   # Load Screen 1
+   observeEvent(input$s1, {
+     req(input$s1)
+     tryCatch({
+       data <- load_ranked_screen(input$s1$datapath)
+       s1_data(data)
+       showNotification("Screen 1 loaded!", type = "message")
+     }, error = function(e) {
+       showNotification(paste("Error:", e$message), type = "error")
+     })
+   })
+   
+   output$s1_loaded <- reactive({ !is.null(s1_data()) })
+   outputOptions(output, "s1_loaded", suspendWhenHidden = FALSE)
+   
+   # Load Screen 2
+   observeEvent(input$s2, {
+     req(input$s2)
+     tryCatch({
+       data <- load_two_sheet_excel(input$s2$datapath)
+       s2_data(data)
+       showNotification("Screen 2 loaded!", type = "message")
+     }, error = function(e) {
+       showNotification(paste("Error:", e$message), type = "error")
+     })
+   })
+   
+   # Load Screen 3
+   observeEvent(input$s3, {
+     req(input$s3)
+     tryCatch({
+       data <- load_mageck(input$s3$datapath)
+       s3_data(data)
+       showNotification("Screen 3 loaded!", type = "message")
+     }, error = function(e) {
+       showNotification(paste("Error:", e$message), type = "error")
+     })
+   })
+   
+   output$s3_loaded <- reactive({ !is.null(s3_data()) })
+   outputOptions(output, "s3_loaded", suspendWhenHidden = FALSE)
+   
+   # Load Screen 4
+   observeEvent(input$s4, {
+     req(input$s4)
+     tryCatch({
+       data <- load_ranked_screen(input$s4$datapath)
+       s4_data(data)
+       showNotification("Screen 4 loaded!", type = "message")
+     }, error = function(e) {
+       showNotification(paste("Error:", e$message), type = "error")
+     })
+   })
+   
+   output$s4_loaded <- reactive({ !is.null(s4_data()) })
+   outputOptions(output, "s4_loaded", suspendWhenHidden = FALSE)
+   
+   # Load Screen 5
+   observeEvent(input$s5, {
+     req(input$s5)
+     tryCatch({
+       data <- load_ranked_screen(input$s5$datapath)
+       s5_data(data)
+       showNotification("Screen 5 loaded!", type = "message")
+     }, error = function(e) {
+       showNotification(paste("Error:", e$message), type = "error")
+     })
+   })
+   
+   output$s5_loaded <- reactive({ !is.null(s5_data()) })
+   outputOptions(output, "s5_loaded", suspendWhenHidden = FALSE)
+   
+   # Load Screen 6
+   observeEvent(input$s6, {
+     req(input$s6)
+     tryCatch({
+       df <- load_gene_file(input$s6$datapath)
+       s6_data(df$gene_id)
+       showNotification("Screen 6 loaded!", type = "message")
+     }, error = function(e) {
+       showNotification(paste("Error:", e$message), type = "error")
+     })
+   })
+   
+   # Process filtered screens
+   screens <- reactive({
+     result <- list()
+     
+     # Screen 1
+     if (!is.null(s1_data())) {
+       s1 <- s1_data()
+       if (!is.null(s1$positive)) {
+         result$s1_pos <- head(s1$positive$gene_id, input$s1n)
+       }
+       if (!is.null(s1$negative)) {
+         result$s1_neg <- head(s1$negative$gene_id, input$s1n)
+       }
+     }
+     
+     # Screen 2
+     if (!is.null(s2_data())) {
+       s2 <- s2_data()
+       result$s2_sheet1 <- s2$s1$gene_id
+       result$s2_sheet2 <- s2$s2$gene_id
+     }
+     
+     # Screen 3
+     if (!is.null(s3_data())) {
+       s3 <- s3_data()
+       filtered <- s3[s3$fdr_mageck <= input$fdr, ]
+       result$s3_fdr <- filtered$gene_id
+     }
+     
+     # Screen 4
+     if (!is.null(s4_data())) {
+       s4 <- s4_data()
+       if (!is.null(s4$positive)) {
+         result$s4_pos <- head(s4$positive$gene_id, input$s4n)
+       }
+       if (!is.null(s4$negative)) {
+         result$s4_neg <- head(s4$negative$gene_id, input$s4n)
+       }
+     }
+     
+     # Screen 5
+     if (!is.null(s5_data())) {
+       s5 <- s5_data()
+       if (!is.null(s5$positive)) {
+         result$s5_pos <- head(s5$positive$gene_id, input$s5n)
+       }
+       if (!is.null(s5$negative)) {
+         result$s5_neg <- head(s5$negative$gene_id, input$s5n)
+       }
+     }
+     
+     # Screen 6
+     if (!is.null(s6_data())) {
+       result$s6_fixed <- s6_data()
+     }
+     
+     return(result)
+   })
+   
+   # Combined dataframe
+   combined_df <- reactive({
+     req(orig_data())
+     
+     orig <- orig_data()
+     scr <- screens()
+     
+     # Get all unique genes
+     all_genes <- unique(c(orig, unlist(scr)))
+     
+     # Build dataframe
+     df <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE)
+     df$in_original <- df$gene_id %in% orig
+     
+     # Add columns for each screen
+     for (screen_name in names(scr)) {
+       df[[screen_name]] <- df$gene_id %in% scr[[screen_name]]
+     }
+     
+     # Count supporting screens
+     screen_cols <- names(scr)
+     df$num_screens_supporting <- rowSums(df[, screen_cols, drop = FALSE])
+     
+     # Apply filters
+     df <- df[df$num_screens_supporting >= input$min_support, ]
+     
+     if (input$intersect) {
+       df <- df[df$in_original, ]
+     }
+     
+     # Sort by support
+     df <- df[order(-df$num_screens_supporting), ]
+     
+     return(df)
+   })
+   
+   # Overview table
+   output$overview <- renderDT({
+     req(orig_data())
+     
+     orig <- orig_data()
+     scr <- screens()
+     
+     screen_names <- names(scr)
+     gene_counts <- sapply(scr, length)
+     overlaps <- sapply(scr, function(x) length(intersect(x, orig)))
+     
+     df <- data.frame(
+       Screen = screen_names,
+       Total_Genes = gene_counts,
+       Overlap_with_Original = overlaps,
+       Overlap_Percentage = round(100 * overlaps / gene_counts, 1),
+       stringsAsFactors = FALSE
+     )
+     
+     datatable(df, options = list(pageLength = 20, dom = 't'), rownames = FALSE)
+   })
+   
+   # Combined results table
+   output$combined <- renderDT({
+     df <- combined_df()
+     datatable(df, 
+               options = list(pageLength = 25, scrollX = TRUE),
+               rownames = FALSE)
+   })
+   
+   # Bar plot
+   output$bar <- renderPlot({
+     req(orig_data())
+     
+     orig <- orig_data()
+     scr <- screens()
+     
+     if (length(scr) == 0) return(NULL)
+     
+     df <- data.frame(
+       Screen = names(scr),
+       Overlap = sapply(scr, function(x) length(intersect(x, orig))),
+       stringsAsFactors = FALSE
+     )
+     
+     ggplot(df, aes(x = reorder(Screen, Overlap), y = Overlap)) +
+       geom_col(fill = "steelblue") +
+       geom_text(aes(label = Overlap), hjust = -0.2, size = 4) +
+       coord_flip() +
+       labs(title = "Overlap Counts with Original List",
+            x = "Screen",
+            y = "Number of Overlapping Genes") +
+       theme_minimal() +
+       theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
+   })
+   
+   # UpSet plot
+   output$upset <- renderPlot({
+     req(orig_data())
+     
+     orig <- orig_data()
+     scr <- screens()
+     
+     if (length(scr) < 2) {
+       plot.new()
+       text(0.5, 0.5, "Need at least 2 screens loaded", cex = 1.5)
+       return()
+     }
+     
+     # Build upset dataframe
+     all_genes <- unique(c(orig, unlist(scr)))
+     upset_df <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE)
+     upset_df$Original <- upset_df$gene_id %in% orig
+     
+     for (screen_name in names(scr)) {
+       upset_df[[screen_name]] <- upset_df$gene_id %in% scr[[screen_name]]
+     }
+     
+     # Remove gene_id column
+     upset_df <- upset_df[, -1]
+     
+     # Create upset plot
+     upset(upset_df, 
+           colnames(upset_df),
+           name = "Gene Set Intersections",
+           width_ratio = 0.2,
+           set_sizes = FALSE)
+   })
+   
+   # Download handler
+   output$download_combined <- downloadHandler(
+     filename = function() {
+       paste0("combined_results_", Sys.Date(), ".csv")
+     },
+     content = function(file) {
+       df <- combined_df()
+       write_csv(df, file)
+     }
+   )
+ }
> 
> # ============================================================================
> # RUN APP
> # ============================================================================
> 
> shinyApp(ui = ui, server = server)

               

