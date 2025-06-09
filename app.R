# LDSeeker Shiny App with dynamic LD_info handling, download button, and cleanup
library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(DT)
library(bslib)

# UI
ui <- navbarPage(
  title = "LDSeeker",
  theme = bs_theme(bootswatch = "flatly", version = 5),
  
  tabPanel("LD Analysis",
           sidebarLayout(
             sidebarPanel(
               fileInput("file", "Upload GWAS File (Tab-separated .txt file)", accept = c('.txt')),
               numericInput("r2threshold", "R² Threshold", value = 0.8, min = 0, max = 1, step = 0.01),
               selectInput("pop", "Population", choices = c("EUR", "EAS", "AFR", "SAS"), selected = "EUR"),
               numericInput("maf", "MAF", value = 0.05, min = 0, max = 1, step = 0.01),
               selectInput("ref", "Reference Panel", choices = c("Pheno_Scanner", "TOP_LD", "Hap_Map", "all_panels")),
               fileInput("imp_list", "Optional SNP List (.txt)", accept = ".txt"),
               actionButton("run", "Run LDSeeker", class = "btn-primary")
             ),
             
             mainPanel(
               h4("LD Analysis Results"),
               withSpinner(DTOutput("result")),
               downloadButton("downloadData", "Download Results"),
               textOutput("runtime")
             )
           )
  ),
  
  tabPanel("Help",
           fluidPage(
             h2("Help"),
             p("Upload your GWAS file (CSV/TSV). Optionally upload a SNP list."),
             p("Set R² threshold, MAF, population, and reference panel."),
             p("Click 'Run LDSeeker' to process the data and view results.")
           )
  ),
  
  tabPanel("About",
           fluidPage(
             h2("About"),
             p("LDSeeker Version 1.0.0 (August 2025)"),
             p("Author: Pantelis Bagos"),
             p("Open-source under GPLv3 License")
           )
  )
)

# Server logic
temp_dir <- tempdir()
server <- function(input, output, session) {
  # Keep track of the current result file
  output_file <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    req(input$file)
    start_time <- Sys.time()
    
    # Copy uploaded GWAS file to temp directory
    gwas_path <- file.path(temp_dir, basename(input$file$name))
    file.copy(input$file$datapath, gwas_path, overwrite = TRUE)
    
    # Handle optional SNP list
    imp_args <- c()
    if (!is.null(input$imp_list)) {
      imp_path <- file.path(temp_dir, basename(input$imp_list$name))
      file.copy(input$imp_list$datapath, imp_path, overwrite = TRUE)
      imp_args <- c("--imp_list", shQuote(imp_path))
    }
    
    # Build and run Python CLI
    args <- c(
      "LDSeeker.py",
      "--file-path", shQuote(gwas_path),
      "--r2threshold", input$r2threshold,
      "--pop", shQuote(input$pop),
      "--maf", input$maf,
      "--ref", shQuote(input$ref),
      imp_args
    )
    try(system2("python", args, stdout = TRUE, stderr = TRUE), silent = TRUE)
    
    # Find output file(s) starting with LD_info
    ld_files <- list.files(
      './',
      pattern    = "^LD_info_[A-Za-z0-9_]+\\.txt$",
      full.names = TRUE
    )
    
    if (length(ld_files) > 0) {
      # Use the first matching file
      ld_path <- ld_files[[1]]
      df <- read.table(ld_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      
      # Render the table
      output$result <- renderDT({
        datatable(df, options = list(pageLength = 10, scrollX = TRUE))
      })
      
      # Show runtime and filename
      output$runtime <- renderText({
        sprintf("File: %s | Total Time: %.2f seconds", 
                basename(ld_path), 
                as.numeric(difftime(Sys.time(), start_time, units = "secs")))
      })
      
      # Store for download and cleanup
      output_file(ld_path)
    } else {
      showNotification("No LD_info file found. Python script may have failed.", type = "error")
      output$result <- renderDT(NULL)
      output$runtime <- renderText("")
      output_file(NULL)
    }
  })
  
  # Download handler for results
  output$downloadData <- downloadHandler(
    filename = function() {
      req(output_file())
      paste0(tools::file_path_sans_ext(basename(output_file())), ".csv")
    },
    content = function(file) {
      df <- read.table(output_file(), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      write.csv(df, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  # Cleanup LD_info files when session ends
  session$onSessionEnded(function() {
    to_remove <- list.files('./',   pattern    = "^LD_info_[A-Za-z0-9_]+\\.txt$",  full.names = TRUE)
    file.remove(to_remove)
  })
}

# Run the app
shinyApp(ui = ui, server = server)


 