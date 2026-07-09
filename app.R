library(shiny)
library(bslib)
library(fontawesome)
library(data.table)
library(DT)
library(plotly)
library(base64enc)
library(gaston)
library(visNetwork)
library(igraph)
library(arrow)

# Fix: Ensure 'p' refers to the Shiny HTML paragraph tag
p <- shiny::p

# 500MB max Input size
options(shiny.maxRequestSize = 500 * 1024^2)

# --- GENERATE DEMO FILES ---
if(!file.exists("LD_patterns_demo.txt")){
  writeLines(c(
    "SNP\tCHR\tBP\tA1\tA2\tBETA\tSE",
    "rs743749\t22\t37398195\tA\tG\t-0.006387442\t9.898344223",
    "rs9306493\t22\t45682425\tA\tG\t-0.015022874\t9.594216875",
    "rs739043\t22\t37645230\tG\tA\t-0.005243055\t9.788226204",
    "rs242885\t22\t34423169\tA\tG\t-0.019996628\t9.449498344",
    "rs5765043\t22\t45231883\tG\tA\t-0.007225636\t9.599864029",
    "rs9625200\t22\t27700318\tA\tG\t0.007320953\t9.914661823",
    "rs17807317\t22\t17680519\tC\tA\t0.005180513\t9.805693943",
    "rs5764726\t22\t19528458\tC\tA\t0.161667412\t0.544457863",
    "rs80571\t22\t49504833\tA\tG\t0.15896526\t0.54073433"
  ), "LD_patterns_demo.txt")
}

if(!file.exists("LD_assoc_demo.txt")){
  writeLines(c(
    "SNP\tCHR\tBP\tAlleles\tMAF\tDistance\tCorrelated_Alleles\tP-value",
    "rs7837688\t8\t128539360\t(T/G)\t0.0909\t0\t\"T=T,G=G\"\t1.02E-06",
    "rs4242384\t8\t128518554\t(C/A)\t0.0808\t-20806\t\"C=T,A=G\"\t2.66E-06",
    "rs4242382\t8\t128517573\t(A/G)\t0.0808\t-21787\t\"A=T,G=G\"\t2.82E-06",
    "rs11988857\t8\t128531873\t(G/A)\t0.1061\t-7487\t\"G=T,A=G\"\t3.78E-06",
    "rs6983267\t8\t128413305\t(G/T)\t0.4848\t-126055\tNA\t1.05E-05",
    "rs7017300\t8\t128525268\t(C/A)\t0.1263\t-14092\t\"C=T,A=G\"\t1.80E-05"
  ), "LD_assoc_demo.txt")
}

if(!file.exists("LDresults_demo.txt")){
  writeLines(c(
    "CHROM_A\tPOS1\trsID1\tMAJ1\tNONMAJ1\tMAF1\tCHROM_B\tPOS2\trsID2\tMAJ2\tNONMAJ2\tMAF2\tR\tR2",
    "22\t16576248\trs5746647\tT\tG\t0.06\t22\t16577670\trs4239851\tG\tA\t0.06\t1\t1",
    "22\t16576248\trs5746647\tT\tG\t0.06\t22\t16578670\trs1987449\tT\tC\t0.06\t1\t1",
    "22\t16576248\trs5746647\tT\tG\t0.06\t22\t16582469\trs2006108\tG\tA\t0.06\t1\t1",
    "22\t16576248\trs5746647\tT\tG\t0.06\t22\t16583531\trs4819768\tA\tC\t0.06\t1\t1",
    "22\t16576248\trs5746647\tT\tG\t0.06\t22\t16584133\trs5747963\tC\tT\t0.06\t1\t1",
    "22\t16576248\trs5746647\tT\tG\t0.06\t22\t16584520\trs2006338\tT\tC\t0.06\t1\t1"
  ), "LDresults_demo.txt")
}

# --- CONFIGURATION & DATA ---
default_pops <- c("EUR", "AFR", "AMR", "EAS", "SAS")
hapmap_pops <- c("YRI", "CHB", "JPT", "CEU", "MKK", "LWK", "CHD", "GIH", "TSI", "MXL", "ASW")
HGDP_pops <- c("OCEANIA", "EUROPE", "AFRICA", "EAST_ASIA", "CENTRAL_SOUTH_ASIA", "MIDDLE_EAST", "AMERICA")
UKBB_pops <- c("EUR", "AFR", "EAS", "CSA", "MID", "AMR")
LASI_DAD_pops <- c("IND")

# --- ADVANCED BIOINFORMATIC THEME ---
bio_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly", 
  primary = "#00d1b2",   # Bio-Teal
  secondary = "#2c3e50", # Deep Navy
  success = "#23d160",
  warning = "#ffb703",
  base_font = font_google("Inter"),
  heading_font = font_google("Rajdhani"), 
  code_font = font_google("Fira Code")
) %>%
  bs_add_rules("
    /* Fluid Layout & Background */
    body { background-color: #fcfdfe; }

    /* Navbar Customization */
    .navbar {
      background: linear-gradient(90deg, #1a252f 0%, #2c3e50 100%) !important;
      border-bottom: 4px solid #00d1b2;
      padding: 0.8rem 2rem;
      box-shadow: 0 4px 12px rgba(0,0,0,0.1);
    }
    .navbar-brand {
      font-weight: 800;
      letter-spacing: 1px;
      text-transform: uppercase;
      font-family: 'Rajdhani', sans-serif;
    }
     
    /* Nav Icons Alignment */
    .nav-link svg { margin-right: 8px; }

    /* Home Page Cards Styling - Matplotlib Colors */
    .feature-card {
      border: none !important;
      transition: all 0.4s cubic-bezier(0.175, 0.885, 0.32, 1.275);
      overflow: hidden;
      color: white !important;
      height: 100%;
    }
    .feature-card:hover {
      transform: translateY(-10px);
      box-shadow: 0 15px 30px rgba(0,0,0,0.15) !important;
    }
    .feature-card .card-header {
      background: transparent !important;
      border: none;
      font-size: 1.25rem;
      padding-top: 1.5rem;
      font-weight: 700;
      color: white !important;
    }
    .feature-card p { opacity: 0.95; font-size: 0.95rem; }
    .feature-card .btn {
      background: rgba(255,255,255,0.25);
      border: 1px solid rgba(255,255,255,0.5);
      color: white;
      font-weight: 600;
      backdrop-filter: blur(5px);
    }
    .feature-card .btn:hover { background: white; color: #333; }

    /* Input Cards Styling */
    .input-group-card {
      background-color: #ffffff;
      border: 1px solid #e2e8f0;
      border-radius: 8px;
      box-shadow: 0 1px 3px rgba(0,0,0,0.05);
      height: 100%;
    }
    .input-group-header {
      background-color: #f8fafc;
      border-bottom: 1px solid #e2e8f0;
      padding: 10px 15px;
      font-weight: 700;
      color: #2c3e50;
      font-size: 0.95rem;
      border-radius: 8px 8px 0 0;
    }
    .input-group-body { padding: 15px; }

    .about-logo {
      max-width: 200px;
      width: 100%;
      height: auto;
      margin-bottom: 10px;
    }
  ")

# --- UI ---
ui <- page_navbar(
  id = "navbar_id",
  title = div(
    style = "display: flex; align-items: center;",
    span(fa("magnifying-glass", fill = "#00d1b2"), style = "margin-right: 12px; font-size: 1.5rem;"),
    span("LDSeeker", style = "color: white; font-size: 1.4rem;")
  ),
  theme = bio_theme,
  
  nav_spacer(),
  
  # ---------- Home Tab ----------
  nav_panel(
    title = span(fa("house"), " Home"),
    value = "home_tab",
    fluidPage(
      div(
        style = "text-align: center; padding: 60px 0 40px 0;",
        tags$img(
          src = if(file.exists("logo.png")) base64enc::dataURI(file = "logo.png", mime = "image/png") else "",
          alt = "Logo",
          class = "about-logo",
          style = if(!file.exists("logo.png")) "display:none;" else ""
        ),
        p("Exploring linkage disequilibrium in GWAS using multiple panels", style = "color: #718096; font-size: 1.25rem; font-weight: 300; letter-spacing: 1px;"),
        hr(style = "width: 100px; margin: 20px auto; border-top: 3px solid #00d1b2;")
      ),
      layout_column_wrap(
        width = 1/2, 
        style = "max-width: 1200px; margin: 0 auto; gap: 2rem;",
        # Matplotlib Blue: #1f77b4
        card(
          class = "feature-card",
          style = "background-color: #1f77b4 !important;",
          card_header(span(fa("project-diagram"), " LDpatterns")),
          card_body(
            p("Explore Linkage Disequilibrium patterns given the desired settings in a large collection of LD reference panels and multiple populations"),
            div(style="flex-grow: 1;"),
            actionButton("btn_goto_analysis", "Go to LDpatterns", class = "w-100 mt-3")
          )
        ),
        # Matplotlib Orange: #ff7f0e
        card(
          class = "feature-card",
          style = "background-color: #ff7f0e !important;",
          card_header(span(fa("chart-line"), " LDassoc")),
          card_body(
            p("Interactively visualize association p-value results and linkage disequilibrium patterns for a genomic region of interest."),
            div(style="flex-grow: 1;"),
            actionButton("btn_goto_ldassoc", "Go to LDassoc", class = "w-100 mt-3")
          )
        ),
        # Matplotlib Green: #2ca02c
        card(
          class = "feature-card",
          style = "background-color: #2ca02c !important;",
          card_header(span(fa("scissors"), " LDPruning")),
          card_body(
            p("Identify independent genetic variants and prune redundant datasets based on configurable R² and p-value thresholds."),
            div(style="flex-grow: 1;"),
            actionButton("btn_goto_pruning", "Go to LDpruning", class = "w-100 mt-3")
          )
        ),
        # Matplotlib Red: #d62728
        card(
          class = "feature-card",
          style = "background-color: #d62728 !important;",
          card_header(span(fa("th"), " LD Visualization")),
          card_body(
            p("Visualize LD patterns with interactive heatmaps and LD plots"),
            div(style="flex-grow: 1;"),
            actionButton("btn_goto_heatmap", "Visualize Data", class = "w-100 mt-3")
          )
        )
      )
    )
  ),
  
  # ---------- LDpatterns Tab ----------
  nav_panel(
    title = span(fa("project-diagram"), " LDpatterns"),
    value = "analysis_tab",
    fluidPage(
      # --- SETTINGS ROW ---
      card(
        card_header(span(fa("cogs"), " Analysis Settings")),
        card_body(
          layout_column_wrap(
            width = 1/3,
            # Col 1: Input
            div(class = "input-group-card",
                div(class = "input-group-header", "1. Input Data"),
                div(class = "input-group-body",
                    fileInput("input_file", "Insert GWAS file or rsIDs", accept = c(".txt", ".tsv", ".csv")),
                    textAreaInput("paste_rsids", "Or Paste rsIDs:", placeholder = "rs123\nrs456\n...", rows = 3),
                    actionLink("load_example", "Load Example Data", icon = icon("arrow-circle-down"))
                )
            ),
            # Col 2: Reference
            div(class = "input-group-card",
                div(class = "input-group-header", "2. Reference & Mode"),
                div(class = "input-group-body",
                    selectInput("ref_panel", "Reference Panel", 
                                choices = c(
                                  "1000 Genomes Project (hg38) - High Coverage" = "1000G_hg38_high_cov",
                                  "1000 Genomes Project (hg38)" = "1000G_hg38", 
                                  "Pheno Scanner" = "Pheno_Scanner", 
                                  "TOP-LD" = "TOP_LD", 
                                  "HapMap" = "Hap_Map",
                                  "LASI-DAD" = "LASI_DAD",
                                  "Human Genome Diversity Project (HGDP)" = "HGDP",
                                  "UK Biobank (UKBB)" = "UKBB"
                                ), 
                                selected = "1000G_hg38"),
                    selectInput("population", "Population Ancestry", choices = default_pops),
                    radioButtons("pairwise", "LD Mode", choices = c("Pairwise" = "YES", "Standard" = "NO"), inline = TRUE)
                )
            ),
            # Col 3: Filters
            div(class = "input-group-card",
                div(class = "input-group-header", "3. Filters & Run"),
                div(class = "input-group-body",
                    layout_column_wrap(width=1/2,
                                       numericInput("r2threshold", "R² Threshold", value = 0.8, min = 0, max = 1, step = 0.01),
                                       numericInput("maf", "MAF Threshold", value = 0.01, min = 0, max = 1, step = 0.01)
                    ),
                    actionButton("run_ld", "EXECUTE LDpatterns", icon = icon("play"), class = "btn-primary w-100 mt-3")
                )
            )
          )
        )
      ),
      br(),
      uiOutput("preview_section"),
      br(),
      layout_column_wrap(
        width = 1/3,
        value_box(title = "Unique Variant Count", value = uiOutput("val_variants"), showcase = fa("microscope"), theme = "primary"),
        value_box(title = "Average R²", value = uiOutput("val_mean_r2"), showcase = fa("chart-bar"), theme = "info"),
        value_box(title = "High LD Links", value = uiOutput("val_high_ld"), showcase = fa("link"), theme = "success")
      ),
      br(),
      card(
        card_header(
          div(
            class = "d-flex justify-content-between align-items-center",
            span("LD Result Matrix"),
            downloadButton("download_ld_matrix", "Download NxN Matrix (R²)", class = "btn-sm btn-light text-dark")
          )
        ),
        full_screen = TRUE,
        DTOutput("ld_table")
      )
    )
  ),
  
  # ---------- LDpruning Tab (NEW) ----------
  nav_panel(
    title = span(fa("scissors"), " LDpruning"),
    value = "pruning_tab",
    fluidPage(
      card(
        card_header(span(fa("scissors"), " Pruning Configuration")),
        card_body(
          layout_column_wrap(
            width = 1/3,
            # Col 1: Input Data
            div(class = "input-group-card",
                div(class = "input-group-header", "1. Input Data"),
                div(class = "input-group-body",
                    fileInput("prune_file", "Insert GWAS file (Tab Delimited)", accept = c(".txt", ".tsv", ".csv")),
                    uiOutput("prune_cols_ui") # Dynamic Column Selectors
                )
            ),
            # Col 2: Pruning & LD Parameters
            div(class = "input-group-card",
                div(class = "input-group-header", "2. Pruning & LD Parameters"),
                div(class = "input-group-body",
                    layout_column_wrap(width = 1/2,
                                       numericInput("prune_r2_thresh", "R² Pruning Cutoff", value = 0.2, min = 0, max = 1, step = 0.01),
                                       numericInput("prune_maf", "MAF Threshold", value = 0.01, min = 0, max = 1, step = 0.01)
                    ),
                    layout_column_wrap(width = 1/2,
                                       selectInput("prune_metric", "Pruning Priority", choices = c("Smallest P-value" = "P", "Largest |Z-score|" = "Z"), selected = "P"),
                                       numericInput("prune_z_thresh", "Pruning |Z| Threshold", value = 2.0, min = 0, step = 0.1)
                    ),
                    layout_column_wrap(width = 1/2,
                                       selectInput("prune_pval_dir", "Pruning P Filter", choices = c("P < Threshold" = "below", "P > Threshold" = "above", "No P Filter" = "none"), selected = "none"),
                                       numericInput("prune_pval_thresh", "Pruning P Threshold", value = 0.05, min = 0, max = 1, step = 1e-4)
                    ),
                    hr(),
                    h6("Significance pre-filter", style = "margin-top: 5px; color: #2c3e50; font-weight: 700;"),
                    layout_column_wrap(width = 1/2,
                                       selectInput("prune_significance", "Keep", choices = c("No significance filter" = "none", "Significant only" = "significant", "Non-significant only" = "nonsignificant"), selected = "none"),
                                       selectInput("prune_sig_metric", "Metric", choices = c("P-value" = "P", "|Z-score|" = "Z"), selected = "P")
                    ),
                    numericInput("prune_sig_thresh", "Significance Threshold", value = 5e-8, min = 0, step = 1e-8)
                )
            ),
            # Col 3: Reference & Run
            div(class = "input-group-card",
                div(class = "input-group-header", "3. Reference & Execute"),
                div(class = "input-group-body",
                    selectInput("prune_ref_panel", "Reference Panel", 
                                choices = c(
                                  "1000 Genomes (hg38)"="1000G_hg38", 
                                  "1000G High Cov"="1000G_hg38_high_cov",
                                  "Pheno Scanner"="Pheno_Scanner", 
                                  "TOP-LD"="TOP_LD", 
                                  "HapMap"="Hap_Map",
                                  "LASI-DAD" = "LASI_DAD",
                                  "HGDP" = "HGDP",
                                  "UK Biobank" = "UKBB"
                                ), selected = "UKBB"),
                    selectInput("prune_population", "Population", choices = default_pops),
                    actionButton("run_pruning", "EXECUTE LDpruning", icon = icon("cut"), class = "btn-success w-100 mt-3")
                )
            )
          )
        )
      ),
      br(),
      # NEW: Input Preview for Pruning
      uiOutput("prune_preview_section"),
      br(),
      
      # Results Area Stats
      layout_column_wrap(
        width = 1/3,
        value_box(title = "Initial Variants", value = uiOutput("val_prune_initial"), showcase = fa("list-ol"), theme = "secondary"),
        value_box(title = "Variants Removed", value = uiOutput("val_prune_removed"), showcase = fa("trash-alt"), theme = "danger"),
        value_box(title = "Final Variants", value = uiOutput("val_prune_final"), showcase = fa("check-circle"), theme = "success")
      ),
      br(),
      
      # Pruned Dataset
      card(
        card_header(
          div(
            class = "d-flex justify-content-between align-items-center",
            span(fa("table"), " Final Pruned Dataset"),
            downloadButton("download_pruned", "Download Pruned Data", class = "btn-sm btn-light text-dark")
          )
        ),
        full_screen = TRUE,
        DTOutput("pruned_table")
      )
    )
  ),
  
  # ---------- LDassoc Tab ----------
  nav_panel(
    title = span(fa("chart-line"), " LDassoc"),
    value = "ldassoc_tab",
    fluidPage(
      # --- CONFIGURATION ROW ---
      card(
        card_header(span(fa("sliders-h"), " Association Analysis Configuration")),
        card_body(
          layout_column_wrap(
            width = 1/3,
            # Col 1: Data Inputs
            div(class = "input-group-card",
                div(class = "input-group-header", "1. Data Source"),
                div(class = "input-group-body",
                    fileInput("assoc_file", "Upload Association Data", accept = c(".txt", ".tsv", ".csv")),
                    actionLink("load_example_assoc", "Load Example Data", icon = icon("arrow-circle-down")),
                    hr(),
                    textInput("assoc_index_snp", "Index rsID", placeholder = "e.g., rs7837688"),
                    textInput("assoc_target_chr", "Chromosome", placeholder = "e.g., 8")
                )
            ),
            # Col 2: Mapping (Dynamic)
            div(class = "input-group-card",
                div(class = "input-group-header", "2. Map Columns"),
                div(class = "input-group-body",
                    uiOutput("assoc_columns_ui")
                )
            ),
            # Col 3: Parameters
            div(class = "input-group-card",
                div(class = "input-group-header", "3. Parameters & Run"),
                div(class = "input-group-body",
                    layout_column_wrap(width = 1/2,
                                       selectInput("assoc_ref_panel", "Reference Panel", 
                                                   choices = c(
                                                     "1000 Genomes Project (hg38) - High Coverage" = "1000G_hg38_high_cov",
                                                     "1000 Genomes (hg38)"="1000G_hg38", 
                                                     "Pheno Scanner"="Pheno_Scanner", 
                                                     "TOP-LD"="TOP_LD", 
                                                     "HapMap"="Hap_Map",
                                                     "LASI-DAD" = "LASI_DAD",
                                                     "Human Genome Diversity Project (HGDP)" = "HGDP",
                                                     "UK Biobank (UKBB)" = "UKBB"
                                                   )),
                                       selectInput("assoc_population", "Population", choices = default_pops)
                    ),
                    layout_column_wrap(width = 1/2,
                                       numericInput("assoc_r2", "R² Cutoff", value = 0.1, step=0.01),
                                       numericInput("assoc_maf", "MAF Cutoff", value = 0.01, step=0.01)
                    ),
                    selectInput("ldassoc_color_scale", "Plot Palette", choices = c("Viridis", "Cividis", "Plasma", "Inferno", "Magma", "Blues", "RdBu", "YlGnBu", "Greys"), selected = "Viridis"),
                    actionButton("run_ldassoc", "EXECUTE LDassoc", icon = icon("chart-line"), class = "btn-primary w-100 mt-2")
                )
            )
          )
        )
      ),
      br(),
      
      
      # --- RESULTS ROW ---
      
      
      
      
      layout_column_wrap(
        width = 1/2,
        card(
          card_header(span(fa("table"), " Association Data Preview")),
          DTOutput("assoc_preview_table")
        ),
        card(
          card_header(span(fa("dna"), " LDassoc Results")),
          full_screen = TRUE,
          DTOutput("ldassoc_result_table")
        )
      ),
      
      br(),
      card(
        card_header(span(fa("chart-area"), " Regional Association & LD Plot")),
        plotlyOutput("ldassoc_plot", height = "500px")
      )
      
    )
  ),
  
  
  # ---------- Visualization Tab (Updated) ----------
  nav_panel(
    title = span(fa("eye"), " LD Visualization"),
    value = "viz_tab",
    layout_sidebar(
      sidebar = sidebar(
        title = "Settings",
        fileInput("ld_heat_file", "Upload Results", accept = c(".txt", ".tsv")),
        actionLink("load_example_heatmap", "Load Example Data", icon = icon("arrow-circle-down")),
        hr(),
        uiOutput("heatmap_chr_ui"),
        numericInput("n_rows", "Rows per Batch", value = 100, min = 10, step = 10),
        uiOutput("heatmap_slider_ui"),
        hr(),
        uiOutput("ld_metric_ui"),
        sliderInput("r2_range", "Intensity Range", min = 0, max = 1, value = c(0, 1), step = 0.01),
        selectInput("color_scale", "Palette", choices = c("Viridis", "Cividis", "Plasma", "Inferno", "Magma", "Blues", "RdBu", "YlGnBu", "Greys"), selected = "RdBu"),
        
        hr(),
        h5("Network Customization"),
        selectInput("net_layout", "Network Layout", 
                    choices = c(
                      "Fruchterman-Reingold" = "layout_with_fr",
                      "Kamada-Kawai" = "layout_with_kk",
                      "Nicely (Auto)" = "layout_nicely",
                      "Circle" = "layout_in_circle",
                      "Grid" = "layout_on_grid",
                      "Sphere" = "layout_on_sphere",
                      "Random" = "layout_randomly",
                      "Davidson-Harel" = "layout_with_dh",
                      "GEM" = "layout_with_gem",
                      "Graphopt" = "layout_with_graphopt",
                      "MDS" = "layout_with_mds"
                    ), 
                    selected = "layout_with_fr"),
        selectInput("net_node_shape", "Node Shape", choices = c("Dot"="dot", "Square"="square", "Triangle"="triangle", "Diamond"="diamond", "Star"="star"), selected = "dot"),
        selectInput("net_node_color", "Node Color", 
                    choices = c(
                      "Default"="#97C2FC",
                      "Bio-Teal"="#00d1b2", 
                      "Deep Navy"="#2c3e50", 
                      "Red"="#e74c3c", 
                      "Orange"="#e67e22", 
                      "Green"="#27ae60",
                      "Purple"="#8e44ad",
                      "Grey"="#95a5a6", 
                      "Black"="#000000"
                    ), 
                    selected = "#97C2FC"),
        
        checkboxInput("show_values", "Show Labels (Heatmap)", value = FALSE)
      ),
      navset_card_tab(
        full_screen = TRUE,
        nav_panel("Interactive Heatmap", plotlyOutput("heatmap", height = "750px")),
        nav_panel("Interactive Network", visNetworkOutput("ld_network", height = "750px")),
        nav_panel("LD Plot", plotOutput("gaston_plot", height = "750px"))
      )
    )
  ),
  
  
  # ---------- SNP to Gene Tab ----------
  nav_panel(
    title = "SNP to Gene",
    value = "snp_gene_tab",
    fluidPage(
      card(
        card_header(span(fa("dna"), " SNP to Gene Mapping")),
        card_body(
          layout_column_wrap(
            width = 1/3,
            div(class = "input-group-card",
                div(class = "input-group-header", "1. Input Data"),
                div(class = "input-group-body",
                    fileInput("snp_gene_file", "Upload GWAS, rsID list, or LDSeeker results", accept = c(".txt", ".tsv", ".csv")),
                    textAreaInput(
                      "snp_gene_paste_rsids",
                      "Or paste rsIDs / SNPs:",
                      placeholder = "rs123\nrs456\nrs789",
                      rows = 5
                    )
                )
            ),
            div(class = "input-group-card",
                div(class = "input-group-header", "2. SNP Identifier Columns"),
                div(class = "input-group-body",
                    uiOutput("snp_gene_cols_ui"),
                    p(
                      em("Detected examples include rsID1, rsID_1, rsID2, rsid1, rsid2, rsid, rsID_A, rsID_B and the same patterns using SNP/snp."),
                      style = "font-size: 0.9rem; color: #718096;"
                    ),
                    uiOutput("snp_gene_chr_ui")
                )
            ),
            div(class = "input-group-card",
                div(class = "input-group-header", "3. Reference & Run"),
                div(class = "input-group-body",
                    uiOutput("snp_gene_ref_status"),
                    actionButton("run_snp_gene", "MAP SNPs TO GENES", icon = icon("dna"), class = "btn-primary w-100 mt-3")
                )
            )
          )
        )
      ),
      br(),
      uiOutput("snp_gene_preview_section"),
      br(),
      layout_column_wrap(
        width = 1/3,
        value_box(title = "Unique Input SNPs", value = uiOutput("val_snp_gene_input"), showcase = fa("list"), theme = "primary"),
        value_box(title = "Mapped SNPs", value = uiOutput("val_snp_gene_mapped"), showcase = fa("link"), theme = "success"),
        value_box(title = "Unique Genes", value = uiOutput("val_snp_gene_genes"), showcase = fa("dna"), theme = "info")
      ),
      br(),
      card(
        card_header(
          div(
            class = "d-flex justify-content-between align-items-center",
            span(fa("table"), " SNP to Gene Results"),
            downloadButton("download_snp_gene", "Download Mapped Results", class = "btn-sm btn-light text-dark")
          )
        ),
        full_screen = TRUE,
        DTOutput("snp_gene_result_table")
      )
    )
  ),
  
  # ---------- Documentation ----------
  nav_panel(
    title = span(fa("book"), " Documentation"),
    value = "doc_tab",
    fluidPage(
      div(class = "container", style = "max-width: 1200px;",
          h2("Documentation", style = "color: #2c3e50; font-family: 'Rajdhani'; font-weight: 800; border-bottom: 2px solid #00d1b2; padding-bottom: 10px; margin-bottom: 30px;"),
          
          accordion(
            id = "doc_accordion",
            open = "Software Overview",
            
            accordion_panel(
              title = "1. Software Overview & Methodology",
              icon = fa("circle-info"),
              div(
                p(strong("LDSeeker"), " is an open-source tool designed to query and compute Linkage Disequilibrium (LD) metrics across diverse human populations using multiple reference panels."),
                tags$ul(
                  tags$li(strong("LDpatterns:"), " Explore LD given a set of rsIDs (from GWAS or a list). In 'Standard' mode, it queries LD relative to the input set. In 'Pairwise' mode, it computes an NxN correlation matrix for all input variants."),
                  tags$li(strong("LDassoc:"), " It highlights an 'Index SNP' and maps all other regional SNPs relative to its LD and genomic position."),
                  tags$li(strong("LDpruning:"), " Pruning is used to identify independent genetic signals. Useful for obtaining independent sets from a GWAS."),
                  tags$li(strong("SNP to Gene:"), " Annotate uploaded GWAS files, rsID/SNP lists, or LDSeeker outputs by matching SNP identifiers to gene symbols from local chromosome-wise parquet mapping files.")
                )
              )
            ),
            
            accordion_panel(
              title = "2. Input File Specifications",
              icon = fa("file-import"),
              div(
                p("All tools support Tab-Delimited (.txt, .tsv). Maximum file size is 500MB."),
                h6("Required Columns per Module:", style = "color: #2c3e50; margin-top: 15px;"),
                tags$table(class = "table table-sm table-hover",
                           tags$thead(tags$tr(tags$th("Module"), tags$th("Required Fields"), tags$th("Optional Fields"))),
                           tags$tbody(
                             tags$tr(tags$td("LDpatterns"), tags$td("SNP (rsID), CHR"), tags$td("BP (Pos)")),
                             tags$tr(tags$td("LDassoc"), tags$td("SNP, CHR, BP, P-value"), tags$td("Alleles, BETA")),
                             tags$tr(tags$td("LDpruning"), tags$td("SNP, CHR"), tags$td("P-value (for significance-based pruning)")),
                             tags$tr(tags$td("SNP to Gene"), tags$td("At least one rsID/SNP identifier column"), tags$td("CHR/CHROM columns help restrict the chromosome-wise parquet lookup"))
                           )
                ),
                p(em("Note: Genetic coordinates should be consistent with the selected Reference Panel build (e.g., hg38 for 1000G High Cov)."), style = "font-size: 0.9rem;")
              )
            ),
            
            accordion_panel(
              title = "3. LD Metrics & Thresholds",
              icon = fa("chart-simple"),
              layout_column_wrap(
                width = 1/3,
                div(h6("Linkage Disequilibrium (r²)"), p("Measures the squared correlation of allele frequencies. Values range from 0 to 1, where 1 indicates perfect correlation.")),
                div(h6("D-Prime (D')"), p("Normalized linkage disequilibrium. Values near 1 indicate no historical recombination between variants, even if r² is low.")),
                div(h6("Minor Allele Frequency (MAF)"), p("Filters out variants where the minor allele occurs below a specific frequency in the population (e.g., < 0.01)."))
              )
            ),
            
            accordion_panel(
              title = "4. SNP to Gene Mapping",
              icon = fa("dna"),
              div(
                p(strong("SNP to Gene"), " maps SNP/rsID identifiers from uploaded GWAS files, plain rsID lists, or LDSeeker result tables to genes."),
                tags$ul(
                  tags$li("The local reference folder must be named ", code("snps_genes_ref"), " and must be located in the app working directory."),
                  tags$li("The reference files are chromosome-wise parquet files such as ", code("snp_gene_map_chr1.parquet"), ", ", code("snp_gene_map_chr2.parquet"), ", ..., read with the ", code("arrow"), " library."),
                  tags$li("Each parquet file must contain exactly the mapping columns ", code("Chromosome"), ", ", code("SNP"), " and ", code("Gene"), "."),
                  tags$li("The app automatically detects SNP identifier columns such as ", code("rsID1"), ", ", code("rsID_1"), ", ", code("rsID2"), ", ", code("rsid1"), ", ", code("rsid2"), ", ", code("rsid"), ", ", code("rsID_A"), ", ", code("rsID_B"), " and equivalent ", code("SNP"), "/", code("snp"), " names."),
                  tags$li("For each selected identifier column, the final table receives a corresponding gene annotation column, for example ", code("rsID1_Gene"), " and ", code("rsID2_Gene"), ". If one SNP maps to multiple genes, the genes are collapsed with semicolons.")
                ),
                p(em("Tip: If your input contains chromosome columns such as CHR, CHROM, CHROM_A or CHROM_B, the app uses them to read only the required chromosome parquet files."), style = "font-size: 0.9rem;")
              )
            ),
            
            accordion_panel(
              title = "5. Reference Panel Statistics",
              icon = fa("database"),
              div(
                p("Detailed population counts and available samples for each reference panel."),
                tableOutput("ref_pop_table"),
                hr(),
                div(
                  style = "display: grid; grid-template-columns: 1fr; gap: 1rem;",
                  card(
                    full_screen = TRUE,
                    card_header(span(fa("users"), " Total samples per reference panel")),
                    plotlyOutput("doc_samples_plot", height = "420px")
                  ),
                  card(
                    full_screen = TRUE,
                    card_header(span(fa("layer-group"), " Population groups per reference panel")),
                    plotlyOutput("doc_pops_count_plot", height = "420px")
                  ),
                  card(
                    full_screen = TRUE,
                    card_header(span(fa("dna"), " Unique rsIDs / variants by reference panel")),
                    plotlyOutput("doc_variants_plot", height = "520px")
                  )
                )
              )
            )
          )
      )
    )
  ),
  
  # ---------- About Tab (New) ----------
  nav_panel(
    title = span(fa("info-circle"), " About"),
    value = "about_tab",
    fluidPage(
      br(),
      fluidRow(
        column(8, offset = 2,
               card(
                 class = "shadow-lg",
                 card_header(span(fa("users"), " Credits & Contact")),
                 card_body(
                   div(style = "text-align: center; padding: 20px;",
                       h4("LDSeeker", style = "font-family: 'Rajdhani', sans-serif; font-weight: 700; color: #2c3e50; margin-bottom: 30px;"),
                       
                       h5("Main Developer"),
                       p("Georgios A. Manios", style="font-size: 1.1rem; font-weight: 500;"),
                       p("Computational Genetics Group", style="color: #718096;"),
                       tags$a(href = "mailto:gmanios@uth.gr", "gmanios@uth.gr", style="color: #00d1b2; font-weight: 600; text-decoration: none;"),
                       
                       hr(style="width: 60%; margin: 20px auto; border-top: 1px solid #e2e8f0;"),
                       
                       h5("Correspondence"),
                       p("Dr. Pantelis G. Bagos"),
                       tags$a(href = "mailto:pbagos@compgen.org", "pbagos@compgen.org", style="color: #00d1b2; font-weight: 600; text-decoration: none;"),
                       
                       hr(style="width: 60%; margin: 20px auto; border-top: 1px solid #e2e8f0;"),
                       
                       p(tags$a(href = "https://github.com/pbagos/LDSeeker/", 
                                span(fa("github"), " Visit GitHub Repository"), 
                                target = "_blank", 
                                class="btn btn-outline-dark"))
                   )
                 )
               ),
               br(),
               card(
                 card_header(span(fa("download"), " Download Reference panels and resources")),
                 card_body(
                   tags$div(
                     class = "table-responsive",
                     tags$table(
                       class = "table table-sm table-hover align-middle",
                       tags$thead(
                         tags$tr(
                           tags$th("Resource"),
                           tags$th("Web address")
                         )
                       ),
                       tags$tbody(
                         tags$tr(tags$td("LDSeeker GitHub repository"), tags$td(tags$a(href = "https://github.com/pbagos/LDSeeker/", "https://github.com/pbagos/LDSeeker/", target = "_blank"))),
                         tags$tr(tags$td("LDSeeker web tool"), tags$td(tags$a(href = "https://compgen.dib.uth.gr/LDSeeker/", "https://compgen.dib.uth.gr/LDSeeker/", target = "_blank"))),
                         tags$tr(tags$td("UKBB LD reference panel"), tags$td(tags$a(href = "https://zenodo.org/records/18847278", "https://zenodo.org/records/18847278", target = "_blank"))),
                         tags$tr(tags$td("1000 Genomes Project – high coverage LD reference panel"), tags$td(tags$a(href = "https://zenodo.org/records/18849097", "https://zenodo.org/records/18849097", target = "_blank"))),
                         tags$tr(tags$td("1000 Genomes Project LD reference panel"), tags$td(tags$a(href = "https://zenodo.org/records/18861527", "https://zenodo.org/records/18861527", target = "_blank"))),
                         tags$tr(tags$td("TOP-LD LD reference panel"), tags$td(tags$a(href = "http://195.251.108.185/ref_panels/TOP_LD", "http://195.251.108.185/ref_panels/TOP_LD", target = "_blank"))),
                         tags$tr(tags$td("HapMap LD reference panel"), tags$td(tags$a(href = "https://zenodo.org/records/20213914", "https://zenodo.org/records/20213914", target = "_blank"))),
                         tags$tr(tags$td("LASI-DAD LD reference panel"), tags$td(tags$a(href = "https://zenodo.org/records/20256884", "https://zenodo.org/records/20256884", target = "_blank"))),
                         tags$tr(tags$td("HGDP LD reference panel"), tags$td(tags$a(href = "http://195.251.108.185/HGDP/", "http://195.251.108.185/HGDP/", target = "_blank"))),
                         tags$tr(tags$td("SNP to Gene mapping reference"), tags$td(tags$a(href = "https://zenodo.org/records/20161132", "https://zenodo.org/records/20161132", target = "_blank")))
                       )
                     )
                   )
                 )
               )
        )
      )
    )
  )
)

# --- SERVER ---
server <- function(input, output, session) {
  rv <- reactiveValues(
    ld = NULL, 
    input_preview = NULL, 
    assoc_data = NULL, 
    ldassoc_res = NULL,
    prune_data = NULL,
    pruned_result = NULL,
    snp_gene_data = NULL,
    snp_gene_result = NULL,
    snp_gene_map = NULL,
    snp_gene_selected_cols = NULL,
    analysis_file_path = NULL, # Tracks if demo or uploaded file is active in Analysis tab
    heatmap_data = NULL # Stores loaded heatmap data
  )
  
  # Navigation Logic
  observeEvent(input$btn_goto_analysis, { nav_select("navbar_id", "analysis_tab") })
  observeEvent(input$btn_goto_pruning, { nav_select("navbar_id", "pruning_tab") })
  observeEvent(input$btn_goto_heatmap, { nav_select("navbar_id", "viz_tab") })
  observeEvent(input$btn_goto_ldassoc, { nav_select("navbar_id", "ldassoc_tab") })
  
  # Dynamic Population Picker Logic (General)
  observeEvent(input$ref_panel, {
    if (input$ref_panel == "Hap_Map") {
      updateSelectInput(session, "population", label = "Population (HapMap)", choices = hapmap_pops, selected = hapmap_pops[1])
    } else if (input$ref_panel == "HGDP") {
      updateSelectInput(session, "population", label = "Population (HGDP)", choices = HGDP_pops, selected = HGDP_pops[1])
    } else if (input$ref_panel == "UKBB") {
      updateSelectInput(session, "population", label = "Population (UKBB)", choices = UKBB_pops, selected = UKBB_pops[1])
    } else if (input$ref_panel == "LASI_DAD") {
      updateSelectInput(session, "population", label = "Population (LASI-DAD)", choices = LASI_DAD_pops, selected = LASI_DAD_pops[1])
    } else {
      updateSelectInput(session, "population", label = "Population Ancestry", choices = default_pops, selected = default_pops[1])
    }
  })
  
  # Dynamic Population Picker Logic (LDassoc)
  observeEvent(input$assoc_ref_panel, {
    if (input$assoc_ref_panel == "Hap_Map") {
      updateSelectInput(session, "assoc_population", label = "Population (HapMap)", choices = hapmap_pops, selected = hapmap_pops[1])
    } else if (input$assoc_ref_panel == "HGDP") {
      updateSelectInput(session, "assoc_population", label = "Population (HGDP)", choices = HGDP_pops, selected = HGDP_pops[1])
    } else if (input$assoc_ref_panel == "UKBB") {
      updateSelectInput(session, "assoc_population", label = "Population (UKBB)", choices = UKBB_pops, selected = UKBB_pops[1])
    } else if (input$assoc_ref_panel == "LASI_DAD") {
      updateSelectInput(session, "assoc_population", label = "Population (LASI-DAD)", choices = LASI_DAD_pops, selected = LASI_DAD_pops[1])
    } else {
      updateSelectInput(session, "assoc_population", label = "Population Ancestry", choices = default_pops, selected = default_pops[1])
    }
  })
  
  # Dynamic Population Picker Logic (Pruning)
  observeEvent(input$prune_ref_panel, {
    if (input$prune_ref_panel == "Hap_Map") {
      updateSelectInput(session, "prune_population", label = "Population (HapMap)", choices = hapmap_pops, selected = hapmap_pops[1])
    } else if (input$prune_ref_panel == "HGDP") {
      updateSelectInput(session, "prune_population", label = "Population (HGDP)", choices = HGDP_pops, selected = HGDP_pops[1])
    } else if (input$prune_ref_panel == "UKBB") {
      updateSelectInput(session, "prune_population", label = "Population (UKBB)", choices = UKBB_pops, selected = UKBB_pops[1])
    } else if (input$prune_ref_panel == "LASI_DAD") {
      updateSelectInput(session, "prune_population", label = "Population (LASI-DAD)", choices = LASI_DAD_pops, selected = LASI_DAD_pops[1])
    } else {
      updateSelectInput(session, "prune_population", label = "Population Ancestry", choices = default_pops, selected = default_pops[1])
    }
  })
  
  # --- LDpatterns LOGIC ---
  
  # Input File Preview Logic (File Upload)
  observeEvent(input$input_file, {
    req(input$input_file)
    updateTextAreaInput(session, "paste_rsids", value = "")
    tryCatch({
      rv$input_preview <- fread(input$input_file$datapath)
      rv$analysis_file_path <- NULL # Reset custom path so input_file path is used by default in execution logic
    }, error = function(e) {
      rv$input_preview <- data.frame(Error = paste("Could not read file:", e$message))
    })
  })
  
  # Load Example Data Logic (LDpatterns)
  observeEvent(input$load_example, {
    # Load from the local file written at startup
    if(file.exists("LD_patterns_demo.txt")){
      rv$input_preview <- fread("LD_patterns_demo.txt")
      rv$analysis_file_path <- "LD_patterns_demo.txt" # Flag to use this file
      # Clear paste input to avoid confusion
      updateTextAreaInput(session, "paste_rsids", value = "")
    } else {
      showNotification("Demo file not found.", type = "error")
    }
  })
  
  # Input File Preview Logic (Text Paste)
  observeEvent(input$paste_rsids, {
    req(input$paste_rsids)
    txt <- input$paste_rsids
    
    # Only process if there is actual content
    if (grepl("\\w", txt)) {
      tryCatch({
        clean_txt <- gsub("[, ]+", "\n", txt)
        df <- fread(text = clean_txt, header = FALSE)
        if (ncol(df) == 1) names(df) <- "SNP"
        rv$input_preview <- df
        rv$analysis_file_path <- NULL
      }, error = function(e) {
        rv$input_preview <- data.frame(Error = paste("Could not parse text:", e$message))
      })
    }
  })
  
  output$preview_section <- renderUI({
    req(rv$input_preview)
    card(
      card_header(span(fa("file-alt"), " Input File Preview")),
      DTOutput("preview_table")
    )
  })
  
  output$preview_table <- renderDT({
    req(rv$input_preview)
    datatable(rv$input_preview, style = "bootstrap5", options = list(pageLength = 5, scrollX = TRUE))
  })
  
  # Statistics Calculations
  output$val_variants <- renderUI({
    if (is.null(rv$ld) || nrow(rv$ld) == 0) return("0")
    res <- rv$ld
    cols <- names(res)
    unique_snps <- 0
    if (all(c("rsID1", "rsID2") %in% cols)) {
      unique_snps <- length(unique(c(res$rsID1, res$rsID2)))
    } else if (all(c("SNP_A", "SNP_B") %in% cols)) {
      unique_snps <- length(unique(c(res$SNP_A, res$SNP_B)))
    } else {
      unique_snps <- nrow(res)
    }
    formatC(unique_snps, big.mark = ",")
  })
  
  output$val_mean_r2 <- renderUI({
    if (is.null(rv$ld)) return("0.000")
    r2_col <- grep("^R2$", names(rv$ld), value = TRUE, ignore.case = TRUE)
    if (length(r2_col) > 0) sprintf("%.3f", mean(rv$ld[[r2_col[1]]], na.rm = TRUE)) else "0.000"
  })
  
  output$val_high_ld <- renderUI({
    if (is.null(rv$ld)) return("0")
    r2_col <- grep("^R2$", names(rv$ld), value = TRUE, ignore.case = TRUE)
    if (length(r2_col) > 0) formatC(sum(rv$ld[[r2_col[1]]] >= 0.9, na.rm = TRUE), big.mark = ",") else "0"
  })
  
  # LDpatterns Execution
  observeEvent(input$run_ld, {
    # Check if we have a file uploaded, a demo file loaded, or text pasted
    has_file <- !is.null(input$input_file)
    has_demo <- !is.null(rv$analysis_file_path)
    has_text <- nzchar(input$paste_rsids)
    
    if (!has_file && !has_demo && !has_text) {
      showNotification("Please upload a file, load example data, or paste rsIDs first.", type = "error")
      return()
    }
    
    # rv$input_preview <- NULL # Don't clear preview, useful context
    file_path <- NULL
    
    withProgress(message = "Executing", value = 0, {
      incProgress(0.1, "Validating inputs...")
      
      if (has_demo) {
        # Use the demo file path directly
        file_path <- rv$analysis_file_path
      } else if (has_text) {
        tmp_file <- tempfile(fileext = ".txt")
        clean_txt <- gsub("[, ]+", "\n", input$paste_rsids)
        df <- fread(text = clean_txt, header = FALSE)
        if (ncol(df) == 1) names(df) <- "SNP"
        fwrite(df, tmp_file, sep = "\t")
        file_path <- tmp_file
      } else if (has_file) {
        file_path <- input$input_file$datapath
      }
      
      cmd <- paste("python LDSeeker.py", "--file-path", shQuote(file_path), "--r2threshold", input$r2threshold, "--pop", input$population, "--maf", input$maf, "--ref", input$ref_panel, "--pairwise", input$pairwise)
      print(cmd)
      incProgress(0.3, "Computing Linkage Disequilibrium...")
      system(cmd)
      
      incProgress(0.4, "Parsing result matrices...")
      file_name <- if (input$pairwise == "YES") "LD_info_chr_all_pairwise.txt" else "LD_info_chr_all.txt"
      if(file.exists(file_name)) rv$ld <- fread(file_name)
      incProgress(0.2, "Rendering view...")
    })
  })
  
  output$ld_table <- renderDT({
    req(rv$ld)
    datatable(rv$ld, extensions = 'Buttons', style = 'bootstrap5',
              options = list(pageLength = 10, scrollX = TRUE, dom = 'Bfrtip',
                             buttons = list('copy', 'csv', 'excel')
              )
    )
  }, server = FALSE)
  
  # --- LD MATRIX DOWNLOAD HANDLER ---
  output$download_ld_matrix <- downloadHandler(
    filename = function() {
      paste0("LD_Matrix_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(rv$ld)
      
      # Extract Data
      df <- copy(rv$ld)
      
      # Detect Columns
      col_s1 <- NULL
      col_s2 <- NULL
      col_val <- NULL
      
      if (all(c("rsID1", "rsID2") %in% names(df))) {
        col_s1 <- "rsID1"; col_s2 <- "rsID2"
      } else if (all(c("SNP_A", "SNP_B") %in% names(df))) {
        col_s1 <- "SNP_A"; col_s2 <- "SNP_B"
      }
      
      # Find R2 column (case insensitive)
      r2_candidates <- grep("^R2$", names(df), value = TRUE, ignore.case = TRUE)
      if (length(r2_candidates) > 0) col_val <- r2_candidates[1]
      
      if (is.null(col_s1) || is.null(col_s2) || is.null(col_val)) {
        write.csv(data.frame(Error = "Data format not suitable for matrix conversion (Missing SNP pairs or R2)."), file, row.names = FALSE)
        return()
      }
      
      # Select and Rename
      dat <- df[, c(col_s1, col_s2, col_val), with = FALSE]
      names(dat) <- c("S1", "S2", "Val")
      
      # Symmetrize (ensure A-B and B-A exist)
      dat_rev <- dat[, .(S1 = S2, S2 = S1, Val = Val)]
      dat_full <- unique(rbind(dat, dat_rev))
      
      # Cast to Matrix using data.table's dcast
      mat_dt <- dcast(dat_full, S1 ~ S2, value.var = "Val", fill = 0)
      
      # Fill Diagonal with 1 if missing
      snps <- names(mat_dt)[-1] # First col is S1
      for (s in snps) {
        if (s %in% mat_dt$S1) {
          row_idx <- which(mat_dt$S1 == s)
          set(mat_dt, i = row_idx, j = s, value = 1.0)
        }
      }
      
      write.csv(mat_dt, file, row.names = FALSE)
    }
  )
  
  # --- LDpruning LOGIC ---
  
  # 1. Read Pruning Input File and UI
  observeEvent(input$prune_file, {
    req(input$prune_file)
    tryCatch({
      df <- fread(input$prune_file$datapath)
      rv$prune_data <- df
      col_names <- names(df)
      
      first_matching_col <- function(pattern, fallback = NA_character_) {
        hits <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
        if (length(hits) > 0) hits[1] else fallback
      }
      none_if_missing <- function(x) {
        if (length(x) == 0 || is.na(x) || is.null(x)) "None" else x
      }
      
      # Robust defaults for common GWAS column names
      default_snp  <- first_matching_col("^(SNP|RSID|RS_ID)$|rs|snp|variant", col_names[1])
      default_chr  <- first_matching_col("^(CHR|CHROM|CHROMOSOME)$|chr|chrom", col_names[1])
      default_pval <- first_matching_col("^(P|PVALUE|P_VALUE|PVAL)$|p.?val")
      default_z    <- first_matching_col("^(Z|ZSCORE|Z_SCORE)$|z.?score")
      default_beta <- first_matching_col("^(BETA|EFFECT|B)$|beta|effect")
      default_se   <- first_matching_col("^(SE|STDERR|STD_ERROR|STANDARD_ERROR)$|std.?err|standard.?error")
      
      output$prune_cols_ui <- renderUI({
        tagList(
          selectInput("prune_col_snp", "SNP/rsID Column", choices = col_names, selected = default_snp),
          selectInput("prune_col_chr", "Chromosome Column", choices = col_names, selected = default_chr),
          selectInput("prune_col_pval", "P-value Column (Optional)", choices = c("None", col_names), selected = none_if_missing(default_pval)),
          selectInput("prune_col_z", "Z-score Column (Optional)", choices = c("None", col_names), selected = none_if_missing(default_z)),
          selectInput("prune_col_beta", "BETA Column (Optional)", choices = c("None", col_names), selected = none_if_missing(default_beta)),
          selectInput("prune_col_se", "SE Column (Optional)", choices = c("None", col_names), selected = none_if_missing(default_se))
        )
      })
      
      # Preview Data
      output$prune_preview_section <- renderUI({
        card(
          card_header(span(fa("file-alt"), " Raw Input Preview")),
          DTOutput("prune_preview_table")
        )
      })
      
      output$prune_preview_table <- renderDT({
        req(rv$prune_data)
        datatable(rv$prune_data, style = "bootstrap5", options = list(pageLength = 5, scrollX = TRUE))
      })
      
    }, error = function(e) {
      showNotification(paste("Error reading pruning file:", e$message), type = "error")
    })
  })
  
  # 2. Execute LDSeeker pruning through the Python backend
  observeEvent(input$run_pruning, {
    req(input$prune_file, rv$prune_data, input$prune_col_snp, input$prune_col_chr)
    
    withProgress(message = "Running LDSeeker LD pruning", value = 0, {
      tryCatch({
        incProgress(0.1, "Preparing GWAS input for LDSeeker...")
        
        df_python <- copy(rv$prune_data)
        
        # Standardize the columns required by LDSeeker.py.
        # The original columns are preserved, while canonical SNP/CHR/P/Z/BETA/SE
        # columns are added or overwritten for the Python argument names.
        df_python[, SNP := as.character(get(input$prune_col_snp))]
        df_python[, CHR := get(input$prune_col_chr)]
        
        p_col_for_python <- "__NO_P_COLUMN__"
        z_col_for_python <- "__NO_Z_COLUMN__"
        beta_col_for_python <- "__NO_BETA_COLUMN__"
        se_col_for_python <- "__NO_SE_COLUMN__"
        
        if (!is.null(input$prune_col_pval) && input$prune_col_pval != "None") {
          df_python[, P := as.numeric(as.character(get(input$prune_col_pval)))]
          p_col_for_python <- "P"
        }
        if (!is.null(input$prune_col_z) && input$prune_col_z != "None") {
          df_python[, Z := as.numeric(as.character(get(input$prune_col_z)))]
          z_col_for_python <- "Z"
        }
        if (!is.null(input$prune_col_beta) && input$prune_col_beta != "None") {
          df_python[, BETA := as.numeric(as.character(get(input$prune_col_beta)))]
          beta_col_for_python <- "BETA"
        }
        if (!is.null(input$prune_col_se) && input$prune_col_se != "None") {
          df_python[, SE := as.numeric(as.character(get(input$prune_col_se)))]
          se_col_for_python <- "SE"
        }
        
        # Remove rows with missing essential identifiers before calling LDSeeker.
        df_python <- df_python[!is.na(SNP) & SNP != "" & !is.na(CHR) & CHR != ""]
        if (nrow(df_python) == 0) {
          showNotification("No variants with non-missing SNP and CHR columns were found.", type = "warning")
          return()
        }
        
        tmp_prune_input <- tempfile(fileext = ".txt")
        fwrite(df_python, tmp_prune_input, sep = "\t", quote = FALSE, na = "NA")
        
        # Clear old LDSeeker outputs so the UI never displays stale pruning results.
        unlink(c("LD_info_chr_all_pairwise.txt", "LD_pruned_kept.txt", "LD_pruned_tmp_batch_*.txt"))
        
        incProgress(0.25, "Executing LDSeeker.py with --ld-prune YES...")
        
        python_cmd <- Sys.which("python")
        if (python_cmd == "") python_cmd <- Sys.which("python3")
        if (python_cmd == "") stop("Could not find python or python3 in PATH.")
        
        cmd_args <- c(
          "LDSeeker.py",
          "--file-path", tmp_prune_input,
          "--r2threshold", as.character(input$prune_r2_thresh),
          "--pop", input$prune_population,
          "--maf", as.character(input$prune_maf),
          "--ref", input$prune_ref_panel,
          "--pairwise", "YES",
          "--ld-prune", "YES",
          "--ld-prune-prefix", "LD_pruned",
          "--ld-prune-metric", input$prune_metric,
          "--ld-prune-p-col", p_col_for_python,
          "--ld-prune-z-col", z_col_for_python
        )
        
        if (!is.null(input$prune_pval_dir) && input$prune_pval_dir != "none" && p_col_for_python == "P") {
          cmd_args <- c(
            cmd_args,
            "--ld-prune-p-threshold", as.character(input$prune_pval_thresh),
            "--ld-prune-p-threshold-mode", input$prune_pval_dir
          )
        }
        
        if (z_col_for_python == "Z") {
          cmd_args <- c(cmd_args, "--ld-prune-z-threshold", as.character(input$prune_z_thresh))
        }
        
        if (!is.null(input$prune_significance) && input$prune_significance != "none") {
          cmd_args <- c(
            cmd_args,
            "--significance", input$prune_significance,
            "--significance-metric", input$prune_sig_metric,
            "--significance-threshold", as.character(input$prune_sig_thresh),
            "--sig-p-col", p_col_for_python,
            "--sig-z-col", z_col_for_python,
            "--beta-col", beta_col_for_python,
            "--se-col", se_col_for_python
          )
        }
        
        message("Running command: ", paste(c(python_cmd, shQuote(cmd_args)), collapse = " "))
        cmd_output <- system2(python_cmd, args = cmd_args, stdout = TRUE, stderr = TRUE)
        message(paste(cmd_output, collapse = "\n"))
        
        incProgress(0.55, "Reading LD_pruned_kept.txt...")
        pruned_file <- "LD_pruned_kept.txt"
        if (!file.exists(pruned_file)) {
          rv$pruned_result <- NULL
          showNotification("LDSeeker finished, but LD_pruned_kept.txt was not found.", type = "error")
          return()
        }
        
        rv$pruned_result <- fread(pruned_file)
        
        
        initial_variants <- length(unique(df_python$SNP))
        final_snp_col <- if ("SNP" %in% names(rv$pruned_result)) "SNP" else names(rv$pruned_result)[1]
        final_variants <- length(unique(rv$pruned_result[[final_snp_col]]))
        removed_variants <- max(initial_variants - final_variants, 0)
        
        output$val_prune_initial <- renderUI(formatC(initial_variants, big.mark = ","))
        output$val_prune_removed <- renderUI(formatC(removed_variants, big.mark = ","))
        output$val_prune_final <- renderUI(formatC(final_variants, big.mark = ","))
        
        incProgress(1, "Done!")
        showNotification("LD pruning completed. LD_pruned_kept.txt is displayed in the pruning results.", type = "message")
        
      }, error = function(e) {
        showNotification(paste("LD pruning failed:", e$message), type = "error")
      })
    })
  })
  
  output$pruned_table <- renderDT({
    req(rv$pruned_result)
    datatable(rv$pruned_result, extensions = 'Buttons', style = 'bootstrap5',
              options = list(pageLength = 15, scrollX = TRUE, dom = 'Bfrtip',
                             buttons = list('copy', 'csv', 'excel')))
  })
  
  output$download_pruned <- downloadHandler(
    filename = function() { paste0("LD_pruned_kept_", format(Sys.Date(), "%Y%m%d"), ".txt") },
    content = function(file) {
      req(rv$pruned_result)
      fwrite(rv$pruned_result, file, sep = "\t")
    }
  )
  
  # --- LDASSOC LOGIC ---
  
  # 1. Read Assoc File and Generate Dynamic Selectors
  observeEvent(input$assoc_file, {
    req(input$assoc_file)
    tryCatch({
      df <- fread(input$assoc_file$datapath)
      rv$assoc_data <- df
      col_names <- names(df)
      
      # Heuristics for defaults
      default_chr <- col_names[grep("chr", col_names, ignore.case = TRUE)[1]]
      default_pval <- col_names[grep("p.?val", col_names, ignore.case = TRUE)[1]]
      default_snp <- col_names[grep("rs|snp|variant", col_names, ignore.case = TRUE)[1]]
      default_bp <- col_names[grep("bp|pos|location", col_names, ignore.case = TRUE)[1]]
      
      output$assoc_columns_ui <- renderUI({
        tagList(
          selectInput("assoc_col_chr", "Chromosome Column", choices = col_names, selected = default_chr),
          selectInput("assoc_col_pval", "P-value Column", choices = col_names, selected = default_pval),
          selectInput("assoc_col_snp", "SNP/rsID Column", choices = col_names, selected = default_snp),
          selectInput("assoc_col_bp", "BP/Position Column", choices = col_names, selected = default_bp)
        )
      })
    }, error = function(e) {
      showNotification(paste("Error reading association file:", e$message), type = "error")
    })
  })
  
  # Load Example for LDassoc
  observeEvent(input$load_example_assoc, {
    if(file.exists("LD_assoc_demo.txt")){
      df <- fread("LD_assoc_demo.txt")
      rv$assoc_data <- df
      col_names <- names(df)
      
      # Pre-fill inputs based on demo data
      updateTextInput(session, "assoc_index_snp", value = "rs7837688")
      updateTextInput(session, "assoc_target_chr", value = "8")
      
      # Manually set column selectors for the demo file structure
      output$assoc_columns_ui <- renderUI({
        tagList(
          selectInput("assoc_col_chr", "Chromosome Column", choices = col_names, selected = "CHR"),
          selectInput("assoc_col_pval", "P-value Column", choices = col_names, selected = "P-value"),
          selectInput("assoc_col_snp", "SNP/rsID Column", choices = col_names, selected = "SNP"),
          selectInput("assoc_col_bp", "BP/Position Column", choices = col_names, selected = "BP")
        )
      })
    } else {
      showNotification("Demo file not found.", type = "error")
    }
  })
  
  output$assoc_preview_table <- renderDT({
    req(rv$assoc_data)
    datatable(rv$assoc_data, style = "bootstrap5", options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # 2. Execute LDassoc
  observeEvent(input$run_ldassoc, {
    req(input$assoc_index_snp, input$assoc_target_chr, rv$assoc_data)
    # Note: Removed direct req(input$assoc_file) to allow running with demo data
    req(input$assoc_col_chr, input$assoc_col_pval, input$assoc_col_snp, input$assoc_col_bp)
    
    withProgress(message = "Processing LDassoc", value = 0, {
      
      # A. Create Temp Input File using FULL Association Data
      incProgress(0.1, "Preparing input file...")
      tmp_index_file <- tempfile(fileext = ".txt")
      
      # Copy data and prepare specific columns for LDSeeker (SNP and CHR are required)
      input_df <- copy(rv$assoc_data)
      ld_input <- input_df[, c(input$assoc_col_snp, input$assoc_col_chr), with = FALSE]
      setnames(ld_input, c(input$assoc_col_snp, input$assoc_col_chr), c("SNP", "CHR"))
      
      fwrite(ld_input, tmp_index_file, sep = "\t")
      
      # B. Run LDSeeker (Pairwise=YES on the FULL file)
      incProgress(0.3, "Computing Pairwise LD...")
      cmd <- paste("python LDSeeker.py", 
                   "--file-path", shQuote(tmp_index_file), 
                   "--r2threshold", input$assoc_r2, 
                   "--pop", input$assoc_population, 
                   "--maf", input$assoc_maf, 
                   "--ref", input$assoc_ref_panel, 
                   "--pairwise", "YES") 
      
      system(cmd)
      
      # C. Read Results, Filter for Index, and Merge
      incProgress(0.5, "Filtering and Merging...")
      ld_result_file <- "LD_info_chr_all_pairwise.txt" 
      
      if(file.exists(ld_result_file)) {
        ld_res <- fread(ld_result_file)
        
        # --- NORMALIZE MAF COLUMN NAME ---
        maf_col <- grep("^MAF$", names(ld_res), value = TRUE, ignore.case = TRUE)
        if(length(maf_col) > 0) setnames(ld_res, maf_col[1], "MAF")
        
        index_snp <- input$assoc_index_snp
        
        # Determine column names for SNP pairs (usually rsID1/rsID2 or SNP_A/SNP_B)
        c1 <- if("rsID1" %in% names(ld_res)) "rsID1" else "SNP_A"
        c2 <- if("rsID2" %in% names(ld_res)) "rsID2" else "SNP_B"
        
        if(c1 %in% names(ld_res) && c2 %in% names(ld_res)) {
          # FILTER: Keep rows where Index SNP matches either column
          ld_res <- ld_res[get(c1) == index_snp | get(c2) == index_snp]
          
          # Identify the "Other" SNP to merge with Association Data
          ld_res[, RS_Number := ifelse(get(c1) == index_snp, get(c2), get(c1))]
        } else {
          showNotification("Could not identify SNP columns in LD output", type = "error")
        }
        
        # D. Merge with Association Data (P-values & BP)
        assoc_df <- rv$assoc_data
        
        # Prepare subset of Association Data (Added BP column)
        keep_cols <- c(input$assoc_col_snp, input$assoc_col_pval, input$assoc_col_bp)
        assoc_subset <- assoc_df[, ..keep_cols]
        names(assoc_subset) <- c("RS_Number", "P_value_Input", "BP")
        
        # Merge
        merged_df <- merge(ld_res, assoc_subset, by = "RS_Number", all.x = TRUE)
        
        # --- ENSURE INDEX SNP IS INCLUDED ---
        if (!index_snp %in% merged_df$RS_Number) {
          # Find index SNP in original association data
          index_row <- assoc_subset[RS_Number == index_snp]
          
          if(nrow(index_row) > 0) {
            # Construct a row for the index SNP
            # We use fill=TRUE in rbind to handle missing columns
            new_row <- data.table(
              RS_Number = index_snp,
              P_value_Input = index_row$P_value_Input,
              BP = index_row$BP,
              R2 = 1.0, # Index SNP has perfect LD with itself
              D_prime = 1.0,
              MAF = NA, # Unknown from just association data, unless in ld_res
              Distance = 0
            )
            # Append to results
            merged_df <- rbind(merged_df, new_row, fill=TRUE)
          }
        }
        
        # --- DEDUPLICATE RESULTS ---
        merged_df <- unique(merged_df, by = "RS_Number")
        
        # --- CALCULATE DISTANCE (Abs Diff from Index SNP) ---
        # Get Index SNP Position from the merged data
        index_pos_val <- merged_df[RS_Number == index_snp, BP]
        
        if (length(index_pos_val) > 0 && !all(is.na(index_pos_val))) {
          # Use first occurrence if duplicates exist and ensure numeric
          ref_bp <- as.numeric(index_pos_val[1])
          merged_df[, Distance := abs(as.numeric(BP) - ref_bp)]
        } else {
          merged_df[, Distance := NA]
        }
        
        # Standardize Dprime column name if needed
        if("D_prime" %in% names(merged_df)) setnames(merged_df, "D_prime", "Dprime")
        
        # E. Formatting Final Table
        setnames(merged_df, "P_value_Input", "P-value")
        
        # Ensure standard columns exist for display
        if(!"Coord" %in% names(merged_df) && "BP" %in% names(merged_df)) merged_df[, Coord := paste0("chr", input$assoc_target_chr, ":", BP)]
        if(!"Coord" %in% names(merged_df)) merged_df[, Coord := "NA"]
        
        # Define preferred order, including MAF if present
        desired_order <- c("RS_Number", "BP", "Coord", "Alleles", "MAF", "Distance", "Dprime", "R2", "P-value")
        
        # Intersect will keep MAF if it exists in merged_df, and drop it if it doesn't
        final_cols <- intersect(desired_order, names(merged_df))
        rv$ldassoc_res <- merged_df[, ..final_cols]
        
      } else {
        showNotification("LD Analysis failed or returned no results.", type = "warning")
        rv$ldassoc_res <- data.frame(Status = "No LD results found.")
      }
      
      incProgress(0.2, "Done.")
    })
  })
  
  # 3. Interactive Plot Logic
  output$ldassoc_plot <- renderPlotly({
    req(rv$ldassoc_res)
    df <- rv$ldassoc_res
    
    # Check if necessary columns exist
    if(!all(c("BP", "P-value", "R2", "RS_Number") %in% names(df))) {
      return(plot_ly() %>% layout(title = "Data incomplete for plotting (Missing BP, P-value, or R2)"))
    }
    
    # Ensure numeric
    df$BP <- as.numeric(df$BP)
    df$`P-value` <- as.numeric(df$`P-value`)
    df$R2 <- as.numeric(df$R2)
    df$logP <- -log10(df$`P-value`)
    
    # Assign marker size based on Index SNP
    df$marker_size <- ifelse(df$RS_Number == input$assoc_index_snp, 15, 10)
    
    # Tooltip
    df$hover <- paste0(
      "SNP: ", df$RS_Number, "<br>",
      "BP: ", df$BP, "<br>",
      "P-value: ", formatC(df$`P-value`, format="e", digits=3), "<br>",
      "R2: ", round(df$R2, 3), "<br>",
      if("Dprime" %in% names(df)) paste0("D': ", round(df$Dprime, 3)) else ""
    )
    
    plot_ly(df, x = ~BP, y = ~logP, 
            type = 'scatter', mode = 'markers',
            marker = list(
              size = ~marker_size, 
              color = ~R2, 
              colorscale = input$ldassoc_color_scale, # Use user selection
              colorbar = list(title = "R2"),
              showscale = TRUE,
              line = list(color = 'grey', width = 1)
            ),
            text = ~hover, hoverinfo = "text"
    ) %>%
      layout(
        title = paste("Regional Association & LD (Index:", input$assoc_index_snp, ")"),
        xaxis = list(title = "Chromosome Position (BP)"),
        yaxis = list(title = "-log10(P-value)"),
        hovermode = "closest"
      )
  })
  
  output$ldassoc_result_table <- renderDT({
    req(rv$ldassoc_res)
    datatable(rv$ldassoc_res, extensions = 'Buttons', style = 'bootstrap5',
              options = list(pageLength = 10, scrollX = TRUE, dom = 'Bfrtip',
                             buttons = list('copy', 'csv', 'excel')))
  })
  
  
  # --- HEATMAP LOGIC ---
  
  # Load Uploaded File
  observeEvent(input$ld_heat_file, {
    req(input$ld_heat_file)
    tryCatch({
      df <- fread(input$ld_heat_file$datapath, header = TRUE, sep = "\t", data.table = FALSE)
      rv$heatmap_data <- df
    }, error = function(e) {
      showNotification(paste("Error reading heatmap file:", e$message), type = "error")
    })
  })
  
  # Load Demo File
  observeEvent(input$load_example_heatmap, {
    if(file.exists("LDresults_demo.txt")){
      df <- fread("LDresults_demo.txt", header = TRUE, sep = "\t", data.table = FALSE)
      rv$heatmap_data <- df
    } else {
      showNotification("Demo file not found.", type = "error")
    }
  })
  
  # Process Heatmap Data (from reactive value)
  ld_heat_data <- reactive({
    req(rv$heatmap_data)
    df <- rv$heatmap_data
    
    # Normalization of column names
    if ("SNP_A" %in% names(df)) names(df)[names(df) == "SNP_A"] <- "rsID1"
    if ("SNP_B" %in% names(df)) names(df)[names(df) == "SNP_B"] <- "rsID2"
    if ("D_prime" %in% names(df)) names(df)[names(df) == "D_prime"] <- "Dprime"
    if ("r" %in% names(df)) names(df)[names(df) == "r"] <- "R"
    if ("d" %in% names(df)) names(df)[names(df) == "d"] <- "D"
    
    # Standardize CHR
    chr_col <- grep("^CHR$|^Chromosome$", names(df), value = TRUE, ignore.case = TRUE)
    if (length(chr_col) > 0) names(df)[names(df) == chr_col[1]] <- "CHR"
    
    validate(need(all(c("rsID1", "rsID2") %in% names(df)), "Missing rsID1 or rsID2 columns."))
    df
  })
  
  # Filtered Heatmap Data based on Chromosome selection
  ld_heat_data_filtered <- reactive({
    req(ld_heat_data())
    df <- ld_heat_data()
    
    if ("CHR" %in% names(df)) {
      # If multiple chromosomes exist, select based on input
      if (!is.null(input$heatmap_filter_chr) && input$heatmap_filter_chr %in% df$CHR) {
        df <- df[df$CHR == input$heatmap_filter_chr, ]
      } else {
        # Default to first available chromosome if input is not ready/valid
        chrs <- unique(df$CHR)
        if (length(chrs) > 1) {
          df <- df[df$CHR == chrs[1], ]
        }
      }
    }
    df
  })
  
  # Shared Batched/Subsampled Reactive for Visualizations
  ld_viz_subset <- reactive({
    req(input$ld_metric, ld_heat_data_filtered())
    df_full <- ld_heat_data_filtered()
    
    limit <- input$n_rows
    page <- if(!is.null(input$heatmap_page)) input$heatmap_page else 1
    start_idx <- (page - 1) * limit + 1
    end_idx <- min(start_idx + limit - 1, nrow(df_full))
    
    df_full[start_idx:end_idx, ]
  })
  
  # Dynamic UI for Chromosome Selection
  output$heatmap_chr_ui <- renderUI({
    req(ld_heat_data())
    df <- ld_heat_data()
    if ("CHR" %in% names(df)) {
      chrs <- unique(df$CHR)
      # Try numeric sort
      if(is.numeric(chrs)) chrs <- sort(chrs) else chrs <- sort(chrs)
      
      if (length(chrs) > 1) {
        selectInput("heatmap_filter_chr", "Select Chromosome", choices = chrs, selected = chrs[1])
      }
    }
  })
  
  output$heatmap_slider_ui <- renderUI({
    req(ld_heat_data_filtered())
    n_batches <- ceiling(nrow(ld_heat_data_filtered()) / input$n_rows)
    if (n_batches > 1) sliderInput("heatmap_page", "Select Batch (Page)", min = 1, max = n_batches, value = 1, step = 1)
  })
  
  # Dynamic Metric Selector UI
  output$ld_metric_ui <- renderUI({
    req(ld_heat_data_filtered())
    df <- ld_heat_data_filtered()
    choices <- c()
    if("R2" %in% names(df)) choices <- c(choices, "R2")
    if("Dprime" %in% names(df)) choices <- c(choices, "Dprime")
    if("R" %in% names(df)) choices <- c(choices, "R")
    if("D" %in% names(df)) choices <- c(choices, "D")
    
    if(length(choices) == 0) choices <- c("R2") # Fallback
    
    selectInput("ld_metric", "Metric", choices = choices, selected = choices[1])
  })
  
  # Auto-adjust slider based on metric selection
  observeEvent(input$ld_metric, {
    req(input$ld_metric)
    if (input$ld_metric %in% c("R", "D")) {
      updateSliderInput(session, "r2_range", label = "Intensity Range", min = -1, max = 1, value = c(-1, 1), step = 0.01)
    } else {
      updateSliderInput(session, "r2_range", label = "Intensity Range", min = 0, max = 1, value = c(0, 1), step = 0.01)
    }
  })
  
  heatmap_data_matrix <- reactive({
    req(ld_viz_subset())
    df <- ld_viz_subset()
    metric <- input$ld_metric
    
    validate(need(metric %in% names(df), paste("Column", metric, "not found in uploaded file.")))
    
    snps <- sort(unique(c(df$rsID1, df$rsID2)))
    mat <- matrix(NA, nrow = length(snps), ncol = length(snps), dimnames = list(snps, snps))
    
    for (i in 1:nrow(df)) {
      a <- as.character(df$rsID1[i]); b <- as.character(df$rsID2[i])
      val <- df[[metric]][i]
      mat[a, b] <- val; mat[b, a] <- val
    }
    mat[is.na(mat)] <- 0
    diag(mat) <- 1.00
    
    # Return full symmetric matrix for gaston and raw processing
    list(mat = mat, snps = snps)
  })
  
  # Plotly Heatmap
  output$heatmap <- renderPlotly({
    req(input$ld_metric)
    h_data <- heatmap_data_matrix()
    mat <- h_data$mat; snps <- h_data$snps
    metric <- input$ld_metric
    
    # Apply upper triangle masking LOCALLY for Plotly
    mat_masked <- mat
    for (i in 1:nrow(mat_masked)) { 
      for (j in 1:ncol(mat_masked)) { 
        if (j > i) mat_masked[i, j] <- NA 
      } 
    }
    
    hover_text <- matrix("", nrow = nrow(mat_masked), ncol = ncol(mat_masked))
    for (i in seq_along(snps)) {
      for (j in seq_along(snps)) {
        if (!is.na(mat_masked[i, j])) hover_text[i, j] <- paste0("rsID1: ", snps[i], "<br>rsID2: ", snps[j], "<br>", metric, ": ", round(mat_masked[i, j], 3))
      }
    }
    
    p <- plot_ly(x = snps, y = snps, z = mat_masked, type = "heatmap", colorscale = input$color_scale,
                 zmin = input$r2_range[1], zmax = input$r2_range[2], text = hover_text, hoverinfo = "text")
    
    if (input$show_values) {
      annots <- list()
      for (i in seq_along(snps)) {
        for (j in seq_along(snps)) {
          if (!is.na(mat_masked[i, j])) {
            annots[[length(annots) + 1]] <- list(x = snps[j], y = snps[i], text = sprintf("%.2f", mat_masked[i, j]),
                                                 showarrow = FALSE, font = list(color = if(mat_masked[i, j] > 0.5) "black" else "white", size = 10))
          }
        }
      }
      p <- layout(p, annotations = annots)
    }
    
    layout(p, xaxis = list(title = "", tickangle = -90, scaleanchor = "y"),
           yaxis = list(title = "", autorange = "reversed"),
           title = list(text = paste("Linkage Disequilibrium (", metric, ")"), font = list(size = 18, color = "#00d1b2")))
  })
  
  # Gaston LD Plot
  output$gaston_plot <- renderPlot({
    req(heatmap_data_matrix())
    h_data <- heatmap_data_matrix()
    mat <- h_data$mat # Symmetric matrix required by gaston
    
    # Check if matrix is valid for gaston
    validate(need(nrow(mat) > 1, "Not enough SNPs to plot."))
    
    # Gaston plot
    # Using dummy positions since simple matrix input doesn't carry BP
    gaston::LD.plot(mat, snp.positions = 1:nrow(mat))
  })
  
  
  
  # --- INTERACTIVE NETWORK LOGIC ---
  output$ld_network <- renderVisNetwork({
    req(input$ld_metric, ld_viz_subset())
    df_batch <- ld_viz_subset()
    metric <- input$ld_metric
    range <- input$r2_range
    
    # Palette Mapping: Extract colors for logic
    palette_colors <- switch(input$color_scale,
                             "Viridis" = c("#440154", "#21908C", "#FDE725"),
                             "Cividis" = c("#00204D", "#7C7B78", "#FFEA46"),
                             "Plasma"  = c("#0D0887", "#CC4678", "#F0F921"),
                             "Inferno" = c("#000004", "#BB3754", "#FCFDBF"),
                             "Magma"   = c("#000004", "#B4367A", "#FCFDBF"),
                             "Blues"   = c("#EFF3FF", "#6BAED6", "#084594"),
                             "RdBu"    = c("#053061", "#F7F7F7", "#67001F"),
                             "YlGnBu"  = c("#FFFFD9", "#41B6C4", "#081D58"),
                             "Greys"   = c("#FFFFFF", "#969696", "#000000"),
                             c("#00d1b2", "#2c3e50") # Default
    )
    
    # Filter current batch by user intensity range
    links_df <- as.data.table(df_batch)
    links <- links_df[get(metric) >= range[1] & get(metric) <= range[2]]
    
    validate(need(nrow(links) > 0, "No LD correlations found in the selected range for this batch."))
    
    # 1. Colorscale logic for Edges
    vals <- as.numeric(links[[metric]])
    
    # Define domain based on metric
    is_signed <- input$ld_metric %in% c("R", "D")
    domain_min <- if(is_signed) -1 else 0
    domain_max <- 1
    
    # Generate gradient
    col_func <- colorRampPalette(palette_colors)
    cols_vec <- col_func(100)
    
    # Map values to color indices
    norm_vals <- (vals - domain_min) / (domain_max - domain_min)
    norm_vals <- pmax(0, pmin(1, norm_vals)) # Clamp to [0,1]
    indices <- floor(norm_vals * 99) + 1
    edge_hex <- cols_vec[indices]
    
    # Nodes: Compute node size based on connectivity count within the current batch
    snps_all <- c(links$rsID1, links$rsID2)
    node_counts <- table(snps_all)
    unique_snps <- names(node_counts)
    
    nodes <- data.frame(
      id = unique_snps,
      label = unique_snps,
      title = paste0("SNP: ", unique_snps, "<br>Connectivity (Current Batch): ", as.numeric(node_counts)),
      value = as.numeric(node_counts) * 2 # Scale node size
      # Removed group="SNP" so global visNodes color applies
    )
    
    # Edge logic: normal size and gradient coloring
    edges <- data.frame(
      from = links$rsID1,
      to = links$rsID2,
      width = 2, # Normal edge size
      title = paste0(metric, ": ", round(as.numeric(links[[metric]]), 3)),
      color = edge_hex
    )
    
    visNetwork(nodes, edges) %>%
      visNodes(
        font = list(size = 14, color = "#34495e"), # Fixed dark color for text readability
        shape = input$net_node_shape,
        color = list(
          background = input$net_node_color,
          border = input$net_node_color,
          highlight = list(background = "#ffb703", border = "#ffb703")
        )
      ) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visInteraction(hover = TRUE, navigationButtons = TRUE) %>%
      visPhysics(enabled = FALSE) %>%
      visIgraphLayout(layout = input$net_layout, randomSeed = 123)
  })
  
  
  
  
  
  
  # --- SNP TO GENE LOGIC ---
  
  normalize_chr_values <- function(x) {
    x <- trimws(as.character(x))
    x <- toupper(gsub("^CHR", "", x, ignore.case = TRUE))
    x <- gsub("\\.0$", "", x)
    x[x %in% c("", "NA", "NAN", "NULL")] <- NA_character_
    x
  }
  
  detect_snp_identifier_cols <- function(col_names) {
    if (length(col_names) == 0) return(character(0))
    clean <- tolower(gsub("[^a-z0-9]", "", col_names))
    exact <- grepl("^(rsid|rs|snp)([0-9]+|[ab])?$", clean)
    common <- clean %in% c(
      "rsnumber", "rsidnumber", "rsidentifier",
      "snpid", "snpname", "snpidentifier",
      "variant", "variantid", "variantname",
      "marker", "markername", "markerid"
    )
    contains <- grepl("rsid|snp", clean)
    unique(col_names[exact | common | contains])
  }
  
  detect_chr_cols <- function(col_names) {
    if (length(col_names) == 0) return(character(0))
    clean <- tolower(gsub("[^a-z0-9]", "", col_names))
    # Matches CHR, CHROM, #CHROM, CHROMOSOME, CHROM_A/B, chr22, and descriptive
    # variants like chr_name, chr_no, chr_id, chrom_num, chromosome_number, etc.
    keep <- grepl("^(chr|chrom|chromosome)([0-9]+|[ab]|no|num|number|name|id|pos)*$", clean)
    unique(col_names[keep])
  }
  
  get_unique_snps_from_cols <- function(df, snp_cols) {
    snps <- unlist(lapply(snp_cols, function(cc) as.character(df[[cc]])), use.names = FALSE)
    snps <- trimws(snps)
    snps <- snps[!is.na(snps) & snps != ""]
    unique(snps)
  }
  
  get_chr_values_from_input <- function(df) {
    chr_cols <- detect_chr_cols(names(df))
    if (length(chr_cols) == 0) return(character(0))
    chrs <- unlist(lapply(chr_cols, function(cc) as.character(df[[cc]])), use.names = FALSE)
    chrs <- normalize_chr_values(chrs)
    unique(chrs[!is.na(chrs) & chrs != ""])
  }
  
  read_snp_gene_reference <- function(snps, chrs = character(0), ref_dir = "snps_genes_ref") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("The arrow package is required. Install it with install.packages('arrow').")
    }
    if (!dir.exists(ref_dir)) {
      stop(paste0("Reference folder not found: ", ref_dir, ". Place snp_gene_map_chr*.parquet files there."))
    }
    
    snps <- unique(trimws(as.character(snps)))
    snps <- snps[!is.na(snps) & snps != ""]
    if (length(snps) == 0) stop("No valid SNP/rsID identifiers were found in the selected columns.")
    snp_keys <- unique(tolower(snps))
    
    all_files <- list.files(ref_dir, pattern = "^snp_gene_map_chr.*\\.parquet$", full.names = TRUE, ignore.case = TRUE)
    if (length(all_files) == 0) {
      stop(paste0("No files matching snp_gene_map_chr*.parquet were found in ", ref_dir, "."))
    }
    
    # Chromosome label carried by each reference file (e.g. snp_gene_map_chr22.parquet -> "22").
    file_chrs <- normalize_chr_values(
      sub("^snp_gene_map_chr", "", tools::file_path_sans_ext(basename(all_files)), ignore.case = TRUE)
    )
    
    # Chromosomes actually present in the uploaded GWAS.
    chrs <- normalize_chr_values(chrs)
    chrs <- unique(chrs[!is.na(chrs) & chrs != ""])
    
    if (length(chrs) > 0) {
      # STRICT chromosome-wise mapping: read ONLY the parquet files whose chromosome
      # is present in the GWAS. No silent fall-back to the full genome.
      keep         <- file_chrs %in% chrs
      files        <- all_files[keep]
      used_chrs    <- sort(unique(file_chrs[keep]))
      missing_chrs <- sort(setdiff(chrs, file_chrs))  # GWAS chromosomes with no reference file
      
      if (length(files) == 0) {
        stop(paste0(
          "No SNP-gene reference files were found for the chromosome(s) in the input: ",
          paste(sort(chrs), collapse = ", "),
          ". Expected files like snp_gene_map_chr", sort(chrs)[1], ".parquet in ", ref_dir, "."
        ))
      }
    } else {
      # No chromosome column detected (e.g. a plain rsID list): scan every reference file.
      files        <- all_files
      used_chrs    <- sort(unique(file_chrs))
      missing_chrs <- character(0)
    }
    
    mapped_list <- lapply(files, function(f) {
      out <- NULL
      if (requireNamespace("dplyr", quietly = TRUE)) {
        out <- tryCatch({
          ds <- arrow::open_dataset(f, format = "parquet")
          dplyr::collect(
            dplyr::select(
              dplyr::filter(ds, SNP %in% snps),
              Chromosome, SNP, Gene
            )
          )
        }, error = function(e) NULL)
      }
      
      if (is.null(out)) {
        out <- arrow::read_parquet(f, col_select = c("Chromosome", "SNP", "Gene"))
      }
      
      out <- as.data.table(out)
      missing_cols <- setdiff(c("Chromosome", "SNP", "Gene"), names(out))
      if (length(missing_cols) > 0) {
        stop(paste0("Reference file ", basename(f), " is missing columns: ", paste(missing_cols, collapse = ", ")))
      }
      out[, SNP_key := tolower(trimws(as.character(SNP)))]
      out[SNP_key %in% snp_keys, .(Chromosome, SNP, Gene, SNP_key)]
    })
    
    mapped <- rbindlist(mapped_list, fill = TRUE)
    if (nrow(mapped) == 0) {
      mapped <- data.table(Chromosome = character(), SNP = character(), Gene = character(), SNP_key = character())
    } else {
      mapped <- unique(mapped[!is.na(SNP_key) & SNP_key != ""])
    }
    
    # Report which chromosomes / files were used (read by the caller for user feedback).
    attr(mapped, "chr_requested") <- sort(chrs)
    attr(mapped, "chr_used")      <- used_chrs
    attr(mapped, "chr_missing")   <- missing_chrs
    attr(mapped, "files_used")    <- basename(files)
    mapped
  }  
  output$snp_gene_ref_status <- renderUI({
    ref_dir <- "snps_genes_ref"
    if (!dir.exists(ref_dir)) {
      div(
        class = "alert alert-warning mb-0",
        span(fa("triangle-exclamation"), " Reference folder not found: ", code(ref_dir))
      )
    } else {
      n_files <- length(list.files(ref_dir, pattern = "^snp_gene_map_chr.*\\.parquet$", ignore.case = TRUE))
      div(
        class = "alert alert-success mb-0",
        span(fa("folder-open"), " Found ", n_files, " chromosome-wise SNP-gene parquet files in ", code(ref_dir), ".")
      )
    }
  })
  
  observeEvent(input$snp_gene_file, {
    req(input$snp_gene_file)
    updateTextAreaInput(session, "snp_gene_paste_rsids", value = "")
    tryCatch({
      rv$snp_gene_data <- fread(input$snp_gene_file$datapath)
      rv$snp_gene_result <- NULL
      rv$snp_gene_map <- NULL
    }, error = function(e) {
      rv$snp_gene_data <- NULL
      rv$snp_gene_result <- NULL
      rv$snp_gene_map <- NULL
      showNotification(paste("Could not read SNP to Gene input file:", e$message), type = "error")
    })
  })
  
  observeEvent(input$snp_gene_paste_rsids, {
    txt <- input$snp_gene_paste_rsids
    if (!is.null(txt) && grepl("\\w", txt)) {
      tryCatch({
        clean_txt <- gsub("[,; ]+", "\n", txt)
        df <- fread(text = clean_txt, header = FALSE)
        if (ncol(df) == 1) names(df) <- "SNP"
        rv$snp_gene_data <- df
        rv$snp_gene_result <- NULL
        rv$snp_gene_map <- NULL
      }, error = function(e) {
        showNotification(paste("Could not parse pasted rsIDs/SNPs:", e$message), type = "error")
      })
    }
  })
  
  output$snp_gene_cols_ui <- renderUI({
    if (is.null(rv$snp_gene_data)) {
      return(helpText("Upload a file or paste rsIDs to detect SNP identifier columns."))
    }
    col_names <- names(rv$snp_gene_data)
    detected <- detect_snp_identifier_cols(col_names)
    selected <- if (length(detected) > 0) detected else col_names[1]
    tagList(
      checkboxGroupInput(
        "snp_gene_cols",
        "Columns to map",
        choices = col_names,
        selected = selected
      ),
      if (length(detected) == 0) {
        div(class = "alert alert-warning", "No obvious rsID/SNP columns were detected. Please select the correct column manually.")
      } else {
        div(class = "alert alert-info", paste("Auto-detected:", paste(detected, collapse = ", ")))
      }
    )
  })
  output$snp_gene_chr_ui <- renderUI({
    if (is.null(rv$snp_gene_data)) return(NULL)
    col_names <- names(rv$snp_gene_data)
    detected  <- detect_chr_cols(col_names)
    tagList(
      selectInput(
        "snp_gene_chr_col",
        "Chromosome column (restricts parquet files)",
        choices  = c("Auto-detect" = "__AUTO__", "All chromosomes" = "__ALL__", col_names),
        selected = if (length(detected) > 0) "__AUTO__" else "__ALL__"
      ),
      if (length(detected) > 0) {
        div(class = "alert alert-info", style = "font-size: 0.85rem;",
            paste0("Auto-detected chromosome column(s): ", paste(detected, collapse = ", ")))
      } else {
        div(class = "alert alert-warning", style = "font-size: 0.85rem;",
            "No chromosome column auto-detected \u2014 pick one to restrict, or leave 'All chromosomes'.")
      }
    )
  })
  output$snp_gene_preview_section <- renderUI({
    req(rv$snp_gene_data)
    card(
      card_header(span(fa("file-alt"), " SNP to Gene Input Preview")),
      DTOutput("snp_gene_preview_table")
    )
  })
  
  output$snp_gene_preview_table <- renderDT({
    req(rv$snp_gene_data)
    datatable(rv$snp_gene_data, style = "bootstrap5", options = list(pageLength = 5, scrollX = TRUE))
  })
  
  observeEvent(input$run_snp_gene, {
    req(rv$snp_gene_data)
    snp_cols <- input$snp_gene_cols
    if (is.null(snp_cols) || length(snp_cols) == 0) {
      showNotification("Please select at least one rsID/SNP identifier column.", type = "warning")
      return()
    }
    
    withProgress(message = "Mapping SNPs to genes", value = 0, {
      tryCatch({
        incProgress(0.15, "Collecting SNP identifiers...")
        df   <- copy(rv$snp_gene_data)
        snps <- get_unique_snps_from_cols(df, snp_cols)
        chr_choice <- input$snp_gene_chr_col
        if (is.null(chr_choice) || chr_choice == "__AUTO__") {
          chrs <- get_chr_values_from_input(df)
        } else if (chr_choice == "__ALL__") {
          chrs <- character(0)
        } else if (chr_choice %in% names(df)) {
          vals <- normalize_chr_values(as.character(df[[chr_choice]]))
          chrs <- unique(vals[!is.na(vals) & vals != ""])
        } else {
          chrs <- get_chr_values_from_input(df)
        }
        
        if (length(chrs) > 0) {
          incProgress(0, paste0("Restricting to chromosome(s): ", paste(sort(chrs), collapse = ", ")))
        } else {
          incProgress(0, "No chromosome column detected; scanning all reference files.")
        }
        
        incProgress(0.45, "Reading chromosome-wise parquet reference files...")
        # Chromosome-aware lookup: if CHR/CHROM/CHROM_A/CHROM_B columns exist,
        # read only the snp_gene_map_chr*.parquet files for chromosomes present in the input.
        snp_gene_map <- read_snp_gene_reference(snps, chrs = chrs, ref_dir = "snps_genes_ref")
        rv$snp_gene_map <- snp_gene_map
        
        used_chr    <- attr(snp_gene_map, "chr_used")
        missing_chr <- attr(snp_gene_map, "chr_missing")
        if (length(missing_chr) > 0) {
          showNotification(
            paste0("No reference file for chromosome(s): ", paste(missing_chr, collapse = ", "),
                   ". SNPs on these chromosomes were not mapped."),
            type = "warning", duration = 8
          )
        }
        
        incProgress(0.25, "Adding gene columns to the input table...")
        if (nrow(snp_gene_map) > 0) {
          collapsed_map <- snp_gene_map[!is.na(Gene) & Gene != "", .(
            Gene = paste(sort(unique(as.character(Gene))), collapse = ";")
          ), by = SNP_key]
        } else {
          collapsed_map <- data.table(SNP_key = character(), Gene = character())
        }
        
        for (cc in snp_cols) {
          key <- tolower(trimws(as.character(df[[cc]])))
          gene_values <- collapsed_map$Gene[match(key, collapsed_map$SNP_key)]
          df[, (paste0(cc, "_Gene")) := gene_values]
        }
        
        # Place each generated Gene column immediately after its corresponding SNP/rsID column.
        current_cols <- names(df)
        ordered_cols <- character(0)
        for (cc in current_cols) {
          ordered_cols <- c(ordered_cols, cc)
          gene_cc <- paste0(cc, "_Gene")
          if (cc %in% snp_cols && gene_cc %in% current_cols) ordered_cols <- c(ordered_cols, gene_cc)
        }
        ordered_cols <- unique(c(ordered_cols, setdiff(current_cols, ordered_cols)))
        setcolorder(df, ordered_cols)
        
        rv$snp_gene_result <- df
        rv$snp_gene_selected_cols <- snp_cols
        
        incProgress(0.15, "Rendering results...")
        showNotification(
          paste0("SNP to Gene mapping completed",
                 if (length(used_chr) > 0) paste0(" for chromosome(s): ", paste(used_chr, collapse = ", ")) else "",
                 "."),
          type = "message"
        )
      }, error = function(e) {
        rv$snp_gene_result <- NULL
        showNotification(paste("SNP to Gene mapping failed:", e$message), type = "error")
      })
    })
  })
  
  output$val_snp_gene_input <- renderUI({
    if (is.null(rv$snp_gene_data)) return("0")
    snp_cols <- input$snp_gene_cols
    if (is.null(snp_cols) || length(snp_cols) == 0) snp_cols <- detect_snp_identifier_cols(names(rv$snp_gene_data))
    if (length(snp_cols) == 0) return("0")
    formatC(length(get_unique_snps_from_cols(rv$snp_gene_data, snp_cols)), big.mark = ",")
  })
  
  output$val_snp_gene_mapped <- renderUI({
    if (is.null(rv$snp_gene_map) || nrow(rv$snp_gene_map) == 0) return("0")
    mapped <- unique(rv$snp_gene_map[!is.na(Gene) & Gene != "", SNP_key])
    formatC(length(mapped), big.mark = ",")
  })
  
  output$val_snp_gene_genes <- renderUI({
    if (is.null(rv$snp_gene_map) || nrow(rv$snp_gene_map) == 0) return("0")
    genes <- unique(rv$snp_gene_map[!is.na(Gene) & Gene != "", as.character(Gene)])
    formatC(length(genes), big.mark = ",")
  })
  
  output$snp_gene_result_table <- renderDT({
    req(rv$snp_gene_result)
    datatable(rv$snp_gene_result, extensions = 'Buttons', style = 'bootstrap5',
              options = list(pageLength = 15, scrollX = TRUE, dom = 'Bfrtip',
                             buttons = list('copy', 'csv', 'excel')))
  }, server = FALSE)
  
  output$download_snp_gene <- downloadHandler(
    filename = function() { paste0("SNP_to_Gene_results_", format(Sys.Date(), "%Y%m%d"), ".txt") },
    content = function(file) {
      req(rv$snp_gene_result)
      fwrite(rv$snp_gene_result, file, sep = "\t", quote = FALSE, na = "NA")
    }
  )
  
  # --- DOCUMENTATION DATA & TABLES ---
  
  # Centralized Data for the reference-panel table and documentation plots
  ref_data_raw <- reactive({
    df <- data.frame(
      Panel_Raw = c(
        "TOP-LD (hg38)", "", "", "",
        "PhenoScanner (Phase 3, hg19 & hg38)", "", "", "", "",
        "HapMap (hg16)", "", "", "", "", "", "", "", "", "", "",
        "1000 Genomes Project Phase 3 (hg38)", "", "", "", "",
        "1000 Genomes Project Phase 3 (hg38), high coverage", "", "", "", "",
        "UKBB (hg19)", "", "", "", "", "",
        "LASI-DAD (hg38)",
        "HGDP", "", "", "", "", "", ""
      ),
      Label = c(
        "EUR", "AFR", "SAS", "EAS",
        "EUR", "AFR", "SAS", "EAS", "AMR",
        "ASW", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "TSI", "YRI",
        "EUR", "AFR", "SAS", "EAS", "AMR",
        "EUR", "AFR", "SAS", "EAS", "AMR",
        "EUR", "AFR", "CSA", "EAS", "AMR", "MID",
        "IND",
        "EUR", "AFR", "CSA", "EAS", "AMR", "MID", "OCN"
      ),
      Sample_Name = c(
        "European", "African", "South Asian", "East Asian",
        "European", "African", "South Asian", "East Asian", "American",
        "African ancestry in Southwest USA",
        "Utah residents with Northern and Western European ancestry from the CEPH collection",
        "Han Chinese in Beijing, China",
        "Chinese in Metropolitan Denver, Colorado",
        "Gujarati Indians in Houston, Texas",
        "Japanese in Tokyo, Japan",
        "Luhya in Webuye, Kenya",
        "Mexican ancestry in Los Angeles, California",
        "Maasai in Kinyawa, Kenya",
        "Toscans in Italy",
        "Yoruba in Ibadan, Nigeria",
        "European", "African", "South Asian", "East Asian", "Admixed American",
        "European", "African", "South Asian", "East Asian", "Admixed American",
        "European", "African", "Central and South Asian", "East Asian", "Admixed American", "Middle Eastern",
        "Indian",
        "European", "African", "Central South Asian", "East Asian", "Admixed American", "Middle Eastern", "Oceanian"
      ),
      Count_Str = c(
        "13,160", "1,335", "239", "844",
        "503", "661", "489", "504", "347",
        "90", "180", "90", "100", "100", "91", "100", "90", "180", "100", "180",
        "632", "893", "601", "585", "490",
        "632", "893", "601", "585", "490",
        "362,446", "6,255", "8,284", "2,700", "987", "1,567",
        "2,680",
        "155", "104", "197", "223", "61", "161", "28"
      ),
      Variants_Str = c(
        "69,524,944", "60,392,677", "22,309,649", "35,538,656",
        "11,159,862", "11,159,862", "11,159,862", "11,159,862", "11,159,862",
        "1,561,113", "1,412,161", "1,328,283", "1,305,880", "1,407,540", "1,296,969", "1,529,438", "1,409,947", "1,419,626", "1,419,920", "1,501,085",
        "8,193,280", "13,876,891", "8,579,150", "7,245,426", "9,347,814",
        "9,263,406", "15,952,015", "9,439,730", "7,989,898", "10,106,451",
        "1,431,634", "1,259,175", "1,379,991", "1,146,910", "1,286,943", "1,307,269",
        "119,948,583",
        "10,144,721", "18,062,644", "11,281,764", "9,631,805", "7,592,946", "11,731,317", "8,147,809"
      ),
      stringsAsFactors = FALSE
    )
    
    # Fill down panel names for grouping while keeping Panel_Raw blank rows for display.
    filled_panel <- as.character(df$Panel_Raw)
    for (i in seq_along(filled_panel)) {
      if (i > 1 && (filled_panel[i] == "" || is.na(filled_panel[i]))) {
        filled_panel[i] <- filled_panel[i - 1]
      }
    }
    df$Panel <- filled_panel
    df$Count <- as.numeric(gsub(",", "", df$Count_Str))
    df$Variants <- as.numeric(gsub(",", "", df$Variants_Str))
    df
  })
  
  # Render Table
  output$ref_pop_table <- renderTable({
    df <- ref_data_raw()
    # Format for display (keep original columns but rename nicely)
    disp_df <- df[, c("Panel_Raw", "Label", "Sample_Name", "Count_Str", "Variants_Str")]
    names(disp_df) <- c("Reference Panel", "Label", "Population Sample", "Number of Samples", "Number of Variants")
    disp_df
  }, striped = TRUE, hover = TRUE, bordered = TRUE, width = "100%")
  
  # Chart 1: Total Samples per Panel
  output$doc_samples_plot <- renderPlotly({
    df <- ref_data_raw()
    
    # Aggregation
    summ_df <- aggregate(Count ~ Panel, data = df, sum)
    summ_df <- summ_df[order(summ_df$Count, decreasing = TRUE), ]
    summ_df$Panel <- factor(summ_df$Panel, levels = summ_df$Panel) # Lock order
    
    plot_ly(summ_df, x = ~Panel, y = ~Count, type = 'bar', 
            marker = list(color = '#00d1b2', line = list(color = '#2c3e50', width = 1))) %>%
      layout(
        title = "Total Samples per Reference Panel",
        xaxis = list(title = "", tickangle = -15),
        yaxis = list(title = "Total Samples"),
        font = list(family = "Inter"),
        margin = list(b = 60)
      )
  })
  
  # Chart 2: Unique Populations per Panel
  output$doc_pops_count_plot <- renderPlotly({
    df <- ref_data_raw()
    
    # Count rows per panel
    count_df <- aggregate(Label ~ Panel, data = df, length)
    count_df <- count_df[order(count_df$Label, decreasing = TRUE), ]
    count_df$Panel <- factor(count_df$Panel, levels = count_df$Panel)
    
    plot_ly(count_df, x = ~Panel, y = ~Label, type = 'bar',
            marker = list(color = '#2c3e50', line = list(color = '#00d1b2', width = 1))) %>%
      layout(
        title = "Sub-Populations per Reference Panel",
        xaxis = list(title = "", tickangle = -15),
        yaxis = list(title = "Number of Populations"),
        font = list(family = "Inter"),
        margin = list(b = 60)
      )
  })
  
  # Chart 3: Unique rsIDs / variants per reference-panel population
  output$doc_variants_plot <- renderPlotly({
    df <- ref_data_raw()
    df$Panel_Label <- paste(df$Panel, df$Label, sep = " / ")
    df <- df[order(df$Variants, decreasing = TRUE), ]
    df$Panel_Label <- factor(df$Panel_Label, levels = df$Panel_Label)
    
    plot_ly(
      df,
      x = ~Panel_Label,
      y = ~Variants,
      type = 'bar',
      text = ~Variants_Str,
      hovertemplate = paste(
        '<b>%{x}</b><br>',
        'Unique rsIDs / variants: %{text}<extra></extra>'
      ),
      marker = list(color = '#ffb703', line = list(color = '#2c3e50', width = 1))
    ) %>%
      layout(
        title = "Unique rsIDs / Variants by Reference Panel",
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "Number of Variants", type = "log"),
        font = list(family = "Inter"),
        margin = list(b = 150)
      )
  })
}

shinyApp(ui, server)