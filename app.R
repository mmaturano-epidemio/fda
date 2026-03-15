library(shiny)
library(shinydashboard)
library(httr)
library(jsonlite)
library(data.table)
library(plotly)
library(DT) 
library(shinyjs)

# ==============================================================================
# 1. CORE API ENGINE & UTILITIES
# ==============================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b

get_fda_data <- function(query_string, limit = 100, count_field = NULL) {
  base_url <- "https://api.fda.gov/drug/event.json"
  params <- list(search = query_string)
  
  if (!is.null(count_field)) {
    params$count <- count_field
  } else {
    params$limit <- limit
  }
  
  resp <- GET(base_url, query = params)
  if (status_code(resp) != 200) return(NULL)
  
  res <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
  return(res)
}

calculate_ror <- function(drug, reaction) {
  q_a <- paste0('patient.drug.medicinalproduct.exact:"', drug, '" AND patient.reaction.reactionmeddrapt.exact:"', reaction, '"')
  res_a <- get_fda_data(q_a, limit = 1)
  
  a <- as.numeric(res_a$meta$results$total %||% 0)
  if (a == 0) return(NULL)
  
  n_drug  <- as.numeric(get_fda_data(paste0('patient.drug.medicinalproduct.exact:"', drug, '"'), limit = 1)$meta$results$total %||% 0)
  n_event <- as.numeric(get_fda_data(paste0('patient.reaction.reactionmeddrapt.exact:"', reaction, '"'), limit = 1)$meta$results$total %||% 0)
  n_total <- as.numeric(get_fda_data('_exists_:safetyreportid', limit = 1)$meta$results$total %||% 0)
  
  b <- n_drug - a
  c <- n_event - a
  d <- n_total - (a + b + c)
  
  ror_val <- (a * d) / (b * c)
  se_ln   <- sqrt(1/a + 1/b + 1/c + 1/d)
  
  return(data.table(
    term = drug, 
    cases = a, 
    ror = round(ror_val, 2),
    ci_lower = round(exp(log(ror_val) - 1.96 * se_ln), 2),
    ci_upper = round(exp(log(ror_val) + 1.96 * se_ln), 2)
  ))
}

# ==============================================================================
# 2. ANALYTICAL FUNCTIONS (Panel A & B)
# ==============================================================================

analyze_drug_safety <- function(drug_name, serious_filter, sample_size = 500) {
  query <- paste0('patient.drug.medicinalproduct:', drug_name)
  data_res <- get_fda_data(query, limit = sample_size)
  
  if (is.null(data_res)) stop("No data found for: ", drug_name)
  
  dt_ev <- as.data.table(data_res$results)
  cols <- names(dt_ev)
  safe_extract <- function(col_name) if(col_name %in% cols) dt_ev[[col_name]] else rep(NA, nrow(dt_ev))
  
  # Extraer todos los flags de seriedad
  dt_demo <- data.table(
    id = dt_ev$safetyreportid,
    hosp = !is.na(safe_extract("seriousnesshospitalization")) & safe_extract("seriousnesshospitalization") == "1",
    death = !is.na(safe_extract("seriousnessdeath")) & safe_extract("seriousnessdeath") == "1",
    lifethreat = !is.na(safe_extract("seriousnesslifethreatening")) & safe_extract("seriousnesslifethreatening") == "1",
    disabling = !is.na(safe_extract("seriousnessdisabling")) & safe_extract("seriousnessdisabling") == "1",
    congenital = !is.na(safe_extract("seriousnesscongenitalanomaly")) & safe_extract("seriousnesscongenitalanomaly") == "1",
    other = !is.na(safe_extract("seriousnessother")) & safe_extract("seriousnessother") == "1",
    sex = fcase(
      safe_extract("patient.patientsex") == "1", "Male",
      safe_extract("patient.patientsex") == "2", "Female",
      default = "Unknown"
    ),
    age = as.numeric(safe_extract("patient.patientonsetage")),
    unit = safe_extract("patient.patientonsetageunit")
  )
  
  # Create serious_any column
  dt_demo[, serious_any := hosp | death | lifethreat | disabling | congenital | other]
  
  dt_demo[, age_y := fcase(unit=="800", age*10, unit=="801", age, unit=="802", age/12, unit=="804", age/365, default=NA_real_)]
  dt_demo[, age_grp := fcase(age_y < 18, "Pediatric", age_y < 65, "Adult", age_y >= 65, "Elderly", default="Unknown")]
  
  dt_rx <- rbindlist(lapply(1:nrow(dt_ev), function(i) {
    rx_list <- dt_ev$patient.reaction[[i]]
    if (is.null(rx_list)) return(NULL)
    res <- as.data.table(rx_list)
    res[, id := dt_ev$safetyreportid[i]]
    return(res[, .(id, reactionmeddrapt)])
  }))
  
  dt_full <- merge(dt_demo, dt_rx, by = "id")
  
  # Aplicar filtro de seriedad
  if (serious_filter == "All reports") {
    dt_filtered <- dt_full
  } else if (serious_filter == "Serious (any)") {
    dt_filtered <- dt_full[serious_any == TRUE]
  } else if (serious_filter == "Hospitalization") {
    dt_filtered <- dt_full[hosp == TRUE]
  } else if (serious_filter == "Death") {
    dt_filtered <- dt_full[death == TRUE]
  } else if (serious_filter == "Life threatening") {
    dt_filtered <- dt_full[lifethreat == TRUE]
  } else if (serious_filter == "Disability") {
    dt_filtered <- dt_full[disabling == TRUE]
  } else if (serious_filter == "Congenital anomaly") {
    dt_filtered <- dt_full[congenital == TRUE]
  } else if (serious_filter == "Other serious") {
    dt_filtered <- dt_full[other == TRUE]
  } else {
    dt_filtered <- dt_full  # fallback
  }
  
  # Filter based summaries 
  hosp_summary <- dt_filtered[, .(cases = .N), by = .(sex, age_grp)]
  top_causes <- dt_filtered[, .(cases = .N), by = reactionmeddrapt][order(-cases)][1:10]
  
  return(list(
    demo = hosp_summary, 
    top = top_causes,
    filter_name = serious_filter,
    drug_name = drug_name
  ))
}

analyze_syndrome_panel_b <- function(syndrome_name, top_n = 10) {
  query <- paste0('patient.reaction.reactionmeddrapt.exact:"', syndrome_name, '"')
  data_res <- get_fda_data(query, count_field = "patient.drug.medicinalproduct.exact")
  
  if (is.null(data_res)) stop("Syndrome not found.")
  
  top_drugs <- as.data.table(data_res$results)[1:(top_n + 5)]
  top_drugs[, term := gsub("\\.", "", term)]
  top_drugs <- top_drugs[, .(count = sum(count)), by = term][order(-count)][1:top_n]
  
  results <- rbindlist(lapply(top_drugs$term, function(d) {
    Sys.sleep(0.1) 
    calculate_ror(d, syndrome_name)
  }))
  
  return(results[order(-ror)])
}

# ==============================================================================
# 3. SHINY USER INTERFACE (UI)
# ==============================================================================

ui <- dashboardPage(
  dashboardHeader(title = "openFDA Explorer"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Welcome / Instructions", tabName = "welcome", icon = icon("info-circle")),
      menuItem("Panel A: Drug Safety", tabName = "panel_a", icon = icon("pills")),
      menuItem("Panel B: Syndrome Analysis", tabName = "panel_b", icon = icon("heartbeat")),
      menuItem("About / Contact", tabName = "about", icon = icon("user-md"))
    )
  ),
  
  dashboardBody(
    tags$head(tags$style(HTML("
      .shiny-output-error-validation { color: #d9534f; font-weight: bold; }
      .progress-bar { background-color: #337ab7; }
    "))),
    
    tabItems(
      # WELCOME / INSTRUCTIONS TAB
      tabItem(tabName = "welcome",
              fluidRow(
                box(width = 12, title = "About this application", status = "primary", solidHeader = TRUE,
                    p("This tool allows exploration of pharmacovigilance data from the openFDA database."),
                    p("It provides two complementary analyses:"),
                    tags$ul(
                      tags$li(strong("Panel A (Drug Safety):"), " Enter a drug name to see demographic characteristics and top adverse reactions. You can filter by seriousness criteria (e.g., hospitalization, death, life-threatening)."),
                      tags$li(strong("Panel B (Syndrome Analysis):"), " Enter an adverse reaction (syndrome) to calculate the Reporting Odds Ratio (ROR) for the most frequently reported drugs.")
                    ),
                    br(),
                    
                    # ROR explanation 
                    h4("Understanding the Reporting Odds Ratio (ROR)"),
                    p("The ROR is a measure of disproportionate reporting used in pharmacovigilance to detect potential safety signals. It compares the odds of exposure to a specific drug among reports of a specific reaction versus all other reactions."),
                    p(strong("Interpretation:")),
                    tags$ul(
                      tags$li("ROR = 1: The reaction is reported as often with the drug as with all other drugs (no signal)."),
                      tags$li("ROR > 1: The reaction is reported more frequently with the drug (potential signal)."),
                      tags$li("ROR > 2 is often considered a signal of interest, especially if the lower bound of the 95% confidence interval exceeds 1.")
                    ),
                    p(strong("Why does it take time?")),
                    p("For each candidate drug, the application makes four separate queries to the FDA API (to obtain the numbers needed for the 2x2 table). If a syndrome is associated with many drugs (e.g., common reactions like diarrhoea), the calculation may take up to a minute. The progress bar shows the current step."),
                    br(),
                    
                    h4("Understanding the data: Sample vs. Full database"),
                    p(strong("Panel A uses a sample of 500 reports.")),
                    tags$ul(
                      tags$li("The numbers shown are counts within this sample."),
                      tags$li("The seriousness filter applies to the sample, so absolute numbers reflect the filtered subset.")
                    ),
                    p(strong("Panel B uses global counts but limits to top 10 drugs.")),
                    br(),
                    
                    h4("Performance notes"),
                    p("Typical waiting time: 15-30 seconds. Progress bar shows current step."),
                    br(),
                    
                    h4("Example search terms"),
                    tags$ul(
                      tags$li(strong("Panel A:"), " ROSIGLITAZONE, IBUPROFEN, ATORVASTATIN"),
                      tags$li(strong("Panel B:"), " STEVENS-JOHNSON SYNDROME, HEPATIC FAILURE, CARDIAC ARREST")
                    ),
                    p(a("MedDRA Web Search", href="https://www.meddra.org/how-to-use/tools/web-based-browser", target="_blank"))
                )
              )
      ),
      
      # PANEL A TAB 
      tabItem(tabName = "panel_a",
              fluidRow(
                box(width = 4, title = "Search parameters", status = "primary", solidHeader = TRUE,
                    textInput("drug_input", "Generic drug name:", 
                              value = "ROSIGLITAZONE", placeholder = "e.g., ROSIGLITAZONE"),
                    selectInput("serious_filter", "Seriousness filter:",
                                choices = c("All reports", 
                                            "Serious (any)", 
                                            "Hospitalization", 
                                            "Death", 
                                            "Life threatening", 
                                            "Disability", 
                                            "Congenital anomaly", 
                                            "Other serious"),
                                selected = "Hospitalization"),
                    actionButton("run_drug", "Analyze drug", class = "btn-primary btn-block"),
                    br(),
                    helpText("Examples: ROSIGLITAZONE, IBUPROFEN, ATORVASTATIN"),
                    helpText(strong("Note:"), "Analysis based on 500 most recent reports. Filter applied after download.")
                ),
                box(width = 8, title = "Hospitalization demographics", status = "info", solidHeader = TRUE,
                    plotlyOutput("demo_plot")
                )
              ),
              fluidRow(
                box(width = 12, title = "Top 10 causes (based on filtered sample)", status = "warning", solidHeader = TRUE,
                    DTOutput("top_rx_table")
                )
              )
      ),
      
      # PANEL B TAB 
      tabItem(tabName = "panel_b",
              fluidRow(
                box(width = 4, title = "Search parameters", status = "danger", solidHeader = TRUE,
                    textInput("syndrome_input", "Syndrome / adverse reaction:", 
                              value = "STEVENS-JOHNSON SYNDROME", placeholder = "e.g., HEPATIC FAILURE"),
                    actionButton("run_syndrome", "Identify culprit drugs", class = "btn-danger btn-block"),
                    hr(),
                    helpText("Use British English (e.g., DIARRHOEA, not DIARRHEA)."),
                    helpText("Examples: STEVENS-JOHNSON SYNDROME, HEPATIC FAILURE, CARDIAC ARREST"),
                    helpText(strong("Note:"), "Analysis uses global counts but limits to top 10 drugs."),
                    uiOutput("meddra_link")
                ),
                box(width = 8, title = "Reporting Odds Ratio (ROR) signals", status = "danger", solidHeader = TRUE,
                    DTOutput("ror_table")
                )
              )
      ),
      
      # ABOUT / CONTACT TAB
      tabItem(tabName = "about",
              fluidRow(
                box(width = 12, title = "About this application", status = "primary", solidHeader = TRUE,
                    h4("OpenFDA Pharmacovigilance Explorer"),
                    p("This interactive dashboard was developed to facilitate the exploration and analysis of adverse event reports from the FDA's openFDA platform. It is intended for researchers, regulators, and the general public interested in drug safety signals."),
                    br(),
                    h4("Data source"),
                    p("All data are obtained in real-time from the ", 
                      a("openFDA Drug Adverse Event endpoint", href="https://open.fda.gov/apis/drug/event/", target="_blank"), 
                      ". The application queries the API directly; no local copy of the data is stored."),
                    p(em("Note: This dashboard is an independent visualization project and is not affiliated with or endorsed by the U.S. Food and Drug Administration (FDA).")),
                    p(em("Results are intended for exploratory research purposes only and do not constitute medical or regulatory advice. Disproportionality signals do not imply causality.")),
                    br(),
                    h4("About the author"),
                    p("Manuel Maturano is a physician and epidemiologist specialized in data analysis and clinical evidence evaluation for public health decision-making."),
                    p("He currently works in the Health Technology Assessment Committee of the Neuquén Province Ministry of Health, Argentina, at the intersection of clinical epidemiology, health economics, and data science."),
                    p("Contact: ", a("mmaturano.epidemio@gmail.com", href="mailto:mmaturano.epidemio@gmail.com")),
                    br(),
                    h4("Technical notes"),
                    p("This project was built entirely in R, prioritizing computational efficiency through the use of the ", code("data.table"), " package. The interactive interface uses ", code("shinydashboard"), ", with visualizations created using ", code("plotly"), " and ", code("DT"), ". Progress indicators and button locking are implemented with ", code("shinyjs"), "."),
                    p("The application makes multiple calls to the openFDA API; response times depend on the complexity of the query and the API's own performance."),
                    br(),
                    h4("License and attribution"),
                    p("This application and its source code are made available under the ", 
                      a("MIT License", href="https://opensource.org/licenses/MIT", target="_blank"), 
                      ". You are free to use, modify, and distribute it, provided that you retain the following attribution:"),
                    p(strong("If you use this code or a modified version of it in your own work, please include the following credit:")),
                    p("\"Based on the openFDA Pharmacovigilance Explorer by Manuel Maturano (mmaturano.epidemio@gmail.com). \n Copyright (c) 2025 Manuel Maturano\""),
                    p("The original source code is available at: ", 
                      a("github.com/mmaturano-epidemio/fda", href="https://github.com/mmaturano-epidemio/fda", target="_blank")),
                    h4("Version"),
                    p("1.0.0 - March 2025")
                )
              )
      )
    )
  )
)

# ==============================================================================
# 4. SHINY SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {
  
  output$meddra_link <- renderUI({
    tagList("Check terms at:", a("MedDRA Web Search", href="https://www.meddra.org/how-to-use/tools/web-based-browser", target="_blank"))
  })
  
  # Reactive execution for Panel A with button disabling
  drug_data <- eventReactive(input$run_drug, {
    shinyjs::disable("run_drug")
    
    withProgress(message = 'Processing request...', value = 0, {
      incProgress(0.2, detail = "Contacting FDA API...")
      Sys.sleep(0.5)
      incProgress(0.4, detail = "Downloading drug safety data...")
      result <- analyze_drug_safety(
        drug_name = toupper(trimws(input$drug_input)),
        serious_filter = input$serious_filter
      )
      incProgress(0.8, detail = "Processing demographics and reactions...")
      Sys.sleep(0.5)
      incProgress(1, detail = "Done")
    })
    
    shinyjs::enable("run_drug")
    result
  })
  
  # Render Demographics Plot (Plotly)
  output$demo_plot <- renderPlotly({
    req(drug_data())
    dt <- drug_data()$demo
    
    dt_clean <- dt[sex != "Unknown" & age_grp != "Unknown"]
    
    plot_ly(dt_clean, x = ~age_grp, y = ~cases, color = ~sex, type = 'bar', barmode = 'group') %>%
      layout(title = paste0("Patients by age and sex (", drug_data()$filter_name, ")"),
             xaxis = list(title = "Age group"),
             yaxis = list(title = "Number of cases in sample"))
  })
  
  # Render Top Reactions Table (DT)
  output$top_rx_table <- renderDT({
    req(drug_data())
    datatable(drug_data()$top, 
              options = list(pageLength = 10, dom = 't'),
              colnames = c("Reaction (MedDRA PT)", "Cases in sample"))
  })
  
  # Panel B
  syndrome_data <- eventReactive(input$run_syndrome, {
    shinyjs::disable("run_syndrome")
    
    withProgress(message = 'Processing request...', value = 0, {
      incProgress(0.1, detail = "Contacting FDA API...")
      Sys.sleep(0.5)
      incProgress(0.3, detail = "Identifying top drugs for this syndrome...")
      incProgress(0.6, detail = "Calculating ROR for each drug...")
      result <- analyze_syndrome_panel_b(toupper(trimws(input$syndrome_input)))
      incProgress(0.9, detail = "Formatting results...")
      Sys.sleep(0.5)
      incProgress(1, detail = "Done")
    })
    
    shinyjs::enable("run_syndrome")
    result
  })
  
  output$ror_table <- renderDT({
    req(syndrome_data())
    dt <- syndrome_data()
    datatable(dt,
              rownames = FALSE,
              options = list(pageLength = 10, dom = 'tp'),
              colnames = c("Suspect drug", "Global cases", "ROR", "Lower 95% CI", "Upper 95% CI")) %>%
      formatStyle('ror', 
                  backgroundColor = styleInterval(c(2, 5), c('white', 'lightyellow', '#ffcccc')))
  })
}

options(shiny.error = NULL)

shinyApp(ui = ui, server = server)
