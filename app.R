library(shiny)
library(shinydashboard)
library(httr)
library(jsonlite)
library(data.table)
library(plotly)
library(DT) 

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

analyze_drug_safety <- function(drug_name, sample_size = 500) {
  query <- paste0('patient.drug.medicinalproduct:', drug_name)
  data_res <- get_fda_data(query, limit = sample_size)
  
  if (is.null(data_res)) stop("No data found for: ", drug_name)
  
  dt_ev <- as.data.table(data_res$results)
  cols <- names(dt_ev)
  safe_extract <- function(col_name) if(col_name %in% cols) dt_ev[[col_name]] else rep(NA, nrow(dt_ev))
  
  dt_demo <- data.table(
    id = dt_ev$safetyreportid,
    hosp = !is.na(safe_extract("seriousnesshospitalization")) & safe_extract("seriousnesshospitalization") == "1",
    sex = fcase(
      safe_extract("patient.patientsex") == "1", "Male",
      safe_extract("patient.patientsex") == "2", "Female",
      default = "Unknown"
    ),
    age = as.numeric(safe_extract("patient.patientonsetage")),
    unit = safe_extract("patient.patientonsetageunit")
  )
  
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
  
  hosp_summary <- dt_full[hosp == TRUE, .(cases = .N), by = .(sex, age_grp)]
  top_causes <- dt_full[hosp == TRUE, .(cases = .N), by = reactionmeddrapt][order(-cases)][1:10]
  
  return(list(demo = hosp_summary, top = top_causes))
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
  dashboardHeader(title = "OpenFDA Pharmacovigilance"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Panel A: Drug Safety", tabName = "panel_a", icon = icon("pills")),
      menuItem("Panel B: Syndrome Analysis", tabName = "panel_b", icon = icon("heartbeat"))
    )
  ),
  
  dashboardBody(
    # Custom CSS for validation errors
    tags$head(tags$style(HTML(".shiny-output-error-validation { color: #d9534f; font-weight: bold; }"))),
    
    tabItems(
      # PANEL A TAB
      tabItem(tabName = "panel_a",
              fluidRow(
                box(width = 4, title = "Search Parameter", status = "primary", solidHeader = TRUE,
                    textInput("drug_input", "Enter Generic Drug Name:", value = "IBUPROFEN"),
                    actionButton("run_drug", "Analyze Drug", class = "btn-primary btn-block")
                ),
                box(width = 8, title = "Hospitalization Demographics", status = "info", solidHeader = TRUE,
                    plotlyOutput("demo_plot")
                )
              ),
              fluidRow(
                box(width = 12, title = "Top 10 Causes of Hospitalization", status = "warning", solidHeader = TRUE,
                    DTOutput("top_rx_table")
                )
              )
      ),
      
      # PANEL B TAB
      tabItem(tabName = "panel_b",
              fluidRow(
                box(width = 4, title = "Search Parameter", status = "danger", solidHeader = TRUE,
                    textInput("syndrome_input", "Enter Syndrome/Reaction:", value = "STEVENS-JOHNSON SYNDROME"),
                    actionButton("run_syndrome", "Identify Culprits", class = "btn-danger btn-block"),
                    hr(),
                    helpText("Hint: Use British English (e.g., Diarrhoea)."),
                    uiOutput("meddra_link")
                ),
                box(width = 8, title = "Reporting Odds Ratio (ROR) Signals", status = "danger", solidHeader = TRUE,
                    DTOutput("ror_table")
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
  
  # External link UI
  output$meddra_link <- renderUI({
    tagList("Check terms at:", a("MedDRA Web Search", href="https://www.meddra.org/how-to-use/tools/web-based-browser", target="_blank"))
  })
  
  # Reactive execution for Panel A
  drug_data <- eventReactive(input$run_drug, {
    withProgress(message = 'Fetching FDA Data...', value = 0.5, {
      analyze_drug_safety(toupper(trimws(input$drug_input)))
    })
  })
  
  # Render Demographics Plot (Plotly)
  output$demo_plot <- renderPlotly({
    req(drug_data())
    dt <- drug_data()$demo
    
    dt_clean <- dt[sex != "Unknown" & age_grp != "Unknown"]
    
    plot_ly(dt_clean, x = ~age_grp, y = ~cases, color = ~sex, type = 'bar', barmode = 'group') %>%
      layout(title = "Hospitalized Patients by Age and Sex",
             xaxis = list(title = "Age Group"),
             yaxis = list(title = "Number of Cases"))
  })
  
  # Render Top Reactions Table (DT)
  output$top_rx_table <- renderDT({
    req(drug_data())
    datatable(drug_data()$top, 
              options = list(pageLength = 5, dom = 't'),
              colnames = c("Reaction (MedDRA PT)", "Local Sample Cases"))
  })
  
  # Reactive execution for Panel B
  syndrome_data <- eventReactive(input$run_syndrome, {
    withProgress(message = 'Calculating ROR Signals...', value = 0.5, {
      analyze_syndrome_panel_b(toupper(trimws(input$syndrome_input)))
    })
  })
  
  # Render Culprit ROR Table (DT)
  output$ror_table <- renderDT({
    req(syndrome_data())
    dt <- syndrome_data()
    
    datatable(dt,
              rownames = FALSE,
              options = list(pageLength = 10, dom = 'tp'),
              colnames = c("Suspect Drug", "Global Cases", "ROR", "Lower 95% CI", "Upper 95% CI")) %>%
      formatStyle('ror', 
                  backgroundColor = styleInterval(c(2, 5), c('white', 'lightyellow', '#ffcccc')))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
