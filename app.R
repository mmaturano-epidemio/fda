library(shiny)
library(shinydashboard)
library(httr)
library(jsonlite)
library(data.table)
library(plotly)
library(DT)

# ==============================================================================
# 1. CORE ENGINE
# ==============================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b

get_fda_data <- function(query_string, limit = 100, count_field = NULL) {
  base_url <- "https://api.fda.gov/drug/event.json"
  params <- list(search = query_string)
  if (!is.null(count_field)) params$count <- count_field else params$limit <- limit
  
  resp <- GET(base_url, query = params)
  if (status_code(resp) != 200) return(NULL)
  fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
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
  
  data.table(term = drug, cases = a, ror = round(ror_val, 2),
             ci_lower = round(exp(log(ror_val) - 1.96 * se_ln), 2),
             ci_upper = round(exp(log(ror_val) + 1.96 * se_ln), 2))
}

# ==============================================================================
# 2. ANALYTICS
# ==============================================================================

analyze_drug_safety <- function(drug_name, sample_size = 500) {
  query <- paste0('patient.drug.medicinalproduct:', drug_name)
  data_res <- get_fda_data(query, limit = sample_size)
  if (is.null(data_res)) return(NULL)
  
  dt_ev <- as.data.table(data_res$results)
  cols  <- names(dt_ev)
  safe_extract <- function(col_name) if(col_name %in% cols) dt_ev[[col_name]] else rep(NA, nrow(dt_ev))
  
  dt_demo <- data.table(
    id = dt_ev$safetyreportid,
    hosp = !is.na(safe_extract("seriousnesshospitalization")) & safe_extract("seriousnesshospitalization") == "1",
    sex = fcase(safe_extract("patient.patientsex") == "1", "Male", safe_extract("patient.patientsex") == "2", "Female", default = "Unknown"),
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
    res[, .(id, reactionmeddrapt)]
  }))
  
  dt_full <- merge(dt_demo, dt_rx, by = "id")
  hosp_summary <- dt_full[hosp == TRUE, .(cases = .N), by = .(sex, age_grp)]
  top_causes   <- dt_full[hosp == TRUE, .(local_cases = .N), by = reactionmeddrapt][order(-local_cases)][1:10]
  
  signals <- rbindlist(lapply(top_causes$reactionmeddrapt, function(rx) {
    Sys.sleep(0.1) 
    res <- calculate_ror(drug_name, rx)
    if(!is.null(res)) res[, reactionmeddrapt := rx]
    res
  }))
  
  if(nrow(signals) > 0) {
    top_final <- merge(top_causes, signals, by = "reactionmeddrapt")
    top_final <- top_final[, .(Reaction = reactionmeddrapt, Local_N = local_cases, Global_N = cases, ROR = ror, `95% CI` = paste0("[", ci_lower, "-", ci_upper, "]"))]
    top_final <- top_final[order(-Local_N)]
  } else {
    top_final <- top_causes
  }
  
  list(demo = hosp_summary, top = top_final)
}

analyze_syndrome_panel_b <- function(syndrome_name, top_n = 10) {
  query <- paste0('patient.reaction.reactionmeddrapt.exact:"', syndrome_name, '"')
  data_res <- get_fda_data(query, count_field = "patient.drug.medicinalproduct.exact")
  if (is.null(data_res)) return(NULL)
  
  top_drugs <- as.data.table(data_res$results)[1:(min(nrow(data_res$results), top_n + 5))]
  top_drugs[, term := gsub("\\.", "", term)]
  top_drugs <- top_drugs[, .(count = sum(count)), by = term][order(-count)][1:min(.N, top_n)]
  
  results <- rbindlist(lapply(top_drugs$term, function(d) {
    Sys.sleep(0.1)
    calculate_ror(d, syndrome_name)
  }))
  results
}

# ==============================================================================
# 3. UI & SERVER
# ==============================================================================

ui <- dashboardPage(
  dashboardHeader(title = "FDA Pharmacovigilance"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Drug Safety (Panel A)", tabName = "panel_a", icon = icon("pills")),
      menuItem("Syndrome Analysis (Panel B)", tabName = "panel_b", icon = icon("vial"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML(".shiny-output-error-validation { color: #d9534f; font-weight: bold; }"))),
    tabItems(
      tabItem(tabName = "panel_a",
              fluidRow(
                box(width = 4, status = "primary", solidHeader = TRUE, title = "Settings",
                    textInput("drug_input", "Generic Drug Name:", "IBUPROFEN"),
                    actionButton("run_drug", "Run Analysis", class = "btn-block btn-primary")),
                box(width = 8, status = "info", title = "Hospitalization Demographics", plotlyOutput("demo_plot"))
              ),
              fluidRow(
                box(width = 12, status = "warning", title = "Disproportionality Analysis (ROR)", DTOutput("top_rx_table"))
              )),
      tabItem(tabName = "panel_b",
              fluidRow(
                box(width = 4, status = "danger", solidHeader = TRUE, title = "Settings",
                    textInput("syndrome_input", "MedDRA Condition:", "STEVENS-JOHNSON SYNDROME"),
                    actionButton("run_syndrome", "Find Culprits", class = "btn-block btn-danger"),
                    hr(),
                    helpText("Hint: Use British English (e.g. Diarrhoea)."),
                    uiOutput("meddra_link")),
                box(width = 8, status = "danger", title = "Reporting Odds Ratio Signals", DTOutput("ror_table"))
              ))
    )
  )
)

server <- function(input, output, session) {
  
  output$meddra_link <- renderUI({
    tagList("Check terms at:", a("MedDRA Web Search", href="https://www.meddra.org/how-to-use/tools/web-based-browser", target="_blank"))
  })
  
  drug_data <- eventReactive(input$run_drug, {
    clean_drug <- toupper(trimws(input$drug_input))
    validate(need(clean_drug != "", "Please enter a drug name."))
    withProgress(message = 'Calculating Signals...', {
      res <- analyze_drug_safety(clean_drug)
      validate(need(!is.null(res), "No data found. Use generic names (e.g. ASPIRIN)."))
      res
    })
  })
  
  output$demo_plot <- renderPlotly({
    req(drug_data())
    plot_ly(drug_data()$demo[sex != "Unknown"], x = ~age_grp, y = ~cases, color = ~sex, type = 'bar') %>%
      layout(barmode = 'group', yaxis = list(title = "Hospitalizations"))
  })
  
  output$top_rx_table <- renderDT({
    req(drug_data())
    datatable(drug_data()$top, options = list(dom = 't'), rownames = FALSE) %>%
      formatStyle('ROR', backgroundColor = styleInterval(c(2, 5), c('white', 'lightyellow', '#ffcccc')))
  })
  
  syndrome_data <- eventReactive(input$run_syndrome, {
    clean_syn <- toupper(trimws(input$syndrome_input))
    validate(need(clean_syn != "", "Please enter a condition."))
    withProgress(message = 'Mining database...', {
      res <- analyze_syndrome_panel_b(clean_syn)
      validate(need(!is.null(res) && nrow(res) > 0, "Condition not found. Check spelling (British English) or try a more general term."))
      res
    })
  })
  
  output$ror_table <- renderDT({
    req(syndrome_data())
    datatable(syndrome_data(), rownames = FALSE, colnames = c("Drug", "Cases", "ROR", "L95", "U95")) %>%
      formatStyle('ror', backgroundColor = styleInterval(c(2, 5), c('white', 'lightyellow', '#ffcccc')))
  })
}

shinyApp(ui, server)
