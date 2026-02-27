library(shiny)
library(surveillance)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)

load("kerala_fire_model.RData")

ui <- page_sidebar(
  title = "Kerala Thermal Anomaly Intelligence",
  theme = bs_theme(version = 5, bootswatch = "darkly"),
  
  sidebar = sidebar(
    title = "Analysis Controls",
    selectInput(
      inputId = "selected_unit",
      label = "Select Area:",
      choices = c("Kerala", "Alappuzha", "Ernakulam", "Idukki", "Kannur", "Kasaragod", 
                  "Kollam", "Kottayam", "Kozhikode", "Malappuram", "Palakkad", 
                  "Pathanamthitta", "Thiruvananthapuram", "Thrissur", "Wayanad"),
      selected = "Kerala"
    ),
    sliderInput(
      inputId = "year_range",
      label = "Year Range:",
      min = 2021, max = 2026,
      value = c(2021, 2026),
      sep = ""
    ),
    hr(),
    uiOutput("stats_display"),
    hr(),
    
    # Expandable Information Section to replace the report
    accordion(
      accordion_panel(
        "1. The 1km Grid Rule (Data Source)",
        p("Data is sourced from NASA FIRMS (VIIRS/MODIS) thermal anomaly sensors. The y-axis represents the raw count of 1km by 1km grid squares that emitted extreme heat during that week. The satellite does not measure the physical size of the fire. A high number means many distinct grids were triggered, which could be dozens of small fires or a few massive ones stretching across multiple kilometers.")
      ),
      accordion_panel(
        "2. Agriculture vs. Wildfire",
        p("These detections represent extreme thermal anomalies, not strictly forest fires. The context changes by geography. In heavily forested districts like Wayanad or Idukki, spikes generally indicate true vegetation wildfires. However, in wetland and coastal districts with zero forest cover like Alappuzha, these hotspots represent seasonal agricultural clearing, specifically post-harvest paddy stubble burning in the Kuttanad region.")
      ),
      accordion_panel(
        "3. Mathematical Risk Layers",
        p("The model deconstructs risk into three sources. Endemic (Grey) is the baseline environmental risk driven by constant seasonal factors. Autoregressive (Blue) represents local heat recurrence, meaning a fire in a district increases the odds of more fires in that exact same district. Spatiotemporal (Orange) is neighborhood spread, representing the risk of fire activity spilling over the border from directly adjacent districts.")
      )
    )
  ),
  
  card(
    full_screen = TRUE,
    card_header(span(bsicons::bs_icon("layers"), " Deconstructed Thermal Risk Analysis")),
    plotOutput("distPlot", height = "600px")
  )
)

server <- function(input, output) {
  
  plot_data <- reactive({
    is_total <- (input$selected_unit == "Kerala")
    
    fit_end <- predict(fire_model_FINAL, type = "endemic")
    fit_ar  <- predict(fire_model_FINAL, type = "epi.own")
    fit_ne  <- predict(fire_model_FINAL, type = "epi.neighbours")
    
    raw_obs <- tail(fire_sts@observed, nrow(fit_end))
    exact_dates <- 2021 + (seq_len(nrow(fit_end))) / 52
    
    if (is_total) {
      df <- data.frame(
        Date = exact_dates,
        Endemic = rowSums(fit_end),
        Autoregressive = rowSums(fit_ar),
        Spatiotemporal = rowSums(fit_ne),
        Observed = rowSums(raw_obs)
      )
    } else {
      df <- data.frame(
        Date = exact_dates,
        Endemic = fit_end[, input$selected_unit],
        Autoregressive = fit_ar[, input$selected_unit],
        Spatiotemporal = fit_ne[, input$selected_unit],
        Observed = raw_obs[, input$selected_unit]
      )
    }
    
    df %>% filter(Date >= input$year_range[1] & Date <= input$year_range[2])
  })
  
  output$stats_display <- renderUI({
    df <- plot_data()
    tagList(
      p(strong("Region: "), input$selected_unit),
      p(strong("Thermal Anomalies in View: "), sum(df$Observed, na.rm = TRUE)),
      p(strong("Avg Weekly Risk Score: "), round(mean(df$Endemic + df$Autoregressive + df$Spatiotemporal, na.rm = TRUE), 2))
    )
  })
  
  output$distPlot <- renderPlot({
    df <- plot_data()
    
    df_long <- df %>%
      pivot_longer(
        cols = c(Endemic, Autoregressive, Spatiotemporal),
        names_to = "Component",
        values_to = "Risk"
      ) %>%
      mutate(Component = factor(Component, levels = c("Spatiotemporal", "Autoregressive", "Endemic")))
    
    ggplot() +
      geom_area(data = df_long, aes(x = Date, y = Risk, fill = Component), alpha = 0.8) +
      geom_point(data = df, aes(x = Date, y = Observed), color = "white", size = 1.2, alpha = 0.7) +
      scale_fill_manual(values = c(
        "Endemic" = "grey50", 
        "Autoregressive" = "#0055ff", 
        "Spatiotemporal" = "#ff8c00"
      )) +
      labs(y = "Thermal Anomalies Detected", x = "Year", fill = "Risk Source") +
      scale_x_continuous(breaks = seq(2021, 2026, 1)) +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_line(color = "gray30"),
        panel.grid.minor = element_blank(),
        text = element_text(color = "white"),
        axis.text = element_text(color = "gray80"),
        legend.position = "top",
        legend.background = element_rect(fill = "#212529", color = NA),
        legend.text = element_text(color = "white"),
        plot.background = element_rect(fill = "#212529", color = NA),
        panel.background = element_rect(fill = "#212529", color = NA)
      )
  })
}

shinyApp(ui = ui, server = server)