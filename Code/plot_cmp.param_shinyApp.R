library(tidyverse)
library(lubridate)
library(shiny)

analyze_APEX_runs <- function(base_dir = "D:/APEX1905_New/PEST", num_runs = 30) {
  # File location for CROP parameter default values
  crop_file <- file.path(base_dir, "CROP1905.DAT")
  
  # Import parameter names
  param_cols <- read.delim(crop_file, sep = "", dec = ".", nrows = 1) %>%
    select(-CROP) %>%
    pivot_longer(values_to = "paramName", names_to = "index", everything()) %>%
    mutate(paramID = as.numeric(gsub(index, pattern = "X", replacement = ""))) %>%
    select(-index)
  
  # Import functional group names (CPNM)
  crop_names <- read.delim(crop_file, sep = "", dec = ".", skip = 167, nrows = 14,
                           header = F) %>%
    select(V1:V2) %>%
    rename(plantID = V1, plantName = V2)
  
  # Initialize empty lists to store data 
  cmp_data_list <- list()
  param_data_list <- list()
  rmse_list <- list()
  
  # Loop through the runs and import files
  for (i in seq_len(num_runs)) {
    cmp_file <- file.path(base_dir, paste0("CONUNN_AGM", i, ".cmp"))
    # Read the lines of the file into a character vector
    lines <- readLines(cmp_file)
    
    par2par_file <- file.path(base_dir, paste0("cropcase", i, ".par"))
    
    ## Import data from cmp_file (cage bm) and par2par_file (parameters used for each run)
    line_index <- grep("RMSE", lines)[1]
    
    cmp_data <- read.delim(skip = 1, text = lines[1:(line_index-1)], 
                           sep = "", dec = ".") %>%
      mutate(APEXrun = i, # ID column for run 
             date = ymd(paste(Year, Month, Day, sep = "-"))) # date format column
    
    # Each list component corresponds to 1 run
    cmp_data_list[[i]] <- cmp_data
    
    # Parameters used for each run
    par2par_lines <- read.delim(par2par_file, header = F, skip = 1)
    
    ## Parameter data
    # Split and convert par2par_data to a data frame
    par2par_values_ID <- strsplit(as.character(par2par_lines$V1), "\\s+")  # Split by whitespace
    param_data <- do.call(rbind, par2par_values_ID) %>% 
      as.data.frame() %>% # Convert to data frame
      mutate(plantID = as.numeric(substr(V1, 1, 3)),
             paramID = as.numeric(substr(V1, 4, 5)),
             paramValue = as.numeric(V2),
             APEXrun = i) %>%
      select(-c(V1:V4))
    
    # Add in parameter names
    param_data <- merge(param_data, param_cols)
    # Save values of each run in list components
    param_data_list[[i]] <- param_data
    
    ## Import RMSE data of each APEX run
    # Find the index of the last text string before the RMSE value
    target_line_index <- grep("T-5% CONFIDENCE", lines)[1]
    
    # Read in summary stats based on index
    stats <- read.delim(text = lines[target_line_index:length(lines)], 
                        sep = "", dec = ".")
    
    # Extract RMSE value
    rmse_value <- stats$RMSE[1] %>% as.numeric()
    rmse_list[[i]] <- rmse_value
  }
  
  # All of the experimental and simulated data (n = num_runs)
  combined_data <- bind_rows(cmp_data_list)
  
  # Parameters inputs used in each run
  param_input <- bind_rows(param_data_list) %>% 
    merge(crop_names)
  
  # Combining RMSE of all APEX runs
  rmse_data <- unlist(rmse_list) %>% 
    as.data.frame() %>%
    rename(RMSE = ".") %>%
    mutate(APEXrun = as.numeric(1:num_runs)) %>%
    merge(param_input, by = "APEXrun") %>%
    select(-paramID) %>%
    pivot_wider(names_from = paramName, values_from = paramValue)
  
  list(sim_vs_exp_data = combined_data, run_input = param_input, rmse_output = rmse_data)
}

# Shiny App UI
ui <- fluidPage(
  titlePanel("APEX ParaView"),
  sidebarLayout(
    sidebarPanel(
      textInput("base_dir", "Base Directory:", base_dir),
      numericInput("num_runs", "Number of Runs:", 30),
      textInput("crop_plot", "Plant group to plot:", ""),
      textInput("param_plot", "Parameter to plot:", "")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Visualization", plotOutput("apex_plot")),
        tabPanel("RMSE Data", 
                 dataTableOutput("rmse_table"),  # Use dataTableOutput for better table rendering
                 downloadButton("downloadTable", "Download Table")
        ),
        tabPanel("RMSE Parameter Plots", 
                 plotOutput("rmse_scatter_plot", height = 400),
                 tabPanel("Lowest RMSE Plot: 1 Plant Group", plotOutput("lowest_rmse_plot")),
                 tabPanel("Lowest RMSE Plot: All Plant Groups", plotOutput("lowest_rmse_plot_inter"))
        )
      )
    )
  )
)

# Shiny App Server
server <- function(input, output) {
  APEX.results <- reactive({
    analyze_APEX_runs(base_dir = input$base_dir,
                      num_runs = input$num_runs)
  })
  
  # Visualization of simulated vs experimental (field) data for all APEX runs
  output$apex_plot <- renderPlot({
    # Pulling out the cage biomass results
    simp_vs_exp_data <- APEX.results()[[1]]
    
    # Subsetting based on plant functional group
    filtered_data <- simp_vs_exp_data %>%
      filter(CPNM2 == input$crop_plot)
    
    # Plotting all APEX LHS runs vs experimental (field) data (black)
    ggplot(filtered_data) +
      geom_line(aes(x = date, y = S_BIOMAS, color = as.factor(APEXrun))) +
      geom_line(aes(x = date, y = Biomass), color = "black",
                linewidth = 1.5) +
      # geom_smooth(aes(x = date, y = S_BIOMAS)) +
      theme_bw() +
      theme(text = element_text(size = 15, family = "serif")) +
      xlab("Year") +
      ylab(paste("Cage Biomass (kg/ha):", input$crop_plot)) +
      labs(color = "APEX LHS Run") +
      scale_y_continuous(limits = c(0, NA)) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  })
  
  # Output table of all APEX runs' RMSE values with their corresponding parameter inputs
  output$rmse_table <- renderDataTable({
    rmse_df <- APEX.results()[[3]] %>%
      filter(plantName == input$crop_plot)
  })
  
  rmse_download <- reactive({
    APEX.results()[[3]]  %>%
      filter(plantName == input$crop_plot)
  })
  
  # Download handler to download the RMSE and parameter input data
  output$downloadTable <- downloadHandler(
    
    filename = function() {
      paste("rmse_table_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(rmse_download(), file)
    }
  )
  
  # Point plot of parameter values vs RMSE for ONE plant group
  # For interpreting influence of a parameter on model fitness
  output$rmse_scatter_plot <- renderPlot({
    crop_plot <- input$crop_plot
    param_plot <- input$param_plot
    
    rmse_df <- APEX.results()[[3]] %>% 
      filter(plantName == crop_plot) %>%
      select(APEXrun, RMSE, plantName, all_of(param_plot))
    
    ggplot() +
      geom_line(aes(x = rmse_df[[4]], y = rmse_df$RMSE)) +
      geom_point(aes(x = rmse_df[[4]], y = rmse_df$RMSE)) +
      theme_bw() +
      theme(text = element_text(size = 15, family = "serif")) +
      xlab(paste("All", param_plot, "Parameter Values:", crop_plot)) +
      ylab("RMSE") 
    
  })
  
  # Visualization of simulated vs experimental (field) data for ONE plant group
  # for ONLY the APEX run w/ the lowest RMSE
  output$lowest_rmse_plot <- renderPlot({
    simp_vs_exp_data <- APEX.results()[[1]]
    rmse_df <- APEX.results()[[3]]
    
    filtered_data <- simp_vs_exp_data %>% 
      filter(APEXrun == as.numeric(rmse_df[which.min(rmse_df$RMSE), "APEXrun"]),
             CPNM2 == input$crop_plot)
    
    
    lowest_rmse_value <- rmse_df[which.min(rmse_df$RMSE), "RMSE"] %>% as.numeric()
    
    ggplot(filtered_data) +
      geom_line(aes(x = date, y = S_BIOMAS, color = as.factor(APEXrun))) +
      geom_point(aes(x = date, y = S_BIOMAS, color = as.factor(APEXrun))) +
      geom_line(aes(x = date, y = Biomass), color = "black") +
      geom_point(aes(x = date, y = Biomass), color = "black") +
      theme_bw() +
      theme(text = element_text(size = 15, family = "serif")) +
      xlab("Year") +
      ylab(paste("Cage Biomass (kg/ha):", input$crop_plot)) +
      labs(color = "APEX LHS Run") +
      scale_y_continuous(limits = c(0, NA)) +
      geom_text(x = max(filtered_data$date), y = (max(filtered_data$Biomass) - 50), 
                label = paste("RMSE:", round(lowest_rmse_value, 2)), 
                hjust = 1, vjust = 1, size = 5,  family = "serif")
    
    
  })
  
  # Visualization of simulated data for ALL plants groups that were parameterized
  # for ONLY the APEX run w/ the lowest RMSE
  output$lowest_rmse_plot_inter <- renderPlot({
    simp_vs_exp_data <- APEX.results()[[1]]
    rmse_df <- APEX.results()[[3]]
    
    filtered_data <- simp_vs_exp_data %>% 
      filter(APEXrun == as.numeric(rmse_df[which.min(rmse_df$RMSE), "APEXrun"]))
    
    
    lowest_rmse_value <- rmse_df[which.min(rmse_df$RMSE), "RMSE"] %>% as.numeric()
    
    ggplot(filtered_data) +
      geom_line(aes(x = date, y = S_BIOMAS, color = CPNM2)) +
      geom_point(aes(x = date, y = S_BIOMAS, color = CPNM2)) +
      theme_bw() +
      theme(text = element_text(size = 15, family = "serif")) +
      xlab("Year") +
      ylab(paste("Cage Biomass (kg/ha)")) +
      labs(color = "Plant Group") +
      scale_y_continuous(limits = c(0, NA)) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      geom_text(x = max(filtered_data$date), y = (max(filtered_data$S_BIOMAS) - 50), 
                label = paste("RMSE:", round(lowest_rmse_value, 2)), 
                hjust = 1, vjust = 1, size = 5,  family = "serif")
    
    
  })
}

shinyApp(ui = ui, server = server)
