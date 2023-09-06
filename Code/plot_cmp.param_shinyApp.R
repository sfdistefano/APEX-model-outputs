library(tidyverse)
library(lubridate)
library(shiny)

base_dir = "D:/APEX1905_New/PEST"

# Define the analyze_APEX_runs function
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
  
  ## Import RMSE data of each APEX run
  # Read in summary stats based on index
  rmse_file <- file.path(base_dir, "CONUNN_AGM_ALL.cmp")
  stats <- read.delim(rmse_file, 
                      sep = "", dec = ".")
  
  # Extract RMSE value
  rmse_values <- stats %>% 
    select(Case., RMSECAGE) %>%
    rename(APEXrun = Case., RMSE = RMSECAGE)
  
  # Initialize empty lists to store data 
  cmp_data_list <- list()
  param_data_list <- list()
  
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
    
    # Print a message when the run is completed
    cat("Run", i, "completed.\n\n")
    
  }
  
  # All of the experimental and simulated data (n = num_runs)
  combined_data <- bind_rows(cmp_data_list)
  
  # Parameters inputs used in each run
  param_input <- bind_rows(param_data_list) %>% 
    merge(crop_names)
  
  # Combining RMSE of all APEX runs
  rmse_data <- rmse_values %>%
    merge(param_input, by = "APEXrun") %>%
    select(-paramID) %>%
    pivot_wider(names_from = paramName, values_from = paramValue)
  
  # Load pastID_ecosite data
  pastID_ecosite <- read.csv("D:/APEX data and scripts/Data/PastureID_ecosite_34subareas.csv")
  
  ## original cage biomass used (day*pasture*plot*functional group)
  biomass.plot <- read.csv("D:/APEX data and scripts/Data/CPER Biomass/AGM_Biomass_Widecln_attr_2023-07-18.csv") %>% 
    rename(Year = YearSampled) %>%
    filter(Year >= 2014)  %>% # filter for year CARM started
    pivot_longer(cols = AG:NA., names_to = "FGCode", values_to = "kgPerHa") %>%
    group_by(Year, Pasture, Ecosite, Treatment, Block, Plot, FGCode) %>%
    summarize(kgPerHa = mean(kgPerHa))%>% # summarize by plot of each pasture
    filter(!(Plot == 5 & Pasture == "18S")) %>% # prescribed burn plots
    filter(!(Plot == 6 & Pasture == "18S")) %>%
    filter(!(Plot == 5 & Pasture == "19N")) %>%
    filter(!(Plot == 6 & Pasture == "19N"))
  
  # preparing dataframe for plotting SE across plots
  biomass.plot_info <- merge(biomass.plot, pastID_ecosite,
                             by = c("Pasture", "Ecosite"), 
                             all.x = TRUE) %>%
    mutate(Month = 8, Day = 12, date = ymd(paste(Year, Month, Day, sep = "-"))
    ) %>% # adding date info, when collected in field
    rename(CPNM = FGCode, Y = Year)
  
  biomass.plot_info$CPNM2 <- biomass.plot_info$CPNM %>% 
    gsub(pattern = "C3PG", replacement = "CSPG") %>%
    gsub(pattern = "FORB", replacement = "FRB3") %>%
    gsub(pattern = "SS", replacement = "SSHB") %>%
    gsub(pattern = "AG", replacement = "VUOC")
  
  # biomass data is summarized to day*pasture*functional group
  biomass.pasture <-  biomass.plot %>% 
    mutate(Month = 8, Day = 12,
           date = ymd(paste(Year, Month, Day, sep = "-"))
    ) %>% # adding date info, when collected in field
    group_by(date, Year, Pasture, Ecosite, Treatment, FGCode) %>%
    summarize(Biomass = mean(kgPerHa)) %>% # summarized by pasture
    rename(CPNM = FGCode, Y = Year) %>%
    filter(!(Pasture %in% c("1W", "28N", "32W"))) %>% # pastures not included in CARM or TRM
    filter(!(CPNM == "SD")) # removing standing dead (growth from previous season)
  
  ## adding pasture ID to match with APEX output files
  biomass.pasture <- merge(biomass.pasture, pastID_ecosite, 
                           by = c("Pasture", "Ecosite"), all.x = TRUE) 
  
  ## fixing functional groups names so they match APEX names
  biomass.pasture$CPNM2 <- biomass.pasture$CPNM %>% 
    gsub(pattern = "C3PG", replacement = "CSPG") %>%
    gsub(pattern = "FORB", replacement = "FRB3") %>%
    gsub(pattern = "SS", replacement = "SSHB") %>%
    gsub(pattern = "AG", replacement = "VUOC")
  
  biomass.ecosite <- biomass.pasture %>% 
    group_by(date, Y, Ecosite, CPNM2) %>%
    summarize(Biomass_es = mean(Biomass))
  
  
  list(sim_vs_exp_data = combined_data, 
       run_input = param_input, 
       rmse_output = rmse_data,
       pastID_ecosite = pastID_ecosite,
       biomass_plot = biomass.plot,
       biomass_ecosite = biomass.ecosite,
       biomass_plot_info = biomass.plot_info)
}

# Shiny App UI
ui <- fluidPage(
  titlePanel("APEX ParamaView"),
  sidebarLayout(
    sidebarPanel(
      textInput("base_dir", "Base Directory:", base_dir),
      numericInput("num_runs", "Number of Runs:", 30),
      textInput("crop_plot", "Plant group to plot:", ""),
      textInput("param_plot", "Parameter to plot:", "")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("RMSE Data", 
                 dataTableOutput("rmse_table"),  # Use dataTableOutput for better table rendering
                 downloadButton("downloadTable", "Download Table")
        ),
        tabPanel("RMSE Parameter Plots", 
                 plotOutput("rmse_scatter_plot"),
                 plotOutput("lowest_rmse_plot"),
                 plotOutput("lowest_rmse_plot_inter")
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
  
  # Output table of all APEX runs' RMSE values with their corresponding parameter inputs
  output$rmse_table <- renderDataTable({
    rmse_df <- APEX.results()[[3]] %>%
      filter(plantName == input$crop_plot)
  })
  
  # Download handler to download the RMSE and parameter input data
  rmse_download <- reactive({
    APEX.results()[[3]]  %>%
      filter(plantName == input$crop_plot)
  })
  
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
    rmse_df <- APEX.results() [[3]]
    pastID_ecosite <- APEX.results()[[4]]
    biomass.plot <- APEX.results()[[5]]
    biomass.ecosite <- APEX.results()[[6]]
    biomass.plot_info <- APEX.results()[[7]]
    
    # Filter data to the run of lowest RMSE value
    filtered_data <- simp_vs_exp_data %>% 
      filter(APEXrun == as.numeric(rmse_df[which.min(rmse_df$RMSE), "APEXrun"]),
             CPNM2 == input$crop_plot) %>%
      merge(pastID_ecosite, by = "PastureID")
    lowest_rmse_value <- rmse_df[which.min(rmse_df$RMSE), "RMSE"] %>% as.numeric()
    
    # Pull out the identifiers to filter the experimental data
    es.name <- filtered_data$Ecosite %>% unique()
    func.group <- filtered_data$CPNM2 %>% unique()
    prop.num <- filtered_data$Proportion %>% unique()
    
    # Filter the experimental data
    # mean calculated at ecosite level
    field01 <- biomass.ecosite %>% 
      filter(Ecosite %in% es.name & CPNM2 %in% func.group) %>%
      mutate(Biomass_wt = Biomass_es*prop.num) %>%
      group_by(Y, date, CPNM2) %>%
      summarize(Biomass_fg = sum(Biomass_wt, na.rm = T)) # sum of 1 func.group across ecosites w/in pasture
    
    # SE calculated from all plots of ecosite
    field02 <- biomass.plot_info %>%
      filter(Ecosite %in% es.name & CPNM2 %in% func.group) %>%
      group_by(date, Y, CPNM2) %>%
      summarize(se = sd(kgPerHa, na.rm = T)/sqrt(length(kgPerHa)))
    
    # Plotting results
    ggplot(filtered_data) +
      geom_line(data = filtered_data, 
                aes(x = date, y = S_BIOMAS, color = paste0(input$crop_plot, "-APEX"))) +
      geom_point(data = filtered_data, 
                 aes(x = date, y = S_BIOMAS, color = paste0(input$crop_plot, "-APEX"))) +
      geom_point(data = field01, 
                 aes(x = date, y = Biomass_fg, color = paste0(input$crop_plot, "-CPER"))) +
      geom_line(data = field01, 
                aes(x = date, y = Biomass_fg, color = paste0(input$crop_plot, "-CPER"))) +
      geom_errorbar(data = field02, 
                    aes(x = date, ymin = field01$Biomass_fg - se, ymax = field01$Biomass_fg + se, 
                        color = paste0(input$crop_plot, "-CPER")), 
                    width = 50) +
      theme_bw() +
      theme(text = element_text(size = 15, family = "serif")) +
      xlab("Year") +
      ylab(paste("Cage Biomass (kg/ha):", input$crop_plot)) +
      labs(color = "Data Source") +
      scale_y_continuous(limits = c(0, NA)) +
      geom_text(x = max(filtered_data$date), y = (max(filtered_data$Biomass) - 50), 
                label = paste("RMSE:", round(lowest_rmse_value, 2)), 
                hjust = 1, vjust = 1, size = 5,  family = "serif") +
      scale_x_date(date_breaks = "1 year",date_labels = "%Y") +
      scale_color_brewer(palette = "Dark2")
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
