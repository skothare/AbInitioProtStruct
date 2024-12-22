# auto_Validation.R
# Group 2: Ab Initio Secondary Structure Prediction of Proteins
# Date: 12th December, 2024
# Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare

library(shiny)
library(ggplot2)
library(reshape2)

############################################
# HELPER FUNCTIONS
############################################

# Function to calculate overall accuracy (percentage of correctly predicted residues)
# Accuracy = (Number of correct predictions / Total residues) * 100
calculate_accuracy <- function(pred_seq, true_seq) {
  matches <- sum(strsplit(pred_seq, NULL)[[1]] == strsplit(true_seq, NULL)[[1]])
  return(matches / nchar(pred_seq) * 100)
}

# Function to calculate per-class precision, recall, and F1-score
# Precision = TP / (TP + FP)
# Recall = TP / (TP + FN)
# F1 = 2 * (Precision * Recall) / (Precision + Recall)
#
# TP: True Positives - correct identification of a class
# FP: False Positives - model predicted a class incorrectly
# FN: False Negatives - model missed instances of a class
#
# "types" is a vector of classes you consider (e.g., c("H","E","C","T"))
calculate_metrics <- function(pred_seq, true_seq, types = c("H", "E", "C", "T")) {
  pred_chars <- strsplit(pred_seq, NULL)[[1]]
  true_chars <- strsplit(true_seq, NULL)[[1]]
  
  # Initialize a data frame to store metrics for each type
  metrics <- data.frame(Type = types, Precision = NA, Recall = NA, F1_Score = NA)
  
  for (i in 1:length(types)) {
    type <- types[i]
    
    # Calculate True Positives (TP), False Positives (FP), and False Negatives (FN) for the current type
    TP <- sum(pred_chars == type & true_chars == type)
    FP <- sum(pred_chars == type & true_chars != type)
    FN <- sum(pred_chars != type & true_chars == type)
    
    # Calculate Precision and Recall, handling division by zero
    Precision <- ifelse((TP + FP) > 0, (TP / (TP + FP)) * 100, NA)
    Recall <- ifelse((TP + FN) > 0, (TP / (TP + FN)) * 100, NA)
    
    # Calculate F1-Score only if both Precision and Recall are defined
    F1_Score <- if (!is.na(Precision) && !is.na(Recall) && (Precision + Recall) > 0) {
      2 * (Precision * Recall) / (Precision + Recall)
    } else {
      NA
    }
    
    # Assign calculated metrics to the data frame
    metrics$Precision[i] <- Precision
    metrics$Recall[i] <- Recall
    metrics$F1_Score[i] <- F1_Score
  }
  
  return(metrics)
}

############################################
# SHINY UI
############################################

ui <- fluidPage(
  titlePanel("Validation of Secondary Structure Prediction Models"),
  sidebarLayout(
    sidebarPanel(
      fileInput("csv_file", "Upload CSV File:", accept = ".csv"),
      actionButton("predict", "Predict")
    ),
    mainPanel(
      # Plot 1: Accuracy per Protein and Model
      plotOutput("accuracyPlot"),
      
      # Table 1: Average Accuracy and Standard Deviation per Model
      tableOutput("averageAccuracyTable"),
      
      # Table 2: Average Precision, Recall, and F1-Score per Model
      tableOutput("modelMetricsTable"),
      
      # Table 3: Average Precision, Recall, and F1-Score per Structure Type and Model
      tableOutput("structureMetricsTable"),
      
      # Plot 2: Average Accuracy by Model
      plotOutput("averageAccuracyPlot"),
      
      # Plot 3: Accuracy vs Protein Length
      plotOutput("lengthAccuracyPlot"),
      
      # Plot 4: Precision per Structure Type and Model
      plotOutput("precisionPlot"),
      
      # Plot 5: Recall per Structure Type and Model
      plotOutput("recallPlot"),
      
      # Plot 6: F1-Score by Structure Type and Model
      plotOutput("f1ScorePlot")
    )
  )
)

############################################
# SHINY SERVER
############################################

server <- function(input, output) {
  
  observeEvent(input$predict, {
    
    req(input$csv_file)
    
    # Read CSV file: expected columns - ProteinName, ProteinSequence, DSSPSequence
    csv_data <- read.csv(input$csv_file$datapath, stringsAsFactors = FALSE)
    
    # Initialize lists to store metrics for each protein/model
    metrics_list <- list()
    
    # Models we are analyzing (adjust as needed)
    models <- c("Chou-Fasman", "GOR", "HMM")
    
    for (i in 1:nrow(csv_data)) {
      protein_name <- csv_data[i, "ProteinName"]
      protein_sequence <- gsub("\\s+", "", csv_data[i, "ProteinSequence"])  # Clean whitespace
      dssp_sequence <- gsub("\\s+", "", csv_data[i, "DSSPSequence"])        # Clean whitespace
      
      # Ensure sequences match in length
      if (nchar(protein_sequence) != nchar(dssp_sequence)) {
        warning(paste("Sequence length mismatch for", protein_name))
        next
      }
      
      # Call the Go executable with the input sequence (adjust command as necessary)
      # Ensure that the Go executable is accessible and has execute permissions
      go_command <- sprintf("./AbInitioPS %s", protein_sequence)
      prediction_result <- system(go_command, intern = TRUE)
      
      # Extract predictions for each model from the command output
      predictions <- list()
      for (model in models) {
        line <- grep(model, prediction_result, value = TRUE)
        predictions[[model]] <- sub("^.*?:\\s*", "", line)
      }
      
      # For each model, calculate accuracy and per-class metrics (precision, recall, F1)
      for (model in models) {
        pred_seq <- predictions[[model]]
        acc <- calculate_accuracy(pred_seq, dssp_sequence)
        per_class_metrics <- calculate_metrics(pred_seq, dssp_sequence, types = c("H", "E", "C", "T"))
        
        # Add identifying info to the metrics data frame (protein, model, accuracy, length)
        per_class_metrics$Protein <- protein_name
        per_class_metrics$Model <- model
        per_class_metrics$Accuracy <- acc
        per_class_metrics$Length <- nchar(protein_sequence)
        
        # Store results in a list for later aggregation
        metrics_list[[length(metrics_list) + 1]] <- per_class_metrics
      }
    }
    
    # Combine all per-protein, per-model metrics into one data frame
    metrics_data <- do.call(rbind, metrics_list)
    
    # Convert F1_Score column to numeric (should already be numeric, but just to be safe)
    metrics_data$F1_Score <- as.numeric(metrics_data$F1_Score)
    
    ############################################
    # PLOT 1: Accuracy per Protein and Model
    # Shows how well each model performed on each protein overall
    # Error bars represent Standard Error of the Mean (SEM) per Protein and Model
    ############################################
    output$accuracyPlot <- renderPlot({
      # Calculate mean accuracy per Protein and Model
      protein_model_acc <- aggregate(Accuracy ~ Protein + Model, data = metrics_data, mean)
      
      # Calculate Standard Error (SEM) per Protein and Model
      sem_data <- aggregate(Accuracy ~ Protein + Model, data = metrics_data, function(x) {
        if (length(x) > 1) {
          return(sd(x) / sqrt(length(x)))
        } else {
          return(0)  # SEM is zero if only one observation exists
        }
      })
      colnames(sem_data) <- c("Protein", "Model", "SEM")
      
      # Merge SEM with mean accuracy data
      protein_model_acc <- merge(protein_model_acc, sem_data, by = c("Protein", "Model"))
      
      ggplot(protein_model_acc, aes(x = Protein, y = Accuracy, fill = Model)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        # Add error bars per Protein and Model
        geom_errorbar(aes(ymin = Accuracy - SEM, ymax = Accuracy + SEM),
                      width = 0.2, position = position_dodge(0.9)) +
        theme_minimal() +
        labs(title = "Accuracy of Structure Prediction Across Models",
             x = "Protein",
             y = "Accuracy (%)",
             fill = "Model") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    })
    
    ############################################
    # TABLE 1: Average Accuracy and Standard Deviation per Model
    # Displays the numerical values alongside the plot
    ############################################
    output$averageAccuracyTable <- renderTable({
      # Calculate average accuracy per model
      average_accuracy_data <- aggregate(Accuracy ~ Model, data = metrics_data, mean)
      # Calculate standard deviation for each model
      accuracy_sd <- aggregate(Accuracy ~ Model, data = metrics_data, sd)
      # Merge mean and SD data frames
      avg_sd_data <- merge(average_accuracy_data, accuracy_sd, by = "Model")
      # Rename columns for clarity
      colnames(avg_sd_data) <- c("Model", "Average Accuracy (%)", "Standard Deviation (%)")
      return(avg_sd_data)
    }, striped = TRUE, bordered = TRUE, hover = TRUE)
    
    ############################################
    # TABLE 2: Average Precision, Recall, and F1-Score per Model
    # Displays the numerical values alongside the Precision, Recall, and F1-Score plots
    ############################################
    output$modelMetricsTable <- renderTable({
      # Calculate average Precision, Recall, and F1-Score per model
      avg_metrics <- aggregate(cbind(Precision, Recall, F1_Score) ~ Model, data = metrics_data, mean, na.rm = TRUE)
      # Rename columns for clarity
      colnames(avg_metrics) <- c("Model", "Average Precision (%)", "Average Recall (%)", "Average F1-Score (%)")
      return(avg_metrics)
    }, striped = TRUE, bordered = TRUE, hover = TRUE)
    
    ############################################
    # TABLE 3: Average Precision, Recall, and F1-Score per Structure Type and Model
    # Displays the numerical values alongside the Precision, Recall, and F1-Score plots
    ############################################
    output$structureMetricsTable <- renderTable({
      # Calculate average Precision, Recall, and F1-Score per Model and Structure Type
      structure_metrics <- aggregate(cbind(Precision, Recall, F1_Score) ~ Model + Type, data = metrics_data, mean, na.rm = TRUE)
      # Rename columns for clarity
      colnames(structure_metrics) <- c("Model", "Structure Type", "Average Precision (%)", "Average Recall (%)", "Average F1-Score (%)")
      return(structure_metrics)
    }, striped = TRUE, bordered = TRUE, hover = TRUE)
    
    ############################################
    # PLOT 2: Average Accuracy by Model
    # Shows mean accuracy per model with error bars representing standard deviation
    ############################################
    # Calculate average accuracy per model
    average_accuracy_data <- aggregate(Accuracy ~ Model, data = metrics_data, mean)
    # Calculate standard deviation for error bars
    accuracy_sd <- aggregate(Accuracy ~ Model, data = metrics_data, sd)
    average_accuracy_data$Error <- accuracy_sd$Accuracy
    
    output$averageAccuracyPlot <- renderPlot({
      ggplot(average_accuracy_data, aes(x = Model, y = Accuracy, fill = Model)) +
        geom_bar(stat = "identity") +
        geom_errorbar(aes(ymin = Accuracy - Error, ymax = Accuracy + Error),
                      width = 0.2) +
        theme_minimal() +
        labs(title = "Average Accuracy by Model",
             x = "Model",
             y = "Average Accuracy (%)")
    })
    
    ############################################
    # PLOT 3: Accuracy vs Protein Length
    # Shows if there's any correlation between protein length and prediction accuracy
    ############################################
    # Prepare data: average accuracy per Protein and Model
    protein_model_acc <- aggregate(Accuracy ~ Protein + Model + Length, data = metrics_data, mean)
    
    output$lengthAccuracyPlot <- renderPlot({
      ggplot(protein_model_acc, aes(x = Length, y = Accuracy, color = Model, shape = Model)) +
        geom_point(size = 2) +
        theme_minimal() +
        labs(title = "Accuracy vs Protein Length",
             x = "Protein Length (Residues)",
             y = "Accuracy (%)")
    })
    
    ############################################
    # PLOT 4: Precision per Structure Type and Model
    # Shows how precise each model is for each structure type
    ############################################
    # Rearrange the 'Type' factor to desired order: H, E, C, T
    metrics_data$Type <- factor(metrics_data$Type, levels = c("H", "E", "C", "T"))
    
    # Aggregate precision by Model and Type
    precision_data <- aggregate(Precision ~ Model + Type, data = metrics_data, mean, na.rm = TRUE)
    
    output$precisionPlot <- renderPlot({
      ggplot(precision_data, aes(x = Type, y = Precision, fill = Model)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        theme_minimal() +
        labs(
          title = "Average Precision by Structure Type and Model",
          x = "Structure Type",
          y = "Precision (%)",
          fill = "Model"
        ) +
        scale_x_discrete(limits = c("H", "E", "C", "T"))  # Ensure order is H, E, C, T
    })
    
    ############################################
    # PLOT 5: Recall per Structure Type and Model
    # Shows how well each model recalls each structure type
    ############################################
    # Aggregate recall by Model and Type
    recall_data <- aggregate(Recall ~ Model + Type, data = metrics_data, mean, na.rm = TRUE)
    
    output$recallPlot <- renderPlot({
      ggplot(recall_data, aes(x = Type, y = Recall, fill = Model)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        theme_minimal() +
        labs(
          title = "Average Recall by Structure Type and Model",
          x = "Structure Type",
          y = "Recall (%)",
          fill = "Model"
        ) +
        scale_x_discrete(limits = c("H", "E", "C", "T"))  # Ensure order is H, E, C, T
    })
    
    ############################################
    # PLOT 6: F1-Score by Structure Type and Model
    # Combines precision and recall into a single metric
    ############################################
    # Aggregate F1_Score by Model and Type
    f1_data <- aggregate(F1_Score ~ Model + Type, data = metrics_data, mean, na.rm = TRUE)
    
    output$f1ScorePlot <- renderPlot({
      ggplot(f1_data, aes(x = Type, y = F1_Score, fill = Model)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        theme_minimal() +
        labs(
          title = "Average F1-Score by Structure Type and Model",
          x = "Structure Type",
          y = "F1-Score (%)",
          fill = "Model"
        ) +
        scale_x_discrete(limits = c("H", "E", "C", "T"))  # Ensure order is H, E, C, T
    })
    
  })
}

# Run the Shiny application
shinyApp(ui = ui, server = server)
