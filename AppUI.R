# AppUI.R
# Group 2: Ab Initio Secondary Structure Prediction of Proteins
# Date: 12th December, 2024
# Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare

# Load required libraries
library(shiny)       # For creating interactive web applications
library(ggplot2)     # For creating visualizations
library(reshape2)    # For data manipulation and reshaping

# Function to calculate the match percentage between two sequences
# Input: Two character sequences (seq1 and seq2)
# Output: Percentage of matching characters
calculate_match_percentage <- function(seq1, seq2) {
  matches <- sum(strsplit(seq1, NULL)[[1]] == strsplit(seq2, NULL)[[1]])
  return(matches / nchar(seq1) * 100)  # Match percentage
}

# Function to prepare data for graph visualization
# Input: A list of predictions and selected models
# Output: Data frame suitable for plotting comparisons
prepare_comparison_data <- function(predictions, options) {
  data <- data.frame()
  for (option in options) {
    pred_seq <- predictions[[option]]
    data <- rbind(data, data.frame(
      Model = option,  # Model name
      Position = 1:nchar(pred_seq),  # Sequence positions
      Structure = strsplit(pred_seq, NULL)[[1]]  # Secondary structure
    ))
  }
  return(data)
}

# Function to read a FASTA file and extract the sequence
# Input: File path to the FASTA file
# Output: Protein sequence as a single string
read_fasta <- function(file_path) {
  lines <- readLines(file_path)  # Read all lines from the file
  sequence <- paste(lines[-1], collapse = "")  # Combine all lines except the header
  return(sequence)
}

# Function to calculate the character composition of a sequence
# Input: A protein secondary structure sequence
# Output: Percentage of each structure type ('H', 'E', 'C', 'T')
calculate_character_composition <- function(sequence) {
  total_length <- nchar(sequence)  # Total sequence length
  composition <- sapply(c('H', 'E', 'C', 'T'), function(char) {
    sum(strsplit(sequence, NULL)[[1]] == char) / total_length * 100  # Percentage calculation
  })
  return(composition)
}

# Define the User Interface (UI) for the application
ui <- fluidPage(
  titlePanel("Protein Secondary Structure Prediction"),
  sidebarLayout(
    sidebarPanel(
      # Input options: Enter sequence or upload FASTA file
      radioButtons("input_type", "Input Type:",
                   choices = c("Enter Protein Sequence" = "text",
                               "Upload FASTA File" = "file")),
      # Conditional input for text sequence
      conditionalPanel(
        condition = "input.input_type == 'text'",
        textInput("input_string", "Enter Protein Sequence:")
      ),
      # Conditional input for file upload
      conditionalPanel(
        condition = "input.input_type == 'file'",
        fileInput("fasta_file", "Upload FASTA File:")
      ),
      # Checkbox for selecting prediction models
      checkboxGroupInput("model_choices", "Select Models:",
                         choices = c("Chou-Fasman", "GOR", "HMM")),
      # Options to show graphs
      checkboxInput("show_graph", "Show Comparison Graph", value = FALSE),
      checkboxInput("show_composition_graph", "Show Character Composition Graph", value = FALSE),
      # Predict button
      actionButton("predict", "Predict")
    ),
    mainPanel(
      # Outputs: Predicted sequences, match percentages, and plots
      verbatimTextOutput("predictedSequences"),
      verbatimTextOutput("matchPercentage"),
      conditionalPanel(
        condition = "input.show_graph == true",
        plotOutput("comparisonPlot")
      ),
      conditionalPanel(
        condition = "input.show_composition_graph == true",
        plotOutput("compositionPlot")
      )
    )
  )
)

# Define the server logic for the application
server <- function(input, output) {
  
  observeEvent(input$predict, {
    req(input$model_choices)  # Ensure at least one model is selected
    
    # Determine input type and get sequence
    if (input$input_type == "text") {
      req(input$input_string)  # Ensure input string is not empty
      sequence <- input$input_string
    } else if (input$input_type == "file") {
      req(input$fasta_file)  # Ensure a file is uploaded
      sequence <- read_fasta(input$fasta_file$datapath)
    }
    
    # Call the Go executable with the input sequence
    go_command <- sprintf("./AbInitioPS %s", sequence)
    prediction_result <- system(go_command, intern = TRUE)  # Execute and capture output
    
    # Extract predictions based on selected models
    predictions <- list()
    for (choice in input$model_choices) {
      line <- grep(choice, prediction_result, value = TRUE)  # Find relevant lines
      predictions[[choice]] <- sub("^.*?:\\s*", "", line)  # Extract predictions
    }
    
    # Display predicted sequences
    output$predictedSequences <- renderText({
      paste(sapply(input$model_choices, function(choice) {
        paste(choice, "Prediction:", predictions[[choice]])
      }), collapse = "\n")
    })
    
    # Calculate and display match percentages between selected models
    match_texts <- c()
    if ("Chou-Fasman" %in% input$model_choices && "GOR" %in% input$model_choices) {
      cf_gor_percent <- calculate_match_percentage(predictions[["Chou-Fasman"]], predictions[["GOR"]])
      match_texts <- c(match_texts, paste("Percentage Matching Between Chou-Fasman and GOR:", round(cf_gor_percent, 2), "%"))
    }
    if ("Chou-Fasman" %in% input$model_choices && "HMM" %in% input$model_choices) {
      cf_hmm_percent <- calculate_match_percentage(predictions[["Chou-Fasman"]], predictions[["HMM"]])
      match_texts <- c(match_texts, paste("Percentage Matching Between Chou-Fasman and HMM:", round(cf_hmm_percent, 2), "%"))
    }
    if ("GOR" %in% input$model_choices && "HMM" %in% input$model_choices) {
      gor_hmm_percent <- calculate_match_percentage(predictions[["GOR"]], predictions[["HMM"]])
      match_texts <- c(match_texts, paste("Percentage Matching Between GOR and HMM:", round(gor_hmm_percent, 2), "%"))
    }
    
    output$matchPercentage <- renderText({
      paste(match_texts, collapse = "\n")
    })
    
    # Plot comparison graph if selected
    if (input$show_graph && length(input$model_choices) > 1) {
      comparison_data <- prepare_comparison_data(predictions, input$model_choices)
      
      output$comparisonPlot <- renderPlot({
        ggplot(comparison_data, aes(x = Position, y = Model)) +
          geom_tile(aes(fill = Structure), color = "white") +
          scale_fill_manual(values = c("H" = "red", "E" = "green", "T" = "blue", "C" = "purple")) +
          theme_minimal() +
          labs(title = "Predicted Secondary Structure Comparisons",
               x = "Position in Sequence",
               y = "",
               fill = "Structure")
      })
    }
    
    # Calculate character composition for each selected model
    composition_data <- data.frame()
    for (choice in input$model_choices) {
      comp <- calculate_character_composition(predictions[[choice]])
      composition_data <- rbind(composition_data, data.frame(
        Model = choice,
        Character = names(comp),
        Percentage = unlist(comp)
      ))
    }
    
    # Plot character composition graph if selected
    if (input$show_composition_graph) {
      output$compositionPlot <- renderPlot({
        ggplot(composition_data, aes(x = Character, y = Percentage, fill = Model)) +
          geom_bar(stat = "identity", position = "dodge") +
          theme_minimal() +
          labs(title = "Character Composition of Predicted Sequences",
               x = "Character Type",
               y = "Percentage (%)")
      })
    }
  })
}

# Run the Shiny application
shinyApp(ui = ui, server = server)
