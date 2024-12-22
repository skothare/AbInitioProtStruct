# Group 2: Protein Secondary Structure Prediction
# Date: 12th December, 2024
# Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare

# Demo link: https://drive.google.com/file/d/1X8fdxQFILHfyiLvSTlYzYjZ6NzGQ1mI8/view?usp=drive_link 

## Project Motivation
In this project, we were motivated by a central theme: predicting the secondary structure of proteins from their amino acid sequences using straightforward, pre-defined heuristic methods. Rather than relying on modern, data-intensive machine learning approaches, we chose to revisit and implement earlier, “classic” algorithms—such as Chou-Fasman, GOR, and simple Hidden Markov Models (HMMs)—to assess their predictive accuracy against established “ground-truth” or industry-standard predictions (e.g., derived from DSSP or reliable online prediction tools). This historical and heuristic-driven perspective allows us to understand how early computational biology tackled a complex prediction problem without the luxury of advanced machine learning frameworks.

## Computational Problem
Given an input amino acid sequence (a string of single-letter amino acid codes), produce an output string indicating the secondary structure state (e.g., H, E, T, C for helix, strand, turn, coil) for each residue. The objective is to maximize the percentage match between this predicted structure and the actual experimentally validated structure. We operate under minimal computational overhead and use models that rely on relatively simple statistical insights or state-based frameworks, without large-scale model training.

## Folder Structure

```
AbInitioPS/
├── Figures/
│   ├── Miscellaneous ".png" files
├── GOR_InfoVals/
│   ├── InfoVal_aHelix.csv
│   ├── InfoVal_bStrand.csv
│   ├── InfoVal_bTurn.csv
│   ├── InfoVal_Coil.csv
├── GORTests/
│   ├── GORPredict/
│   │   ├── Input/
│   │   ├── Output/
│   ├── SlideWindow/
│       ├── Input/
│       ├── Output/
├── .RData
├── .Rhistory
├── AccuracyTestDataset_50.csv
├── AppUI.R
├── auto_Validation.R
├── CF_functions_test.go
├── CF_functions.go
├── datatypes.go
├── EM_main.go
├── GOR_functions_test.go
├── GOR_functions.go
├── AbInitioPS
├── hmm_functions_test.go
├── HMM_functions.go
├── main.go
├── README.md
```

## Description of Files

- **`CF_functions.go`**: Implements the Chou-Fasman algorithm for secondary structure prediction.
- **`CF_functions_test.go`**: Unit tests for `CF_functions.go`.
- **`datatypes.go`**: Contains shared data types used across different modules.
- **`EM_main.go`**: Runs the Expectation-Maximization (EM) algorithm for HMM training.
- **`GOR_functions.go`**: Implements the GOR method for secondary structure prediction.
- **`GOR_functions_test.go`**: Unit tests for `GOR_functions.go`.
- **`HMM_functions.go`**: Implements the HMM model and algorithms for secondary structure prediction.
- **`hmm_functions_test.go`**: Unit tests for `HMM_functions.go`.
- **`main.go`**: Main entry point to the application. Integrates and executes different models.
- **`AppUI.R`**: R Shiny application for running the prediction algorithms via a user interface.
- **`auto_Validation.R`**: R Shiny application for validating model performance using metrics like precision, recall, and F1-score.
- **`AccuracyTestDataset_50.csv`**: Example dataset used for testing.
- **`GOR_InfoVals/`**: Contains CSV files with informational values for different secondary structure types (e.g., Helix, Strand, Turn, Coil).
- **`GORTests/GORPredict/Input/`**: Input test files for GOR prediction module.
- **`GORTests/GORPredict/Output/`**: Output test files for GOR prediction module.
- **`GORTests/SlideWindow/Input/`**: Input test files for sliding window functions.
- **`GORTests/SlideWindow/Output/`**: Output test files for sliding window functions.

## Instructions for Building and Running the Go Code

### Prerequisites
- **Go** (1.18 or higher)
- **Operating System**: macOS or Windows

### Build the Code
1. Open a terminal (macOS) or Command Prompt/PowerShell (Windows).
2. Navigate to the folder containing `main.go`.
3. Run the following command to build the application:
   ```sh
   go build -o Group2 main.go
   ```

### Run the Code
1. Execute the built binary with the input sequence as an argument:
   For macOS/Linux:
   ```sh
   ./Group2 "ETGTVPAINYLGAGYDHVRGNPVGDPSSMGDPGIRPPVLRF"
   ```
   For Windows:
   ```sh
   Group2.exe "ETGTVPAINYLGAGYDHVRGNPVGDPSSMGDPGIRPPVLRF"
   ```
   Replace the input sequence with your desired sequence.

## Running Tests

### Run all tests in the package:
```sh
go test ./...
```

### Run tests for a specific file (e.g., `CF_functions_test.go`):
```sh
go test -v CF_functions_test.go
```

## Using VS Code to Load and Run the Group2 Package

### Prerequisites
- Install **Visual Studio Code**.
- Install the **Go extension** for VS Code.
- Ensure **Go** is properly installed and added to your system PATH.

### Steps to Load the Group2 Package in VS Code
1. Open **Visual Studio Code**.
2. Click **File > Open Folder** and select the `Group2` folder.
3. Ensure the Go extension detects the project and prompts for setting up the environment.
4. Verify the `main.go` file is present in the project root.

### Build the Code in VS Code
1. Open the terminal in VS Code (**View > Terminal** or `Ctrl + ``).
2. Run the following command to build the project:
   ```sh
   go build -o Group2 main.go
   ```

### Run the Program in VS Code
1. In the terminal, execute the built binary with the input sequence:
   ```sh
   ./Group2 "ETGTVPAINYLGAGYDHVRGNPVGDPSSMGDPGIRPPVLRF"
   ```
   Replace the input sequence as needed.

### Run Tests in VS Code
1. Open the terminal.
2. Run the following command to execute all tests:
   ```sh
   go test ./...
   ```
3. To test a specific file, run:
   ```sh
   go test -v CF_functions_test.go
   ```

## Running R Shiny Applications

### Prerequisites
- **R** (version 4.1 or higher)
- **R Shiny package** installed (`install.packages("shiny")`)
- Additional R libraries: `ggplot2`, `reshape2`.

### Run `AppUI.R`
1. Open RStudio.
2. Load `AppUI.R`.
3. Click **Run App**.
4. Use the interface to:
   - Enter a protein sequence or upload a FASTA file.
   - Choose models (`Chou-Fasman`, `GOR`, or `HMM`).
   - View predicted secondary structures and comparison graphs.

#### Inputs
- **Protein sequence** (text or FASTA file).
- **Model selection** (`Chou-Fasman`, `GOR`, `HMM`).

#### Outputs
- Predicted secondary structures.
- Comparison graphs and character composition graphs.

### Run `auto_Validation.R`
1. Open RStudio.
2. Load `auto_Validation.R`.
3. Click **Run App**.
4. Upload a CSV file with the following columns:
   - **ProteinName**: Name of the protein.
   - **ProteinSequence**: Original sequence.
   - **DSSPSequence**: True secondary structure sequence.

#### Outputs
- **Per-protein accuracy**.
- **Average accuracy per model**.
- **Graphs for precision, recall, and F1-score**.

## Infographic for Chou-Fasman and GOR Model Mechanisms
![CF](Figures/CF.png)

**Figure 1:** A figure depicting the overall Chou Fasman algorithm (note: steps have been skipped). a. Shows the favourable helix nucleation site, b. Shows the final extended helix window, c. Shows the favourable sheet nucleation site, d. Shows the sliding window looking for more nucleation sites (in this case it was unfavourable), e. Looking for turn regions (propensity score not added for simplicity) f. Shows found favourable turn regions, g. Shows updating predicted structures onto the final prediction structure, h. Shows overlaps being handled based on higher propensity scores.

![GOR1](Figures/GOR1.png)

**Figure 2:** Calculating the information values using probabilities of finding an amino acid in a particular structural segment.


![GOR2](Figures/GOR2.png)

**Figure 3:** Stepwise calculation of structural information scores and final structure prediction.
