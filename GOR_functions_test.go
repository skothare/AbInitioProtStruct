// Group 2: Ab Initio Secondary Structure Prediction of Proteins
// Date: 12th December, 2024
// Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare
// GOR_functions_test.go
package main

import (
	"bufio"
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

/******************************************************************************************************
TEST SLIDEWINDOW()
******************************************************************************************************/

// parseParamTable parses parameter lines into an InfoValTable.
// It skips any lines that start with "Position".
func parseParamTable(lines []string) (InfoValTable, error) {
	params := make(InfoValTable)
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line == "" || strings.HasPrefix(line, "#") || strings.HasPrefix(line, "Position") {
			continue // Skip empty lines, comments, and Position headers
		}
		parts := strings.Split(line, ",")
		if len(parts) != 18 { // 1 AA code + 17 values
			return nil, fmt.Errorf("invalid parameter line: %s", line)
		}
		aa := parts[0]
		values := make([]float64, 17)
		for i, valStr := range parts[1:] {
			val, err := strconv.ParseFloat(valStr, 64)
			if err != nil {
				return nil, fmt.Errorf("invalid float value '%s' in AA '%s': %v", valStr, aa, err)
			}
			values[i] = val
		}
		params[aa] = values
	}
	return params, nil
}

// readSlideWindowTestCase reads and parses a SlideWindow test case from input_#.txt
func readSlideWindowTestCase(t *testing.T, testCase int) (sequence string, windowSize int, index int,
	alphaParams, betaParams, turnParams, coilParams InfoValTable) {

	inputPath := filepath.Join("GORTests/SlideWindow", "Input", fmt.Sprintf("input_%d.txt", testCase))
	file, err := os.Open(inputPath)
	if err != nil {
		t.Fatalf("Failed to open input file %d: %v", testCase, err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Initialize slices to hold parameter lines
	var alphaLines, betaLines, turnLines, coilLines []string
	sequence = ""
	windowSize = 17 // default, will override if specified
	index = 0       // default, will override if specified

	currentSection := ""

	for scanner.Scan() {
		line := scanner.Text()
		line = strings.TrimSpace(line)

		if strings.HasPrefix(line, "#") {
			// Identify the current section
			if strings.Contains(line, "Alpha Parameters") {
				currentSection = "alpha"
				continue
			} else if strings.Contains(line, "Beta Parameters") {
				currentSection = "beta"
				continue
			} else if strings.Contains(line, "Turn Parameters") {
				currentSection = "turn"
				continue
			} else if strings.Contains(line, "Coil Parameters") {
				currentSection = "coil"
				continue
			} else if strings.Contains(line, "Sequence") {
				currentSection = "sequence"
				continue
			} else if strings.Contains(line, "WindowSize") {
				currentSection = "windowSize"
				continue
			} else if strings.Contains(line, "Index") {
				currentSection = "index"
				continue
			} else {
				currentSection = ""
				continue
			}
		}

		switch currentSection {
		case "alpha":
			alphaLines = append(alphaLines, line)
		case "beta":
			betaLines = append(betaLines, line)
		case "turn":
			turnLines = append(turnLines, line)
		case "coil":
			coilLines = append(coilLines, line)
		case "sequence":
			sequence = line
		case "windowSize":
			val, err := strconv.Atoi(line)
			if err != nil {
				t.Fatalf("Invalid WindowSize in input file %d: %v", testCase, err)
			}
			windowSize = val
		case "index":
			val, err := strconv.Atoi(line)
			if err != nil {
				t.Fatalf("Invalid Index in input file %d: %v", testCase, err)
			}
			index = val
		}
	}

	if err := scanner.Err(); err != nil {
		t.Fatalf("Error reading input file %d: %v", testCase, err)
	}

	// Parse parameter tables
	alphaParams, err = parseParamTable(alphaLines)
	if err != nil {
		t.Fatalf("Error parsing Alpha Parameters in test case %d: %v", testCase, err)
	}

	betaParams, err = parseParamTable(betaLines)
	if err != nil {
		t.Fatalf("Error parsing Beta Parameters in test case %d: %v", testCase, err)
	}

	turnParams, err = parseParamTable(turnLines)
	if err != nil {
		t.Fatalf("Error parsing Turn Parameters in test case %d: %v", testCase, err)
	}

	coilParams, err = parseParamTable(coilLines)
	if err != nil {
		t.Fatalf("Error parsing Coil Parameters in test case %d: %v", testCase, err)
	}

	return
}

// readSlideWindowExpectedOutput reads expected scores from output_#.txt
func readSlideWindowExpectedOutput(t *testing.T, testCase int) (expectedAlpha, expectedBeta, expectedTurn, expectedCoil float64) {
	outputPath := filepath.Join("GORTests/SlideWindow", "Output", fmt.Sprintf("output_%d.txt", testCase))
	file, err := os.Open(outputPath)
	if err != nil {
		t.Fatalf("Failed to open output file %d: %v", testCase, err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()
		line = strings.TrimSpace(line)
		if strings.HasPrefix(line, "scoreAlpha:") {
			valStr := strings.TrimSpace(strings.TrimPrefix(line, "scoreAlpha:"))
			val, err := strconv.ParseFloat(valStr, 64)
			if err != nil {
				t.Fatalf("Invalid scoreAlpha value in output file %d: %v", testCase, err)
			}
			expectedAlpha = val
		} else if strings.HasPrefix(line, "scoreBeta:") {
			valStr := strings.TrimSpace(strings.TrimPrefix(line, "scoreBeta:"))
			val, err := strconv.ParseFloat(valStr, 64)
			if err != nil {
				t.Fatalf("Invalid scoreBeta value in output file %d: %v", testCase, err)
			}
			expectedBeta = val
		} else if strings.HasPrefix(line, "scoreTurn:") {
			valStr := strings.TrimSpace(strings.TrimPrefix(line, "scoreTurn:"))
			val, err := strconv.ParseFloat(valStr, 64)
			if err != nil {
				t.Fatalf("Invalid scoreTurn value in output file %d: %v", testCase, err)
			}
			expectedTurn = val
		} else if strings.HasPrefix(line, "scoreCoil:") {
			valStr := strings.TrimSpace(strings.TrimPrefix(line, "scoreCoil:"))
			val, err := strconv.ParseFloat(valStr, 64)
			if err != nil {
				t.Fatalf("Invalid scoreCoil value in output file %d: %v", testCase, err)
			}
			expectedCoil = val
		}
	}

	if err := scanner.Err(); err != nil {
		t.Fatalf("Error reading output file %d: %v", testCase, err)
	}

	return
}

/*
	floatEquals checks if two float64 numbers are equal within a small tolerance.
	Input: a, b (float64), epsilon (float64)
	Output: bool
*/

func floatEquals(a, b, epsilon float64) bool {
	return (a-b) < epsilon && (b-a) < epsilon
}

// TestSlideWindow tests the SlideWindow function using predefined test cases
func TestSlideWindow(t *testing.T) {
	// Define the number of test cases
	numTestCases := 5

	for testCase := 1; testCase <= numTestCases; testCase++ {
		// Read and parse the input file
		sequence, windowSize, index, alphaParams, betaParams, turnParams, coilParams := readSlideWindowTestCase(t, testCase)

		// Read and parse the expected output
		expectedAlpha, expectedBeta, expectedTurn, expectedCoil := readSlideWindowExpectedOutput(t, testCase)

		// Call the SlideWindow function
		actualAlpha, actualBeta, actualTurn, actualCoil := SlideWindow(sequence, len(sequence), windowSize, index, alphaParams, betaParams, turnParams, coilParams)

		// Define a small epsilon for floating-point comparison
		epsilon := 1e-6

		// Compare the actual scores with expected scores
		if !floatEquals(actualAlpha, expectedAlpha, epsilon) {
			t.Errorf("Test case %d: scoreAlpha mismatch. Expected: %.6f, Got: %.6f", testCase, expectedAlpha, actualAlpha)
		}
		if !floatEquals(actualBeta, expectedBeta, epsilon) {
			t.Errorf("Test case %d: scoreBeta mismatch. Expected: %.6f, Got: %.6f", testCase, expectedBeta, actualBeta)
		}
		if !floatEquals(actualTurn, expectedTurn, epsilon) {
			t.Errorf("Test case %d: scoreTurn mismatch. Expected: %.6f, Got: %.6f", testCase, expectedTurn, actualTurn)
		}
		if !floatEquals(actualCoil, expectedCoil, epsilon) {
			t.Errorf("Test case %d: scoreCoil mismatch. Expected: %.6f, Got: %.6f", testCase, expectedCoil, actualCoil)
		}
	}
}

/******************************************************************************************************
TEST GORPREDICT()
******************************************************************************************************/

/*
parseInputFile parses the input_#.txt file and returns the sequence and parameter tables.
Input: filePath (string) of the input file.
Output: sequence (string), alphaParams, betaParams, turnParams, coilParams (InfoValTable), error
*/
func parseInputFile(filePath string) (string, InfoValTable, InfoValTable, InfoValTable, InfoValTable, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return "", nil, nil, nil, nil, fmt.Errorf("failed to open input file %s: %v", filePath, err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var alphaParams, betaParams, turnParams, coilParams InfoValTable
	var sequence string
	var currentSection string
	var sequenceCount int

	// Buffer to store lines of current section
	var buffer []string

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		if line == "" {
			continue // Skip empty lines
		}

		if strings.HasPrefix(line, "#") {
			// Process previous section
			if currentSection != "" && len(buffer) > 0 {
				if currentSection == "Sequence" {
					if sequenceCount >= 1 {
						return "", nil, nil, nil, nil, errors.New("multiple sequences detected in input file")
					}
					sequence = strings.TrimSpace(buffer[0])
					sequenceCount++
				} else {
					reader := strings.NewReader(strings.Join(buffer, "\n"))
					table, err := parseParameterTable(reader)
					if err != nil {
						return "", nil, nil, nil, nil, fmt.Errorf("error parsing %s parameters in file %s: %v", currentSection, filePath, err)
					}
					switch currentSection {
					case "Alpha":
						alphaParams = table
					case "Beta":
						betaParams = table
					case "Turn":
						turnParams = table
					case "Coil":
						coilParams = table
					}
				}
			}

			// Identify new section
			if strings.HasPrefix(line, "# Alpha Parameters") {
				currentSection = "Alpha"
			} else if strings.HasPrefix(line, "# Beta Parameters") {
				currentSection = "Beta"
			} else if strings.HasPrefix(line, "# Turn Parameters") {
				currentSection = "Turn"
			} else if strings.HasPrefix(line, "# Coil Parameters") {
				currentSection = "Coil"
			} else if strings.HasPrefix(line, "# Sequence") {
				currentSection = "Sequence"
			} else {
				currentSection = "" // Unknown or irrelevant section
			}

			buffer = []string{} // Reset buffer
			continue
		}

		if currentSection != "" {
			buffer = append(buffer, line)
		}
	}

	// Process the last section
	if currentSection != "" && len(buffer) > 0 {
		if currentSection == "Sequence" {
			if sequenceCount >= 1 {
				return "", nil, nil, nil, nil, errors.New("multiple sequences detected in input file")
			}
			sequence = strings.TrimSpace(buffer[0])
			sequenceCount++
		} else {
			reader := strings.NewReader(strings.Join(buffer, "\n"))
			table, err := parseParameterTable(reader)
			if err != nil {
				return "", nil, nil, nil, nil, fmt.Errorf("error parsing %s parameters in file %s: %v", currentSection, filePath, err)
			}
			switch currentSection {
			case "Alpha":
				alphaParams = table
			case "Beta":
				betaParams = table
			case "Turn":
				turnParams = table
			case "Coil":
				coilParams = table
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return "", nil, nil, nil, nil, fmt.Errorf("error scanning input file %s: %v", filePath, err)
	}

	// Validate that all parameter tables and sequence are present
	if alphaParams == nil {
		return "", nil, nil, nil, nil, errors.New("alpha parameters missing in input file")
	}
	if betaParams == nil {
		return "", nil, nil, nil, nil, errors.New("beta parameters missing in input file")
	}
	if turnParams == nil {
		return "", nil, nil, nil, nil, errors.New("turn parameters missing in input file")
	}
	if coilParams == nil {
		return "", nil, nil, nil, nil, errors.New("coil parameters missing in input file")
	}
	if sequence == "" {
		return "", nil, nil, nil, nil, errors.New("sequence missing in input file")
	}

	return sequence, alphaParams, betaParams, turnParams, coilParams, nil
}

// parseParameterTable parses a parameter table section into an InfoValTable.
func parseParameterTable(reader io.Reader) (InfoValTable, error) {
	csvReader := csv.NewReader(reader)
	csvReader.TrimLeadingSpace = true

	// Read the header line
	header, err := csvReader.Read()
	if err != nil {
		return nil, fmt.Errorf("failed to read header: %v", err)
	}

	if len(header) < 2 {
		return nil, errors.New("invalid header in parameter table")
	}

	positions := header[1:]
	expectedWindowSize := len(positions)
	if expectedWindowSize != 17 { // For window size 17 (-8 to +8)
		return nil, fmt.Errorf("expected 17 positions, got %d", expectedWindowSize)
	}

	// Parse positions to ensure they are -8 to +8
	for i, pos := range positions {
		expectedPos := strconv.Itoa(-8 + i)
		if pos != expectedPos {
			return nil, fmt.Errorf("expected position %d, got %s", -8+i, pos)
		}
	}

	table := make(InfoValTable)

	// Read each amino acid line
	for {
		record, err := csvReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("error reading parameter table: %v", err)
		}

		if len(record) != 18 { // Amino acid + 17 scores
			return nil, fmt.Errorf("invalid record length: expected 18, got %d", len(record))
		}

		aa := record[0]
		scores := make([]float64, 17)
		for i := 1; i < 18; i++ {
			score, err := strconv.ParseFloat(record[i], 64)
			if err != nil {
				return nil, fmt.Errorf("invalid score for amino acid %s at position %s: %v", aa, header[i], err)
			}
			scores[i-1] = score
		}
		table[aa] = scores
	}

	return table, nil
}

/*
OutputGORSequence2 returns the predicted secondary structure as a single string.
Input: Slice of GORPredictionResult
Output: structureString (string)
*/
func OutputGORSequence2(predictions []GORPredictionResult) string {
	result := make([]string, len(predictions))
	for i, pred := range predictions {
		result[i] = pred.PredictedStructure
	}

	resultString := strings.Join(result, "")
	return resultString
}

/*******************GORPredict TESTING FUNCTIONS*********************/

/*
parseExpectedOutput parses the output_#.txt file and returns the expected structure string.
Input: filePath (string) of the output file.
Output: expectedStructure (string), error
*/
func parseExpectedOutput(filePath string) (string, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return "", fmt.Errorf("failed to open output file %s: %v", filePath, err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var expectedStructure string

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		// Assuming the output file contains only the structure string
		expectedStructure = line
		break // Read only the first non-empty line
	}

	if err := scanner.Err(); err != nil {
		return "", fmt.Errorf("error reading output file %s: %v", filePath, err)
	}

	if expectedStructure == "" {
		return "", fmt.Errorf("no structure string found in output file %s", filePath)
	}

	return expectedStructure, nil
}

/*
TestGORPredict tests the GORPredict function by comparing the predicted structure string against the expected structure string from output_#.txt files.
*/
func TestGORPredict(t *testing.T) {
	inputDir := "GORTests/GORPredict/Input"
	outputDir := "GORTests/GORPredict/Output"

	// Determine the number of test cases by counting input files
	files, err := os.ReadDir(inputDir)
	if err != nil {
		t.Fatalf("Failed to read input directory: %v", err)
	}

	// Filter input files matching the pattern input_#.txt
	var inputFiles []os.DirEntry
	for _, file := range files {
		if strings.HasPrefix(file.Name(), "input_") && strings.HasSuffix(file.Name(), ".txt") {
			inputFiles = append(inputFiles, file)
		}
	}

	if len(inputFiles) == 0 {
		t.Fatalf("No input files found in %s", inputDir)
	}

	for _, inputFile := range inputFiles {
		testCaseName := strings.TrimSuffix(inputFile.Name(), ".txt")
		testCaseNumber := strings.TrimPrefix(testCaseName, "input_")
		outputFileName := fmt.Sprintf("output_%s.txt", testCaseNumber)
		outputFilePath := filepath.Join(outputDir, outputFileName)

		inputFilePath := filepath.Join(inputDir, inputFile.Name())

		// Parse input file
		sequence, alphaParams, betaParams, turnParams, coilParams, err := parseInputFile(inputFilePath)
		if err != nil {
			t.Errorf("Test case %s: failed to parse input file: %v", testCaseNumber, err)
			continue
		}

		// Parse expected output
		expectedStructure, err := parseExpectedOutput(outputFilePath)
		if err != nil {
			t.Errorf("Test case %s: failed to parse output file: %v", testCaseNumber, err)
			continue
		}

		// Run GORPredict
		predictions, err := GORPredict(sequence, alphaParams, betaParams, turnParams, coilParams)
		if err != nil {
			t.Errorf("Test case %s: GORPredict returned error: %v", testCaseNumber, err)
			continue
		}

		// Get the predicted structure string
		predictedStructure := OutputGORSequence2(predictions)

		// Compare predicted structure with expected structure
		if predictedStructure != expectedStructure {
			t.Errorf("Test case %s: structure mismatch.\nExpected: %s\nGot:      %s", testCaseNumber, expectedStructure, predictedStructure)
		} else {
			t.Logf("Test case %s: PASS", testCaseNumber)
		}
	}
}
