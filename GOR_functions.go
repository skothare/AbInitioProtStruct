// Group 2: Ab Initio Secondary Structure Prediction of Proteins
// Date: 12th December, 2024
// Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare

package main

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

/*******************FUNCTIONS*********************/

/*
	ReadGORParameters reads GOR parameters from a CSV file.

Input: Filename (string) of the .csv file.
Output: Returns a map with key strings and float64 slice. Amino acid are single-letter codes and values of floats represent the information values at positions -8 to +8 relative to the central residue.
*/
func ReadGORParameters(filename string) (InfoValTable, error) {
	// Open the CSV file
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("failed to open file %s: %v", filename, err)
	}
	defer file.Close()

	// Create a new CSV reader
	reader := csv.NewReader(file)
	reader.TrimLeadingSpace = true

	// Read the header line to skip it
	_, err = reader.Read()
	if err != nil {
		return nil, fmt.Errorf("failed to read header: %v", err)
	}

	// Initialize the parameters map
	params := make(map[string][]float64)

	// Read the rest of the CSV file
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break // End of file reached
		}
		if err != nil {
			return nil, fmt.Errorf("failed to read record: %v", err)
		}

		// The first field is the amino acid code
		aa := record[0]

		// The rest are the parameter values
		values := make([]float64, 0, 17) // There are 17 positions from -8 to +8
		for _, valStr := range record[1:] {
			val, err := strconv.ParseFloat(valStr, 64)
			if err != nil {
				return nil, fmt.Errorf("failed to parse float value %s: %v", valStr, err)
			}
			values = append(values, val)
		}

		// Store in the map
		params[aa] = values
	}

	return params, nil
}

/*
	GORPredict predicts the secondary structure of a protein sequence using the GOR method and the provided parameters for each structure type. It returns a slice of predicted structures corresponding to each residue.

Input: Protein sequence of type string, alpha-helix, beta sheet, turn, and coil information value tables.
Output: Returns a slice of predicted structures.
*/
func GORPredict(sequence string, alphaParams, betaParams, turnParams, coilParams InfoValTable) ([]GORPredictionResult, error) {

	// Set the window size for this GOR method.
	windowSize := 17 // Positions from -8 to +8

	seqLen := len(sequence) // Length of the given protein sequence
	// Create a slice of GORPredictionResult objects that is equal to the length of the sequence.
	predictions := make([]GORPredictionResult, seqLen) // Returned by this function.

	// Iterate over each residue in the sequence
	for i := 0; i < seqLen; i++ {
		// For the ith index of the corresponding residue in the sequence:

		// Calculate scores for each structure within a sliding window
		scoreAlpha, scoreBeta, scoreTurn, scoreCoil := SlideWindow(sequence, seqLen, windowSize, i, alphaParams, betaParams, turnParams, coilParams)

		// Determine the structure with the highest score
		maxScore := scoreAlpha
		structure := "H" // Helix by default

		if scoreBeta > maxScore {
			maxScore = scoreBeta
			structure = "E" // Beta-strand
		}
		if scoreTurn > maxScore {
			maxScore = scoreTurn
			structure = "T" // Turn
		}
		if scoreCoil > maxScore {
			maxScore = scoreCoil
			structure = "C" // Coil
		}
		// Store the prediction result
		predictions[i] = GORPredictionResult{
			Position:           i + 1, // Positions starting from 1
			Residue:            string(sequence[i]),
			ScoreAlpha:         scoreAlpha,
			ScoreBeta:          scoreBeta,
			ScoreTurn:          scoreTurn,
			ScoreCoil:          scoreCoil,
			PredictedStructure: structure,
		}
	}

	return predictions, nil
}

/*
	SlideWindow(): Returns four float64 values that represent four information scores for the alpha helix, beta sheet, turn, and coil structures.

Input: A string object of the entire protein sequence, three integers representing the length of the protein sequence, size of the reading window, index of the central residue of the window, and four InfoValTable (map) objects that contain the slices floats as the values to the amino acid single-letter keys.
Output: Four float64 values representing the calculated informations scores for the four structures.
*/
func SlideWindow(sequence string, seqLen, windowSize, i int, alphaParams, betaParams, turnParams, coilParams InfoValTable) (float64, float64, float64, float64) {
	halfWindow := windowSize / 2 // half the size of a given reading window.

	// Initialize scores for each structure type
	scoreAlpha, scoreBeta, scoreTurn, scoreCoil := 0.0, 0.0, 0.0, 0.0

	// Slide the reading window over the sequence
	for pos := -halfWindow; pos <= halfWindow; pos++ { // Ranges from -8 through 8 indices within the window
		windowIndex := i + pos // Calculate the current index in the window based on the central residue.

		// Check if windowIndex is within bounds
		if windowIndex < 0 || windowIndex >= seqLen {
			continue // Skip positions outside the sequence
		}

		// Get the amino acid at the window position (window position is the true index in protein seq.)
		aa := string(sequence[windowIndex])

		// Get the position index for parameters (0 to 16 indices in the slice of information values (column) in the information value table.)
		paramIndex := pos + halfWindow // Calculates the parameter index i.e. information value in a given row.

		// Get the information values for the amino acid
		alphaVals, ok := alphaParams[aa] // Finds the row for this amino acid at the index for the loop.

		var betaVals, turnVals, coilVals []float64
		// If the amino acid, "aa", does not appear in the map, use the last row's values (row X) from the InfoValTable (zero values).
		if !ok {
			alphaVals = alphaParams["X"]
			betaVals = betaParams["X"]
			turnVals = turnParams["X"]
			coilVals = coilParams["X"]

		} else { // Set other structure's values
			betaVals = betaParams[aa]
			turnVals = turnParams[aa]
			coilVals = coilParams[aa]
		}

		// Accumulate scores
		scoreAlpha += alphaVals[paramIndex]
		scoreBeta += betaVals[paramIndex]
		scoreTurn += turnVals[paramIndex]
		scoreCoil += coilVals[paramIndex]
	}

	return scoreAlpha, scoreBeta, scoreTurn, scoreCoil
}

/*
	OutputDetailedTable()

Input: A slice of the GORPredictionResult datatype that includes the structure prediction string, structure information content scores, amino acids and their corresponding positions.
Output: None. Function prints out the table to the console.
Prints a formatted table to the console, showing:
Position
Amino acid residue
Scores for alpha-helix (Ia), beta-sheet (Ib), turn (It), and coil (Ic)
The predicted structure type (St)
*/
func OutputDetailedTable(predictions []GORPredictionResult) {
	// Print the headers of the table.
	fmt.Printf("   AA    Ia    Ib    It    Ic    St\n")
	fmt.Println("        ---------------------------")
	// Print out the predictions.
	for _, pred := range predictions {
		// Determine the maximum score and mark it with an asterisk
		scores := []float64{pred.ScoreAlpha, pred.ScoreBeta, pred.ScoreTurn, pred.ScoreCoil}
		maxScore := scores[0]
		maxIndex := 0
		for idx, score := range scores {
			if score > maxScore {
				maxScore = score
				maxIndex = idx
			}
		}

		// Format scores and add asterisk to the maximum score
		scoreStrs := make([]string, 4)
		for idx, score := range scores {
			formattedScore := fmt.Sprintf("%5.0f", score)
			if idx == maxIndex {
				formattedScore += "*" // Add asterisk to the highest score
			}
			scoreStrs[idx] = formattedScore
		}

		// Print the line
		fmt.Printf("%3d %-1s %s %s %s %s   %s\n", pred.Position, pred.Residue,
			scoreStrs[0], scoreStrs[1], scoreStrs[2], scoreStrs[3], pred.PredictedStructure)
	}

	// Output the predicted secondary structure as a string
	fmt.Print("Predicted GOR secondary structure: ")
	for _, pred := range predictions {
		fmt.Print(pred.PredictedStructure)
	}
	fmt.Println()
}

/*
OutputGORSequence(): This is an output visualization function.
Input: Slice containing the GORPredictionResult objects. OutputGORSequence takes as input the same as above.
Output: Returns the sequence, does not include table.
*/
func OutputGORSequence(predictions []GORPredictionResult) string {
	// Print the headers of the table.
	//fmt.Printf("   AA    Ia    Ib    It    Ic    St\n")
	//fmt.Println("        ---------------------------")
	// Print out the predictions.
	for _, pred := range predictions {
		// Determine the maximum score and mark it with an asterisk
		scores := []float64{pred.ScoreAlpha, pred.ScoreBeta, pred.ScoreTurn, pred.ScoreCoil}
		maxScore := scores[0]
		maxIndex := 0
		for idx, score := range scores {
			if score > maxScore {
				maxScore = score
				maxIndex = idx
			}
		}

		// Format scores and add asterisk to the maximum score
		scoreStrs := make([]string, 4)
		for idx, score := range scores {
			formattedScore := fmt.Sprintf("%5.0f", score)
			if idx == maxIndex {
				formattedScore += "*" // Add asterisk to the highest score
			}
			scoreStrs[idx] = formattedScore
		}
	}

	result := make([]string, len(predictions))
	for i, pred := range predictions {
		result[i] = pred.PredictedStructure
	}

	resultString := strings.Join(result, "")
	return resultString
}
