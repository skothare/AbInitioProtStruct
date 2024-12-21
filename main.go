// Group 2: Ab Initio Secondary Structure Prediction of Proteins
// Date: 12th December, 2024
// Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare
package main

import (
	"fmt"
	"os"
	"strings"
)

// Main function to execute the secondary structure prediction
func main() {
	/** GOR & CF **/
	// Load GOR parameter files
	alphaParams, err := ReadGORParameters("GOR_InfoVals/InfoVal_aHelix.csv")
	if err != nil {
		fmt.Printf("Error reading alpha parameters: %v\n", err)
		return
	}

	betaParams, err := ReadGORParameters("GOR_InfoVals/InfoVal_bStrand.csv")
	if err != nil {
		fmt.Printf("Error reading beta parameters: %v\n", err)
		return
	}

	turnParams, err := ReadGORParameters("GOR_InfoVals/InfoVal_bTurn.csv")
	if err != nil {
		fmt.Printf("Error reading turn parameters: %v\n", err)
		return
	}

	coilParams, err := ReadGORParameters("GOR_InfoVals/InfoVal_Coil.csv")
	if err != nil {
		fmt.Printf("Error reading coil parameters: %v\n", err)
		return
	}

	// Check for input sequence argument
	if len(os.Args) < 2 {
		fmt.Println("Please provide a string argument")
		return
	}
	input := os.Args[1]

	// Convert input sequence to uppercase
	sequence := strings.ToUpper(input)

	if sequence == "" {
		fmt.Println("No sequence provided.")
		return
	}

	// Validate the input sequence for valid amino acids
	if !isValidSequence(sequence) {
		fmt.Println("Invalid sequence: sequence must only contain valid amino acid codes (A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V)")
		return
	}

	runeSequence := []rune(sequence) // Convert string to runes for CF prediction

	// Predict the secondary structure using GOR method
	gorPredictions, err := GORPredict(sequence, alphaParams, betaParams, turnParams, coilParams)
	if err != nil {
		fmt.Printf("Error in GOR prediction: %v\n", err)
		return
	}

	/** HMM **/
	// Define states and symbols for HMM
	states := []string{"Helix", "ESheet", "Coil", "Turn"}
	symbols := []string{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}

	hmm := NewHMM(states, symbols) // Create the HMM object

	hmm.Initial = []float64{0.29, 0.29, 0.33, 0.09} // Initial probabilities for each state

	hmm.Transition = [][]float64{ // Transition probabilities between states
		{0.4, 0.3, 0.2, 0.1},     // Helix transitions
		{0.3, 0.39, 0.21, 0.1},   // Sheet transitions
		{0.25, 0.2, 0.4, 0.15},   // Coil transitions
		{0.22, 0.22, 0.22, 0.34}, // Turn transitions
	}

	hmm.Emission = [][]float64{ // Emission probabilities for each state and symbol
		// Helix emissions
		{0.09600544711756695, 0.010667271901951884, 0.04312301407172038, 0.08556513844757149,
			0.04153427144802542, 0.03313663186563777, 0.01702224239673173, 0.06672719019518839,
			0.06967771221062188, 0.12210621879255561, 0.02973218338629142, 0.028370403994552883,
			0.02292328642759873, 0.060372219700408535, 0.056513844757149344, 0.061053109396277803,
			0.04312301407172038, 0.06740807989105765, 0.01702224239673173, 0.027916477530640037},
		// Sheet emissions
		{0.061148086522462564, 0.018718801996672214, 0.024126455906821963, 0.04492512479201331,
			0.05698835274542429, 0.04159733777038269, 0.02454242928452579, 0.08527454242928452,
			0.04076539101497504, 0.11314475873544093, 0.0262063227953411, 0.022462562396006656,
			0.014143094841930116, 0.03410981697171381, 0.042429284525790346, 0.044509151414309486,
			0.08153078202995008, 0.14101497504159735, 0.025374376039933443, 0.05698835274542429},
		// Coil emissions
		{0.06551990722844994, 0.014302280633938926, 0.08233475067645922, 0.054310011596443754,
			0.029764205643602628, 0.08909934286818709, 0.03710862002319289, 0.031310398144569,
			0.04580595284112872, 0.06378044066486277, 0.020487050637804406, 0.057209122535755705,
			0.08658678005411674, 0.03691534596057209, 0.040780827212988015, 0.09161190568225744,
			0.06358716660224198, 0.049671434093544645, 0.01256281407035176, 0.027251642829532276},

		// Turn emissions
		{0.054404145077720206, 0.007772020725388601, 0.10233160621761658, 0.06347150259067358,
			0.012953367875647668, 0.20725388601036268, 0.019430051813471502, 0.006476683937823834,
			0.05569948186528497, 0.03238341968911917, 0.011658031088082901, 0.10880829015544041,
			0.09844559585492228, 0.04404145077720207, 0.05181347150259067, 0.06347150259067358,
			0.031088082901554404, 0.015544041450777202, 0.0025906735751295338, 0.010362694300518135},
	}

	hmmOutput := hmm.Viterbi(sequence) // Predict using Viterbi algorithm

	// Print results of predictions from different methods
	fmt.Printf("Input sequence: %s\n", sequence)
	fmt.Printf("Predicted Chou-Fasman secondary structure: %s\n", ChouFasmanPredictSS(runeSequence))
	fmt.Printf("Predicted GOR secondary structure: %s\n", OutputGORSequence(gorPredictions))
	fmt.Printf("Predicted HMM secondary structure: %s\n", hmmOutput) // H for Helix etc.
}

// Function to validate if the sequence contains valid amino acids only
func isValidSequence(sequence string) bool {
	for _, aa := range sequence {
		_, exists := propensities[aa] // Check if amino acid is valid based on propensities map
		if !exists {
			return false
		}
	}
	return true
}

/*package main

import (
	"fmt"
	"os"
	"strings"
)

func main() {
	// Read the parameter files for GOR
	alphaParams, err := ReadGORParameters("GOR_InfoVals/InfoVal_aHelix.csv")
	if err != nil {
		fmt.Printf("Error reading alpha parameters: %v\n", err)
		return
	}

	betaParams, err := ReadGORParameters("GOR_InfoVals/InfoVal_bStrand.csv")
	if err != nil {
		fmt.Printf("Error reading beta parameters: %v\n", err)
		return
	}

	turnParams, err := ReadGORParameters("GOR_InfoVals/InfoVal_bTurn.csv")
	if err != nil {
		fmt.Printf("Error reading turn parameters: %v\n", err)
		return
	}

	coilParams, err := ReadGORParameters("GOR_InfoVals/InfoVal_Coil.csv")
	if err != nil {
		fmt.Printf("Error reading coil parameters: %v\n", err)
		return
	}

	if len(os.Args) < 2 {
		fmt.Println("Please provide a string argument")
		return
	}
	input := os.Args[1]

	sequence := strings.ToUpper(input)

	if sequence == "" {
		fmt.Println("No sequence provided.")
		return
	}

	// Validate sequence
	if !isValidSequence(sequence) {
		fmt.Println("Invalid sequence: sequence must only contain valid amino acid codes (A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V)")
		return
	}

	// Convert string to runes for CF prediction
	runeSequence := []rune(sequence)

	// Predict the secondary structure using GOR
	gorPredictions, err := GORPredict(sequence, alphaParams, betaParams, turnParams, coilParams)
	if err != nil {
		fmt.Printf("Error in GOR prediction: %v\n", err)
		return
	}

	// Define states (secondary structures) and symbols (amino acids) for HMM
	states := []string{"Helix", "ESheet", "Coil", "Turn"}
	symbols := []string{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}

	// Create the HMM object
	hmm := NewHMM(states, symbols)

	// Manually set initial probabilities for each state
	hmm.Initial = []float64{0.29, 0.29, 0.33, 0.09} // Example probabilities

	// Manually set transition probabilities between states
	hmm.Transition = [][]float64{
		{0.4, 0.3, 0.2, 0.1},     // Helix transitions
		{0.3, 0.39, 0.21, 0.1},   // Sheet transitions
		{0.25, 0.2, 0.4, 0.15},   // Coil transitions
		{0.22, 0.22, 0.22, 0.34}, // Turn transitions
	}

	// Manually set emission probabilities (likelihood of symbols given each state)
	hmm.Emission = [][]float64{
		// Helix emissions
		{0.09600544711756695, 0.010667271901951884, 0.04312301407172038, 0.08556513844757149,
			0.04153427144802542, 0.03313663186563777, 0.01702224239673173, 0.06672719019518839,
			0.06967771221062188, 0.12210621879255561, 0.02973218338629142, 0.028370403994552883,
			0.02292328642759873, 0.060372219700408535, 0.056513844757149344, 0.061053109396277803,
			0.04312301407172038, 0.06740807989105765, 0.01702224239673173, 0.027916477530640037},
		// Sheet emissions
		{0.061148086522462564, 0.018718801996672214, 0.024126455906821963, 0.04492512479201331,
			0.05698835274542429, 0.04159733777038269, 0.02454242928452579, 0.08527454242928452,
			0.04076539101497504, 0.11314475873544093, 0.0262063227953411, 0.022462562396006656,
			0.014143094841930116, 0.03410981697171381, 0.042429284525790346, 0.044509151414309486,
			0.08153078202995008, 0.14101497504159735, 0.025374376039933443, 0.05698835274542429},
		// Coil emissions
		{0.06551990722844994, 0.014302280633938926, 0.08233475067645922, 0.054310011596443754,
			0.029764205643602628, 0.08909934286818709, 0.03710862002319289, 0.031310398144569,
			0.04580595284112872, 0.06378044066486277, 0.020487050637804406, 0.057209122535755705,
			0.08658678005411674, 0.03691534596057209, 0.040780827212988015, 0.09161190568225744,
			0.06358716660224198, 0.049671434093544645, 0.01256281407035176, 0.027251642829532276},

		// Turn emissions
		{0.054404145077720206, 0.007772020725388601, 0.10233160621761658, 0.06347150259067358,
			0.012953367875647668, 0.20725388601036268, 0.019430051813471502, 0.006476683937823834,
			0.05569948186528497, 0.03238341968911917, 0.011658031088082901, 0.10880829015544041,
			0.09844559585492228, 0.04404145077720207, 0.05181347150259067, 0.06347150259067358,
			0.031088082901554404, 0.015544041450777202, 0.0025906735751295338, 0.010362694300518135},
	}

	// Predict the secondary structure using the Viterbi algorithm for HMM
	hmmOutput := hmm.Viterbi(sequence)

	// Print the input and the predicted outputs
	fmt.Printf("Input sequence: %s\n", sequence)
	fmt.Printf("Predicted Chou-Fasman secondary structure: %s\n", ChouFasmanPredictSS(runeSequence))
	fmt.Printf("Predicted GOR secondary structure: %s\n", OutputGORSequence(gorPredictions))
	fmt.Printf("Predicted HMM secondary structure: %s\n", hmmOutput) // H for Helix, E for Sheet, etc.
}

func isValidSequence(sequence string) bool {
	for _, aa := range sequence {
		// Check if amino acid exists in our propensities map
		_, exists := propensities[aa]
		if !exists {
			return false
		}
	}
	return true
}
*/
