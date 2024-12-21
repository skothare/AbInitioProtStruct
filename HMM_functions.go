// Group 2: Ab Initio Secondary Structure Prediction of Proteins
// Date: 12th December, 2024
// Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare

package main

import (
	"log"
)

// NewHMM initializes an HMM with the given states and symbols.
// It prepares empty probability matrices and mappings for states and symbols.
func NewHMM(states, symbols []string) *HMM {
	stateMapping := make(map[string]int)
	for i, state := range states {
		stateMapping[state] = i
	}
	symbolMapping := make(map[string]int)
	for i, symbol := range symbols {
		symbolMapping[symbol] = i
	}
	numStates := len(states)
	//numSymbols := len(symbols)

	return &HMM{
		States:        states,
		Symbols:       symbols,
		Transition:    make([][]float64, numStates), // Empty matrix for state transitions
		Emission:      make([][]float64, numStates), // Empty matrix for emissions
		Initial:       make([]float64, numStates),   // Empty array for initial probabilities
		StateMapping:  stateMapping,                 // Map state names to indices
		SymbolMapping: symbolMapping,                // Map symbol names to indices
	}
}

// Viterbi computes the most likely sequence of states for a given observation sequence.
// It uses dynamic programming to find the optimal path through the HMM.
func (hmm *HMM) Viterbi(sequence string) string {
	T := len(sequence)   // Length of the observation sequence
	N := len(hmm.States) // Number of states in the HMM

	// v[t][j] holds the highest probability of any path that ends in state j at time t.
	v := make([][]float64, T)
	// backpointer[t][j] stores the previous state that led to state j at time t.
	backpointer := make([][]int, T)

	// Initialize the v and backpointer arrays for each time step.
	for t := 0; t < T; t++ {
		v[t] = make([]float64, N)
		backpointer[t] = make([]int, N)
	}

	// Initialization step: Compute probabilities for the first observation.
	for i := 0; i < N; i++ {
		symbol := string(sequence[0]) // First symbol in the sequence
		emissionIdx, ok := hmm.SymbolMapping[symbol]
		if !ok {
			log.Fatalf("Invalid symbol: %s", symbol) // Handle invalid symbols
		}
		v[0][i] = hmm.Initial[i] + hmm.Emission[i][emissionIdx] // Log probabilities
		backpointer[0][i] = -1                                  // No previous state at the first step
	}

	// Recursion step: Compute probabilities for subsequent observations.
	for t := 1; t < T; t++ {
		for j := 0; j < N; j++ {
			maxVal := -1e9 // Initialize to a very low probability
			maxIdx := -1   // Index of the best previous state
			for i := 0; i < N; i++ {
				val := v[t-1][i] + hmm.Transition[i][j] // Transition from state i to state j
				if val > maxVal {                       // Keep track of the best probability
					maxVal = val
					maxIdx = i
				}
			}
			symbol := string(sequence[t]) // Current symbol in the sequence
			emissionIdx, ok := hmm.SymbolMapping[symbol]
			if !ok {
				log.Fatalf("Invalid symbol: %s", symbol) // Handle invalid symbols
			}
			v[t][j] = maxVal + hmm.Emission[j][emissionIdx] // Update probability for state j
			backpointer[t][j] = maxIdx                      // Record the state that led to state j
		}
	}

	// Termination step: Backtrack to find the most likely sequence of states.
	bestPath := make([]int, T) // Array to store the optimal path
	lastState := 0             // Index of the last state in the best path
	maxVal := -1e9             // Maximum probability of any path ending at the last time step
	for i := 0; i < N; i++ {
		if v[T-1][i] > maxVal {
			maxVal = v[T-1][i]
			lastState = i
		}
	}
	bestPath[T-1] = lastState // Start backtracking from the last state

	// Backtrack through the backpointer array to reconstruct the path.
	for t := T - 2; t >= 0; t-- {
		bestPath[t] = backpointer[t+1][bestPath[t+1]]
	}

	// Convert the sequence of state indices to state labels (e.g., H, E, T, C).
	result := ""
	for _, stateIdx := range bestPath {
		result += string(hmm.States[stateIdx][0]) // Use the first letter of each state
	}

	return result // Return the predicted sequence of states
}
