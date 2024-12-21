// Group 2: Ab Initio Secondary Structure Prediction of Proteins
// Date: 12th December, 2024
// Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare
package main

import (
	"reflect"
	"testing"
)

func TestNewHMM(t *testing.T) {
	states := []string{"Helix", "Sheet", "Coil", "Turn"}
	symbols := []string{"A", "C", "D"}

	hmm := NewHMM(states, symbols)

	tests := []struct {
		name            string
		actualStates    []string
		expectedStates  []string
		actualSymbols   []string
		expectedSymbols []string
	}{
		{
			name:           "Check states initialization",
			actualStates:   hmm.States,
			expectedStates: states,
		},
		{
			name:            "Check symbols initialization",
			actualSymbols:   hmm.Symbols,
			expectedSymbols: symbols,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if !reflect.DeepEqual(tt.actualStates, tt.expectedStates) && !reflect.DeepEqual(tt.actualSymbols, tt.expectedSymbols) {
				t.Errorf("Failed %s: got %v and %v, expected %v and %v", tt.name, tt.actualStates, tt.actualSymbols, tt.expectedStates, tt.expectedSymbols)
			}
		})
	}
}

func TestViterbi(t *testing.T) {
	hmm := NewHMM(
		[]string{"Helix", "Sheet", "Coil"},
		[]string{"A", "C", "G", "T"},
	)

	hmm.Transition = [][]float64{
		{0.6, 0.3, 0.1},
		{0.1, 0.7, 0.2},
		{0.3, 0.3, 0.4},
	}
	hmm.Emission = [][]float64{
		{0.2, 0.4, 0.3, 0.1},
		{0.1, 0.3, 0.4, 0.2},
		{0.4, 0.1, 0.2, 0.3},
	}
	hmm.Initial = []float64{0.5, 0.3, 0.2}

	tests := []struct {
		name         string
		sequence     string
		expectedPath string
	}{
		{
			name:         "Simple observation sequence",
			sequence:     "ACG",
			expectedPath: "HHH",
		},
		{
			name:         "Repeated observation sequence",
			sequence:     "AAA",
			expectedPath: "HHH",
		},
		{
			name:         "Mixed observation sequence",
			sequence:     "GTA",
			expectedPath: "SSS",
		},
		{
			name:         "Complex observation sequence",
			sequence:     "TCGA",
			expectedPath: "SSSS",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			resultPath := hmm.Viterbi(tt.sequence)
			if resultPath != tt.expectedPath {
				t.Errorf("%s: got %v, expected %v", tt.name, resultPath, tt.expectedPath)
			}
		})
	}
}
