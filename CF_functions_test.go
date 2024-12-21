// Group 2: Ab Initio Secondary Structure Prediction of Proteins
// Date: 12th December, 2024
// Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare
package main

import (
	"math"
	"reflect"
	"testing"
)

func TestChouFasmanPredictSS(t *testing.T) {
	// Define test cases
	tests := []struct {
		name     string // Name of the test case
		sequence []rune // Input sequence
		expected string // Expected secondary structure prediction
	}{
		{
			name:     "Beta Sheet and Coil",
			sequence: []rune{'V', 'V', 'V', 'I', 'I', 'I', 'E', 'E', 'E', 'L', 'L', 'L', 'P', 'G', 'G'},
			expected: "EEEEEHHHHHCTTTT",
		},
		{
			name:     "Helix Region",
			sequence: []rune{'A', 'A', 'A', 'A', 'A', 'A'},
			expected: "HHHHHH", // Alanine forms a strong helix
		},
		{
			name:     "Turn Region",
			sequence: []rune{'S', 'P', 'N', 'D'},
			expected: "TTTT",
		},
		{
			name:     "Mixed Helix and Coil",
			sequence: []rune{'E', 'E', 'K', 'K', 'L'},
			expected: "CCCCC",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := ChouFasmanPredictSS(tt.sequence)
			if result != tt.expected {
				t.Errorf("Test %s failed. Expected %s but got %s", tt.name, tt.expected, result)
			}
		})
	}
}

func normalizeRegions(regions []Region) []Region {
	if regions == nil {
		return []Region{}
	}
	return regions
}

func TestPredictHelix(t *testing.T) {
	tests := []struct {
		name     string   // Name of the test case
		sequence []rune   // Input sequence
		expected []Region // Expected output regions
	}{
		{
			name:     "Strong Helix",
			sequence: []rune{'A', 'A', 'A', 'A', 'A', 'A'},
			expected: []Region{
				{start: 0, end: 6, score: 1.42, structure: 'H'},
			},
		},
		{
			name:     "Mixed Helix and Coil",
			sequence: []rune{'E', 'E', 'K', 'K', 'L', 'L'},
			expected: []Region{
				{start: 0, end: 6, score: 1.2933333333333332, structure: 'H'},
			},
		},
		{
			name:     "No Helix",
			sequence: []rune{'G', 'P', 'S', 'T', 'G', 'P'},
			expected: []Region{}, // No helix regions expected
		},
		{
			name:     "Extended Helix",
			sequence: []rune{'M', 'A', 'L', 'E', 'K', 'A', 'R', 'L'},
			expected: []Region{
				{start: 0, end: 8, score: 1.295, structure: 'H'}, // Example score
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := PredictHelix(tt.sequence)

			if !reflect.DeepEqual(normalizeRegions(result), normalizeRegions(tt.expected)) {
				t.Errorf("Test %s failed. Expected %v but got %v", tt.name, tt.expected, result)
			}
		})
	}
}

func TestIsHelix(t *testing.T) {
	tests := []struct {
		name     string // Name of the test case
		window   []rune // Input window of amino acids
		expected bool   // Expected output: true if it's a helix, false otherwise
	}{
		{
			name:     "Strong Helix",
			window:   []rune{'A', 'A', 'A', 'A', 'A', 'A'}, // All alanine
			expected: true,                                 // High helix propensity
		},
		{
			name:     "Weak Helix with Breaker",
			window:   []rune{'K', 'E', 'L', 'P', 'A', 'A'}, // Includes proline (P) as a breaker
			expected: false,                                // Proline breaks helix formation
		},
		{
			name:     "Very Low Propensity",
			window:   []rune{'G', 'G', 'G', 'S', 'S', 'T'}, // Multiple residues with low propensities
			expected: false,                                // Does not meet criteria for forming a helix
		},
		{
			name:     "Borderline Helix",
			window:   []rune{'M', 'Q', 'L', 'H', 'K', 'E'}, // Mixed strong and weak formers
			expected: true,                                 // Meets criteria for helix formation
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := IsHelix(tt.window)
			if result != tt.expected {
				t.Errorf("Test %s failed. Expected %v but got %v", tt.name, tt.expected, result)
			}
		})
	}
}

func TestPredictSheet(t *testing.T) {
	tests := []struct {
		name     string   // Name of the test case
		sequence []rune   // Input sequence
		expected []Region // Expected output regions
	}{
		{
			name:     "Strong Beta Sheet",
			sequence: []rune{'V', 'I', 'Y', 'F', 'L'}, // All strong sheet formers
			expected: []Region{
				{start: 0, end: 5, score: 1.4899999999999998, structure: 'E'}, // Example score
			},
		},
		{
			name:     "Mixed Strong and Weak",
			sequence: []rune{'V', 'I', 'G', 'Y', 'F'},
			expected: []Region{
				{start: 0, end: 5, score: 1.38, structure: 'E'}, // Example score
			},
		},
		{
			name:     "Extended Sheet",
			sequence: []rune{'V', 'V', 'I', 'I', 'Y', 'F', 'L'}, // Continuous strong formers
			expected: []Region{
				{start: 0, end: 7, score: 1.5357142857142858, structure: 'E'}, // Example score
			},
		},
		{
			name:     "No Sheet",
			sequence: []rune{'A', 'R', 'N', 'D', 'C'}, // Low sheet propensity amino acids
			expected: []Region{},                      // No sheet regions expected
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := PredictSheet(tt.sequence)

			if !reflect.DeepEqual(normalizeRegions(result), normalizeRegions(tt.expected)) {
				t.Errorf("Test %s failed. Expected %v but got %v", tt.name, tt.expected, result)
			}
		})
	}
}

func TestIsSheet(t *testing.T) {
	tests := []struct {
		name     string // Name of the test case
		window   []rune // Input window of amino acids
		expected bool   // Expected output: true if it's a sheet, false otherwise
	}{
		{
			name:     "Strong Beta Sheet",
			window:   []rune{'V', 'I', 'Y', 'F', 'L'}, // All strong sheet formers
			expected: true,                            // High sheet propensity
		},
		{
			name:     "Low Propensity",
			window:   []rune{'G', 'P', 'S', 'T', 'G'}, // Very low sheet propensity amino acids
			expected: false,                           // Does not meet criteria for forming a sheet
		},
		{
			name:     "Low Propensity",
			window:   []rune{'A', 'R', 'N', 'D', 'C'}, // Low sheet propensity amino acids
			expected: false,                           // Does not meet criteria for forming a sheet
		},
		{
			name:     "Borderline Sheet",
			window:   []rune{'T', 'Y', 'V', 'L', 'I'}, // Mixed formers and low breakers
			expected: true,                            // Meets criteria for sheet formation
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := IsSheet(tt.window)
			if result != tt.expected {
				t.Errorf("Test %s failed. Expected %v but got %v", tt.name, tt.expected, result)
			}
		})
	}
}

func TestPredictTurn(t *testing.T) {
	tests := []struct {
		name     string   // Name of the test case
		sequence []rune   // Input sequence
		expected []Region // Expected output regions
	}{
		{
			name:     "Strong Turn",
			sequence: []rune{'S', 'P', 'N', 'D'}, // High turn propensity
			expected: []Region{
				{start: 0, end: 4, score: 1.4925, structure: 'T'}, // Example score
			},
		},
		{
			name:     "Mixed Turn",
			sequence: []rune{'S', 'P', 'G', 'D'}, // Contains residues with moderate turn propensity
			expected: []Region{
				{start: 0, end: 4, score: 1.5125, structure: 'T'}, // Example score
			},
		},
		{
			name:     "Longer Sequence with Multiple Turns",
			sequence: []rune{'S', 'P', 'N', 'D', 'G'}, // Contains multiple turn regions
			expected: []Region{
				{start: 0, end: 4, score: 1.4925, structure: 'T'}, // Turn at residues 0–3
				{start: 1, end: 5, score: 1.545, structure: 'T'},  // Turn at residues 3–6
			},
		},
		{
			name:     "No Turn",
			sequence: []rune{'V', 'I', 'L', 'F'}, // Residues with very low turn propensity
			expected: []Region{},                 // No turn regions expected
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := PredictTurn(tt.sequence)

			if !reflect.DeepEqual(normalizeRegions(result), normalizeRegions(tt.expected)) {
				t.Errorf("Test %s failed. Expected %v but got %v", tt.name, tt.expected, result)
			}
		})
	}
}

func TestIsTurn(t *testing.T) {
	tests := []struct {
		name     string // Name of the test case
		window   []rune // Input window of amino acids
		expected bool   // Expected output: true if it's a turn, false otherwise
	}{
		{
			name:     "Strong Turn",
			window:   []rune{'S', 'P', 'N', 'D'}, // High turn propensity and positional probability
			expected: true,
		},
		{
			name:     "No Turn due to low positional probability",
			window:   []rune{'A', 'N', 'D', 'A'}, // Low positional probability
			expected: false,
		},
		{
			name:     "Mixed Turn",
			window:   []rune{'S', 'P', 'G', 'D'}, // Moderate propensity, high positional probability
			expected: true,
		},
		{
			name:     "No Turn due to Low propensity",
			window:   []rune{'V', 'I', 'L', 'F'}, // Low propensity
			expected: false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := IsTurn(tt.window)
			if result != tt.expected {
				t.Errorf("Test %s failed. Expected %v but got %v", tt.name, tt.expected, result)
			}
		})
	}
}

func TestCalculateAveragePropensity(t *testing.T) {
	tests := []struct {
		name          string  // Name of the test case
		window        []rune  // Input sequence (window)
		structureType rune    // Structure type ('H', 'E', 'T')
		expected      float64 // Expected average propensity
	}{
		{
			name:          "Single Residue",
			window:        []rune{'A'},
			structureType: 'H',
			expected:      1.42,
		},
		{
			name:          "Multiple Residues Beta-Sheet",
			window:        []rune{'B', 'E', 'F', 'G'},
			structureType: 'E',
			expected:      0.625,
		},
		{
			name:          "Long Sequence Alpha-Helix",
			window:        []rune{'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'},
			structureType: 'H',
			expected:      1.42,
		},
		{
			name:          "Empty Window",
			window:        []rune{},
			structureType: 'H', // No residues to calculate propensity for
			expected:      0.0, // Expect 0 as the average when window is empty
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := CalculateAveragePropensity(tt.window, tt.structureType)

			if result != tt.expected {
				t.Errorf("Test %s failed. Expected %v but got %v", tt.name, tt.expected, result)
			}
		})
	}
}

// Test function for ClassifyOverlap
func TestClassifyOverlap(t *testing.T) {
	// Helper function to compare slices
	compareSlices := func(a, b []rune) bool {
		if len(a) != len(b) {
			return false
		}
		for i := 0; i < len(a); i++ {
			if a[i] != b[i] {
				return false
			}
		}
		return true
	}

	// Define test cases
	tests := []struct {
		name         string
		sequence     []rune
		helixRegions []Region
		sheetRegions []Region
		turnRegions  []Region
		expected     []rune
	}{
		{
			name:     "No overlap",
			sequence: []rune{'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'},
			helixRegions: []Region{
				{start: 0, end: 2, structure: 'H', score: 10},
			},
			sheetRegions: []Region{
				{start: 3, end: 5, structure: 'E', score: 20},
			},
			turnRegions: []Region{
				{start: 6, end: 8, structure: 'T', score: 30},
			},
			expected: []rune{'H', 'H', 0, 'E', 'E', 0, 'T', 'T'},
		},
		{
			name:     "Two regions overlap",
			sequence: []rune{'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'},
			helixRegions: []Region{
				{start: 0, end: 4, structure: 'H', score: 10},
			},
			sheetRegions: []Region{
				{start: 2, end: 6, structure: 'E', score: 20},
			},
			turnRegions: []Region{
				{start: 4, end: 8, structure: 'T', score: 30},
			},
			expected: []rune{'H', 'H', 'E', 'E', 'T', 'T', 'T', 'T'},
		},
		{
			name:     "Three regions overlap",
			sequence: []rune{'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'},
			helixRegions: []Region{
				{start: 0, end: 3, structure: 'H', score: 10},
			},
			sheetRegions: []Region{
				{start: 2, end: 5, structure: 'E', score: 20},
			},
			turnRegions: []Region{
				{start: 2, end: 7, structure: 'T', score: 30},
			},
			expected: []rune{'H', 'H', 'T', 'T', 'T', 'T', 'T', 0},
		},
		{
			name:     "Full overlap with same score",
			sequence: []rune{'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'},
			helixRegions: []Region{
				{start: 0, end: 4, structure: 'H', score: 10},
			},
			sheetRegions: []Region{
				{start: 0, end: 4, structure: 'E', score: 10},
			},
			turnRegions: []Region{
				{start: 0, end: 4, structure: 'T', score: 10},
			},
			expected: []rune{'E', 'E', 'E', 'E', 0, 0, 0, 0},
		},
	}

	// Run tests
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := make([]rune, len(tt.sequence))
			// Call the ClassifyOverlap function
			ClassifyOverlap(tt.sequence, tt.helixRegions, tt.sheetRegions, tt.turnRegions, result)

			// Compare result with expected
			if !compareSlices(result, tt.expected) {
				t.Errorf("Expected: %v, got: %v", tt.expected, result)
			}
		})
	}
}

func TestExtendHelix(t *testing.T) {
	tests := []struct {
		name          string
		sequence      []rune
		start         int
		expectedStart int
		expectedEnd   int
		expectedScore float64
	}{
		{
			name:          "Helix extends only forward",
			sequence:      []rune{'V', 'V', 'V', 'I', 'I', 'I', 'E', 'E', 'E', 'L', 'L', 'L', 'P', 'G', 'G'},
			start:         0,
			expectedStart: 0,
			expectedEnd:   10,
			expectedScore: 1.216,
		},
		{
			name:          "Helix extends only backwards",
			sequence:      []rune{'G', 'G', 'G', 'G', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'},
			start:         8,
			expectedStart: 5,
			expectedEnd:   14,
			expectedScore: 1.42,
		},
		{
			name:          "Helix extends both ways",
			sequence:      []rune{'G', 'G', 'G', 'G', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'},
			start:         8,
			expectedStart: 5,
			expectedEnd:   16,
			expectedScore: 1.42,
		},
		{
			name:          "Helix extends both ways",
			sequence:      []rune{'G', 'G', 'G', 'G', 'A', 'A', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G'},
			start:         4,
			expectedStart: 4,
			expectedEnd:   10,
			expectedScore: 1.42,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			start, end, score := ExtendHelix(tt.sequence, tt.start)

			// Validate the start, end, and score
			if start != tt.expectedStart {
				t.Errorf("expected start %d, got %d", tt.expectedStart, start)
			}
			if end != tt.expectedEnd {
				t.Errorf("expected end %d, got %d", tt.expectedEnd, end)
			}
			if score != tt.expectedScore {
				t.Errorf("expected score %f, got %f", tt.expectedScore, score)
			}
		})
	}
}
func TestExtendSheet(t *testing.T) {
	tests := []struct {
		name          string
		sequence      []rune
		start         int
		expectedStart int
		expectedEnd   int
		expectedScore float64
	}{
		{
			name:          "Sheet extends only forward",
			sequence:      []rune{'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'D', 'D', 'D', 'D'},
			start:         0,
			expectedStart: 0,
			expectedEnd:   6,
			expectedScore: 1.47,
		},
		{
			name:          "Sheet extends only backward",
			sequence:      []rune{'D', 'D', 'D', 'D', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y'},
			start:         6,
			expectedStart: 4,
			expectedEnd:   11,
			expectedScore: 1.47,
		},
		{
			name:          "No extension possible (max extension)",
			sequence:      []rune{'D', 'D', 'D', 'D', 'Y', 'Y', 'Y', 'Y', 'Y', 'D', 'D', 'D', 'D'},
			start:         4,
			expectedStart: 4,
			expectedEnd:   9,
			expectedScore: 1.47,
		},
		{
			name:          "Extends both ways",
			sequence:      []rune{'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y'},
			start:         4,
			expectedStart: 0,
			expectedEnd:   15,
			expectedScore: 1.47,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Run ExtendSheet and capture results
			start, end, score := ExtendSheet(tt.sequence, tt.start)

			// Test start position
			if start != tt.expectedStart {
				t.Errorf("Expected start %d, but got %d", tt.expectedStart, start)
			}

			// Test end position
			if end != tt.expectedEnd {
				t.Errorf("Expected end %d, but got %d", tt.expectedEnd, end)
			}

			// Test score with some tolerance
			const epsilon = 1e-6
			if math.Abs(tt.expectedScore-score) > epsilon {
				t.Errorf("Expected score %v, but got %v", tt.expectedScore, score)
			}
		})
	}
}

func TestOverlap(t *testing.T) {
	tests := []struct {
		name     string
		region1  Region
		region2  Region
		expected bool
	}{
		{
			name:     "Regions overlap completely",
			region1:  Region{start: 1, end: 5},
			region2:  Region{start: 3, end: 6},
			expected: true, // Overlap from 3 to 5
		},
		{
			name:     "Regions are completely disjoint",
			region1:  Region{start: 1, end: 5},
			region2:  Region{start: 6, end: 10},
			expected: false, // No overlap
		},
		{
			name:     "Region1 is completely within Region2",
			region1:  Region{start: 3, end: 5},
			region2:  Region{start: 1, end: 6},
			expected: true, // Overlap from 3 to 5
		},
		{
			name:     "Regions touch but do not overlap",
			region1:  Region{start: 1, end: 5},
			region2:  Region{start: 5, end: 10},
			expected: false, // They touch at 5 but do not overlap
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Overlap(tt.region1, tt.region2)
			if result != tt.expected {
				t.Errorf("expected %v, got %v", tt.expected, result)
			}
		})
	}
}

func TestMax(t *testing.T) {
	tests := []struct {
		name     string
		a        int
		b        int
		expected int
	}{
		{
			name:     "a is greater than b",
			a:        5,
			b:        3,
			expected: 5, // 5 is greater than 3
		},
		{
			name:     "b is greater than a",
			a:        2,
			b:        8,
			expected: 8, // 8 is greater than 2
		},
		{
			name:     "a is equal to b",
			a:        7,
			b:        7,
			expected: 7, // Both are equal, so return either
		},
		{
			name:     "a is negative and b is positive",
			a:        -3,
			b:        4,
			expected: 4, // 4 is greater than -3
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Max(tt.a, tt.b)
			if result != tt.expected {
				t.Errorf("expected %d, got %d", tt.expected, result)
			}
		})
	}
}

func TestMin(t *testing.T) {
	tests := []struct {
		name     string
		a        int
		b        int
		expected int
	}{
		{
			name:     "a is smaller than b",
			a:        3,
			b:        5,
			expected: 3, // 3 is smaller than 5
		},
		{
			name:     "b is smaller than a",
			a:        8,
			b:        2,
			expected: 2, // 2 is smaller than 8
		},
		{
			name:     "a is equal to b",
			a:        7,
			b:        7,
			expected: 7, // Both are equal, so return either
		},
		{
			name:     "a is negative and b is positive",
			a:        -4,
			b:        3,
			expected: -4, // -4 is smaller than 3
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Min(tt.a, tt.b)
			if result != tt.expected {
				t.Errorf("expected %d, got %d", tt.expected, result)
			}
		})
	}
}
