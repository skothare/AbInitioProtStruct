// Group 2: Ab Initio Secondary Structure Prediction of Proteins
// Date: 12th December, 2024
// Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare

package main

// AminoAcidPropensities represents the propensity of an amino acid to form different secondary structures
// Each float64 value indicates the likelihood of forming:
// - alpha Helix
// - beta Sheet
// - Turn
type AminoAcidPropensities struct {
	alphaHelix float64 // Propensity to form an alpha helix
	betaSheet  float64 // Propensity to form a beta sheet
	turn       float64 // Propensity to form a turn
}

// BendProbabilities represents the positional probabilities of an amino acid causing a bend in a protein's secondary structure
// p1, p2, p3, p4 correspond to different positions in a tetrapeptide window
type BendProbabilities struct {
	p1, p2, p3, p4 float64
}

// Region represents a section of a protein sequence with a specific secondary structure
type Region struct {
	start     int     // Starting index of the region
	end       int     // Ending index of the region
	score     float64 // Propensity score for the region's structure
	structure rune    // Type of secondary structure (H: Helix, E: Sheet, T: Turn, C: Coil)
}

// propensities is a lookup table of secondary structure propensities for each amino acid
// The values have been taken from the Chou-Fasman Paper
// Keys are amino acid one-letter codes
// Values indicate propensities for alpha helix, beta sheet, and turn formations
var propensities = map[rune]AminoAcidPropensities{
	'A': {1.42, 0.83, 0.66}, // Alanine
	'R': {0.98, 0.93, 0.95}, // Arginine
	'N': {0.67, 0.89, 1.56}, // Asparagine
	'D': {1.01, 0.54, 1.46}, // Aspartic Acid
	'C': {0.70, 1.19, 1.19}, // Cysteine
	'Q': {1.11, 1.10, 0.98}, // Glutamine
	'E': {1.51, 0.37, 0.74}, // Glutamic Acid
	'G': {0.57, 0.75, 1.64}, // Glycine
	'H': {1.00, 0.87, 0.95}, // Histidine
	'I': {1.08, 1.60, 0.47}, // Isoleucine
	'L': {1.21, 1.30, 0.59}, // Leucine
	'K': {1.16, 0.74, 1.01}, // Lysine
	'M': {1.45, 1.05, 0.60}, // Methionine
	'F': {1.13, 1.38, 0.60}, // Phenylalanine
	'P': {0.57, 0.55, 1.52}, // Proline
	'S': {0.77, 0.75, 1.43}, // Serine
	'T': {0.83, 1.19, 0.96}, // Threonine
	'W': {1.08, 1.37, 0.96}, // Tryptophan
	'Y': {0.69, 1.47, 1.14}, // Tyrosine
	'V': {1.06, 1.70, 0.50}, // Valine
}

// bendProbabilitiesTable is a lookup table of bend probabilities for each amino acid
// Provides positional probabilities (p1, p2, p3, p4) for causing a bend in a tetrapeptide
// These probabilities are used in turn prediction algorithms
var bendProbabilitiesTable = map[rune]BendProbabilities{
	'A': {0.060, 0.076, 0.035, 0.058},
	'R': {0.070, 0.106, 0.099, 0.085},
	'N': {0.161, 0.083, 0.191, 0.091},
	'D': {0.147, 0.110, 0.179, 0.081},
	'C': {0.149, 0.050, 0.117, 0.128},
	'Q': {0.074, 0.098, 0.037, 0.098},
	'E': {0.056, 0.060, 0.077, 0.064},
	'G': {0.102, 0.085, 0.190, 0.152},
	'H': {0.140, 0.047, 0.093, 0.054},
	'I': {0.043, 0.034, 0.013, 0.056},
	'L': {0.061, 0.025, 0.036, 0.070},
	'K': {0.055, 0.115, 0.072, 0.095},
	'M': {0.068, 0.082, 0.014, 0.055},
	'F': {0.059, 0.041, 0.065, 0.065},
	'P': {0.102, 0.301, 0.034, 0.068},
	'S': {0.120, 0.139, 0.125, 0.106},
	'T': {0.086, 0.108, 0.065, 0.079},
	'W': {0.089, 0.073, 0.064, 0.167},
	'Y': {0.082, 0.065, 0.114, 0.125},
	'V': {0.062, 0.048, 0.028, 0.053},
}

// Following type is a map that contains the information values for a given structure.
type InfoValTable map[string][]float64

// GORPredictionResult holds the scores and predicted structure for a residue
type GORPredictionResult struct {
	Position int    // Position of the residue in the sequence
	Residue  string // The amino acid at that position
	// Total scores (tallies) per structure
	ScoreAlpha         float64
	ScoreBeta          float64
	ScoreTurn          float64
	ScoreCoil          float64
	PredictedStructure string // Predicted structure string
}

// GORPropensities represents context-based propensities for secondary structures
type GORPropensities struct {
	helixMatrix map[rune]float64
	sheetMatrix map[rune]float64
	turnMatrix  map[rune]float64
}

// GORModelPredictor stores the probability matrices used in the GOR method
var GORModel = GORPropensities{
	helixMatrix: map[rune]float64{
		// Simplified example values (actual values would be based on a statistical model)
		'A': 1.2, 'C': 0.5, 'D': 0.8, 'E': 1.3, 'F': 0.9,
		// Add more amino acids here...
	},
	sheetMatrix: map[rune]float64{
		'A': 0.8, 'C': 1.3, 'D': 0.7, 'E': 0.6, 'F': 1.4,
		// Add more amino acids here...
	},
	turnMatrix: map[rune]float64{
		'A': 0.9, 'C': 0.7, 'D': 1.4, 'E': 0.8, 'F': 0.5,
		// Add more amino acids here...
	},
}

// HMM represents a Hidden Markov Model.
type HMM struct {
	States        []string       // List of possible states (e.g., Helix, Sheet, Turn, Coil)
	Symbols       []string       // List of possible observation symbols (e.g., amino acids)
	Transition    [][]float64    // Transition probabilities between states
	Emission      [][]float64    // Emission probabilities for symbols given a state
	Initial       []float64      // Initial probabilities of each state
	StateMapping  map[string]int // Maps state names to indices for easier lookup
	SymbolMapping map[string]int // Maps symbol names to indices for easier lookup
}
