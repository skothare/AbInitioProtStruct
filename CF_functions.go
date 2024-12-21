// Group 2: Ab Initio Secondary Structure Prediction of Proteins
// Date: 12th December, 2024
// Members: Arth Banka, Riti Bhatia, Sanchitha Kuthethoor, Sumeet Kothare

package main

// ChouFasmanPredictSS()
// Input: a slice of runes, each elements corresponds to an amino acid residue
// Output: a string that predicts the secondary structure of a protein by employing the Chou-Fasman model
func ChouFasmanPredictSS(sequence []rune) string {
	length := len(sequence)
	result := make([]rune, length)

	// Default for all residues is coil
	for i := range result {
		result[i] = 'C'
	}

	// First, predict for each structure type respectively
	helixRegions := PredictHelix(sequence)

	sheetRegions := PredictSheet(sequence)

	turnRegions := PredictTurn(sequence)

	// Resolve the overlapping regions and assign structure types to the result slice accordingly
	ClassifyOverlap(sequence, helixRegions, sheetRegions, turnRegions, result)

	return string(result)
}

// PredictHelix()
// Input: a slice of runes sequence
// Output: a slice of Region datatypes of all the regions with high propensities of being an alpha helix
func PredictHelix(sequence []rune) []Region {
	var regions []Region

	// Slide a 6-residue window across the sequence
	for i := 0; i <= len(sequence)-6; i++ {
		if IsHelix(sequence[i : i+6]) { // Check 6-residue window
			// If it's a potential helix nucleation site, extend it from either sides
			start, end, score := ExtendHelix(sequence, i)
			regions = append(regions, Region{start, end, score, 'H'}) //Append the helix region
			i = end - 1                                               // Skip ahead to avoid overlapping predictions
		}
	}

	return regions
}

// IsHelix()
// Input: a window of length 6 (slice of runes)
// Output: a boolean after checking if the window (nucleation site) has a high propensity of forming an alpha helix
func IsHelix(window []rune) bool {
	if len(window) < 6 {
		return false
	}

	formerCount := 0.0
	breakerCount := 0
	avgPropensity := 0.0

	// Analyze each amino acid in the window
	for _, aa := range window {
		// Immediately return false if Proline is found
		if aa == 'P' {
			return false
		}
		prop := propensities[aa]
		avgPropensity += prop.alphaHelix

		// Count helix formers and breakers
		if prop.alphaHelix >= 1.05 {
			formerCount += 1.0 // Strong former
		} else if prop.alphaHelix >= 1.0 {
			formerCount += 0.5 // Weak former
		} else if prop.alphaHelix <= 0.69 {
			breakerCount++
		} else if aa == 'P' {
			return false
		}
	}

	// Calculate the average propensity
	avgPropensity /= float64(len(window))

	// Check if the window meets helix formation criteria
	return formerCount >= 4.0 &&
		breakerCount < 2 &&
		avgPropensity >= 1.03
}

// PredictSheet()
// Input: a slice of runes sequence
// Output: a slice of Region datatypes of all the regions with high propensities of being a beta sheet
func PredictSheet(sequence []rune) []Region {
	var regions []Region

	// Slide a 5-residue window across the sequence
	for i := 0; i <= len(sequence)-5; i++ { // 5 residue window
		if IsSheet(sequence[i : i+5]) {
			// If it's a potential sheet nucleation site, extend the region
			start, end, score := ExtendSheet(sequence, i)
			regions = append(regions, Region{start, end, score, 'E'}) // Append sheet region
			i = end - 1
		}
	}

	return regions
}

// IsSheet()
// Input: a window of length 5 (slice of runes)
// Output: a boolean after checking if the window (nucleation site) has a high propensity of forming a beta sheet
func IsSheet(window []rune) bool {
	if len(window) < 5 {
		return false
	}

	formerCount := 0
	breakerCount := 0
	avgPropensity := 0.0

	// Analyze each amino acid in the window
	for _, aa := range window {
		prop := propensities[aa]
		avgPropensity += prop.betaSheet

		// Count sheet formers and breakers
		if prop.betaSheet >= 1.0 {
			formerCount++ // Former
		} else if prop.betaSheet <= 0.75 {
			breakerCount++ // Breaker
		}
	}

	// Calculate average propensity
	avgPropensity /= float64(len(window))

	// Check if the window meets sheet formation criteria
	return formerCount >= 3 &&
		breakerCount <= 1 &&
		avgPropensity >= 1.05
}

// ExtendHelix()
// Input: the sequence which is a slice of runes, and a int start corresponding to an index value in the slice
// Output: the start and end ints of the helix, and a float corresponding to the score (avg propensity)
func ExtendHelix(sequence []rune, start int) (int, int, float64) {
	// Start with initial 6-residue nucleation site
	currentStart := start
	currentEnd := start + 6

	// Continue extending as long as possible in either direction
	for {
		canStillExtend := false

		// Try extending forward if not at sequence end
		if currentEnd < len(sequence) {
			forwardRegion := sequence[currentStart : currentEnd+1]
			forwardProp := CalculateAveragePropensity(forwardRegion, 'H')

			// Check for tetrapeptide breakers in the forward direction
			hasBreaker := false
			if currentEnd+4 <= len(sequence) {
				for i := currentEnd - 3; i <= currentEnd; i++ {
					tetrapeptide := sequence[i : i+4]
					tetraProp := CalculateAveragePropensity(tetrapeptide, 'H')
					if tetraProp < 1.00 {
						hasBreaker = true
						break // If set of tetrapeptide breakers identified, then break the loop
					}
				}
			}

			// Extend if meets propensity and no breakers
			if forwardProp >= 1.03 && !hasBreaker {
				currentEnd++
				canStillExtend = true
			}
		}

		// Try extending backward if not at sequence start
		if currentStart > 0 {
			backwardRegion := sequence[currentStart-1 : currentEnd]
			backwardProp := CalculateAveragePropensity(backwardRegion, 'H')

			// Check for tetrapeptide breakers in the backward direction
			hasBreaker := false
			if currentStart >= 3 {
				for i := currentStart - 3; i <= currentStart; i++ {
					if i+4 <= currentEnd {
						tetrapeptide := sequence[i : i+4]
						tetraProp := CalculateAveragePropensity(tetrapeptide, 'H')
						if tetraProp < 1.00 {
							hasBreaker = true
							break
						}
					}
				}
			}
			// Extend if meets propensity and no breakers
			if backwardProp >= 1.03 && !hasBreaker {
				currentStart--
				canStillExtend = true
			}
		}

		// Stop if no more extension possible
		if !canStillExtend {
			break
		}
	}

	// Calculate final region score
	finalRegion := sequence[currentStart:currentEnd]
	totalScore := CalculateAveragePropensity(finalRegion, 'H')

	return currentStart, currentEnd, totalScore
}

// ExtendSheet()
// Input: the sequence which is a slice of runes, and a int start corresponding to an index value in the slice
// Output: the start and end ints of the sheet, and a float corresponding to the score (avg propensity)
func ExtendSheet(sequence []rune, start int) (int, int, float64) {
	// Start with initial 5-residue nucleus
	currentStart := start
	currentEnd := start + 5

	// Continue extending as long as possible in either direction
	for {
		canStillExtend := false

		// Try extending forward if not at sequence end
		if currentEnd < len(sequence) {
			forwardRegion := sequence[currentStart : currentEnd+1]
			forwardProp := CalculateAveragePropensity(forwardRegion, 'E')

			// Check for tetrapeptide breakers in the forward direction
			hasBreaker := false
			if currentEnd+4 <= len(sequence) {
				for i := currentEnd - 3; i <= currentEnd; i++ {
					tetrapeptide := sequence[i : i+4]
					tetraProp := CalculateAveragePropensity(tetrapeptide, 'E')
					if tetraProp < 1.00 {
						hasBreaker = true
						break
					}
				}
			}
			// Extend if meets propensity and no breakers
			if forwardProp >= 1.05 && !hasBreaker {
				currentEnd++
				canStillExtend = true
			}
		}

		// Try extending backward if not at sequence start
		if currentStart > 0 {
			backwardRegion := sequence[currentStart-1 : currentEnd]
			backwardProp := CalculateAveragePropensity(backwardRegion, 'E')

			// Check for tetrapeptide breakers in the backward direction
			hasBreaker := false
			if currentStart >= 3 {
				for i := currentStart - 3; i <= currentStart; i++ {
					if i+4 <= currentEnd {
						tetrapeptide := sequence[i : i+4]
						tetraProp := CalculateAveragePropensity(tetrapeptide, 'E')
						if tetraProp < 1.00 {
							hasBreaker = true
							break
						}
					}
				}
			}
			// Extend if meets propensity and no breakers
			if backwardProp >= 1.05 && !hasBreaker {
				currentStart--
				canStillExtend = true
			}
		}

		// Stop if no more extension possible
		if !canStillExtend {
			break
		}
	}
	// Calculate final region score
	finalRegion := sequence[currentStart:currentEnd]
	totalScore := CalculateAveragePropensity(finalRegion, 'E')

	return currentStart, currentEnd, totalScore
}

// PredictTurn()
// Input: sequence which is a slice of runes
// Output: a slice of Region datatypes of all the regions likely of being turns
func PredictTurn(sequence []rune) []Region {
	var regions []Region

	// Slide a 4-residue window across the sequence
	for i := 0; i <= len(sequence)-4; i++ {
		if IsTurn(sequence[i : i+4]) { // 4 residue window
			regions = append(regions, Region{i, i + 4, CalculateAveragePropensity(sequence[i:i+4], 'T'), 'T'}) //Append the calculated score
		}
	}

	return regions
}

// IsTurn()
// Input: a window of length 4 (slice of runes)
// Output: a boolean after checking if the window (nucleation site) has a high propensity of being a turn
func IsTurn(window []rune) bool {
	if len(window) < 4 {
		return false
	}

	// Calculate positional probability of tetrapeptide being a turn
	pt := bendProbabilitiesTable[window[0]].p1 *
		bendProbabilitiesTable[window[1]].p2 *
		bendProbabilitiesTable[window[2]].p3 *
		bendProbabilitiesTable[window[3]].p4

	// Ensure the two middle residues meet a minimum bend probability
	if propensities[window[1]].turn < 0.5 || propensities[window[2]].turn < 0.5 {
		return false
	}

	// Calculate the average turn propensity
	avgPropensity := CalculateAveragePropensity(window, 'T')

	// Check if the window meets turn formation criteria
	return avgPropensity >= 1.0 && pt >= 0.000075
}

// ClassifyOverlap()
// Input: sequence as a slice of runes, the helix Regions, sheet Regions, and turn Regions which are all slices of Region datatypes
// and result which is a slice of rines as well
// Output: it does not output anything, but assigns the overlapping helix/sheet/turn regions based on which has the highest propensity scores
func ClassifyOverlap(sequence []rune, helixRegions, sheetRegions, turnRegions []Region, result []rune) {
	// Combine all regions into a single slice
	allRegions := make([]Region, 0)
	allRegions = append(allRegions, helixRegions...)
	allRegions = append(allRegions, sheetRegions...)
	allRegions = append(allRegions, turnRegions...)

	// First mark all regions with their respective structures
	for _, region := range allRegions {
		for i := region.start; i < region.end; i++ {
			result[i] = region.structure
		}
	}

	// Handle overlaps by comparing regions and selecting based on score
	for i, region1 := range allRegions {
		// Only compare with regions that come after region1 to avoid redundant comparisons
		for _, region2 := range allRegions[i+1:] {
			if Overlap(region1, region2) {
				overlapStart := Max(region1.start, region2.start)
				overlapEnd := Min(region1.end, region2.end)

				// Choose structure with higher score
				structure := region1.structure
				if region2.score > region1.score {
					structure = region2.structure
				}

				// Update the overlapping region
				for i := overlapStart; i < overlapEnd; i++ {
					result[i] = structure
				}
			}
		}
	}
}

// CalculateAveragePropensity()
// Input: window which is a slice of runes, and a structureType rune
// Output: the avergage propensities of all the residues in the window for being that structure type as a float64
func CalculateAveragePropensity(window []rune, structureType rune) float64 {
	if len(window) == 0 {
		return 0.0
	}

	sum := 0.0
	for _, aa := range window {
		prop := propensities[aa]
		switch structureType {
		case 'H':
			sum += prop.alphaHelix
		case 'E':
			sum += prop.betaSheet
		case 'T':
			sum += prop.turn
		}
	}
	return sum / float64(len(window))
}

// Overlap()
// Input: two Regions
// Output: a boolean corresponding to whether the two regions overlap with each other
func Overlap(r1, r2 Region) bool {
	return r1.start < r2.end && r2.start < r1.end
}

// Max()
// Input: two integers a and b
// Ouput: the largest integer between the two
func Max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// Min()
// Input: two integers a and b
// Ouput: the smallest integer between the two
func Min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
