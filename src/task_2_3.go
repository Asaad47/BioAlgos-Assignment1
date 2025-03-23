package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"sort"
	"strings"
)

// MinimizerIndex stores selected k-mers as minimizers
// mapped to their occurrences in reference genomes
type MinimizerIndex struct {
	minimizers map[string]map[string]int // minimizer -> genome -> count
	k          int
	w          int
}

// SequenceReadMatch tracks matches for a single sequence read
type SequenceReadMatch struct {
	readID       string
	matchedOrgs  map[string]int // maps organism to number of minimizer matches
	totalMatches int            // total number of minimizer matches across all organisms
}

// getMinimizer selects the representative k-mer from a window,
// which is the lexicographically smallest k-mer in this case
func getMinimizer(kmers []string) string {
	sort.Strings(kmers) // Lexicographic order
	return kmers[0]     // choose the smallest k-mer as minimizer
}

// buildMinimizerIndex creates an index storing only minimizers
func buildMinimizerIndex(genomeFiles []string, k, w int) *MinimizerIndex {
	minimizerIndex := &MinimizerIndex{minimizers: make(map[string]map[string]int), k: k, w: w}

	for _, genomeFile := range genomeFiles {
		file, err := os.Open(genomeFile)
		if err != nil {
			log.Fatalf("Failed to open file: %v", err)
		}
		defer file.Close()

		genomeName := getOrganismShortName(genomeFile)
		fmt.Printf("Processing genome: %s\n", genomeName)
		genomeLength := 0

		scanner := bufio.NewScanner(file)
		// buffer to store the last w+k-2 characters from previous line
		var prevChars []byte
		prevCharsLen := 0

		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, ">") {
				continue
			}
			line = strings.ToLower(strings.TrimSpace(line))
			genomeLength += len(line)
			lineBytes := []byte(line)

			// process minimizers that span the boundary between prevChars and current line
			if prevCharsLen > 0 {
				// create a window that spans the boundary
				boundaryWindow := make([]byte, w+k-1)
				copy(boundaryWindow[:], prevChars[:prevCharsLen])
				copy(boundaryWindow[prevCharsLen:], lineBytes[:w+k-1-prevCharsLen])

				// extract k-mers from the boundary window
				kmers := make([]string, 0, w)
				for j := 0; j < w; j++ {
					kmers = append(kmers, string(boundaryWindow[j:j+k]))
				}
				minimizer := getMinimizer(kmers)
				if _, exists := minimizerIndex.minimizers[minimizer]; !exists {
					minimizerIndex.minimizers[minimizer] = make(map[string]int)
				}
				minimizerIndex.minimizers[minimizer][genomeName]++
			}

			// process minimizers entirely within the current line
			for i := 0; i <= len(lineBytes)-w-k+1; i++ {
				window := lineBytes[i : i+w+k-1]
				kmers := make([]string, 0, w)
				for j := 0; j < w; j++ {
					kmers = append(kmers, string(window[j:j+k]))
				}
				minimizer := getMinimizer(kmers)
				if _, exists := minimizerIndex.minimizers[minimizer]; !exists {
					minimizerIndex.minimizers[minimizer] = make(map[string]int)
				}
				minimizerIndex.minimizers[minimizer][genomeName]++
			}

			// store the last w+k-2 characters for the next line
			if len(line) >= w+k-2 {
				prevChars = make([]byte, w+k-2)
				copy(prevChars, lineBytes[len(lineBytes)-(w+k-2):])
				prevCharsLen = w + k - 2
			} else {
				// if line is shorter than w+k-2, shift existing chars and prepend new ones
				newPrevChars := make([]byte, w+k-2)
				if prevCharsLen > 0 {
					shift := len(line)
					copy(newPrevChars[:w+k-2-shift], prevChars[shift:])
				}
				copy(newPrevChars[w+k-2-len(line):], lineBytes)
				prevChars = newPrevChars
				prevCharsLen = w + k - 2
			}
		}
		fmt.Printf("Genome length (%s): %d\n", genomeName, genomeLength)
	}
	return minimizerIndex
}

// classifyReadsMinimizer classifies reads using the minimizer index
func classifyReadsMinimizer(readFiles []string, index *MinimizerIndex) map[string]*SequenceReadMatch {
	readMatches := make(map[string]*SequenceReadMatch)

	for _, readFile := range readFiles {
		file, err := os.Open(readFile)
		if err != nil {
			log.Fatalf("Failed to open file: %v", err)
		}
		defer file.Close()

		scanner := bufio.NewScanner(file)
		var readID, sequence string
		lineNum := 0

		for scanner.Scan() {
			line := scanner.Text()
			switch lineNum % 4 {
			case 0: // Header line
				if strings.HasPrefix(line, "@") {
					readID = strings.TrimSpace(line[1:]) + "_" + readFile
					readMatches[readID] = &SequenceReadMatch{
						readID:      readID,
						matchedOrgs: make(map[string]int),
					}
				}
			case 1: // Sequence line
				sequence = strings.ToLower(line)
				for i := 0; i <= len(sequence)-index.w-index.k+1; i++ {
					window := sequence[i : i+index.w+index.k-1]
					kmers := make([]string, 0, index.w)
					for j := 0; j < index.w; j++ {
						kmers = append(kmers, window[j:j+index.k])
					}
					minimizer := getMinimizer(kmers)

					if genomeMatches, exists := index.minimizers[minimizer]; exists {
						for genome, count := range genomeMatches {
							readMatches[readID].matchedOrgs[genome] += count
							readMatches[readID].totalMatches += count
						}
					}
				}
			}
			lineNum++
		}
	}
	return readMatches
}

func getOrganismShortName(path string) string {
	filename := path
	if strings.Contains(filename, "GCF_000005845") {
		return "E. coli"
	} else if strings.Contains(filename, "GCF_000009045") {
		return "B. subtilis"
	} else if strings.Contains(filename, "GCF_000006765") {
		return "P. aeruginosa"
	} else if strings.Contains(filename, "GCF_000013425") {
		return "S. aureus"
	} else if strings.Contains(filename, "GCF_000195955") {
		return "M. tuberculosis"
	}
	return filename
}

func main() {
	k := 31
	w := 10
	genomeFiles := []string{
		"../data/1_ecol_ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna",
		"../data/2_bsub_ncbi_dataset/ncbi_dataset/data/GCF_000009045.1/GCF_000009045.1_ASM904v1_genomic.fna",
		"../data/3_paer_ncbi_dataset/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna",
		"../data/4_saur_ncbi_dataset/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna",
		"../data/5_mtub_ncbi_dataset/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna",
	}

	readFiles := []string{
		"../data/sequence_reads/simulated_reads_no_errors_10k_R1.fastq",
		"../data/sequence_reads/simulated_reads_no_errors_10k_R2.fastq",
		"../data/sequence_reads/simulated_reads_miseq_10k_R1.fastq",
		"../data/sequence_reads/simulated_reads_miseq_10k_R2.fastq",
	}

	fmt.Println("\n" + strings.Repeat("=", 80))
	fmt.Println("Minimizer-Based Classification Report")
	fmt.Println(strings.Repeat("=", 80))

	fmt.Printf("\nBuilding minimizer index with k=%d and w=%d\n", k, w)
	index := buildMinimizerIndex(genomeFiles, k, w)

	fmt.Printf("\nClassifying reads using minimizers.\n")
	readMatches := classifyReadsMinimizer(readFiles, index)

	// Calculate statistics
	orgReadCounts := make(map[string]int)
	orgMinimizerCounts := make(map[string]int)
	multipleMatches := 0
	uniqueMatches := 0
	noMatches := 0

	for _, match := range readMatches {
		if len(match.matchedOrgs) > 1 {
			multipleMatches++
		} else if len(match.matchedOrgs) == 1 {
			uniqueMatches++
		} else {
			noMatches++
		}
		for org := range match.matchedOrgs {
			orgReadCounts[org]++
			orgMinimizerCounts[org] += match.matchedOrgs[org]
		}
	}

	// Report results
	fmt.Printf("\n1. Classification Results:\n")
	fmt.Printf("    Total sequence reads processed: %d\n", len(readMatches))
	fmt.Printf("    Matched (minimizer) reads per organism:\n")
	for _, orgName := range []string{"E. coli", "B. subtilis", "P. aeruginosa", "S. aureus", "M. tuberculosis"} {
		fmt.Printf("     %-15s: %d reads, %d minimizer matches\n", orgName, orgReadCounts[orgName], orgMinimizerCounts[orgName])
	}

	fmt.Printf("\n2. Match Statistics:\n")
	fmt.Printf("    Reads with unique matches: %d (%.2f%%)\n",
		uniqueMatches, float64(uniqueMatches)*100/float64(len(readMatches)))
	fmt.Printf("    Reads with multiple matches: %d (%.2f%%)\n",
		multipleMatches, float64(multipleMatches)*100/float64(len(readMatches)))
	fmt.Printf("    Reads with no matches: %d (%.2f%%)\n",
		noMatches, float64(noMatches)*100/float64(len(readMatches)))

	fmt.Println("\n" + strings.Repeat("=", 80))
}
