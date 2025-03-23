package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strings"
)

// SequenceReadMatch tracks matches for a single sequence read
type SequenceReadMatch struct {
	readID       string
	matchedOrgs  map[string]int // maps organism to number of k-mer matches
	totalMatches int            // total number of k-mer matches across all organisms
}

// KmerStats tracks the occurrence of a k-mer across genomes
type KmerStats struct {
	occurrences map[string]int // maps genome name to count
	totalCount  int            // total occurrences across all genomes
}

func buildKmerIndex(genomeFiles []string, k int) (map[string]*KmerStats, int) {
	kmerIndex := make(map[string]*KmerStats)
	totalGenomeLength := 0

	for _, genomeFile := range genomeFiles {
		fastaFile, err := os.Open(genomeFile)
		if err != nil {
			log.Fatalf("Failed to open FASTA file: %v", err)
		}
		defer fastaFile.Close()

		orgName := getOrganismShortName(genomeFile)
		fmt.Printf("Processing genome: %s\n", orgName)
		genomeLength := 0

		scanner := bufio.NewScanner(fastaFile)
		// buffer to store the last k-1 characters from previous line
		var prevChars [31]byte // fixed size array for k-1 characters (k=31)
		prevCharsLen := 0

		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, ">") {
				continue
			}
			line = strings.ToLower(strings.TrimSpace(line))
			genomeLength += len(line)
			lineBytes := []byte(line)

			// process k-mers that span the boundary between prevChars and current line
			if prevCharsLen > 0 {
				var boundaryKmer [31]byte
				copy(boundaryKmer[:], prevChars[:prevCharsLen])
				copy(boundaryKmer[prevCharsLen:], lineBytes[:k-prevCharsLen])

				kmer := string(boundaryKmer[:])
				if _, exists := kmerIndex[kmer]; !exists {
					kmerIndex[kmer] = &KmerStats{
						occurrences: make(map[string]int),
					}
				}
				kmerIndex[kmer].occurrences[orgName]++
				kmerIndex[kmer].totalCount++
			}

			// process k-mers entirely within the current line
			for i := 0; i <= len(lineBytes)-k; i++ {
				kmer := string(lineBytes[i : i+k])

				if _, exists := kmerIndex[kmer]; !exists {
					kmerIndex[kmer] = &KmerStats{
						occurrences: make(map[string]int),
					}
				}

				kmerIndex[kmer].occurrences[orgName]++
				kmerIndex[kmer].totalCount++
			}

			// store the last k-1 characters for the next line
			if len(line) >= k-1 {
				copy(prevChars[:], lineBytes[len(lineBytes)-(k-1):])
				prevCharsLen = k - 1
			} else {
				// if line is shorter than k-1, shift existing chars and prepend new ones
				shift := len(line)
				copy(prevChars[:k-1-shift], prevChars[shift:])
				copy(prevChars[k-1-shift:], lineBytes)
				prevCharsLen = k - 1
			}
		}
		fmt.Printf("Genome length (%s): %d\n", orgName, genomeLength)
		totalGenomeLength += genomeLength
	}

	return kmerIndex, totalGenomeLength
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

func extractKmers(sequence string, k int) []string {
	sequence = strings.ToLower(strings.TrimSpace(sequence))
	kmers := make([]string, 0, len(sequence)-k+1)
	for i := 0; i <= len(sequence)-k; i++ {
		kmers = append(kmers, sequence[i:i+k])
	}
	return kmers
}

func classifyReads(readFiles []string, kmerIndex map[string]*KmerStats, k int) map[string]*SequenceReadMatch {
	readMatches := make(map[string]*SequenceReadMatch)

	// Process each read file
	for _, readFile := range readFiles {
		file, err := os.Open(readFile)
		if err != nil {
			log.Fatalf("Failed to open read file: %v", err)
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
				sequence = line
				kmers := extractKmers(sequence, k)

				// check each k-mer against the index
				for _, kmer := range kmers {
					if stats, exists := kmerIndex[kmer]; exists {
						// add matches for each organism
						for organism, count := range stats.occurrences {
							readMatches[readID].matchedOrgs[organism] += count
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

func main() {
	k := 31
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
	fmt.Println("K-mer Based Classification Report")
	fmt.Println(strings.Repeat("=", 80))

	fmt.Printf("\nBuilding k-mer index with k = %d:\n", k)
	kmerIndex, _ := buildKmerIndex(genomeFiles, k)

	fmt.Printf("\nClassifying reads.\n")
	readMatches := classifyReads(readFiles, kmerIndex, k)

	orgReadCounts := make(map[string]int)
	orgKmerCounts := make(map[string]int)
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
			orgKmerCounts[org] += match.matchedOrgs[org]
		}
	}

	// Report:
	fmt.Printf("\n1. Classification Results:\n")
	fmt.Printf("	Total sequence reads processed: %d\n", len(readMatches))
	fmt.Printf("	Matched (k-mer) reads per organism:\n")
	for _, orgName := range []string{"E. coli", "B. subtilis", "P. aeruginosa", "S. aureus", "M. tuberculosis"} {
		fmt.Printf("     %-15s: %d reads, %d k-mer matches\n", orgName, orgReadCounts[orgName], orgKmerCounts[orgName])
	}

	fmt.Printf("\n2. Match Statistics:\n")
	fmt.Printf("	Reads with unique matches: %d (%.2f%%)\n",
		uniqueMatches, float64(uniqueMatches)*100/float64(len(readMatches)))
	fmt.Printf("	Reads with multiple matches: %d (%.2f%%)\n",
		multipleMatches, float64(multipleMatches)*100/float64(len(readMatches)))
	fmt.Printf("	Reads with no matches: %d (%.2f%%)\n",
		noMatches, float64(noMatches)*100/float64(len(readMatches)))

	fmt.Println("\n" + strings.Repeat("=", 80))
}
