package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strings"
)

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

func main() {
	k := 31
	genomeFiles := []string{
		"../data/1_ecol_ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna",
		"../data/2_bsub_ncbi_dataset/ncbi_dataset/data/GCF_000009045.1/GCF_000009045.1_ASM904v1_genomic.fna",
		"../data/3_paer_ncbi_dataset/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna",
		"../data/4_saur_ncbi_dataset/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna",
		"../data/5_mtub_ncbi_dataset/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna",
	}

	fmt.Println("\n" + strings.Repeat("=", 80))
	fmt.Println("K-mer Index Analysis Report")
	fmt.Println(strings.Repeat("=", 80))

	fmt.Printf("\nBuilding k-mer index with k = %d:\n", k)
	kmerIndex, totalGenomeLength := buildKmerIndex(genomeFiles, k)

	totalKmers := len(kmerIndex)
	maxTheoreticalKmers := int(math.Pow(4, float64(k))) // 4^k possible k-mers for DNA
	theoreticalKmers := totalGenomeLength - (k-1)*len(genomeFiles)

	kmersPerOrg := make(map[string]int)
	for _, stats := range kmerIndex {
		for org, count := range stats.occurrences {
			if count > 0 {
				kmersPerOrg[org]++
			}
		}
	}

	// Report:
	fmt.Printf("\n1. Data Structure Description:\n")
	fmt.Printf("	Used a map[string]*KmerStats structure.\n")

	fmt.Printf("\n2. K-mer Index Statistics:\n")
	fmt.Printf("	Total unique k-mers in index: %d\n", totalKmers)
	fmt.Printf("	K-mers per organism:\n")
	for _, orgName := range []string{"E. coli", "B. subtilis", "P. aeruginosa", "S. aureus", "M. tuberculosis"} {
		fmt.Printf("	- %-20s: %d unique k-mers\n", orgName, kmersPerOrg[orgName])
	}

	fmt.Printf("\n3. Theoretical and Discrepancy Analysis:\n")
	fmt.Printf("	The actual number of k-mers: %d\n", totalKmers)
	fmt.Printf("	The max theoretical number of k-mers (4^k): %d\n", maxTheoreticalKmers)
	fmt.Printf("	The theoretical number of k-mers (total genome length - (k-1) * number of genomes): %d\n", theoreticalKmers)

	fmt.Println("\n" + strings.Repeat("=", 80))
}
