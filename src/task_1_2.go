package main

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

type AhoCorasick struct {
	root *Node
}

type Node struct {
	children map[rune]*Node
	failLink *Node
	pattern  string
	readID   string
}

// ReadMatch holds information about which organisms a read matches
type ReadMatch struct {
	readID       string
	matchedOrgs  map[string]bool
	matchedCount map[string]int
}

func NewAhoCorasick() *AhoCorasick {
	return &AhoCorasick{root: &Node{children: make(map[rune]*Node)}}
}

func (ac *AhoCorasick) AddPattern(pattern string, readID string) {
	current := ac.root
	for _, char := range pattern {
		if current.children[char] == nil {
			current.children[char] = &Node{children: make(map[rune]*Node)}
		}
		current = current.children[char]
	}
	current.pattern = pattern
	current.readID = readID
}

func (ac *AhoCorasick) ComputeFailureLinks() {
	queue := []*Node{ac.root}
	ac.root.failLink = nil

	for len(queue) > 0 {
		current := queue[0]
		queue = queue[1:]

		for char, child := range current.children {
			queue = append(queue, child)

			if current == ac.root {
				child.failLink = ac.root
				continue
			}

			failLink := current.failLink
			for failLink != nil && failLink.children[char] == nil {
				failLink = failLink.failLink
			}

			if failLink == nil {
				child.failLink = ac.root
			} else {
				child.failLink = failLink.children[char]
			}
		}
	}
}

func (ac *AhoCorasick) Search(textFilePath string, organismName string) map[string]int {
	textFile, err := os.Open(textFilePath)
	if err != nil {
		fmt.Println("Error opening text file:", err)
		return nil
	}
	defer textFile.Close()

	scanner := bufio.NewScanner(textFile)

	matches := make(map[string]int) // maps readID to count of matches

	textIndex := 0
	current := ac.root

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			continue // skip header lines in FASTA format
		}
		line = strings.ToLower(strings.TrimSpace(line))

		for _, char := range line {
			textIndex++
			for current != nil && current.children[char] == nil {
				current = current.failLink
			}

			if current == nil {
				current = ac.root
				continue
			}

			current = current.children[char]

			temp := current
			for temp != nil {
				if temp.pattern != "" {
					// record match, incrementing counter if exists
					matches[temp.readID] = matches[temp.readID] + 1
				}
				temp = temp.failLink
			}
		}
	}

	if err := scanner.Err(); err != nil {
		fmt.Println("Error reading text file:", err)
	}

	return matches
}

func (ac *AhoCorasick) BuildTrieFromFastq(patternFilePath string) (int, []string) {
	numPatterns := 0
	readIDs := []string{}

	patternFile, err := os.Open(patternFilePath)
	if err != nil {
		fmt.Println("Error opening pattern file:", err)
		return numPatterns, readIDs
	}
	defer patternFile.Close()

	scanner := bufio.NewScanner(patternFile)
	var readID, seq string
	lineNum := 0

	for scanner.Scan() {
		line := scanner.Text()
		switch lineNum % 4 {
		case 0: // header line
			if strings.HasPrefix(line, "@") {
				readID = strings.TrimSpace(line[1:]) // remove @ and whitespace
				readIDs = append(readIDs, readID)
			}
		case 1: // sequence line
			seq = strings.ToLower(strings.TrimSpace(line))
			ac.AddPattern(seq, readID)
			numPatterns++
		}
		lineNum++
	}

	if err := scanner.Err(); err != nil {
		fmt.Println("Error reading pattern file:", err)
	}

	ac.ComputeFailureLinks()
	return numPatterns, readIDs
}

func getOrganismShortName(path string) string {
	filename := filepath.Base(path)

	// map filenames to short names
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
	readFiles := []string{
		"../data/sequence_reads/simulated_reads_no_errors_10k_R1.fastq",
		"../data/sequence_reads/simulated_reads_no_errors_10k_R2.fastq",
		"../data/sequence_reads/simulated_reads_miseq_10k_R1.fastq",
		"../data/sequence_reads/simulated_reads_miseq_10k_R2.fastq",
	}

	genomeFiles := []string{
		"../data/1_ecol_ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna",
		"../data/2_bsub_ncbi_dataset/ncbi_dataset/data/GCF_000009045.1/GCF_000009045.1_ASM904v1_genomic.fna",
		"../data/3_paer_ncbi_dataset/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna",
		"../data/4_saur_ncbi_dataset/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna",
		"../data/5_mtub_ncbi_dataset/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna",
	}

	for _, readFile := range readFiles {
		fmt.Printf("\n=== Processing %s ===\n", readFile)

		readMatches := make(map[string]*ReadMatch)

		ac := NewAhoCorasick()
		numPatterns, readIDs := ac.BuildTrieFromFastq(readFile)
		fmt.Printf("Added %d patterns from %d reads\n", numPatterns, len(readIDs))

		for _, readID := range readIDs {
			readMatches[readID] = &ReadMatch{
				readID:       readID,
				matchedOrgs:  make(map[string]bool),
				matchedCount: make(map[string]int),
			}
		}

		for _, genomeFile := range genomeFiles {
			orgName := getOrganismShortName(genomeFile)
			fmt.Printf("Searching against %s genome...\n", orgName)

			matches := ac.Search(genomeFile, orgName)

			for readID, matchInfo := range matches {
				if _, exists := readMatches[readID]; exists {
					readMatches[readID].matchedOrgs[orgName] = true
					readMatches[readID].matchedCount[orgName] = matchInfo
				}
			}
		}

		fmt.Printf("\n--- Classification Report for %s ---\n", readFile)

		orgReadCounts := make(map[string]int) // maps organism name to count of reads
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
			}
		}

		// print organism-specific counts
		fmt.Println("\nReads matching each organism:")
		for _, orgName := range []string{"E. coli", "B. subtilis", "P. aeruginosa", "S. aureus", "M. tuberculosis"} {
			fmt.Printf("%s: %d reads\n", orgName, orgReadCounts[orgName])
		}

		// print summary
		fmt.Printf("\nTotal reads: %d\n", len(readMatches))
		fmt.Printf("Reads matching exactly one organism: %d (%.2f%%)\n",
			uniqueMatches, float64(uniqueMatches)*100/float64(len(readMatches)))
		fmt.Printf("Reads matching multiple organisms: %d (%.2f%%)\n",
			multipleMatches, float64(multipleMatches)*100/float64(len(readMatches)))
		fmt.Printf("Reads with no matches: %d (%.2f%%)\n",
			noMatches, float64(noMatches)*100/float64(len(readMatches)))
	}
}
