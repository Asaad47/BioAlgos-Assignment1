package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

type AhoCorasick struct {
	root *Node
}

type Node struct {
	children map[rune]*Node
	failLink *Node
	pattern  string
}

func NewAhoCorasick() *AhoCorasick {
	return &AhoCorasick{root: &Node{children: make(map[rune]*Node)}}
}

func (ac *AhoCorasick) AddPattern(pattern string) {
	current := ac.root
	for _, char := range pattern {
		if current.children[char] == nil {
			current.children[char] = &Node{children: make(map[rune]*Node)}
		}
		current = current.children[char]
	}
	current.pattern = pattern
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

func (ac *AhoCorasick) Search(textFilePath string) []string {
	textFile, err := os.Open(textFilePath)
	if err != nil {
		fmt.Println("Error opening text file:", err)
		return []string{}
	}
	defer textFile.Close()

	scanner := bufio.NewScanner(textFile)

	matches := []string{}
	textIndex := 0
	current := ac.root
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			continue
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
					matches = append(matches, fmt.Sprintf("(%s, %d, %d)", temp.pattern, textIndex-len(temp.pattern)+1, textIndex))
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

func (ac *AhoCorasick) BuildTrie(patternFilePath string) {
	numPatterns := 0
	patternFile, err := os.Open(patternFilePath)
	if err != nil {
		fmt.Println("Error opening pattern file:", err)
		return
	}
	defer patternFile.Close()

	scanner := bufio.NewScanner(patternFile)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "@") || strings.HasPrefix(line, "+") || strings.HasPrefix(line, "I") { // TODO: check if this is correct
			continue
		}
		line = strings.ToLower(strings.TrimSpace(line))
		ac.AddPattern(line)
		numPatterns++
	}

	if err := scanner.Err(); err != nil {
		fmt.Println("Error reading pattern file:", err)
	}
	fmt.Println("Number of patterns added:", numPatterns)

	ac.ComputeFailureLinks()
	fmt.Println("Failure links computed!")
}

// TODO: implement main function and test
func main() {
	if len(os.Args) == 2 && os.Args[1] == "test" {
		fmt.Println("Ran: go run ahocorasick.go test")
		ac := NewAhoCorasick()
		ac.BuildTrie("test/pattern.txt")
		matches := ac.Search("test/text.txt")
		fmt.Println("Found matches:")
		for _, match := range matches {
			fmt.Println(match)
		}
		fmt.Println("Number of matches:", len(matches))
		return
	} else if len(os.Args) == 2 && os.Args[1] == "all" {
		patternFilePaths := []string{
			"../data/sequence-reads/simulated_reads_no_errors_10k_R1.fastq",
			"../data/sequence-reads/simulated_reads_no_errors_10k_R2.fastq",
			// "../data/sequence-reads/simulated_reads_miseq_10k_R1.fastq",
			// "../data/sequence-reads/simulated_reads_miseq_10k_R2.fastq",
		}
		textFilePaths := []string{
			"../data/1_ecol_ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna",
			"../data/2_bsub_ncbi_dataset/ncbi_dataset/data/GCF_000009045.1/GCF_000009045.1_ASM904v1_genomic.fna",
			"../data/3_paer_ncbi_dataset/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna",
			"../data/4_saur_ncbi_dataset/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna",
			"../data/5_mtub_ncbi_dataset/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna",
		}
		for _, patternFilePath := range patternFilePaths {
			ac := NewAhoCorasick()
			ac.BuildTrie(patternFilePath)
			for _, textFilePath := range textFilePaths {
				matches := ac.Search(textFilePath)
				fmt.Println("Pattern file:", patternFilePath)
				fmt.Println("Text file:", textFilePath)
				fmt.Println("Found matches:")
				// for _, match := range matches {
				// 	fmt.Println(match)
				// }
				fmt.Println("Number of matches:", len(matches))
			}
		}
	} else if len(os.Args) == 2 && os.Args[1] == "pattern-test" {
		patternFilePath := "../data/sequence-reads/simulated_reads_no_errors_10k_R1.fastq"
		// textFilePath := "../data/1_ecol_ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
		textFilePath := "test/simulated_no_errors_text.txt" // outputs two matches (2nd and 3rd patterns)

		fmt.Println("Pattern file:", patternFilePath)
		fmt.Println("Text file:", textFilePath)
		ac := NewAhoCorasick()
		ac.BuildTrie(patternFilePath)
		matches := ac.Search(textFilePath)
		fmt.Println("Found matches:")
		for _, match := range matches {
			fmt.Println(match)
		}
		fmt.Println("Number of matches:", len(matches))
	} else if len(os.Args) == 2 && os.Args[1] == "text-test" {
		patternFilePath := "test/ecoli_pattern.txt"
		textFilePath := "../data/1_ecol_ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
		ac := NewAhoCorasick()
		ac.BuildTrie(patternFilePath)
		matches := ac.Search(textFilePath)
		fmt.Println("Found matches:")
		for _, match := range matches {
			fmt.Println(match)
		}
		fmt.Println("Number of matches:", len(matches))
	}
}
