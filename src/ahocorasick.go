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
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			continue
		}
		line = strings.ToLower(strings.TrimSpace(line))

		current := ac.root
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

func (ac *AhoCorasick) BuildTrie(patternFilePaths []string) {
	for _, patternFilePath := range patternFilePaths {
		patternFile, err := os.Open(patternFilePath)
		if err != nil {
			fmt.Println("Error opening pattern file:", err)
			continue
		}
		defer patternFile.Close()

		scanner := bufio.NewScanner(patternFile)
		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, "@") || strings.HasPrefix(line, "+") { // TODO: check if this is correct
				continue
			}
			line = strings.ToLower(strings.TrimSpace(line))
			ac.AddPattern(line)
			fmt.Println("Added pattern:", line)
		}

		if err := scanner.Err(); err != nil {
			fmt.Println("Error reading pattern file:", err)
		}
	}

	ac.ComputeFailureLinks()
	fmt.Println("Failure links computed")
}

// TODO: implement main function and test
func main() {
	if len(os.Args) == 2 && os.Args[1] == "test" {
		fmt.Println("Ran: go run ahocorasick.go test")
		ac := NewAhoCorasick()
		ac.BuildTrie([]string{
			"test/pattern1.txt",
			"test/pattern2.txt",
			"test/pattern3.txt",
		})
		matches := ac.Search("test/text.txt")
		fmt.Println("Found matches:")
		for _, match := range matches {
			fmt.Println(match)
		}
		return
	}
}
