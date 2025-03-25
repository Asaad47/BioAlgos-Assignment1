# Introduction
- [Introduction](#introduction)
  - [Setup](#setup)
- [Task 1: Metagenome Classification by String Matching](#task-1-metagenome-classification-by-string-matching)
  - [Task 1.1 (Multiple Matches)](#task-11-multiple-matches)
  - [Task 1.2 (Exact Matching)](#task-12-exact-matching)
  - [Task 1.4 (Comparison with Bioinformatics tools)](#task-14-comparison-with-bioinformatics-tools)
- [Task 2: Metagenomic classification by k-mer index](#task-2-metagenomic-classification-by-k-mer-index)
  - [Task 2.1 (Build the k-mer Index)](#task-21-build-the-k-mer-index)
  - [Task 2.2 (Implement Classification)](#task-22-implement-classification)
  - [Task 2.3 (Minimizers)](#task-23-minimizers)
- [Task 3: Real-world data and tools](#task-3-real-world-data-and-tools)
  - [Task 3.1 (Comparison)](#task-31-comparison)
    - [Brief Analysis](#brief-analysis)
  - [Task 3.2 (Real-world use case)](#task-32-real-world-use-case)
    - [Brief Analysis](#brief-analysis-1)

In this assignment, I had the opportunity to work with DNAs of 5 different organisms and sequence reads to classify them into different organisms, and also to use bioinformatics tools to classify the sequence reads. To do this classification, Aho-Corasick algorithm is used in [Task 1.2](#task-12-exact-matching) to do exact matching and compared with BLAST tool in [Task 1.4](#task-14-comparison-with-bioinformatics-tools). K-mer index and minimizer based classification are implemented in [Task 2.1](#task-21-build-the-k-mer-index), [Task 2.2](#task-22-implement-classification), and [Task 2.3](#task-23-minimizers). Finally, Kraken2 tool is used in [Task 3.1](#task-31-comparison) to classify the same sequence reads and compared to the results of the previous tasks, and also used to classify other sequence reads of human gut and wastewater samples in [Task 3.2](#task-32-real-world-use-case).

## Setup

To reproduce the results for all tasks, install the genome and reads data to the `data/` directory.

Save given sequence reads to `data/sequence_reads/` directory. Download genome data from NCBI and rename each genome ncbi dataset directory to its corresponding name of the following:
- `1_ecol_ncbi_dataset/`
- `2_bsub_ncbi_dataset/`
- `3_paer_ncbi_dataset/`
- `4_saur_ncbi_dataset/`
- `5_mtub_ncbi_dataset/`

In the end, the `data/` directory should look like this:

```bash
data/
├── 1_ecol_ncbi_dataset/
├── 2_bsub_ncbi_dataset/
├── 3_paer_ncbi_dataset/
├── 4_saur_ncbi_dataset/
├── 5_mtub_ncbi_dataset/
├── sequence_reads/
│   ├── simulated_reads_no_errors_10k_R1.fastq
│   ├── simulated_reads_no_errors_10k_R2.fastq
│   ├── simulated_reads_miseq_10k_R1.fastq
│   ├── simulated_reads_miseq_10k_R2.fastq
```

For each task, a command is described to reproduce the results by running the program in the `src/` directory.

# Task 1: Metagenome Classification by String Matching

For this task, Aho-Corasick algorithm is implemented in `task_1_2.go` and BLAST tool is used to compare the results in [Task 1.4](#task-14-comparison-with-bioinformatics-tools).

The relevant files for this task are:
- `src/task_1_2.go`
- `src/run_blast_analysis.sh`

## Task 1.1 (Multiple Matches)

- Proposed method for handling multiple matches
  - Track all matches with counts:
    - For each read, maintain a map of organisms to match counts
    - Record both:
      - Number of distinct matches (how many organisms matched)
      - Number of matching positions for each organism
  - Classification rules:
    - If a read matches exactly one organism: classify as that organism
    - If a read matches multiple organisms:
      - Compare the number of matching positions for each organism
      - If one organism has significantly more matches: classify as that organism
      - Otherwise: mark as ambiguous/multiple match
  - Biological justification:
    - Multiple matches often indicate conserved regions between organisms
    - Some bacterial species share significant portions of their genomes due to, for example, common ancestry or preservation of essential genes.
  - Computational efficiency
    - Space efficiency:
      - Using a map structure for tracking matches requires O(m) space per read, where m is the number of matching organisms
      - From the implementation results in [Task 1.2](#task-12-exact-matching), we see that most reads match 0 or 1 organism, making this space overhead minimal
    - Time efficiency:
      - The approach requires only one pass through the matches
      - No need for complex post-processing or multiple comparisons
      - Can be integrated directly into the existing Aho-Corasick implementation shown in `task_1_2.go`

## Task 1.2 (Exact Matching)

Used Aho-Corasick algorithm to find exact matches.

- Total reads processed: 20,000
- Total reads matching each organism:
  - E. coli: 3,269
  - B. subtilis: 553
  - P. aeruginosa: 538
  - S. aureus: 549
  - M. tuberculosis: 548
- Reads matching exactly one organism: 5,437 (27.19%)
- Reads matching multiple organisms: 10 (0.05%)
- Reads with no matches: 14,553 (72.76%)

<br>

- Run `/usr/bin/time -l go run task_1_2.go` to get the results for all files.


Output for this task:

```bash
=== Processing ../data/sequence_reads/simulated_reads_no_errors_10k_R1.fastq ===
Added 5000 patterns from 5000 reads
Searching against E. coli genome...
Searching against B. subtilis genome...
Searching against P. aeruginosa genome...
Searching against S. aureus genome...
Searching against M. tuberculosis genome...

--- Classification Report for ../data/sequence_reads/simulated_reads_no_errors_10k_R1.fastq ---

Reads matching each organism:
E. coli: 3000 reads
B. subtilis: 500 reads
P. aeruginosa: 500 reads
S. aureus: 504 reads
M. tuberculosis: 500 reads

Total reads: 5000
Reads matching exactly one organism: 4996 (99.92%)
Reads matching multiple organisms: 4 (0.08%)
Reads with no matches: 0 (0.00%)

=== Processing ../data/sequence_reads/simulated_reads_no_errors_10k_R2.fastq ===
Added 5000 patterns from 5000 reads
Searching against E. coli genome...
Searching against B. subtilis genome...
Searching against P. aeruginosa genome...
Searching against S. aureus genome...
Searching against M. tuberculosis genome...

--- Classification Report for ../data/sequence_reads/simulated_reads_no_errors_10k_R2.fastq ---

Reads matching each organism:
E. coli: 53 reads
B. subtilis: 5 reads
P. aeruginosa: 3 reads
S. aureus: 8 reads
M. tuberculosis: 6 reads

Total reads: 5000
Reads matching exactly one organism: 65 (1.30%)
Reads matching multiple organisms: 5 (0.10%)
Reads with no matches: 4930 (98.60%)

=== Processing ../data/sequence_reads/simulated_reads_miseq_10k_R1.fastq ===
Added 5000 patterns from 5000 reads
Searching against E. coli genome...
Searching against B. subtilis genome...
Searching against P. aeruginosa genome...
Searching against S. aureus genome...
Searching against M. tuberculosis genome...

--- Classification Report for ../data/sequence_reads/simulated_reads_miseq_10k_R1.fastq ---

Reads matching each organism:
E. coli: 215 reads
B. subtilis: 48 reads
P. aeruginosa: 35 reads
S. aureus: 36 reads
M. tuberculosis: 42 reads

Total reads: 5000
Reads matching exactly one organism: 376 (7.52%)
Reads matching multiple organisms: 0 (0.00%)
Reads with no matches: 4624 (92.48%)

=== Processing ../data/sequence_reads/simulated_reads_miseq_10k_R2.fastq ===
Added 5000 patterns from 5000 reads
Searching against E. coli genome...
Searching against B. subtilis genome...
Searching against P. aeruginosa genome...
Searching against S. aureus genome...
Searching against M. tuberculosis genome...

--- Classification Report for ../data/sequence_reads/simulated_reads_miseq_10k_R2.fastq ---

Reads matching each organism:
E. coli: 1 reads
B. subtilis: 0 reads
P. aeruginosa: 0 reads
S. aureus: 1 reads
M. tuberculosis: 0 reads

Total reads: 5000
Reads matching exactly one organism: 0 (0.00%)
Reads matching multiple organisms: 1 (0.02%)
Reads with no matches: 4999 (99.98%)
        5.24 real         5.43 user         0.23 sys
           449642496  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
               46141  page reclaims
                 102  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                 500  signals received
                 263  voluntary context switches
                7103  involuntary context switches
           819593799  instructions retired
           319308165  cycles elapsed
            13616384  peak memory footprint
```


## Task 1.4 (Comparison with Bioinformatics tools)

- Execution time:
  - Not much difference in execution time between my implementation in Task 1.2 and `blastn` tool. They both took around 5 seconds to run. However, my implementation only considers exact matches while BLAST considers approximate matches.
- Memory usage:
  - My implementation's maximum resident set size was around 429MB while BLAST's was around 184MB.
- The resulting counts are quite different.
  - The BLAST tool uses a more sophisticated algorithm to find approximate matches and returns a larger number of matches.
  - BLAST considers statistical significance of matches and requires parameters of `-max_target_seqs` and `-evalue` to control the number of matches returned.

<br>

`run_blast_analysis.sh` is a script that runs the BLAST tool on the sequence reads and generates a summary of the results. It first concatenates the genome files into one file and builds a BLAST database from it. Then, it converts the sequence reads to FASTA format and runs BLAST on each read file. Finally, it generates a summary of the results.

`run_blast_analysis.sh` assumes that the genome files are in the `data/` directory and to be run in the `src/` directory. Also, it uses `seqtk` to convert the FASTQ files to FASTA format, so a total of 3 commands are needed to be installed: `makeblastdb`, `blastn`, and, `seqtk`.

- Run `sh run_blast_analysis.sh` to get the results for this task.

Output for this task:

```bash
Concatenating genome files into five_genomes.fna...
Building BLAST database...


Building a new DB, current time: 03/24/2025 03:29:37
New DB name:   ./src/reference_db
New DB title:  ../data/five_genomes.fna
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 3000000000B
Adding sequences from FASTA; added 5 sequences in 0.076303 seconds.


Converting FASTQ files to FASTA format...
Running BLAST for each read file...
Output file ../data/sequence_reads/simulated_reads_no_errors_10k_R1_blast_results.txt not found. Creating empty file...
Processing ../data/sequence_reads/simulated_reads_no_errors_10k_R1.fasta...
        1.34 real         0.35 user         0.04 sys
           192512000  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
               11576  page reclaims
                 887  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                  11  messages sent
                  23  messages received
                   0  signals received
                   6  voluntary context switches
                 376  involuntary context switches
          4843612200  instructions retired
          1380976406  cycles elapsed
           118376064  peak memory footprint
Output file ../data/sequence_reads/simulated_reads_no_errors_10k_R2_blast_results.txt not found. Creating empty file...
Processing ../data/sequence_reads/simulated_reads_no_errors_10k_R2.fasta...
        1.31 real         0.38 user         0.03 sys
           182747136  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
               11711  page reclaims
                  50  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                  11  messages sent
                  23  messages received
                   0  signals received
                   7  voluntary context switches
                  88  involuntary context switches
          4907632152  instructions retired
          1357133345  cycles elapsed
           116295232  peak memory footprint
Output file ../data/sequence_reads/simulated_reads_miseq_10k_R1_blast_results.txt not found. Creating empty file...
Processing ../data/sequence_reads/simulated_reads_miseq_10k_R1.fasta...
        1.51 real         0.51 user         0.04 sys
           187269120  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
               12385  page reclaims
                  37  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                  11  messages sent
                  23  messages received
                   0  signals received
                   6  voluntary context switches
                 109  involuntary context switches
          6427546323  instructions retired
          1824478843  cycles elapsed
           112985600  peak memory footprint
Output file ../data/sequence_reads/simulated_reads_miseq_10k_R2_blast_results.txt not found. Creating empty file...
Processing ../data/sequence_reads/simulated_reads_miseq_10k_R2.fasta...
        1.52 real         0.50 user         0.03 sys
           188432384  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
               11983  page reclaims
                  37  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                  11  messages sent
                  23  messages received
                   0  signals received
                  11  voluntary context switches
                 130  involuntary context switches
          6303387720  instructions retired
          1793370325  cycles elapsed
           113887104  peak memory footprint
Generating classification summary...
Summary for ../data/sequence_reads/simulated_reads_no_errors_10k_R1_blast_results.txt:
3997 NC_000913.3
 683 NC_000962.3
 635 NC_000964.3
 625 NC_002516.2
 656 NC_007795.1
------------------------------------------------
Summary for ../data/sequence_reads/simulated_reads_no_errors_10k_R2_blast_results.txt:
4095 NC_000913.3
 652 NC_000962.3
 691 NC_000964.3
 626 NC_002516.2
 653 NC_007795.1
------------------------------------------------
Summary for ../data/sequence_reads/simulated_reads_miseq_10k_R1_blast_results.txt:
5286 NC_000913.3
 634 NC_000962.3
 696 NC_000964.3
 656 NC_002516.2
 740 NC_007795.1
------------------------------------------------
Summary for ../data/sequence_reads/simulated_reads_miseq_10k_R2_blast_results.txt:
4461 NC_000913.3
 637 NC_000962.3
 750 NC_000964.3
 683 NC_002516.2
 756 NC_007795.1
------------------------------------------------
Analysis complete!
```

# Task 2: Metagenomic classification by k-mer index

For this task, k-mer index is implemented in [Task 2.1](#task-21-build-the-k-mer-index), k-mer based classification is implemented in [Task 2.2](#task-22-implement-classification), and minimizer based classification is implemented in [Task 2.3](#task-23-minimizers).

The relevant files for this task are:
- `src/task_2_1.go`
- `src/task_2_2.go`
- `src/task_2_3.go`

## Task 2.1 (Build the k-mer Index)

1. Data Structure Description:
     - Used a `map[string]*KmerStats` structure where:
       - Key: k-mer sequence
       - Value: `KmerStats` struct containing:
         - `occurrences`: map of organism names to their counts
         - `totalCount`: total occurrences across all genomes
     - Justification: This structure allows O(1) lookup of k-mer occurrences, efficient tracking across multiple genomes, and works effeciently in our case with small genome sizes.

2. K-mer Index Statistics:
     - Total unique k-mers in index: 14,128,469
     - K-mers per organism:
       - E. coli        : 2,923,616 unique k-mers
       - B. subtilis    : 2,662,309 unique k-mers
       - P. aeruginosa  : 3,975,633 unique k-mers
       - S. aureus      : 1,784,492 unique k-mers
       - M. tuberculosis: 2,786,106 unique k-mers

3. Theoretical Analysis:
     - Max theoretical number of possible k-mers generally (k = 31): 4,611,686,018,427,387,904
       - This is calculated as 4^k, where 4 represents the possible DNA bases (A,T,C,G).
     - Theoretical number of possible k-mers (k = 31): 22,354,405
       - This is calculated as the total genome length minus (k-1) times the number of genomes.

4. Discrepancy Analysis:
     - The actual number of k-mers (14,128,469) is much smaller than the max theoretical (4,611,686,018,427,387,904) because:
       - DNA sequences are not random and they follow biological patterns.
       - Many theoretical combinations never appear in real DNA due to biological constraints.
       - Genome sequences have significant redundancy and patterns.
     - Also, the actual number of k-mers (14,128,469) is significantly smaller than the possible theoretical limit (22,354,405):
       - This is because the genome sequences have significant redundancy and patterns.

<br>

- Run `/usr/bin/time -l go run task_2_1.go` to get the results for this task.

Output for this task:

```bash
================================================================================
K-mer Index Analysis Report
================================================================================

Building k-mer index with k = 31:
Processing genome: E. coli
Genome length (E. coli): 4641652
Processing genome: B. subtilis
Genome length (B. subtilis): 4215606
Processing genome: P. aeruginosa
Genome length (P. aeruginosa): 6264404
Processing genome: S. aureus
Genome length (S. aureus): 2821361
Processing genome: M. tuberculosis
Genome length (M. tuberculosis): 4411532

1. Data Structure Description:
        Used a map[string]*KmerStats structure.

2. K-mer Index Statistics:
        Total unique k-mers in index: 14128469
        K-mers per organism:
        E. coli        : 2923616 unique k-mers
        B. subtilis    : 2662309 unique k-mers
        P. aeruginosa  : 3975633 unique k-mers
        S. aureus      : 1784492 unique k-mers
        M. tuberculosis: 2786106 unique k-mers

3. Theoretical and Discrepancy Analysis:
        The actual number of k-mers: 14128469
        The max theoretical number of k-mers (4^k): 4611686018427387904
        The theoretical number of k-mers (total genome length - (k-1) * number of genomes): 22354405

================================================================================
       18.02 real        35.57 user         6.33 sys
          5141397504  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
             1206375  page reclaims
                2216  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                2576  signals received
                1437  voluntary context switches
              157227  involuntary context switches
           813768512  instructions retired
           321841346  cycles elapsed
            13714688  peak memory footprint
```

## Task 2.2 (Implement Classification)

1. Classification Results:
    - Total reads processed: 20,000
    - Matched reads per organism:
      - E. coli        : 6,190 reads, 668,055 k-mer matches
      - B. subtilis    : 1,056 reads, 113,283 k-mer matches
      - P. aeruginosa  : 1,100 reads, 105,376 k-mer matches
      - S. aureus      : 1,091 reads, 111,586 k-mer matches
      - M. tuberculosis: 1,032 reads, 111,353 k-mer matches
    - Reads with unique matches: 10,110 (50.55%)
    - Reads with multiple matches: 144 (0.72%)
    - Reads with no matches: 9,746 (48.73%)

2. Comparison with String Matching Approach:
     - The k-mer based approach shows different results because:
       - A single read can match multiple times in different positions
       - K-mer matches can overlap and be counted multiple times
       - String matching looks for exact sequence matches, while k-mer matching allows for partial matches

3. Handling Multiple Matches:
    - Two numbers are reported for each organism:
      - The first number is the number of reads that have a k-mer match to the organism.
      - The second number is the number of k-mer matches for the organism.
      - The second number is higher than the first number because a single read can have multiple k-mer matches to the same organism.
    - The reported unique matches are based on the first number. This is because read matches are counted as a match if they have a k-mer match to the organism.


<br>

- Run `/usr/bin/time -l go run task_2_2.go` to get the results for this task.

Output for this task:

```bash
================================================================================
K-mer Based Classification Report
================================================================================

Building k-mer index with k = 31:
Processing genome: E. coli
Genome length (E. coli): 4641652
Processing genome: B. subtilis
Genome length (B. subtilis): 4215606
Processing genome: P. aeruginosa
Genome length (P. aeruginosa): 6264404
Processing genome: S. aureus
Genome length (S. aureus): 2821361
Processing genome: M. tuberculosis
Genome length (M. tuberculosis): 4411532

Classifying reads.

1. Classification Results:
        Total sequence reads processed: 20000
        Matched (k-mer) reads per organism:
     E. coli        : 6190 reads, 668055 k-mer matches
     B. subtilis    : 1056 reads, 113283 k-mer matches
     P. aeruginosa  : 1100 reads, 105376 k-mer matches
     S. aureus      : 1091 reads, 111586 k-mer matches
     M. tuberculosis: 1032 reads, 111353 k-mer matches

2. Match Statistics:
        Reads with unique matches: 10110 (50.55%)
        Reads with multiple matches: 144 (0.72%)
        Reads with no matches: 9746 (48.73%)

================================================================================
       11.39 real        32.05 user         3.51 sys
          5394333696  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
              786747  page reclaims
                2205  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                2155  signals received
                1326  voluntary context switches
              107546  involuntary context switches
           809544219  instructions retired
           317889284  cycles elapsed
            13534528  peak memory footprint
```

## Task 2.3 (Minimizers)

1. Classification Results:
    - Total sequence reads processed: 20,000
    - Matched (minimizer) reads per organism:
      - E. coli        : 6,177 reads, 3,972,502 minimizer matches
      - B. subtilis    : 1,041 reads, 672,715 minimizer matches
      - P. aeruginosa  : 1,084 reads, 632,656 minimizer matches
      - S. aureus      : 1,084 reads, 669,813 minimizer matches
      - M. tuberculosis: 1,024 reads, 664,123 minimizer matches

2. Match Statistics:
    - Reads with unique matches: 10,106 (50.53%)
    - Reads with multiple matches: 126 (0.63%)
    - Reads with no matches: 9,768 (48.84%)

3. Comparison with K-mer Based Approach:
    - The minimizer based approach shows different results because:
      - A single read can match multiple times in different positions
      - Minimizer matching allows for partial matches
    - The classification accuracy results are not significantly different between the two approaches.
      - Minimzer based appraoch has a slightly lower matching rate (51.16%) than the k-mer based approach (51.27%).

4. Memory Consumption Comparison:
    - The minimizer based approach has a significant reduction in memory consumption compared to the k-mer based approach.
      - Minimzer based appraoch has a maximum resident set size of 1,753,595,904 bytes (~1.63GB) compared to the k-mer based approach's 5,394,333,696 bytes (~5.02GB).


<br>

- Run `/usr/bin/time -l go run task_2_3.go` to get the results for this task.

Output for this task:

```bash
/usr/bin/time -l go run task_2_3.go

================================================================================
Minimizer-Based Classification Report
================================================================================

Building minimizer index with k=31 and w=10
Processing genome: E. coli
Genome length (E. coli): 4641652
Processing genome: B. subtilis
Genome length (B. subtilis): 4215606
Processing genome: P. aeruginosa
Genome length (P. aeruginosa): 6264404
Processing genome: S. aureus
Genome length (S. aureus): 2821361
Processing genome: M. tuberculosis
Genome length (M. tuberculosis): 4411532

Classifying reads using minimizers.

1. Classification Results:
    Total sequence reads processed: 20000
    Matched (minimizer) reads per organism:
     E. coli        : 6177 reads, 3972502 minimizer matches
     B. subtilis    : 1041 reads, 672715 minimizer matches
     P. aeruginosa  : 1084 reads, 632656 minimizer matches
     S. aureus      : 1084 reads, 669813 minimizer matches
     M. tuberculosis: 1024 reads, 664123 minimizer matches

2. Match Statistics:
    Reads with unique matches: 10106 (50.53%)
    Reads with multiple matches: 126 (0.63%)
    Reads with no matches: 9768 (48.84%)

================================================================================
        7.80 real        19.21 user         0.40 sys
          1753595904  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
              124256  page reclaims
                 105  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                1917  signals received
                 951  voluntary context switches
               14127  involuntary context switches
           818889034  instructions retired
           314528087  cycles elapsed
            13845824  peak memory footprint
```


# Task 3: Real-world data and tools

For this task, Kraken2 tool is used in Ibex cluster to classify the given sequence reads in [Task 3.1](#task-31-comparison) and used it to classify the gut and wastewater sequence reads in [Task 3.2](#task-32-real-world-use-case).

The relevant files for this task are:
- `sra_download.sh`
- `src/group_species.py`
- `results/combined_summary_grouped.txt`
- `results/combined_summary.txt`
- `results/simulated_reads_miseq_10k_R1_report.txt`
- `results/simulated_reads_miseq_10k_R2_report.txt`
- `results/simulated_reads_no_errors_10k_R1_report.txt`
- `results/simulated_reads_no_errors_10k_R2_report.txt`

<br>

<!-- Downloaded and installed Kraken2 using the following commands:
```bash
wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.14.tar.gz

tar -xzf v2.14.tar.gz

cd kraken2-2.14

./install_kraken2.sh ../kraken2-commands

cd ../kraken2-commands

echo 'export PATH=$PWD:$PATH' >> ~/.bashrc
source ~/.bashrc
``` -->

Kraken2 was installed in Ibex cluster. To check if it is installed correctly, run the following command:
```bash
kraken2 --version
```
it returned the following output:
```bash
Kraken version 2.1.3
Copyright 2013-2023, Derrick Wood (dwood@cs.jhu.edu)
```

## Task 3.1 (Comparison)

Sequence of commands to prepare custom database for Kraken2:
```bash
mkdir -p kraken2_custom_db/library

pip install ncbi-genome-download

echo -e "GCF_000005845.2\nGCF_000009045.1\nGCF_000006765.1\nGCF_000013425.1\nGCF_000195955.2" > accessions.txt

ncbi-genome-download bacteria -A accessions.txt --formats fasta --output-folder library/

mv library/refseq/bacteria/*/*.fna.gz library/
gunzip library/*.fna.gz

mv library/GCF_000005845.2_ASM584v2_genomic.fna library/511145.fna
mv library/GCF_000009045.1_ASM904v1_genomic.fna library/224308.fna
mv library/GCF_000006765.1_ASM676v1_genomic.fna library/208964.fna
mv library/GCF_000013425.1_ASM1342v1_genomic.fna library/93061.fna
mv library/GCF_000195955.2_ASM19595v2_genomic.fna library/83332.fna

kraken2-build --download-taxonomy --db kraken2_custom_db
for file in library/*.fna; do
    kraken2-build --add-to-library "$file" --db kraken2_custom_db --no-masking  # --no-masking is used to avoid masking the reads and to finish the process faster
done

kraken2-build --build --db kraken2_custom_db
```

- Run the following command to verify the database:
```bash
>>> kraken2-inspect --db kraken2_custom_db

# Database options: nucleotide db, k = 35, l = 31
# Spaced mask = 11111111111111111111111111111111110011001100110011001100110011
# Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101
# Total taxonomy nodes: 41
# Table size: 7176774
# Table capacity: 10324845
# Min clear hash value = 0
100.00  7176774 0       R       1       root
100.00  7176774 0       R1      131567    cellular organisms
100.00  7176774 652     R2      2           Bacteria
 51.33  3684078 40      K       1783272       Bacillati
 31.66  2272323 0       P       1239            Bacillota
 31.66  2272323 0       C       91061             Bacilli
 31.66  2272323 482     O       1385                Bacillales
 18.95  1359651 0       F       186817                Bacillaceae
 18.95  1359651 0       G       1386                    Bacillus
 18.95  1359651 0       G1      653685                    Bacillus subtilis group
 18.95  1359651 0       S       1423                        Bacillus subtilis
 18.95  1359651 0       S1      135461                        Bacillus subtilis subsp. subtilis
 18.95  1359651 1359651 S2      224308                          Bacillus subtilis subsp. subtilis str. 168
 12.71  912190  0       F       90964                 Staphylococcaceae
 12.71  912190  0       G       1279                    Staphylococcus
 12.71  912190  0       S       1280                      Staphylococcus aureus
 12.71  912190  912190  S1      93061                       Staphylococcus aureus subsp. aureus NCTC 8325
 19.67  1411715 0       P       201174          Actinomycetota
 19.67  1411715 0       C       1760              Actinomycetes
 19.67  1411715 0       O       85007               Mycobacteriales
 19.67  1411715 0       F       1762                  Mycobacteriaceae
 19.67  1411715 0       G       1763                    Mycobacterium
 19.67  1411715 0       G1      77643                     Mycobacterium tuberculosis complex
 19.67  1411715 0       S       1773                        Mycobacterium tuberculosis
 19.67  1411715 1411715 S1      83332                         Mycobacterium tuberculosis H37Rv
 48.66  3492044 0       K       3379134       Pseudomonadati
 48.66  3492044 0       P       1224            Pseudomonadota
 48.66  3492044 506     C       1236              Gammaproteobacteria
 27.95  2005625 0       O       72274               Pseudomonadales
 27.95  2005625 0       F       135621                Pseudomonadaceae
 27.95  2005625 0       G       286                     Pseudomonas
 27.95  2005625 0       G1      136841                    Pseudomonas aeruginosa group
 27.95  2005625 0       S       287                         Pseudomonas aeruginosa
 27.95  2005625 2005625 S1      208964                        Pseudomonas aeruginosa PAO1
 20.70  1485913 0       O       91347               Enterobacterales
 20.70  1485913 0       F       543                   Enterobacteriaceae
 20.70  1485913 0       G       561                     Escherichia
 20.70  1485913 0       S       562                       Escherichia coli
 20.70  1485913 0       S1      83333                       Escherichia coli K-12
 20.70  1485913 1485913 S2      511145                        Escherichia coli str. K-12 substr. MG1655
```
- Have sequence reads stored in `seq_reads/` directory and run the following command to classify the reads (copies of the report files can be found in `results/` directory):
```bash
for file in seq_reads/*.fastq; do
    kraken2 --db kraken2_custom_db/kraken2_custom_db --report "${file%.fastq}_report.txt" --output "${file%.fastq}_output.txt" --threads 4 "$file"
done
```
Output:
```bash
Loading database information... done.
5000 sequences (1.50 Mbp) processed in 0.155s (1938.8 Kseq/m, 583.57 Mbp/m).
  5000 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
Loading database information... done.
5000 sequences (1.50 Mbp) processed in 0.236s (1270.0 Kseq/m, 382.27 Mbp/m).
  5000 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
Loading database information... done.
5000 sequences (0.62 Mbp) processed in 0.080s (3739.0 Kseq/m, 467.38 Mbp/m).
  5000 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
Loading database information... done.
5000 sequences (0.62 Mbp) processed in 0.062s (4860.6 Kseq/m, 607.57 Mbp/m).
  5000 sequences classified (100.00%)
  0 sequences unclassified (0.00%)
```

- Run the following command to generate the species abundance summary:
```bash
for file in seq_reads/*_report.txt; do
  echo "=== Species Abundance Summary for $file ==="
  awk '$4 ~ /^S/ {printf "%s\t%s\t%s%%\n", $6, $2, $1}' "$file"
  echo ""
done
```
Output:
```bash
=== Species Abundance Summary for seq_reads/simulated_reads_miseq_10k_R1_report.txt ===
Escherichia     2991    59.82%
Escherichia     2991    59.82%
Escherichia     2991    59.82%
Pseudomonas     500     10.00%
Pseudomonas     500     10.00%
Staphylococcus  500     10.00%
Staphylococcus  500     10.00%
Bacillus        500     10.00%
Bacillus        500     10.00%
Bacillus        500     10.00%
Mycobacterium   500     10.00%
Mycobacterium   500     10.00%

=== Species Abundance Summary for seq_reads/simulated_reads_miseq_10k_R2_report.txt ===
Escherichia     2993    59.86%
Escherichia     2993    59.86%
Escherichia     2993    59.86%
Pseudomonas     500     10.00%
Pseudomonas     500     10.00%
Staphylococcus  500     10.00%
Staphylococcus  500     10.00%
Bacillus        500     10.00%
Bacillus        500     10.00%
Bacillus        500     10.00%
Mycobacterium   500     10.00%
Mycobacterium   500     10.00%

=== Species Abundance Summary for seq_reads/simulated_reads_no_errors_10k_R1_report.txt ===
Escherichia     2990    59.80%
Escherichia     2990    59.80%
Escherichia     2990    59.80%
Pseudomonas     500     10.00%
Pseudomonas     500     10.00%
Staphylococcus  500     10.00%
Staphylococcus  500     10.00%
Bacillus        499     9.98%
Bacillus        499     9.98%
Bacillus        499     9.98%
Mycobacterium   500     10.00%
Mycobacterium   500     10.00%

=== Species Abundance Summary for seq_reads/simulated_reads_no_errors_10k_R2_report.txt ===
Escherichia     2993    59.86%
Escherichia     2993    59.86%
Escherichia     2993    59.86%
Pseudomonas     500     10.00%
Pseudomonas     500     10.00%
Staphylococcus  500     10.00%
Staphylococcus  500     10.00%
Bacillus        500     10.00%
Bacillus        500     10.00%
Bacillus        500     10.00%
Mycobacterium   500     10.00%
Mycobacterium   500     10.00%
```

### Brief Analysis

- The abundance estimates are identical in each file, with no apparent variation between samples.
- All reads are classified — no unclassified or ambiguous reads appear in the report.
- The species abundance summary shows different results to other approaches of Task 1 and 2, where here the abundance estimates are similar between all reads and only matching the `simulated_reads_no_errors_10K_R1.fastq` classification. But in previous tasks, the abundance estimates were different between different reads.
- The Kraken2 classification is faster than the other approaches of Task 1 and 2 and uses less memory. 
- Considering `simulated_reads_no_errors_10K_R1.fastq` file, the impact of the minimizer scheme did not affect the classification results, and this might be due to the fact that the reads are simulated and do not contain any errors. Also, the genomes used are not closely related to each other. 
- In general, using full k-mer index would be more accurate than using minimizer index as the full index captures more low frequency k-mers.
- So, the two methods differ in a trade-off between computational resources and accuracy. The minimizer index is faster and uses less memory but at the cost of accuracy while the full index is more accurate but slower and has higher memory usage.

## Task 3.2 (Real-world use case)

Sequence of commands to prepare custom database for Kraken2:

- Download and extract the Standard-8 database:
```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20241228.tar.gz
tar -xvzf k2_standard_08gb_20241228.tar.gz
```
- Move all extracted files to `standard8/` directory.
  - To verify that the database is correctly configured, run `kraken2-inspect --db standard8` command.

- Download the sequence reads from the SRA database using `sra_download.sh` script.
  - The script downloads the reads from `ebi.ac.uk` api as I was not able to download the reads using `fastq-dump` command.
- Move files to `SRR_reads/` directory, and unzip them using `gunzip` command.

- Run the following command to classify the reads:
```bash
for fq in SRR_reads/*.fastq; do
  kraken2 --db standard8 --report "${fq%.fastq}_report.txt" --output "${fq%.fastq}_output.txt" --threads 4 "$fq"
done
```
Output:
```bash
Loading database information... done.
19515475 sequences (2926.03 Mbp) processed in 72.718s (16102.2 Kseq/m, 2414.27 Mbp/m).
  8619706 sequences classified (44.17%)
  10895769 sequences unclassified (55.83%)
Loading database information... done.
19515475 sequences (2924.44 Mbp) processed in 72.638s (16120.0 Kseq/m, 2415.62 Mbp/m).
  8479389 sequences classified (43.45%)
  11036086 sequences unclassified (56.55%)
Loading database information... done.
16558897 sequences (2482.83 Mbp) processed in 62.079s (16004.3 Kseq/m, 2399.68 Mbp/m).
  11953045 sequences classified (72.19%)
  4605852 sequences unclassified (27.81%)
Loading database information... done.
16558897 sequences (2480.85 Mbp) processed in 61.863s (16060.3 Kseq/m, 2406.15 Mbp/m).
  11616117 sequences classified (70.15%)
  4942780 sequences unclassified (29.85%)
Loading database information... done.
19795127 sequences (2968.08 Mbp) processed in 74.023s (16045.1 Kseq/m, 2405.79 Mbp/m).
  10503634 sequences classified (53.06%)
  9291493 sequences unclassified (46.94%)
Loading database information... done.
19795127 sequences (2965.80 Mbp) processed in 73.697s (16116.1 Kseq/m, 2414.58 Mbp/m).
  10214819 sequences classified (51.60%)
  9580308 sequences unclassified (48.40%)
Loading database information... done.
21457804 sequences (3212.91 Mbp) processed in 78.021s (16501.6 Kseq/m, 2470.81 Mbp/m).
  6372221 sequences classified (29.70%)
  15085583 sequences unclassified (70.30%)
Loading database information... done.
21457804 sequences (3205.17 Mbp) processed in 77.202s (16676.6 Kseq/m, 2491.00 Mbp/m).
  6165002 sequences classified (28.73%)
  15292802 sequences unclassified (71.27%)
Loading database information... done.
19866883 sequences (2976.50 Mbp) processed in 73.742s (16164.7 Kseq/m, 2421.84 Mbp/m).
  9853808 sequences classified (49.60%)
  10013075 sequences unclassified (50.40%)
Loading database information... done.
19866883 sequences (2970.16 Mbp) processed in 74.103s (16085.9 Kseq/m, 2404.89 Mbp/m).
  9535786 sequences classified (48.00%)
  10331097 sequences unclassified (52.00%)
Loading database information... done.
877209 sequences (132.46 Mbp) processed in 2.562s (20543.7 Kseq/m, 3102.10 Mbp/m).
  479877 sequences classified (54.70%)
  397332 sequences unclassified (45.30%)
Loading database information... done.
877209 sequences (132.46 Mbp) processed in 2.634s (19984.6 Kseq/m, 3017.68 Mbp/m).
  547789 sequences classified (62.45%)
  329420 sequences unclassified (37.55%)
Loading database information... done.
2 sequences (0.00 Mbp) processed in 0.026s (4.6 Kseq/m, 0.70 Mbp/m).
  1 sequences classified (50.00%)
  1 sequences unclassified (50.00%)
Loading database information... done.
1113095 sequences (168.08 Mbp) processed in 3.142s (21257.3 Kseq/m, 3209.86 Mbp/m).
  441313 sequences classified (39.65%)
  671782 sequences unclassified (60.35%)
Loading database information... done.
1113095 sequences (168.08 Mbp) processed in 3.349s (19944.6 Kseq/m, 3011.64 Mbp/m).
  435173 sequences classified (39.10%)
  677922 sequences unclassified (60.90%)
Loading database information... done.
4 sequences (0.00 Mbp) processed in 0.047s (5.2 Kseq/m, 0.78 Mbp/m).
  0 sequences classified (0.00%)
  4 sequences unclassified (100.00%)
Loading database information... done.
523931 sequences (79.11 Mbp) processed in 1.553s (20246.5 Kseq/m, 3057.23 Mbp/m).
  158405 sequences classified (30.23%)
  365526 sequences unclassified (69.77%)
Loading database information... done.
523931 sequences (79.11 Mbp) processed in 1.520s (20682.6 Kseq/m, 3123.07 Mbp/m).
  176571 sequences classified (33.70%)
  347360 sequences unclassified (66.30%)
Loading database information... done.
181195 sequences (27.36 Mbp) processed in 0.654s (16625.6 Kseq/m, 2510.47 Mbp/m).
  115033 sequences classified (63.49%)
  66162 sequences unclassified (36.51%)
Loading database information... done.
181195 sequences (27.36 Mbp) processed in 0.648s (16786.7 Kseq/m, 2534.79 Mbp/m).
  127914 sequences classified (70.59%)
  53281 sequences unclassified (29.41%)
Loading database information... done.
546485 sequences (79.50 Mbp) processed in 1.735s (18896.9 Kseq/m, 2749.09 Mbp/m).
  368317 sequences classified (67.40%)
  178168 sequences unclassified (32.60%)
Loading database information... done.
546485 sequences (79.56 Mbp) processed in 1.850s (17728.6 Kseq/m, 2580.96 Mbp/m).
  414350 sequences classified (75.82%)
  132135 sequences unclassified (24.18%)
```

- Run the following command to generate the species abundance summary (A copy of this summary can be found in `results/combined_summary.txt`):
```bash
( for report in SRR_reads/*_report.txt; do echo "=== Summary for $report ==="; awk '$4 ~ /^S/ {printf "%s\t%s\t%s%%\n", $6, $2, $1}' "$report"; echo ""; done ) > combined_summary.txt

```

- Run the following command to group the species and print summary results (A copy of this summary can be found in `results/combined_summary_grouped.txt`):
  - The script combines the species values of the same species and prints the top species for each sample.
```bash
python group_species.py
```
  The output is:
```bash
SRR_reads/SRR11412973_1_report.txt
  Total species: 1118
  Top species: Bacteroides (17.29%)

SRR_reads/SRR11412973_2_report.txt
  Total species: 1169
  Top species: Bacteroides (17.32%)

SRR_reads/SRR11412976_1_report.txt
  Total species: 745
  Top species: Bacteroides (29.57%)

SRR_reads/SRR11412976_2_report.txt
  Total species: 885
  Top species: Bacteroides (29.57%)

SRR_reads/SRR11412979_1_report.txt
  Total species: 1223
  Top species: Segatella (74.85%)

SRR_reads/SRR11412979_2_report.txt
  Total species: 1329
  Top species: Segatella (74.78%)

SRR_reads/SRR11412980_1_report.txt
  Total species: 1478
  Top species: Bacteroides (21.18%)

SRR_reads/SRR11412980_2_report.txt
  Total species: 1544
  Top species: Bacteroides (21.26%)

SRR_reads/SRR11412984_1_report.txt
  Total species: 1353
  Top species: Segatella (43.74%)

SRR_reads/SRR11412984_2_report.txt
  Total species: 1410
  Top species: Segatella (43.76%)

SRR_reads/SRR21907296_1_report.txt
  Total species: 26
  Top species: Severe (99.87%)

SRR_reads/SRR21907296_2_report.txt
  Total species: 62
  Top species: Severe (98.7%)

SRR_reads/SRR21907296_report.txt
  Total species: 1
  Top species: Severe (100.0%)

SRR_reads/SRR21907303_1_report.txt
  Total species: 255
  Top species: Severe (98.59%)

SRR_reads/SRR21907303_2_report.txt
  Total species: 99
  Top species: Severe (97.75%)

SRR_reads/SRR21907303_report.txt
  Total species: 0
  Top species:  (0.0%)

SRR_reads/SRR21907307_1_report.txt
  Total species: 27
  Top species: Severe (99.63%)

SRR_reads/SRR21907307_2_report.txt
  Total species: 44
  Top species: Severe (94.12%)

SRR_reads/SRR21907330_1_report.txt
  Total species: 10
  Top species: Severe (99.97%)

SRR_reads/SRR21907330_2_report.txt
  Total species: 14
  Top species: Xanthomonas (72.06%)

SRR_reads/SRR21907332_1_report.txt
  Total species: 18
  Top species: Severe (99.99%)

SRR_reads/SRR21907332_2_report.txt
  Total species: 19
  Top species: Severe (99.98%)
```

### Brief Analysis

- The number of species classified for gut samples (between 745 and 1455) is generally higher than the number of species classified for wastewater samples (between 10 and 255).
- The top species classified in the gut samples are `Bacteroides` and `Segatella` while the top species classified in the wastewater samples are `Severe` and `Xanthomonas`.
  - The former are related to gut microbiomes while the latter are related to plant pathogens that can occur in wastewater samples.
- The two differences mentioned above can be good classifiers for gut and wastewater samples. Further analysis can be generally done, but I think the two differences mentioned above are already good indicators of the sample type in this case.