# Introduction
- [Introduction](#introduction)
- [Task 1: Metagenome Classification by String Matching](#task-1-metagenome-classification-by-string-matching)
  - [Task 1.1 (Multiple Matches)](#task-11-multiple-matches)
  - [Task 1.2 (Exact Matching)](#task-12-exact-matching)
  - [Task 1.4 (Comparison with Bioinformatics tools)](#task-14-comparison-with-bioinformatics-tools)
- [Task 2: Metagenomic classification by k-mer index](#task-2-metagenomic-classification-by-k-mer-index)
  - [Task 2.1 (Build the k-mer Index)](#task-21-build-the-k-mer-index)
  - [Task 2.2 (Implement Classification)](#task-22-implement-classification)
  - [Task 2.3 (Minimzers)](#task-23-minimzers)
- [Task 3: Real-world data and tools](#task-3-real-world-data-and-tools)
  - [Task 3.1 (Comparison)](#task-31-comparison)
  - [Task 3.2 (Real-world use case)](#task-32-real-world-use-case)

# Task 1: Metagenome Classification by String Matching

## Task 1.1 (Multiple Matches)

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

Run `/usr/bin/time -l go run task_1_2.go all` to get the results for all files.

<details>
<summary>Output for this task:</summary>
<br>

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

</details>


<details>
<summary>To reproduce the results for this task:</summary>
<br>

Install the genome and reads data to the `data/` directory.

Save sequence reads to `data/sequence_reads/` directory. Download genome data from NCBI and rename each genome ncbi dataset directory to its corresponding name of the following:
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

Run the following command (in `src/` directory) to get the results for all files:

```bash
/usr/bin/time -l go run task_1_2.go all
```


</details>

## Task 1.4 (Comparison with Bioinformatics tools)

- Execution time:
  - Not much difference in execution time between my implementation in Task 1.2 and `blastn` tool. They both took around 5 seconds to run. However, my implementation only considers exact matches while BLAST considers approximate matches.
- Memory usage:
  - My implementation's maximum resident set size was around 429MB while BLAST's was around 184MB.
- The resulting counts are quite different.
  - The BLAST tool uses a more sophisticated algorithm to find approximate matches and returns a larger number of matches.
  - BLAST considers statistical significance of matches and requires parameters of `-max_target_seqs` and `-evalue` to control the number of matches returned.


Run `sh run_blast_analysis.sh` to get the results for this task.

`run_blast_analysis.sh` assumes that the genome files are in the `data/` directory and to be run in the `src/` directory. Also, it uses `seqtk` to convert the FASTQ files to FASTA format, so a total of 3 commands are needed to be installed: `makeblastdb`, `blastn`, and, `seqtk`.

<details>
<summary>Output for this task:</summary>
<br>

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

</details>

# Task 2: Metagenomic classification by k-mer index

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

Run `/usr/bin/time -l go run task_2_1.go` to get the results for this task.

<details>
<summary>Output for this task:</summary>
<br>

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

</details>

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

Run `/usr/bin/time -l go run task_2_2.go` to get the results for this task.

<details>
<summary>Output for this task:</summary>
<br>

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

</details>

## Task 2.3 (Minimzers)

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

Run `/usr/bin/time -l go run task_2_3.go` to get the results for this task.

<details>
<summary>Output for this task:</summary>
<br>

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

</details>


# Task 3: Real-world data and tools

## Task 3.1 (Comparison)

## Task 3.2 (Real-world use case)