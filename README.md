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

## Task 2.3 (Minimzers)

# Task 3: Real-world data and tools

## Task 3.1 (Comparison)

## Task 3.2 (Real-world use case)