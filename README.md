# Introduction

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

## Task 2.2 (Implement Classification)

## Task 2.3 (Minimzers)

# Task 3: Real-world data and tools

## Task 3.1 (Comparison)

## Task 3.2 (Real-world use case)