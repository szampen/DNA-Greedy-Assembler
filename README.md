# DNA Greedy Assembler

A Python-based tool for assembling genomic sequences from random fragments using **Greedy Overlap-Layout-Consensus** approach.
Developed as a part of Bioinformatics course.

## Project Overview

The goal of this tool is to reconstruct a continuous DNA sequence (contig) from a set of short, overlapping "reads".
It simulates the process of genomic assembly by:
1. **Generating random reads** from a reference sequence with a specified coverage.
2. **Finding overlap** using a suffix-prefix matching algorithm.
3. **Merging sequences** iteratively based on the best available overlap

## Key Features

* **Iterative Assembly Algorithm:** Uses an optimized loop-based approach to merge fragments, ensuring stability for larger datasets.
* **Bioinformatics Integration:** Uses the `Biopython` library for handling standard FASTA file formats.
* **Customizable Parameters:** Easily adjust fragment length and expected coverage via command-line arguments.
* **Real-world Testing:** Validated on 10,000 nucleotide fragment of the Homo sapiens chromosome Y.

## Tech stack

* **Language:** Python 3.x
* **Libraries:** Biopython, Random
* **Data formats:** FASTA (.fa, .fasta)

## Performance & Complexity

The algorithm uses a greedy heuristic with a computational complexity of approximately $O(n^2 \cdot m)$, where $n$ is the number of reads and $m$ is the average read length.

> [!TIP]
> This project is ideal for understanding the fundamentals of assembly. For gigabase-scale genomes, technologies like De Bruijn Graphs or Suffix Trees are typically preferred for better memory efficiency.

## Installation & Usage

1. **Clone the repository:**
   ```bash
   git clone https://github.com/szampen/DNA-Greedy-Assembler.git
   
2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   
3. **Run the assembler:**
   ```bash
   python assembler.py [your_filename.fa/.fasta]
   
## Example output

```text
Testing on sequence length: 10000
Successfully generated 250 reads to reads.fasta
Iteration 10: 240 conting remaining, best overlap 198.
Iteration 20: 230 conting remaining, best overlap 196.
Iteration 30: 220 conting remaining, best overlap 194.
...
Assembly finished. No more overlaps found.

==============================
ASSEMBLY RESULTS
==============================
Contig 1: Length 4968
Top Contig: CTCCCCTCGGGACCACCCCAGACCCCC... [truncated]
Contig 2: Length 4504
Contig 3: Length 436

Total length of all contigs: 9908
