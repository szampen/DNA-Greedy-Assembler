import random
import logging
import sys
from typing import List, Tuple
from Bio import SeqIO
from Bio import Seq

# Configuration of basic logging for better console visualization
logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stdout)

def generate_random_reads(sequence : str, fragments_len : int = 200, expected_coverage : int = 5, output_file : str = "reads.fasta") -> None:
    """
    Simulates random DNA sequencing by generating reads from a reference sequence.
    """
    number_of_reads = int((len(sequence) / fragments_len) * expected_coverage)
    reads = []
    max_start = len(sequence) - fragments_len

    if max_start <= 0:
        logging.warning(f"Sequence too short for specified fragment length")
        return

    for i in range(number_of_reads):
        start_pos = random.randint(0, max_start)
        read = sequence[start_pos:start_pos + fragments_len]
        reads.append(SeqIO.SeqRecord(Seq.Seq(read),f"read_{i + 1}",description=""))

    SeqIO.write(reads,output_file,"fasta")
    logging.info(f"Successfully generated {number_of_reads} reads to {output_file}")

def find_overlap(seq1: str, seq2: str, min_overlap: int = 1) -> int:
    """
    Finds the maximum suffix-prefix overlap between two sequences.
    """
    max_possible_overlap = min(len(seq1), len(seq2))

    for overlap_len in range(max_possible_overlap, min_overlap, -1):
        if seq1[-overlap_len:] == seq2[:overlap_len]:
            return overlap_len

    return 0

def load_reads_fasta(fasta_filename : str) -> List[str]:
    """
    Loads sequence from a FASTA file into a list of strings.
    """
    return [str(record.seq) for record in SeqIO.parse(fasta_filename,"fasta")]


def merge_reads_into_sequence(reads: List[str]) -> List[str]:
    """
    Assembles reads into contigs using an iterative greedy approach.
    """
    iteration = 0
    while True:
        best_overlap = 0
        best_pairs: List[Tuple[int, int]] = []

        # Searching for the best match between pairs
        for i in range(len(reads)):
            for j in range(len(reads)):
                if i == j:
                    continue

                overlap = find_overlap(reads[i], reads[j])
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_pairs = [(i,j)]
                elif overlap == best_overlap and overlap > 0:
                    best_pairs.append((i,j))

        # If a match is not found, then finish assembling
        if best_overlap == 0:
            logging.info("Assembly finished. No more overlaps found.")
            break

        # Random pick from the pair with the best score
        best_i, best_j = random.choice(best_pairs)
        merged_seq = reads[best_i] + reads[best_j][best_overlap:]

        # Creating new list without merged fragments, adding new contig
        new_reads = [read for k, read in enumerate(reads) if k not in (best_i,best_j)]
        new_reads.append(merged_seq)
        reads = new_reads

        iteration += 1
        if iteration % 10 == 0:
            logging.info(f"Iteration {iteration}: {len(reads)} contig remaining, best overlap {best_overlap}.")

    return sorted(reads, key=len, reverse=True)


def run_test(sequence_file: str, sample_size: int = 10000) -> None:
    """
    Main test pipeline: Loads sequence, generate reads and assembles them.
    """
    try:
        records = list(SeqIO.parse(sequence_file, "fasta"))
        if not records:
            logging.error("Source file is empty or invalid.")
            return

        original_sequence = str(records[0].seq)[:sample_size]
        logging.info(f"Testing on sequence length: {len(original_sequence)}")

        generate_random_reads(original_sequence)

        reads = load_reads_fasta("reads.fasta")
        contigs = merge_reads_into_sequence(reads)

        print("\n" + "="*30)
        print("ASSEMBLY RESULTS")
        print("="*30)
        for i,contig in enumerate(contigs):
            print(f"Contig {i+1}: Length {len(contig)}")
            if i == 0: # Printing only the longest contig for a clear result
                print(f"Top Contig: {contig[:100]}")

        total_len = sum(len(c) for c in contigs)
        print(f"\nTotal length of all contigs: {total_len}")

    except FileNotFoundError:
        logging.error(f"File {sequence_file} not found. Please provide a valid FASTA file.")


if __name__ == "__main__":
    # Test is computed by default on >gi|224514821|ref|NT_167199.1| Homo sapiens chromosome Y genomic contig, GRCh37.p5 Primary Assembly
    input_file = sys.argv[1] if len(sys.argv) > 1 else "hs_ref_GRCh37.p5_chrY.fa"

    run_test(input_file)