import random
from Bio import SeqIO
from Bio import Seq

def generate_random_reads(sequence, fragments_len = 200, expected_cover_len = 5):
    number_of_reads = int((len(sequence) / fragments_len) * expected_cover_len)
    reads = []
    for i in range(number_of_reads):
        max_start = len(sequence) - fragments_len
        if max_start > 0:
            start_pos = random.randint(0, max_start)
            read = sequence[start_pos:start_pos + fragments_len]
            reads.append(SeqIO.SeqRecord(Seq.Seq(read),f"read_{i + 1}",description=""))
    SeqIO.write(reads,"reads.fasta","fasta")

def find_overlap(seq1, seq2, min_overlap = 1):
    max_overlap = min(len(seq1), len(seq2))

    for overlap_len in range(max_overlap, min_overlap, -1):
        if seq1[-overlap_len:] == seq2[:overlap_len]:
            return overlap_len

    return 0

def load_reads_fasta(fasta_filename):
    reads = []
    for record in SeqIO.parse(fasta_filename,"fasta"):
        reads.append(str(record.seq))
    return reads

def merge_reads_into_sequence(reads,iteration = 0):
    best_overlap = 0
    best_pairs = []
    for i in range(len(reads)):
        for j in range(len(reads)):
            if i != j:
                overlap = find_overlap(reads[i], reads[j])
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_pairs = [(i,j)]
                elif overlap == best_overlap and overlap > 0:
                    best_pairs.append((i,j))

    if best_overlap == 0:
        return sorted(reads, key=len, reverse=True)

    best_i,best_j = random.choice(best_pairs)

    merged_seq = reads[best_i] + reads[best_j][best_overlap:]

    new_reads = []
    for k, read in enumerate(reads):
        if k != best_i and k != best_j:
            new_reads.append(read)
    new_reads.append(merged_seq)

    if iteration % 10 == 0:
        print(f"   Iteracja {iteration}: {len(new_reads)} kontigów, najlepszy overlap: {best_overlap}")

    return merge_reads_into_sequence(new_reads, iteration + 1)

def test(sequence_file):
    records = list(SeqIO.parse(sequence_file, "fasta"))
    original_sequence = str(records[0].seq)
    fragment = original_sequence[:10000]
    print(len(fragment))
    print(fragment)

    generate_random_reads(fragment)

    contigs = merge_reads_into_sequence(load_reads_fasta("reads.fasta"))

    sum = 0
    for i, contig in enumerate(contigs):
        sum += len(contig)
        print(f"\nKontig {i+1}:")
        print(f"  Długość: {len(contig)} nukleotydów")
        print(contig)
    print(f"\nSum: {sum}")

#test wykonywany jest na >gi|224514821|ref|NT_167199.1| Homo sapiens chromosome Y genomic contig, GRCh37.p5 Primary Assembly
test("hs_ref_GRCh37.p5_chrY.fa")