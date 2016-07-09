from Bio.SubsMat import MatrixInfo
from Bio import SeqIO


def fasta_to_list(path_to_alignment):
    """Parse through a fasta formatted file and store its contents in a
    list. Uses SeqIO from Biopython.

    :param path_to_alignment: a string representation of the path to the
    fasta formatted file
    :return: a list with all the aligned sequences
    """
    sequence_collection = []
    for seq_record in SeqIO.parse(path_to_alignment, """fasta"""):
        sequence_collection.append(seq_record.seq)
    return sequence_collection


def pair_score(pair, matrix):
    """Return the score of a pair from a defined a substitution matrix.

    :param pair: a tuple representation of the pair to check
    :param matrix: the substitution matrix used to score pairs
    :return: the score pertaining to the pair
    """
    # since A-C and C-A will have the same score, sometimes the substitution
    # matrix will only contain one of these two.
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]


def alignment_score(sequence_1, sequence_2, matrix, gap_open, gap_extension):
    """Compare two sequence and generate a pairwise score.

    :precondition: both alignments must contain the same amount of elements
    :param sequence_1: one of the sequence to compare
    :param sequence_2: the other sequence to compare
    :param matrix: the substitution matrix used to score pairs
    :param gap_open: the penalty for extending a gap
    :param gap_extension: the sum of pair scores for this alignment
    :return: the pairwise score between these two sequences
    """
    score = 0.0
    in_gap = False
    for position in range(len(sequence_1)):
        pair = (sequence_1[position], sequence_2[position])
        if not in_gap:
            if '-' in pair:
                in_gap = True
                score -= gap_open
            else:
                score += pair_score(pair, matrix)
        else:
            if '-' not in pair:
                in_gap = False
                score += pair_score(pair, matrix)
            else:
                score -= gap_extension
    return score


def sum_of_pairs(alignments, matrix, gap_open, gap_extension):
    """Generate all possible pairs and get their alignment scores.

    :param alignments: the collection of all sequences
    :param matrix: the substitution matrix used to score pairs
    :param gap_open: the penalty for extending a gap
    :param gap_extension: the sum of pair scores for this alignment
    :return:
    """
    sequence_total = len(alignments)
    score = 0.0
    for initial_sequence in range(sequence_total):
        for next_sequence in range(initial_sequence + 1, sequence_total):
            score += alignment_score(alignments[initial_sequence],
                                     alignments[next_sequence],
                                     matrix, gap_open, gap_extension)
    return score


if __name__ == "__main__":
    sample = fasta_to_list("./sample.fasta")
    sample_score = sum_of_pairs(sample, MatrixInfo.blosum62, 10, 0.2)
    print(sample_score)
