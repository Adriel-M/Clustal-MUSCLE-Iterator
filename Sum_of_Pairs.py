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
    first_in_gap = False
    second_in_gap = False
    for position in range(len(sequence_1)):
        pair = (sequence_1[position], sequence_2[position])
        # Both in gap
        if first_in_gap and second_in_gap:
            # Both gaps are extended
            if pair == ("-", "-"):
                score -= gap_extension
                score -= gap_extension
            # Only the first has a gap
            elif pair[0] == "-":
                score -= gap_extension
                second_in_gap = False
            # Only the second has a gap
            elif pair[1] == "-":
                score -= gap_extension
                first_in_gap = False
            # No more gaps
            else:
                first_in_gap = False
                second_in_gap = False
                score += pair_score(pair, matrix)
        elif first_in_gap:
            # First gap is extended, open new gap in second
            if pair == ("-", "-"):
                score -= gap_open
                score -= gap_extension
                second_in_gap = True
            # Only the first has a gap, first is extended
            elif pair[0] == "-":
                score -= gap_extension
            # Only second has a gap, first gap closes
            elif pair[1] == "-":
                score -= gap_open
                second_in_gap = True
                first_in_gap = False
            # No more gaps
            else:
                first_in_gap = False
                score += pair_score(pair, matrix)
        elif second_in_gap:
            # Second gap is extended, open new gap in first
            if pair == ("-", "-"):
                score -= gap_open
                score -= gap_extension
                first_in_gap = True
            # Only the first has a gap, second gap closes
            elif pair[0] == "-":
                score -= gap_open
                first_in_gap = True
                second_in_gap = False
            # Only second has a gap, second is extended
            elif pair[1] == "-":
                score -= gap_extension
            # No more gaps
            else:
                second_in_gap = False
                score += pair_score(pair, matrix)
        else:
            if '-' in pair:
                if pair[0] == "-":
                    first_in_gap = True
                    score -= gap_open
                if pair[1] == "-":
                    second_in_gap = True
                    score -= gap_extension
            else:
                score += pair_score(pair, matrix)
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
