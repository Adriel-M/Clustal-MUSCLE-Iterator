from Sum_of_Pairs import *
from Bio.SubsMat import MatrixInfo
import os
import subprocess


def iterative_clustal(path_to_megacc, path_to_unaligned, min_gap_open,
                      max_gap_open, gap_open_interval, min_gap_extension,
                      max_gap_extension, gap_extension_interval,
                      substitution_matrix):
    best_gap_config = []
    best_score = 0
    current_gap_open = min_gap_open
    while current_gap_open <= max_gap_open:
        formatted_gap_open = "%.2f" % current_gap_open
        current_gap_extension = min_gap_extension
        while current_gap_extension <= max_gap_extension:
            formatted_gap_extension = "%.2f" % current_gap_extension
            config_name = "clustal_open" + formatted_gap_open + "_extension" \
                          + formatted_gap_extension
            generate_config("configuration/" + config_name + ".mao",
                            formatted_gap_open, formatted_gap_extension)
            run_clustal(path_to_megacc, "configuration/" + config_name +
                        ".mao", path_to_unaligned, "results/" +
                        config_name + ".fasta")
            alignments = fasta_to_list("results/" + config_name + ".fasta")
            # default blastp is gap_open: 11, gap_extension: 1
            score = sum_of_pairs(alignments, substitution_matrix, 11.0, 1.0)
            current_config = [current_gap_open, current_gap_extension]
            if score == best_score:
                best_gap_config.append(current_config)
            elif score > best_score:
                best_gap_config.clear()
                best_score = score
                best_gap_config.append(current_config)
            current_gap_extension += gap_extension_interval
        current_gap_open += gap_open_interval
    return best_gap_config


def generate_config(config_path, gap_open, gap_extension):
    """Generate configuration of megacc

    :param config_path: output config path
    :param gap_open: the penalty for opening a gap
    :param gap_extension: the penalty for extending a gap
    """
    # check if file exists, if it does, don't generate a new config
    if not os.path.isfile(config_path):
        gap_open = float(gap_open)
        gap_extension = float(gap_extension)
        gap_open = "%.2f" % gap_open
        gap_extension = "%.2f" % gap_extension
        # Hardcoded for now.
        template = open("configuration/template/clustal_default.mao", "r")
        config = open(config_path, "w")
        while True:
            line = template.readline()
            if not line:
                break
            if "Penalty" in line:
                if "ProteinPWGapOpeningPenalty" in line:
                    config.write("ProteinPWGapOpeningPenalty   = " + gap_open)
                elif "ProteinMAGapOpeningPenalty" in line:
                    config.write("ProteinMAGapOpeningPenalty   = " + gap_open)
                elif "ProteinPWGapExtensionPenalty" in line:
                    config.write("ProteinPWGapExtensionPenalty = " +
                                 gap_extension)
                elif "ProteinMAGapExtensionPenalty" in line:
                    config.write("ProteinMAGapExtensionPenalty = " +
                                 gap_extension)
                config.write("\n")
            else:
                config.write(line)


def run_clustal(path_to_megacc, configuration_path, unaligned_sequences,
                output_path):
    """Runs clustal

    :param path_to_megacc: file path to megacc program
    :param configuration_path: file path to megacc alignment configuration
    :param unaligned_sequences: a collection of unaligned sequences
    :param output_path: file path to output the results
    """
    command_to_run = [path_to_megacc, "-a", configuration_path, "-d",
                      unaligned_sequences, "-f", "Fasta", "-o", output_path]
    subprocess.call(command_to_run)

if __name__ == "__main__":
    best = iterative_clustal("mega/megacc", "example/lol.fas", 1, 10, 0.5,
                             0.10, 1, 0.10, MatrixInfo.blosum62)
    print(best)
