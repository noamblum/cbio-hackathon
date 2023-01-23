import argparse
from itertools import groupby
from hmmlearn import hmm
import numpy as np
import gzip as gz

SIZE_MOTIF = 10
NUM_MOTIF = 1
NUM_HIDDEN_STATES = 1 + SIZE_MOTIF * NUM_MOTIF
NUM_EMISSIONS = 4
MAX_NUM_ITER = 5
NUM_TRAIN = 100
A = 0
C = 1
G = 2
T = 3


# HIDDEN_STATES = [B, M1 , M2 , M3, M4, M5, M6, M7, M8, M9, M10]
def initialize_p0_array():
    p = 0.99
    p0_array = [p, 1 - p, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    return p0_array


# HIDDEN_STATES = [B, M1 , M2 , M3, M4, M5, M6, M7, M8, M9, M10]
def initialize_emission_matrix():
    motif = np.loadtxt("data/motifs/Nanog/homerResults/motif1.motif",
                       comments=">")

    emission_matrix = np.zeros(shape=(NUM_HIDDEN_STATES, NUM_EMISSIONS))
    emission_matrix[0] = 1 / NUM_EMISSIONS
    emission_matrix[1:] = motif
    return emission_matrix


# HIDDEN_STATES = [B, M1 , M2 , M3, M4, M5, M6, M7, M8, M9, M10]
def initialize_transition_matrix():
    p = 0.5
    transition_matrix = np.array(
        [[p, 1 - p, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])


def initialize_nanog_motif_hmm():
    nanog_motif_hmm = hmm.CategoricalHMM(n_components=NUM_HIDDEN_STATES,
                                         n_iter=MAX_NUM_ITER,
                                         verbose=True,
                                         params='st', init_params='')
    nanog_motif_hmm.n_features = 4
    nanog_motif_hmm.emissionprob_ = initialize_emission_matrix()
    nanog_motif_hmm.startprob_ = initialize_p0_array()
    nanog_motif_hmm.transmat_ = initialize_transition_matrix()
    return nanog_motif_hmm


def fastaread(fasta_name):
    """
    Read a fasta file. For each sequence in the file, yield the header and the
     actual sequence. You may keep this function, edit it, or delete it and
     implement your own reader.
    """
    f = gz.open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line:
    line.decode().startswith(">")))
    for header in faiter:
        header = next(header)[1:].decode().strip()
        seq = "".join(s.decode().strip() for s in next(faiter))
        yield header, seq


def load_and_preprocess_fasta(filename):
    train_seqs = []
    train_len_seqs = []

    for header, orig_seq in fastaread(filename):
        seq = []
        for char in orig_seq:
            if char == "A" or "a":
                seq.append([A])
            elif char == "C" or "c":
                seq.append([C])
            elif char == "G" or "g":
                seq.append([G])
            elif char == "T" or "t":
                seq.append([T])
        train_seqs += seq
        train_len_seqs.append(len(seq))

    return np.array(train_seqs), train_len_seqs


def train_hmm(nanog_motif_hmm, len_seqs, seqs):
    for i in range(NUM_TRAIN):
        nanog_motif_hmm.fit(seqs, len_seqs)
        nanog_motif_hmm.transmat_[1:] = initialize_transition_matrix()[1:]
        nanog_motif_hmm.transmat_[0, 1] = np.sum(
            nanog_motif_hmm.transmat_[0, 1:])
        nanog_motif_hmm.transmat_[0, 2:] = 0

def load_fasta_to_seqs_array(filename):
    seqs_array = []
    for header, orig_seq in fastaread(filename):
        seq = []
        for char in orig_seq:
            if char == "A" or "a":
                seq.append([A])
            elif char == "C" or "c":
                seq.append([C])
            elif char == "G" or "g":
                seq.append([G])
            elif char == "T" or "t":
                seq.append([T])

        seqs_array.append(seq)
    return seqs_array

def write_decode_to_files(nanog_motif_hmm, likelihood_output_file,
                          prediction_output_file, filename):
    seqs_array = load_fasta_to_seqs_array(filename)
    f_likelihood_output_file = open(likelihood_output_file, 'w')
    f_prediction_output = open(prediction_output_file, 'w')

    for seq in seqs_array:
        f_likelihood_output_file.write(
            str(nanog_motif_hmm.decode(seq)[0]) + '\n')
        hidden_states_array = nanog_motif_hmm.predict(seq)
        for state in hidden_states_array:
            if 0 <= state <= 3:
                f_prediction_output.write('I')
            else:
                f_prediction_output.write('N')
        f_prediction_output.write('\n')


def main():
    nanog_motif_hmm = initialize_nanog_motif_hmm()

    train_seqs, train_len_seqs = \
        load_and_preprocess_fasta("data\fasta\Nanog.fa.txt")

    train_hmm(nanog_motif_hmm, train_len_seqs, train_seqs)


if __name__ == '__main__':
    main()
