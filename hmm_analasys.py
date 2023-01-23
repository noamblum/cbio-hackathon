import numpy as np
import pandas as pd
import nanog_motif_hmm

if __name__ == '__main__':

    motif1 = "CCATTARCAT"  # advanced known data, should be  known when running hmm
    motif2 = "CCAATCARNG"
    motif3 = "TTAATAGCCC"
    motifs = (motif1, motif2, motif3)
    data = {"motifs": None, "start": '', "stop": ""}
    result_sequences = []  # get results from hmm file
    motif_num = 1

    for i in range(1, 4):
        result_sequences = []  # get results from hmm file
        with open(f"nanog_motif_predictions{i}.txt") as file:
            lines = file.readlines()
            for line in lines:
                result_sequences.append(line.strip())

        has_motif = ['M' in seq for seq in result_sequences]
        end = [seq.rfind('M') for seq in result_sequences]
        start = [seq.rfind('B', 0, e) for seq, e in zip(result_sequences, end)]
        df = pd.DataFrame({"has_motif": has_motif, "start": start, "end": end}, index=result_sequences)

        print(df)
        df.to_csv(f"find_motif_in_seq-motif{i}-{motifs[i - 1]}.csv")
