import numpy as np
import pandas as pd
import nanog_motif_hmm



if __name__ == '__main__':

    motif1 = "CCATTARCAT"  # advanced known data, should be  known when running hmm
    data = {"motifs":None,"start":'',"stop":""}
    result_sequences = nanog_motif_hmm.list_prediction #get results from hmm file
    index = result_sequences  # seq to enter / regions
    df = pd.DataFrame(data,index = result_sequences)
    motif_num = 1



    for seq in result_sequences:
        i = 0
        while i < len(seq):
            if (seq[i] == 'I'):
                for j in range(6):
                    if (seq[j] != 'I'):
                        break
                df[seq][motif] = True
                df[seq]["start"] = i
                df[seq]["end"] = i+5
                i = i+6
            else:
                i = i+1

    print(df)
