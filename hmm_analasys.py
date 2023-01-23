import numpy as np
import pandas as pd


data = {"motifs": ["1", "2", "3"],"start":'',"stop":""}
result_sequences = [] #get results from hmm file
index = result_sequences  # seq to enter / regions
df = pd.DataFrame(data, index=index)


motif = 'm'  # advanced known data, should be  known when running hmm
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
            
