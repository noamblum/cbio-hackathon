import csv
import numpy as np
import pandas as pd

DISTANCE = 40

NanogData = np.zeros((27986, 4), dtype=np.float64)
PouData = np.zeros((27985, 4), dtype=np.float64)
SoxData = np.zeros((27985, 4), dtype=np.float64)


def ReadFile(path):
    df = pd.read_table(path, sep='\t', names=["chrom", "start", "end", "peak", "mean", "strand"])
    return df[["chrom", "start", "end", "mean"]]
    # with open(path) as f:
    #     file = csv.reader(f, delimiter="\t")
    #     for i, line in enumerate(file):
    #         array[i][0] = float(line[0][3:])
    #         array[i][1] = line[1]
    #         array[i][2] = line[2]
    #         array[i][3] = line[4]


def Match(array_main, array_ref, PorS):
    return pd.merge(array_main, array_ref, on=["chrom", "start", "end"], how="inner")
    # x = np.zeros((27986, 6), dtype=np.float64)
    # x[:, 0] = NanogData[:, 0]
    # x[:, 1] = NanogData[:, 1]
    # x[:, 2] = NanogData[:, 2]
    # result = pd.DataFrame(x, columns=['name', 'startRange', 'endRange', 'pou', 'sox', 'count'])
    # for i in range(27985):
    #     for j in range(27985):
    #         if array_main[i][0] == array_ref[j][0]:
    #             # ---------------         Nanog
    #             #        ---------------- PorS
    #             if((array_main[i][1]<array_ref[j][1]) & (array_main[i][2]<array_ref[j][2]) & (array_main[i][2]>array_ref[j][1])):
    #                 result[PorS][j] = 1
    #             #         --------------- Nanog
    #             #---------------          PorS
    #             elif((array_main[i][1]>array_ref[j][1]) & (array_main[i][2]>array_ref[j][2]) & (array_main[i][1]<array_ref[j][2])):
    #                 result[PorS][j] = 1
    #             # ------------------------- Nanog
    #             #      --------------       PorS
    #             elif((array_main[i][1]<array_ref[j][1]) & (array_main[i][2]>array_ref[j][2])):
    #                 result[PorS][j] = 1
    #             #       --------------      Nanog
    #             # ------------------------- PorS
    #             elif ((array_main[i][1]>array_ref[j][1]) & (array_main[i][2]<array_ref[j][2])):
    #                 result[PorS][j] = 1
    #     print(i)
    # result.to_csv("data.csv")


def main():
    nanog_data = ReadFile("data/bed/Nanog.bed")
    pouf_data = ReadFile("data/bed/Pou5f3.bed")
    sox_data = ReadFile("data/bed/Sox19b.bed")
    Match(nanog_data, pouf_data, 'pou').to_csv("Nanog_Pou5f3_overlap.csv", index=False)
    Match(nanog_data, sox_data, 'sox').to_csv("Nanog_Sox19b_overlap.csv", index=False)






# df1 = np.zeros((27986, 5), dtype=np.float64)
# with open("data\\bed\\Nanog.bed")as f:
#     file = csv.reader(f, delimiter="\t")
#     for i, line in enumerate(file):
#         df1[i][0] = float(line[0][3:])
#         df1[i][1] = line[1]
#         df1[i][2] = line[2]
#         df1[i][3] = line[4]
#
# df2 = np.zeros((27985, 4), dtype=np.float64)
# with open("data\\bed\\Pou5f3.bed")as f:
#     file = csv.reader(f, delimiter="\t")
#     for i, line in enumerate(file):
#         df2[i, 0] = float(line[0][3:])
#         df2[i, 1] = line[1]
#         df2[i, 2] = line[2]
#         df2[i, 3] = line[4]
#
# df3 = np.zeros((27985, 4), dtype=np.float64)
# with open("data\\bed\\Sox19b.bed")as f:
#     file = csv.reader(f, delimiter="\t")
#     for i, line in enumerate(file):
#         df2[i, 0] = float(line[0][3:])
#         df2[i, 1] = line[1]
#         df2[i, 2] = line[2]
#         df2[i, 3] = line[4]

    #
    # for i in range(27985):
    #     for j in range(27985):
    #         if array1[i][0] == array2[j][0]:
    #             if ((array1[i][1] - DISTANCE) >= array2[j][1]) & ((array1[i][1] + DISTANCE) > array2[j][2]) & ((array1[i][1] - DISTANCE) < array2[j][2]):
    #                 result[PorS][j] = 1
    #                 print("TT")
    #             elif ((array1[i][1] - DISTANCE) < array2[j][1]) & ((array1[i][1] + DISTANCE) < array2[j][2]) & ((array1[i][1] + DISTANCE) < array2[j][2]):
    #                 result[PorS][j] = 1
    #
    #             elif ((array1[i][1] - DISTANCE) < array2[j][1]) & ((array1[i][1] + DISTANCE) > array2[j][2]) & ((array1[i][1] - DISTANCE) < array2[j][2]) & (
    #                     (array1[i][1] + DISTANCE) < array2[j][2]):
    #                 result[PorS][j] = 1
    #                 print("TT")



# df1 = pd.read_table("data\\bed\\Nanog.bed", header=None, names=["name", "start", "end", "id", "value", "strand"])
# df2 = pd.read_table("data\\bed\\Pou5f3.bed", header=None, names=["name", "start", "end", "id", "value", "strand"])


# df1['motif'] = ''
# for i in range(1786):
#     for j in range(1949):
#         if df1[i][0] == df2[j][0]:
#             if(((df1[i][1] - 40) > df2[j][1]) & ((df1[i][1] + 40) > df2[j][2]) & ((df1[i][1] - 40) < df2[j][2])):
#                 df1[i][6] = "yes"
#             elif ((df1[i][1] - 40) < df2[j][1]) & ((df1[i][1] + 40) < df2[j][2]) & ((df1[i][1] + 40) < df2[j][2]):
#                 df1[i][6] = "yes"
#             elif ((df1[i][1] - 40) < df2[j][1]) & ((df1[i][1] + 40) > df2[j][2]) & ((df1[i][1] - 40) < df2[j][2]) & ((df1[i][1] + 40)< df2[j][2]):
#                 df1[i][6] = "yes"


# df1['motif'] = ''
# df1.loc[((df1.start - 40) > df2.start) & ((df1.start + 40) > df2.end) & ((df1.start - 40) < df2.end), 'motif'] = 'yes'
# df1.loc[((df1.start - 40) < df2.start) & ((df1.start + 40) < df2.end) & ((df1.start + 40) < df2.end), 'motif'] = 'yes'
# df1.loc[((df1.start - 40) < df2.start) & ((df1.start + 40) > df2.end) & ((df1.start - 40) < df2.end) & ((df1.start + 40)< df2.end), 'motif'] = 'yes'

if __name__ == '__main__':
    main()