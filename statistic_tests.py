from scipy.stats import chi2_contingency
import pandas as pd
import numpy as np
# Create the contingency table
if __name__ == '__main__':
    df = pd.read_csv("processed_data.csv")
    print(df.columns)
    df1 = df[df['has_motif_CCATTARCAT'] == True]
    df2 = df[df['has_motif_CCAATCARNG'] == True]
    df3 = df[df['has_motif_TTAATAGCCC'] == True]

    motifA = df1['count'].value_counts(sort=True, ascending=True).to_numpy()
    motifB = df2['count'].value_counts(sort=True, ascending=True).to_numpy()
    motifC = df3['count'].value_counts(sort=True, ascending=True).to_numpy()
    # Perform the chi-squared test


    chi2, p, dof, expected = chi2_contingency(np.array([motifA,motifB,motifC]))
    table = pd.DataFrame({"CCATTARCAT":motifA, "CCAATCARNG" :motifB,"TTAATAGCCC":motifC}, index = [0,1,2])
    print(table)
    print("Chi-squared test statistic:", chi2)
    print("p-value:", p)
