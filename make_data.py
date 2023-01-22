import numpy as np
import pandas as pd

SOX_FILE = "data/raw/sox.tsv"
POUF_FILE = "data/raw/pouf.tsv"
NANOG_FILE = "data/raw/nanog.tsv"
ATAC_FILE = "data/raw/all_ATAC_peak_category_and_location.csv"

TSV_COLUMNS = ["peak", "length", "covered", "sum", "mean0", "mean", "min", "max"]


atac_df = pd.read_csv(ATAC_FILE)
nanog_df = pd.read_table(NANOG_FILE, sep='\t', names=TSV_COLUMNS)
pouf_df = pd.read_table(POUF_FILE, sep='\t', names=TSV_COLUMNS)
sox_df = pd.read_table(SOX_FILE, sep='\t', names=TSV_COLUMNS)

for df in (nanog_df, pouf_df, sox_df):
    dtmp = pd.merge(df[['peak']], atac_df.drop_duplicates('peak', keep='first'), on='peak', how='left')
    df['start'] = dtmp['start'].values
    df['end'] = dtmp['end'].values
    df['chr'] = dtmp['chr'].values

for df, name in zip((nanog_df, pouf_df, sox_df),("Nanog", "Pou5f3", "Sox19b")):
    cutoff = np.median(df["mean"])
    above_cutoff_df = df[df["mean"] >= cutoff]
    with open(f"data/bed/{name}.bed", "w") as of:
        for _, r in above_cutoff_df.iterrows():
            of.write(f"chr{r['chr']}\t{r['start']}\t{r['end']}\t{r['peak']}\t{r['mean']}\t+\n")
