import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import scikit_posthocs as sp

BLUE = "#3399FF"
RED = "#FF3399"
PURPLE = "#CC99FF"
if __name__ == '__main__':

    # Prepare the data
    data = {
           'CCATTARCAT': [True, True, False, True, True, True],
            'CCAATCARNG': [False, True, False, False, True, True],
            'TTAATAGCCC': [True, False, False, False, False, True],
            'sox': [False, False, False, True, False, True],
            'pouf': [True, False, True, False, False, True],
            'value': [15, 30, 8, 12, 5, 11],
            'count': [1,2,3,2,2,1]}
    # df = pd.DataFrame(data)
    df = pd.read_csv("processed_data.csv")
    print(df.columns)
    df.rename(columns = {'has_motif_CCATTARCAT': 'CCATTARCAT',
                         'has_motif_CCAATCARNG': 'CCAATCARNG',
                         'has_motif_TTAATAGCCC' : 'TTAATAGCCC'}, inplace=True)
    # create a dataframe with the regions that have a motif
    df['motif'] = df[['CCATTARCAT', 'CCAATCARNG', 'TTAATAGCCC']].apply(
        lambda x: ''.join(x.index[x]), axis=1)

    ###############bar plot################################
    m1 = df.loc[df["CCATTARCAT"], "value"]
    m2 = df.loc[df["CCAATCARNG"], "value"]
    m3 = df.loc[df["TTAATAGCCC"], "value"]
    _, pv = kruskal(m1,m2,m3)
    post_hoc_pvalues = sp.posthoc_dunn([m1,m2,m3]) > 0.5
    post_hoc_pvalues.to_html("posthoc.html")
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(16,12), sharey=True, sharex=True)
    fig.suptitle(f"ChIP-Seq Strength Differs Between Motifs\nKruskal-Wallis test p-value: {pv:.2E}", size=20)
    fig.supxlabel("ChIP-Seq strength", fontsize=14)
    fig.supylabel("# peaks density", fontsize=14)
    ax0.set_title("Motif CCATTARCAT", fontdict={"size": 16})
    ax1.set_title("Motif CCAATCARNG", fontdict={"size": 16})
    ax2.set_title("Motif TTAATAGCCC", fontdict={"size": 16})
    ax0.hist(m1, bins=100, density=True, range=(0,15))
    ax1.hist(m2, bins=100, density=True, range=(0,15))
    ax2.hist(m3, bins=100, density=True, range=(0,15))
    fig.savefig("strength_distribution_by_motif.png")

    title = "Barplot of pick region with motif and without motif"
    motif_regions = df[df['motif'] != ""]
    count_motif_regions = len(motif_regions)

    # create a dataframe with the regions that do not have a motif
    no_motif_regions = df[df['motif'] == ""]
    count_no_motif_regions = len(no_motif_regions)

    fig = plt.subplots(figsize=(10, 7))
    # create a bar plot
    bar1 = plt.bar('Motif Regions',
            count_motif_regions,
            color=BLUE)
    bar2 = plt.bar('No Motif Regions',
             count_no_motif_regions,
            color= RED)
    plt.legend(['Motif Regions', 'No Motif Regions'])
    plt.xlabel("Region Type")
    plt.ylabel("Count")
    plt.title(title)
    #plt.show()
    plt.savefig("bar_plot_has_motif.png")


    ###############box plot################################


    title = "Box plot"
    sns.boxplot(x='motif', y='value',
                data=df, palette='pastel')
    sns.stripplot(x='motif', y='value', hue='count',
                  data=df, jitter=True, dodge=True, palette='pastel', size=5,
                  marker='o')
    plt.xlabel("Motifs")
    plt.title(title)
    #plt.show()
    plt.savefig("box_plot.png")


    ###############bar plot2################################

    # Create a bar for the total number of true values
    motifs = ["CCATTARCAT", "CCAATCARNG", "TTAATAGCCC"]
    proteins = ["sox", "pouf"]

    motif_total_count = [] # number of picks for each motif [2,0]
    motif_vs_motif_total_count = [] # number of picks for each motif [2,2]
    motif_vs_all_motifs_total_count = (df['CCATTARCAT'].astype(int) + df['CCAATCARNG'].astype(int) + df['TTAATAGCCC'].astype(int)).count()
    motif_both_proteins_counts = []

    for motif in motifs:
        motif_total_count.append(df.query(f'{motif} == True').shape[0])
        motif_both_proteins_counts.append(df.query(f'{motif} == True & sox == True & pouf == True').shape[0])
        motif_vs_motif = []
        for motif in motifs:
            motif_vs_motif.append(df.query(f'{motif} == True & {motif} == True').shape[0] - motif_vs_all_motifs_total_count)
        motif_vs_motif_total_count.append(motif_vs_motif)

    prob_count = []
    for protein in proteins:
        motif_counts = []
        for i,motif in enumerate(motifs):
            motif_counts.append(df.query(f'{motif} == True & {protein} == True').shape[0] - motif_both_proteins_counts[i])
            prob_count.append(motif_counts)
    ind = np.arange(len(motifs))
    width = 0.35

    fig = plt.subplots(figsize=(10, 7))
    p_sox = plt.bar(ind, prob_count[0], width, color=RED)
    p_pouf = plt.bar(ind, prob_count[1], width, bottom=prob_count[0], color=BLUE)
    p_both = plt.bar(ind, motif_both_proteins_counts, width, bottom=[i+j for i, j in zip(prob_count[0], prob_count[1])], color=PURPLE)
    plt.ylabel('Number of picks')
    plt.xlabel('Motif seq')
    plt.title('number of picks with sox and pouf  for each motif of nanog')
    plt.xticks(ind, ('CCATTARCAT', 'CCAATCARNG', 'TTAATAGCCC'))
    plt.legend((p_sox,p_pouf,p_both), ('sox', 'pouf', 'both'))

    #plt.show()
    plt.savefig("bar_plot_count.png")




    fig = plt.subplots(figsize=(10, 7))
    CCATTARCAT = plt.bar(ind, motif_vs_motif[0], width, color=RED)
    CCAATCARNG = plt.bar(ind, motif_vs_motif[1], width, bottom=motif_vs_motif[0], color=BLUE)
    # TTAATAGCCC = plt.bar(ind,  motif_vs_motif[2], width, bottom=[i+j for i, j in zip(motif_vs_motif[0], motif_vs_motif[1])], color=PURPLE)
    plt.ylabel('number of picks')
    plt.title('number of picks for each motif')
    plt.xticks(ind, ('CCATTARCAT', 'CCAATCARNG', 'TTAATAGCCC'))
    plt.legend((p_sox,p_pouf,p_both), ('sox', 'pouf', 'both'))

    #plt.show()
    plt.savefig("bar_plot_count.png")


    #
    #
    # motif_counts = df[df['CCATTARCAT'] == True].count()
    # motif_counts = df[df == True].count()
    #
    # motif_protein_counts = df.groupby(["CCATTARCAT", "CCAATCARNG", "TTAATAGCCC"])[
    #     "number_of_pro"].value_counts().unstack(fill_value=0)
    # colors = cm.rainbow(np.linspace(0, 1, 4))
    # for motif in ["CCATTARCAT", "CCAATCARNG", "TTAATAGCCC"]:
    #     plt.bar(motif, motif_counts[motif], color='gray')
    #     # Create bars for the number of true values for each protein count
    #     for i in range(4):
    #         plt.bar(motif, motif_protein_counts[motif][i],
    #                 bottom=sum(motif_protein_counts[motif][:i]),
    #                 color=colors[i])
    # # Show the plot
    # #plt.show()



    #
    # # filter the dataframe to only include rows where the CCATTARCAT is equal to True
    # df_CCATTARCAT = df[df["CCATTARCAT"] == True]
    #
    # # filter the dataframe to only include rows where the CCAATCARNG is equal to True
    # df_CCAATCARNG = df[df["CCAATCARNG"] == True]
    #
    # # filter the dataframe to only include rows where the TTAATAGCCC is equal to True
    # df_TTAATAGCCC = df[df["TTAATAGCCC"] == True]
    #
    # # concatenate the three dataframe
    # df_motifs = pd.concat([df_CCATTARCAT, df_CCAATCARNG, df_TTAATAGCCC])
    #
    # # create the boxplot
    # bp = sns.boxplot(x="motif", y="value", hue="number_of_pro",
    #                  data=df_motifs, showfliers=True, fliersize=2)
    #
    # # loop through the box, whiskers and caps of the boxplot and adjust the color and transparency
    # for patch in bp.artists:
    #     r, g, b, a = patch.get_facecolor()
    #     patch.set_facecolor((r, g, b, .3))
    #
    # # add title and axis labels
    # plt.title("Motifs")
    # plt.xlabel("Motifs")
    # plt.ylabel("Value")
    #
    # # show the plot
    # #plt.show()
    #
    # # fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))
    # # for i,df in enumerate([df_CCATTARCAT, df_CCAATCARNG, df_TTAATAGCCC]):
    # #     sns.boxplot(y='value',
    # #                 data=df, palette='pastel', ax=ax[i])
    # #     sns.stripplot(y='value', hue='number_of_pro',
    # #                   data=df, jitter=True, dodge=True, palette='pastel', size=5,
    # #                   marker='o', ax=ax[i])
    # # # plt.xlabel(["CCATTARCAT", "CCAATCARNG", "TTAATAGCCC"])
    # # plt.ylabel("Value")
    # # #plt.show()
    #
    #
    #
    #
    # # title = "Box plot"
    # #
    # # motifs = ['CCATTARCAT', 'CCAATCARNG', 'TTAATAGCCC']
    # # for motif in motifs:
    # #     sns.boxplot(x=motif, y='value', hue='number_of_pro',
    # #                 data=df, showfliers=True)
    # #     sns.stripplot(x=motif, y='value',
    # #                   hue='number_of_pro', data=df,
    # #                   jitter=True, alpha=0.5)
    # # # sns.boxplot(x='motif', y='value',
    # # #             data=df, palette='pastel')
    # # # sns.stripplot(x='motif', y='value', hue='number_of_pro',
    # # #               data=df, jitter=True, dodge=True, palette='pastel', size=5,
    # # #               marker='o')
    # # plt.title(title)
    # # #plt.show()
    # # # filter the dataframe to only include rows where the CCATTARCAT is equal to True
    # # df_CCATTARCAT = df[df["CCATTARCAT"] == True]
    # #
    # # # create the boxplot for CCATTARCAT
    # # bp = sns.boxplot(x="CCATTARCAT", y="value",
    # #                  hue="number_of_pro", data=df_CCATTARCAT,
    # #                  showfliers=True, fliersize=2)
    # #
    # # # loop through the box, whiskers and caps of the boxplot and adjust the color and transparency
    # # for patch in bp.artists:
    # #     r, g, b, a = patch.get_facecolor()
    # #     patch.set_facecolor((r, g, b, .3))
    # #
    # # # add title and axis labels
    # # plt.title("CCATTARCAT")
    # # plt.xlabel("CCATTARCAT")
    # # plt.ylabel("Value")
    # #
    # # # show the plot
    # # #plt.show()
    # #
    #
    #
    # # df1 = pd.readcsv('data/bed/Nanog.bed', sep='/t')
    # # df2 = pd.readcsv('data/bed/Nanog.bed', sep='/t')
    # #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #

