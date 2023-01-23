import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

BLUE = "#3399FF"
RED = "#FF3399"
PURPLE = "#CC99FF"
if __name__ == '__main__':

    # Prepare the data
    data = {
           'motifA': [True, True, False, True, True, True],
            'motifB': [False, True, False, False, True, True],
            'motifC': [True, False, False, False, False, True],
            'sox': [False, False, False, True, False, True],
            'pouf': [True, False, True, False, False, True],
            'value': [15, 30, 8, 12, 5, 11],
            'number_of_pro': [1,2,3,2,2,1]}
    df = pd.DataFrame(data)

    # create a dataframe with the regions that have a motif
    df['motif'] = df[['motifA', 'motifB', 'motifC']].apply(
        lambda x: ''.join(x.index[x]), axis=1)

    ###############bar plot################################
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
    plt.show()
    plt.savefig("bar_plot_has_motif.png")


    ###############box plot################################


    title = "Box plot"
    sns.boxplot(x='motif', y='value',
                data=df, palette='pastel')
    sns.stripplot(x='motif', y='value', hue='number_of_pro',
                  data=df, jitter=True, dodge=True, palette='pastel', size=5,
                  marker='o')
    plt.xlabel("Motifs")
    plt.title(title)
    plt.show()
    plt.savefig("box_plot.png")


    ###############bar plot2################################

    # Create a bar for the total number of true values
    motifs = ["motifA", "motifB", "motifC"]
    proteins = ["sox", "pouf"]

    motif_total_count = []
    motif_both_proteins_counts = []
    for motif in motifs:
        motif_total_count.append(df.query(f'{motif} == True').shape[0])
        motif_both_proteins_counts.append(df.query(f'{motif} == True & sox == True & pouf == True').shape[0])

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
    plt.ylabel('number of picks')
    plt.title('number of picks with sox and pouf  for each motif of nanog')
    plt.xticks(ind, ('motifA', 'motifB', 'motifC'))
    plt.legend((p_sox,p_pouf,p_both), ('sox', 'pouf', 'both'))

    plt.show()
    plt.savefig("bar_plot_count.png")


    #
    #
    # motif_counts = df[df['motifA'] == True].count()
    # motif_counts = df[df == True].count()
    #
    # motif_protein_counts = df.groupby(["motifA", "motifB", "motifC"])[
    #     "number_of_pro"].value_counts().unstack(fill_value=0)
    # colors = cm.rainbow(np.linspace(0, 1, 4))
    # for motif in ["motifA", "motifB", "motifC"]:
    #     plt.bar(motif, motif_counts[motif], color='gray')
    #     # Create bars for the number of true values for each protein count
    #     for i in range(4):
    #         plt.bar(motif, motif_protein_counts[motif][i],
    #                 bottom=sum(motif_protein_counts[motif][:i]),
    #                 color=colors[i])
    # # Show the plot
    # plt.show()



    #
    # # filter the dataframe to only include rows where the motifA is equal to True
    # df_motifA = df[df["motifA"] == True]
    #
    # # filter the dataframe to only include rows where the motifB is equal to True
    # df_motifB = df[df["motifB"] == True]
    #
    # # filter the dataframe to only include rows where the motifC is equal to True
    # df_motifC = df[df["motifC"] == True]
    #
    # # concatenate the three dataframe
    # df_motifs = pd.concat([df_motifA, df_motifB, df_motifC])
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
    # plt.show()
    #
    # # fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))
    # # for i,df in enumerate([df_motifA, df_motifB, df_motifC]):
    # #     sns.boxplot(y='value',
    # #                 data=df, palette='pastel', ax=ax[i])
    # #     sns.stripplot(y='value', hue='number_of_pro',
    # #                   data=df, jitter=True, dodge=True, palette='pastel', size=5,
    # #                   marker='o', ax=ax[i])
    # # # plt.xlabel(["MotifA", "MotifB", "MotifC"])
    # # plt.ylabel("Value")
    # # plt.show()
    #
    #
    #
    #
    # # title = "Box plot"
    # #
    # # motifs = ['motifA', 'motifB', 'motifC']
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
    # # plt.show()
    # # # filter the dataframe to only include rows where the motifA is equal to True
    # # df_motifA = df[df["motifA"] == True]
    # #
    # # # create the boxplot for motifA
    # # bp = sns.boxplot(x="motifA", y="value",
    # #                  hue="number_of_pro", data=df_motifA,
    # #                  showfliers=True, fliersize=2)
    # #
    # # # loop through the box, whiskers and caps of the boxplot and adjust the color and transparency
    # # for patch in bp.artists:
    # #     r, g, b, a = patch.get_facecolor()
    # #     patch.set_facecolor((r, g, b, .3))
    # #
    # # # add title and axis labels
    # # plt.title("MotifA")
    # # plt.xlabel("MotifA")
    # # plt.ylabel("Value")
    # #
    # # # show the plot
    # # plt.show()
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

