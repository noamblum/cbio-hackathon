import pandas as pd
import matplotlib.pyplot as plt
if __name__ == '__main__':

    # Prepare the data
    data = {'motif': ['A', 'A', None, 'B', None, 'C'],
           'binary_var': ['true', 'false', 'true', 'false', 'true', 'false'],
            'value': [15, 10, 8, 12, 7, 11],
            'number_of_pro': [1,2,3,2,2,1]}
    df = pd.DataFrame(data)

    # create a dataframe with the regions that have a motif

    title = "Barplot of pick region with motif and without motif"
    motif_regions = df[df['motif'].isnull() == False]
    count_motif_regions = len(motif_regions)

    # create a dataframe with the regions that do not have a motif
    no_motif_regions = df[df['motif'].isnull()]
    count_no_motif_regions = len(no_motif_regions)

    # create a bar plot
    plt.bar(['Motif Regions', 'No Motif Regions'],
            [count_motif_regions, count_no_motif_regions],
            color=['blue', 'red'])
    plt.legend(["Motif", "No Motif"])
    plt.xlabel("Region Type")
    plt.ylabel("Count")
    plt.title(title)
    plt.show()

    #
    #
    # df_true = df[df['binary_var'] == 'true']
    # df_false = df[df['binary_var'] == 'false']
    #
    # # Create the bar plot
    # fig, ax = plt.subplots()
    # ax.bar(df_true['group'], df_true['count'], width=-0.2, align='edge',
    #        color='blue', label = "True")
    # ax.bar(df_false['group'], df_false['count'], width=0.2, align='edge',
    #        color='red', label = "False")
    # ax.set_xlabel('Motif')
    # ax.set_ylabel('Count')
    # ax.legend(title = "")
    # plt.show()
    #
    #
    #
    # # Create the box plot
    # fig, ax = plt.subplots()
    # colors = df['number_of_pro']
    # motif = df['group']
    # pick_val = df['value']
    # scatter = ax.scatter(motif, pick_val, c=colors, cmap='Blues',
    #            alpha=0.5)
    #
    # # Add a colorbar
    # cbar = plt.colorbar(scatter)
    # cbar.set_label('Number of proteins')
    #
    # # Add x-axis labels
    # ax.set_xticklabels(df['group'].unique())
    #
    # # Show the plot
    # plt.show()







