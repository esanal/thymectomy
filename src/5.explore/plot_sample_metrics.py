'''
Fix UMI/Cell, possibly via copying methods (or maybe the data modified in my other script where)
'''


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib import ticker
import numpy as np

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
font = {
    "family": "Myriad Pro",
    #'weight' : 'bold',
    "size": 8,
}
mpl.rc("font", **font)
plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
plt.rcParams['axes.xmargin'] = 0.001
plt.rcParams['axes.ymargin'] = 0.001
mpl.rcParams.update({'figure.autolayout': True})
cm = 1 / 2.54

# data read
pickles_dir_cleaned = "../../results/mergedGSIDs_mixcr_th1_collapsed_clones_naives_cleaned.pkl"
data = pd.read_pickle(pickles_dir_cleaned)

# UMI/Cell calculation
data["UMI/Cell(UMI)"] = data["UMI count"] / data["UMIorCell"]

subsets = data.subset.unique()

# Colors
color_list = {"Young": "#e97131", "Old": "#9c0054", "ThymectomizedYoung": "#3957bd"}
color_list_cell = {
    "CD4CM": "#006d2c",    # Forest Green
    "CD4EM": "#41ae76",    # Seafoam
    "CD4NCD31": "#88419d", # Purple
    "CD4NCD31-": "#8c96c6",# Lavender
    "CD8EM": "#c51b7d",    # Magenta
    "CD8N": "#2b8cbe",     # Teal
    "CD4Treg": "black"   # Graphite
}

# Ind order
data.individual = pd.Categorical(
    data.individual,
    categories=[
        "Y-Tx2",
        "Y-Tx4",
        "Y-Tx5",
        "Y-Tx6",
        "Y-Tx9",
        "Y-Tx10",
        "Y1",
        "Y3",
        "Y4",
        "Y5",
        "Y6",
        "Y7",
        "O5",
        "O6",
        "O9",
        "O10",
        "O12",
        "O13",
    ],
    ordered=True,
)
order_ind = [
    "Y-Tx2",
    "Y-Tx4",
    "Y-Tx5",
    "Y-Tx6",
    "Y-Tx9",
    "Y-Tx10",
    "Y1",
    "Y3",
    "Y4",
    "Y5",
    "Y6",
    "Y7",
    "O5",
    "O6",
    "O9",
    "O10",
    "O12",
    "O13",
]
order_ind_treg = ["Y-Tx2", "Y-Tx4", "Y-Tx5", "Y-Tx6", "Y-Tx9", "Y-Tx10", "Y3", "Y4", "Y5", "Y6", "Y7", "O5", "O6", "O9", "O10", "O12", "O13"]

# sep. dataframe into Treg and all
data_treg = data[data.subset == "CD4Treg"].copy()
data = data[data.subset != "CD4Treg"].copy()
data[
    [
        "individual",
        "subset",
        "chain",
        "Group",
        "Cell Number",
        "Cell Number Estimated",
        "UMIorCell_code",
        "Read count",
        "UMI count",
        "clonetype",
        "Read/Cell(UMI)",
        "UMI/Cell(UMI)",
        "Clonetype/Cell(UMI)"
   ]
].to_csv("../../results/explore_thymectomy_overview.csv")

def plot_overview(subset_data, plotName = "all", order_ind = order_ind):
    for c in ["TRA", "TRB"]:
        # for subset in subsets:
        subset_data_cur = subset_data[(subset_data.chain == c)].copy()
        subset_data_cur["subset"] = subset_data_cur["subset"].cat.remove_unused_categories().copy()
        print(f"Plotting = {list(subset_data_cur.subset.unique())}")
        print(subset_data_cur.columns)
        # plot arguments
        ## plots Cell number (A), Reads (B), UMI count (C), clonetype (D), Reads/Cell (E), UMI/Cell (F)
        plots_dict = {"A": "Cell Number",
                      "B": "Read count",
                      "C": "UMI count",
                      "D": "clonetype",
                      "E": "Read/Cell(UMI)",
                      "F": "UMI/Cell(UMI)",
                      "G": "Clonetype/Cell(UMI)"}
        color_list_current = color_list_cell
        current_color_group = "subset"
        f_width = 7.3
        f_height = 9
        fig, axs = plt.subplot_mosaic("AAA;BBB;CCC;DDD;EEE;FFF;GGG", figsize=(f_width, f_height))
        if plotName == "Treg":
            order_ind = order_ind_treg
        else:
            order_ind = order_ind
        # start plotting
        for plot in plots_dict.keys():
            sns.barplot(
                data=subset_data_cur,
                x="individual",
                y=plots_dict[plot],
                hue="subset",
                ax=axs[plot],
                palette=color_list_current,
                order=order_ind
            )
            axs[plot].set(xlabel='Individual')
            axs[plot].spines['top'].set_visible(False)
            axs[plot].spines['right'].set_visible(False)
            if plot in ["A", "B", "C", "D", "E"]:
                axs[plot].set_yscale("log")
                axs[plot].yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=10))
                # For log scale, set minimum to small value
                axs[plot].set_ylim(bottom=0.1)
                # Disable minor ticks for log plots
                axs[plot].yaxis.set_minor_formatter(ticker.NullFormatter())
                axs[plot].yaxis.set_minor_locator(ticker.NullLocator())
            elif plot in ["G"]:
                axs[plot].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
                axs[plot].yaxis.get_offset_text().set_fontsize(6)  # Make the scientific notation smaller
            else:
                axs[plot].yaxis.set_major_locator(ticker.MaxNLocator(integer=True, steps=[1, 2, 5, 10]))
                # For linear scale, start from 0
                axs[plot].set_ylim(bottom=0)
            # Reduce y-axis tick label size
            axs[plot].tick_params(axis='y', labelsize=6)
            if plot in ["E", "F"]:
                axs[plot].axhline(y=1, color="black", linestyle="dotted", linewidth=0.7, zorder=0)
            if plot != "G":
                axs[plot].legend_.remove()
                axs[plot].set_xticks([])
            axs[plot].set_xlabel('')
            axs["A"].set_ylabel('# of Cells (x$10^5$)')
            axs["B"].set_ylabel('# of Reads')
            axs["C"].set_ylabel('# of UMIs')
            axs["D"].set_ylabel('# of Clonotypes')
            axs["E"].set_ylabel('Reads/Cell(or UMI)')
            axs["F"].set_ylabel('UMIs/Cell(or UMI)')
            axs["G"].set_ylabel('Clonetypes/Cell(or UMI)')

            axs["G"].set_ylim(0, 1)

            # change x-axis label rotation    
            axs[plot].tick_params(axis='x', rotation=45)
            # Add faint grid
            axs[plot].set_axisbelow(True)
            axs[plot].grid(True, which='major', axis='y', linestyle='--', alpha=0.2, color='gray')
            # add label outside of the plot
            axs[plot].text(
                -0.1,  # adjust left position
                0.85,   # adjust top position
                plot,
                transform=axs[plot].transAxes,
                fontsize=12,
                fontweight="bold",
                va="bottom",
                ha="left",
            )

        # Create a separate figure for the legend
        figlegend = plt.figure(figsize=(f_width, 1))
        # Get the legend from the last plot and add it to the new figure
        legend_elements = axs["G"].get_legend_handles_labels()
        figlegend.legend(legend_elements[0], legend_elements[1], 
                        loc='center', frameon=False, 
                        fontsize=5, title_fontsize=5,
                        ncol=2)  # adjust ncol as needed
        figlegend.tight_layout()
        
        # Save the legend separately
        figlegend.savefig(f'../../figures/sample_metrics_legend_{plotName}.pdf', bbox_inches='tight', dpi=300)
        plt.close(figlegend)
        
        # Legend
        handles, labels = axs["G"].get_legend_handles_labels()

        fig.legend(
            handles, 
            labels,
            loc='upper center',
            bbox_to_anchor=(0.5, 0.02),  # centered horizontally, near bottom
            ncol=len(labels),            # all items in one row
            fontsize=5,
            frameon=False
        )

        # Make room for the legend at the bottom
        fig.subplots_adjust(bottom=0.08)

        # Remove the legend from the subplot itself
        axs["G"].legend_.remove()
        
        fig.savefig(
            f"../../figures/sample_metrics_{c}_{plotName}_naive_cleaned.pdf",
            dpi=600
        )
        # Report UMI>1 samples
        n_umi_larger_1 = len(subset_data[subset_data["UMI/Cell(UMI)"] > 1])
        n_umi_less_1 = len(subset_data[subset_data["UMI/Cell(UMI)"] < 1])
        print(
            f"total size = {len(subset_data)}, UMI/Cell(UMI)>1 = {n_umi_larger_1}, UMI/Cell(UMI)<1 = {n_umi_less_1}"
        )
        mpl.pyplot.close()
plot_overview(data, plotName="thymectomy")
plot_overview(data_treg, plotName="Treg", order_ind=order_ind_treg)
