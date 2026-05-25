import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib import ticker

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
data.individual = data.individual.str.replace("A","O").str.replace("T","Y-Tx")
data.individual

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
data[
    [
        "individual",
        "subset",
        "chain",
        "clonetype",
        "Cell Number",
        "Group",
        "Read count",
        "UMI count",
        'Read/Cell(UMI)',
        "UMI/Cell(UMI)",
    ]
].to_csv("/home/erdem/Dropbox/Research/thymectomy/results_overview_naive_cleaned.csv")

# sep. dataframe into Treg and all
data_treg = data[data.subset == "CD4Treg"].copy()
data_f = data[data.subset != "CD4Treg"].copy()

def plot_overview(subset_data, plotName = "all"):
    for c in ["TRA", "TRB"]:
        # for subset in subsets:
        subset_data_cur = subset_data[(subset_data.chain == c)].copy()
        subset_data_cur["subset"] = subset_data_cur["subset"].cat.remove_unused_categories().copy()
        print(f"Plotting = {list(subset_data_cur.subset.unique())}")
        print(subset_data_cur.columns)
        # plot arguments
        ## plots Reads (A), UMI count (B), Reads/Cell (C), UMI/Cell (D)
        plots_dict = {"A": ["UMIorCell","Read count"],
                      "B": ["UMIorCell","UMI count"],
                      "C": ["UMIorCell","Read/Cell(UMI)"],
                      "D": ["UMIorCell","UMI/Cell(UMI)"],
                      "E": ["UMI/Cell(UMI)","Read/Cell(UMI)"]}
        color_list_current = color_list_cell
        current_color_group = "subset"
        f_width = 5.58
        f_height = 9.9409449
        fig, axs = plt.subplot_mosaic("ABE;CDE", figsize=(f_height, f_width))
        
        # Set negative vertical spacing to reduce height between plots
        plt.subplots_adjust(wspace = 0.4, hspace=0.4)#, left=0.15, right=0.95, top=0.95, bottom=0.05)
        
        scatter_args = {
            "s": 26,
            "alpha": 0.7,
            "clip_on": False,
            "zorder": 10,
            "linewidth": 0,
        }
        # start plotting
        for plot in plots_dict.keys():
            sns.scatterplot(
                data=subset_data_cur,
                x=plots_dict[plot][0],
                y=plots_dict[plot][1],
                ax=axs[plot],
                hue=current_color_group,
                style="individual",
                legend=1,
                palette=color_list_current,
                **scatter_args
                )
                

            if plot in ["A", "C"]:
                axs[plot].set_yscale("log")
                axs[plot].ticklabel_format(axis='x', style='sci', scilimits=(0,0), useMathText=True)
                axs[plot].xaxis.get_offset_text().set_visible(False)
                # Disable minor ticks for log plots
                axs[plot].yaxis.set_minor_formatter(ticker.NullFormatter())
                axs[plot].yaxis.set_minor_locator(ticker.NullLocator())
            elif plot in ["B"]:
                axs[plot].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
                axs[plot].xaxis.get_offset_text().set_visible(False)
                axs[plot].yaxis.get_offset_text().set_visible(False)  # Make the scientific notation smaller
                axs[plot].ticklabel_format(axis='x', style='sci', scilimits=(0,0), useMathText=True)
                axs[plot].yaxis.set_minor_formatter(ticker.NullFormatter())
                axs[plot].yaxis.set_minor_locator(ticker.NullLocator())
            else:
                axs[plot].yaxis.set_major_locator(ticker.MaxNLocator(integer=True, steps=[1, 2, 5, 10]))
                axs[plot].ticklabel_format(axis='x', style='sci', scilimits=(0,0), useMathText=True)
                # For linear scale, start from 0
                axs[plot].set_ylim(bottom=1)

            # Reduce x and y-axis tick label size
            axs[plot].tick_params(axis='both', labelsize=8)

            if plot in ["D", "E"]:
                axs[plot].axhline(y=1, color="black", linestyle="dotted", linewidth=0.7, zorder=0)
                axs[plot].axhline(y=1, color="black", linestyle="dotted", linewidth=0.7, zorder=0)
            if plot != "D":
                axs[plot].legend_.remove()
            axs[plot].set_xlabel('')
            # specific requirements for plots
            axs["C"].set_ylim(bottom=1, top = 1000)
            #axs["C"].axhline(y=1, color="black", linestyle="dotted", linewidth=0.7, zorder=0)
            #axs["E"].axvline(x=1, color="black", linestyle="dotted", linewidth=0.7, zorder=0)

            axs["E"].set_yscale("log")
            # Add faint grid
            axs[plot].grid(True, which='major', linestyle='--', alpha=0.2, color='gray')
            # add label outside of the plot
            axs[plot].text(
                -0.2,  # adjust left position
                0.99,   # adjust top position
                plot,
                transform=axs[plot].transAxes,
                fontsize=12,
                fontweight="bold",
                va="bottom",
                ha="right",
            )
        # add x axis label
        axs["A"].set_ylabel('# of Reads')
        axs["B"].set_ylabel('# of UMIs (x$10^5$)')
        axs["C"].set_ylabel('Reads/Cell (or UMI)')
        axs["D"].set_ylabel('UMIs/Cell (or UMI)')
        axs["E"].set_ylabel('Reads/Cell (or UMI)')
        axs["C"].set_xlabel('Cell Number (or UMI) (x$10^5$)')
        axs["D"].set_xlabel('Cell Number (or UMI) (x$10^5$)')
        axs["E"].set_xlabel('UMIs/Cell (or UMI)')

        # x = y
        xlim = axs["B"].get_xlim()
        ylim = axs["B"].get_ylim()

        axs["B"].axline((0, 0), slope=1, ls="--", c="gray", alpha=0.7, zorder=0)

        axs["B"].set_xlim(xlim)
        axs["B"].set_ylim(ylim)
        # legend
        handles, labels = axs["D"].get_legend_handles_labels()

        legend = fig.legend(
            handles,
            labels,
            loc='upper left',           # anchor point of the legend box
            bbox_to_anchor=(0.9, 0.75),  # (x, y) position in figure coordinates (0-1)
            fontsize=6,
            frameon=False,
            ncol=1
        )

        rename_map = {"subset": "Subset", "individual": "Individual"}

        for text in legend.get_texts():
            old_label = text.get_text()
            if old_label in rename_map:
                text.set_text(rename_map[old_label])
                text.set_ha('center')  # center the group header
                #text.set_position((text.get_position()[0] + 10, text.get_position()[1]))  # adjust x offset

        # Make room on the right for the legend
        #fig.subplots_adjust(right=0.99)  # adjust as needed

        # Remove legend from subplot
        axs["D"].legend_.remove()
        fig.subplots_adjust(wspace = 0.3)
        fig.savefig(
            f"../../figures/sample_metrics_scatter_{c}_{plotName}_naive_cleaned.pdf",
            dpi=600,
            #pad_inches=0.15,
            bbox_inches='tight'
        )
        mpl.pyplot.close()
plot_overview(data_f, plotName="thymectomy")
plot_overview(data_treg, plotName="Treg")
