import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib import ticker
import math
import numpy as np
from brokenaxes import brokenaxes
from matplotlib.ticker import MultipleLocator
from scipy import stats

mpl.rcParams["pdf.fonttype"] = 42
plt.rcParams.update({
    'font.size': 7,
    'font.family': 'Arial',  # or 'Helvetica', 'DejaVu Sans'
    'axes.linewidth': 1.0,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1,
    'xtick.major.size': 0,
    'ytick.major.size': 4,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'axes.spines.right': False,
    'axes.spines.top': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'xtick.labelbottom': False
})

plt.rcParams['path.simplify'] = True
plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
plt.rcParams['axes.xmargin'] = 0.001
plt.rcParams['axes.ymargin'] = 0.001
mpl.rcParams.update({'figure.autolayout': True})
cm = 1 / 2.54

# data read (collapsed clones already!)
pickles_dir = "~/Hosts/vacuole1/former-NOBINFBACKUP/thymectomy/results/pickles/mixcr_th1_raw_collapsed_clones.pkl"
data_f = pd.read_pickle(pickles_dir)

data_old = data_f.copy()

subsets = data_f.subset.unique()

def clean_by_lineage(group):
    # {Target_Subsets: Source_Subsets_to_Subtract}
    rules = {("CD4NCD31", "CD4NCD31-"): ["CD4CM", "CD4EM"], ("CD8N",): ["CD8EM"]}

    # Clean in the loop
    for targets, sources in rules.items():
        # Collect sequences to remove for the specific rule
        to_remove = set()
        source_rows = group[group["subset"].isin(sources)]
        for seq_list in source_rows["nSeqCDR3"]:
            to_remove.update(seq_list)

        if not to_remove:
            continue

        # Identify which rows in the group match the targets of this rule
        target_indices = group[group["subset"].isin(targets)].index

        for idx in target_indices:
            row = group.loc[idx]

            # Filter indices
            keep_idx = [
                i for i, seq in enumerate(row["nSeqCDR3"]) if seq not in to_remove
            ]

            # If no changes needed, skip to save time
            if len(keep_idx) == len(row["nSeqCDR3"]):
                continue

            # Update list columns
            list_cols = [
                "nSeqCDR3",
                "umi_counts",
                "umi_tag_read",
                "readCount_list",
                "v",
                "d",
                "j",
                "c",
            ]
            for col in list_cols:
                group.at[idx, col] = [row[col][i] for i in keep_idx]

            # Update cell number based on removed UMIs (estimated!!!)
            group.at[idx, "Cell Number"] = round(
                (group.at[idx, "Cell Number"] / group.at[idx, "UMI count"])
                * sum(group.at[idx, "umi_counts"])
            )
            group.at[idx, "Cell Number Estimated"] = True
            # Update metrics
            new_nseq = group.at[idx, "nSeqCDR3"]
            group.at[idx, "clonetype"] = len(new_nseq)
            group.at[idx, "UMI count"] = sum(group.at[idx, "umi_counts"])
            group.at[idx, "Read count"] = sum(group.at[idx, "readCount_list"])
    return group


# --- Run the update ---
data_f = data_f.groupby(["individual", "chain"], group_keys=False).apply(
    clean_by_lineage
)

# Colors
# color_list = {"Young": "#e96f00", "Aged": "#b90000", "Thymectomized": "#0032ff"}
color_list = {"Young": "#e97131", "Aged": "#9c0054", "Thymectomized": "#3957bd"}
color_list_cell = {
    "CD4CM": "#223470",
    "CD4EM": "#3957bd",
    "CD4NCD31": "#f0a931",
    "CD4NCD31-": "#e97131",
    "CD8EM": "#904700",
    "CD8N": "#9c0054",
    "CD4Treg": "black",
}

# Individual order
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

# titles to print
titles = ["CD4CM", "CD4EM", "CD4NCD31", "CD4NCD31-", "CD4Treg", "CD8EM", "CD8N"]
titles_to_print = {
    "CD4CM": r"$\mathregular{CD4^+ T_{CM}}$",
    "CD4EM": r"$\mathregular{CD4^+ T_{EM}}$",
    "CD4NCD31": r"$\mathregular{CD4^+ CD31^+ Naive}$",
    "CD4NCD31-": r"$\mathregular{CD4^+ CD31^- Naive}$",
    "CD4Treg": r"$\mathregular{CD4^+ T_{reg}}$",
    "CD8EM": r"$\mathregular{CD8^+ T_{EM}}$",
    "CD8N": r"$\mathregular{CD8^+ Naive}$"
}
titles = [titles[i] for i in [6, 2, 3, 5, 1, 0]]

data_f.individual = data_f.individual.str.replace("A","O").str.replace("T","Y-Tx")

data_f.individual = pd.Categorical(
    data_f.individual,
    categories=order_ind,
    ordered=True,
)


# Calculate clone per cell
data_f["Clonetype/Cell"] = data_f["clonetype"] / data_f["Cell Number"]
data_f["Clonetype/Cell(UMI)"] = np.where(
    data_f["UMI count"] < data_f["Cell Number"],
    data_f["clonetype"] / data_f["UMI count"],
    data_f["clonetype"] / data_f["Cell Number"],
)

# Read quality: checks for if there are not enough read per cell
data_f["Read/Cell"] = data_f["Read count"] / data_f["Cell Number"]
data_f["read_quality"] = np.where(data_f["Read/Cell"] >= 1, 1, 0)
data_f["Read/Cell(UMI)"] = np.where(
    data_f["UMI count"] < data_f["Cell Number"],
    data_f["Read count"] / data_f["UMI count"],
    data_f["Read count"] / data_f["Cell Number"],
)
data_f["read_quality_umi_included"] = np.where(data_f["Read/Cell(UMI)"] >= 1, 1, 0)

# Proxy for the cell numbers: Cell number or UMI
data_f["UMIorCell"] = np.where(
    data_f["UMI count"] < data_f["Cell Number"],
    data_f["UMI count"],
    data_f["Cell Number"],
)
data_f["UMIorCell_code"] = np.where(
    data_f["UMI count"] >= data_f["Cell Number"],
    "cell",
    "umi"
)

# Calculate % of abundance top 10 clones
## sort clones
data_f['umi_counts'].apply(lambda x: sorted(x, reverse=True))
## find total abundance
totals = data_f['umi_counts'].apply(sum)
## calculate percentages
data_f['umi_count_percent'] = [[100*x/total for x in lst] for lst,total in zip(data_f['umi_counts'], totals)]
data_f['umi_count_percent_10th'] = [np.cumsum(lst)[9] for lst in data_f['umi_count_percent']]

data_f[["individual", "chain", "subset", "umi_count_percent", "umi_counts"]].to_csv('~/Sync/thymectomy_clone_dist_data_collapsed_naive_cleanup.csv', index=False)

# export data  with collapsed clones and cleaned up naives!!!
pickles_dir_cleaned = "~/Hosts/vacuole1/former-NOBINFBACKUP/thymectomy/results/pickles/mixcr_th1_raw_collapsed_clones_naives_cleaned.pkl"
data_f.to_pickle(pickles_dir_cleaned)

# Remove samples with less than 500 umis!
data_f_removed =data_f[data_f['UMI count'] <= 500]
print(f"Removed samples with less than 500 UMIs: {data_f_removed.shape[0]}")
data_f = data_f[data_f['UMI count'] >= 500].copy()




# sep. dataframe into Treg and all
data_treg = data_f[data_f.subset == "CD4Treg"].copy()
data_f = data_f[data_f.subset != "CD4Treg"].copy()

# export data
data_f[["individual", "chain", "subset", 'Clonetype/Cell(UMI)', 'umi_count_percent_10th', 'UMIorCell']].to_csv('thymectomy_clone_dist_data_collapsed_naive_cleanup.csv', index=False)

order = ["Thymectomized", "Young", "Aged"]
x_positions = [0, 1, 2]

# Define brokenaxes limits
broken_y_lims = {
    # "CD8N":      {"ylim_bottom": (0, 5+1), "ylim_top": (55-1, 60), "tick_step_bottom": 2.5,   "tick_step_top": 2.5},
    # "CD4NCD31":  {"ylim_bottom": (0, 5+1), "ylim_top": (7.5-1, 10), "tick_step_bottom": 2.5,  "tick_step_top": 2.5},
    "CD4CM":     {"ylim_bottom": (0, 4+1), "ylim_top": (58-1, 66), "tick_step_bottom": 2,   "tick_step_top": 2},
}


def plot_points(data, column_to_plot, ax_to_plot, filled, clipping = False):
    if data.empty:
        return
    
    for group_idx, group in enumerate(order):
        group_data = data[data["Group"] == group]
        if group_data.empty:
            continue
        
        y_values = group_data[column_to_plot].values
        #breakpoint()
        individuals = group_data.individual.values
        # Sort by y-value to process from bottom to top
        sorted_indices = np.argsort(y_values)
        sorted_y = y_values[sorted_indices]
        
        # Sort individuals
        sorted_individuals = individuals[sorted_indices]
        
        # Calculate minimum y-distance that requires x-adjustment
        # Points closer than this in y need different x positions
        y_range = ax_to_plot.get_ylim()
        min_y_distance = (y_range[1] - y_range[0]) * 0.001  # scale
        
        # Generate x-positions that avoid overlap
        x_positions = np.zeros(len(sorted_y))
        
        for i in range(len(sorted_y)):
            x_offset = 0
            collision = True
            max_attempts = 50
            attempt = 0
            
            while collision and attempt < max_attempts:
                collision = False
                test_x = group_idx + x_offset
                
                # Check for collisions with previously placed points
                for j in range(i):
                    y_diff = abs(sorted_y[i] - sorted_y[j])
                    x_diff = abs(test_x - x_positions[j])
                    
                    # If points are close in y and x, they overlap
                    if x_diff < 0.08:
                        collision = True
                        break
                
                if collision:
                    # Alternate between left and right, increasing offset
                    if attempt % 2 == 0:
                        x_offset = (attempt // 2 + 1) * 0.04
                    else:
                        x_offset = -(attempt // 2 + 1) * 0.04
                    attempt += 1
                else:
                    x_positions[i] = test_x
                    break
            
            # Fallback if no position found
            if attempt >= max_attempts:
                x_positions[i] = group_idx + np.random.normal(0, 0.04)
        
        # Marker for Y-Tx 10
        marker_list = ["*" if ind == "Y-Tx10" else "o" for ind in sorted_individuals] 

        # Plot points with original y-values and adjusted x-positions
        for xi, yi, mkr in zip(x_positions, sorted_y, marker_list):
            # Adjust size: stars ('*') look smaller than circles ('o')
            current_s = 30 if mkr == "*" else 9        
            if filled:
                ax_to_plot.scatter(
                    x=xi,
                    y=yi,
                    marker=mkr,
                    color=color_list[group],
                    edgecolors="black",
                    linewidths=0.5,
                    s=current_s,
                    alpha=1,
                    clip_on=clipping,
                    zorder=100
                )
            else:
                ax_to_plot.scatter(
                    x=xi,
                    y=yi,
                    marker=mkr,
                    facecolors='white',
                    edgecolors="black",
                    s=current_s,
                    linewidths=0.5,
                    alpha=1,
                    clip_on=clipping,
                    zorder=100
                )

# Broken axes plotter
def plot_broken_axes(
    original_ax,
    subset_selected_data,
    ylim_bottom=(0, 10),
    ylim_top=(50, 60),
    x_positions=x_positions,
    order=order,
    color_list=color_list,
    y_label="",
    show_y_tick_labels=True,
    show_x_tick_labels = False,
    tick_step_bottom=None,
    tick_step_top=None,
    rng=None):
    data = subset_selected_data.copy()
    filled_data = data[data["UMIorCell_code"] != "umi"].copy()
    hollow_data = data[data["UMIorCell_code"] == "umi"].copy()

    # 1. Reset the original axis but KEEP it in the figure
    original_ax.clear()
    original_ax.axis("off") 
    
    # 2. Use its space to create the split
    ## calculate the range (span) of each axis
    range_top = ylim_top[1] - ylim_top[0]
    range_bot = ylim_bottom[1] - ylim_bottom[0]

    gs = original_ax.get_subplotspec().subgridspec(
        2, 1, 
        hspace=0.15, 
        height_ratios=[range_top, range_bot]
    )

    ax_top = original_ax.figure.add_subplot(gs[0])
    ax_bot = original_ax.figure.add_subplot(gs[1])

    def plot_data_on(ax_in, df, df_filled, df_hollow):
        
        # Calculate medians 
        medians = df.groupby("Group", observed=True)["umi_count_percent_10th"].median()
        
        # Ensure we only plot if the index exists to avoid 'KeyError' or plotting old data
        available_medians = [medians[g] if g in medians.index else 0 for g in order]

        ax_in.bar(
            x_positions,
            available_medians,
            color=[color_list[group] for group in order],
            width=0.45,
            edgecolor="none",
            alpha=0.7,
        )
        
        plot_points(df_filled, "umi_count_percent_10th", ax_in, filled=True, clipping=True)
        plot_points(df_hollow, "umi_count_percent_10th", ax_in, filled=False, clipping=True)

    # plot
    plot_data_on(ax_top, data, filled_data, hollow_data)
    plot_data_on(ax_bot, data, filled_data, hollow_data)

    # configure the break
    ylim_bottom
    ax_bot.set_ylim(*ylim_bottom)  # bottom range
    ax_top.set_ylim(*ylim_top)  # top range

    # tick density
    if tick_step_bottom is not None:
        ax_bot.yaxis.set_major_locator(MultipleLocator(tick_step_bottom))
        #ax_bot.yaxis.set_major_locator(MaxNLocator(nbins='auto', prune='upper'))
    if tick_step_top is not None:
        ax_top.yaxis.set_major_locator(MultipleLocator(tick_step_top))

    # hide spines between
    ax_bot.spines["top"].set_visible(False)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.xaxis.tick_top()
    ax_top.tick_params(labeltop=False)  # don't put tick labels at the top
    ax_bot.xaxis.tick_bottom()

    # diagonal cut lines 
    d = 0.015  # size of cut line
    kwargs = dict(transform=ax_top.transAxes, color="k", clip_on=False, linewidth=1)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal

    kwargs.update(transform=ax_bot.transAxes)  # switch to the bottom axes
    ax_bot.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

    if show_y_tick_labels == False:
        ax_top.tick_params(axis="y", labelleft=False)
        ax_bot.tick_params(axis="y", labelleft=False)

    if show_x_tick_labels == True:
        ax_top.tick_params(axis="x", labelbottom=False)
        ax_bot.tick_params(axis="x", labelbottom=True)
        ax_bot.set_xticks(x_positions)
        ax_bot.set_xticklabels(["Y-Tx", "Y", "O"], rotation=0, weight="bold")

    ax_top.set_ylabel("")
    if y_label != "":
        ax_top.set_ylabel(titles_to_print[y_label])
        ax_top.yaxis.set_label_coords(-0.25, 0.5, transform=original_ax.transAxes)

    ax_top.set_xlim(-0.5, 2.5)
    ax_bot.set_xlim(-0.5, 2.5)

    original_ax.set_ylabel("Centered Label", labelpad=20)
    original_ax.yaxis.set_label_position("left")
    original_ax.get_yaxis().set_visible(True) # Make just the label visible

    ax_top.yaxis.get_major_ticks()[0].tick1line.set_visible(False)
    ax_bot.yaxis.get_major_ticks()[-1].tick1line.set_visible(False)

    original_ax.top_split = ax_top
    original_ax.bot_split = ax_bot

def stats_groups(dfs, column = "Clonetype/Cell(UMI)"):
    # Statistics
    combinations = [
        ("Young", "Aged"),
        ("Young", "Thymectomized"),
        ("Thymectomized", "Aged"),
    ]
    pvalues = []
    pvalues_toplot = []
    combinations_toplot = []
    pvalues_toplot_corrected = []
    combinations_toplot_corrected = []
    for combination in combinations:
        # q=0
        data_1 = dfs[(dfs.Group == combination[0])][column].to_list()
        data_2 = dfs[(dfs.Group == combination[1])][column].to_list()

        # mannwhitneyu
        s_0, p = stats.mannwhitneyu(data_1, data_2)
        bonferroni_threshold = 0.05 / 3
        pvalues.append(p)
        pCutOff = 0.05
        if p <= pCutOff:
            pvalues_toplot.append(p)
            combinations_toplot.append(combination)
        if p <= bonferroni_threshold:
            pvalues_toplot_corrected.append(get_sig_code(p))
            combinations_toplot_corrected.append(combination)
    return(pvalues, pvalues_toplot, combinations_toplot, pvalues_toplot_corrected, combinations_toplot_corrected)

def get_sig_code(p):
    if p <= 0.0001:
        return '****'
    elif p <= 0.001:
        return '***'
    elif p <= 0.01:
        return '**'
    elif p <= 0.05:
        return '*'
    else:
        return 'ns'



def add_significance(ax, x1, x2, y, h, text):
    """
    ax: the axes object
    x1, x2: the x-positions of the two bars
    y: the height where the line starts
    h: the height of the 'tips' of the bracket
    text: the significance string (e.g., '***' or 'p < 0.01')
    """
    # Draw the bracket: left tip -> top line -> right tip
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.75, c='black')
    # Add the text/stars
    ax.text((x1+x2)*.5, y+h*0.5, text, ha='center', va='bottom', color='black', fontsize=7)

# Modified plot function
def plot_overview(subset_data, plotName="all"):
    # filter data
    subset_data_cur_alpha = subset_data[(subset_data.chain == "TRA")].copy()
    subset_data_cur_beta = subset_data[(subset_data.chain == "TRB")].copy()
    subset_data_cur_alpha["subset"] = subset_data_cur_alpha["subset"].cat.remove_unused_categories()
    subset_data_cur_beta["subset"] = subset_data_cur_beta["subset"].cat.remove_unused_categories()

    # Get unique groups
    unique_groups = subset_data_cur_alpha.subset.unique()
    n_groups = len(unique_groups)
    if n_groups > 1:
        unique_groups = titles
    
    n_rows = n_groups
    if n_rows == 6:
        height_ratios = [1, 1, 1, 0.001, 1, 1, 1]
        n_rows = 7
    else:
        height_ratios = [1]
    cm = 1/2.54
    
    # Create a single figure with 4 columns
    fig_width = 16 * cm  # Width for 4 columns
    fig_height = n_rows * 3.1 * cm
    
    # Create figure
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # Create GridSpec with more space between columns 1 and 2 (0-indexed: 1 is between col1 and col2)
    left_width = 1.0
    right_width = 1.0
    spacer_width = 0.5  # Adjust this to control space between column groups
    
    gs = fig.add_gridspec(
        n_rows, 5,  # 5 columns: 2 left, 1 spacer, 2 right
        width_ratios=[left_width, left_width, spacer_width, right_width, right_width],
        height_ratios=height_ratios,
        hspace=0.3,
        wspace=0.1  # Small spacing within each group
    )
    
    # Create axes - skip the spacer column (index 2)
    axes = []
    for i in range(n_rows):
        if i == 3:
            continue
        row_axes = []
        # Left group columns (0, 1)
        row_axes.append(fig.add_subplot(gs[i, 0]))
        row_axes.append(fig.add_subplot(gs[i, 1]))
        # Skip the spacer column (gs[i, 2])
        # Right group columns (3, 4)
        row_axes.append(fig.add_subplot(gs[i, 3]))
        row_axes.append(fig.add_subplot(gs[i, 4]))
        axes.append(row_axes)
    axes = np.array(axes)
    
    # Order age groups
    order = ["Thymectomized", "Young", "Aged"]
    x_positions = [0, 1, 2]
    # Order of subsets
    order_subsets = {
        "CD8N": {"richness": [0, 0], "pool": [3, 0]},
        "CD4NCD31": {"richness": [1, 0], "pool": [4, 0]},
        "CD4NCD31-": {"richness": [2, 0], "pool": [5, 0]},
        "CD8EM": {"richness": [0, 2], "pool": [3, 2]},
        "CD4EM": {"richness": [1, 2], "pool": [4, 2]},
        "CD4CM": {"richness": [2, 2], "pool": [5, 2]},
    }

    for i, subset in enumerate(unique_groups):
        # Select subset and read_quality
        subset_selected_alpha = subset_data_cur_alpha[
            (subset_data_cur_alpha["subset"] == subset) &
            (subset_data_cur_alpha["read_quality_umi_included"] == 1)
        ]
        subset_selected_beta = subset_data_cur_beta[
            (subset_data_cur_beta["subset"] == subset) &
            (subset_data_cur_beta["read_quality_umi_included"] == 1)
        ]

        # Stats
        stat_configs = [
        ("alpha", subset_selected_alpha, None),
        ("beta",  subset_selected_beta,  None),
        ("alpha_pool", subset_selected_alpha, "umi_count_percent_10th"),
        ("beta_pool",  subset_selected_beta,  "umi_count_percent_10th")
        ]

        results = {}

        for label, data, col in stat_configs:
            # Call the function (passing 'column' only if it exists)
            args = {"column": col} if col else {}
            res = stats_groups(data, **args)
            # Store the result and print the 5th element (the 'corrected' pvals)
            results[label] = res
            #print(f"{subset} {'Pool ' if 'pool' in label else 'Richness '}Pvals ({label.split('_')[0]}): {list(zip(res[2], get_sig_code(res[3])))}")
            print(f"{subset} {'Pool ' if 'pool' in label else 'Richness '}Pvals ({label.split('_')[0]}): {res[4]}")
            print(f"{subset} {'Pool ' if 'pool' in label else 'Richness '}Pvals ({label.split('_')[0]}): {res[3]}")

        # Subset alpha and beta subsets
        filled_data_alpha = subset_selected_alpha[subset_selected_alpha["UMIorCell_code"] != "umi"]
        filled_data_beta = subset_selected_beta[subset_selected_beta["UMIorCell_code"] != "umi"]
        hollow_data_alpha = subset_selected_alpha[subset_selected_alpha["UMIorCell_code"] == "umi"]
        hollow_data_beta = subset_selected_beta[subset_selected_beta["UMIorCell_code"] == "umi"]
        #breakpoint()
        
        # COLUMN 0: Alpha chain - Richness per cell (0-1)
        #ax0 = axes[i, 0]
        a, b = order_subsets[subset]["richness"]
        ax0 = axes[a, b]
        medians_alpha = subset_selected_alpha.groupby("Group", observed=True)["Clonetype/Cell(UMI)"].median()
        
        # Plot barsorder_subsets[subset]["richness"]
        ax0.bar(x_positions, medians_alpha[order].values,
            color=[color_list[group] for group in order],
            width=0.45, edgecolor='none', alpha=0.7)
        
        plot_points(filled_data_alpha, "Clonetype/Cell(UMI)", ax0, filled=True)
        plot_points(hollow_data_alpha, "Clonetype/Cell(UMI)", ax0, filled=False)

        # Set y-axis limits for left columns (0-1)
        ax0.set_ylim(0, 1.0)
        ax0.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax0.set_ylabel(titles_to_print[subset])  # Empty y-label
        ax0.set_xlim(-0.5, 2.5)
        
        # COLUMN 1: Beta chain - Richness per cell (0-1)
        ax1 = axes[a, b+1]
        medians_beta = subset_selected_beta.groupby("Group", observed=True)["Clonetype/Cell(UMI)"].median()

        ax1.bar(
            x_positions,
            medians_beta[order].values,
            color=[color_list[group] for group in order],
            width=0.45,
            edgecolor="none",
            alpha=0.7,
        )

        plot_points(filled_data_beta, "Clonetype/Cell(UMI)", ax1, filled=True)
        plot_points(hollow_data_beta, "Clonetype/Cell(UMI)", ax1, filled=False)

        
        # Set y-axis limits for left columns (0-1)
        ax1.set_ylim(0, 1.0)
        ax1.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax1.set_ylabel("")  # Empty y-label
        ax1.set_xlim(-0.5, 2.5)

        # COLUMN 2: Alpha chain - Top 10 clones (0-100) - Note: this is now index 2 in axes array
        a, b = order_subsets[subset]["pool"]
        ax2 = axes[a, b]
        if subset in broken_y_lims.keys():
            print(subset)
            plot_broken_axes(axes[i, 2], subset_selected_alpha, **broken_y_lims[subset], y_label=subset)
        else:
            medians_alpha_top = subset_selected_alpha.groupby("Group", observed=True)["umi_count_percent_10th"].median()            
            
            ax2.bar(x_positions, medians_alpha_top[order].values,
                color=[color_list[group] for group in order],
                width=0.45, edgecolor='none', alpha=0.7)
            #max_y = max([hollow_data_alpha["umi_count_percent_10th"].max(), filled_data_alpha["umi_count_percent_10th"].max()])

            #add_significance(ax2, 1, 2, max_y*1.1, 0.1, "**")
            #add_significance(ax2, 0, 2, max_y*1.25, 0.1, "*")
            plot_points(filled_data_alpha, "umi_count_percent_10th", ax2, filled=True)
            plot_points(hollow_data_alpha, "umi_count_percent_10th", ax2, filled=False)

        ax2.set_xlim(-0.5, 2.5)
        ax2.set_ylabel(titles_to_print[subset])
        
        # COLUMN 3: Beta chain - Top 10 clones (0-100) - Note: this is now index 3 in axes array
        ax3 = axes[a, b + 1]
        if subset in broken_y_lims.keys():
            plot_broken_axes(axes[i, 3], subset_selected_beta, show_y_tick_labels=False, **broken_y_lims[subset])
        else:
            medians_beta_top = subset_selected_beta.groupby("Group", observed=True)["umi_count_percent_10th"].median()
            
            ax3.bar(x_positions, medians_beta_top[order].values,
                color=[color_list[group] for group in order],
                width=0.45, edgecolor='none', alpha=0.7)

            #max_y = max([hollow_data_alpha["umi_count_percent_10th"].max(), filled_data_alpha["umi_count_percent_10th"].max()])
            #add_significance(ax3, 1, 2, max_y*1.1, 0.1, "**")
            #add_significance(ax3, 0, 2, max_y*1.25, 0.1, "*")


            plot_points(filled_data_beta, "umi_count_percent_10th", ax3, filled=True)
            plot_points(hollow_data_beta, "umi_count_percent_10th", ax3, filled=False)

        ax3.set_xlim(-0.5, 2.5)
            
        # Set y-axis limits for right columns (0-100) - same as column 2
        if i in [1]:
            axes[i, 2].set_ylim(0, 3)
            axes[i, 3].set_ylim(0, 3)
            axes[i, 2].set_yticks([0, 1, 2, 3])
            axes[i, 3].set_yticks([0, 1, 2, 3])
        elif i in [3]:
            axes[i, 2].set_ylim(0, 100)
            axes[i, 3].set_ylim(0, 100)
            axes[i, 2].set_yticks([0, 25, 50, 75, 100])
            axes[i, 3].set_yticks([0, 25, 50, 75, 100])
        elif i in [4]:
            axes[i, 2].set_ylim(0, 75)
            axes[i, 3].set_ylim(0, 75)
            axes[i, 2].set_yticks([0, 25, 50, 75])
            axes[i, 3].set_yticks([0, 25, 50, 75])    
        
        ax3.set_ylabel("")  # Empty y-label
        
        # Set cell type labels for leftmost column only
        # axes[i, 0].set_ylabel(titles_to_print[subset])
        
        # Hide y-ticks for columns 1 and 3 (right side of each group)
        axes[i, 1].tick_params(axis='y', labelleft=False)
        axes[i, 3].tick_params(axis='y', labelleft=False)
        
        # Set x-ticks only for bottom row
        if i == n_rows - 1:
            for j in range(4):
                target_ax = axes[i, j]

                # Check if this axis was processed by plot_broken_axes
                if hasattr(target_ax, 'bot_split'):
                    target_ax = target_ax.bot_split

                target_ax.set_xticks(x_positions)
                target_ax.set_xticklabels(["Y-Tx", "Y", "O"], rotation=45, weight="bold")
                target_ax.tick_params(axis='x', labelbottom=True)

        else:
            for j in range(4):
                axes[i, j].set_xticks([])
                axes[i, j].set_xticklabels([])
                axes[i, j].tick_params(axis='x', labelbottom=False)
    
    # Set column titles
    axes[0, 0].set_title('TCR alpha', pad=10)
    axes[0, 1].set_title('TCR beta', pad=10)
    axes[0, 2].set_title('TCR alpha', pad=10)
    axes[0, 3].set_title('TCR beta', pad=10)
    
    # Add left y-axis label (Richness per cell) on far left
    fig.text(0.02, 0.70, 'Richness per cell (or UMI)', 
            va='center', ha='center', rotation='vertical', fontsize=7)
    
    # Add right y-axis label (% pool size) in the spacer area between columns
    fig.text(0.02, 0.3, '% pool size of top 10 clones', 
            va='center', rotation='vertical', fontsize=7)
    
    # Add Panel Labels (A, B, C, D)
    # Adjust these x and y coordinates slightly based on your specific margins
    labels = [
        ('A', 0.02, 0.9), # Top Left (Richness Alpha/Beta)
        ('B', 0.48, 0.9), # Top Right (Placeholder/Other)
        ('C', 0.02, 0.48), # Bottom Left (Pool Size Alpha/Beta)
        ('D', 0.48, 0.48)  # Bottom Right
    ]

    for label, x, y in labels:
        fig.text(x, y, label, fontsize=12, fontweight='bold', va='top', ha='left')

    # Save figure
    fig.savefig(
        f"/home/erdem/Dropbox/Research/thymectomy/src/4.clone_per_cell_or_UMI/{plotName}_collapsed_cleanedupNaives_dev.pdf",
        dpi=300, bbox_inches='tight'
    )
    plt.close(fig)

# Plot both datasets
plot_overview(data_f, plotName="thymectomy_cellVSClone")
plot_overview(data_treg, plotName="Treg_cellVsClone")

