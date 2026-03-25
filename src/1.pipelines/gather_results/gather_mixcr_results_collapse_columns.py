import glob
import pandas as pd
import numpy as np
import sys
import csv

csv.field_size_limit(sys.maxsize)

mixcr_results_dir = "/home/erdem/NOBINFBACKUP/thymectomy/results/mixcr_from_ip_tra_trb_re_run/*/threshold_1/*/*/*_result_*.tsv"
mixcr_results = glob.glob(mixcr_results_dir)
print(f"File paths to start with: {mixcr_results[:5]}")

cells = []; individuals =[]; chains = []; clonetype_count = []; UMIs = []; UMICountTotal = []; UMIs_str = []; readCount_list = []; readCountTotal = []
v = []; d = []; j = []; c = []; nSeqCDR3 = []
done_cells = 0

for cell_dir in mixcr_results:
    # cell_dir = mixcr_results[0]
    print(f'{done_cells} is done out of {len(mixcr_results)}\r')
    done_cells += 1
    cell_cur = cell_dir.split("/")[11]
    cell_ind_chain = cell_cur.split("_")
    cell = cell_ind_chain[0]
    if cell == "CD4CD31N":
        cell = "CD4NCD31"
    if cell == "CD4CD31-N":
        cell = "CD4NCD31-"
    individual = cell_ind_chain[1]
    chain_input = cell_ind_chain[2]
    chain = cell_ind_chain[4]
    if chain_input != chain[:3]:
        continue
    
    #print current sample
    print(cell_ind_chain)
    # read sample
    cur_res = pd.read_csv(cell_dir, sep = "\t", header = 0, engine='python')
    #check if empty
    #subset_df = mixcr_results_df.loc[(mixcr_results_df.individual == individual) & (mixcr_results_df.subset == cell) & (mixcr_results_df.chain == chain_input)]
    if len(cur_res) == 0:
        print(cell_ind_chain, " is empty!")

    # Fix UMI count data by converting format to list
    umi_as_list = []
    umi_tags = cur_res["tagCounts"].tolist()
    for cloneID, umi_list in enumerate(umi_tags):
        # umi_list = cur_res["tagCounts"][0]
        # umis = cur_res["tagCounts"][8000]
        if ":" not in umi_list:
            umis = umi_list[1:-1].split(",")
            umis = [[umi.split("=")[0], int(float(umi.split("=")[1]))] for umi in umis]
        else:
            umis = [[umi_list.split(": ")[0], int(float(umi_list.split(": ")[1]))]]
        umi_as_list.append(umis)
    cur_res["tagCounts"] = umi_as_list


    # Group by clonetypes
    cur_res = cur_res.groupby("nSeqCDR3", as_index=False).agg("sum")

    # sort by abundance
    cur_res = cur_res.sort_values('uniqueUMICount', ascending = False)
    
    #Grouped umis re-calculated
    umi_as_list_clones_collapsed = cur_res["tagCounts"].tolist()

    # populate lists with data
    ## Read data
    readCount_list.append(cur_res["readCount"].tolist())
    readCountTotal.append(sum(cur_res["readCount"].tolist()))
    ## nucleotide sequence CDR3
    nSeqCDR3.append(cur_res["nSeqCDR3"].tolist())
    # Subset&Individual
    cells.append(cell)
    individuals.append(individual)
    ## UMI data
    UMIs.append(cur_res["uniqueUMICount"].tolist())
    UMICountTotal.append(sum(cur_res["uniqueUMICount"].tolist()))
    
    UMIs_str.append(umi_as_list_clones_collapsed)
    chains.append(chain_input)
    #add VDJC
    v.append(cur_res["allVHitsWithScore"].tolist())
    d.append(cur_res["allDHitsWithScore"].tolist())
    j.append(cur_res["allJHitsWithScore"].tolist())
    c.append(cur_res["allCHitsWithScore"].tolist())
    clonetype_count.append(len(cur_res))

#combine the lists of data into a dataframe
data_f = pd.DataFrame({"individual": individuals, "subset": cells, "chain": chains, "nSeqCDR3": nSeqCDR3,
                        "clonetype": clonetype_count,
                        "umi_counts": UMIs, "UMI count": UMICountTotal, "umi_tag_read": UMIs_str,
                        "readCount_list": readCount_list, "Read count": readCountTotal,
                        "v":v, "d":d, "j":j, "c":c})

#Get number of cells
samples_all = pd.read_pickle("/home/erdem/NOBINFBACKUP/thymectomy/results/pickles/samples_all.dat")

#combine cell number data
merged = data_f.merge(samples_all[["individual", "subset", "Cell Number"]], how = "left", on = ["individual", "subset"])
#check if any cell numbers left empty...
print(merged[merged["Cell Number"].isna()])

# Correct group names
merged.loc[merged.individual.str.startswith("T"), "Group"] = "Thymectomized"
merged.loc[merged.individual.str.startswith("A"), "Group"] = "Aged"
merged.loc[merged.individual.str.startswith("Y"), "Group"] = "Young"
merged = merged.astype({"individual": "category", "subset": "category", "chain": "category", "Group": "category"})
print(merged.info())
#save file
#merged.to_hdf('mixcr_th1_raw.h5', key='data', mode='w', format = "table")
merged.to_pickle("/home/erdem/NOBINFBACKUP/thymectomy/results/pickles/mixcr_th1_raw_collapsed_clones.pkl")
