import glob
import pandas as pd
import numpy as np
import sys
import csv

csv.field_size_limit(sys.maxsize)

mixcr_results_dir = "../../../results/results/combinedGSIDs/*.1.clones.txt"

mixcr_results = glob.glob(mixcr_results_dir)
print(f"File paths to start with: {mixcr_results[:5]}")
cells = []; individuals =[]; chains = []; clonetype_count = []; UMIs = []; UMICountTotal = []; UMIs_str = []; readCount_list = []; readCountTotal = []; cell_count = []
v = []; d = []; j = []; c = []; nSeqCDR3 = []
done_samples = 0
dfs_sample = []

# load metadata
metadata = pd.read_csv("../../../data/metadata/metadata.csv")

for cell_dir in mixcr_results:
    gsid_cur = cell_dir.split("/")[-1].split(".")[0]
    if gsid_cur not in metadata["Sample number Genomescan"].unique():
        continue
    done_samples += 1
    print(f'{done_samples} is done out of {len(mixcr_results)}\r')
    
    meta_row = metadata[metadata["Sample number Genomescan"] == gsid_cur]
    
    cell = meta_row["subset"].item()
    if cell == "CD4CD31N":
        cell = "CD4NCD31"
    if cell == "CD4CD31-N":
        cell = "CD4NCD31-"
    cur_cell_count = meta_row["Cell Number"].item()
    individual = meta_row["individual"].item() 
    
    # read sample
    cur_res = pd.read_csv(cell_dir, sep = "\t", header = 0, engine='python')
    if len(cur_res) == 0:
        print(cell_ind_chain, " is empty!")
    
    # determine chain
    target_cols = ["allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore", "allCHitsWithScore"]
    mask_tra = cur_res[target_cols].map(lambda x: "TRA" in str(x)).any(axis=1)
    mask_trb = cur_res[target_cols].map(lambda x: "TRB" in str(x)).any(axis=1)
    cur_res["chain"] = np.select([mask_tra, mask_trb], ["TRA", "TRB"], default="other")
    cur_res_tra = cur_res[mask_tra]
    cur_res_trb = cur_res[mask_trb]

    cur_res_tra = cur_res_tra.groupby("nSeqCDR3", as_index=False).agg("sum")
    cur_res_trb = cur_res_trb.groupby("nSeqCDR3", as_index=False).agg("sum")

    # sort by abundance
    cur_res_tra = cur_res_tra.sort_values("uniqueMoleculeCount", ascending = False)
    cur_res_trb = cur_res_trb.sort_values("uniqueMoleculeCount", ascending = False)

    # populate lists with data
    ## Read data
    readCount_list.append(cur_res_tra["readCount"].tolist())
    readCount_list.append(cur_res_trb["readCount"].tolist())
    readCountTotal.append(sum(cur_res_tra["readCount"].tolist()))
    readCountTotal.append(sum(cur_res_trb["readCount"].tolist()))
    ## nucleotide sequence CDR3
    nSeqCDR3.append(cur_res_tra["nSeqCDR3"].tolist())
    nSeqCDR3.append(cur_res_trb["nSeqCDR3"].tolist())
    # Subset&Individual
    cells.append(cell)
    cells.append(cell)
    individuals.append(individual)
    individuals.append(individual)
    cell_count.append(cur_cell_count)
    cell_count.append(cur_cell_count)
    
    ## UMI data
    UMIs.append(cur_res_tra["uniqueMoleculeCount"].tolist())
    UMIs.append(cur_res_trb["uniqueMoleculeCount"].tolist())
    UMICountTotal.append(sum(cur_res_tra["uniqueMoleculeCount"].tolist()))
    UMICountTotal.append(sum(cur_res_trb["uniqueMoleculeCount"].tolist()))
    ## chains
    chains.append("TRA")
    chains.append("TRB")
    #add VDJC
    v.append(cur_res_tra["allVHitsWithScore"].tolist())
    v.append(cur_res_trb["allVHitsWithScore"].tolist())
    d.append(cur_res_tra["allDHitsWithScore"].tolist())
    d.append(cur_res_trb["allDHitsWithScore"].tolist())
    j.append(cur_res_tra["allJHitsWithScore"].tolist())
    j.append(cur_res_trb["allJHitsWithScore"].tolist())
    c.append(cur_res_tra["allCHitsWithScore"].tolist())
    c.append(cur_res_trb["allCHitsWithScore"].tolist())
    clonetype_count.append(len(cur_res_tra))
    clonetype_count.append(len(cur_res_trb))

#combine the lists of data into a dataframe
merged = pd.DataFrame({"individual": individuals, "subset": cells, "chain": chains, "nSeqCDR3": nSeqCDR3,
                        "clonetype": clonetype_count,
                        "UMI count": UMICountTotal, 
                        "umi_counts": UMIs,
                        "Cell Number": cell_count,
                        "readCount_list": readCount_list, "Read count": readCountTotal,
                        "v":v, "d":d, "j":j, "c":c})
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
merged.to_pickle("../../../results/mergedGSIDs_mixcr_th1_collapsed_clones.pkl")
