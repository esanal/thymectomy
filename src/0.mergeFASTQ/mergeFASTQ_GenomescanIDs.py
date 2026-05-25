import pandas as pd
import os
import ast
import subprocess
from snakemake.io import glob_wildcards
from pathlib import Path

FASTQ_DIR = Path("/home/erdem/Hosts/vacuole1/thymectomy/data/fastq")
output_dir = "/home/erdem/Hosts/vacuole1/thymectomy/data/fastq_combined_GSIDs"

SAMPLES = [
    f.name.split('_')[1]
    for f in FASTQ_DIR.glob("*_R1.fastq.gz")
]

# Index the directory once for speed
all_files = os.listdir(FASTQ_DIR)

processed = []
# Iterate through rows
for Sid in SAMPLES:
    if Sid not in processed:
        print(f"Processing {Sid}")
        for direction in ['R1', 'R2']:
            matching_files = []
            # Find all files matching the IDs for this direction
            #for sample_id in []:
            matches = [f for f in FASTQ_DIR.glob(f"*{Sid}*{direction}*.fastq.gz")]
        
            # Deduplicate matches
            matching_files = sorted([str(a) for a in list(set(matches))])
             
            if not matching_files:
                print(f"  [!] No {direction} files found for IDs: {ids}")
                continue

            output_file = os.path.join(output_dir, f"{Sid}_{direction}.fastq.gz")

            # --- LOGIC: Concatenate or Copy ---
            if len(matching_files) > 1:
                # Multiple IDs: Concatenate
                inputs = " ".join(f"'{f}'" for f in matching_files) # Quoting paths for safety
                print(f"  [Merge]\n {"\n".join(matching_files)} files\n ->\n {output_file}")
                subprocess.run(f"cat {inputs} > '{output_file}'", shell=True, check=True)
        
            else:
                # Only one ID: Copy
                source_file = matching_files[0]
                print(f"  [Copy] Single file -> {output_file}")
                subprocess.run(f"cp '{source_file}' '{output_file}'", shell=True, check=True)
        processed.append(Sid) 

print("\nProcess complete. Merged files are in:", output_dir)
