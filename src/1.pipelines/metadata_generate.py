
labels = {'HD1': 'Y1', 'HD11': 'Y5', 'HD2': 'Y2', 'HD6': 'Y3', 'HD8': 'Y4',
                      '609330500': 'A1', '610570345': 'A2', '617350339': 'A3', '620750446': 'A4',
                                    '63001509': 'A5', '63074408': 'A6'}
# print metadata file
sampleID = []
read1_file_name = []
read2_file_name = []
for i,sample in samples.iterrows():
    individual = sample["Donor"]
    cell_pop = sample["Cell population"]
    sample_number = sample["Sample number Genomescan"]
    if not glob.glob(fastq_dir+"*"+sample_number+"*R1*"):
        continue
    r1 = glob.glob(fastq_dir+"*"+sample_number+"*R1*")[0].split("/")[-1]
    r2 = glob.glob(fastq_dir+"*"+sample_number+"*R2*")[0].split("/")[-1]
    cell_pop = cell_pop.replace(" ", "").replace("+","").replace("naive", "N")
    individual = individual.replace("-","")
    individual = labels[individual]
    sample_id = "-".join([individual,cell_pop])
    #gather
    sampleID.append(sample_id)
    read1_file_name.append(r1)
    read2_file_name.append(r2)
    combine the file and write
    metadata = {"sampleID": sampleID, "read1_file_name": read1_file_name, "read2_file_name": read2_file_name}
    metadata = pd.DataFrame(metadata)
    metadata.to_csv(meta_dir+'metadata.csv', index = False)
