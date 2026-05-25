import pandas as pd
import glob
import pickle

#new labels of individuals
labels = {'HD1': 'Y1', 'HD11': 'Y5', 'HD2': 'Y2', 'HD6': 'Y3', 'HD8': 'Y4', 'HD15': "Y6", "HD18": "Y7", '609330500': 'A1', '610570345': 'A2', '617350339': 'A3', '620750446': 'A4', 'UU030509': 'A5', 'UU030408': 'A6', '0507':"A7", '0384': "A8", 'UU030379': "A9", 'UU030314': "A10", '446': "A11", "UU030316": "A12", "UU030393": "A13", 'RTHYM5':"T2", 'RTHYM10': "T4" ,  'RTHYM17':"T8", 'RTHYM12':"T5", 'RTHYM20': "T10", 'RTHYM14': "T6", 'RTHYM16': "T7", 'RTHYM9': "T3", 'RTHYM1':"T1", 'RTHYM18': "T9"}


samples_all = pd.read_excel("../data/metadata/samples_all.xlsx", usecols = "B,C,E,F,H:K", sheet_name = "Overview donors")

samples_all.drop([13,40,41,82,83,98], inplace = True) #remove empty rows depending on the excel format!
samples_all = samples_all.fillna(method="ffill")
samples_all.Donor = samples_all.Donor.str.replace("-","")
samples_all["Population"] = samples_all["Population"].str.replace(" ", "").str.replace("+","",regex=False).str.replace("naive", "N")
samples_all["Donor"] = samples_all["Donor"].str.strip()
samples_all["individual"] = samples_all["Donor"].map(labels)
samples_all.rename(columns = {"Population": "subset", "Cell number": "Cell Number"}, inplace = True)
samples_all.reset_index(inplace = True, drop = True)

#write metadata to a file
meta_dir = "../data/metadata/"
#with open(pickles_dir+'samples_all.dat', 'wb') as handle:
#    pickle.dump(samples_all, handle, protocol=pickle.HIGHEST_PROTOCOL)
#print('Sample data saved.')
samples_all.to_csv(meta_dir+"metadata.csv")

#fastq_dir = "/home/erdem/Hosts/vacuole1/thymectomy/data/fastq/"
fastq_dir = "/home/erdem/Hosts/vacuole1/thymectomy/data/fastq_combined_GSIDs/"

#gsids for all
sampleID = []
read1_file_name = []
read2_file_name = []
for i,sample in samples_all.iterrows():
    individual = sample["Donor"]
    individual_label = sample["individual"]
    cell_pop = sample["subset"]
    sample_number = sample["Sample number Genomescan"]
    if not glob.glob(fastq_dir+"*"+sample_number+"*R1*"):
        print(f"{individual} {sample_number} \t is not found")
        continue
    #print(f"{individual} \t is not found")
    r1 = glob.glob(fastq_dir+"*"+sample_number+"*R1*")[0].split("/")[-1]
    r2 = glob.glob(fastq_dir+"*"+sample_number+"*R2*")[0].split("/")[-1]
    sample_id = "-".join([individual_label,cell_pop])
    #add to lists
    sampleID.append(sample_id)
    read1_file_name.append(r1)
    read2_file_name.append(r2)
    #combine the file and write
metadata_gsid = {"sampleID": sampleID, "read1_file_name": read1_file_name, "read2_file_name": read2_file_name}
metadata_gsid = pd.DataFrame(metadata_gsid)
metadata_gsid.to_csv(meta_dir+'metadata_gsid.csv', index = False)
