import glob
import os
import pandas as pd
import re

force = ""
threshold = 1

def alignment(r1, r2, vdjca_output, report_name):
    '''Step 1: Alignment'''
    command = ['$MIXCR', "align", '--preset', 'takara-human-tcr-V2-cdr3', '--rna', '-OsaveOriginalReads=true', '--report', report_name, r1, r2, vdjca_output, force]
    command_run = " ".join(command)
    print(f'STEP 1: alignment running with the command: \n {command_run} \n')
    os.system(command_run)
    return(print("Step 1 is done.\n\n\n"))

def refineTagsAndSort(vdjca, refined_vdjca, memory = '2950'):
    '''Step 2: UMI correction'''
    command = ['$MIXCR' , "refineTagsAndSort", "--memory-budget", memory, vdjca, refined_vdjca, force]
    command_run = " ".join(command)  
    print(f'STEP 2: UMI correction is running with the command: \n {command_run} \n')
    os.system(command_run)
    return(print("Step 2 is done.\n\n\n"))

def assemble(vdjca, clna, threshold = -1):
    '''Step 3: Assembly'''
    command = ['$MIXCR' , "assemble",  vdjca, clna, '--write-alignments', force]
    if threshold != -1:
        command = ['$MIXCR' , "assemble",  vdjca, clna, '--write-alignments', '--dont-split-clones-by V', '--dont-split-clones-by J', '--dont-split-clones-by C', force, "-P", f"assembler.minRecordsPerConsensus={str(threshold)}"] 
    else:
        command = ['$MIXCR' , "assemble",  vdjca, clna, '--write-alignments', '--dont-split-clones-by V', '--dont-split-clones-by J', '--dont-split-clones-by C', force]
    command_run = " ".join(command)
    print(f'STEP 3: Clone assembly is running with the command: \n {command_run} \n')
    os.system(command_run)
    return(print("Step 3 is done.\n\n\n"))

def exportClones(clns, result):
    '''Step 4: Export clones'''
    command = ['$MIXCR' , "exportClones", "-tagCounts", clns, result, force]
    command_run = " ".join(command)
    print(f'STEP 4: Clone export is running with the command: \n {command_run} \n')
    os.system(command_run)
    return(print("Step 4 is done.\n\n\n"))

def exportReads(clna):
    '''Step 5: Export clone reads'''
    command = ['$MIXCR' , "exportReadsForClones",  clna, fastq, '-s', force]
    command_run = " ".join(command) 
    os.system(command_run)
    print(f'STEP 5: Exporting original reads with the command: \n {command_run} \n')
    #os.system(command_run)
    return(print("Step 5 is done. Analysis ends here.\n\n\n\n\n"))
def make_dir(direct):
    #make individual folder
    try:
           os.makedirs(direct)
    except FileExistsError:
           # directory already exists
              pass
    os.chdir(direct)
    return()
    
fastq_dir_r1 = "/home/erdem/NOBINFBACKUP/thymectomy/results/all_preprocess/combine/"
fastq_dir_r2 = "/home/erdem/NOBINFBACKUP/thymectomy/results/all_preprocess/combine/r2_mixcr/"
results_dir = "/home/erdem/NOBINFBACKUP/thymectomy/results/mixcr_from_ip_tra_trb_re_run/"

#samples
samples = pd.read_pickle("/home/erdem/NOBINFBACKUP/thymectomy/results/pickles/samples_all.dat")

#quit()
print(samples.individual.unique())
for individual in samples.individual.unique()[::-1]:#labels.keys():
    individual_label =individual 
    ind_data = samples[samples['individual']==individual]
    #make_dir(results_dir + individual_label)
    for i, row in ind_data.iterrows():
        #print(row)
        cell = row.loc["subset"]
        #if (individual_label != "A10"):
        #    continue
        #if cell != "CD4NCD31-":
        #    print(cell)
        #    continue
        #if cell != "CD4NCD31-":
        #    print(cell)
        #    continue
        #print(individual_label, cell) 
        sample_number = row.loc["Sample number Genomescan"]
        sample_number = re.sub('-',  '', sample_number)
        #print(sample_number)
        #get A B fastq's
        r1_TRA = glob.glob(fastq_dir_r1 + individual_label + "-" + cell + "_TRA*" + "R1*")
        r2_TRA = glob.glob(fastq_dir_r2 + individual_label + "-" + cell + "_TRA*" + "R2*")
        r1_TRB = glob.glob(fastq_dir_r1 + individual_label + "-" + cell + "_TRB*" + "R1*")
        r2_TRB = glob.glob(fastq_dir_r2 + individual_label + "-" + cell + "_TRB*" + "R2*")
        #if ((len(r1_TRA) > 1) or (len(r1_TRA) > 1) or (len(r1_TRA) > 1)   or (len(r1_TRA) > 1)):
        #print(individual_label, cell) 
        #    print(r1_TRA)
         #   print(r2_TRA)
          #  print(r1_TRB)
           # print(r2_TRB)
        if not r1_TRA:
            print(f"{individual} {cell}\t not found.")
            continue
        r1_TRA = r1_TRA[0]
        r2_TRA = r2_TRA[0]
        r1_TRB = r1_TRB[0]
        r2_TRB = r2_TRB[0]

        #### make dir for TRA ####
        make_dir(results_dir + "/" + cell + "/" + individual_label + "/TRA") 
        print("[TRA] Analyzing... : ", r1_TRA, r2_TRA)
        run_label = "_".join([cell, individual_label, "TRA"])
        
        #step 1: alignment
        alignment(r1_TRA, r2_TRA, run_label+'.vdjca', run_label+'_report.txt')
        
        #step 2:umicorrection
        refineTagsAndSort(run_label+'.vdjca', run_label+'_refined.vdjca')
        

        #### make dir for thresholds ####
        make_dir(results_dir + "/" + cell + "/" + "threshold_" + str(threshold) + "_no_split" + "/" + individual_label + "/" + "/TRA") 
        print("[TRA] Analyzing... : ", r1_TRA, r2_TRA)
        vdjca_loc_TRA = results_dir + "/" + cell + "/" + individual_label + "/TRA/" 
        
        #step 3: assemble
        assemble(vdjca_loc_TRA + run_label+'_refined.vdjca', run_label+'_refined.clna', threshold = threshold)

        #step 4
        exportClones(run_label+'_refined.clna', run_label+'_result.tsv')
        print('[TRA] Done!/n')        



        ##########################
        #### make dir for TRB ####
        ##########################
        make_dir(results_dir + "/" + cell + "/" + individual_label + "/TRB") 
        print("[TRB] Analyzing... : ", r1_TRB, r2_TRB)
        run_label = "_".join([cell, individual_label, "TRB"])
        
        #step 1: alignment
        alignment(r1_TRB, r2_TRB, run_label+'.vdjca', run_label+'_report.txt')
        
        #step 2:umicorrection
        refineTagsAndSort(run_label+'.vdjca', run_label+'_refined.vdjca')

        #### make dir for thresholds ####
        make_dir(results_dir + "/" + cell + "/" + "threshold_" + str(threshold) + "_no_split"+ "/" + individual_label + "/" + "/TRB") 
        print("[TRB] Analyzing... : ", r1_TRB, r2_TRB)
        vdjca_loc_TRB = results_dir + "/" + cell + "/" + individual_label + "/TRB/" 
 
        #step 3: assemble
        assemble(vdjca_loc_TRB + run_label + '_refined.vdjca', run_label+'_refined.clna', threshold = threshold)

        #step 4
        exportClones(run_label+'_refined.clna', run_label+'_result.tsv')
        print("[TRB] Done!/n/n")          

