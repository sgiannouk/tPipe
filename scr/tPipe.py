#usr/bin/python
###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.2.0"

from multiprocessing import Pool
import multiprocessing as mp
import ConfigParser, argparse
from datetime import datetime
import os, subprocess
import commands, shutil
import numpy as np
import xlsxwriter, sys, xlrd
import natsort as ns
from operator import itemgetter

libs = ["m","p"]
reftLib = {}
secstructLib = {}
nrm_em = {}

# Lists and other requirements. 
Dflex = []
Dfl1 = []
Dfl2 = []
DD5ex = []
DD51 = []
DD52 = []
DD3ex = []
DD31 = []
DD32 = []
DD5Hex = []
DD5H1 = []
DD5H2 = []
DD3Hex = []
DD3H1 = []
DD3H2 = []
DMidex = []
DMid1 = []
DMid2 = []

dics0 = [Dflex, DD5Hex, DD3Hex, DD5ex, DD3ex, DMidex]
dics1 = [Dfl1, DD5H1, DD3H1, DD51, DD31, DMid1]
dics2 = [Dfl2, DD5H2, DD3H2, DD52, DD32, DMid2]
group_dics = [dics0, dics1 ,dics2]

start_time = datetime.now()

# Parsing input arguments
def parse_commandline():

    global outputDir
    global aligned_folder
    global processed_data_folder
    global reports_folder
    global log_folder
    global qc_folder
    global results_folder
    global implementation
    global thrds



    usage = "tpipe <configuration_file.txt>"
    
    description = "DESCRIPTION\
                   \n-----------\
                \ntPipe was designed to preprocess and analyse small RNA-Seq datasets.\
                \n\n\tThe analysis contains the following steps:\
                \n\t1. Quality Control\
                \n\t2. Preprocessing\
                \n\t3. Aligning reads against the reference tLibraries\
                \n\t4. Classification\
                \n\t5. Expression Matrix\
                \n\t6. Basic Statistics"
    
    epilog = " -- September 2016 | Stavros Giannoukakos -- "
    
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, version=__file__.split(".")[0]+" version "+__version__, epilog=epilog)

    parser.add_argument ("configuration_file", type=str, nargs=1, help="\nInput configuration file to parse and load necessary info and options")
    
    # Number of threads to be used
    parser.add_argument("--threads", dest="threads", default=1, 
                        help="\nNumber of threads to be used in the analysis")
    
    # Adapter sequence
    parser.add_argument("--adapter", dest="adapter",
                        help="\nSequence of the 3' end adapter. (default: Illumina Small RNA 3' Adapter: \"TGGAATTCTCGGGTGCCAAGGAACTC\")")
    
    # Skip Quality Control
    parser.add_argument("-no_QC", dest="no_QualityControl", action="store_true", 
                        help="\nIf this option is activated, the Quality Control step will be skipped.\
                        \nThe results will contain no QC report.\nATTENTION: By activating it, the time of the analysis will decrease!")

    # Skip preprocessing
    parser.add_argument("-no_Preprocessing", dest="no_Preprocessing", action='store_true', 
                        help="\nIf this option is activated, the preprocessing step will be skipped.\
                        \nATTENTION: By activating it, the time of the analysis will decrease!")
 

    # Skip Expression Matrix
    parser.add_argument("-no_stats", dest="no_stats", action='store_true', 
                        help="\nIf this option is activated, statistics analysis will be skipped.\
                        \nATTENTION: By activating it, the time of the analysis will decrease!")    

    # input folder option
    parser.add_argument("-tLibs", dest="libraries", nargs=1, help="\nPath of the input directory that contains the Bowtie1 indexes of the tLibraries.")

    # input folder option
    parser.add_argument("-i", dest="input_dir", nargs=1, help="\nPath of the input directory that contains the data.")
    
    # output folder options
    parser.add_argument("-o", dest="output_dir", nargs=1, help="\nPath of the output directory that analysis will be stored. (default: current directory)")


    # Get the options and return them
    args = parser.parse_args()

    # Reading configuration file
    config = ConfigParser.ConfigParser()
    config.read(args.configuration_file)


    # Output folders
    folders = ["Aligned_files", "Processed_data", "Reports", "Reports/Log_reports", "Reports/QCreports", "Results"]

    # Error control and file creation from conf.file
    for section_name in config.sections():
        for name, value in config.items(section_name):
            if not value:
                print "No input value in %s: %s" %(section_name,name)
            else:
                value = value.strip("\"")
                if name == "threads":
                    thrds = value
                elif name == "output_folder":
                    try:
                        if not (value == "-"):
                            if not os.path.exists(value):
                                os.makedirs(value)
                                outputDir = value
                            else:
                                outputDir = value
                        else:
                            outputDir = os.path.join(os.getcwd(),"tPipe_analysis")
                        for i, fol in enumerate(folders):
                            if not os.path.exists(os.path.join(outputDir,fol)):
                                os.makedirs(os.path.join(outputDir,fol))
                    except Exception as e: print "FATAL ERROR - ", e  
    
    # Creating the above folders
    for fol in folders:
        if not os.path.exists(os.path.join(outputDir, fol)):
            os.makedirs(os.path.join(outputDir, fol))
    
    #Define global output folder
    aligned_folder = os.path.join(outputDir,"Aligned_files")
    processed_data_folder = os.path.join(outputDir,"Processed_data")
    reports_folder = os.path.join(outputDir,"Reports")
    log_folder = os.path.join(outputDir,"Reports/Log_reports")
    qc_folder = os.path.join(outputDir,"Reports/QCreports")
    results_folder = os.path.join(outputDir,"Results")


    (infiles, qcfiles) = check_files("".join(config.get("IO", "input_folder")))

    # Output log file
    sys.stdout=open(os.path.join(reports_folder,"tPipeLog.log"),"w", 0)
    print str("tpipe is about to begin - " + start_time.strftime("%d.%m.%Y %H:%M:%S"))
    print "\n\nConfiguration File - Input Options".upper()
    print "------------------------------------------------------------"
    for each_section in config.sections():
        print "\n"+each_section.upper()
        for (each_key, each_val) in config.items(each_section):
            print each_key + ": "+each_val
    print "\n------------------------------------------------------------"
    print "------------------------------------------------------------\n\n"
    return (config, infiles, qcfiles)

# Checking input arguments
def check_files(input_folder):

    input_files = []

    for paths, subdirs, files in os.walk(input_folder):
        for name in files:
            if name.startswith("."):
                pass
            elif not name.endswith((".fastq",".fastq.gz",".fq",".fq.gz")):
                print "Unidentified input format in file: %s\nCheck help option" %fileR2.split("/")[-1]
            else:
                input_files.append(os.path.join(paths, name))

    # Converting the list of files into string - needed for fastQC input 
    qcfiles = " ".join(input_files)

    return (input_files, qcfiles)


############# QC #############

# Quality Control of the raw reads
def quality_control(fastqc_files):

    run_qc = " ".join(["fastqc", "--threads", thrds, "-q", "-o", qc_folder, fastqc_files])
    subprocess.call(run_qc, shell=True) 

    return


############# PREPROCESSING #############

# Preprocessing of the data
def preprocessing(args):

    nameoffile = args[0].split("/")[-1].split(".")[0]
    print "Preprocessing of: %s" %(args[0].split("/")[-1].split(".")[0])
    # Removing 3' adapter
    run_adaptrm = " ".join(["cutadapt", "-a", args[1], "-m 15", "-q 15","--trim-n", 
        "--too-short-output", os.path.join(processed_data_folder, nameoffile+".short.fastq"),
        "-o", os.path.join(processed_data_folder, nameoffile+".tr.fastq"), args[0]]) #,">",os.path.join(log_folder, nameoffile+"_cutadapt.log" ,
    with open(os.path.join(log_folder,nameoffile+"_cutadapt.log"),"w") as out: subprocess.call(run_adaptrm,shell=True,stdout=out) 

    # Cleaning
    for path, subdirs, files in os.walk(processed_data_folder):
        for name in files:
            # Rename according to the size-depth of the file and collapsing
            if name.endswith(".tr.fastq") and nameoffile in name:
                ffiles = os.path.join(path, name)
                
                rt = int(commands.getoutput('grep + -c '+ ffiles))
                result = chformat(rt)
                collapsing_fastq = " ".join(["fastx_collapser", "-i", ffiles, "-o", ffiles.split(".")[0]+"#"+str(result)+".fasta"])
                subprocess.call(collapsing_fastq, shell=True)
                os.remove(ffiles)
            if args[2] == True:
                # Compressing short reads
                if name.endswith(".short.fastq") and nameoffile in name:
                    compr = " ".join(["gzip -1", os.path.join(path,name)])
                    subprocess.call(compr, shell=True)     
  
    return

# Detecting the depth of the input file.
def chformat(num):

    mg = 0

    while abs(num) >= 1000:
        mg += 1
        num /= 1000.0

    return '%.2f%s' % (num, ['', 'K', 'M', 'G', 'T'][mg])


############# MAPPING #############

# Aligning reads against reference libraries and filtering
def align(inputf, ind, p_ind, mismatch, refgenome, tlibs):

    index = ind.split("/")[-1] # tLibrary 
    prev_ind = p_ind.split("/")[-1] # Previous tLibrary
    
    frmat = ".fasta" if inputf is processed_data_folder else "_" + prev_ind + ".un.reads.fa"
    for path, subdirs, files in os.walk(inputf):
        for name in files:
            if name.endswith(frmat):
                
                data = os.path.join(path, name)
                
                naming = name[:-6] if inputf is processed_data_folder else name[:-16]
                file_name = naming+"_"+str(mismatch)+"_"+index
                    
                # Creating output folder
                rfolder = os.path.join(aligned_folder, file_name.split("#")[0])
                if not os.path.exists(rfolder): os.makedirs(rfolder) 
                
                # Naming output files.  
                unaligned = os.path.join(rfolder, file_name + ".un.reads.fa")
                aligned = os.path.join(rfolder, file_name + ".al.fa")
                btresult = os.path.join(rfolder, file_name + ".out")
                
                # Running bowtie1 / printing stats 
                # print "Mapping: %s with %s mismatches" %(data.split("/")[-1].split(".")[0], str(mismatch))
                run = " ".join(["bowtie", "--threads", thrds, "-v", str(mismatch) , "-a", "--best", "--strata",
                    "--norc", ind, "--un", unaligned,"--al", aligned, "-f", data, btresult])
                subprocess.call(run, shell=True)
                if inputf is aligned_folder and not mismatch == "2": os.remove(data)
                print "\n"

    
    
    ### Filtering Step ###
    # Rerunning reads, found with 1 and 2 that are potentially tRNAs, against the reference genome (hg19) this time.
    if index == "m" and mismatch == 2:
        # Remapping 1mm detected reads against ref genome allowing no mm
        ## Remapping 2mm detected reads against ref genome allowing up to 1 mm
        for i in range(1,3):
            for path, subdirs, files in os.walk(aligned_folder):
                for name in files:
                    if name.endswith(str(i)+"_m.al.fa"):
                        # Obtaining data to realign
                        data = os.path.join(path, name)
                        aligned = os.path.join(path, name[:-10]+str(i)+"control.al.fa")
                        # Run aligned reads against the reference genome - need at least one hit for the control
                        filter_aligned = " ".join(["bowtie", "--quiet", "--threads", thrds, "-v 0" , "-k 1", "--best", refgenome, "--al", aligned, "-f", data])
                        subprocess.call(filter_aligned, shell=True)
                        print "Applying filtering step in the %s sample\n" %name[:-10]
                        filter_out(name[:-10], str(i))


        # Cleaning the data
        for path, subdirs, files in os.walk(aligned_folder):
            for name in files:
                if name.endswith(("_m.al.fa", "_p.al.fa", "_1_m.out", "_2_m.out")):
                    os.remove(os.path.join(path, name))
                elif name.endswith(".un.reads.fa"):
                    file = os.path.join(path, name)
                    os.rename(file, file.split("#")[0]+".unassigned.fasta")
                elif name.endswith("control.al.fa"):
                    ctr = os.path.join(path, name)
                    with open(os.path.join(path, name.split("#")[0]+".discarded.fasta"), 'a+') as disc_out, open(ctr,"r") as ctr_in:
                        for line in ctr_in: disc_out.write(line)
                        os.remove(ctr)
                elif name.endswith(".out") and "_p." in name or "_m." in name:
                    exact = os.path.join(path, name)
                    with open(os.path.join(path, name[:-8]+"_0.out"), 'a+') as ex_out, open(exact,"r") as ex_in:
                        for line in ex_in: ex_out.write(line)
                        os.remove(exact) 

    return

# Removing reads that were mapped against the genome with less mismatches. 
def filter_out(filename, mm):

    seqs = {}
    btoutput = []
    remove_reads = []
    
    # Read and save reads that did not pass the control.
    for path, subdirs, files in os.walk(aligned_folder):
        for name in files:
            if filename in name and  name.endswith(str(mm)+"control.al.fa"):
                to_remove = os.path.join(path, name)
                with open(to_remove) as fin:
                    for line in fin:
                        if line.startswith(">"): remove_reads.append(line.strip()[1:])

    for path, subdirs, files in os.walk(aligned_folder):
        for name in files:
            # Editing bowtie output files
            if filename in name and  name.endswith(str(mm)+"_m.out"):
                bt_reads = os.path.join(path, name)
                with open(bt_reads, "r") as oin, open(bt_reads.replace("_m", ""), "w+") as oout:
                    for line in oin:
                        if not line.split("\t")[0].strip() in remove_reads: oout.write(line)
    
    return 


############# CLASSIFICATION #############

# Classification of the mapped data
def main_analyis(keep_ss, reference_tLib, reference_secstr):

    global book

    ## Extensive EM sheets.
    # Output file 1: Expression Matrix
    book = xlsxwriter.Workbook(os.path.join(results_folder, "Expression_Matrix(extensive).xlsx"))
    headerf = book.add_format({'align' : 'center', 'valign' : 'vcenter'})
    
    fl0 = book.add_worksheet("FL|exact")
    fl1 = book.add_worksheet("FL|1")
    fl2 = book.add_worksheet("FL|2")

    DH50 = book.add_worksheet("5h|exact")
    DH51 = book.add_worksheet("5h|1")
    DH52 = book.add_worksheet("5h|2")

    DH30 = book.add_worksheet("3h|exact")
    DH31 = book.add_worksheet("3h|1")
    DH32 = book.add_worksheet("3h|2")

    D50 = book.add_worksheet("5p|exact")
    D51 = book.add_worksheet("5p|1")
    D52 = book.add_worksheet("5p|2")

    D30 = book.add_worksheet("3p|exact")
    D31 = book.add_worksheet("3p|1")
    D32 = book.add_worksheet("3p|2")

    D0mid = book.add_worksheet("D2T|exact")
    D1mid = book.add_worksheet("D2T|1")
    D2mid = book.add_worksheet("D2T|2")

    sheets0 = [fl0,D50,D30,DH50,DH30,D0mid]
    sheets1 = [fl1,D51,D31,DH51,DH31,D1mid]
    sheets2 = [fl2,D52,D32,DH52,DH32,D2mid]

    sheets = sheets0+sheets1+sheets2
    group_sh = [sheets0, sheets1, sheets2]

    # Creating each sheet with the appropriate annotation etc.
    for sh in sheets: sh.write(0,0,"tRNA_Annotation", headerf), sh.write(0,1, "Sequence", headerf)

    # Creating a list with the reference tLibraries.
    with open(reference_tLib) as refin:
        for line in refin:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                seq = line.strip().upper()
                reftLib[name] = name,seq

    # Creating a dictionary with the reference secondary str. library.
    with open(reference_secstr) as fin:
        for lines in fin:
            # Reading secondary structure reference library.
            if lines.startswith(">"):
                name = lines[1:].split(" ")[0].strip()
            elif lines.startswith("SS_info"):
                ss_info = lines.split(" ")[1].strip()
            elif lines.startswith("Seq"):
                seq = lines.split(" ")[1].strip()
            elif lines.startswith("Str"):
                str_info = lines.split(" ")[1].strip()
                secstructLib[name]=(ss_info,str_info,seq) 
    
    
    # Initialise variables
    num = []

    number = 0
    count = 0
    # Obtaining the number of samples for analysis.
    for path, subdirs, files in os.walk(aligned_folder):
        for name in files:
            if name.endswith("0.out"):
                number += 1
            
            newname = name[:-6]
            f = os.path.join(path, name)
            
            # Appending the input bowtie files and several info.
            if name.endswith("0.out") and newname in f:
                num.append((f, count, newname, 0))
            elif name.endswith("1.out") and newname in f:
                num.append((f, count, newname, 1))
            elif name.endswith("2.out") and newname in f:
                num.append((f, count, newname, 2))
        count += 1

    # Classifying the input data.
    for r_val in num:
        classify(r_val[0], r_val[1], r_val[2], r_val[3],keep_ss)
        for sh in sheets:
            sh.write(0, (r_val[1]+1), r_val[2],headerf)

    # final_data = []
    # Creating Extensive Expression Matrix 
    for sh_idx, sh_sheets in enumerate(group_sh):
        for sh_idx1, sh_sheets1 in enumerate(sh_sheets):
            for d_idx, dval in enumerate(group_dics):
                for d_idx1, dval1 in enumerate(dval):
                    if sh_idx == d_idx and sh_idx1 == d_idx1:
                        write_excel((dval1, number, sh_sheets1))
                        # final_data.append((dval1, number, sh_idx, sh_idx1))


    # Multiprocessing
    # pool = mp.Pool(processes=int(thrds))
    # results = pool.map(create_excel, final_data)
    # pool.close()

    book.close()

    return 

# Main classification function
def classify(input_file, indx, name, mismatch, keep_ss):


    # Depending on the mismatch, obtaining the proper dictionary-list.           
    if mismatch == int(0):D = dics0
    if mismatch == int(1):D = dics1
    if mismatch == int(2):D = dics2

    sample = indx+1
    # Reading the input bowtie data.
    data = reading_bowtie_file(input_file)

    # Analysing the inout data and collect the appropriate info for the classification. 
    for d_data in data:
        if d_data[0] in reftLib:

            # Obtaining the ss info for each read. 
            lref = len(reftLib[d_data[0]][1])
            llib = len(d_data[1])
            sa, sb = secondary_structure_info(d_data[0].strip())

            # Example: In the secondary structure file we obtain for each read this info:
            #   >chrX-X-X.trna6-7-5-IleGAT.pre
            #   SS_info: D:10-26;A:28-44;T:50-66
            #   Seq: GGCCGGTTAGCTCAGTTGGTAAGAGCGTGGTGCTGATAACACCAAGGTCGCGGGCTCGACTCCCGCACCGGCCA
            #   Str: >>>>>>>..>>>>.........<<<<.>>>>>.......<<<<<.....>>>>>.......<<<<<<<<<<<<.

            # Were SS_info is the start and end of the D-loop, the Anticodon-loop and the T-loop
            # and Str in the secondary structure itself. 

            # Here by in "sa" variable we obtain all the secondary structure info (SS_info) (so, here: sa = D:10-26;A:28-44;T:50-66)
            # and by sb, we obtain the secondary structure itself (so, sb = >>>>>>>..>>>>.........<<<<.>>>>>.......<<<<<.....>>>>>.......<<<<<<<<<<<<.)
            

            ending_point = int(d_data[3]+llib) # Where: d_data[3] = (starting point), llib = (length of the tRNA of tRF that we found=llib)
            Dl_end = int(sa.split("-")[1].split(";")[0].strip())    # End of the D-loop 
            Al_start = int(sa.split("-")[1].split(":")[1].strip())  # Beginning of Anticodon-loop 
            Al_end = int(sa.split("-")[2].split(";")[0].strip())    # End of Anitcodon-loop
            Tl_start = int(sa.split("-")[2].split(":")[1].strip())  # Beginning of T-loop 


          

            ### Classification Algorithm ###
            # FULL length
            # Here the "d_data[3]" (this is directly obtained from bowtie output) is the starting point of the tRNA or tRF. If the starting point is 
            # 0 and also the ending point of the tRNA is equal to the reference tRNA length, then we have full length molecule 
            if (int(d_data[3]) == int(0) and ending_point == lref):
                nn0 = str(d_data[0]+"_"+str(d_data[3]+1)+":"+str(d_data[3]+llib)+"-FL")
                seq2 = final_seq(d_data[1], d_data[3], d_data[5], lref, llib, d_data[4], name, nn0, sb, d_data[2],keep_ss)
                D[0].append((nn0, seq2, d_data[2], sample, d_data[5], d_data[4], name.split("#")[0]))
                
            # 5 - HALFs (5p - somewhere in A-Loop)
            # Here the starting point should be 0 - so starting from the beginning and the ending point shoule be between the boundaries of the Anticodon-loop 
            elif (int(d_data[3]) == int(0) and  (Al_start <= ending_point <= Al_end)):
                nn0 = str(d_data[0]+"_"+str(d_data[3]+1)+":"+str(d_data[3]+llib)+"-5H")
                seq2 = final_seq(d_data[1], d_data[3], d_data[5], lref, llib, d_data[4], name, nn0, sb, d_data[2],keep_ss)
                D[3].append((nn0, seq2, d_data[2], sample, d_data[5], d_data[4], name.split("#")[0]))
            
            #3 - HALFS (somewhere in A-loop - 3p)
            # Here the starting point can be somewhere between the boundaries of the Anticodon-loop and end at the last base of the full length molecule (CCA). 
            elif (ending_point == lref and Al_start <= int(d_data[3]) <= Al_end):
                nn0 = str(d_data[0]+"_"+str(d_data[3]+1)+":"+str(d_data[3]+llib)+"-3H")
                seq2 = final_seq(d_data[1], d_data[3], d_data[5], lref, llib, d_data[4], name, nn0, sb, d_data[2],keep_ss)
                D[4].append((nn0, seq2, d_data[2], sample, d_data[5], d_data[4], name.split("#")[0]))
            
            #5p - tRFs
            # Here the starting point is at the beginning and the end anywhere in the D-loop
            elif (int(d_data[3]) == int(0) and ending_point <= Dl_end):
                nn0 = str(d_data[0]+"_"+str(d_data[3]+1)+":"+str(d_data[3]+llib)+"-5P")
                seq2 = final_seq(d_data[1], d_data[3], d_data[5], lref, llib, d_data[4], name, nn0, sb, d_data[2],keep_ss)
                D[1].append((nn0, seq2, d_data[2], sample, d_data[5], d_data[4], name.split("#")[0]))

            #3p - tRFs
            # Here the starting point should be anywhere in the T-loop and ends at the last base of the refence molecule (CCA)
            elif (ending_point == lref and int(d_data[3]) >= Tl_start):
                nn0 = str(d_data[0]+"_"+str(d_data[3]+1)+":"+str(d_data[3]+llib)+"-3P")       
                seq2 = final_seq(d_data[1], d_data[3], d_data[5], lref, llib, d_data[4], name, nn0, sb, d_data[2],keep_ss)
                D[2].append((nn0, seq2, d_data[2], sample, d_data[5], d_data[4], name.split("#")[0])) 

            #Mid - tRFs
            # Everything that is left 
            else:
                nn0 = str(d_data[0]+"_"+str(d_data[3]+1)+":"+str(d_data[3]+llib)+"-D2T")
                seq2 = final_seq(d_data[1], d_data[3], d_data[5], lref, llib, d_data[4], name, nn0, sb, d_data[2],keep_ss)
                D[5].append((nn0, seq2, d_data[2], sample, d_data[5], d_data[4], name.split("#")[0]))

    return 

# Reading bowtie output files and saving the necessary information.
def reading_bowtie_file(bt_file):
    
    num_mut = 0
    data = []
    with open(bt_file,"r") as fp:
        for line in fp:

            # Obtaining basic information for each read
            counts = int(line.split("\t")[0].split("-")[1].strip())
            annotation = line.split("\t")[2].strip()
            sequence = line.split("\t")[4].strip().upper()
            s_point = int(line.split("\t")[3].strip())
            
            # Saving mutations 
            if len(line.split()) == 7: # 0 mismatches
                mut = " "
                num_mut = 0
            else:
                mut = line.split("\t")[7].strip()
                if line.split("\t")[1].strip() == "+":
                    if len(mut.split(":")) == 2: # 1 mismatch 
                        num_mut = 1
                    if len(mut.split(":")) == 3: # 2 mismatches
                        num_mut = 2
                elif line.split("\t")[1].strip() == "-": # In case of a reverse (-) strand           
                    if len(mut.split(":")) == 2: # 1 mismatch
                        mut = str(len(sequence)-int(mut.split(":")[0])-1) +":"+ mut.split(":")[1]
                        num_mut = 1
                    elif len(mut.split(":")) == 3: # 2 mismatches
                        mut1 = str(len(sequence)-int(mut.split(",")[0].split(":")[0])-1) +":"+ mut.split(",")[0].split(":")[1]
                        mut2 = str(len(sequence)-int(mut.split(",")[1].split(":")[0])-1) +":"+ mut.split(",")[1].split(":")[1]
                        mut = mut1+","+mut2
                        num_mut = 2
            data.append((annotation, sequence, counts, s_point, mut, num_mut))


    return data

# Obtaining secondary structure information for each read. This file looks like a .fastq file (kind of)
def secondary_structure_info(infile):

    # Initialise variables
    a = 0
    b = 0

    # Returns secondary structure info (a) and secondary structure (b) for each read.
    for key, vals in secstructLib.iteritems():
        if infile == key:
            a = vals[0]
            b = vals[1]
    

    return (a,b)

# Filling gaps at the beginning and at the end, where it is necessary. 
def final_seq(inseq, start, num_mut, lref, llib, mut, name, nn0, sb, counts, keep_ss):
    
    # Output file: fasta file containing the modified names, but the raw sequences.
    # Creating a fasta file with  modified names, but the raw sequences.
    # with open(os.path.join(results_folder,name.split("#")[0]+".fasta"),"a+") as f: f.write(">"+nn0+" "+str(counts)+"\n"+inseq+"\n")

    seqq=0
    seq_in = 0
    
    if num_mut == 0:
        seq_in = inseq
    else:
        seq_in = creating_mut(inseq,mut, num_mut)

    diff = int(lref - (start+llib))

    if start == 0 and diff == 0:
        seqq = seq_in
        seq_ss = sb
    elif start == 0 and diff !=0:
        seqq = seq_in+("-"*(lref - (start+llib)))
        seq_ss = sb[:-diff]
    elif start != 0 and diff ==0:
        seqq = ("-"*start)+seq_in
        seq_ss = sb[start: ]
    elif start != 0 and diff !=0:
        seqq = ("-"*start)+seq_in+("-"*diff)
        seq_ss = sb[start: -diff]

    # Output file 2: Secondary Structure of each read
    if keep_ss == True: 
        with open(os.path.join(results_folder,name.split("#")[0]+".ss"),"a") as aout: aout.write(">"+nn0+"\n"+seq_ss+"\n")
    
    # Output file 3: fasta file containing the modified names and sequences.
    with open(os.path.join(results_folder,name.split("#")[0]+".mod"),"a+") as aout: aout.write(">"+nn0+" "+str(counts)+"\n"+seqq+"\n")
    
    return seqq

# Find the position of the mutations, change to the correct mutation and make the letter lower. 
def creating_mut(sequence, mut, number_of_mut):
    
    seqq1 = 0
    if number_of_mut == 1:
        seqq1 = sequence[ :int(mut.split(",")[0].split(":",1)[0])] + sequence[int(mut.split(",")[0].split(":",1)[0])].lower() + sequence[int(mut.split(",")[0].split(":",1)[0])+1: ]
    elif number_of_mut == 2:
        seq2 = sequence[ :int(mut.split(",")[0].split(":",1)[0])] + sequence[int(mut.split(",")[0].split(":",1)[0])].lower() + sequence[int(mut.split(",")[0].split(":",1)[0])+1: ]
        seqq1 = seq2[ :int(mut.split(",")[1].split(":",1)[0])] + seq2[int(mut.split(",")[1].split(":",1)[0])].lower() + seq2[int(mut.split(",")[1].split(":",1)[0])+1: ]
    
    return seqq1

# Output file: summarised excel file.
def write_excel(argment):
    
    # Multiprocessing
    # names = {0:("FL|exact","5p|exact","3p|exact","5h|exact","3h|exact","D2T|exact"),1:("FL|1","5p|1","3p|1","5h|1","3h|1","D2T|1"),2:("FL|2","5p|2","3p|2","5h|2","3h|2","D2T|2")}
    # wbpage = book.get_worksheet_by_name(names[argment[2]][argment[3]])
    
    wbpage = argment[2]
    

    newlist = []
    st = [(item[0],item[1]) for item in argment[0]]

    # Create tRNA annotation and sequences for the excel. 
    for i in st:
        if i not in newlist:
            newlist.append(i)

    # # Creating an empty matrix that will host each table-sheet.
    Matrix = np.zeros((len(newlist), argment[1]),dtype=int)   

    for idx, val in enumerate(newlist):
        wbpage.write((idx+1), 0, val[0])
        wbpage.write((idx+1), 1, val[1])

        for val1 in argment[0]:
            if (val1[0],val1[1]) == val:
                Matrix[newlist.index((val1[0],val1[1])), (val1[3]-2)] = int(val1[2])

    for (x,y), value in np.ndenumerate(Matrix):
        wbpage.write(x+1, y+2, value)

    return 


############# STATS #############

# Production of basic and compact analysis
def statistics(keep_mod):
    
    #DONE##SUMMARY OF EACH SAMPLE#############################################################################################################
    for path, subdirs, files in os.walk(results_folder):
        for name in files:
            if name.endswith(".mod"):
                totalCounts = sum_per_tRNA(name, os.path.join(path,name))
                mm_stats(name[:-4], totalCounts)
                nrm_em[name[:-4]] = sum(totalCounts.values())
                # with open(os.path.join(results_folder,"total_counts.csv"),"a") as aout: aout.write(name[:-4]+"\t"+str(sum(totalCounts.values()))+"\n")
                if keep_mod == False: os.remove(os.path.join(path,name))

    ##COMPACT EXPRESSION MATRIX##########################################################################################################
    # Creating a more compact excel this time by category of tRNA and tRF and not taking into account the mismatches.
    book2 = xlsxwriter.Workbook(os.path.join(results_folder, "Expression_Matrix(compact).xlsx"))
    sh1 = book2.add_worksheet("Full Length")
    sh2 = book2.add_worksheet("5 halfs")
    sh3 = book2.add_worksheet("3 halfs")
    sh4 = book2.add_worksheet("5 prime")
    sh5 = book2.add_worksheet("3 prime")
    sh6 = book2.add_worksheet("D2T")

    b2sheets = [sh1,sh2,sh3,sh4,sh5,sh6]
    
    wb = xlrd.open_workbook(os.path.join(results_folder,"Expression_Matrix(extensive).xlsx"))
    sheets_comp = wb.sheets()

    # Reading the tRNA-tRF.xlsx and write the new excel file.
    rowidx = 0
    for idx, sh1 in enumerate(sheets_comp):
        for idx2, sh2 in enumerate(b2sheets):
            if(idx/3) == idx2:
                if idx in [0,3,6,9,12,15]:
                    for row in range(sh1.nrows):
                        for col in range(sh1.ncols):
                            sh2.write(row, col, sh1.cell(row, col).value)
                    rowidx += sh1.nrows
                elif idx in [1,4,7,10,13,16]:
                    for row in range(sh1.nrows):
                        if row == 0: pass
                        else:
                            for col in range(sh1.ncols):
                                sh2.write(rowidx+row-1, col, sh1.cell(row, col).value)
                    rowidx += sh1.nrows                
                elif idx in [2,5,8,11,14,17]:
                    for row in range(sh1.nrows):
                        if row == 0: pass
                        else:
                            for col in range(sh1.ncols):
                                sh2.write(rowidx+row-2, col, sh1.cell(row, col).value)
                    rowidx = 0    
	

    return

# Creating the normilised version of the compact expression matrix
def normalised_matrix():
    ##NORMALISED COMPACT EXPRESSION MATRIX##########################################################################################################
    book3 = xlsxwriter.Workbook(os.path.join(results_folder, "Expression_Matrix(normalised).xlsx"))
    sh1 = book3.add_worksheet("Full Length")
    sh2 = book3.add_worksheet("5 halfs")
    sh3 = book3.add_worksheet("3 halfs")
    sh4 = book3.add_worksheet("5 prime")
    sh5 = book3.add_worksheet("3 prime")
    sh6 = book3.add_worksheet("D2T")

    b3sheets = [sh1,sh2,sh3,sh4,sh5,sh6]
    
    wb = xlrd.open_workbook(os.path.join(results_folder, "Expression_Matrix(compact).xlsx"))
    sheets_comp = wb.sheets()

    for idx, sh11 in enumerate(sheets_comp):
		for idx2, sh22 in enumerate(b3sheets):
			if idx == idx2:
				for row in range(sh11.nrows):
					for col in range(sh11.ncols):
						if col == 0 or col == 1:
							sh22.write(row, col, sh11.cell(row, col).value)
						else:
							for idx, keys in enumerate(nrm_em):
								if sh11.cell(0, col).value.split("#")[0] == keys:
									sh22.write(row, idx+2, sh11.cell(row, col).value)
									if row != 0:
										sh22.write(row, idx+2, float(sh11.cell(row, col).value)/(float(nrm_em[keys])/1e06))

    return

# Statistical analysis of mm  
def mm_stats(filenm, totalCounts):
        
    
    categorised = {}
    # Categorising detected molecules with 1mm
    for vals in dics1:
        for elm in vals:
            if filenm in elm:
            	counts = int(elm[2])
                annotation = elm[0].split("_")[1].strip()
                trna = elm[0].split("_")[0].strip()
                mut = str(int(elm[5].split(":")[0])+int(elm[0].split("_")[1].split(":")[0]))+":"+elm[5].split(":")[1]
            	
            	if (trna, mut) in categorised:
                    categorised[(trna, mut)].append((counts, annotation))
                else:
                    categorised[(trna, mut)] = [(counts, annotation)]
	
	# Categorising detected molecules with 2mm
    for vals2 in dics2:
        for elm2 in vals2:
            if filenm in elm2:
                counts2 = int(elm2[2])
                annotation2 = elm2[0].split("_")[1].strip()
                trna2 = elm2[0].split("_")[0].strip()
                mut1 = str(int(elm2[5].split(",")[0].split(":")[0])+int(elm2[0].split("_")[1].split(":")[0]))+":"+elm2[5].split(",")[0].split(":")[1]
                mut2 = str(int(elm2[5].split(",")[1].split(":")[0])+int(elm2[0].split("_")[1].split(":")[0]))+":"+elm2[5].split(",")[1].split(":")[1]

                if (trna2, mut1) in categorised:
                    categorised[(trna2, mut1)].append((counts2, annotation2))
                else:
                    categorised[(trna2, mut1)] = [(counts2, annotation2)]
                
                if (trna2, mut2) in categorised:
                    categorised[(trna2, mut2)].append((counts2, annotation2))
                else:
                    categorised[(trna2, mut2)] = [(counts2, annotation2)]     
    
    # Saving the data           
    with open(os.path.join(results_folder, filenm+".modstat.csv"),"w") as stat_out:
        stat_out.write("tRNA_name\tMM_position\ttMM_type\ttotalRC_MM\ttotalRC_tRNA\ttpercentage\n")
        for a, b in ns.natsorted(categorised.iteritems()):
        	total_mol = len(b)
        	resulted_counts = sum([pair[0] for pair in b])
       		if a[0] in totalCounts:
       			for elms in b:
        			stat_out.write(a[0]+"_"+elms[1]+"\t"+a[1].split(":")[0]+"\t"+a[1].split(":")[1].split(">")[0]+"->"+
                                   a[1].split(":")[1].split(">")[1]+"\t"+str(resulted_counts)+"\t"+str(totalCounts[a[0]])+"\t"+
                                   str(round((float(resulted_counts)/float(totalCounts[a[0]])),5))+"\n")
    return

# Categorising tRNAs and tRFs per tRNA molecule. 
def sum_per_tRNA(nameoffile, pathoffile):
    
    # Reading the modified files. 
    allseqs = {}
    totalCounts = {}

    with open(pathoffile) as fin:
        for line in fin:
            
            if line.startswith(">"):
                seqname = line[1:].split("_")[0].strip()
                num = line.split(" ")[1].strip()
            else:
                seq = line.strip()
                if seqname in allseqs:
                    allseqs[seqname].append((seq, int(num)))
                else:
                    allseqs[seqname] = [(seq, int(num))]                 

    # Categorising and summarising the reads.
    with open(os.path.join(results_folder, nameoffile[:-4]+".summary"), "a+") as fout:
        for key, val in ns.natsorted(allseqs.iteritems()):
            val = sorted(val, key=itemgetter(1), reverse=True)
            for idx, dict_val in enumerate(val,1):
                result = sum(n for _, n in val)
                countofmol = len(allseqs[key])
                if key in secstructLib:
                    if countofmol == 1:
                        fout.write(key+"\n"+secstructLib[key][2].upper()+" "+str(result)+"\n"+secstructLib[key][1]+"\n"+dict_val[0]+" "+str(dict_val[1])+(2*"\n"))
                    elif (countofmol != 1 and idx == 1):
                        fout.write(key+"\n"+secstructLib[key][2].upper()+" "+str(result)+"\n"+secstructLib[key][1]+"\n"+dict_val[0]+" "+str(dict_val[1]))
                    elif (countofmol != 1 and idx != countofmol):
                        fout.write("\n"+dict_val[0]+" "+str(dict_val[1]))
                    elif (countofmol != 1 and idx == countofmol):
                        fout.write("\n"+dict_val[0]+" "+str(dict_val[1])+("\n"*2))  
                
                totalCounts[key]=result
    
    return totalCounts

############# MAIN #############

def main():

    (argms, infiles, qcfiles) = parse_commandline()
    numoffiles = len(infiles)

    ### 1. Quality_Control

    if (argms.getboolean("Quality_Control", "no_QC")) == True:
        print "\n1. NO Quality Control reports will be generated. This step has been skipped!"
        shutil.rmtree(qc_folder)
    else:
        OC_stime = datetime.now()
        print "\n%s | 1. Quality Control reports for the data are being generated: in progress .." %OC_stime.strftime("%d.%m.%Y %H:%M:%S")
        print "------------------------------------------------------------"
        print "Number of files: %s | Number of processors used: %s" %(numoffiles,thrds)
        quality_control(qcfiles)
        QC_etime = datetime.now() 
        print "\tQuality Control reports: Completed Successfully in %s" %(QC_etime - OC_stime)

    ### 2. Preprocessing
    if (argms.getboolean("Preprocessing", "no_Preprocessing")) == True: print "\n2. NO Preprocessing of the data will be performed. This step has been skipped!"
    else:
        # Obtaining the adapter sequence
        adapter = "TGGAATTCTCGGGTGCCAAGGAACTC" if argms.get("Preprocessing", "adapter_sequence") is "-" else "".join(argms.get("Preprocessing", "adapter_sequence"))
        # preprocessing in parallel using multiple cores
        Preprocessing_stime = datetime.now()
        print "\n%s | 2. Preprocessing the data: in progress .." %str(Preprocessing_stime.strftime("%d.%m.%Y %H:%M:%S"))
        print "------------------------------------------------------------"
        prepro_args = [(file, adapter, argms.getboolean("Preprocessing", "keep_short_reads"))for file in infiles]
        pool = mp.Pool(processes=int(thrds))
        pool.map(preprocessing, prepro_args)
        Preprocessing_etime = datetime.now()
        print "\tPreprocessing the data: Completed Successfully in %s" %(Preprocessing_etime - Preprocessing_stime)
    
    

    ### 3. Align against the reference libraries
    aligning_stime = datetime.now()
    print "\n%s | 3. Aligning against the reference libraries: in progress .." %(aligning_stime.strftime("%d.%m.%Y %H:%M:%S"))
    print "------------------------------------------------------------"

    libraries = []
    for i in libs:
        bowtie_ind = os.path.join(argms.get("Aligning", "tLibs"), i)
        libraries.append(bowtie_ind)

    for mism in range(3):
        for i in range(len(libraries)):
            # tLib: mature | mm: 0
            if (i==0 and mism==0):
                # print libraries[i], mism
                align(processed_data_folder,libraries[i], libraries[i], mism, argms.get("Aligning", "Ref_Genome"),argms.get("Aligning", "tLibs"))
            # tLib: premature | mm: 0    
            elif (i==1 and mism==0):
                # print libraries[i], mism
                align(aligned_folder,libraries[i], libraries[i-1], mism, argms.get("Aligning", "Ref_Genome"),argms.get("Aligning", "tLibs"))
            # tLib: mature | mm: 1
            elif (i == 0 and mism == 1):
                # print libraries[i], mism
                align(aligned_folder,libraries[i], libraries[i-1], mism, argms.get("Aligning", "Ref_Genome"),argms.get("Aligning", "tLibs"))
            # tLib: mature | mm: 2
            elif (i==0 and mism==2):
                # print libraries[i], mism
                align(aligned_folder,libraries[i], libraries[i], mism, argms.get("Aligning", "Ref_Genome"),argms.get("Aligning", "tLibs")) 
    
    aligning_etime = datetime.now()
    print "\tAligning: Completed Successfully in %s" %(aligning_etime - aligning_stime)

    ref_tLib= os.path.join(argms.get("Aligning", "tLibs"),"reflibs_mu.fa")
    ref_secstr =  os.path.join(argms.get("Aligning", "tLibs"), "sslib_mu.fa")

    ### 4. Classification
    clas_stime = datetime.now()
    print "\n%s | 4. Classification of the data: in progress .." %(clas_stime.strftime("%d.%m.%Y %H:%M:%S"))
    print "------------------------------------------------------------"
    main_analyis(argms.getboolean("Classification", "keep_ss_files"),ref_tLib, ref_secstr)
    clas_etime = datetime.now()
    print "\tClassification: Completed Successfully in %s" %(clas_etime - clas_stime)
    print "\n5. Expression Matrix has been generated"
    
    ### 6. Stats - Output files 
    if (argms.getboolean("Basic_Statistics", "no_stats")) == True: print "\n4. No statistics will be generated. This step has been skipped!"
    else: 
        stat_stime = datetime.now()
        print "\n%s | 6. Basic statistical analysis of the data: in progress .." %(clas_stime.strftime("%d.%m.%Y %H:%M:%S"))
        print "------------------------------------------------------------"
        statistics(argms.getboolean("Basic_Statistics", "keep_modified"))
        if argms.getboolean("Basic_Statistics", "keep_normalised_em") == True: normalised_matrix()
        stat_etime = datetime.now()
        print "\tStatistical analysis: Completed Successfully in %s" %(stat_etime - stat_stime)


    # Final report
    try:
        # Running MultiQC to get a summed Quality Control, summed evaluation of featureCounts and STAR aligner.
        run_mc = " ".join(["multiqc", "-q", "-o", reports_folder, "-n", "FinalReport", reports_folder])
        subprocess.call(run_mc, shell=True)
        #Deleting .zip files from FastQC reports.
        for path, subdirs, files in os.walk(reports_folder):
            for name in files:
                if name.endswith("fastqc.zip"):
                    zipped_reports = os.path.join(path, name)
                    os.remove(zipped_reports)
        # Removing byproducts produced by MultiQC.           
        shutil.rmtree(os.path.join(reports_folder,"FinalReport_data"))
        if (argms.getboolean("Aligning", "keep_aligned")) == False: shutil.rmtree(aligned_folder)
        
    except:
        pass

    end_time = datetime.now()

    print "\n-------------------------------------------"
    print "-------------------------------------------"
    print "Number of analysed files: %d" %numoffiles
    print "Analysis completed in: %s" %(end_time - start_time)

if __name__ == "__main__":
    main()
