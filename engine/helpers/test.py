
import re
import os
import subprocess
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo 
import json
import Bio.Align.Applications as apps
import tracemalloc as malloc
import timeit
import csv
import time
import numpy as np

def settings():
    file = open("settings.json")
    settings = json.load(file)
    return settings 

def generate_alignments():
    #intitials sequences statistics
     #min length
     
     #max length
    programs = settings().get("program_commands")
    src_files = os.listdir("sequences")
    #data collection into csv file
    fields=["Program","Seq","Exec Time(ms)","Mem Usage(Kb)","Similarity Score","Seq Mean length","Mean PID(%)", "Total Seq Count"]
    rows = []
    print("generating alignments")
    for program in programs:
         print (program)
         #the commands
         if program == "linsi":
            for src_file in src_files:
                    out_file = "{}_{}.fasta".format(program,src_file.split(".")[0])
                    file = open("./Data/{}".format(out_file),"w+")
                    malloc.start()
                    start =  timeit.default_timer()
                    proc = subprocess.Popen(["linsi","sequences/{}".format(src_file)],stdout=file,encoding="utf-8")
                    mem= malloc.get_tracemalloc_memory()
                    end = timeit.default_timer()
                    exec_time = end-start
                    malloc.stop()
                    time.sleep(20)
                    data = parse_file("./Data/{}".format(out_file),"fasta")
                    res = score_aln(data["sequences"])
                    mPID = res["identities"]/data["total_lengths"]*100
                    rows.append([program,src_file.split(".")[0],"{:.4f}".format(exec_time),\
                       mem,res["Sim Score"],data["mean_length"],"{:.4f}".format(mPID),res["Seq Count"]])
                    

         if program == "muscle":
             for src_file in src_files:
                    out_file = "{}_{}.fasta".format(program,src_file.split(".")[0])
                    file = open("./Data/{}".format(out_file),"w+")
                    infile = "./sequences/{}".format(src_file)
                    malloc.start()
                    start = timeit.default_timer()
                    proc = subprocess.Popen(["muscle","-in",infile],stdout=file,encoding="utf-8") 
                    mem = malloc.get_tracemalloc_memory()
                    end = timeit.default_timer()
                    exec_time = end-start
                    malloc.stop()
                    time.sleep(20)
                    data = parse_file("./Data/{}".format(out_file),"fasta")
                    res = score_aln(data["sequences"])
                    mPID = res["identities"]/data["total_lengths"]*100
                    rows.append([program,src_file.split(".")[0],"{:.4f}".format(exec_time),\
                       mem,res["Sim Score"],data["mean_length"],"{:.4f}".format(mPID),res["Seq Count"]])

         if program == "t_coffee":
              for src_file in src_files:
                   out_file = "./Data/{}_{}.fasta".format(program,src_file.split(".")[0])
                   infile = "./sequences/{}".format(src_file)
                   tcoffee_cline = apps.TCoffeeCommandline(infile=infile,output= "fasta_aln",outfile= out_file)
                   malloc.start()
                   start = timeit.default_timer()
                   tcoffee_cline()
                   mem = malloc.get_tracemalloc_memory()
                   end = timeit.default_timer()
                   exec_time = end-start
                   malloc.stop()
                   data = parse_file(out_file,"fasta")
                   res = score_aln(data["sequences"])
                   mPID = res["identities"]/data["total_lengths"]*100
                   rows.append([program,src_file.split(".")[0],"{:.4f}".format(exec_time),\
                       mem,res["Sim Score"],data["mean_length"],"{:.4f}".format(mPID),res["Seq Count"]])

         if program ==  "kalign":
              for src_file in src_files:
                   out_file = "./Data/{}_{}.fasta".format(program,src_file.split(".")[0])
                   file = open("{}".format(out_file),"w+")
                   infile = "./sequences/{}".format(src_file)
                   malloc.start()
                   start = timeit.default_timer()
                   proc = subprocess.Popen(["kalign",infile],stdout=file,encoding="utf-8")
                   mem = malloc.get_tracemalloc_memory()
                   end = timeit.default_timer()
                   exec_time = end-start
                   malloc.stop()
                   time.sleep(40)
                   file.close()
                   data = parse_file(out_file ,"fasta")
                   res = score_aln(data["sequences"])
                   mPID = res["identities"]/data["total_lengths"]*100
                   rows.append([program,src_file.split(".")[0],"{:.4f}".format(exec_time),\
                       mem,res["Sim Score"],data["mean_length"],"{:.4f}".format(mPID),res["Seq Count"]])

         if  program == "dialign2-2": 
              for src_file in src_files: 
                    out_file = "./Data/{}_{}".format(program,src_file.split(".")[0]) 
                    infile = "./sequences/{}".format(src_file)
                    cmd = apps.DialignCommandline(input=infile,fn=out_file,fa=True)
                    malloc.start()
                    start = timeit.default_timer()
                    cmd()
                    end = timeit.default_timer()
                    exec_time = end-start
                    mem = malloc.get_tracemalloc_memory()
                    malloc.stop()
                    data = parse_file(out_file+".fa","fasta")
                    res = score_aln(data["sequences"])
                    mPID = res["identities"]/data["total_lengths"]*100
                    rows.append([program,src_file.split(".")[0],"{:.4f}".format(exec_time),\
                       mem,res["Sim Score"],data["mean_length"],"{:.4f}".format(mPID),res["Seq Count"]])

         if program == "probcons":
              for src_file in src_files:
                   out_file = "./Data/{}_{}.fasta".format(program,src_file.split(".")[0])
                   infile = "./sequences/{}".format(src_file)
                   cmd = apps.ProbconsCommandline(input=infile)
                   malloc.start()
                   start = timeit.default_timer()
                   stdout,stderr = cmd()
                   end = timeit.default_timer()
                   exec_time = end-start
                   mem = malloc.get_tracemalloc_memory()
                   malloc.stop()
                   with open(out_file, "w+" ) as handle:
                        handle.write(stdout)

                   data = parse_file(out_file,"fasta")
                   res = score_aln(data["sequences"])
                   mPID = res["identities"]/data["total_lengths"]*100
                   rows.append([program,src_file.split(".")[0],"{:.4f}".format(exec_time),\
                       mem,res["Sim Score"],data["mean_length"],"{:.4f}".format(mPID),res["Seq Count"]])     
         # Prank not able to generate alignment files
         if program == "clustalo":
              for src_file in src_files:
                   out_file = "./Data/{}_{}".format(program,src_file.split(".")[0])
                   infile = "./sequences/{}".format(src_file)
                   clustalo = apps.ClustalOmegaCommandline(infile=infile, outfile=out_file, verbose=True, auto=True)
                   start = timeit.default_timer()
                   clustalo()
                   end = timeit.default_timer()
                   exec_time = end-start
                   mem = malloc.get_tracemalloc_memory()
                   malloc.stop()
                   data = parse_file(out_file,"fasta")
                   res = score_aln(data["sequences"])
                   mPID = res["identities"]/data["total_lengths"]*100
                   rows.append([program,src_file.split(".")[0],"{:.4f}".format(exec_time),\
                       mem,res["Sim Score"],data["mean_length"],"{:.4f}".format(mPID),res["Seq Count"]])
                   
         with open("perfomance_data.csv","w+") as csvfile:
              CSVwriter = csv.writer(csvfile)
              CSVwriter.writerow(fields)
              CSVwriter.writerows(rows)

                      

def score_aln(sequences):
          count = 0
          scores = []
          pids = []
          results = dict()
          while count < len(sequences):
               #remove the sequence from list
               p_score = 0
               p_identity = 0
               sequence = sequences.pop(count)
               #now loop in the remaining sequences and compare the residuals
               for this_seq in sequences:
                    i=0
                    while i < len(sequence):
                         #genetate a tuple of pairwise residuals
                         pairwise = (sequence[i].upper(),this_seq[i].upper())
                         if "-" in pairwise:
                              print("processing gaps ...")
                         else:
                              if sub_matrix(62).get(pairwise) is not None:
                                    p_score +=  sub_matrix(62).get(pairwise)
                                    if pairwise[0]== pairwise[1]:
                                         p_identity +=1
                              else:
                                    #print(pairwise)
                                    p_score += sub_matrix(62).get(rev_tuple(pairwise))
                                    if pairwise[0]== pairwise[1]:
                                         p_identity +=1
                         #end the loop
                         i += 1
               #push the pairwise into the list for scores
                    scores.append(p_score)
                    pids.append((p_identity/len(sequences))*100)          
               #replace the deleted sequences
               sequences.insert(count,sequence)          
               count += 1
          scaled_dwn_scores = list(map(lambda x: x/1,scores))  
          results.update({"Seq Count": len(sequences)})   
          results.update({"Sim Score": np.mean(scaled_dwn_scores)})
          results.update({"identities": np.mean(pids)})
          return results
            
"""
with open("./Data/BBA0011_Mafft_linsi_aligned.fasta") as file:
     content = str(file.read())
     sequences = re.split(">seq[0-9]*",content)
     individual_seqs = []
     count = 1
     del(sequences[0])
     print("SeqID,          Length")
     for sequence in sequences:
          psequence = sequence.strip().replace("\n","")
          individual_seqs.append(psequence)
          #print the length of each seq
          print("{}         {}".format(count,len(psequence)))
          count += 1
     print(individual_seqs)

     score_aln(individual_seqs)
 
     #Loop through the list of individual sequences
"""  
def parse_file(file,format):
     sequences = []
     seq_IDs = []
     seq_lengths = []
     seq_details = dict()
     for seq_record in SeqIO.parse(file, format):
          seq_IDs.append(seq_record.id)
          sequences.append(seq_record.seq)
          seq_lengths.append(len(seq_record.seq))

     seq_details.update({"seq_IDS":seq_IDs})
     seq_details.update({"sequences":sequences})
     seq_details.update({"total_lengths":sum(seq_lengths)})
     seq_details.update({"mean_length":np.mean(seq_lengths)})
     #call another function
     return seq_details

# method that returns two substitition matrices dictionaries
def sub_matrix(percent):
    if percent==50:
        return MatrixInfo.blosum50
    else:
        return MatrixInfo.blosum62

# a method that reverses a tuple        
def rev_tuple(t):
    return tuple(reversed(t))


#generate_alignments()
#data = parse_file("./Data/BBA0011_Mafft_linsi_aligned.fasta","fasta")
#score = score_aln(data["sequences"])     
#print(sub_matrix(50))
#print("===========================================================================================")
#print(sub_matrix(62))
