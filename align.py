#from turtle import pd

from utils import Config
from utils import Utils
import os
import re
import subprocess
from Bio import SeqIO
#from Bio.SubsMat import MatrixInfo 
import json
import Bio.Align.Applications as apps
import tracemalloc as malloc
import timeit
import time
import csv
import time
import numpy as np
from stats import Stat
import pandas as pd
import asyncio
import threading
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import uuid

class Align:
    def __init__(self,aln_uuid,file_name):
        self.programs = Config.settings()["programs"]
        self.input_files_path = Config.settings()["raw_seqs_path"]
        self.output_files_path = Config.settings()["aligned_seqs_path"]
        self.result_path =  Config.settings()["results_path"]
        self.all_scores = {}
        self.fields=["Program","Seq","Exec Time(ms)","Mem Usage(Kb)","Similarity Score","Seq Mean length","Mean PID(%)", "Total Seq Count"]
        self.rows = []
        self.exec_times = {}
        self.mem_use = {}
        self.aln_uuid = aln_uuid
        self.file_name = file_name
        #remove old alignments

    def generate_alignment(self,program,file):
        src_file = file
        aln_id = self.aln_uuid
        #aln_file = Config.settings()['aligned_seqs_path']+"program_"+src_file.split(".")[0]+".fasta"
        score_data = None
        print(file)
        #data collection into csv file
        print("generating alignments")
        #for program in programs:
        print (program)
        #the commands
        if program == "linsi":
                out_file = "{}_{}.fasta".format(program,aln_id)
                file = open("./Data/{}".format(out_file),"w+")
                malloc.start()
                start = timeit.default_timer()
                proc = subprocess.Popen(["linsi",src_file],stdout=file,encoding="utf-8")
                mem = malloc.get_tracemalloc_memory()
                end = timeit.default_timer()
                exec_time = end-start
                malloc.stop()
                self.exec_times.update({program:exec_time})
                self.mem_use.update({program:mem})
                    
        if program == "muscle":
         
                out_file = "{}_{}.fasta".format(program,aln_id)
                file = open("./Data/{}".format(out_file),"w+")
                infile = src_file
                malloc.start()
                start = timeit.default_timer()
                proc = subprocess.Popen(["muscle","-in",infile],stdout=file,encoding="utf-8") 
                mem = malloc.get_tracemalloc_memory()
                end = timeit.default_timer()
                exec_time = end-start
                malloc.stop()
                self.exec_times.update({program:exec_time})
                self.mem_use.update({program:mem})

        if program == "t_coffee":
            
            out_file = "./Data/{}_{}.fasta".format(program,aln_id)
            infile = src_file
            tcoffee_cline = apps.TCoffeeCommandline(infile=infile,output= "fasta_aln",outfile= out_file)
            malloc.start()
            start = timeit.default_timer()
            tcoffee_cline()
            mem = malloc.get_tracemalloc_memory()
            end = timeit.default_timer()
            exec_time = end-start
            malloc.stop()
            self.exec_times.update({program:exec_time})
            self.mem_use.update({program:mem})   

        if program ==  "kalign":
      
            out_file = "./Data/{}_{}.fasta".format(program,aln_id)
            file = open("{}".format(out_file),"w+")
            infile = src_file
            malloc.start()
            start = timeit.default_timer()
            proc = subprocess.Popen(["kalign",infile],stdout=file,encoding="utf-8")
            mem = malloc.get_tracemalloc_memory()
            end = timeit.default_timer()
            exec_time = end-start
            malloc.stop()
            self.exec_times.update({program:exec_time})
            self.mem_use.update({program:mem})
            
        if  program == "dialign2-2": 
           
                out_file = "./Data/{}_{}.fasta".format(program,aln_id) 
                infile = src_file
                cmd = apps.DialignCommandline(input=infile,fn=out_file,fa=True)
                malloc.start()
                start = timeit.default_timer()
                cmd()
                end = timeit.default_timer()
                exec_time = end-start
                mem = malloc.get_tracemalloc_memory()
                malloc.stop()
                self.exec_times.update({program:exec_time})
                self.mem_use.update({program:mem})
                #remove the unwanted file
                if os.path.isfile(out_file):
                        os.remove(out_file)
                        print("File has been deleted")
                else:
                    print("File does not exist")

        if program == "probcons":
                out_file = "./Data/{}_{}.fasta".format(program,aln_id)
                infile = src_file
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
                self.exec_times.update({program:exec_time})
                self.mem_use.update({program:mem})

        # Prank not able to generate alignment files
        if program == "clustalo":
                out_file = "./Data/{}_{}.fasta".format(program,aln_id)
                infile = src_file
                print("***********************************alng with clustalo************************")
                clustalo = apps.ClustalOmegaCommandline(infile=infile, outfile=out_file, verbose=True, auto=True)
                start = timeit.default_timer()
                clustalo()
                end = timeit.default_timer()
                exec_time = end-start
                mem = malloc.get_tracemalloc_memory()
                malloc.stop()
                self.exec_times.update({program:exec_time})
                self.mem_use.update({program:mem})
                #result = list(filter(lambda x: "." in x, files))
        #remove the file after the alignments
        #if os.path.isfile(src_file):
            #os.remove(src_file)       
        return None

    def generate_scores(self,program):
         print(f"+++++++++++++++++This is the program: {program}++++++++++++++++")
         item = re.compile(program)
         files = os.listdir(Config.settings()["aligned_seqs_path"])
         print(f'here are the file+++++++++ {files}')
         result = list(filter(lambda x: "." in x, files))
         file_name = list(filter(item.match,result))[0] 
         file_path = os.path.join(Config.settings()["aligned_seqs_path"],file_name)
         data = Utils(file_path,"fasta").parse_file()
         res = Stat(data["sequences"]).score_align()
         mPID = res["identities"]
         self.rows.append([program,self.file_name,"{:.4f}".format(self.exec_times.get(program)),\
         self.mem_use.get(program),res["Sim Score"],data["mean_length"],"{:.4f}".format(mPID),res["Seq Count"]])
         self.all_scores.update({program: res["scores"]}) 

    def generate_stats(self):
        with open(self.result_path,"w+") as csvfile:
            CSVwriter = csv.writer(csvfile)
            CSVwriter.writerow(self.fields)
            CSVwriter.writerows(self.rows)

        #score_data = pd.DataFrame.from_dict(self.all_scores)
        #score_data.transpose() 
        #print(score_data)
        os.system(f"rm {Config.settings()['aligned_seqs_path']}*_{self.aln_uuid }*")
        stat_result = Stat.generate_stats(self.all_scores)
        Stat.sort_results()  
        with open("sorted_metrics.csv","a+") as csvfile:
               CSVwriter = csv.writer(csvfile)
               CSVwriter.writerow(stat_result)
              
