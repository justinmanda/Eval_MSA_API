from flask import Flask, jsonify, request
from flask_cors import CORS
import uuid
import tracemalloc as malloc
import timeit
import time
import csv
import time
import numpy as np
from stats import Stat
import threading
import os
from align import Align
from utils import Config
from multiprocessing import Process
import time
from flask_mail import Mail

app = Flask(__name__)
app.config['MAIL_SERVER']='smtp.gmail.com'
app.config['MAIL_PORT'] = 465
app.config['MAIL_USERNAME'] = 'justinmandah@gmail.com'
app.config['MAIL_PASSWORD'] = 'justin4686'
app.config['MAIL_USE_TLS'] = False
app.config['MAIL_USE_SSL'] = True
mail = Mail(app)
CORS(app)

@app.route("/",methods=["GET"])
def welcome():
    return "<p> This is a program that evaluates the accuracy of the multiple sequence alignment programs<br/>\
            based on the blosum62 substitution matrix. In essence, this program assesses how the alignment<br>programs\
            uphold the concept of homology during the genaration of the sequence alignment. This program is <br> \
            only suited to work with protein sequences. Please input the required information below to give it a try</p> "

            

@app.route("/api/v1/align/sequences",methods=["POST"])
def generate_alignment():
    data = request.get_json()
    file_uuid = uuid.uuid4().hex
    file_name = data.get("filename")
    file = "./tmp/{}_sequence.fasta".format(file_uuid)
    with open(file,'w',encoding = 'utf-8') as f:
        f.write(data.get("sequence"))

    print(data)
    start_time = timeit.default_timer() 
    programs = data.get("programs") 
    task_cb = Process(target=run_data_processor, args=(file_uuid, programs,file,file_name))
    task_cb.start()
    end_time = timeit.default_timer()
    print("{} minutes taken to align with {} aligners".format((end_time-start_time)/60,len(Config.settings()["programs"])))
    return '',200     

def run_data_processor(file_uuid,programs,file,file_name):
    align = Align(file_uuid,file_name) 
    threads = []
    for program in programs:
        threads.append(threading.Thread(target=align.generate_alignment, args=(program.lower(),file,)))
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join() 
    #Now generate the scores 
    #check that all files have been generated
    count = len(os.listdir(Config.settings()["aligned_seqs_path"]))
    while count != 7:
        time.sleep(2.4)
        count = len(os.listdir(Config.settings()["aligned_seqs_path"]))
    #end of task blocking code    
    scoring_threads = []
    for program in programs:
        scoring_threads.append(threading.Thread(target=align.generate_scores, args=(program.lower(),)))
    for thred in scoring_threads:
        thred.start()
    for thred in scoring_threads:
        thred.join()  
    #align.generate_scores(program)
    stats = align.generate_stats()
#Send email of results       
   
app.run()            
