import json
from Bio import SeqIO
from Bio import AlignIO
import numpy as np
class Config(object):
    @staticmethod
    def settings():
        file = open("settings.json")
        settings = json.load(file)
        return settings  
    
class Builtins:
    def __init__(self):
        pass

    def rev_tuple(self,t):
        return tuple(reversed(t))

class Utils: 
    def __init__(self,path,format): 
        self.path = path
        self.format = format

    def parse_file(self):
        sequences = []
        seq_IDs = []
        seq_lengths = []
        seq_details = dict()
        for seq_record in SeqIO.parse(self.path, self.format):
            seq_IDs.append(seq_record.id)
            sequences.append(seq_record.seq)
            seq_lengths.append(len(seq_record.seq))

        seq_details.update({"seq_IDS":seq_IDs})
        seq_details.update({"sequences":sequences})
        seq_details.update({"total_lengths":sum(seq_lengths)})
        seq_details.update({"mean_length":np.mean(seq_lengths)})
        #call another function
        return seq_details
    
    def get_alignment_obj(path):
        try:
            alignment = AlignIO.read(path, "fasta")
        except Exception as e:
            print(f"An erreor occurred: {e}")
        return alignment
