from Bio import SeqIO
from Bio.SubsMat import MatrixInfo 
import json
import Bio.Align.Applications as apps
import tracemalloc as malloc
import timeit
import csv
import time
import numpy as np
from utils import *
import pandas as pd
from scipy import stats
from utils import Config
import statsmodels.api as sm
from statsmodels.formula.api import ols
from bioinfokit.analys import stat as biostat
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import math
from Bio import AlignIO
class Stat:
    def __init__(self,seqs):
        self.sequences = seqs

    def score_align(self):
        count = 0
        scores = []
        pids = []
        results = dict()
        while count < len(self.sequences):
            #remove the sequence from list
            p_score = 0
            p_identity = 0
            sequence = self.sequences.pop(count)
            ref_seq = ""
            #now loop in the remaining self.sequences and compare the residuals
            for this_seq in self.sequences:
                i=0
                while i < len(sequence):
                        #genetate a tuple of pairwise residuals
                        pairwise = (sequence[i].upper(),this_seq[i].upper())
                        if "-" in pairwise:
                            print("processing gaps ...")
                        else:
                            if self.sub_matrix(62).get(pairwise) is not None:
                                p_score +=  self.sub_matrix(62).get(pairwise)
                                ref_seq +=  sequence[i].upper()
                                if pairwise[0]== pairwise[1]:
                                        p_identity +=1
                            else:
                                #print(pairwise)
                                p_score += self.sub_matrix(62).get(Builtins().rev_tuple(pairwise))
                                ref_seq +=  sequence[i].upper()
                                if pairwise[0]== pairwise[1]:
                                        p_identity +=1
                        #end the loop
                        i += 1
            #push the pairwise into the list for scores
                scores.append(p_score)
                pids.append(p_identity/len(ref_seq)*100)  
                #print(f"Here is are the ids: {pids}") 
                #print(f"The reference sequence is: {ref_seq}")       
            #replace the deleted self.sequences
            self.sequences.insert(count,sequence)          
            count += 1
        scaled_dwn_scores = list(map(lambda x: x/1,scores))  
        results.update({"Seq Count": len(self.sequences)})   
        results.update({"Sim Score": np.mean(scaled_dwn_scores)})
        results.update({"identities": np.mean(pids)})
        results.update({"scores": scores})
        return results
    
    def sub_matrix(self,percent):
        if percent==50:
            return MatrixInfo.blosum50
        else:
            return MatrixInfo.blosum62
    def entropy_score(self):
        alignment_length = len(self.sequences[0])
        entropy_values = []

        for i in range(alignment_length):
            column = [sequence[i] for sequence in self.sequences if sequence[i] != '-']
            num_residues = len(column)
            residue_counts = {}
            
            for residue in column:
                if residue in residue_counts:
                    residue_counts[residue] += 1
                else:
                    residue_counts[residue] = 1
            
            entropy = 0
            for count in residue_counts.values():
                frequency = count / num_residues
                entropy -= frequency * math.log2(frequency)

            entropy_values.append(entropy)

        return (sum(entropy_values)/len(entropy_values)) #return mean antropy
    from Bio import SeqIO
    #An MSA with higher gap content value has decreased quality
    def gap_content(self):
        alignment_length = len(self.sequences[0])
        gap_counts = [0] * alignment_length

        for i in range(alignment_length):
            gap_counts[i] = sum(1 for sequence in self.sequences if sequence[i] == '-')

        gap_percentage = [count / len(self.sequences) for count in gap_counts]

        return (sum(gap_percentage)/len(gap_percentage)) #mean gap content
    
    def sum_of_pairs_score(self):
        sp_score = 0
        alignment_length = len(self.sequences[0])  # Assuming all sequences have the same length
        total_residuals = 0
        for i in range(alignment_length):
            for j in range(len(self.sequences)):
                residue_i = self.sequences[j][i]
                for k in range(j + 1, len(self.sequences)):
                    residue_k = self.sequences[k][i]
                    if residue_i != '-' and residue_k != '-':
                        total_residuals += 1
                        if residue_i == residue_k:
                            sp_score += 1

        return (sp_score/total_residuals)
    
    def calculate_identity_score(msa_file): #pass the file here
        alignment = AlignIO.read(msa_file, "fasta")
        num_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()

        identical_count = 0
        for i in range(alignment_length):
            column = alignment[:, i]
            if len(set(column)) == 1:
                identical_count += 1

        identity_score = identical_count / alignment_length * 100

        return identity_score
    @staticmethod
    def sort_results():        
      df = pd.read_csv(Config.settings()["results_path"])
      df.sort_values(by=["Similarity Score"],ascending=False,inplace=True)
      print(df)

    @staticmethod
    def generate_stats(data):
        print("Results of Statistical data analysis")
        print("============================================================================================")
        #reshape the data frame for statmodels package
        #df_melt = pd.melt(data.reset_index(), id_vars=['index'], value_vars=Config.settings()["programs"])
        # replace column names
        #df_melt.columns = ['index', 'programs', 'blosum_scores']
        #model = ols('blosum_scores ~ C(programs)', data=df_melt).fit()
        #anova_table = sm.stats.anova_lm(model, typ=2)
        #p_value = float(anova_table["PR(>F)"][0])
        #anova with biostat
        #res = biostat()
        #res.anova_stat(df=df_melt, res_var='value', anova_model='blosum_scores ~ C(program)')
        #res.anova_summary
        #res.tukey_hsd(df=df_melt, res_var='blosum_scores', xfac_var='programs', anova_model='blosum_scores ~ C(programs)')
        #print(res.tukey_summary)
        #print(anova_table)
        #from statsmodels.stats.multicomp import pairwise_tukeyhsd

        # perform multiple pairwise comparison (Tukey HSD)
        #m_comp = pairwise_tukeyhsd(endog=df_melt['blosum_scores'], groups=df_melt['programs'], alpha=0.05)
        #print(m_comp)
        print("Do a Kruskal Wallis Non parametric test")
        args = [data[program] for program in data]
        sig,pvalue = stats.kruskal(*args)
        print(f"Kruskal Wallis pvalue is: {pvalue}")
        if pvalue < 0.05:
            print("Pvalue is less than 0.05. This means the quality of alignment varies significantly for the programs")
        #genarate some graphs
        #ax = sns.violinplot(x='programs', y='blosum_scores', data=df_melt)
        #ax = sns.swarmplot(x="programs", y="blosum_scores", data=df_melt, size = 20,color='#7d0013')
        #plt.show()
        




                   