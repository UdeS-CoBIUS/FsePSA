#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``fse_main.py`` **module description**:

This module allows to run a pairwise alignment algorithm on every pair of 
a set of coding sequences given in a multifasta file and write the resulting 
alignments in a srspair alignment file. The alignment algorithm is described 
in manuscript:
    * F. Bélanger, A. Rachati, A. Ouangraoua. Aligning protein-coding sequences 
with frameshift extension penalties. (2016).

.. moduleauthor:: François Bélanger and Aida Ouangraoua

February 2016

"""

import argparse
from Bio import SeqIO
import time

from scoring_matrix import *
from fse import fse

LENGTH_ALN_LINE = 50

def build_arg_parser():
    parser = argparse.ArgumentParser(description='Align sequence')
    parser.add_argument('-go', '--gapopen', type=float, default=-10)
    parser.add_argument('-ge', '--gapextend', type=float, default=-0.2)
    parser.add_argument('-fso', '--fsopen', type=float, default=-10)
    parser.add_argument('-fse', '--fsextend', type=float, default=-0.2)
    parser.add_argument('-aa', '--aminoacidmatrix', default='resources/BLOSUM62.txt')
    parser.add_argument('-d', '--datasequence', nargs=2, default=['examples/example_data.fasta', 'fasta'])
    parser.add_argument('-o', '--outfile', default='examples/result_example_data.srspair')
    parser.add_argument('-of', '--outformat', default='srspair')
    return parser

def compute_fslength(sequence1,sequence2):
    total_fsopen = total_fsextend = 0
    sequence1_ = sequence1.tomutable()
    sequence2_ = sequence2.tomutable()
    position1 = 0
    position2 = 0
    inFrameshift = False
    for i in range(len(sequence1)):
        if sequence1[i] != "-":
            position1 += 1
        if sequence2[i] != "-":
            position2 += 1
        if(sequence1[i] != "-" and sequence2[i] != "-"):
            if(inFrameshift == True):
                total_fsextend += 1
            elif(((position1 % 3 != 0 and  position2 % 3 == 0) or
            (position1 % 3 == 0 and  position2 % 3 != 0)) and
            sequence1[i-2] != "-" and sequence1[i-1] != "-" and
            sequence2[i-2] != "-" and sequence2[i-1] != "-"):
                inFrameshift = True
                total_fsopen += 1
                total_fsextend += 3
                if (sequence1[i-3] == "-"):
                    sequence1_[i-3] = "!"
                if (sequence2[i-3] == "-"):
                    sequence2_[i-3] = "!"
        else:
            inFrameshift = False                
    return total_fsopen, total_fsextend, sequence1_, sequence2_

def compute_alignment_composition(sequence1,sequence2):
    nb_identity = 0
    nb_gap1 = 0
    nb_gap2 = 0
    markup_line = ""
    for i in range(len(sequence1)):
        if sequence1[i] == "-" or sequence1[i] == "!":
            nb_gap1 += 1
            if sequence1[i] == "-":
                markup_line += " "
            if sequence1[i] == "!":
                markup_line += "!"
                sequence1[i] = "-"
        elif sequence2[i] == "-" or  sequence2[i] == "!":
            nb_gap2 += 1
            if sequence2[i] == "-":
                markup_line += " "
            if sequence2[i] == "!":
                markup_line += "!"
                sequence2[i] = "-"
        elif sequence1[i] == sequence2[i]:
            nb_identity += 1
            markup_line += "|"
        else:
            markup_line += "."
    return nb_identity, nb_gap1, nb_gap2, markup_line,sequence1,sequence2

def print_alignment_header(sequence1_id,sequence2_id,arg, outputfile):
    outputfile.write("\n")
    outputfile.write("#=======================================\n")
    outputfile.write("#\n")
    outputfile.write("# Aligned_sequences: 2\n")
    outputfile.write("# 1: " + sequence1_id + "\n")
    outputfile.write("# 2: " + sequence2_id + "\n")
    outputfile.write("# Matrix: " + arg.aminoacidmatrix + "\n")
    outputfile.write("# GapOpen_penalty: " + str(arg.gapopen) + "\n")
    outputfile.write("# GapExtend_penalty: " + str(arg.gapextend) + "\n")
    outputfile.write("# FSopen_penalty: " + str(arg.fsopen) + "\n")
    outputfile.write("# FSextend_penalty: " + str(arg.fsextend) + "\n")
    outputfile.write("#\n")

def print_alignment(total_length,total_identity,total_gap,total_fsopen,total_fsextend,alignment,score,outputfile):
    outputfile.write("# Length: " + str(total_length) + "\n")
    outputfile.write("# Identity:       " + str(total_identity)+"/"+str(total_length) + " " + str(round(100.0*total_identity/total_length,1)) + "%\n")
    outputfile.write("# Similarity:     " + str(total_identity)+"/"+str(total_length) + " " + str(round(100.0*total_identity/total_length,1)) + "%\n")
    outputfile.write("# Gaps:           " + str(total_gap)+"/"+str(total_length) + " " + str(round(100.0*total_gap/total_length,1)) + "%\n")
    outputfile.write("# FSopen:         " + str(total_fsopen)+"/"+str(total_length) + " " + str(round(100.0*total_fsopen/total_length,1)) + "%\n")
    outputfile.write("# FrameShift:     " + str(total_fsextend)+"/"+str(total_length) + " " + str(round(100.0*total_fsextend/total_length,1)) + "%\n")    
    outputfile.write("# Score: " + str(score) + "\n")
    outputfile.write("#\n")
    outputfile.write("#\n")
    outputfile.write("#=======================================\n")
    outputfile.write("\n")
    outputfile.write(str(alignment))    
    outputfile.write("#---------------------------------------\n")
    outputfile.write("#---------------------------------------\n")

def format_alignment(sequence1_id,sequence2_id,sequence1,sequence2,outformat):
    total_identity = total_fsopen = total_fsextend = total_gap1 = total_gap2 = 0
    alignment = ""
    aln_srspair = ""
    aln1_fasta = aln2_fasta = ""
    start = 0
    while start < len(sequence1):
        line_length = LENGTH_ALN_LINE
        if(start+line_length > len(sequence1)):
            line_length = len(sequence1)-start
        subsequence1 = sequence1[start:start+line_length]
        subsequence2 = sequence2[start:start+line_length]
        aln1_fasta += subsequence1 + "\n"
        aln2_fasta += subsequence2 + "\n"
        nb_identity,nb_gap1,nb_gap2,markup_line,subsequence1,subsequence2 = compute_alignment_composition(subsequence1,subsequence2)

        prefix1 = sequence1_id
        for i in range (len(sequence1_id),20):
            prefix1 += " "
        end_prefix1 = " "+ str(start-total_gap1+1) + " "
        prefix1 = prefix1[:20-len(end_prefix1)] + end_prefix1
            
        prefix2 = sequence2_id
        for i in range (len(sequence2_id),20):
            prefix2 += " "
        end_prefix2 = " "+ str(start-total_gap2+1) + " "
        prefix2 = prefix2[:20-len(end_prefix2)] + end_prefix2
            
        aln_srspair +=  prefix1 + subsequence1 + "     " + str(start+line_length - total_gap1 - nb_gap1) + "\n"
        prefix_markup_line = ""
        for i in range (20):
            prefix_markup_line += " "
        aln_srspair +=  prefix_markup_line + markup_line + "\n"
        aln_srspair += prefix2 +subsequence2 + "     " + str(start+line_length - total_gap2 - nb_gap2) + "\n\n"

        total_identity += nb_identity
        total_gap1 += nb_gap1
        total_gap2 += nb_gap2
        start += line_length

    if(outformat == "srspair"):
        alignment = aln_srspair
    if(outformat == "fasta"):
        alignment = ">" + sequence1_id + "\n" + aln1_fasta + "\n" + ">" + sequence2_id + "\n" + aln2_fasta

    total_length = len(sequence1)
    total_gap = total_gap1 + total_gap2
    return alignment,total_length,total_identity,total_gap
    
def print_result(sequence1_id,sequence2_id,score,sequence1, sequence2, arg, outputfile):
    total_fsopen, total_fsextend, sequence1, sequence2 = compute_fslength(sequence1,sequence2)
    alignment,total_length,total_identity,total_gap = format_alignment(sequence1_id,sequence2_id,sequence1,sequence2,arg.outformat)

    print_alignment_header(sequence1_id,sequence2_id,arg, outputfile)
    print_alignment(total_length,total_identity,total_gap,total_fsopen,total_fsextend,alignment,score,outputfile)
    
def print_file_header(arg, outputfile):
    outputfile.write("########################################\n")
    outputfile.write("# Program: fse\n")
    outputfile.write("# Rundate: " + time.strftime("%c") + "\n")
    outputfile.write("# Commandline: python fse_main.py\n")
    outputfile.write("#    --datasequence " + arg.datasequence[0] + " " + arg.datasequence[1] + "\n")
    outputfile.write("#    --gapopen " + str(arg.gapopen) + "\n")
    outputfile.write("#    --gapextend " + str(arg.gapextend) + "\n")
    outputfile.write("#    --fsopen " + str(arg.fsopen) + "\n")
    outputfile.write("#    --fsextend " + str(arg.fsextend) + "\n")    
    outputfile.write("#    --aminoacidmatrix " + arg.aminoacidmatrix + "\n")  
    outputfile.write("#    --outfile " + arg.outfile + "\n")
    outputfile.write("#    --outformat " + arg.outformat + "\n")
    outputfile.write("# Report_file: " + arg.outfile + "\n")
    outputfile.write("########################################\n")

def print_score_matrix(score_matrix,outputfile):
    outputfile.write("\n")
    outputfile.write("#=======================================\n")
    outputfile.write("#\n")
    outputfile.write("# Score_matrix\n")
    outputfile.write("#\n")
    outputfile.write("#=======================================\n")
    outputfile.write("\n")

    outputfile.write(score_matrix)

    outputfile.write("#---------------------------------------\n")
    outputfile.write("#---------------------------------------\n")

def main():
    parser = build_arg_parser()
    arg = parser.parse_args()

    outputfile = open(arg.outfile,"w")
    print_file_header(arg, outputfile)

    saa = ScoringMatrix(arg.aminoacidmatrix)
    saa.load()
    san = ScoringMatrix()
    san.init_similarity()

    seq_id_table = [] 
    seq_table = []
    seq_file = arg.datasequence

    for record in SeqIO.parse(seq_file[0], seq_file[1]) :
        seq_id_table.append(record.id)
        seq_table.append(record.seq)

    score_matrix = ""
    for i in range(len(seq_table)):
        score_matrix += seq_id_table[i] + "\t"
    score_matrix += "\n" 

    for i in range(len(seq_table)):
        for j in range(i+1):
            score_matrix += "\t"
        for j in range(i+1,len(seq_table)):                
            score, sequence1, sequence2 = fse(seq_table[i], seq_table[j], arg, saa, san)
            print_result(seq_id_table[i],seq_id_table[j],score,sequence1, sequence2, arg,outputfile)
            score_matrix += str(score) + "\t"
        score_matrix += "\n"
        
    print_score_matrix(score_matrix,outputfile)
    outputfile.close()

if __name__ == '__main__':
    main()
