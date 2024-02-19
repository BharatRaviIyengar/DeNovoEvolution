import os
os.chdir("path")
import random
from Bio import SeqIO


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def make_pete_graph(Orthogroups):
    f_pete = open("data_for_pete_model","w")
    f_anna_ortho = open("data_rst_anna_ortho", "w")
    f_anna_ORF = open("data_rst_anna_MISSING_ORF", "w")
    nb_ortho_strand_0 = 0
    nb_ORF_lost_strand_0 = 0
    nb_ortho_strand_0_total = 0
    nb_ORF_expected_strand_0 = 0
    
    nb_ortho_strand_1 = 0
    nb_ORF_lost_strand_1 = 0
    nb_ortho_strand_1_total = 0
    nb_ORF_expected_strand_1 = 0
    
    nb_ortho_strand_2 = 0
    nb_ORF_lost_strand_2 = 0
    nb_ortho_strand_2_total = 0
    nb_ORF_expected_strand_2 = 0
    
    strand_0_7pop = 0
    strand_0_6pop = 0
    strand_0_5pop = 0
    strand_0_4pop = 0
    strand_0_3pop = 0
    strand_0_2pop = 0
    strand_0_1pop = 0
    strand_1_7pop = 0
    strand_1_6pop = 0
    strand_1_5pop = 0
    strand_1_4pop = 0
    strand_1_3pop = 0
    strand_1_2pop = 0
    strand_1_1pop = 0
    strand_2_7pop = 0
    strand_2_6pop = 0
    strand_2_5pop = 0
    strand_2_4pop = 0
    strand_2_3pop = 0
    strand_2_2pop = 0
    strand_2_1pop = 0
    for line in Orthogroups:
        dico_line = {}
        proper_line = line.split("\n")[0]
        list_elts = proper_line.split(",")[1:]
        pos = 0
        list_strands = []
        while pos < len(list_elts):
            pop = list_elts[pos].split("_")[0]
            strand = list_elts[pos+1]
            if strand not in list_strands:
                list_strands.append(strand)
            dico_line[pop] = strand
            pos+=2
        print (line)
        print (dico_line)
        print ("************")
        if len(list_strands) == 1:
            if strand == "0":
                if len(dico_line) == 7:
                    strand_0_7pop += 1
                elif len(dico_line) == 6:
                    strand_0_6pop += 1
                elif len(dico_line) == 5:
                    strand_0_5pop += 1
                elif len(dico_line) == 4:
                    strand_0_4pop += 1
                elif len(dico_line) == 3:
                    strand_0_3pop += 1
                elif len(dico_line) == 2:
                    strand_0_2pop += 1
                elif len(dico_line) == 1:
                    strand_0_1pop += 1
                if len(dico_line)>1 and "ZB" in dico_line.keys():
                    if len(dico_line)<7: 
                        nb_ortho_strand_0 += 1
                        nb_ORF_lost_strand_0 += (7-len(dico_line))
                    nb_ortho_strand_0_total += 1
                    nb_ORF_expected_strand_0 += 7
            elif strand == "1":
                if len(dico_line) == 7:
                    strand_1_7pop += 1
                elif len(dico_line) == 6:
                    strand_1_6pop += 1
                elif len(dico_line) == 5:
                    strand_1_5pop += 1
                elif len(dico_line) == 4:
                    strand_1_4pop += 1
                elif len(dico_line) == 3:
                    strand_1_3pop += 1
                elif len(dico_line) == 2:
                    strand_1_2pop += 1
                elif len(dico_line) == 1:
                    strand_1_1pop += 1
                if len(dico_line)>1 and "ZB" in dico_line.keys():
                    if len(dico_line)<7: 
                        nb_ortho_strand_1 += 1
                        nb_ORF_lost_strand_1 += (7-len(dico_line))
                    nb_ortho_strand_1_total += 1
                    nb_ORF_expected_strand_1 += 7
            elif strand == "2":
                if len(dico_line) == 7:
                    strand_2_7pop += 1
                elif len(dico_line) == 6:
                    strand_2_6pop += 1
                elif len(dico_line) == 5:
                    strand_2_5pop += 1
                elif len(dico_line) == 4:
                    strand_2_4pop += 1    
                elif len(dico_line) == 3:
                    strand_2_3pop += 1
                elif len(dico_line) == 2:
                    strand_2_2pop += 1
                elif len(dico_line) == 1:
                    strand_2_1pop += 1
                if len(dico_line)>1 and "ZB" in dico_line.keys():
                    if len(dico_line)<7: 
                        nb_ortho_strand_2 += 1
                        nb_ORF_lost_strand_2 += (7-len(dico_line))
                    nb_ortho_strand_2_total += 1
                    nb_ORF_expected_strand_2 += 7
    f_pete.write("strand"+","+"nb_line"+","+"nb_ORFs"+"\n")
    f_pete.write("0"+","+"1"+","+str(strand_0_1pop)+"\n")
    f_pete.write("0"+","+"2"+","+str(strand_0_2pop)+"\n")
    f_pete.write("0"+","+"3"+","+str(strand_0_3pop)+"\n")
    f_pete.write("0"+","+"4"+","+str(strand_0_4pop)+"\n")
    f_pete.write("0"+","+"5"+","+str(strand_0_5pop)+"\n")
    f_pete.write("0"+","+"6"+","+str(strand_0_6pop)+"\n")
    f_pete.write("0"+","+"7"+","+str(strand_0_7pop)+"\n")
    f_pete.write("1"+","+"1"+","+str(strand_1_1pop)+"\n")
    f_pete.write("1"+","+"2"+","+str(strand_1_2pop)+"\n")
    f_pete.write("1"+","+"3"+","+str(strand_1_3pop)+"\n")
    f_pete.write("1"+","+"4"+","+str(strand_1_4pop)+"\n")
    f_pete.write("1"+","+"5"+","+str(strand_1_5pop)+"\n")
    f_pete.write("1"+","+"6"+","+str(strand_1_6pop)+"\n")
    f_pete.write("1"+","+"7"+","+str(strand_1_7pop)+"\n")
    f_pete.write("2"+","+"1"+","+str(strand_2_1pop)+"\n")
    f_pete.write("2"+","+"2"+","+str(strand_2_2pop)+"\n")
    f_pete.write("2"+","+"3"+","+str(strand_2_3pop)+"\n")
    f_pete.write("2"+","+"4"+","+str(strand_2_4pop)+"\n")
    f_pete.write("2"+","+"5"+","+str(strand_2_5pop)+"\n")
    f_pete.write("2"+","+"6"+","+str(strand_2_6pop)+"\n")
    f_pete.write("2"+","+"7"+","+str(strand_2_7pop)+"\n")
    
    percOrthoLost_F0 = float(100)*float(nb_ortho_strand_0)/float(nb_ortho_strand_0_total)
    percORFlost_F0 = float(100)*float(nb_ORF_lost_strand_0)/float(nb_ORF_expected_strand_0)
    
    percOrthoLost_F1 = float(100)*float(nb_ortho_strand_1)/float(nb_ortho_strand_1_total)
    percORFlost_F1 = float(100)*float(nb_ORF_lost_strand_1)/float(nb_ORF_expected_strand_1)

    percOrthoLost_F2 = float(100)*float(nb_ortho_strand_2)/float(nb_ortho_strand_2_total)
    percORFlost_F2 = float(100)*float(nb_ORF_lost_strand_2)/float(nb_ORF_expected_strand_2)
    
    f_anna_ortho.write("strand,nbOrtho,percOrtho"+"\n")
    f_anna_ortho.write("F0,"+str(nb_ortho_strand_0)+","+str(percOrthoLost_F0)+"\n")
    f_anna_ortho.write("F1,"+str(nb_ortho_strand_1)+","+str(percOrthoLost_F1)+"\n")
    f_anna_ortho.write("F2,"+str(nb_ortho_strand_2)+","+str(percOrthoLost_F2)+"\n")
    
    f_anna_ORF.write("strand,nbMissingORFs,percMissing"+"\n")
    f_anna_ORF.write("F0,"+str(nb_ORF_lost_strand_0)+","+str(percORFlost_F0)+"\n")
    f_anna_ORF.write("F1,"+str(nb_ORF_lost_strand_1)+","+str(percORFlost_F1)+"\n")
    f_anna_ORF.write("F2,"+str(nb_ORF_lost_strand_2)+","+str(percORFlost_F2)+"\n")
    
    f_pete.close()
    f_anna_ortho.close()
    f_anna_ORF.close()


            

Orthogroups = openFile("orthogroups_correct.txt")
make_pete_graph(Orthogroups)




