import os
os.chdir("path")
import random
from Bio import SeqIO


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def make_dico_orfs(my_file):
    dico = {}
    for line in my_file:
        elts = line.split(",")
        name_ORF = elts[0]
        name_gene = elts[5]
        pos_in_gene = elts[8]
        dico[name_ORF] = [name_gene,pos_in_gene]
    return dico


def deal_blast(my_file,dico):
    final_list = []
    my_query = ""
    sublist = []
    for line in my_file:
        if line[0]!="#":
            elts = line.split()
            query = elts[0]
            target = elts[1]
            start_query = int(elts[3])
            end_query = int(elts[4])
            start_target = int(elts[5])
            end_target = int(elts[6])
            if query!= my_query:
                if my_query == "":
                    my_query = query
                    sublist = [query]
                else:
                    final_list.append(sublist)
                    sublist = [query]
                    my_query = query
            if query!=target:
                gene_query = dico[query][0]
                gene_target = dico[target][0]
                if gene_query == gene_target and start_query == start_target and end_query == end_target:
                    sublist.append(target)
    return final_list


def verify_list_similarity(l1,l2):
    same = False
    if len(l1) == len(l2):
        same = True
        for elt in l1:
            if elt in l2:
                continue
            else:
                same = False
                break
    return same
                  
    
def marge_orthogroups(List_AK,List_DK,List_GI,List_SW,List_UM,List_YE,List_ZB):
    final_orthogroups = []
    for sublist in List_AK:
        #if "AK_CDS_49_2_2" in sublist:
            #print ("lala")
        for orther_sublist in List_DK:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_DK.remove(orther_sublist)
                break
        for orther_sublist in List_GI:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_GI.remove(orther_sublist)
                break
        for orther_sublist in List_SW:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_SW.remove(orther_sublist)
                break
        for orther_sublist in List_UM:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_UM.remove(orther_sublist)
                break
        for orther_sublist in List_YE:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_YE.remove(orther_sublist)
                break
        for orther_sublist in List_ZB:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_ZB.remove(orther_sublist)
                break
        final_orthogroups.append(sublist)
        
    for sublist in List_DK:
        for orther_sublist in List_GI:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_GI.remove(orther_sublist)
                break
        for orther_sublist in List_SW:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_SW.remove(orther_sublist)
                break
        for orther_sublist in List_UM:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_UM.remove(orther_sublist)
                break
        for orther_sublist in List_YE:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_YE.remove(orther_sublist)
                break
        for orther_sublist in List_ZB:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_ZB.remove(orther_sublist)
                break
        final_orthogroups.append(sublist)
        
    for sublist in List_GI:
        for orther_sublist in List_SW:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_SW.remove(orther_sublist)
                break
        for orther_sublist in List_SW:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_SW.remove(orther_sublist)
                break
        for orther_sublist in List_UM:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_UM.remove(orther_sublist)
                break
        for orther_sublist in List_YE:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_YE.remove(orther_sublist)
                break
        for orther_sublist in List_ZB:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_ZB.remove(orther_sublist)
                break
        final_orthogroups.append(sublist)
        
    for sublist in List_SW:
        for orther_sublist in List_UM:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_UM.remove(orther_sublist)
                break
        for orther_sublist in List_YE:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_YE.remove(orther_sublist)
                break
        for orther_sublist in List_ZB:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_ZB.remove(orther_sublist)
                break
        final_orthogroups.append(sublist)
        
    for sublist in List_UM:
        for orther_sublist in List_YE:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_YE.remove(orther_sublist)
                break
        for orther_sublist in List_ZB:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_ZB.remove(orther_sublist)
                break
        final_orthogroups.append(sublist)
        
    for sublist in List_YE:
        for orther_sublist in List_ZB:
            similarity = verify_list_similarity(sublist,orther_sublist)
            if similarity == True:
                List_ZB.remove(orther_sublist)
                break
        final_orthogroups.append(sublist)
    for sublist in List_ZB:
        final_orthogroups.append(sublist)
    #for sublist in final_orthogroups:
        #if "AK_CDS_49_2_2" in sublist:
            #print ("lili")
    return final_orthogroups
    

def make_final_file(final_orthogroups):
    c = 1
    f = open("final_orthogroups.txt", "w")
    for sublist in final_orthogroups:
        f.write("ort_"+str(c))
        for orf in sublist:
            f.write(","+orf)
        f.write("\n")
        c+=1
    f.close()


def refine_orthogroups(orthogroups):
    new  =[]
    while len(orthogroups)>0:
        my_sublist = orthogroups[0]
        new.append(my_sublist)
        orthogroups.remove(my_sublist)
        for elt in orthogroups:
            verif = verify_list_similarity(my_sublist,elt)
            if verif == True:
                orthogroups.remove(elt)
    return (new)
                      

all_info = openFile("All_datafile")
dico_orf_gene = make_dico_orfs(all_info)

AK_output = openFile("AK_output.txt")
List_AK = deal_blast(AK_output,dico_orf_gene)

DK_output = openFile("DK_output.txt")
List_DK = deal_blast(DK_output,dico_orf_gene)

GI_output = openFile("GI_output.txt")
List_GI = deal_blast(GI_output,dico_orf_gene)

UM_output = openFile("UM_output.txt")
List_UM = deal_blast(UM_output,dico_orf_gene)

YE_output = openFile("YE_output.txt")
List_YE = deal_blast(YE_output,dico_orf_gene)

SW_output = openFile("SW_output.txt")
List_SW = deal_blast(SW_output,dico_orf_gene)

ZB_output = openFile("ZB_output.txt")
List_ZB = deal_blast(ZB_output,dico_orf_gene)

orthogroups_round1 = marge_orthogroups(List_AK,List_DK,List_GI,List_SW,List_UM,List_YE,List_ZB)
print (len(orthogroups_round1))
final_orthogroups = refine_orthogroups(orthogroups_round1)
print (len(final_orthogroups))
make_final_file(final_orthogroups)











