import os
os.chdir("path")
import random
from Bio import SeqIO


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def detectdoublons(myfile):
    dic = {}
    l_doublons = []
    for line in myfile:
        mylist = line.split("\n")[0].split(",")
        name_orthogroup = mylist[0]
        for elt in mylist[1:]:
            if elt not in dic.keys():
                dic[elt] = name_orthogroup
            else:
                if name_orthogroup not in l_doublons:
                    l_doublons.append(name_orthogroup)
                if dic[elt] not in l_doublons:
                    l_doublons.append(dic[elt])
    return l_doublons


def extract_orthogroups(myfile):
    d = {}
    for line in myfile:
        mylist = line.split("\n")[0].split(",")
        name_orthogroup = mylist[0]
        new_list = []
        for elt in mylist[1:]:
            new_list.append(elt)
        d[name_orthogroup] = new_list
    return d


def redo_orthogroups(doublons, dico):
    subdico = {}
    final_joined = []
    print (len(dico))
    for orthogroups in dico.keys():
        if orthogroups in doublons:
            subdico[orthogroups] = dico[orthogroups]
    for orthogroup in subdico.keys():
        present = False
        for k in final_joined:
            if orthogroup in k:
                present = True
        if present == False:
            list_seq = subdico[orthogroup]
            list_similar_ortho = [orthogroup]
            for seq in list_seq:
                for other_ortho in subdico.keys():
                    if seq in subdico[other_ortho]:
                        if other_ortho not in list_similar_ortho:
                            list_similar_ortho.append(other_ortho)
            final_joined.append(list_similar_ortho)
    new_sub_dico = {}
    for joined_groups in final_joined:
        liste_new_seqs = []
        for ortho in joined_groups:
            for seq in dico[ortho]:
                if seq not in liste_new_seqs:
                    liste_new_seqs.append(seq)
        new_sub_dico[joined_groups[0]] = liste_new_seqs
    for joined_groups in final_joined:
        for ortho in joined_groups:
            if ortho in dico.keys():
                del dico[ortho]
    print (len(new_sub_dico))
    print (len(dico))
    for new_ortho in new_sub_dico.keys():
        dico[new_ortho] = new_sub_dico[new_ortho]
    print (len(dico))
    #print (final_joined)
    #c = 0
    #l = []
    #d = ""
    #for elt in final_joined:
        #for j in elt:
            #c+=1
            #if j not in l:
                #l.append(j)
            #else:
                #d = j
    #print (c)
    #print (d)
    #for joined_groups in final_joined:
        #dico_lala = {}
        #for ortho in joined_groups:
            #for seq in dico[ortho]:
                #if seq in dico_lala.keys():
                    #dico_lala[seq]+=1
                #else:
                    #dico_lala[seq]=1
        #nbseq = 0
        #nbseqalone = 0
        #print (dico_lala)
        #for seqs in dico_lala.keys():
            #if dico_lala[seq] == 1:
                #nbseqalone+=1
            #nbseq+=1
        #print (joined_groups)
        #print (nbseq)
        #print (nbseqalone)
        #print ("*********")
   
        
def sort_list(my_list):
    new_list = []
    #print (my_list)
    if len(my_list)<2:
        new_list = my_list
    else:
        my_min = 1000
        seq = ""
        while len(my_list)>0:
            for my_seq in my_list:
                nb = int(my_seq.split("_")[3])
                if nb < my_min:
                    my_min = nb
                    seq = my_seq
            new_list.append(seq)
            my_list.remove(seq)
            my_min = 1000
            seq = ""
    #print (new_list)
    #print ("***********")
    return new_list


#def shuffle_combs(all_dico):
    #matrix = []
    #while len(all_dico)>0:
        #num_min = 1000
        #name_line = ""
        #for i in all_dico.keys():
            #if len(all_dico[i])<num_min:
                #name_line = i
                #num_min = len(all_dico[i])
        #new_l = []
        #for index in all_dico[name_line]:
            #new_l.append({index:all_dico[name_line][index]})
        #matrix.append(new_l)
        #Max_comb = len(new_l)
        
        
    #print (matrix)
                
    
def algo_attrib(dico_filled,dico_strand):
    new_ortho = []
    nb_lines = len(dico_filled)
    taille = 0
    equal = True
    for pop in dico_filled.keys():
        taille = len(dico_filled[pop])
        break
    for pop in dico_filled.keys():
        if len(dico_filled[pop])!= taille:
            equal = False
            break   
    if equal == True:
        #print (dico_filled)
        for i in range(taille):
            sublist = []
            for line in dico_filled.keys():
                orf = dico_filled[line][i]
                sublist.append(orf)
            new_ortho.append(sublist)
        #print (new_ortho)
    else:
        all_dico = {}
        for line in dico_filled.keys():
            list_orf = dico_filled[line]
            new_dico = {}
            c = 1
            for orf in list_orf:
                new_dico[line+"_"+str(c)] = dico_strand[orf]
                c+=1
            all_dico["D_"+line] = new_dico
        list_frames = []
        for pop in all_dico.keys():
            for pos in all_dico[pop].keys():
                fr = all_dico[pop][pos]
                if fr not in list_frames:
                    list_frames.append(fr)
        if len(list_frames) == 1:
            #print (all_dico)
            taille_max = 0
            for i in all_dico.keys():
                if len(all_dico[i])>taille_max:
                    taille_max = len(all_dico[i])
            for i in range(taille_max):
                sublist = []
                for line in dico_filled.keys():
                    if len(dico_filled[line]) >= (i+1):
                        orf = dico_filled[line][i]
                        sublist.append(orf)
                new_ortho.append(sublist)
            #print (all_dico)
            #print (dico_filled)
            #print (new_ortho)
            #print ("**************")
        #else:
            #print (all_dico)
            #print ("lqlq")
            #continue
    return new_ortho


    #fake_dico = {"D_AK":{"AK_1":"1","AK_2":"0","AK_3":"2","AK_4":"1","AK_5":"0","AK_6":"0","AK_7":"1"},
                 #"D_DK":{"DK_1":"2","DK_2":"1","DK_3":"0","DK_4":"2","DK_5":"1"},
                 #"D_GI":{"GI_1":"1","GI_2":"0","GI_3":"0","GI_4":"2","GI_5":"1","GI_6":"2"},
                 #"D_SW":{"SW_1":"1","SW_2":"1","SW_3":"0","SW_4":"2","SW_5":"1","SW_6":"1"},
                 #"D_UM":{"UM_1":"0","UM_2":"1","UM_3":"2","UM_4":"1"},
                 #"D_YE":{"YE_1":"0","YE_2":"1","YE_3":"2","YE_4":"1","YE_5":"2"},
                 #"D_ZB":{"ZB_1":"1","ZB_2":"2","ZB_3":"0","ZB_4":"1","ZB_5":"1","ZB_6":"0"}}
    #shuffle_combs(fake_dico)
    

def split_the_ortho(list_AK,list_DK,list_GI,list_SW,list_UM,list_YE,list_ZB,dico_strand,dico_orthogroups,orthogroup):
    new_ortho = []
    prop = ""
    c = 1
    if len(list_AK)<2 and len(list_DK)<2 and len(list_GI)<2 and len(list_SW)<2 and len(list_UM)<2 and len(list_YE)<2 and len(list_ZB)<2:
        new_ortho = [dico_orthogroups[orthogroup]]
        prop = "no_initial_dupli"
    else:
        sorted_AK = sort_list(list_AK)
        sorted_DK = sort_list(list_DK)
        sorted_SW = sort_list(list_SW)
        sorted_GI = sort_list(list_GI)
        sorted_UM = sort_list(list_UM)
        sorted_YE = sort_list(list_YE)
        sorted_ZB = sort_list(list_ZB)
        nb_list_filed = 0
        dico_filled = {}
        len_AK = len(sorted_AK)
        len_DK = len(sorted_DK)
        len_GI = len(sorted_GI)
        len_SW = len(sorted_SW)
        len_UM = len(sorted_UM)
        len_YE = len(sorted_YE)
        len_ZB = len(sorted_ZB)
        if len_AK > 0:
            nb_list_filed += 1
            dico_filled["AK"] = sorted_AK
        if len_DK > 0:
            nb_list_filed += 1
            dico_filled["DK"] = sorted_DK
        if len_GI > 0:
            nb_list_filed += 1
            dico_filled["GI"] = sorted_GI
        if len_SW > 0:
            nb_list_filed += 1
            dico_filled["SW"] = sorted_SW
        if len_UM > 0:
            nb_list_filed += 1
            dico_filled["UM"] = sorted_UM
        if len_YE > 0:
            nb_list_filed += 1
            dico_filled["YE"] = sorted_YE
        if len_ZB > 0:
            nb_list_filed += 1
            dico_filled["ZB"] = sorted_ZB
        if nb_list_filed == 1:
            #print (dico_filled)
            for pop in dico_filled.keys():
                for seqs in dico_filled[pop]:
                    new_ortho.append([seqs])
            #print (new_ortho)
        else:
            new_ortho = algo_attrib(dico_filled,dico_strand)
    return new_ortho,prop

        
def split_correct_orthogroups(dico_orthogroups,dico_strand):
    file_correct = open("orthogroups_correct.txt", "w")
    file_uncorrect = open("orthogroups_incorrect.txt", "w")
    nb_same_frame = 0
    nb_diff_frame = 0
    nb_no_dup = 0
    nb_single = 0
    no_initial_dupli = 0
    nbb = 0
    nbc = 0
    compteur_new_orthogroups = 1
    for orthogroup in dico_orthogroups.keys():
        list_AK = []
        list_DK = []
        list_GI = []
        list_SW = []
        list_UM = []
        list_YE = []
        list_ZB = []
        for seq in dico_orthogroups[orthogroup]:
            pop = seq.split("_")[0]
            if pop == "AK":
                list_AK.append(seq)
            elif pop == "DK":
                list_DK.append(seq)
            elif pop == "GI":
                list_GI.append(seq)
            elif pop == "SW":
                list_SW.append(seq)
            elif pop == "UM":
                list_UM.append(seq)
            elif pop == "YE":
                list_YE.append(seq)
            elif pop == "ZB":
                list_ZB.append(seq)
        new_ortho,prop = split_the_ortho(list_AK,list_DK,list_GI,list_SW,list_UM,list_YE,list_ZB,dico_strand,dico_orthogroups,orthogroup)
        if prop == "no_initial_dupli":
            no_initial_dupli += 1
        #print (new_ortho)
        if len(new_ortho) == 0:
            for i in dico_orthogroups[orthogroup]:
                file_uncorrect.write(orthogroup)
                file_uncorrect.write(","+i)
            file_uncorrect.write("\n")
        else:
            for new_orthogroup in new_ortho:
                list_AK = []
                list_DK = []
                list_GI = []
                list_SW = []
                list_UM = []
                list_YE = []
                list_ZB = []
                print (new_orthogroup)
                for seq in new_orthogroup:
                    pop = seq.split("_")[0]
                    if pop == "AK":
                        list_AK.append(seq)
                    elif pop == "DK":
                        list_DK.append(seq)
                    elif pop == "GI":
                        list_GI.append(seq)
                    elif pop == "SW":
                        list_SW.append(seq)
                    elif pop == "UM":
                        list_UM.append(seq)
                    elif pop == "YE":
                        list_YE.append(seq)
                    elif pop == "ZB":
                        list_ZB.append(seq)
                list_strand = []
                nb_ligne = 0
                if len(list_AK)>0:
                    strand = dico_strand[list_AK[0]]
                    if strand not in list_strand:
                        list_strand.append(strand)
                    nb_ligne+=1
                if len(list_DK)>0:
                    strand = dico_strand[list_DK[0]]
                    if strand not in list_strand:
                        list_strand.append(strand)
                    nb_ligne+=1
                if len(list_GI)>0:
                    strand = dico_strand[list_GI[0]]
                    if strand not in list_strand:
                        list_strand.append(strand)
                    nb_ligne+=1
                if len(list_SW)>0:
                    strand = dico_strand[list_SW[0]]
                    if strand not in list_strand:
                        list_strand.append(strand)
                    nb_ligne+=1
                if len(list_YE)>0:
                    strand = dico_strand[list_YE[0]]
                    if strand not in list_strand:
                        list_strand.append(strand)
                    nb_ligne+=1
                if len(list_UM)>0:
                    strand = dico_strand[list_UM[0]]
                    if strand not in list_strand:
                        list_strand.append(strand)
                    nb_ligne+=1
                if len(list_ZB)>0:
                    strand = dico_strand[list_ZB[0]]
                    if strand not in list_strand:
                        list_strand.append(strand)
                    nb_ligne+=1
                if nb_ligne>1:
                    if len(list_strand) > 1:
                        nb_diff_frame+=1
                    else:
                        nb_same_frame+=1
                else:
                    nb_single +=1
                file_correct.write("orthogroup_"+str(compteur_new_orthogroups))
                for seq in new_orthogroup:
                    file_correct.write(","+seq+","+dico_strand[seq])
                file_correct.write("\n")
                compteur_new_orthogroups+=1
            
    #print ("nb total no dup : " + str(nb_no_dup))
    print (no_initial_dupli)
    print ("nb same only one seq : " + str(nb_single))
    print ("nb diff frame : " + str(nb_diff_frame))
    print ("nb same frame : " + str(nb_same_frame))

        
    def make_dico_strand(my_file):
    dico_strand = {}
    for ligne in my_file:
        nameORF = ligne.split(",")[0]
        strandORF = ligne.split(",")[8]
        dico_strand[nameORF] = strandORF
    #print (dico_strand)
    return dico_strand

orthogroups_file = openFile("final_orthogroups.txt")
l_doublon = detectdoublons(orthogroups_file)
#print (len(l_doublon))
dico_orthogroups = extract_orthogroups(orthogroups_file)
redo_orthogroups(l_doublon, dico_orthogroups)
file_all_data = openFile("All_datafile")
dico_strand = make_dico_strand(file_all_data)
split_correct_orthogroups(dico_orthogroups,dico_strand)





