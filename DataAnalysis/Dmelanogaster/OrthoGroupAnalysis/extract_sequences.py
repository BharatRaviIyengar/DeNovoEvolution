import os
os.chdir("path")
import random
from Bio import SeqIO


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def store_gtf(my_file):
    list_final = []
    for lines in my_file:
        new_line = lines.split("\n")[0]
        l = new_line.split()
        start = int(l[3])-1
        end = int(l[4])
        chrom = l[0]
        transcript = l[8].split(";")[0]
        gene = l[9]
        pos_in_gene = l[7]
        sublist =[start, end, chrom, transcript, gene, pos_in_gene]
        list_final.append(sublist)
    return list_final


def store_prots(my_file):
    dico_prot = {}
    with open(my_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            name_seq = record.id.split("_")[0]+"_"+record.id.split("_")[1]
            #print (name_seq)
            seq = record.seq
            dico_prot[name_seq] = seq
    return dico_prot


def assess_several_orf_per_gene(data):
    new_data = []
    dico_gene = {}
    for sublist in data:
        gene_name = sublist[5]
        if gene_name not in dico_gene.keys():
            dico_gene[gene_name] = [int(sublist[1])]
        else:
            dico_gene[gene_name].append(int(sublist[1]))
    for gene in dico_gene.keys():
        if len(dico_gene[gene])>1:
            dico_gene[gene].sort()
    compteur = 1
    gene_attributed_number = {}
    for sublist in data:
        gene_name = sublist[5]
        if len(dico_gene[gene_name])>1:
            if gene_name in gene_attributed_number.keys():
                number_ORF = gene_attributed_number[gene_name]
            else:
                gene_attributed_number[gene_name] = compteur
                number_ORF = gene_attributed_number[gene_name]
                compteur+=1
                
            start = int(sublist[1])
            total = len(dico_gene[gene_name])
            pos = 1
            for coord in dico_gene[gene_name]:
                if coord == start:
                    break
                pos+=1
            new_orf_name = sublist[0]+str(number_ORF)+"_"+str(pos)+"_"+str(total)
            new_sublist = [new_orf_name]
            for elt in sublist[1:]:
                new_sublist.append(elt)
            new_data.append(new_sublist)
        else:
            new_orf_name = sublist[0]+str(compteur)
            new_sublist = [new_orf_name]
            for elt in sublist[1:]:
                new_sublist.append(elt)
            new_data.append(new_sublist)
            compteur+=1
    return new_data


def get_all_info(fasta, list_gtf,dico_prot,name_pop,dico_gtf):
    #compteur = 1
    final_list = []
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            name_seq = record.id
            prot_seq = dico_prot[name_seq]
            seq = record.seq
            l_name_seq = name_seq.split(":")
            end_seq = l_name_seq[3].split("(")
            start = int(end_seq[0].split("-")[0])
            end = int(end_seq[0].split("-")[1])
            sign = (end_seq[1].split(")")[0])
            chrom = l_name_seq[2]
            name_CDS = name_pop+"_"+"CDS_"#+str(compteur)
            trouve = False
            for sublists in list_gtf:
                if sublists[0] == start and sublists[1] == end and sublists[2] == chrom:
                    name_transcript = sublists[4]
                    name_gene = dico_gtf[name_transcript]
                    liste_my_seq = [name_CDS,start,end,chrom,sign,name_gene,sublists[3],sublists[4],sublists[5],seq,prot_seq]
                    final_list.append(liste_my_seq)
                    break
            #compteur+=1
    final_final_list = assess_several_orf_per_gene(final_list)
    return final_final_list


def edit_final_file(data, name1,name2,name3):
    f = open(name1, "w")
    nb_0 = 0
    nb_1 = 0
    nb_2 = 0
    for sublist in data:
        f.write(sublist[0])
        for elt in sublist[1:]:
            f.write(","+str(elt))
        f.write("\n")
    f.close()
    f = open(name2, "w")
    for sublist in data:
        f.write(">"+sublist[0]+"\n")
        f.write(str(sublist[9])+"\n")
        if sublist[8] == "0":
            nb_0+=1
        elif sublist[8] == "1":
            nb_1+=1
        elif sublist[8] == "2":
            nb_2+=1
    f.close()
    f = open(name3, "w")
    for sublist in data:
        f.write(">"+sublist[0]+"\n")
        f.write(str(sublist[10])+"\n")
    f.close()
    print (nb_0)
    print (nb_1)
    print (nb_2)
    print ("*******")
    
def parse_gtf(gtf):
    dico = {}
    for line in gtf:
        list_elts = line.split()
        if len(list_elts)>2 and list_elts[2] == "mRNA":
            data = list_elts[8]
            list_data = data.split(";")
            transcript = "anti-"+list_data[0].split("-")[1]
            gene = list_data[1].split("-")[2]
            #print (transcript)
            #print (gene)
            dico[transcript] = gene
    return dico
    

print ("AK")
gtf_AK = openFile("AK_realasORFs.gtf")
gtf_genome_AK = openFile("AK5deNovoGenomeAnnotation.gff")
dico_gtf_AK = parse_gtf(gtf_genome_AK)
list_gtf_AK = store_gtf(gtf_AK)
dico_AK = store_prots("AK_realasORF_prot.fa")
final_list_AK = get_all_info("AK_realasORFs.fa", list_gtf_AK,dico_AK,"AK",dico_gtf_AK)
edit_final_file(final_list_AK, "AK_datafile","AK_new_nuc.fa","AK_new_prot.fa")

print ("DK")
gtf_DK = openFile("DK_realasORFs.gtf")
gtf_genome_DK = openFile("DK5deNovoGenomeAnnotation.gff")
dico_gtf_DK = parse_gtf(gtf_genome_DK)
list_gtf_DK  = store_gtf(gtf_DK )
dico_DK = store_prots("DK_realasORF_prot.fa")
final_list_DK  = get_all_info("DK_realasORFs.fa", list_gtf_DK,dico_DK,"DK",dico_gtf_DK)
edit_final_file(final_list_DK, "DK_datafile","DK_new_nuc.fa","DK_new_prot.fa")

print ("GI")
gtf_GI = openFile("GI_realasORFs.gtf")
gtf_genome_GI = openFile("GI5deNovoGenomeAnnotation.gff")
dico_gtf_GI = parse_gtf(gtf_genome_GI)
list_gtf_GI = store_gtf(gtf_GI)
dico_GI = store_prots("GI_realasORF_prot.fa")
final_list_GI = get_all_info("GI_realasORFs.fa", list_gtf_GI,dico_GI,"GI",dico_gtf_GI)
edit_final_file(final_list_GI, "GI_datafile","GI_new_nuc.fa","GI_new_prot.fa")

print ("SW")
gtf_SW = openFile("SW_realasORFs.gtf")
gtf_genome_SW = openFile("SW5deNovoGenomeAnnotation.gff")
dico_gtf_SW = parse_gtf(gtf_genome_SW)
list_gtf_SW = store_gtf(gtf_SW)
dico_SW = store_prots("SW_realasORF_prot.fa")
final_list_SW = get_all_info("SW_realasORFs.fa", list_gtf_SW,dico_SW,"SW",dico_gtf_SW)
edit_final_file(final_list_SW, "SW_datafile","SW_new_nuc.fa","SW_new_prot.fa")

print ("UM")
gtf_UM = openFile("UM_realasORFs.gtf")
gtf_genome_UM = openFile("UMdeNovoGenomeAnnotation.gff")
dico_gtf_UM = parse_gtf(gtf_genome_UM)
list_gtf_UM = store_gtf(gtf_UM)
dico_UM = store_prots("UM_realasORF_prot.fa")
final_list_UM = get_all_info("UM_realasORFs.fa", list_gtf_UM,dico_UM,"UM",dico_gtf_UM)
edit_final_file(final_list_UM, "UM_datafile","UM_new_nuc.fa","UM_new_prot.fa")

print ("YE")
gtf_YE = openFile("YE_realasORFs.gtf")
gtf_genome_YE = openFile("YEdeNovoGenomeAnnotation.gff")
dico_gtf_YE = parse_gtf(gtf_genome_YE)
list_gtf_YE = store_gtf(gtf_YE)
dico_YE = store_prots("YE_realasORF_prot.fa")
final_list_YE = get_all_info("YE_realasORFs.fa", list_gtf_YE,dico_YE,"YE",dico_gtf_YE)
edit_final_file(final_list_YE, "YE_datafile","YE_new_nuc.fa","YE_new_prot.fa")

print ("ZB")
gtf_ZB = openFile("ZB_realasORFs.gtf")
gtf_genome_ZB = openFile("ZambdeNovoGenomeAnnotation.gff")
dico_gtf_ZB = parse_gtf(gtf_genome_ZB)
list_gtf_ZB = store_gtf(gtf_ZB)
dico_ZB = store_prots("ZB_realasORF_prot.fa")
final_list_ZB = get_all_info("ZB_realasORFs.fa", list_gtf_ZB,dico_ZB,"ZB",dico_gtf_ZB)
edit_final_file(final_list_ZB, "ZB_datafile","ZB_new_nuc.fa","ZB_new_prot.fa")


















