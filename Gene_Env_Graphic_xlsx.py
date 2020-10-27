#! /usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import re
import os
from Bio import Entrez
from Bio import SeqIO
import StringIO
import xlsxwriter

files = glob.glob("*.gbff")
print("%i files to index" % len(files))
RefSeq_bct = SeqIO.index_db("RefSeqGB.idx", files, "genbank")

workbook = xlsxwriter.Workbook('Gene_Env_Graphic.xlsx')
worksheet = workbook.add_worksheet()

wrap = workbook.add_format({'text_wrap': True})

Pfam_list = []
IDs_list = []
associated_genes_list = []

filein1 = open("Three_IDs_RedoGluta_to_charge_list.txt", 'r')    # list of IDs 
for line in filein1:
    ligne = line.strip().split()
    IDs_list.append(ligne)
filein1.close()

filein2 = open("Pfam_Seven_prots_list_BP1728_422.txt", 'r')    # list of Pfam predictions 
for line in filein2:
    ligne = line.strip().split("\t")
    Pfam_list.append(ligne)
filein2.close()

print len(IDs_list), IDs_list
print len(Pfam_list), Pfam_list

for i in range (len(IDs_list)):
    ID_prot = IDs_list[i][0]
    print "ID_prot : ", ID_prot
    ID_nuc = IDs_list[i][2]
    print "ID_nuc :", ID_nuc
    if ID_nuc in RefSeq_bct:
        record = RefSeq_bct[ID_nuc]
        organism = record.annotations["source"]
        print "organism :", organism
        Genbank_list = []
        for feature in record.features:
            if (feature.type == "CDS" and feature.qualifiers.get('pseudogene') == None) and (feature.type == "CDS" and feature.qualifiers.get('pseudo') == None):
                CDS_list =[]
                CDS_list.append(feature.strand)
                location = str(feature.location)
                regex3 = re.compile ("[0-9]+")
                locat = regex3.findall(location)
                CDS_list.append(int(locat[0])+1)
                CDS_list.append(int(locat[1]))
                protein_id = feature.qualifiers.get('protein_id')
                if protein_id is None:
                    CDS_list.append('')
                else:
                    protein_id_2 = protein_id[0]
                    CDS_list.append(protein_id_2)
                locus_tag = feature.qualifiers.get('locus_tag')
                if locus_tag is None:
                    CDS_list.append('')
                else:
                    locus_tag_2 = locus_tag[0]
                    CDS_list.append(locus_tag_2)
                gene = feature.qualifiers.get('gene')
                if gene is None:
                    CDS_list.append('')
                else:
                    gene_2 = gene[0]
                    CDS_list.append(gene_2)
                product = feature.qualifiers.get('product')
                if product is None:
                    CDS_list.append('')
                else:
                    product_2 = product[0]
                    CDS_list.append(product_2)
                translation = str(feature.qualifiers.get('translation'))
                regex2 = re.compile("[A-Z]+")
                PROT = regex2.search(translation)
                PROT2 = PROT.group(0)
                CDS_list.append(PROT2)
                Genbank_list.append(CDS_list)
            elif (feature.type == "gene" and feature.qualifiers.get('pseudogene') != None) or (feature.type == "gene" and feature.qualifiers.get('pseudo') != None):
                CDS_list =[]
                CDS_list.append(feature.strand)
                location = str(feature.location)
                regex3 = re.compile ("[0-9]+")
                locat = regex3.findall(location)
                CDS_list.append(int(locat[0])+1)
                CDS_list.append(int(locat[1]))
                CDS_list.append("")
                locus_tag = feature.qualifiers.get('locus_tag')
                locus_tag_2 = locus_tag[0]
                CDS_list.append(locus_tag_2)
                CDS_list.append("")
                CDS_list.append("pseudogene")
                CDS_list.append("")
                Genbank_list.append(CDS_list)
        print "Genbank_list : ", Genbank_list
        for index, cds in enumerate(Genbank_list):
            if ID_prot in cds[3]:
                positionID = index
                locus_list = []
                print "positionID :", positionID
                print "Genbank_list[positionID][3]", Genbank_list[positionID][3]
                for offset in range (7):
                    start = positionID - 3
                    positionR = start + offset
                    print "positionR :", positionR
                    if  positionR < 0:
                        print "forget"
                        continue
                    else:
                        try:
                            sens = int(Genbank_list[positionR][0])
                            try:
                                espacement = int(Genbank_list[positionR + 1][1]) - int(Genbank_list[positionR][2]) - 1
                            except:
                                espacement = "bord"
                            locus = Genbank_list[positionR][3]
                            print "locus :", locus
                            if sens == -1:
                                locus2 = "<<" + locus + "<<"
                                locus_list.append("<<" + locus + "<<")
                            else:
                                locus2 = ">>" + locus + ">>"
                                locus_list.append(">>" + locus + ">>")
                            gene = Genbank_list[positionR][5]
                            locus_list.append(gene)
                            product = Genbank_list[positionR][6]
                            locus_list.append(product)
                            space2 = str(espacement)
                            locus_list.append(espacement)
                            Pfam_Dom_list =[]
                            ID_Pfam = Genbank_list[positionR][4]
                            for j, item in enumerate(Pfam_list):
                                ID_Pfam_file = Pfam_list [j][0]
                                if ID_Pfam_file == ID_Pfam :
                                    dom_def = str(Pfam_list [j][1]) + "|" + str(Pfam_list [j][8])
                                    Pfam_Dom_list.append(dom_def)
                                else:
                                    continue
                            locus_list.append(Pfam_Dom_list)
                            worksheet.write(i, 0, organism , wrap)
                            worksheet.write(i, offset + 1, locus2 + "\n" + ID_Pfam + "\n" + gene + "\n" + product + "\n" + space2 + "\n" + str(Pfam_Dom_list) , wrap)
                        except:
                            continue
                    associated_genes_list.append(locus_list)
            else:
                continue
    else:
        print "This ID_nuc was not found :", ID_nuc
        continue
workbook.close()
