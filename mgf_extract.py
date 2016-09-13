import numpy
import pandas as pd
import csv
import xlrd
import glob
import os
import re
import random
import fileinput


def write_title(df, ptm_type):

    title = "TITLE=sample=1 period=1 cycle="+df['cycle.Experiment'].str.split('.').str.get(0)+" experiment="+df['cycle.Experiment'].str.split('.').str.get(1)
    df['title'] = title
    title_name = "titles"+ ptm_type + ".txt"
    titles = open(title_name, "w")
        #write titles and sequences to file
    titles.write(df.to_csv(sep=':', index=False, header=False, columns=["title","cycle.Experiment", "mostCommonYSsequence", "Precursor_MH+_theoretical","numYSinRow"]))
    titles.close()

        #create the folder for individual mgf files

    mypath = "split_files" + ptm_type
    if not os.path.isdir(mypath):
        os.makedirs(mypath)

        #split the mgf file
    file_match = open("/Users/anamaria/Documents/git/iprg2012/mgf_files/iPRG2012match.mgf", "r")
    ofile = open(mypath + "/0.mgf", "w")

    for line in iter(file_match.readline, b''):

        mine = line.rstrip()

        if mine.startswith("BEGIN IONS"):
            ofile.close()
            title = file_match.readline()
            ofile = open(mypath+ "/"+ title +".mgf", "w")
            ofile.write("BEGIN IONS\n")
            ofile.write(title)

        else:
            ofile.write(line)
    file_match.close()

    #delete the files which don't have a title from titles.txt
    files = glob.glob(mypath +'/*.mgf')
    for f in files:
        if not isOk(f, title_name):
            os.remove(f)

    #merge remaining files in folder
    file_list = glob.glob(mypath + "/*.mgf")

    with open('/Users/anamaria/Documents/git/iprg2012/mgf_files/iPRG2012' + ptm_type +'.mgf', 'w') as file:
        input_lines = fileinput.input(file_list)
        file.writelines(input_lines)

    #delete the files which don't have a title from titles.txt
def isOk (f, titles_file):
    t=f.split("/")[1].split(".")[0].rstrip()
    title_file = open(titles_file, "r")
    for line in iter(title_file.readline, b''):
        titleseq = line.rstrip().split(":")
        #print titleseq
        if t==titleseq[0]:
            with open(f, "a") as bla:
                bla.write("ID="+titleseq[1]+"\n")
                bla.write("SEQ="+titleseq[2]+"\n")
                bla.write("MH+="+titleseq[3]+"\n")
                bla.write("YS="+titleseq[4]+"\nFINISH\n")
            bla.close()
            return True
