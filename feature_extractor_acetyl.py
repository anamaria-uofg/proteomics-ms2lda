import sys
basedir = '/Users/anamaria/Documents/git/iprg2012/file_process'
sys.path.append(basedir)

from document import Document
from peak import Peak
from chem_utils import *
import pylab as plt
from IPython.display import display, HTML
import numpy as np
import pandas as pd
import seaborn as sns
import re
import itertools
import math
from operator import itemgetter



class FeatureExtractor (object):

    def __init__(self, ms1, ms2, start_filter_range,
                 end_filter_range, grouping_tol=20, search_tol = 20):

        self.ms1 = ms1
        self.ms2 = ms2
        self.start_filter_range = start_filter_range
        self.end_filter_range = end_filter_range

        self.grouping_tol = grouping_tol
        self.search_tol = search_tol

        self.peak_list = [] #vector containing peak objects (which contain ms2 peaks)
        self.document_list = [] #vector containing document objects
                                #(which contain ms1 peaks with their associated ms2 peaks)

        print "Starting extraction of peaks"

        for _,row in self.ms2.iterrows():

            peak_id= row['peakID']
            parent_id = row['MSnParentPeakID']
            mz = row['mz']
            rt = row['rt']
            intensity = row['intensity']
            charge = row['charge']

            p = Peak(peak_id, parent_id, mz, rt, intensity, charge)
            self.peak_list.append(p)

        print "Finished ms2 peaks extraction\n Starting ms1 peak extraction"

        for _,row in self.ms1.iterrows():

            peak_id= row['peakID']
            level = row['msLevel']
            mz = row['mz']
            rt = row['rt']
            intensity = row['intensity']
            charge = row['charge']

            ms2_spectra_peak = []
            #add the ms2 peaks to their respective spectra
            for peak in self.peak_list:
                if peak_id == peak.parent_id:
                    ms2_spectra_peak.append(peak)

            d = Document(peak_id, mz, rt, intensity, charge, ms2_spectra_peak)
            self.document_list.append(d)

        print "Adding words to the peak list:"
        #add another attribute to ms1 list containing the words for their respective losses
        self._add_words(self.document_list)

    #add the words to the diff_list object
    def _add_words(self, doc_list):
        words = []
        for doc in doc_list:
            #calculate the differences between ms2 mz
            mz_diff_list = self._create_mz_diff_list(doc)
            #annotate the groups of differences
            words = self._group_mz_diff_list(mz_diff_list)
            #check if there is any group within the given tolerance for
            #HPO3, H3PO4, ser, tyr, thr and annotate them accordingly
            #use this only for phos
            doc.words = self._discretize(words)
            #doc.words = words
            print "    Words added to spectra ", doc.peak_id

    #create the list of mz differences
    def _create_mz_diff_list(self, doc):
        #use itertools to create a list of combinations for the spectra attribute of doc
        #it's faster than double looping
        combo = itertools.combinations(doc.spectra, 2)
        mz_diff_list = []

        for c in combo:
            #calculate the mz difference
            mz_diff = abs(c[1].mz-c[0].mz)
            #determine the lower intensity between the two fragments
            intensity = min(c[1].intensity, c[0].intensity)


            if mz_diff > self.start_filter_range and mz_diff < self.end_filter_range: #because there are a lot of differences around 0
                mz_diff_list.append((mz_diff, intensity))
        #return a sorted list of mz_diff and their respective intensity value
        mz_diff_list_sorted = mz_diff_list.sort()
        return mz_diff_list

    #group the list according to a given grouping_tolerance value
    def _group_mz_diff_list(self, mz_diff_list): #access after _create_mz_diff_list
        #if the group's length is larger than 0
        if len(mz_diff_list) > 0:
            group = {}
            words = []
            group_number = 0

            current_group = []
            current_group.append(mz_diff_list[0])

            group[group_number] = current_group

            for index in range(1,len(mz_diff_list)-1):

                if self._belongs_to(mz_diff_list[index],group[group_number]):
                    current_group.append(mz_diff_list[index])
                    group[group_number] = current_group
                else:
                    group_number += 1
                    current_group = []
                    current_group.append(mz_diff_list[index])
                    group[group_number] = current_group

            words = self._generate_words(group)



        else:
            words = [('none',0.00001,0.000001)]

        return words



    def _generate_words(self, group):

        words=[]

        for key, value in group.items():
            mz_total = 0
            mz_count = 0
            intensity_list = []
            for mz, intensity in value:
                mz_total += mz
                mz_count += 1
                intensity_list.append(intensity)
            mz_mean = mz_total/mz_count
            word_mz_mean = "diff_" + "%.2f" % round(mz_mean, 2)
            final_intensity = max(intensity_list)

            words.append((word_mz_mean, final_intensity, mz_mean))

        return words
        if group is None:
            print "no group"


    def _belongs_to(self, mz_diff, group):

        tol_val = self.grouping_tol
        x=len(group)-1
        current_value = mz_diff[0]

        if x > 0:
            previous_value = group[x][0]
        else:
            previous_value = group[0][0]

        tolerance_max, _ = self._calculate_tolerance(previous_value, tol_val)
        if current_value <= tolerance_max:
            return True
        else:
            return False


    def _discretize(self, words):

        words_new = []
        tol = self.search_tol

        mw = get_aa_dictionary()

        h2o = mw.get("H2O")
        prot = mw.get("proton")
        oh = mw.get("OH")
        o = mw.get("O")
        h = mw.get("H")
        c = mw.get("carbon")

        h3po4 = mw.get("H3PO4")
        hpo3 = mw.get("HPO3")
        hpo3h2o = h2o+hpo3

        ser = mw.get("s")
        tyr = mw.get("y")
        thr = mw.get("t")

        tol_hpo3_max, tol_hpo3_min = self._calculate_tolerance(hpo3, tol)
        tol_h3po4_max, tol_h3po4_min = self._calculate_tolerance(h3po4, tol)
        tol_ser_max, tol_ser_min = self._calculate_tolerance(ser, tol)
        tol_tyr_max, tol_tyr_min = self._calculate_tolerance(tyr, tol)
        tol_thr_max, tol_thr_min = self._calculate_tolerance(thr, tol)

        for word, intensity, mean in words:
            if mean >= tol_hpo3_min and mean <=  tol_hpo3_max:
                word = "hpo3"
            elif mean >= tol_h3po4_min and mean <=  tol_h3po4_max:
                word = "h3po4"
            elif mean >= tol_ser_min and mean <=  tol_ser_max:
                word = "ser"
            elif mean >= tol_tyr_min and mean <=  tol_tyr_max:
                word = "tyr"
            elif mean >= tol_thr_min and mean <=  tol_thr_max:
                word = "thr"
            words_new.append((word, intensity, mean))


        return words_new
    def _discretize(self, words):

        words_new = []
        tol = self.search_tol

        mw = get_aa_dictionary()

        acetyl = mw.get("C2H2O")
       

        tol_acetyl_max, tol_acetyl_min = self._calculate_tolerance(acetyl, tol)
        
        for word, intensity, mean in words:
            if mean >= tol_acetyl_min and mean <=  tol_acetyl_max:
                word = "acetyl"
            words_new.append((word, intensity, mean))


        return words_new

    def _calculate_tolerance(self, number, tolerance):
        upper = number + number * tolerance * 1e-6
        lower = number - number * tolerance * 1e-6
        return (upper, lower)

    def create_corpus(self, intensity_threshold):

        print "Initialising dictionary creation:"
        self.corpus = {}

        for n in range(len(self.document_list)):

            doc = self.document_list[n]
            doc_id = doc.peak_id

            counts = {}

            for word, intensity, mean in doc.words:

                if word in counts:
                    counts[word] += intensity
                else:
                    counts[word] = intensity

            self.corpus[doc_id] = counts

        self.corpus_final = {}

        for doc_id in self.corpus:

            count = self.corpus[doc_id]
            count_final = {}

            for word in count:
                c = count[word]
                if c > intensity_threshold:
                    count_final[word] = c

            self.corpus_final[doc_id] = count_final

        print "Dictionary created."
        return self.corpus_final
        #return self.corpus
