import sys
basedir = '/Users/anamaria/Documents/git/iprg2012/file_process'
sys.path.append(basedir)
from chem_utils import *
from ms_objects import Document, Peak


import numpy as np
import pandas as pd
import re
import math
import itertools


class AnalysisMGF(object):

    def __init__(self, ms1, ms2):
        self.ms1 = ms1
        self.ms2 = ms2
        #self.seq_tol = seq_tol
        #self.diff_tol = diff_tol

        self.peak_list = [] #vector containing peak objects (which contain ms2 peaks)
        self.document_list = [] #vector containing document objects

        mw= get_aa_dictionary()

        self.h3po4 = mw.get("H3PO4")
        self.hpo3 = mw.get("HPO3")

        self.ser = mw.get("s")
        self.tyr = mw.get("y")
        self.thr = mw.get("t")

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

            seq = row['seq']
            #print seq
            d.seq = seq

            if self._check_match(seq):
                d.ms2_theoretical_values = get_submass(seq)
            else:
                d.ms2_theoretical_values = [(None, 0,None)]

            self.document_list.append(d)
        print "Finished ms1 peak extraction"

        #print "Adding theoretical sequences and ion types"
        #self.add_theoretical_seq_ion(self.seq_tol)
        #print "Adding diff list"
        #self.add_mz_diff_list(self.diff_tol)
        #print "Done"

    def add_theoretical_seq_ion(self, seq_tol, diff_tol):

        for ms1 in self.document_list:
            print "Evaluating spectra:", ms1.peak_id

            for ms2_spectra in ms1.spectra:

                experimental_mz = ms2_spectra.mz
                ms2_spectra.poss_seq = None
                ms2_spectra.poss_ion = None

                tol_max, tol_min = self._calculate_tolerance(experimental_mz, seq_tol)

                for seq, mz, ion in ms1.ms2_theoretical_values:

                    if mz > tol_max:
                        break
                    if mz >= tol_min and mz<= tol_max:
                        ms2_spectra.poss_seq = seq
                        ms2_spectra.poss_ion = ion

            self._add_match_list(ms1)

        for ms1 in self.document_list:
            self._calculate_mz_diff_list(ms1)
        self._add_diff_names(diff_tol)
        print "Done."

    def _add_match_list(self, ms1):
        match = 0
        total = 0
        y_ion = 0
        b_ion = 0
        a_ion = 0
        match_list = [0,0,0,0]
        for ms2 in ms1.spectra:
            total += 1

            if ms2.poss_seq is not None:
                match += 1
                if ms2.poss_ion == 'y':
                    y_ion += 1
                elif ms2.poss_ion == 'b':
                    b_ion += 1
                elif ms2.poss_ion == 'a':
                    a_ion += 1
                #print ms1_spectra.peak_id, mass2.mz, mass2.poss_seq, mass2.poss_ion
        match_list = ((self._percent(match, total), self._percent(y_ion, match),
                     self._percent(b_ion, match), self._percent(a_ion, match)))

        ms1.match_seq = match_list

    #create the list of mz differences
    #def _add_mz_diff_list(self, tolerance_val):


    def _calculate_mz_diff_list(self,doc):
        #use itertools to create a list of combinations for the spectra attribute of doc
        #it's faster than double looping
        combo = itertools.combinations(doc.spectra, 2)
        mz_diff_list = []

        for c in combo:
            #calculate the mz difference
            mz_diff = abs(c[1].mz-c[0].mz)
            #determine the lower intensity between the two fragments
            intensity = min(c[1].intensity, c[0].intensity)
            #because there are a lot of differences around 0
            mz_diff_list.append((mz_diff, intensity,
                                 c[1].poss_seq, c[1].poss_ion,
                                c[0].poss_seq, c[0].poss_ion))
            #return a sorted list of mz_diff and their respective intensity value
        mz_diff_list_sorted = mz_diff_list.sort()
        doc.mz_diff_list = mz_diff_list


    def _add_diff_names(self, tolerance_val):

        tol_hpo3_max, tol_hpo3_min = self._calculate_tolerance(self.hpo3, tolerance_val)
        tol_h3po4_max, tol_h3po4_min = self._calculate_tolerance(self.h3po4, tolerance_val)
        tol_ser_max, tol_ser_min = self._calculate_tolerance(self.ser, tolerance_val)
        tol_tyr_max, tol_tyr_min = self._calculate_tolerance(self.tyr, tolerance_val)
        tol_thr_max, tol_thr_min = self._calculate_tolerance(self.thr, tolerance_val)

        for ms1 in self.document_list:
            ms1.HPO3 = ((None, None, None))
            ms1.H3PO4 = ((None, None, None))
            ms1.ser = ((None, None, None))
            ms1.tyr = ((None, None, None))
            ms1.thr = ((None, None, None))
            for diff,intensity, seq1, ion1, seq0, ion0 in ms1.mz_diff_list:
                if diff >= tol_hpo3_min and diff <= tol_hpo3_max:
                    ms1.HPO3 = ((diff, seq1, seq0))
                elif diff >= tol_h3po4_min and diff <= tol_h3po4_max:
                    ms1.H3PO4 = ((diff, seq1, seq0))
                elif diff >= tol_ser_min and diff <= tol_ser_max:
                    ms1.ser = ((diff, seq1, seq0))
                elif diff >= tol_tyr_min and diff <= tol_tyr_max:
                    ms1.tyr = ((diff, seq1, seq0))
                elif diff >= tol_thr_min and diff <= tol_thr_max:
                    ms1.thr = ((diff, seq1, seq0))
                    
    def add_loss(self, tolerance_val):
        tol_hpo3_max, tol_hpo3_min = self._calculate_tolerance(self.hpo3, tolerance_val)
        tol_h3po4_max, tol_h3po4_min = self._calculate_tolerance(self.h3po4, tolerance_val)
        tol_ser_max, tol_ser_min = self._calculate_tolerance(self.ser, tolerance_val)
        tol_tyr_max, tol_tyr_min = self._calculate_tolerance(self.tyr, tolerance_val)
        tol_thr_max, tol_thr_min = self._calculate_tolerance(self.thr, tolerance_val)


        for ms1 in self.document_list:
            print "Evaluating spectra:", ms1.peak_id
            parent_mz = ms1.mz
            loss = []
            ms1.lossHPO3 = None
            ms1.lossH3PO4 = None
            ms1.lossser = None
            ms1.losstyr = None
            ms1.lossthr = None
            
            for ms2_spectra in ms1.spectra:

                experimental_mz = ms2_spectra.mz
                loss = parent_mz - experimental_mz

                if loss >= tol_hpo3_min and loss <= tol_hpo3_max:
                    ms1.lossHPO3 = loss
                elif loss >= tol_h3po4_min and loss <= tol_h3po4_max:
                    ms1.lossH3PO4 = loss
                elif loss >= tol_ser_min and loss <= tol_ser_max:
                    ms1.lossser = loss
                elif loss >= tol_tyr_min and loss <= tol_tyr_max:
                    ms1.losstyr = loss
                elif loss >= tol_thr_min and loss <= tol_thr_max:
                    ms1.lossthr = loss


    def _check_match(self,seq):
        match = re.search('[krmnqpdewc]', seq)
        #continue the analysis only if it doesn't contain other PTMs
        if match == None:
            return True
        else:
            return False

    def _get_tolerance(self,calculated, observed):
        tolerance = abs((calculated - observed)/calculated) * math.pow(10,6)
        return tolerance

    def _calculate_tolerance(self,number, tolerance):
        upper = number + number * tolerance * 1e-6
        lower = number - number * tolerance * 1e-6
        return (upper, lower)

    def _percent(self,a,b):
        if b == 0:
            return 0
        else:
            return float(a)*100/float(b)
