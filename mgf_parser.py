import numpy as np
import pandas as pd
import re
import math
import itertools
import seaborn
import warnings
warnings.filterwarnings('ignore')

def parse_mgf(filename, debug=False):
    ms1_peakids = []
    ms1_peakdata = []
    ms2_peakids = []
    ms2_peakdata = []
    with open(filename, "r") as ins:

        pep_mass = None
        pep_rt = None
        pep_charge = np.nan
        fragments = []
        peak_id = 1
        for line in ins:

            line = line.strip()
            if not line:
                continue # skip empty line

            # split by ' ' or '='
            tokens = re.split(' |=', line)
            tok = tokens[0].upper()

            if tok == 'BEGIN':
                continue
            elif tok == 'TITLE':
                continue
            elif tok == 'RTINSECONDS':
                pep_rt = float(tokens[1])
                ms1_id = peak_id
                peak_id += 1
            elif tok == 'PEPMASS':
                pep_mass = float(tokens[1])
            elif tok == 'CHARGE':
                pep_charge = tokens[1]
            elif tok == 'END':

                if debug:
                    print ms1_id, pep_mass, pep_rt, pep_charge
                ms1_peakdata.append((ms1_id, np.nan, 1, pep_mass, pep_rt, 0, pep_charge))
                ms1_peakids.append(ms1_id)

                for ms2_id, ms2_mass, ms2_intensity in fragments:
                    if debug:
                        print '- %d %f %f' % (ms2_id, ms2_mass, ms2_intensity)
                    ms2_peakdata.append((ms2_id, ms1_id, 2, ms2_mass, 0, ms2_intensity, np.nan))
                    ms2_peakids.append(ms2_id)
                if debug:
                    print

                # reset for the next line
                pep_mass = None
                pep_rt = None
                pep_charge = np.nan
                fragments = []

            else: # read the fragments
                ms2_mass = float(tok)
                ms2_intensity = float(tokens[1])
                fragments.append((peak_id, ms2_mass, ms2_intensity))
                peak_id += 1

    ms1 = pd.DataFrame(ms1_peakdata, index=ms1_peakids,
                       columns=['peakID', 'MSnParentPeakID', 'msLevel', 'mz', 'rt', 'intensity', 'charge'])
    ms2 = pd.DataFrame(ms2_peakdata, index=ms2_peakids,
                       columns=['peakID', 'MSnParentPeakID', 'msLevel', 'mz', 'rt', 'intensity', 'charge'])

    return ms1, ms2

def parse_modified_mgf(filename, debug=False):
    ms1_peakids = []
    ms1_peakdata = []
    ms2_peakids = []
    ms2_peakdata = []
    with open(filename, "r") as ins:

        pep_mass = None
        pep_rt = None
        pep_charge = np.nan
        pep_id = np.nan
        pep_seq = np.nan
        pep_mh = None
        pep_ys = None
        fragments = []
        peak_id = 1

        for line in ins:

            line = line.strip()
            if not line:
                continue # skip empty line

            # split by ' ' or '='
            tokens = re.split(' |=', line)
            tok = tokens[0].upper()

            if tok == 'BEGIN':
                continue
            elif tok == 'TITLE':
                continue
            elif tok == 'RTINSECONDS':
                pep_rt = float(tokens[1])
                ms1_id = peak_id
                peak_id += 1
            elif tok == 'PEPMASS':
                pep_mass = float(tokens[1])
            elif tok == 'CHARGE':
                pep_charge = tokens[1]
            #original end of the .mgf file
            elif tok == 'END':
                continue
            #the extra lines added
            elif tok == 'ID':
                pep_id = tokens[1]
            elif tok == 'SEQ':
                pep_seq = tokens[1]
            elif tok == 'MH+':
                pep_mh = tokens[1]
            elif tok == 'YS':
                pep_ys = tokens[1]
            elif tok == 'FINISH':

                if debug:
                    print pep_mass, pep_rt, pep_charge, pep_id,pep_seq, pep_mh, pep_ys
                ms1_peakdata.append((pep_id, np.nan, 1, pep_mass, pep_rt, 0, pep_charge, pep_seq, pep_mh, pep_ys))
                ms1_peakids.append(pep_id)

                for ms2_id, ms2_mass, ms2_intensity in fragments:
                    if debug:
                        print '- %d %f %f' % (ms2_id, ms2_mass, ms2_intensity)
                    ms2_peakdata.append((ms2_id, pep_id, 2, ms2_mass, 0, ms2_intensity, np.nan))
                    ms2_peakids.append(ms2_id)
                if debug:
                    print

                # reset for the next line
                pep_mass = None
                pep_rt = None
                pep_charge = np.nan
                pep_id = np.nan
                pep_seq = np.nan
                pep_mh = None
                pep_ys = None
                fragments = []

            else: # read the fragments
                ms2_mass = float(tok)
                ms2_intensity = float(tokens[1])
                fragments.append((peak_id, ms2_mass, ms2_intensity))
                peak_id += 1

    ms1 = pd.DataFrame(ms1_peakdata, index=ms1_peakids,
                       columns=['peakID', 'MSnParentPeakID', 'msLevel', 'mz', 'rt', 'intensity', 'charge', 'seq', 'mh', 'ys'])
    ms2 = pd.DataFrame(ms2_peakdata, index=ms2_peakids,
                       columns=['peakID', 'MSnParentPeakID', 'msLevel', 'mz', 'rt', 'intensity', 'charge'])

    return ms1, ms2
