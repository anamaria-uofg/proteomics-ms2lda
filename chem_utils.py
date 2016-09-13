import sys
import numpy as np
import pandas as pd
import re
basedir = '/Users/anamaria/Documents/git/chem'
sys.path.append(basedir)
import formula as f

def get_aa_dictionary():
    mw = {}
    mw["A"] = f.Formula("C3H7NO2").compute_exact_mass()
    mw["R"] = f.Formula("C6H14N4O2").compute_exact_mass()
    mw["N"] = f.Formula("C4H8N2O3").compute_exact_mass()
    mw["D"] = f.Formula("C4H7NO4").compute_exact_mass()
    mw["C"] = f.Formula("C3H7NO2S").compute_exact_mass()
    mw["Q"] = f.Formula("C5H10N2O3").compute_exact_mass()
    mw["E"] = f.Formula("C5H9NO4").compute_exact_mass()
    mw["G"] = f.Formula("C2H5NO2").compute_exact_mass()
    mw["H"] = f.Formula("C6H9N3O2").compute_exact_mass()
    mw["I"] = f.Formula("C6H13NO2").compute_exact_mass()
    mw["L"] = f.Formula("C6H13NO2").compute_exact_mass()
    mw["K"] = f.Formula("C6H14N2O2").compute_exact_mass()
    mw["M"] = f.Formula("C5H11NO2S").compute_exact_mass()
    mw["F"] = f.Formula("C9H11NO2").compute_exact_mass()
    mw["P"] = f.Formula("C5H9NO2").compute_exact_mass()
    mw["S"] = f.Formula("C3H7NO3").compute_exact_mass()
    mw["T"] = f.Formula("C4H9NO3").compute_exact_mass()
    mw["W"] = f.Formula("C11H12N2O2").compute_exact_mass()
    mw["Y"] = f.Formula("C9H11NO3").compute_exact_mass()
    mw["V"] = f.Formula("C5H11NO2").compute_exact_mass()

    mw["HPO3"] = f.Formula("HPO3").compute_exact_mass()
    mw["H3PO4"] = f.Formula("H3PO4").compute_exact_mass()
    mw["C2H2O"] = f.Formula("C2H2O").compute_exact_mass()
    mw["C2H3O"] = f.Formula("C2H3O").compute_exact_mass()

    mw["s"] = mw.get("S") + mw.get("HPO3")
    mw["t"] = mw.get("T") + mw.get("HPO3")
    mw["y"] = mw.get("Y") + mw.get("HPO3")

    mw["carbon"] = f.Formula("C").compute_exact_mass()
    mw["H"] = f.Formula("H").compute_exact_mass()
    mw["O"] = f.Formula("O").compute_exact_mass()
    mw["OH"] = f.Formula("OH").compute_exact_mass()
    mw["H2O"] = f.Formula("H2O").compute_exact_mass()
    mw["proton"] = 1.007276466879

    return mw

#If given a specific peptide sequence, the methods below will calculate
#the masses of a, b and y ions. Generally the ions observed at
#the highest intensities are the b and y ions, which are formed following
#the cleavage of the peptide bond; a ions occur at lower intensities and
#they are formed from b ions following a CO loss.
#The masses were calculated using the following formulas:
#mass_yion = sum(aa residues) + proton + h + oh,
#mass_bion = sum(aa residues) + proton,
#mass_aion = mass_bion - CO,
#where sum(aa residues) = sum of each aa - H20 * no of aa.

def get_ymass(seq):
    mw = get_aa_dictionary()
    ymass = 0
    no_aa = len(seq)
    for l in range(0,no_aa):
        ymass += mw.get(seq[l:l+1])
    ymass += mw.get("OH") + mw.get("H") + mw.get("proton") - (no_aa * mw.get("H2O"))
    return ymass

def get_bmass(seq):
    mw = get_aa_dictionary()
    bmass = 0
    no_aa = len(seq)
    for l in range(0,no_aa):
        bmass += mw.get(seq[l:l+1])
    bmass += mw.get("proton") - (no_aa * mw.get("H2O"))
    return bmass

def get_amass(seq):
    mw = get_aa_dictionary()
    amass = get_bmass(seq) - mw.get("carbon") - mw.get("O")
    return amass

#Next, a method for calculating the mass of all the subfragments,
#i.e. a, b and y ions, when given a peptide sequence.
#The masses of the subfragments stemming from the non-modified peptide are also included.
#This method will return a dataframe containing the mass and the sequence of those fragments.

def get_submass(seq):
    mass = []
    upperseq = seq.upper()

    # obtain y ions masses and sequences
    for y in range(1,len(seq)):
        mass.append((seq[-y:],get_ymass(seq[-y:]), "y"))
        if seq[-y:].strip() != seq[-y:].strip().upper():
            mass.append((seq[-y:].strip().upper(), get_ymass(seq[-y:].strip().upper()), "y"))


    # obtain a and b ions masses and sequences
    for b in range(1,len(seq)+1):
        mass.append((seq[:b],get_bmass(seq[:b]), "b"))
        mass.append((seq[:b],get_amass(seq[:b]), "a"))
        if seq[:b] is not upperseq[:b]:
            mass.append((upperseq[:b],get_bmass(upperseq[:b]), "b"))
            mass.append((upperseq[:b],get_amass(upperseq[:b]), "a"))


    return mass
