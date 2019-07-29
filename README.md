Python 2.7. 

1. mgfextract_notebook.ipynb
• It extracts specific subsets from an MGF file by matching the original MGF
file from against the ground truth (Excel file);
• Uses the method (write_title) from mgf_extract.py to create the required MGF
files.

2. mgf_extract.py
• It splits the original MGF file into multiple folders (one spectrum per folder)
and then chooses the required spectra (based on the title given to each folder
which contains information such as the protein id and the protein sequence),
deletes the unwanted folder and merges the rest of them into one folder.

3. data_pre_analysis.ipynb
• It was used for performing the analysis for the first part of the project.
• It uses mgf_parser.py for parsing the mgf file.
• It uses mgf_analysis.py for assigning and computing all the required values to
each spectrum.
• Creates the histogram of fragment differences found in both phosphorylated
and non-PTM spectra. Takes some time to compute (~7 min), due to v. large
number of differences.
• Calculates the percentage of spectra matching with fragments matching the
calculated values and the percentages of y and b ions for each matching
spectrum.
• Creates the graph for determining the mass tolerance threshold.

4. mgf_parser.py
• The main code for this was developed by PhD student Joe Wandy. I’ve only
adapted it for parsing a slightly modified format of mgf file (with a few extra
lines added to it for 1)protein id; 2)mh+; 3)protein sequence)

5. mgf_analysis.py
• Assigns and computes all the required values to each spectrum as shown in
Figure 1.
• It uses chem_utils.py and the ‘Document’ and ‘Peak’ objects from
ms_objects.py.

6. chem_utils.py
• It creates a dictionary with the masses of all amino acids and calculates the
masses of y and b ions when given a protein sequence.

7. ms_objects.py
• It is only used for creating useful objects to be used in other scripts.

8. new_lda_run_phos_eta_1.ipynb (all three: phos_eta_1, phos_eta_01, acetyl)
• It was used for performing the analysis for the second part of the project.
• It uses feature_extractor.py (feature_extractor_acetyl.py for acetylated
spectra) for preparing the dictionary (‘corpus’) which is going to be used for
the LDA run.
• It creates the PCA plot
• It splits the data to be used for the classifier and performs the evaluation for
each classifier.

9. feature_extractor.py
• It is used for grouping very similar differences (within a certain mass
tolerance) from the list of differences of all spectra. The intensity assigned to
each group within a spectrum is the average of all the intensities assigned to
the differences in the respective group (each difference is assigned the lower
intensity between the two fragments).
• A word is assigned to each group base on the mean value of the differences in
the group, i.e diff_mean-value

. Each MS1 peak is identified by its Peak ID (extracted from the ground truth), m/z, rt, intensity and charge and
MS2 values (extracted from the MGF file). For each MS1 spectra tandem MS values were matched with
theoretical values of fragments and the matching percentage was calculated
