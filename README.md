# UCD_BHLS_Thesis
Repository for all R scripts and results involved in my final year research project for the Biomedical Health and Life Sciences degree in University College Dublin. 

Project title - "Bioinformatic Analysis of Copy Number Alteration Signature Mechanisms in Ovarian Carcinoma".  

There were four method components involved in this project:
1. Confirming copy number alteration (CNA) signature findings from a previously published high grade serous ovarian cancer (HGSOC) study 
   (Macintyre et al 2018 - "Copy number signatures and mutational processes in ovarian carcinoma").
2. Differential gene expression profiling for the CNA signatures identified in 1.
3. RNA signature identification using the results of 2.
4. Survival analysis of CNA and RNA signature classifications.

Scripts for part 1: "Signatures" folder

signature_extrapolation - for CNA signature extrapolation for TCGA HGSOC data, as performed by Macintyre et al 2018.
All_Sig_Functions - all functions required for signature extrapolation and formatting of CNA signature results for downstream analysis.

Scripts for parts 2, 3 and 4: "Post_signature_analysis" folder

DE_RNA_survival - differential gene expression profiling, RNA signature identification and survival analysis.
Post_sig_functions - all extra functions required for the DE_RNA_survival script.

