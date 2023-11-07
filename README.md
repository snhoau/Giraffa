# Giraffa
for lc-ms/ms identification (peptides and protein crosslink)

This script is only for ELinker

We have developed a dual-control crosslinking agent for efficient analysis of multi-species protein interactions. The entire analysis process is based on tandem mass spectrometry, and this program is used to analyze the mass spectrometry results and extract potential protein interaction information. For specific information, including the synthesis of the crosslinking agent, please refer to our published article.

Operating Environment
Giraffa dose not require installation. We provide the environment for reference. Since the program uses process locks, it is recommended to run it on Linux. If an error occurs during execution and prompts for library installation, please use pip to install the required libraries. The libraries that may need to be installed are as follows:

pickle, fcntl, numpy, collections, biopython, pyteomics, multiprocessing

The program includes:

1. A_RATvsSV_ELink_Mpos.py - Pre-computation of interactions, the generated files are used for subsequent programs.

2. A_RATvsSV_ELink_Mpos.py - Calculation of interaction information, outputting the interaction results.

supplymentary/shanghai_sv_mature_200_mix_fltd.fasta - Scorpion venom protein sequences (length <200)
supplymentary/uniprot-proteome-rat-UP000002494_fltd.fasta - Rat protein sequences

supplymentary/2-demo.mgf - Demo mgf file

# Usage

# 1. Pre-computation of interactions

The first step of pre-computation of interactions requires protein sequences that have been identified by conventional experimental methods as input, in fasta format, as well as an integrated mgf file (use a conventional raw to mgf conversion tool, but the mass spectrometry results of multiple runs need to be merged into one mgf file).

Possible modifications in script:
- line 782 - Merged mass spectrometry mgf file (merged)
- line 783 - Identified rat protein sequence fasta file
- line 784 - Identified scorpion venom peptide sequence fasta file
- line 802, 807 - Parameters for virtual digestion
- line 896 - Number of cores used for parallel processing

# 2. Interaction identification
Run after the first program has finished to perform interaction analysis.

Areas that need to be modified:
- line 825 - Identified rat protein sequence fasta file
- line 826 - Identified scorpion venom peptide sequence fasta file
- line 891 - Number of cores used for parallel processing
- line 913 - FDR threshold, or line 668 to select your preferred scoring function to adjust FDR

PS: The output files (PSM and Protein) are two result files. The Protein file contains the results of protein interactions, which can be opened in Excel and divided into columns. The content of each column can be referred to in the article's supplementary materials.
