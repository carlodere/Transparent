##################### TRANSPARENT2.0 ###################

Authors: Carlo Derelitto(1) and Daniele Santoni(1*)

(1)(*) Istituto di Analisi dei Sistemi ed Informatica "Antonio Ruberti" - CNR, Via dei Taurini 19, 00185, Rome, Italy.

(*) Corresponding author Ph.: +39 06 49937119. Fax: +39 06 49937137.
email address: daniele.santoni@iasi.cnr.it

The toolset/pipeline consists of three python scripts (at the bottom of this file you find an example reporting how to run the tool on a test-case contained in the package in input directory):

##################### main.py #################

This script takes in input a list of ENTREZ ID genes and provides a list of TFs potentially regulating genes in the list 

USAGE:
python3 main.py -l input_file_list [-o output_dir] 
-l --list	input gene file
-d --directory	output directory ("output" is the default output directory).

The input gene file (a csv or a text file) must contain a list of ENTREZ gene IDs (one ID for each line).
Each ID gene is parsed, all IDs that are not recognized as valid IDs are reported in the log file.  

The script builds the following output files:
OUTPUT
- "basename".txt (where basename.csv is the input gene file list) providing the table with all 636 TFs and correspondent P-values and adjusted P-values (TFs are ranked in adjusted P-value ascending order).
- "basename".log reporting all performed steps and warnings or errors.


##################### significant.py #################

This script builds a report with all significant identified TFs (P-value or adjusted P-value smaller than a threshold).
For each TF a complete list with all genes showing in their promoter sequences at least one TFBS of the given TFs is also provided.

USAGE
python3 significant.py -f basename.txt -l input_file_list [-p padj|pvalue] [-t P-value significativity threshold] 
INPUT
-f --file	output file of main.py
-l --list	input gene file (as in main.py)
-p --parameter	select TFs based on P-value or adjusted Pvalue (default padj)
-t --threshold	P-value significativity threshold (default 0.05)
-d --directory	output directory ("output" is the default output directory)

OUTPUT
- "basename"_significant_TF_threshold_parametr.csv
- "basename"_significant_TF_threshold_parametr_genes.txt


##################### create_network.py #################

USAGE
python3 create_network.py -f basename.txt [-p padj|pvalue] [-t P-value significativity threshold] [-l input_file_list]

-f --file	output file of main.py
-p --parameter	select TFs based on P-value or adjusted Pvalue (default padj)
-t --threshold	P-value significativity threshold (default 0.05)
-d --directory	output directory ("output" is the default output directory)
-l --list	input gene file

OUTPUT
- "basename"_network_threshold.tsv
- "basename"_network_threshold_genes.tsv (this file is generated when the -l parameter is used)
 
##################### Sample Run #################

1) python3 main.py -l input/Schizophrenia_0.3.csv 
2) python3 significant.py -f output/Schizophrenia_0.3.txt -l input/Schizophrenia_0.3.csv -t 0.01
3) python3 create_network.py -f output/Schizophrenia_0.3.txt -t 0.01
