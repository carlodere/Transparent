import os
import argparse
import logging

import pandas as pd
import numpy as np
from scipy.stats import hypergeom

# PATH
path_TRANSPARENT=os.path.dirname(os.path.abspath(__file__))
path_InputDir=os.path.join(path_TRANSPARENT, 'input')
path_OutputDir=os.path.join(path_TRANSPARENT, 'output')
path_AvarageDir=os.path.join(path_TRANSPARENT, 'averages_TF')

# TERMINAL OPTIONS
def opzioni():
    parser = argparse.ArgumentParser(prog='python3 main.py',
                                 description="Script to find the TFBS of the different TFs in the promoters of the genes provided in input")
    parser.add_argument("-l", "--list", help="Path of the file containing the list of target genes", required=True, action="store")
    parser.add_argument("-d", "--outputdir", help="Output directory")
    args = parser.parse_args()
    return(args)

# INPUT FILE PATH AND NAME
def path_file_name_in(args):
    if os.path.isfile(args.list):
        path_input=args.list
        path, input_filname=os.path.split(args.list)
        input_filname=input_filname.rsplit('.', 1)[0]
    elif os.path.isfile(f'{path_InputDir}/{args.list}'):
        path_input=f'{path_InputDir}/{args.list}'
        input_filname=args.list.rsplit('.', 1)[0]
    else:
        print('File not found')  
    return(path_input, input_filname)

# OUTPUT FILE PATH [-o]
def path_file_name_out(args, input_filname):
    if args.outputdir == None:
        args.outputdir = 'empty'
    if os.path.isdir(args.outputdir): 
        path_output=f'{args.outputdir}/{input_filname}'
    elif args.outputdir=='empty': 
        path_output=f'{path_OutputDir}/{input_filname}'
    elif os.path.isdir(args.outputdir)!=True and os.path.isdir(args.outputdir)!='empty':
        path_output=f'{path_OutputDir}/{args.outputdir}'
    return(path_output)

# REMOVING UNKNOWN GENE IDs 
def delete_no_genes(target_genes_list):
    try:
        target_genes=[]
        tot_genes_df=pd.read_table(f'{path_TRANSPARENT}/tot_transcripts.txt',  header=None, sep=' ')
        total_genes=np.array(tot_genes_df[0])
        target_genes_to_test=open(target_genes_list).readlines()
        for x in target_genes_to_test:
            gene=x.strip()
            gene=int(gene)
            if gene not in total_genes:
                logging.warning(f'ID {gene}: is not a gene')
            else:
                logging.info(f'ID {gene}: analyzed')
                target_genes.append(str(gene))
        return(target_genes)
    except:
        print('Error: Check if the table has the header and delete it')

# IT EXTRACTS FILE LIST OF outputdir
def files_list(path):
    file_list=[]
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            if not file.startswith('.'):
                file_list.append(file)
    return(file_list)

# COMPUTING HYPERGEOMETRIC TEST AND PADJ
def hypergeometric_test(x, M, n, N):
    """
    The hypergeometric distribution models drawing objects from a bin.
    - M is total number of objects
    - n is total number of Type I objects. 
    - x (random variate) represents the number of Type I objects in N drawn without replacement from the total population

    - http://en.wikipedia.org/wiki/Hypergeometric_distribution
    - https://www.biostars.org/p/66729/
    - http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html
    - http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.hypergeometric.html
    - http://stackoverflow.com/questions/6594840/what-are-equivalents-to-rs-phyper-function-in-python
    """

    # COMPUTING THE AVERAGE NUMBER OF TFBSs
    n_TF=0
    for file in os.listdir(path_AvarageDir):
        if file.endswith('.txt'):
            n_TF+=1


    assert n <= M
    assert x <= n
    assert N <= M
    pv_gt = hypergeom.sf(x-1, M, n, N)  # 1-cdf sometimes more accurate
    padj = pv_gt*n_TF
    return pv_gt, padj



def main():
    args=opzioni()
    path_input, output_file_name=path_file_name_in(args)
    path_output = path_file_name_out(args, output_file_name)
    logging.basicConfig(filename=f'{path_output}.log',
                        level=logging.INFO,
                        filemode='w',
                        format='%(asctime)s:%(levelname)s: %(message)s')
    averages_TF_file_list=files_list(path_AvarageDir)
    logging.info(f'Analized the {output_file_name} genes list')
    target_genes= delete_no_genes(path_input)
    TF_name=pd.read_table(f'{path_TRANSPARENT}/human_JASPAR_PWM_list.txt', sep='\t', header=None, index_col=0)
    mat=[]
    for file_name in averages_TF_file_list:
        zero_target_genes=0
        non_zero_target_genes=0
        vett_total_genes=[]
        vett_target_genes=[]
        output_file_name=file_name.replace('PWMTF_', '').strip('_averages.txt')
        print(output_file_name)
        logging.info(f'Processing model: {output_file_name}')  
        check=False
        if len(TF_name.loc[output_file_name])>1:
            check=True
        tot_genes_df=pd.read_table(os.path.join(path_TRANSPARENT, f'averages_TF/{file_name}'), sep=',')
        vett_mean=np.array(tot_genes_df['AVERAGES'])
        logging.info('Count total genes')
        zero_total_genes=np.count_nonzero(vett_mean==0)
        non_zero_total_gene=np.count_nonzero(vett_mean>0)
        vett_total_genes.append(zero_total_genes)
        vett_total_genes.append(non_zero_total_gene)
        logging.info('Count target genes')
        for gene_id in target_genes:
            gene=int(gene_id.strip())
            row=tot_genes_df.loc[(tot_genes_df['GENE_ID']==gene)].to_numpy()
            mean=row[0][3]
            if mean == 0:
                zero_target_genes+=1
            elif mean > 0:
                non_zero_target_genes+=1
        vett_target_genes.append(zero_target_genes)
        vett_target_genes.append(non_zero_target_genes)
        total_genes=non_zero_total_gene+zero_total_genes
        total_genes_target=non_zero_target_genes+zero_target_genes
        logging.info(f'Calculate hypergeometric_test')
        p_value, padj=hypergeometric_test(non_zero_target_genes, total_genes, non_zero_total_gene, total_genes_target)
        if check:
            for TF in list(TF_name.loc[output_file_name][1]):
                vett=[output_file_name, TF, zero_total_genes, non_zero_total_gene, zero_target_genes, non_zero_target_genes, p_value, padj]
                mat.append(vett)
        else:
            vett=[output_file_name, TF_name.loc[output_file_name][1], zero_total_genes, non_zero_total_gene, zero_target_genes, non_zero_target_genes, p_value, padj]
            mat.append(vett)
    df_pro=pd.DataFrame(mat, columns=['Model', 'Name', 'N_total_genes=0', 'N_total_genes>0', 'N_target_genes=0', 'N_target_genes>0', 'p-value', 'padj'])
    df_pro.sort_values(by=['p-value'], inplace=True)
    df_pro.to_csv(f'{path_output}.txt', index=False, sep='\t')
    logging.info('Done.')
    print('DONE')
    
    


if __name__ == "__main__":
    main()
