
"""
Viewing VCF with pandas and parsing and identifying SNPs using cyvcf2
Run the script as follows: python vcf_snp_finder.py PATH_TO_VCF_FILE

Author: Maryam Vazirabad
Date: 02/09/2022
"""
import sys
import io
import os
import cyvcf2
from cyvcf2 import VCF
import numpy as np
from pathlib import Path
import pandas as pd

def vcf_to_df(path): #option of converting vcf to pandas dataframe
     with open(path, "r") as filename:
          lines = [l for l in filename if not l.startswith('##')] #remove header except column names
          head = lines[0].strip().split("\t") #columns
          variant_data = [j.strip().split("\t") for j in lines[1:]] #all variant data
          return pd.DataFrame(variant_data, columns=head) #use info to create dataframe

def find_snps(vcf, quality_check = False): #searching for SNPs and printing them and their basic info
     snp_list = []
     count_snp = 0
     count_transversion = 0
     count_transition = 0
     for variant in vcf: #iterating through variants in file
          if variant.is_snp: #check if variant is SNP
               if quality_check is True: #if we care about quality of variant (default = False)
                    if variant.QUAL < 10: continue #if quality below threshold of 10, continue
               count_snp += 1
               snp_list.append(variant)
               if variant.is_transition: #check if SNP is transition
                    count_transition += 1
               else: #else SNP is transversion
                    count_transversion += 1
     print('\nSNPs found: ')
     print(snp_list)
     print('\nNumber of SNPs: {0}'.format(count_snp))
     print('Number of SNP Transitions: {0}'.format(count_transition))
     print('Number of SNP Transversions: {0}'.format(count_transversion))

if __name__ == "__main__":
     try:
          vcf_path = sys.argv[1] #path to vcf file (vcf or vcf.gz)
          assert vcf_path.endswith('.vcf') or vcf_path.endswith('vcf.gz'), 'Invalid file type'
     except IndexError:
          print("You did not specify a file")
          sys.exit(1) #abort due to error
     vcf = cyvcf2.VCF(Path(vcf_path)) #using cyvcf2 to iterate and query VCF file
     vcf_df = vcf_to_df(Path(vcf_path)) #function for converting vcf to pandas dataframe
     '''
     commenting out this block of code due to possible large input vcf file
     with pd.option_context('display.max_rows', None, 'display.max_columns',None):
          print(vcf_df)
     '''
     find_snps(vcf) #function for finding SNPs in VCF


