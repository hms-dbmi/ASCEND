"""
calc_diag_probs.py

Given:

1. The list of target variants (genome-wide)
2. The list of variants observed in a reference disease cohort
3. Some specific parameters for the run

Calculate the diagnostic probabilities for every variant in list (1)
"""

import numpy as np
import pandas as pd
import gzip
from scipy.stats import invgamma, beta
import time

import argparse

# Likelihood function of the point estimate of the pathogenic fraction per genome
raw_beta_path_frac_params = {
    "AlphaMissense_MT": (571.1010677043221, 153.5712226370307),
    "REVEL_MT": (406.165389810628, 108.31270140607201),
    "PrimateAI-3D_MT": (495.7292995271211, 143.43542373630802)
}

# Beta MLE for between-gene pathogenic fraction distribution
beta_path_frac_params_MLE = {
    'AlphaMissense_MT': (0.7452350231177307, 0.31256095811679524), 
    'PrimateAI-3D_MT': (0.554610415637113, 0.28278307632600785), 
    'REVEL_MT': (0.5423312433175935, 0.26324777945772304)}

# Parameters for per-score ClinVar pathogenicity model fits
# Format: (a_path, b_path, a_benign, b_benign, p_path)
ClinVar_path_params = {'AlphaMissense_MT': (1.771563986206257, 0.97911395626164, 
                                            7.930376747003526, 1.150831170930573, 
                                            0.7461230515823135), 
                       'REVEL_MT' :        (1.179611639965687, 0.80913006574919, 
                                            4.499996982418247, 0.8140249626667113, 
                                            0.753397159038067),
                       'PrimateAI-3D_MT' : (1.394137770284522, 0.87079995439776, 
                                            4.088291033885278, 0.8631403140888494, 
                                            0.846164729591571)}


##### 1. READ IN ANNOTATED VARIANTS FROM THE REFERENCE COHORT #####
# Takes in the output file from ramediesDN2_preprocess2: print_var_info file and summary stats file

def load_reference_variants(var_info_file, diagnostic_genes, verbose=False):
    """
    Load the annotated variants from the reference cohort
    (ramediesDN2_preprocess2/print_var_info output file)

    Args:
        var_info_file (str): Path to the variant information file
        diagnostic_genes (set): Set of diagnostic gene IDs
        verbose (bool): Whether to print verbose output
    """

    var_info_df = pd.read_csv(var_info_file, sep='\t')
    diag_var_info_df = var_info_df[var_info_df['ensembl_gene_id'].isin(diagnostic_genes)].reset_index(drop=True)
    if verbose:
        print(f"Loaded {diag_var_info_df.shape[0]} variants in diagnostic genes from {var_info_file}")
    return diag_var_info_df


def load_reference_sumstats(sum_stats_file, verbose=False):
    """
    Load per-gene summary statistics for the reference cohort
    """
    df = pd.read_csv(sum_stats_file, sep='\t')
    if verbose:
        print(f"Loaded sum stats for {df.shape[0]} genes from {sum_stats_file}")
    return df


##### 2. STATISTICS BLOCK #####

### ENRICHMENT-BASED POSTERIOR CALCULATIONS ###

# Beta mixture PDF
# Used for the construction of per-score pathogenicity P
#   distributions from ClinVar data
# As a vector of arguments, accepts (x, a1, b1, a2, b2, pi) (from ClinVar_path_params)
def beta_mixture_PDF(x, a1, b1, a2, b2, pi):
    return (1 - pi) * beta.pdf(x, a1, b1) + pi * beta.pdf(x, a2, b2)

# Expectation of T_1: the raw inverse-gamma parameter corresponding to 
#   the unscaled fraction of benign cases
def E_T1(mu, n, P_D, p_null = 0.0):
     # mu : mutation rate (# expected)
    # n : observed count
    # P_D : Disease prevalence
    # p_null : probability of misannotation (default 0.0)

    # Mean of the untruncated inverse-gamma distribution
    if n == 0:
        return np.nan
    mean_raw = mu / n
    # Effect of truncation on the mean
    lower_T1_val = P_D / (1 - p_null * (1 - P_D))
    mean_fac = (invgamma.cdf(1, n, scale=mu) - invgamma.cdf(lower_T1_val, n, scale=mu)) / \
               (invgamma.cdf(1, n + 1, scale=mu) - invgamma.cdf(lower_T1_val, n + 1, scale=mu))
    # Mean of the truncated inverse-gamma distribution
    T1 = mean_raw * mean_fac

    # Enforce bounds on T1
    # Case 1: if T1 < P_D, set T1 = P_D (all cases are pathogenic)
    if T1 < P_D:
        T1 = P_D
    # Case 2: if T1 > 1.0, set T1 = 1.0 (all cases are benign)
    elif T1 > 1.0:
        T1 = 1.0

    return mean_raw * mean_fac

# Expectation of T: the scaled fraction of pathogenic cases
def E_T(mu, n, P_D, p_null = 0.0):
    # mu : mutation rate (# expected)
    # n : observed count
    # P_D : Disease prevalence

    # Calculate E[T1]
    E_T1_val = E_T1(mu, n, P_D, p_null)
    # Affine transformation to get E[T]
    return (1 - E_T1_val) / (1 - P_D)


# MAIN ENRICHMENT POSTERIOR CALCULATION FUNCTION
def Enrich_PPost_calc(varN, lambda_param, P_D, p_null):
    # Calculate posterior mean from just the enrichment data
    # Used as the LoF posterior throughout and as Ppost1 for missenses
    # I use the mean estimate from the truncated Inverse Gamma approximation

    # varN : observed count for the variant
    # lambda_param : expected count for the variant
    # P_D : Disease prevalence
    # p_null : probability of misannotation

    PPost = E_T(lambda_param, varN, P_D, p_null)

    return PPost


def PPost_T(mu, n, P_D, T):
    """
    Full posterior PDF calculation for T given mu, n, P_D
    """

    T1 = 1 - T * (1 - P_D)
    Ppost_fac = (1 - P_D) / (invgamma.cdf(1, n + 1, scale=mu) - invgamma.cdf(P_D, n + 1, scale=mu))
    Ppost = Ppost_fac * invgamma.pdf(T1, n + 1, scale=mu)

    return Ppost


##### 3. MAIN FUNCTIONS TO CALCULATE DIAGNOSTIC PROBABILITIES #####
def load_genes(genes_file, verbose=False):
    """
    Load the list of whitelisted genes from a file
    Args:
        genes_file (str): Path to the file with whitelisted genes (one gene (ENS ID) per line)
    Returns:
        set: dedupped list of whitelisted gene IDs
    """
    genes_whitelist = {}
    with open(genes_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            gene_id = line.strip()
            if gene_id:
                genes_whitelist[gene_id] = True
    if verbose:
        print(f"Loaded {len(genes_whitelist)} gene IDs from {genes_file}")

    return genes_whitelist


def precompute_PPost1_Lof(diag_sumstats_df, P_D, pnull, verbose=False, add_one_to_varN=True):
    """
    Precompute LoF posterior probabilities for all genes
    Args:
        diag_sumstats_df (DataFrame): Per-gene summary statistics DataFrame for the reference cohort
        P_D (float): Disease prevalence
        pnull (float): Probability of misannotation
    Returns:
        dict: Dictionary mapping gene IDs to LoF posterior probabilities
    """
    if verbose:
        print("PPost1: Precomputing LoF posterior probabilities for all genes...")

    gene_sumstats_df_genes = diag_sumstats_df['ENS_ID'].unique().tolist()
    PPost1_LoF_dict = {}

    syn_N = diag_sumstats_df['syn_varN'].sum()
    for gene_id in gene_sumstats_df_genes:
        gene_sumstats_df_local = diag_sumstats_df[diag_sumstats_df['ENS_ID'] == gene_id].reset_index(drop=True)
        varN = gene_sumstats_df_local['Lof_varN'].values[0]
        if varN is None or pd.isna(varN):
            varN = 0

        if add_one_to_varN:
            varN = varN + 1 # !!!
        lambda_param = float(gene_sumstats_df_local['Lof_lambda'].values[0]) * syn_N

        PPost1_LoF_dict[gene_id] = Enrich_PPost_calc(varN, lambda_param, P_D, pnull)
        if varN == 0:
            PPost1_LoF_dict[gene_id] = np.nan

        #print(f"Gene {gene_id}: LoF varN = {varN}, lambda = {lambda_param}, PPost1_LoF = {PPost1_LoF_dict[gene_id]}")

    return PPost1_LoF_dict


score2lambda_colname_dict = {
    'AlphaMissense_MT': 'AM_lambda',
    'REVEL_MT': 'REVEL_lambda',
    'PrimateAI-3D_MT': 'PAI_lambda'
}

def precompute_PPost1_missense(diag_sumstats_df, diag_var_info_df, score_type, P_D, verbose=False, add_one_to_varN=True):
    """
    Precompute missense posterior probabilities for all genes and a given score type
    Args:
        diag_sumstats_df (DataFrame): Per-gene summary statistics DataFrame for the reference cohort
        score_type (str): Type of missense score ('AlphaMissense_MT', 'REVEL_MT', 'PrimateAI-3D_MT')
        P_D (float): Disease prevalence
    Returns:
        dict: Dictionary mapping gene IDs to missense posterior probabilities
    """
    if verbose:
        print(f"PPost1: Precomputing {score_type} missense posterior probabilities for all genes...")

    gene_sumstats_df_genes = diag_sumstats_df['ENS_ID'].unique().tolist()
    PPost1_missense_dict = {}

    syn_N = diag_sumstats_df['syn_varN'].sum()

    lambda_column = score2lambda_colname_dict[score_type]
    for gene_id in gene_sumstats_df_genes:
        gene_sumstats_df_local = diag_sumstats_df[diag_sumstats_df['ENS_ID'] == gene_id].reset_index(drop=True)
        lambda_param = float(gene_sumstats_df_local[lambda_column].values[0]) * syn_N

        if pd.isna(lambda_param) or lambda_param == 0.0:
            PPost1_missense_dict[gene_id] = np.nan
            continue

        var_info_df_local = diag_var_info_df[diag_var_info_df['ensembl_gene_id'] == gene_id].reset_index(drop=True)
        var_info_df_local = var_info_df_local[var_info_df_local[score_type].notna()].reset_index(drop=True)
        if add_one_to_varN:
            varN = np.sum(~pd.isna(var_info_df_local[score_type]) * var_info_df_local["SAMPLE_RECURRENCE"]) + 1 #!!!
        else:
            varN = np.sum(~pd.isna(var_info_df_local[score_type]) * var_info_df_local["SAMPLE_RECURRENCE"])

        PPost1_missense_dict[gene_id] = Enrich_PPost_calc(varN, lambda_param, P_D, p_null=0.0)
        if varN == 0:
            PPost1_missense_dict[gene_id] = np.nan

    return PPost1_missense_dict


def precompute_PPost2_missense(verbose=False):
    """
    Precompute score-based missense prior probabilities for all genes with fixed pathogenic fractions
    Args:
        verbose (bool): Whether to print verbose output
    Returns:
        dict: Dictionary mapping score types to missense posterior probabilities
    """
    
    PPost2_missense_dict = {}
    for score_type in ['AlphaMissense_MT', 'REVEL_MT', 'PrimateAI-3D_MT']:
        if verbose:
            print(f"PPost2: Precomputing {score_type} missense posterior probabilities for all genes...")
        
        beta_frac_a, beta_frac_b = raw_beta_path_frac_params[score_type]
        MLE_path_frac = beta_frac_a / (beta_frac_a + beta_frac_b)
        PPost2_missense_dict[score_type] = MLE_path_frac

    return PPost2_missense_dict


def precompute_PPost3_missense(diag_var_info_df, diag_sumstats_df, score_type, P_D, verbose=False):
    """
    Precompute score-based missense posterior probabilities for all genes with per-gene pathogenic fractions
    Args:
        diag_var_info_df (DataFrame): Variant information DataFrame for the reference cohort
        diag_sumstats_df (DataFrame): Per-gene summary statistics DataFrame for the reference cohort
        score_type (str): Type of missense score ('AlphaMissense_MT', 'REVEL_MT', 'PrimateAI-3D_MT')
        verbose (bool): Whether to print verbose output
    Returns:
        tuple: Two dictionaries mapping gene IDs to missense posterior probabilities (PPost3 and PPost4)
    """

    
    if verbose:
        print(f"PPost3&4: Precomputing {score_type} missense posterior probabilities for all genes...")

    gene_sumstats_df_genes = diag_sumstats_df['ENS_ID'].unique().tolist()
    PPost3_dict = {}
    PPost4_dict = {}
    T_range = np.linspace(0.001, 0.999, 1000)

    syn_N = diag_sumstats_df['syn_varN'].sum()
    lambda_column = score2lambda_colname_dict[score_type]

    a1, b1, a2, b2, pi = ClinVar_path_params[score_type]
    i = 0
    for gene_id in gene_sumstats_df_genes:
        i += 1
        if verbose and (i % 1000 == 0):
            print(f"  Processed {i} / {len(gene_sumstats_df_genes)} genes...")

        gene_sumstats_df_local = diag_sumstats_df[diag_sumstats_df['ENS_ID'] == gene_id].reset_index(drop=True)
        lambda_param = float(gene_sumstats_df_local[lambda_column].values[0]) * syn_N

        var_info_df_local = diag_var_info_df[diag_var_info_df['ensembl_gene_id'] == gene_id].reset_index(drop=True)
        var_info_df_local = var_info_df_local[var_info_df_local[score_type].notna()].reset_index(drop=True)

        if pd.isna(lambda_param) or lambda_param == 0.0:
            PPost3_dict[gene_id] = np.nan
            PPost4_dict[gene_id] = np.nan
            continue

        score_vec = var_info_df_local[score_type].values
        varN = np.sum(~pd.isna(var_info_df_local[score_type]) * var_info_df_local["SAMPLE_RECURRENCE"]) + 1
        
        # Score-based per-gene prior
        a_beta, b_beta = beta_path_frac_params_MLE[score_type]
        score_prior_PDF = beta.pdf(T_range, a_beta, b_beta)

        # Enrichment-based likelihood
        enrich_prior_PDF = PPost_T(lambda_param, varN, P_D, T_range)

        T_LL_arr = np.zeros(len(T_range))
        for score in score_vec:
            # Score-based likelihood
            f_P = beta_mixture_PDF(score, a1, b1, a2, b2, pi)
            f_B = 1.0

            T_LL = np.log(T_range * f_P + (1 - T_range) * f_B)
            T_LL_arr += T_LL
        
        # Combine the three components to get posterior PDF
        T_LL_arr = T_LL_arr - np.max(T_LL_arr)  # for numerical stability
        T_LogPost3 = T_LL_arr + np.log(score_prior_PDF + 1e-300)
        T_LogPost4 = T_LL_arr + np.log(enrich_prior_PDF + 1e-300)

        T_Post3_unnorm = np.exp(T_LogPost3 - np.max(T_LogPost3))
        T_Post4_unnorm = np.exp(T_LogPost4 - np.max(T_LogPost4))

        E_T3 = np.sum(T_range * T_Post3_unnorm) / np.sum(T_Post3_unnorm)
        E_T4 = np.sum(T_range * T_Post4_unnorm) / np.sum(T_Post4_unnorm)

        PPost3_dict[gene_id] = E_T3
        PPost4_dict[gene_id] = E_T4

        #print(score_vec)
        #print(f"Gene {gene_id}: {score_type} varN = {varN}, lambda = {lambda_param}, PPost3 = {PPost3_dict[gene_id]}, PPost4 = {PPost4_dict[gene_id]}")

    return PPost3_dict, PPost4_dict


def main_func(gene_list_file, # file containing the whitelist of genes
              ref_var_info_file, # reference cohort variant info file
              ref_sum_stats_file, # reference cohort per-gene summary stats file
              processed_vars_file, # File with variants to be processed (gzipped)
              output_file,        # Path to write the output instead of printing
              colname, # Prefix of the columns to be added
              pnull=0.0, # Assumed probability of misannotation
              verbose=True,
              P_D = 0.01,
              round_until=4,
              add_one_to_varN=True): 
    
    start_time = time.time()
    
    # Load whitelisted ENS IDs
    gene_whitelist_dict =  load_genes(gene_list_file, verbose=verbose)

    # Load reference cohort information
    diag_var_info_df = load_reference_variants(ref_var_info_file, [i for i in gene_whitelist_dict.keys()], verbose=verbose)
    diag_sumstats_df = load_reference_sumstats(ref_sum_stats_file, verbose=verbose)

    if verbose:
        t = time.time()
        print(f"Loaded reference data in {t - start_time:.2f} seconds.")

    # PRECOMPUTE GENOME_WIDE E_Ts here
    # Precompute genome-wide values for the LoF and missense variants
    # PPost1 -- enrichment-based posterior probabilities
    PPost1_LoF_dict = precompute_PPost1_Lof(diag_sumstats_df, P_D, pnull, verbose=verbose, add_one_to_varN=add_one_to_varN)

    PPost1_AlphaMissense_dict = precompute_PPost1_missense(diag_sumstats_df, diag_var_info_df, 'AlphaMissense_MT', P_D, verbose=verbose, add_one_to_varN=add_one_to_varN)
    PPost1_REVEL_dict = precompute_PPost1_missense(diag_sumstats_df, diag_var_info_df, 'REVEL_MT', P_D, verbose=verbose, add_one_to_varN=add_one_to_varN)
    PPost1_PrimateAI_dict = precompute_PPost1_missense(diag_sumstats_df, diag_var_info_df, 'PrimateAI-3D_MT', P_D, verbose=verbose, add_one_to_varN=add_one_to_varN)

    # PPost2 -- score-based posterior probabilities with fixed pathogenic fractions
    PPost2_priors = precompute_PPost2_missense(verbose)

    # PPost3 & PPost4 -- Precompute E_Ts (score-based) for every gene for every missense score
    PPost3_AlphaMissense_dict, PPost4_AlphaMissense_dict = precompute_PPost3_missense(diag_var_info_df, diag_sumstats_df, 'AlphaMissense_MT', P_D, verbose)
    PPost3_REVEL_dict, PPost4_REVEL_dict = precompute_PPost3_missense(diag_var_info_df, diag_sumstats_df, 'REVEL_MT', P_D, verbose)
    PPost3_PrimateAI_dict, PPost4_PrimateAI_dict = precompute_PPost3_missense(diag_var_info_df, diag_sumstats_df, 'PrimateAI-3D_MT', P_D, verbose)

    preprocess_time = time.time()
    if verbose:
        print(f"Preprocessing completed in {preprocess_time - start_time:.2f} seconds.")
        print("Starting main variant processing...")

    # Added columns
    add_header = [f"{colname}_PP_AlphaMissense_MT_1", 
                  f"{colname}_PP_AlphaMissense_MT_2", 
                  f"{colname}_PP_AlphaMissense_MT_3",
                  f"{colname}_PP_AlphaMissense_MT_4",
                  f"{colname}_PP_REVEL_MT_1", 
                  f"{colname}_PP_REVEL_MT_2", 
                  f"{colname}_PP_REVEL_MT_3",
                  f"{colname}_PP_REVEL_MT_4",
                  f"{colname}_PP_PrimateAI-3D_MT_1", 
                  f"{colname}_PP_PrimateAI-3D_MT_2", 
                  f"{colname}_PP_PrimateAI-3D_MT_3",
                  f"{colname}_PP_PrimateAI-3D_MT_4"]
    
    if verbose:
        print(f"Reading {processed_vars_file} and writing results to {output_file}")

    header_dict = None
    if processed_vars_file.endswith('.gz'):
        open_func, mode_in, mode_out = gzip.open, 'rt', 'wt'
    else:
        open_func, mode_in, mode_out = open, 'r', 'w'

    with open_func(processed_vars_file, mode_in) as infile, open_func(output_file, mode_out) as outfile:
        i = 0
        for raw_line in infile:
            i += 1
            if verbose and (i % 1000000 == 0):
                t = time.time()
                print(f"  Processed {i} input variants in {t - preprocess_time:.2f} seconds..., total time: {t - start_time:.2f} seconds")

            raw_line = raw_line.rstrip('\n')
            # header line(s) may start with '#'
            if raw_line.startswith('#'):
                cols = raw_line.lstrip('#').split('\t')
                header_dict = {cols[i]: i for i in range(len(cols))}
                outfile.write(raw_line + '\t' + '\t'.join(add_header) + '\n')
                continue

            parts = raw_line.split('\t')
            # If we haven't seen a header yet, assume the first non-# line is the header
            if header_dict is None:
                cols = parts
                header_dict = {cols[i]: i for i in range(len(cols))}
                outfile.write('\t'.join(parts + add_header) + '\n')
                continue

            syn_flag = parts[header_dict['Synonymous-flag']] == '1'
            if syn_flag:
                PP_vals = ['NA'] * len(add_header)
                outfile.write('\t'.join(parts + [str(x) for x in PP_vals]) + '\n')
                continue

            gene_id = parts[header_dict['ensembl_gene_id']]

            LoF_flag = parts[header_dict['LOFTEE-flag']] == '1'
            # Case of LoF variant: directly compute LoF posterior probability
            if LoF_flag:
                p_enrich = PPost1_LoF_dict.get(gene_id, np.nan)
                pp_LoF_list = [p_enrich] * len(add_header)  # same for all missense scores
                pp_LoF_list = [round(x, round_until) if not pd.isna(x) else 'NA' for x in pp_LoF_list]
                outfile.write('\t'.join(parts + [str(x) for x in pp_LoF_list]) + '\n')
                continue

            missense_flag = parts[header_dict['AlphaMissense_MT']] != "NA" or \
                            parts[header_dict['REVEL_MT']] != "NA" or \
                            parts[header_dict['PrimateAI-3D_MT']] != "NA"

            if missense_flag:
                # AlphaMissense_MT
                score_val = parts[header_dict['AlphaMissense_MT']]
                f_P_val = beta_mixture_PDF(float(score_val), *ClinVar_path_params['AlphaMissense_MT']) if score_val != 'NA' else np.nan
                f_B_val = 1.0

                PP1_AM = PPost1_AlphaMissense_dict.get(gene_id, np.nan)
                #print(f"Gene {gene_id}, AlphaMissense_MT score {score_val}, f_P {f_P_val}, PP1_AM {PP1_AM}")
                p2_AM = PPost2_priors['AlphaMissense_MT']
                PP2_AM = (p2_AM * f_P_val) / (p2_AM * f_P_val + (1 - p2_AM) * f_B_val) if not pd.isna(f_P_val) else np.nan
                p3_AM = PPost3_AlphaMissense_dict.get(gene_id, np.nan)
                PP3_AM = (p3_AM * f_P_val) / (p3_AM * f_P_val + (1 - p3_AM) * f_B_val) if not pd.isna(f_P_val) else np.nan
                p4_AM = PPost4_AlphaMissense_dict.get(gene_id, np.nan)
                PP4_AM = (p4_AM * f_P_val) / (p4_AM * f_P_val + (1 - p4_AM) * f_B_val) if not pd.isna(f_P_val) else np.nan

                # REVEL_MT
                score_val = parts[header_dict['REVEL_MT']]
                f_P_val = beta_mixture_PDF(float(score_val), *ClinVar_path_params['REVEL_MT']) if score_val != 'NA' else np.nan
                f_B_val = 1.0

                PP1_REVEL = PPost1_REVEL_dict.get(gene_id, np.nan)
                #print(f"Gene {gene_id}, REVEL_MT score {score_val}, f_P {f_P_val}, PP1_REVEL {PP1_REVEL}")
                p2_REVEL = PPost2_priors['REVEL_MT']
                PP2_REVEL = (p2_REVEL * f_P_val) / (p2_REVEL * f_P_val + (1 - p2_REVEL) * f_B_val) if not pd.isna(f_P_val) else np.nan
                p3_REVEL = PPost3_REVEL_dict.get(gene_id, np.nan)
                PP3_REVEL = (p3_REVEL * f_P_val) / (p3_REVEL * f_P_val + (1 - p3_REVEL) * f_B_val) if not pd.isna(f_P_val) else np.nan
                p4_REVEL = PPost4_REVEL_dict.get(gene_id, np.nan)
                PP4_REVEL = (p4_REVEL * f_P_val) / (p4_REVEL * f_P_val + (1 - p4_REVEL) * f_B_val) if not pd.isna(f_P_val) else np.nan

                # PrimateAI-3D_MT
                score_val = parts[header_dict['PrimateAI-3D_MT']]
                #print(f"Gene {gene_id}, PrimateAI-3D_MT score {score_val}")
                f_P_val = beta_mixture_PDF(float(score_val), *ClinVar_path_params['PrimateAI-3D_MT']) if score_val != 'NA' else np.nan
                f_B_val = 1.0
                PP1_PAI = PPost1_PrimateAI_dict.get(gene_id, np.nan)
                p2_PAI = PPost2_priors['PrimateAI-3D_MT']
                PP2_PAI = (p2_PAI * f_P_val) / (p2_PAI * f_P_val + (1 - p2_PAI) * f_B_val) if not pd.isna(f_P_val) else np.nan
                p3_PAI = PPost3_PrimateAI_dict.get(gene_id, np.nan)
                PP3_PAI = (p3_PAI * f_P_val) / (p3_PAI * f_P_val + (1 - p3_PAI) * f_B_val) if not pd.isna(f_P_val) else np.nan
                p4_PAI = PPost4_PrimateAI_dict.get(gene_id, np.nan)
                PP4_PAI = (p4_PAI * f_P_val) / (p4_PAI * f_P_val + (1 - p4_PAI) * f_B_val) if not pd.isna(f_P_val) else np.nan

                pp_missense_list = [PP1_AM, PP2_AM, PP3_AM, PP4_AM,
                                    PP1_REVEL, PP2_REVEL, PP3_REVEL, PP4_REVEL,
                                    PP1_PAI, PP2_PAI, PP3_PAI, PP4_PAI]
                
                pp_missense_list = [round(x, round_until) if not pd.isna(x) else 'NA' for x in pp_missense_list]

                outfile.write('\t'.join(parts + [str(x) for x in pp_missense_list]) + '\n')
                continue

            # Non-LoF, non-missense variants
            PP_vals = ["NA"] * len(add_header)
            outfile.write('\t'.join(parts + [str(x) for x in PP_vals]) + '\n')


            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculate diagnostic probabilities for variants using reference cohort data.'
    )

    parser.add_argument('--gene-list-file',    default="gene_list_file.txt", help='Whitelist gene list file (one ENS ID per line), default: gene_list_file.txt')
    parser.add_argument('--ref-var-info-file', required=True, help='Reference cohort variant info (print_var_info output)')
    parser.add_argument('--ref-sum-stats-file',required=True, help='Reference cohort per-gene summary stats file')
    parser.add_argument('--processed-vars-file',required=True, help='File with variants to be processed (can be gzipped)')
    parser.add_argument('--output-file',       required=True, help='Path to write the results (gzipped if .gz)')
    parser.add_argument('--colname',           default='DiagProb', help='Prefix for added columns (default: DiagProb)')
    parser.add_argument('--pnull',             type=float, default=0.25, help='Probability of LoF misannotation (default: 0.25)')
    parser.add_argument('--P_D',               type=float, default=0.01, help='Disease prevalence (default: 0.01)')
    parser.add_argument('--verbose',           action='store_true', help='Enable verbose logging')
    parser.add_argument('--add-one-to-varN',    action='store_true', help='Whether to add one to the observed variant count (varN) in enrichment calculations (default: False)')

    args = parser.parse_args()

    # Print values of all arguments for logging
    print("Running with the following arguments:")
    for arg in vars(args):        print(f"  {arg}: {getattr(args, arg)}")

    main_func(
        gene_list_file=args.gene_list_file,
        ref_var_info_file=args.ref_var_info_file,
        ref_sum_stats_file=args.ref_sum_stats_file,
        processed_vars_file=args.processed_vars_file,
        output_file=args.output_file,
        colname=args.colname,
        pnull=args.pnull,
        verbose=args.verbose,
        P_D=args.P_D,
        add_one_to_varN=args.add_one_to_varN
    )              

