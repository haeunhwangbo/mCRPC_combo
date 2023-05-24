configfile: "config.yaml"

import pandas as pd
import glob
import os
import numpy as np


# DIRECTORIES
FIG_DIR = config['fig_dir']
TABLE_DIR = config['table_dir']
EXPERIMENTAL_DATA_DIR = config['experimental_dir']

def get_trial_list(dataset):
    trials_df = pd.read_csv(config[dataset]['metadata_sheet'], 
                            sep='\t', header=0)
    if dataset == 'placebo':
        trial_list = list(placebo_df['File prefix'])
    else:
        trials_set = set()
        for i in range(trials_df.shape[0]):
            for arm in ['Experimental', 'Control']:  
                t = trials_df.at[i, arm]
                trials_set.add(t)
        trial_list = list(trials_set)
    return trial_list

def get_raw_files(dataset):
    dir = config[dataset]['raw_dir']
    trial_list = get_trial_list(dataset)
    file_list = [f"{dir}/{trial}.csv" for trial in trial_list]
    return file_list

def get_trial_files(dataset):
    dir = config[dataset]['data_dir']
    trial_list = get_trial_list(dataset)
    file_list = [f"{dir}/{trial}.clean.csv" for trial in trial_list]
    return file_list

def get_pred_list(dataset):
    trials_df = pd.read_csv(config[dataset]['metadata_sheet'], 
                            sep='\t', header=0)
    pred_list = []
    for i in range(trials_df.shape[0]):
        name_a = trials_df.at[i, 'Experimental']
        name_b = trials_df.at[i, 'Control']
        pred_list.append(f'{name_a}-{name_b}')
    return pred_list


# FILE LISTS
PFS_TRIALS = expand("{output_file}", output_file=get_trial_files("PFS"))
rPFS_TRIALS = expand("{output_file}", output_file=get_trial_files("rPFS"))
PSA_TRIALS = expand("{output_file}", output_file=get_trial_files("waterfall"))

PFS_PRED_FILES =  expand(f"{config['PFS']['pred_dir']}/{{pred}}_combination_predicted_{{model}}.csv", 
                         pred=get_pred_list("PFS"), model=['ind'])
rPFS_PRED_FILES =  expand(f"{config['rPFS']['pred_dir']}/{{pred}}_combination_predicted_{{model}}.csv", 
                         pred=get_pred_list("rPFS"), model=['ind'])
PSA_PRED_FILES =  expand(f"{config['waterfall']['pred_dir']}/{{pred}}_combination_predicted_{{model}}.csv", 
                         pred=get_pred_list("waterfall"), model=['ind'])

rule all:
    input:
        f"{config['PFS']['table_dir']}/PFS_predictive_power.csv",
        f"{config['rPFS']['table_dir']}/rPFS_predictive_power.csv",
        f"{config['PFS']['fig_dir']}/PFS_survival_plots.pdf",
        f"{config['rPFS']['fig_dir']}/rPFS_survival_plots.pdf",
        f"{config['waterfall']['fig_dir']}/waterfall_survival_plots.pdf",
        f'{FIG_DIR}/CTRPv2_corr_distributions.pdf',
        f'{TABLE_DIR}/experimental_correlation_report.csv'

rule preprocess:
    input:
        f"{config['PFS']['metadata_sheet']}",
        f"{config['rPFS']['metadata_sheet']}",
        f"{config['waterfall']['metadata_sheet']}",
        expand("{input_file}", input_file=get_raw_files("PFS")),
        expand("{input_file}", input_file=get_raw_files("rPFS")),
        expand("{input_file}", input_file=get_raw_files("waterfall")),
        "src/preprocessing.py"
    output:
        PFS_TRIALS,
        rPFS_TRIALS,
        PSA_TRIALS
    conda:
        "env/environment_short.yml"
    shell:
        "python src/preprocessing.py PFS; "
        "python src/preprocessing.py rPFS; "
        "python src/preprocessing.py waterfall"


rule find_seeds:
    input:
        f"{config['PFS']['metadata_sheet']}",
        f"{config['rPFS']['metadata_sheet']}",
        f"{config['waterfall']['metadata_sheet']}",
        PFS_TRIALS,
        rPFS_TRIALS,
        PSA_TRIALS,
        "src/find_median_sim.py"
    output:
        f"{config['PFS']['metadata_sheet_seed']}",
        f"{config['rPFS']['metadata_sheet_seed']}",
        f"{config['waterfall']['metadata_sheet_seed']}",
    shell:
        "python src/find_median_sim.py PFS; "
        "python src/find_median_sim.py rPFS; "
        "python src/find_median_sim.py waterfall"

rule hsa_prediction:
    input:
        f"{config['PFS']['metadata_sheet_seed']}",
        f"{config['rPFS']['metadata_sheet_seed']}",
        f"{config['waterfall']['metadata_sheet_seed']}",
        PFS_TRIALS,
        rPFS_TRIALS,
        PSA_TRIALS,
        "src/hsa_additivity_model.py"
    output:
        PFS_PRED_FILES,
        rPFS_PRED_FILES,
        PSA_PRED_FILES
    shell:
        "python src/hsa_additivity_model.py PFS; "
        "python src/hsa_additivity_model.py rPFS; "
        "python src/hsa_additivity_model.py waterfall"

rule survival_plots:
    input:
        f"{config['PFS']['metadata_sheet_seed']}",
        f"{config['rPFS']['metadata_sheet_seed']}",
        f"{config['waterfall']['metadata_sheet_seed']}",
        PFS_TRIALS,
        rPFS_TRIALS,
        PSA_TRIALS,
        PFS_PRED_FILES,
        rPFS_PRED_FILES,
        PSA_PRED_FILES,
        "src/plotting/plot_survival_curves_suppl.py"
    output:
        f"{config['PFS']['fig_dir']}/PFS_survival_plots.pdf",
        f"{config['rPFS']['fig_dir']}/rPFS_survival_plots.pdf",
        f"{config['waterfall']['fig_dir']}/waterfall_survival_plots.pdf"
    shell:
        "python src/plotting/plot_survival_curves_suppl.py PFS; "
        "python src/plotting/plot_survival_curves_suppl.py rPFS; "
        "python src/plotting/plot_survival_curves_suppl.py waterfall"


rule experimental_correlation:
    input:
        f'{EXPERIMENTAL_DATA_DIR}/CTRPv2_CCL.csv',
        f'{EXPERIMENTAL_DATA_DIR}/CTRPv2_drug.csv',
        f'{EXPERIMENTAL_DATA_DIR}/Recalculated_CTRP_12_21_2018.txt',
        f'{EXPERIMENTAL_DATA_DIR}/CTRPv2_clincal_active_drug_pairwise_corr.csv',
    output:
        f'{FIG_DIR}/CTRPv2_corr_distributions.pdf',
        f'{TABLE_DIR}/experimental_correlation_report.csv'
    script:
        "src/experimental_correlation.py"


rule predictive_power:
    input:
        f"{config['PFS']['metadata_sheet_seed']}",
        f"{config['rPFS']['metadata_sheet_seed']}",
        PFS_TRIALS,
        rPFS_TRIALS,
        PFS_PRED_FILES,
        rPFS_PRED_FILES
    output:
        f"{config['PFS']['table_dir']}/PFS_predictive_power.csv",
        f"{config['rPFS']['table_dir']}/rPFS_predictive_power.csv",
    shell:
        "python src/predictive_power.py PFS; "
        "python src/predictive_power.py rPFS; "
