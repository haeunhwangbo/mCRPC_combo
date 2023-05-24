import pandas as pd
import numpy as np
from multiprocessing import Pool
import argparse
import yaml
import tempfile
import os
from hsa_additivity_model import predict_hsa

with open('config.yaml', 'r') as f:
    CONFIG = yaml.safe_load(f)

NRUN = 100

def make_prediction_for_each_combo(i: int, indf: pd.DataFrame, data_dir: str, pred_dir: str):
    name_a = indf.at[i, 'Experimental']
    name_b = indf.at[i, 'Control']
    corr = indf.at[i, 'Corr']  # experimental spearman correlation value

    df_a = pd.read_csv(f'{data_dir}/{name_a}.clean.csv',
                    header=0, index_col=False)
    df_b = pd.read_csv(f'{data_dir}/{name_b}.clean.csv',
                    header=0, index_col=False)

    for k in range(NRUN):
        seed = k
        ind = predict_hsa(df_a, df_b, name_a, name_b,
                          waterfall=waterfall,
                          rho=corr,
                          seed_ind=seed,
                          save=False)
        ind.to_csv(f'{pred_dir}/{name_a}-{name_b}_combination_predicted_ind_run{seed:02d}.csv')


def make_predictions_diff_seeds(indf: pd.DataFrame, data_dir: str, pred_dir: str, waterfall=False):
    args_list = [(i, indf, data_dir, pred_dir, waterfall) for i in indf.index]
    with Pool(processes=8) as pool:
        pool.starmap(make_prediction_for_each_combo, args_list)


def find_median_sim(indf: pd.DataFrame, pred_dir: str, save=True, outfile=None) -> pd.DataFrame:
    med_df = indf.copy()
    med_df.loc[:, 'ind_median_std'] = np.nan
    med_df.loc[:, 'ind_median_run'] = 0

    for i in range(indf.shape[0]):
        name_a = indf.at[i, 'Experimental']
        name_b = indf.at[i, 'Control']
        ind_arr = np.zeros(NRUN)
        for seed in range(NRUN):
            ind = pd.read_csv(f'{pred_dir}/{name_a}-{name_b}_combination_predicted_ind_run{seed:02d}.csv')
            ind_arr[seed] = ind.loc[2499:2500, 'Time'].mean()
        med_df.loc[i, 'ind_median_std'] = np.std(ind_arr)
        # save run# of the median
        ind_idx = np.argsort(ind_arr)[len(ind_arr)//2]
        med_df.loc[i, 'ind_median_run'] = ind_idx

    if save:
        med_df.to_csv(outfile, index=False, sep='\t')
    
    return med_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset', type=str,
                        help='Dataset to use')
    args = parser.parse_args()

    table_dir = CONFIG['table_dir']
    config_dict = CONFIG[args.dataset]
    sheet = config_dict['metadata_sheet']
    data_dir = config_dict['data_dir']
    outfile = config_dict['metadata_sheet_seed']
    
    indf = pd.read_csv(sheet, sep='\t')

    with tempfile.TemporaryDirectory(dir=table_dir) as temp_dir:
        make_predictions_diff_seeds(indf, data_dir, temp_dir)
        find_median_sim(indf, temp_dir, save=True, outfile=outfile)
