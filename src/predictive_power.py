import numpy as np
import pandas as pd
import yaml
import warnings
import argparse
from multiprocessing import Pool
from coxhazard_test import get_cox_results, create_ipd
warnings.filterwarnings("ignore")

with open('config.yaml', 'r') as f:
    CONFIG = yaml.safe_load(f)

N = 5000
NRUN = 1000
SEED = 0
rng = np.random.default_rng(SEED)


def simulate_one_trial(sampled_patients: np.array, ipd_ori: pd.DataFrame, ipd_control: pd.DataFrame) -> int:
    ipd_sim = ipd_ori.reindex(sampled_patients)
    p, HR, low95, high95 = get_cox_results(ipd_control, ipd_sim)
    success = 0
    if p < 0.05 and high95 < 1:
        success = 1
    return success


def calculate_success_prob(input_df: pd.DataFrame, i: int, data_dir: str, pred_dir: str) -> tuple:
    """Calculate the probability of success for a given trial.

    Args:
        input_df (pd.DataFrame): _description_
        i (int): _description_
        data_dir (str): _description_
        pred_dir (str): _description_

    Returns:
        tuple: _description_
    """
    name_a = input_df.at[i, 'Experimental']
    name_b = input_df.at[i, 'Control']

    n_combo = max(input_df.at[i, 'N_control'].astype(int), 
                  input_df.at[i, 'N_experimental'].astype(int))

    # import prediction
    independent = pd.read_csv(
        f'{pred_dir}/{name_a}-{name_b}_combination_predicted_ind.csv')

    independent = independent[independent['Time']
                              < independent['Time'].max() - 0.1]
    # make ipd using a large N and sample smaller number of patients later
    ipd_ind = create_ipd(independent, n=N)
    
    success_prob_ctrl = 0
    success_prob_exp = 0
    for arm in ['Control', 'Experimental']:
        name_base = input_df.at[i, arm]
        n_base = input_df.at[i, 'N_control'].astype(int)
        try:
            ipd_base = pd.read_csv(f'{data_dir}/{name_base}_indiv.csv')

        except FileNotFoundError:
            df_base = pd.read_csv(f'{data_dir}/{name_base}.clean.csv')
            ipd_base = create_ipd(df_base, n=n_base)
        
        # calculate probability of success using n_combo patients
        success_cnt = 0
        for run in range(NRUN):
            sampled_patients = rng.integers(0, N, n_combo)
            success_cnt += simulate_one_trial(sampled_patients, ipd_ind, ipd_base)
        
        # Cox-PH test usign large N
        p, hr, lower, upper = get_cox_results(ipd_base, ipd_ind)

        if arm == 'Control':
            success_prob_ctrl = success_cnt / NRUN
            p_ctrl, hr_ctrl, lower_ctrl, upper_ctrl = p, hr, lower, upper
        else:
            success_prob_exp = success_cnt / NRUN
            p_exp, hr_exp, lower_exp, upper_exp = p, hr, lower, upper
    
    return i, success_prob_exp, success_prob_ctrl, p_ctrl, hr_ctrl, lower_ctrl, upper_ctrl, p_exp, hr_exp, lower_exp, upper_exp


def predictive_power(metadata, data_dir, pred_dir):
    outdf = metadata.copy()
    ll = []
    with Pool(processes=4) as pool:
        args_list = [(metadata, i, data_dir, pred_dir) for i in range(metadata.shape[0])]
        for result in pool.starmap(calculate_success_prob, args_list):
            print(result)
            ll.append(result)
    tmp = pd.DataFrame(ll, columns=['idx','prob_success_exp', 'prob_success_ctrl', 
                                    'p_ctrl', 'hr_ctrl', 'lower_ctrl', 'upper_ctrl',
                                    'p_exp', 'hr_exp', 'lower_exp', 'upper_exp'])
    tmp = tmp.set_index('idx', drop=True)
    outdf = pd.concat([outdf, tmp], axis=1)
    return outdf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset', type=str, 
                        help='Dataset to use (PFS, rPFS, waterfall)')
    args = parser.parse_args()
    config_dict = CONFIG[args.dataset]
    metadata = pd.read_csv(config_dict['metadata_sheet_seed'], sep='\t')
    data_dir = config_dict['data_dir']
    pred_dir = config_dict['pred_dir']
    table_dir = config_dict['table_dir']

    outdf = predictive_power(metadata, data_dir, pred_dir)
    outdf.to_csv(f'{table_dir}/{args.dataset}_predictive_power.csv', index=False)


if __name__ == '__main__':
    main()
