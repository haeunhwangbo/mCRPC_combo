import numpy as np
import pandas as pd
from coxhazard_test import get_cox_results, create_ipd
from utils import interpolate
from plotting.plot_utils import import_input_data
from pathlib import Path


def hsa_add_diff():
    """Compute difference between additivity and HSA by calculating
    normalized difference between the two curves and performing
    Cox-PH model between the two survival curves for each combination.

    Returns:
        pd.DataFrame: difference between HSA and additivity
    """    
    indir, cox_df = import_input_data()
    diff_df = pd.DataFrame(index=cox_df.index, columns=[
                           'HSA - Control', 'Additivity - HSA',
                           'p', 'HR', 'HRlower', 'HRupper'])
    n = 5000
    for i in range(cox_df.shape[0]):
        name_a = cox_df.at[i, 'Experimental']
        name_b = cox_df.at[i, 'Control']
        path = cox_df.at[i, 'Path'] + '/'

        # import data
        control = pd.read_csv(path + f'{name_b}.clean.csv')
        independent = pd.read_csv(
            indir + f'{name_a}-{name_b}_combination_predicted_ind.csv')
        additive = pd.read_csv(
            indir + f'{name_a}-{name_b}_combination_predicted_add.csv')
        tmax = np.amin(
            [control['Time'].max(), independent['Time'].max(), additive['Time'].max()])

        tmax = np.amin([additive['Time'].max(), independent['Time'].max()])
        independent = independent[independent['Time'] < tmax]
        additive = additive[additive['Time'] < tmax]

        ipd_add = create_ipd(additive, n=n)
        ipd_ind = create_ipd(independent, n=n)

        # set up data for Cox regress
        p, hr, hr_lower, hr_upper = get_cox_results(ipd_add, ipd_ind)
        diff_df.at[i, 'p'] = p
        diff_df.at[i, 'HR'] = hr
        diff_df.at[i, 'HRlower'] = hr_lower
        diff_df.at[i, 'HRupper'] = hr_upper

        f_ctrl = interpolate(control, x='Time', y='Survival')
        f_ind = interpolate(independent, x='Time', y='Survival')
        f_add = interpolate(additive, x='Time', y='Survival')

        timepoints = np.linspace(0, tmax, n)

        diff_df.at[i,
                   'HSA - Control'] = sum(f_ind(timepoints) - f_ctrl(timepoints))/n
        diff_df.at[i,
                   'Additivity - HSA'] = sum(f_add(timepoints) - f_ind(timepoints))/n

    return diff_df.astype(np.float64)


def added_benefit_hsa_add_syn():
    """Compute normalized difference between survival curves. 
    Combo - Control, HSA - Control, Additivity - HSA, Combo - Additivity.

    Returns:
        pd.DataFrame: result dataframe
    """    
    indir, cox_df = import_input_data()
    diff_df = pd.DataFrame(index=cox_df.index, 
                           columns=['Combo - Control', 
                                    'HSA - Control', 
                                    'Additivity - HSA', 
                                    'Combo - Additivity'])
    diff_df.loc[:, 'Model'] = cox_df['Model']
    diff_df.loc[:, 'Combination'] = cox_df['Combination']

    n = 5000
    for i in range(cox_df.shape[0]):
        name_a = cox_df.at[i, 'Experimental']
        name_b = cox_df.at[i, 'Control']
        name_ab = cox_df.at[i, 'Combination']
        path = cox_df.at[i, 'Path'] + '/'

        # import data
        control = pd.read_csv(path + f'{name_b}.clean.csv')
        obs_combo = pd.read_csv(path + f'{name_ab}.clean.csv')
        independent = pd.read_csv(
            indir + f'{name_a}-{name_b}_combination_predicted_ind.csv')
        additive = pd.read_csv(
            indir + f'{name_a}-{name_b}_combination_predicted_add.csv')
        tmax = np.amin(
            [control['Time'].max(), independent['Time'].max(), additive['Time'].max()])
        f_ctrl = interpolate(control, x='Time', y='Survival')
        f_ind = interpolate(independent, x='Time', y='Survival')
        f_add = interpolate(additive, x='Time', y='Survival')
        f_obs = interpolate(obs_combo, x='Time', y='Survival')
        timepoints = np.linspace(0, tmax, n)

        diff_df.at[i,
                'Combo - Control'] = sum(f_obs(timepoints) - f_ctrl(timepoints))/n
        diff_df.at[i,
                'HSA - Control'] = sum(f_ind(timepoints) - f_ctrl(timepoints))/n
        diff_df.at[i,
                'Additivity - HSA'] = sum(f_add(timepoints) - f_ind(timepoints))/n
        diff_df.at[i, 'Combo - Additivity'] = sum(
            f_obs(timepoints) - f_add(timepoints))/n
        
        diff_df.round(5).to_csv('../analysis/additivity_HSA_similarity/difference.csv', index=False)
    return diff_df


if __name__ == '__main__':
    #hsa_add_diff()
    outdir = '../analysis/additivity_HSA_similarity/'
    new_directory = Path(outdir)
    new_directory.mkdir(parents=True, exist_ok=True)
    added_benefit_hsa_add_syn()
