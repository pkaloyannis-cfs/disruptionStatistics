# -*- coding: utf-8 -*-
"""
Created on Wed May 18 10:27:59 2022

@author: yalynskayn
"""

from TCV_CMod_1D2D_new_data import *


def table_disr(machine_param1, machine_param2, X_range, Y_range, jj_50, denom_idx, deltat, name, exp, param_name1, param_name2):
    X, Y = machine_param1, machine_param2
    H1, xedges1, yedges1 = np.histogram2d(
        X[jj_50], Y[jj_50], bins=25, range=[X_range, Y_range])
    # denminator
    H2, xedges2, yedges2 = np.histogram2d(
        X[denom_idx], Y[denom_idx], bins=25, range=[X_range, Y_range])
    total_sum = np.sum(H1) + np.sum(H2)
    weight = (H1 + H2) / total_sum
    # probability = heatmap1/total_sum/deltat
    vps = (H1 / H2) * weight / deltat
    # changing NaNs and infinities to 0.0
    mask = np.where(~np.isfinite(vps))
    vps[mask] = 0.0
    df = pd.DataFrame({f'{name}_{exp}': vps.max()}, [f'{param_name1}, {param_name2}'])
    return df


df = table_disr(all_data.LI, all_data.BETAN, LI_range, BETAN_range, jj_50, denom_idx, deltat, name, exp,
                'Li', 'betan')
df = df.append(table_disr(all_data.KAPPA, all_data.DELTA, KAPPA_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
                          'elongation', 'average_triangularity'))
df = df.append(table_disr(all_data.KAPPA, all_data.DELTA_TOP, KAPPA_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
                          'elongation', 'triangularity_top'))
df = df.append(table_disr(all_data.KAPPA, all_data.DELTA_BOTTOM, KAPPA_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
                          'elongation', 'triangularity_bottom'))
df = df.append(table_disr(all_data.KAPPA, all_data.GAP_in, KAPPA_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
                          'elongation', 'GAP_in'))
df = df.append(table_disr(all_data.KAPPA, all_data.GAP_out, KAPPA_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
                          'elongation', 'GAP_out'))
df = df.append(table_disr(all_data.KAPPA, all_data.Q95, KAPPA_range, Q95_range, jj_50, denom_idx, deltat, name, exp,
                          'elongation', 'Q95'))
df = df.append(table_disr(all_data.KAPPA, all_data.I_P, KAPPA_range, I_P_range, jj_50, denom_idx, deltat, name, exp,
                          'elongation', 'I_P'))
df = df.append(table_disr(all_data.BETAN, all_data.GAP_in, BETAN_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
                          'betan', 'GAP_in'))
df = df.append(table_disr(all_data.BETAN, all_data.GAP_out, BETAN_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
                          'betan', 'GAP_out'))
df = df.append(table_disr(all_data.BETAN, all_data.Q95, BETAN_range, Q95_range, jj_50, denom_idx, deltat, name, exp,
                          'betan', 'Q95'))
df = df.append(table_disr(all_data.BETAN, all_data.I_P, BETAN_range, I_P_range, jj_50, denom_idx, deltat, name, exp,
                          'betan', 'I_P'))
df = df.append(table_disr(all_data.LI, all_data.GAP_in, LI_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
                          'Li', 'GAP_in'))
df = df.append(table_disr(all_data.LI, all_data.GAP_out, LI_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
                          'Li', 'GAP_out'))
df = df.append(table_disr(all_data.LI, all_data.Q95, LI_range, Q95_range, jj_50, denom_idx, deltat, name, exp,
                          'Li', 'Q95'))
df = df.append(table_disr(all_data.LI, all_data.I_P, LI_range, I_P_range, jj_50, denom_idx, deltat, name, exp,
                          'Li', 'I_P'))
df = df.append(table_disr(all_data.GAP_out, all_data.GAP_in, GAP_out_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
                          'GAP_out', 'GAP_in'))
df = df.append(table_disr(all_data.I_P, all_data.Q95, I_P_range, Q95_range, jj_50, denom_idx, deltat, name, exp,
                          'I_P', 'Q95'))
df = df.append(table_disr(all_data.I_P, all_data.GAP_in, I_P_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
                          'I_P', 'GAP_in'))
df = df.append(table_disr(all_data.I_P, all_data.GAP_out, I_P_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
                          'I_P', 'GAP_out'))
df = df.append(table_disr(all_data.Q95, all_data.GAP_in, Q95_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
                          'Q95', 'GAP_in'))
df = df.append(table_disr(all_data.Q95, all_data.GAP_out, Q95_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
                          'Q95', 'GAP_out'))
df = df.append(table_disr(all_data.I_P, all_data.DELTA, I_P_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
                          'I_P', 'average_triangularity'))
df = df.append(table_disr(all_data.I_P, all_data.DELTA_TOP, I_P_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
                          'I_P', 'triangularity_top'))
df = df.append(table_disr(all_data.I_P, all_data.DELTA_BOTTOM, I_P_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
                          'I_P', 'triangularity_bottom'))
df = df.append(table_disr(all_data.Q95, all_data.DELTA, Q95_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
                          'Q95', 'average_triangularity'))
df = df.append(table_disr(all_data.Q95, all_data.DELTA_TOP, Q95_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
                          'Q95', 'triangularity_top'))
df = df.append(table_disr(all_data.Q95, all_data.DELTA_BOTTOM, Q95_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
                          'Q95', 'triangularity_bottom'))
df = df.append(table_disr(all_data.GAP_in, all_data.DELTA, GAP_in_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
                          'GAP_in', 'average_triangularity'))
df = df.append(table_disr(all_data.GAP_in, all_data.DELTA_TOP, GAP_in_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
                          'GAP_in', 'triangularity_top'))
df = df.append(table_disr(all_data.GAP_in, all_data.DELTA_BOTTOM, GAP_in_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
                          'GAP_in', 'triangularity_bottom'))
df = df.append(table_disr(all_data.GAP_out, all_data.DELTA, GAP_out_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
                          'GAP_out', 'average_triangularity'))
df = df.append(table_disr(all_data.GAP_out, all_data.DELTA_TOP, GAP_out_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
                          'GAP_out', 'triangularity_top'))
df = df.append(table_disr(all_data.GAP_out, all_data.DELTA_BOTTOM, GAP_out_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
                          'GAP_out', 'triangularity_bottom'))
df.to_csv(f"disrupt_per_sec{name}_{exp}.csv")

ans = input("Do you have all 3 tables for CMod and TCV? (--> 'no', 'yes) If 'no', then restart the program with inputs for missing data\n")

if ans == 'yes':
    df1 = pd.read_csv("disrupt_per_sec_CMod_VDE.csv")
    df2 = pd.read_csv("disrupt_per_sec_TCV_VD.csv")
    df3 = pd.read_csv("disrupt_per_sec_TCV_VDE.csv")
    df = pd.concat([df1, df2['_TCV_VD'], df3['_TCV_VDE']], axis=1)
    #df = pd.read_csv('DISRUPT_TABLE.csv')
    df = df.rename(columns={"Unnamed: 0": "parametrs"})
    df = df.rename(columns={"_CMod": "CMod", "_TCV_VD": "TCV_VD", "_TCV_VDE": "TCV_VDE"})
    parametrs = df.parametrs

    df.plot.barh(x='parametrs', y={'CMod', 'TCV_VD', 'TCV_VDE'}, figsize=(6, 8))

elif ans == 'no':
    print('you answered no, try with other data')

else:
    print('wrong answer')
