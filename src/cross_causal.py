# pip install causalinference

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# causal inference
from causalinference import CausalModel
import multiprocessing
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import QuantileTransformer


def read_data_TCV(exp):
    # concatenate ALL the Table file, might be useful for the future
    # all_files=sorted(glob('Lite_Table_TCV_*.csv'))
    # defines VD or VDE shotlist from the master table according to the "Flag" parameter
    print('Reading TCV parameters file...')
    mtable = pd.read_excel('TCV_parameters.xlsx', sheet_name='TCV Disr > 50kA')
    if exp == 'VD':
        print('Reading files with VD, other disruptions and non disruptive data...')
        Lite_Table_TCV_Disruptions_VD = pd.read_csv(
            'Lite_Table_TCV_Disruptions_VD_new.csv', index_col=False)
        Lite_Table_TCV_Disruptions_0 = pd.read_csv(
            'Lite_Table_TCV_Disruptions_0_new.csv', index_col=False)
        Lite_Table_TCV_NonDisruptions = pd.read_csv(
            'Lite_Table_TCV_NonDisruptions_new.csv', index_col=False)
        all_files = sorted(('Lite_Table_TCV_Disruptions_VD_new.csv', 'Lite_Table_TCV_Disruptions_0_new.csv',
                            'Lite_Table_TCV_NonDisruptions_new.csv'))
        tcv_all = pd.concat((pd.read_csv(file).assign(filename=file)
                             for file in all_files), ignore_index=True)
        shots_list = mtable[(mtable['Flag']) == 'VD'].shot.values
        np.savetxt("vdshots_tcv.txt", shots_list, fmt="%s")
    if exp == 'VDE':
        print('Reading files with VDE, other disruptions and non disruptive data...')
        all_files = sorted(('Lite_Table_TCV_Disruptions_VDE_new.csv', 'Lite_Table_TCV_Disruptions_0_new.csv',
                            'Lite_Table_TCV_NonDisruptions_new.csv'))
        tcv_all = pd.concat((pd.read_csv(file).assign(filename=file)
                             for file in all_files), ignore_index=True)
        shots_list = mtable[(mtable['Flag']) == 'VDE'].shot.values
    return tcv_all, shots_list


def define_indices_TCV(tcv_all):
    print('Definig indices...')
    # Define indices
    indices_no_disrupt = list(
        tcv_all.loc[tcv_all['filename'] == 'Lite_Table_TCV_NonDisruptions_new.csv'].index)
    indices_disrupt = list(
        tcv_all.loc[tcv_all['filename'] != 'Lite_Table_TCV_NonDisruptions_new.csv'].index)
    indices_flattop = list(tcv_all.index)
    # find indices of the rows when exactly the disruption happend, getting rid of all the dublicated times
    indices_disrupt_time = list(
        tcv_all.Shot[indices_disrupt].index[~tcv_all.Shot[indices_disrupt].duplicated()])
    return indices_no_disrupt, indices_disrupt, indices_flattop, indices_disrupt_time


def read_data_CMod():
    print('Reading CMod data...')
    cmod = pd.read_csv('new_data_CMod_VDE.csv')
    # read shotlist
    VDEshotlist = pd.read_csv('VDE_shotlist.txt')
    VDEshotlist.to_csv('VDE_shotlist.csv')
    vdeshots = pd.read_csv('VDE_shotlist.csv', index_col=None, header=None)
    return cmod, vdeshots


def define_indices_CMod(cmod):
    # Define indices
    # Define indices of missing objects
    print('Definig indices...')
    indices_no_disrupt = list(cmod[cmod.Time2Disr.isnull()].index)
    indices_disrupt = list(cmod[~cmod.Time2Disr.isnull()].index)
    # Define indices of disruption happened
    indices_flattop = list(cmod[abs(cmod.dipprog_dt) <= 1.e3].index)
    indices_disrupt_time = list(cmod[cmod.Time2Disr <= 3e-3].index)
    return indices_no_disrupt, indices_disrupt, indices_flattop, indices_disrupt_time


def indices_operation(machine_param, shots_list, indices_no_disrupt, indices_disrupt, indices_flattop,
                      indices_disrupt_time, name, exp):
    print('Operations with indices are in progress...')
    # find indices of the rows when ecactly the disruption happend in the flattop data, same as indices_disrupt_time
    indices_disrupt_time_in_flattop = list(
        np.intersect1d(indices_disrupt_time, indices_flattop))
    # find indices of the whole shot numbers which end in disruption
    indices_disrupt_in_flattop = list(machine_param.Shot[np.isin(machine_param.Shot, machine_param.Shot[
        indices_disrupt_time_in_flattop])].index)
    # find indices of the whole shot numbers in the database which end in disruption in the flattop data
    indices_flattop_disrupt_in_flattop = list(
        np.intersect1d(indices_flattop, indices_disrupt_in_flattop))
    # find indices of non disruptive shots
    indices_flattop_no_disrupt = list(np.intersect1d(indices_flattop, indices_no_disrupt))

    # for ease of reference we named them in this way
    ii = indices_flattop_no_disrupt
    ff = indices_flattop_disrupt_in_flattop
    if name == '_TCV':
        # reading other disruptions daa file
        other_disruptions = list(
            machine_param.loc[machine_param['filename'] == 'C://Users//nadya//Documents//Python Scripts//causal_inf_cross//Lite_Table_TCV_Disruptions_0_new.csv'].index)
        # all VDE/VD data reading
        if exp == 'VD':
            all_exp_file = machine_param.loc[machine_param['filename'] ==
                                             'C://Users//nadya//Documents//Python Scripts//causal_inf_cross//Lite_Table_TCV_Disruptions_VD_new.csv']
        if exp == 'VDE':
            all_exp_file = machine_param.loc[machine_param['filename'] ==
                                             'C://Users//nadya//Documents//Python Scripts//causal_inf_cross//Lite_Table_TCV_Disruptions_VDE_new.csv']
        # jj - flattop VDE/VD disruptions that are in vd(e)shotlist, idx is choosing Shots from VDE/VD data which are in the vd(e)shots with Ip>50kA
        idx = list(all_exp_file.Shot[np.isin(all_exp_file.Shot, shots_list)].index)
        # subselect flattop phase of disruptive data, part of shotlist
        jj = list(np.intersect1d(ff, idx))
        # all data is flattop, so other_disruptions=other_disruptions_flattop
        other_disruptions_flattop = other_disruptions
        # include data for VDEs flattop where time_until_disrupt<=40ms
        jj_before_40 = list(np.intersect1d(
            jj, list(all_exp_file[all_exp_file.Time2Disr <= 0.04].index)))
        # isolate for VDE shots behavior around 50ms before disruption
        jj_50 = list(np.intersect1d(jj, all_exp_file.loc[
            (all_exp_file['Time2Disr'] >= 0.04) & (all_exp_file['Time2Disr'] <= 0.06)].index))

    if name == '_CMod':
        # find indices in shot variable that coincide with shot numbers in shotlist
        idx = list(machine_param.Shot[np.isin(machine_param.Shot, shots_list)].index)
        # subselect flattop phase of disruptive data, part of VDE shotlist
        jj = list(np.intersect1d(ff, idx))
        # include data for VDEs flattop where time_until_disrupt<=40ms
        jj_before_40 = list(np.intersect1d(
            jj, list(machine_param[machine_param.Time2Disr <= 0.04].index)))
        # isolate for VDE shots behavior around 50ms before disruption
        jj_50 = list(np.intersect1d(jj, machine_param.loc[
            (machine_param['Time2Disr'] >= 0.04) & (machine_param['Time2Disr'] <= 0.06)].index))
        # all flattop disruptions
        disruption = list(machine_param.Shot[ff])
        # shotlist of all other flattop disruptions that are not VDEs
        shotlist_disruptions = np.unique(disruption)
        other_disruptions = list(
            shotlist_disruptions[~np.isin(shotlist_disruptions, shots_list)])
        # indices of all flattop disruptions that are not VDEs
        other_disruptions_flattop = list(
            machine_param.Shot[np.isin(machine_param.Shot, other_disruptions)].index)
        # flattop indices of all flattop disruptions that are not VDEs
        # jj2 = (np.intersect1d(ff,other_disruptions_flattop)).tolist()
        # since the whole data is flattop all_other_disruptions is equal to jj2 which is flattop not vde disruptions
    # concatenate indices in the denominator:
    #denom_idx = np.concatenate((ii, other_disruptions_flattop, jj_before_40), axis=0)
    return jj_50, ii


work_type = input(
    'Which installation do you want to choose?\n - Input 1 if TCV \n - Input 2 if CMod \ -->')
if work_type == '1':
    # function that drops some columns to fasten the data procession and saves the results in the new "Lite" versions
    # fasten_TCV()
    name = '_TCV'
    # takes shots with VD or VDE flag from TCV_parameters.xlsx
    exp = input('Do you want to study VD or VDE? (input VD or VDE) -->').upper()
    all_data, shots_list = read_data_TCV(exp)
    indices_no_disrupt, indices_disrupt, indices_flattop, indices_disrupt_time = define_indices_TCV(
        all_data)
    jj_50, ii = indices_operation(all_data, shots_list, indices_no_disrupt, indices_disrupt,
                                  indices_flattop, indices_disrupt_time, name, exp)

if work_type == '2':
    name = '_CMod'
    exp = "VDE"
    all_data, shots_list = read_data_CMod()
    indices_no_disrupt, indices_disrupt, indices_flattop, indices_disrupt_time = define_indices_CMod(
        all_data)
    jj_50, ii = indices_operation(all_data, shots_list, indices_no_disrupt, indices_disrupt,
                                  indices_flattop, indices_disrupt_time, name, exp)


def norm_func(all_data, params, param, name, exp):
    Y = all_data[param]
    if name == '_CMod':
        if Y.name == 'BETAN':
            range_ = [0, 1.25]
        if Y.name == 'KAPPA':
            range_ = [0.8, 2.0]
        if Y.name == 'GAP_in':
            range_ = [0, 0.175]
        if Y.name == 'GAP_out':
            range_ = [0.025, 0.1]
        if Y.name == 'Q95':
            range_ = [2, 7.5]
        if Y.name == 'I_P':
            range_ = [0, 1.5e6]
        if Y.name == 'DELTA':
            range_ = [-0.2, 0.8]
        if Y.name == 'DELTA_TOP':
            range_ = [-0.2, 0.8]
        if Y.name == 'DELTA_BOTTOM':
            range_ = [-0.2, 0.8]
        if Y.name == 'LI':
            range_ = [0.75, 2.25]
    if name == '_TCV':
        if Y.name == 'BETAN':
            range_ = [0, 1.75]
        if Y.name == 'DELTA':
            range_ = [-2, all_data.DELTA.max()]
        if Y.name == 'DELTA_TOP':
            range_ = [-2.5, 5]
        if Y.name == 'DELTA_BOTTOM':
            range_ = [-2.5, all_data.DELTA_BOTTOM.max()]
        if Y.name == 'KAPPA':
            range_ = [all_data.KAPPA.min(), all_data.KAPPA.max()]
        if exp == 'VDE':
            if Y.name == 'GAP_in':
                range_ = [0, 0.09]
            if Y.name == 'GAP_out':
                range_ = [0, 0.15]
            if Y.name == 'Q95':
                range_ = [2, 8]
            if Y.name == 'I_P':
                range_ = [0, 0.4e6]
            if Y.name == 'LI':
                range_ = [0.5, 2.0]
        if exp == 'VD':
            if Y.name == 'GAP_in':
                range_ = [0, 0.15]
            if Y.name == 'GAP_out':
                range_ = [0, 0.1]
            if Y.name == 'Q95':
                range_ = [1.5, 7.5]
            if Y.name == 'I_P':
                range_ = [0, 4e5]
            if Y.name == 'LI':
                range_ = [0.5, 2.0]
    adn = all_data[params]
    #range_ = [all_data[param].min(), all_data[param].max()]
    adn1 = adn.loc[(adn[param] >= range_[0]) & ((adn[param] <= range_[1]))]
    cols = adn1.columns
    #scaler = RobustScaler()
    scaler = QuantileTransformer(output_distribution="uniform")
    normalized_all_data0 = scaler.fit_transform(adn1)
    normalized_all_data0 = pd.DataFrame(normalized_all_data0, columns=cols)
    return normalized_all_data0


def step_1(params_list, param_to_expl, normalized_all_data, ii, jj_50, D, flag):
    X = pd.concat([normalized_all_data.loc[normalized_all_data.index.isin(ii)],
                   normalized_all_data.loc[normalized_all_data.index.isin(jj_50)]]).sort_index()
    X = X.dropna()
    Y = X[param_to_expl]
    X = X.drop(param_to_expl, axis=1)
    D = D.iloc[D.index.isin(Y.index)]
    cm = CausalModel(Y.values, D.values, X.values)

    if flag == 1:
        print(cm.summary_stats)
        df_ndiff = pd.DataFrame(
            cm.summary_stats['ndiff'], index=X.columns, columns={'normalized_diff'})
        print(abs(df_ndiff).sort_values(by='normalized_diff'))
    return cm, D, X, Y


def propensity_est(cm):
    print(f'Each model propensity - {cm}')
    multiprocessing.Process(cm.est_propensity())
    print(cm.propensity)


def stratification(cm):
    print("Stratifying")
    cm.stratify_s()
    print(cm.strata)

    return cm.strata[len(cm.strata)-1]


def explore_last_strata(param, stratas, X_1, Y_1, D_1, flag):

    X_1[param] = Y_1

    s1 = stratas.summary_stats['Y_t_mean']-stratas.summary_stats['Y_t_sd']
    s2 = stratas.summary_stats['Y_t_mean']+stratas.summary_stats['Y_t_sd']

    X_12 = X_1.loc[(X_1[param] >= s1) & (X_1[param] <= s2)]

    ind_prop1 = X_12.index

    D_12 = D_1[ind_prop1]  # if index is in jj_50 - than 1, if not - 0

    X_12 = X_1[X_1.index.isin(ind_prop1)]
    Y_12 = Y_1[ind_prop1]

    X_12['treatment'] = D_12
    return X_12, Y_12, D_12


D = pd.concat([pd.Series(data=[0]*len(ii), index=ii),
               pd.Series(data=[1]*len(jj_50), index=jj_50)]).sort_index()  # if index is in jj_50 - than 1, if not - 0 - (now just for not disruptive)
params = {'DELTA_BOTTOM', 'DELTA', 'DELTA_TOP', 'GAP_in',
          'I_P', 'LI', 'KAPPA', 'GAP_out', 'BETAN',  'Q95'}
#params = {'BETAN', 'GAP_out', 'Q95', 'DELTA_BOTTOM'}
#param = ['BETAN', 'GAP_out', 'Q95', 'DELTA_BOTTOM']
param = ['BETAN', 'GAP_out', 'DELTA_BOTTOM', 'DELTA',
         'DELTA_TOP', 'GAP_in', 'I_P', 'LI', 'KAPPA', 'Q95']
all_data = all_data[params].dropna()
all_data.I_P.loc[all_data.I_P < 0] = abs(all_data.I_P.loc[all_data.I_P < 0])
all_data = all_data[all_data['BETAN'] < 2]


normalized_all_data = []
data_cm = {}
n = 0

for i in param:
    normalized_all_data.append(pd.DataFrame(norm_func(all_data, params, i, name, exp)))
    # print(normalized_all_data[n].median())
    print('\n', i, '\n')
    #cm[n], X_1[n], X_2[n], Y_1[n], Y_2[n], D_1[n], D_2[n] = step_1(params, param[n], normalized_all_data[n])

    print('Getting parameters')
    if "cm_1" in data_cm:
        data_cm["cm_1"].append(
            step_1(params, param[n], normalized_all_data[n], ii, jj_50, D, flag=0)[0])
    else:
        data_cm["cm_1"] = [
            step_1(params, param[n], normalized_all_data[n], ii, jj_50, D, flag=0)[0]]
    if "X_1" in data_cm:
        data_cm["X_1"].append(
            step_1(params, param[n], normalized_all_data[n], ii, jj_50, D, flag=0)[2])
    else:
        data_cm["X_1"] = [
            step_1(params, param[n], normalized_all_data[n], ii, jj_50, D, flag=0)[2]]
    if "Y_1" in data_cm:
        data_cm["Y_1"].append(
            step_1(params, param[n], normalized_all_data[n], ii, jj_50, D, flag=0)[3])
    else:
        data_cm["Y_1"] = [
            step_1(params, param[n], normalized_all_data[n], ii, jj_50, D, flag=0)[3]]

    if "D_1" in data_cm:
        data_cm["D_1"].append(
            step_1(params, param[n], normalized_all_data[n], ii, jj_50, D, flag=0)[1])
    else:
        data_cm["D_1"] = [
            step_1(params, param[n], normalized_all_data[n], ii, jj_50, D, flag=0)[1]]

    [step_1(params, param[n], normalized_all_data[n], ii, jj_50, D, flag=1)]
    n = n+1

print('Estimating propensity score')
for i in range(len(data_cm['cm_1'])):
    print(data_cm['Y_1'][i].name)
    propensity_est(data_cm['cm_1'][i])

for i in range(len(data_cm['cm_1'])):
    print(data_cm['Y_1'][i].name)
    for j in range(len(data_cm['X_1'][i].columns)):
        print('\t', data_cm['X_1'][i].columns[j], ' - ',
              data_cm['cm_1'][i].propensity['coef'][j+1])


for i in range(len(data_cm['cm_1'])):
    data_cm['X_1'][i]['propensity_score'] = data_cm['cm_1'][i].propensity['fitted']

for i in range(len(data_cm['cm_1'])):
    print(data_cm['Y_1'][i].name)
    if "stratas" in data_cm:
        data_cm["stratas"].append(stratification(data_cm['cm_1'][i]))
    else:
        data_cm["stratas"] = [stratification(data_cm['cm_1'][i])]

for i in range(len(data_cm['cm_1'])):
    print(data_cm['Y_1'][i].name)
    # exploring strata with with the biggest Propensity score
    print(data_cm['stratas'][i].summary_stats)

all_data = all_data[params]

data_cm_final = {}
n = 0
for i in param:
    flag = 0
    print(data_cm['Y_1'][n].name)
    if "X_12" in data_cm_final:
        data_cm_final["X_12"].append(explore_last_strata(
            i, data_cm['stratas'][n], data_cm['X_1'][n],
            data_cm['Y_1'][n], data_cm['D_1'][n], flag)[0])
    else:
        data_cm_final["X_12"] = [explore_last_strata(
            i, data_cm['stratas'][n], data_cm['X_1'][n],
            data_cm['Y_1'][n], data_cm['D_1'][n], flag)[0]]

    if "Y_12" in data_cm_final:
        data_cm_final["Y_12"].append(explore_last_strata(
            i, data_cm['stratas'][n], data_cm['X_1'][n],
            data_cm['Y_1'][n], data_cm['D_1'][n], flag)[1])
    else:
        data_cm_final["Y_12"] = [explore_last_strata(
            i, data_cm['stratas'][n], data_cm['X_1'][n],
            data_cm['Y_1'][n], data_cm['D_1'][n], flag)[1]]

    if "D_12" in data_cm_final:
        data_cm_final["D_12"].append(explore_last_strata(
            i, data_cm['stratas'][n], data_cm['X_1'][n],
            data_cm['Y_1'][n], data_cm['D_1'][n], flag)[2])
    else:
        data_cm_final["D_12"] = [explore_last_strata(
            i, data_cm['stratas'][n], data_cm['X_1'][n],
            data_cm['Y_1'][n], data_cm['D_1'][n], flag)[2]]

    flag = 1
    [explore_last_strata(
        i, data_cm['stratas'][n], data_cm['X_1'][n],
        data_cm['Y_1'][n], data_cm['D_1'][n], flag)]

    n = n+1

for i in range(len(data_cm['cm_1'])):
    if "X_12_orig" in data_cm_final:
        data_cm_final["X_12_orig"].append(
            all_data[params].iloc[data_cm_final['X_12'][i].index])
    else:
        data_cm_final["X_12_orig"] = [
            all_data[params].iloc[data_cm_final['X_12'][i].index]]
    '''data_cm_final['X_12_orig'][i]['treatment'] = data_cm_final['X_12'][i]['treatment']
  data_cm_final['X_12_orig'][i]['propensity_score'] = data_cm['X_1'][i]['propensity_score'].iloc[
      data_cm['X_1'][i].index.isin(data_cm_final['X_12'][i].index)]'''
    #data_cm_final["X_12_orig"][i] = data_cm_final["X_12_orig"][i].dropna()

for i in range(len(data_cm['cm_1'])):
    if "Y_12_orig" in data_cm_final:
        data_cm_final["Y_12_orig"].append(
            all_data[param[i]].iloc[data_cm_final['X_12'][i].index])
    else:
        data_cm_final["Y_12_orig"] = [
            all_data[param[i]].iloc[data_cm_final['X_12'][i].index]]
    if "D_12_orig" in data_cm_final:
        data_cm_final["D_12_orig"].append(D[data_cm_final['X_12'][i].index])
    else:
        data_cm_final["D_12_orig"] = [D[data_cm_final['X_12'][i].index]]

for i in range(len(param)):
    print(data_cm['Y_1'][i].name)
    if "cm_orig" in data_cm_final:
        data_cm_final["cm_orig"].append(CausalModel(data_cm_final['Y_12_orig'][i].values, data_cm_final['D_12_orig'][i].values,
                                                    data_cm_final['X_12_orig'][i][params].drop({param[i]}, axis=1).values))
    else:
        data_cm_final["cm_orig"] = [CausalModel(data_cm_final['Y_12_orig'][i].values, data_cm_final['D_12_orig'][i].values,
                                                data_cm_final['X_12_orig'][i][params].drop({param[i]}, axis=1).values)]

    print(data_cm_final["cm_orig"][i].summary_stats)

ATE = []
ATC = []
ATT = []
for i in data_cm_final["cm_orig"]:
    #print('now we are printing summary stats and estimations for last strata')
    i.est_via_ols()
    print('ATE:', i.estimates['ols']['ate'], 'ATC:',
          i.estimates['ols']['atc'], 'ATT:',
          i.estimates['ols']['att'])
    ATE.append(i.estimates['ols']['ate'])
    ATC.append(i.estimates['ols']['atc'])
    ATT.append(i.estimates['ols']['att'])

if name == '_CMod':
    ATE_CMod = ATE
    ATT_CMod = ATT
    ATC_CMod = ATC
if name == '_TCV' and exp == 'VDE':
    ATE_TCV_VDE = ATE
    ATT_TCV_VDE = ATT
    ATC_TCV_VDE = ATC
if name == '_TCV' and exp == 'VD':
    ATE_TCV_VD = ATE
    ATT_TCV_VD = ATT
    ATC_TCV_VD = ATC


for i in range(len(param)):
    print(param[i])
    plt.title(f"{data_cm_final['Y_12_orig'][i].name} vs Propensity score")

    treated_ind = list(data_cm_final["D_12_orig"]
                       [i].loc[data_cm_final["D_12_orig"][i] == 1].index)
    control_ind = list(data_cm_final["D_12_orig"]
                       [i].loc[data_cm_final["D_12_orig"][i] == 0].index)

    plt.scatter(all_data[param[i]].iloc[control_ind],
                data_cm['X_1'][i]['propensity_score'][control_ind], marker='o', c="b", alpha=0.5, label="Control")
    plt.scatter(all_data[param[i]].iloc[treated_ind],
                data_cm['X_1'][i]['propensity_score'][treated_ind], marker='o', c="r", alpha=0.7, label="Treatment")
    plt.xlabel(data_cm_final['Y_12_orig'][i].name)
    plt.ylabel("Propensity score")
    plt.legend(loc='upper right')

    plt.show()

ind_name = []
for i in data_cm_final['Y_12_orig']:
    ind_name.append(i.name)

if name == '_CMod':
    df_CMod = pd.DataFrame(
        data={'ATE': ATE_CMod, 'ATT': ATT_CMod, 'ATC': ATC_CMod}, index=ind_name, )
    fig, axes = plt.subplots(figsize=(60, 6), nrows=1, ncols=10)
    # ax2=plt.subplot(2,2,2)
    for i in range(len(df_CMod.index)):
        df_CMod.iloc[i].plot.bar(
            ax=axes[i], title=f'Treatment effects CMod - {ind_name[i]}', color={"tab:orange", "tab:green", "tab:blue"})
    plt.show()

if name == '_TCV':
    if exp == 'VDE':
        df_TCV_VDE = pd.DataFrame(
            data={'ATE': ATE_TCV_VDE, 'ATT': ATT_TCV_VDE, 'ATC': ATC_TCV_VDE}, index=ind_name, )
        fig, axes = plt.subplots(figsize=(60, 6), nrows=1, ncols=10)
        # ax2=plt.subplot(2,2,2)
        for i in range(len(df_TCV_VDE.index)):
            df_TCV_VDE.iloc[i].plot.bar(
                ax=axes[i], title=f'Treatment effects CMod - {ind_name[i]}', color={"tab:orange", "tab:green", "tab:blue"})
        plt.show()
    if exp == 'VD':
        df_TCV_VD = pd.DataFrame(
            data={'ATE': ATE_TCV_VD, 'ATT': ATT_TCV_VD, 'ATC': ATC_TCV_VD}, index=ind_name, )
        fig, axes = plt.subplots(figsize=(60, 6), nrows=1, ncols=10)
        # ax2=plt.subplot(2,2,2)
        for i in range(len(df_TCV_VD.index)):
            df_TCV_VD.iloc[i].plot.bar(
                ax=axes[i], title=f'Treatment effects CMod - {ind_name[i]}', color={"tab:orange", "tab:green", "tab:blue"})
        plt.show()
