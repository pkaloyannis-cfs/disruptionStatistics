# -*- coding: utf-8 -*-
"""
@author: yalynskayn
Description of the variables used:
indices_no_disrupt - indices of no disrupted shots
indices_disrupt - indices of disrupted shots
indices_flattop - indices of all flattop data
indices_disrupt_time – indices of disruptions time without dublicates
dipprog_dt - time derivative of Ip(programmed)
indices_disrupt_time_in_flattop - indices of the rows when exactly the disruption happened in the flattop data, same as indices_disrupt_time
indices_disrupt_in_flattop - indices of the whole shot numbers which end in disruption
indices_flattop_disrupt_in_flattop - indices of the whole shot numbers in the database which end in disruption in the flattop data
indices_flattop_no_disrupt - indices of non disruptive shots in flattop
ii (same as indices of non disruptive shots in flattop) and ff (indices of disruptive in the flattop shots) are to simplify the names 
idx - indices in shot variable that coincide with shot numbers in shotlist
jj - subselected flattop phase of disruptive data, part of VDE shotlist
jj_before_40 - data for VDEs flattop where time_until_disrupt<=40ms
jj_50 - isolate for VDE shots behavior around 50ms before disruption
other_disruptions_flattop - shotlist of all other flattop disruptions that are not VDEs
jj2 – flattop indices of all flattop disruptions that are not VDEs, same as other_disruptions_flattop
denom_idx - concatenate indices in the denominator
______________________________________________________________________________

This code “TCV_CMod_1D2D.py” can be used to plot 1D and 2D VDE Disruptivity maps for TCV and CMod tokamaks.
To run the file following documents are needed:
1.	C-Mod: CMod_db_cross.csv, VDE_shotlist.txt.
2.	TCV: 'Lite_Table_TCV_Disruptions_VD.csv', 'Lite_Table_TCV_Disruptions_0.csv', 'Lite_Table_TCV_NonDisruptions.csv'. 'TCV_parameters.xlsx'.
These files can be found on MECA : https://nucleus.iaea.org/sites/fusionportal/MECA/SitePages/Home.aspx
______________________________________________________________________________
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# ignore the NaNs in the future operations for disrupt per sec calculations
np.seterr(divide='ignore', invalid='ignore')


def fasten_TCV():
    file_name = ['Table_TCV_Disruptions_RE.csv', 'Table_TCV_Disruptions_below50kA.csv', 'Table_TCV_Disruptions_VD.csv',
                 'Table_TCV_Disruptions_RD.csv', 'Table_TCV_Disruptions_VDE.csv', 'Table_TCV_Disruptions_0.csv',
                 'Table_TCV_NonDisruptions.csv']
    for file in file_name:
        f = pd.read_csv(file)
        f.drop(['NEavg', 'Iref', 'TEaxis', 'SSXcore', 'Wtot', 'NElid', 'POHM_conf', 'PradTot', 'H98y2calc',
                'TAU_conf_calc', 'NBI', 'ECRH', 'powers_calc', 'OddN', 'EvenN', 'BZERO', 'Vloop', 'BETAP', 'BETAT',
                'QAX', 'TBEO', 'RMAG', 'ZMAG', 'ZCC', 'AREA', 'a_minor', 'R_geom', 'VOL', 'P_LH'], axis=1, inplace=True)
        f.to_csv('Lite_' + file)

# TCV 2022 data


def read_data_TCV(exp):
    # concatenate ALL the Table file, might be useful for the future
    # all_files=sorted(glob('Lite_Table_TCV_*.csv'))
    # defines VD or VDE shotlist from the master table according to the "Flag" parameter
    print('Reading TCV parameters file...')
    mtable = pd.read_excel('TCV_parameters.xlsx', sheet_name='TCV Disr > 50kA')
    if exp == 'VD':
        print('Reading files with VD, other disruptions and non disruptive data...')
        all_files = sorted(('Lite_Table_TCV_Disruptions_VD.csv', 'Lite_Table_TCV_Disruptions_0.csv',
                            'Lite_Table_TCV_NonDisruptions.csv'))
        tcv_all = pd.concat((pd.read_csv(file).assign(filename=file)
                             for file in all_files), ignore_index=True)
        shots_list = mtable[(mtable['Flag']) == 'VD'].shot.values
        np.savetxt("vdshots_tcv.txt", shots_list, fmt="%s")
    if exp == 'VDE':
        print('Reading files with VDE, other disruptions and non disruptive data...')
        all_files = sorted(('Lite_Table_TCV_Disruptions_VDE.csv', 'Lite_Table_TCV_Disruptions_0.csv',
                            'Lite_Table_TCV_NonDisruptions.csv'))
        tcv_all = pd.concat((pd.read_csv(file).assign(filename=file)
                             for file in all_files), ignore_index=True)
        shots_list = mtable[(mtable['Flag']) == 'VDE'].shot.values
        np.savetxt("vdeshots_tcv.txt", shots_list, fmt="%s")
    return tcv_all, shots_list


def read_data_CMod():
    print('Reading CMod data...')
    cmod = pd.read_csv('CMod_db_cross.csv')
    #cmod = pd.read_csv('new_data_CMod_VDE.csv')
    tritop = cmod.DELTA_TOP
    tribot = cmod.DELTA_BOTTOM
    aminor = cmod.aminor
    rmajor = cmod.sx
    aspect_ratio = rmajor / aminor
    overall_tri = (tribot + tritop) / 2
    # add columns
    # add aspect ratio to the 11th column in the dataset after the Rmajor
    cmod.insert(10, 'aspect_ratio', aspect_ratio)
    # add global triangularity to the 9th column in the dataset after the tribot
    cmod.insert(8, 'DELTA', overall_tri)
    # clear from memory anything except the part of interest
    cmod.drop(['n_equal_1_mode', 'n_equal_1_normalized', 'Te_peaking_ECE', 'ne_peaking', 'Mirnov', 'Mirnov_norm_btor',
               'Te_peaking', 'SXR_peaking', 'SXR', 'kappa_area'], axis=1, inplace=True)
    # read shotlist
    VDEshotlist = pd.read_csv('VDE_shotlist.txt')
    VDEshotlist.to_csv('VDE_shotlist.csv')
    vdeshots = pd.read_csv('VDE_shotlist.csv', index_col=None, header=None)
    return cmod, vdeshots


# for TCV data 2022
def define_indices_TCV(tcv_all):
    print('Definig indices...')
    # Define indices
    indices_no_disrupt = list(
        tcv_all.loc[tcv_all['filename'] == 'Lite_Table_TCV_NonDisruptions.csv'].index)
    indices_disrupt = list(
        tcv_all.loc[tcv_all['filename'] != 'Lite_Table_TCV_NonDisruptions.csv'].index)
    indices_flattop = list(tcv_all.index)
    # find indices of the rows when exactly the disruption happend, getting rid of all the dublicated times
    indices_disrupt_time = list(
        tcv_all.Shot[indices_disrupt].index[~tcv_all.Shot[indices_disrupt].duplicated()])
    return indices_no_disrupt, indices_disrupt, indices_flattop, indices_disrupt_time


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
            machine_param.loc[machine_param['filename'] == 'Lite_Table_TCV_Disruptions_0.csv'].index)
        # all VDE/VD data reading
        if exp == 'VD':
            all_exp_file = machine_param.loc[machine_param['filename']
                                             == 'Lite_Table_TCV_Disruptions_VD.csv']
        if exp == 'VDE':
            all_exp_file = machine_param.loc[machine_param['filename']
                                             == 'Lite_Table_TCV_Disruptions_VDE.csv']
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
    denom_idx = np.concatenate((ii, other_disruptions_flattop, jj_before_40), axis=0)
    return jj_50, denom_idx


# plotting 1d and 2d maps
def plot_disruptivity_1D(X, X_range, width, jj, denom_idx, deltat):
    # defining the bin edges
    bins = np.histogram_bin_edges(X, bins=25, range=X_range)
    # numerator
    nn = np.histogram(X[jj], bins)[0]
    # denminator
    dd = np.histogram(X[denom_idx], bins)[0]
    vps = (nn) / (dd * deltat)  # disruptivity
    err = np.sqrt(nn) / (dd * deltat)
    # plot error bars
    plt.figure(figsize=(8, 6), dpi=100, linewidth=0.0)
    plt.bar(bins[0:-1], vps, color='firebrick', width=width,
            align='center', zorder=2, edgecolor=None)
    plt.errorbar(bins[0:-1], vps, yerr=err, fmt='.', markersize=1, barsabove=True, ecolor='black', elinewidth=2,
                 capsize=3, capthick=2,
                 alpha=1, zorder=3)
    plt.xticks(size='large', fontweight='bold')
    plt.yticks(size='large', fontweight='bold')
    # plt.title('1D disruptivity of pure VDEs', fontweight='bold')
    # plt.title('1D disruptivity of pure VDEs',fontweight='bold')
    plt.yscale('log')
    plt.grid(visible=None, which='both', axis='both')


def plot_disruptiviy_2D(X, Y, X_range, Y_range, jj, denom_idx, deltat, X_lab, Y_lab):
    plt.figure(figsize=(8, 6), dpi=100)
    # plt.subplot(1,2,1)
    # numerator
    H1, xedges1, yedges1 = np.histogram2d(X[jj], Y[jj], bins=25, range=[X_range, Y_range])
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
    cax = plt.imshow(vps.T, cmap='viridis', origin='lower', interpolation='spline16', aspect='auto',
                     extent=[xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]])
    plt.xticks(size='large', fontweight='bold')
    plt.yticks(size='large', fontweight='bold')
    plt.title('2D Disruptivity map', fontsize='large', fontweight="bold")
    plt.xlabel(X_lab, fontsize='large', fontweight="bold")
    plt.ylabel(Y_lab, fontsize='large', fontweight="bold")
    plt.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    cbar = plt.colorbar(cax, label='Disruptivity ($s^{-1}$)')
    cbar.ax.tick_params(labelsize='large')
    cbar.set_label(label='Disruptivity ($s^{-1}$)', size='large')


# we are plotting the numerator and denminator counts of the disruptivity equation, as well as
# the 2D disruptivity map, in order to compare them with joint plots of parameters
def plotting_denom_and_nom(X, Y, X_range, Y_range, jj, denom_idx, deltat, X_lab, Y_lab):
    fig = plt.figure(figsize=(21, 6), dpi=100)
    # numerator
    plt.subplot(1, 3, 1)
    H1, xedges1, yedges1 = np.histogram2d(X[jj], Y[jj], bins=25, range=[X_range, Y_range])
    H2, xedges2, yedges2 = np.histogram2d(
        X[denom_idx], Y[denom_idx], bins=25, range=[X_range, Y_range])
    cax = plt.imshow(H1.T, cmap='viridis', origin='lower', interpolation='spline16', aspect='auto',
                     extent=[xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]])
    fig.colorbar(cax, label='counts')
    plt.yticks(size='large', fontweight='bold')
    plt.xticks(size='large', fontweight='bold')
    plt.xlabel(X_lab, fontsize='large', fontweight="bold")
    plt.ylabel(Y_lab, fontsize='large', fontweight="bold")
    plt.title('numerator count', fontweight='bold')
    plt.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    # denominator
    plt.subplot(1, 3, 2)
    cax1 = plt.imshow(H2.T, cmap='viridis', origin='lower', interpolation='spline16', aspect='auto',
                      extent=[xedges2[0], xedges2[-1], yedges2[0], yedges2[-1]])
    fig.colorbar(cax1, label='counts')
    plt.yticks(size='large', fontweight='bold')
    plt.xticks(size='large', fontweight='bold')
    plt.xlabel(X_lab, fontsize='large', fontweight="bold")
    plt.ylabel(Y_lab, fontsize='large', fontweight="bold")
    plt.title('denominator count', fontweight='bold')
    plt.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    # disruptivity map
    plt.subplot(1, 3, 3)
    total_sum = np.sum(H1) + np.sum(H2)
    weight = (H1 + H2) / total_sum
    vpsX_Y = (H1 / H2) * weight / deltat
    mask = np.where(~np.isfinite(vpsX_Y))
    vpsX_Y[mask] = 0.0
    cax2 = plt.imshow(vpsX_Y.T, cmap='viridis', origin='lower', interpolation='spline16', aspect='auto',
                      extent=[xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]])
    fig.colorbar(cax2, label='Disruptivity ($s^{-1}$)')
    plt.yticks(size='large', fontweight='bold')
    plt.xticks(size='large', fontweight='bold')
    plt.xlabel(X_lab, fontsize='large', fontweight="bold")
    plt.ylabel(Y_lab, fontsize='large', fontweight="bold")
    plt.title('disruptivity map', fontweight='bold')
    plt.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))


def plt_1D(machine_param, X_range, jj_50, denom_idx, deltat, name, exp, param_name):
    x = machine_param
    plot_disruptivity_1D(
        x, X_range, (X_range[1] - X_range[0]) / 25, jj_50, denom_idx, deltat)
    plt.ticklabel_format(axis='x', style="sci", scilimits=(0, 0))
    plt.title('disruptivity as a function of ' + param_name, fontweight='bold')
    plt.xlabel(param_name, fontsize='large', fontweight='bold')
    plt.ylabel('disruptivity ($s^{-1}$)', fontsize='large', fontweight='bold')
    im_name = f"{exp}_{param_name}_disruptivity_{name}.png"
    plt.savefig(im_name, dpi=100)
    plt.show()


def plt_2D(machine_param1, machine_param2, X_range, Y_range, jj_50, denom_idx, deltat, name, exp, param_name1, param_name2):
    x, y = machine_param1, machine_param2
    plot_disruptiviy_2D(x, y, X_range, Y_range, jj_50, denom_idx,
                        deltat, param_name1, param_name2)
    im_name = f"2D_{exp}_{param_name1}_{param_name2}_disruptivity_{name}.png"
    plt.savefig(im_name, dpi=100)
    plt.show()
    jp = sns.jointplot(x=x[jj_50], y=y[jj_50])
    jp.ax_marg_x.set_xlim(X_range)
    jp.ax_marg_y.set_ylim(Y_range)
    im_name = f"{exp}_{param_name1}_{param_name2}_joint_plot_{name}.png"
    plt.savefig(im_name, dpi=100)
    plt.show()
    #plotting_denom_and_nom(x, y, X_range, Y_range, jj_50, denom_idx, deltat, param_name1, param_name2)
    #im_name = f"2D_{exp}_{param_name1}_{param_name2}_disruptivity_denom_nom_{name}.png"
    #plt.savefig(im_name, dpi = 100)
    # plt.show()


work_type = input('Which installation do you want to choose?\
\n - Input 1 if TCV \
\n - Input 2 if CMod \ -->')
if work_type == '1':
    # function that drops some columns to fasten the data procession and saves the results in the new "Lite" versions
    # fasten_TCV()
    name = '_TCV'
    # takes shots with VD or VDE flag from TCV_parameters.xlsx
    exp = input('Do you want to study VD or VDE? (input VD or VDE) -->').upper()
    all_data, shots_list = read_data_TCV(exp)
    indices_no_disrupt, indices_disrupt, indices_flattop, indices_disrupt_time = define_indices_TCV(
        all_data)
    jj_50, denom_idx = indices_operation(all_data, shots_list, indices_no_disrupt, indices_disrupt,
                                         indices_flattop, indices_disrupt_time, name, exp)
    # sampling rate
    deltat = 0.002
    # parameters range for plotting
    DELTA_range = [all_data.DELTA.min(), all_data.DELTA.max()]
    DELTA_TOP_range = [all_data.DELTA_TOP.min(), all_data.DELTA_TOP.max()]
    DELTA_BOTTOM_range = [all_data.DELTA_BOTTOM.min(), all_data.DELTA_BOTTOM.max()]
    if exp == 'VDE':
        LI_range = [0.5, 3.0]
        GAP_in_range = [0, 0.09]
        GAP_out_range = [0, 0.15]
        Q95_range = [2, 8]
        I_P_range = [-0.6e6, 0.4e6]
    if exp == 'VD':
        LI_range = [0.5, 2.0]
        GAP_in_range = [0, 0.15]
        GAP_out_range = [0, 0.1]
        Q95_range = [0, 7.5]
        I_P_range = [-4.3e5, 4e5]

    BETAN_range = [0, 2.0]
    KAPPA_range = [all_data.KAPPA.min(), all_data.KAPPA.max()]

if work_type == '2':
    name = '_CMod'
    # exp can be only "VDE" because we have only VDE data for CMod
    exp = "VDE"
    all_data, shots_list = read_data_CMod()
    indices_no_disrupt, indices_disrupt, indices_flattop, indices_disrupt_time = define_indices_CMod(
        all_data)
    jj_50, denom_idx = indices_operation(all_data, shots_list, indices_no_disrupt, indices_disrupt,
                                         indices_flattop, indices_disrupt_time, name, exp)
    # sampling rate
    deltat = 0.005
    # parameters range for plotting
    DELTA_range = [-0.2, 0.8]
    DELTA_TOP_range = [-0.2, 0.8]
    DELTA_BOTTOM_range = [-0.2, 0.8]
    LI_range = [0.75, 2.25]
    BETAN_range = [0, 1.5]
    KAPPA_range = [0.8, 2.0]
    GAP_in_range = [0.05, 0.25]
    GAP_out_range = [0, 0.2]
    Q95_range = [2, 7.5]
    I_P_range = [-1.5e6, 1.5e6]


if __name__ == "__main__":

    flag = input('Do you want to plot 1D or 2D maps? (input 1D or 2D or both)-->').upper()
    if flag == '1D' or flag == 'BOTH':
        # plotting for 1D
        print('Plotting 1D Disruptivity ...')
        plt_1D(all_data.DELTA, DELTA_range, jj_50, denom_idx,
               deltat, name, exp, 'average_triangularity')
        plt_1D(all_data.DELTA_TOP, DELTA_TOP_range, jj_50,
               denom_idx, deltat, name, exp, 'triangularity_top')
        plt_1D(all_data.DELTA_BOTTOM, DELTA_BOTTOM_range, jj_50,
               denom_idx, deltat, name, exp, 'triangularity_bottom')
        plt_1D(all_data.LI, LI_range, jj_50, denom_idx, deltat, name, exp, 'li')
        plt_1D(all_data.BETAN, BETAN_range, jj_50, denom_idx, deltat, name, exp, 'betan')
        plt_1D(all_data.KAPPA, KAPPA_range,  jj_50,
               denom_idx, deltat, name, exp, 'elongation')
        plt_1D(all_data.GAP_in, GAP_in_range,  jj_50,
               denom_idx, deltat, name, exp, 'GAP_in')
        plt_1D(all_data.GAP_out, GAP_out_range,  jj_50,
               denom_idx, deltat, name, exp, 'GAP_out')
        plt_1D(all_data.Q95, Q95_range,  jj_50, denom_idx, deltat, name, exp, 'Q95')
        plt_1D(all_data.I_P, I_P_range,  jj_50, denom_idx, deltat, name, exp, 'I_P')

    if flag == '2D' or flag == 'BOTH':
        # plotting 2D
        print('Plotting 2D Disruptivity ...')
        plt_2D(all_data.LI, all_data.BETAN, LI_range, BETAN_range, jj_50, denom_idx, deltat, name, exp,
               'Li', 'betan')
        plt_2D(all_data.KAPPA, all_data.DELTA, KAPPA_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
               'elongation', 'average_triangularity')
        plt_2D(all_data.KAPPA, all_data.DELTA_TOP, KAPPA_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
               'elongation', 'triangularity_top')
        plt_2D(all_data.KAPPA, all_data.DELTA_BOTTOM, KAPPA_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
               'elongation', 'triangularity_bottom')
        plt_2D(all_data.KAPPA, all_data.GAP_in, KAPPA_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
               'elongation', 'GAP_in')
        plt_2D(all_data.KAPPA, all_data.GAP_out, KAPPA_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
               'elongation', 'GAP_out')
        plt_2D(all_data.KAPPA, all_data.Q95, KAPPA_range, Q95_range, jj_50, denom_idx, deltat, name, exp,
               'elongation', 'Q95')
        plt_2D(all_data.KAPPA, all_data.I_P, KAPPA_range, I_P_range, jj_50, denom_idx, deltat, name, exp,
               'elongation', 'I_P')

        plt_2D(all_data.BETAN, all_data.GAP_in, BETAN_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
               'betan', 'GAP_in')
        plt_2D(all_data.BETAN, all_data.GAP_out, BETAN_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
               'betan', 'GAP_out')
        plt_2D(all_data.BETAN, all_data.Q95, BETAN_range, Q95_range, jj_50, denom_idx, deltat, name, exp,
               'betan', 'Q95')
        plt_2D(all_data.BETAN, all_data.I_P, BETAN_range, I_P_range, jj_50, denom_idx, deltat, name, exp,
               'betan', 'I_P')

        plt_2D(all_data.LI, all_data.GAP_in, LI_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
               'Li', 'GAP_in')
        plt_2D(all_data.LI, all_data.GAP_out, LI_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
               'Li', 'GAP_out')
        plt_2D(all_data.LI, all_data.Q95, LI_range, Q95_range, jj_50, denom_idx, deltat, name, exp,
               'Li', 'Q95')
        plt_2D(all_data.LI, all_data.I_P, LI_range, I_P_range, jj_50, denom_idx, deltat, name, exp,
               'Li', 'I_P')

        plt_2D(all_data.GAP_out, all_data.GAP_in, GAP_out_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
               'GAP_out', 'GAP_in')
        plt_2D(all_data.I_P, all_data.Q95, I_P_range, Q95_range, jj_50, denom_idx, deltat, name, exp,
               'I_P', 'Q95')

        plt_2D(all_data.I_P, all_data.GAP_in, I_P_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
               'I_P', 'GAP_in')
        plt_2D(all_data.I_P, all_data.GAP_out, I_P_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
               'I_P', 'GAP_out')

        plt_2D(all_data.Q95, all_data.GAP_in, Q95_range, GAP_in_range, jj_50, denom_idx, deltat, name, exp,
               'Q95', 'GAP_in')
        plt_2D(all_data.Q95, all_data.GAP_out, Q95_range, GAP_out_range, jj_50, denom_idx, deltat, name, exp,
               'Q95', 'GAP_out')

        plt_2D(all_data.I_P, all_data.DELTA, I_P_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
               'I_P', 'average_triangularity')
        plt_2D(all_data.I_P, all_data.DELTA_TOP, I_P_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
               'I_P', 'triangularity_top')
        plt_2D(all_data.I_P, all_data.DELTA_BOTTOM, I_P_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
               'I_P', 'triangularity_bottom')

        plt_2D(all_data.Q95, all_data.DELTA, Q95_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
               'Q95', 'average_triangularity')
        plt_2D(all_data.Q95, all_data.DELTA_TOP, Q95_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
               'Q95', 'triangularity_top')
        plt_2D(all_data.Q95, all_data.DELTA_BOTTOM, Q95_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
               'Q95', 'triangularity_bottom')

        plt_2D(all_data.GAP_in, all_data.DELTA, GAP_in_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
               'GAP_in', 'average_triangularity')
        plt_2D(all_data.GAP_in, all_data.DELTA_TOP, GAP_in_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
               'GAP_in', 'triangularity_top')
        plt_2D(all_data.GAP_in, all_data.DELTA_BOTTOM, GAP_in_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
               'GAP_in', 'triangularity_bottom')

        plt_2D(all_data.GAP_out, all_data.DELTA, GAP_out_range, DELTA_range, jj_50, denom_idx, deltat, name, exp,
               'GAP_out', 'average_triangularity')
        plt_2D(all_data.GAP_out, all_data.DELTA_TOP, GAP_out_range, DELTA_TOP_range, jj_50, denom_idx, deltat, name, exp,
               'GAP_out', 'triangularity_top')
        plt_2D(all_data.GAP_out, all_data.DELTA_BOTTOM, GAP_out_range, DELTA_BOTTOM_range, jj_50, denom_idx, deltat, name, exp,
               'GAP_out', 'triangularity_bottom')
