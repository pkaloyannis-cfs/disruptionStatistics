from TCV_CMod_1D2D_new_data import *

pd.options.mode.chained_assignment = None  # default='warn'


def plot_interp(df_interpolated):
    if name == '_CMod':
        plt.plot(df_interpolated.Time.loc[df_interpolated.Shot == 1120105029],
                 df_interpolated.BETAN.loc[df_interpolated.Shot == 1120105029])
        plt.title('interpolation betan', fontweight='bold')
        plt.xlabel('Time', fontsize='large', fontweight='bold')
        plt.ylabel('betan', fontsize='large', fontweight='bold')
        plt.show()


def plot_negative(all_data):
    if name == '_CMod':
        plt.plot(all_data.Time.loc[all_data.Shot == 1120105029],
                 all_data.BETAN.loc[all_data.Shot == 1120105029])
        plt.title('old betan data', fontweight='bold')
        plt.xlabel('Time', fontsize='large', fontweight='bold')
        plt.ylabel('betan', fontsize='large', fontweight='bold')
        plt.show()


def interp_data(all_data, name, exp, param_name):
    shot_u = all_data.Shot.loc[all_data[param_name] < 0].unique()
    if len(shot_u) == 0:
        print('No negative data to interpolate')
        return all_data
    else:
        print(f'Interpolation of {param_name} in progress')
        shots_to_remove = []
        for i in shot_u:
            sep_shot = all_data.loc[all_data.Shot == i]
            #print(y.BETAN.loc[y.BETAN<0].count(), y.BETAN.count())
            sep_shot_count = sep_shot.shape[0]
            neg_shot_count = sep_shot[param_name].loc[sep_shot[param_name] < 0].shape[0]
            l = list(sep_shot[param_name])
            if neg_shot_count > sep_shot_count/2 or ((l[-1] < 0) & (l[-2] < 0) & (l[-3] < 0)):
                shots_to_remove.append(i)

        all_data_new = all_data.loc[~all_data.Shot.isin(shots_to_remove)]
        shot_u = all_data_new.Shot.loc[all_data_new[param_name] < 0].unique()
        all_data_new[param_name][all_data_new[param_name] < 0] = np.nan

        data_interp = all_data_new

        for i in shot_u:
            sep_shot_param = all_data_new[param_name].loc[
                all_data_new.Shot == i].interpolate(method='from_derivatives')
            data_interp[param_name][all_data_new.Shot == i] = sep_shot_param
        #data_interp = all_data_new.interpolate(method='pad', limit_direction='forward')
        return data_interp


if name == '_CMod':
    data_interp = interp_data(all_data, name, exp, 'BETAN')
    print(data_interp.shape[0] / all_data.shape[0])
    data_interp = interp_data(data_interp, name, exp, 'Q95')
    print(data_interp.shape[0] / all_data.shape[0])
    data_interp = interp_data(data_interp, name, exp, 'GAP_in')
    data_interp.to_csv(f"new_data{name}_{exp}.csv", index=False)
    print(data_interp.shape[0] / all_data.shape[0])

    plot_interp(data_interp)
    plot_negative(all_data)


if name == '_TCV':
    param = ['BETAN', 'Q95']
    all_data[param[0]].loc[all_data[param[0]] < 0] = abs(
        all_data[param[0]].loc[all_data[param[0]] < 0])
    all_data[param[1]].loc[all_data[param[1]] < 0] = abs(
        all_data[param[1]].loc[all_data[param[1]] < 0])
    all_data.loc[all_data.filename ==
                 'Lite_Table_TCV_Disruptions_0.csv'].to_csv(
        'Lite_Table_TCV_Disruptions_0_new.csv', index=False)
    all_data.loc[all_data.filename ==
                 'Lite_Table_TCV_NonDisruptions.csv'].to_csv(
        'Lite_Table_TCV_NonDisruptions_new.csv', index=False)
    if exp == 'VD':
        all_data.loc[all_data.filename ==
                     'Lite_Table_TCV_Disruptions_VD.csv'].to_csv(
            'Lite_Table_TCV_Disruptions_VD_new.csv', index=False)
    if exp == 'VDE':
        all_data.loc[all_data.filename ==
                     'Lite_Table_TCV_Disruptions_VDE.csv'].to_csv(
            'Lite_Table_TCV_Disruptions_VDE_new.csv', index=False)
