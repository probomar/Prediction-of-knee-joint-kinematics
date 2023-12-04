import matplotlib.pyplot as plt
from input import *
import pandas as pd
import numpy as np


folder_name_lin = 'Results4_linear'
folder_name_nonlin = 'Results4_non-linear'

ligaments = ['ACLa', 'ACLp', 'PCLa', 'PCLp', 'LCL', 'MCLa', 'MCLo', 'MCLd']
dflig_lin = pd.read_csv(folder_name_lin + '/ligaments/' + ligaments[0] + '.csv', sep=',', header=None)
lig_lin = dflig_lin.to_numpy()
dflig_nonlin = pd.read_csv(folder_name_nonlin + '/ligaments/' + ligaments[0] + '.csv', sep=',', header=None)
lig_nonlin = dflig_nonlin.to_numpy()

lig_lin_lenght = np.empty([lig_lin.shape[0], 0])
lig_nonlin_lenght = np.empty([lig_nonlin.shape[0], 0])

for i in range(len(ligaments)):
    dflig_lin = pd.read_csv(folder_name_lin + '/ligaments/' + ligaments[i] + '.csv', sep=',', header=None)
    lig_lin = dflig_lin.to_numpy()
    lig_lin_lenght = np.append(lig_lin_lenght, lig_lin[:, 0].reshape(lig_lin.shape[0], 1), axis=1)

    dflig_nonlin = pd.read_csv(folder_name_nonlin + '/ligaments/' + ligaments[i] + '.csv', sep=',', header=None)
    lig_nonlin = dflig_nonlin.to_numpy()
    lig_nonlin_lenght = np.append(lig_nonlin_lenght, lig_nonlin[:, 0].reshape(lig_nonlin.shape[0], 1), axis=1)

gait_cycle = motion[::step, 2]
gait_cycle_lin = gait_cycle[0:lig_lin_lenght.shape[0]]
gait_cycle_nonlin = gait_cycle[0:lig_nonlin_lenght.shape[0]]

fig, axs = plt.subplots(1, 2, constrained_layout=True, sharey=True)

for i in range(len(ligaments)):
    axs[0].plot(gait_cycle_lin, lig_lin_lenght[:, i], 'x', label=ligaments[i], markersize=3)
    axs[1].plot(gait_cycle_nonlin, lig_nonlin_lenght[:, i], 'x', label=ligaments[i], markersize=3)

axs[0].set_xlabel('% cyklu chůze')
axs[0].set_ylabel('Velikost síly [N]')
axs[0].legend()
axs[0].set_xlim(10, 50)
axs[0].set_ylim(0, 300)
axs[0].set_title('Lineární vazy')

axs[1].set_xlabel('% cyklu chůze')
axs[1].legend()
axs[1].set_xlim(10, 50)
axs[1].set_title('Nelineární vazy')

fig.savefig('results/ligaments.png', dpi=1000)
fig.show()
