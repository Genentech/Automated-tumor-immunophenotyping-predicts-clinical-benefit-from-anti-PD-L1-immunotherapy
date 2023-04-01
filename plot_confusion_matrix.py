
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap

from sklearn.metrics import cohen_kappa_score

viridis = cm.get_cmap('viridis', 4)

newcolors2 = viridis(np.linspace(0, 1, 3))
orange = np.array([0.929, 0.490, 0.192, 1])
blue = np.array([0.267, 0.447, 0.769, 1])
gray = np.array([0.647, 0.647, 0.647, 1])
black = np.array([0, 0, 0, 1])
newcolors2[0, :] = orange
newcolors2[1, :] = blue
newcolors2[2, :] = gray
newcmp2 = ListedColormap(newcolors2)

colors = [orange, blue, gray]


def draw_stats(left_labels, right_labels, fn_out_png):
    print(fn_out_png)
    left_un = np.array(['Inflamed', 'Excluded', 'Desert'])
    right_un = np.array(['Inflamed', 'Excluded', 'Desert'])
    xx = []
    yy = []
    cc = []
    stats = np.zeros((len(left_un), len(right_un)))
    for i in range(len(left_un)):
        #print(left_un[i], left_labels[np.where(left_labels == left_un[i])].shape)
        xxx = []
        yyy = []
        ccc = []
        for j in range(len(right_un)):
            xxx.append(j)
            yyy.append(i)
            ccc.append(np.array(colors[i][:3]))
            (n, ) = left_labels[np.where(
                np.logical_and(left_labels == left_un[i], right_labels == right_un[j]))].shape
            stats[i, j] = n
            #print(left_un[i], right_un[j], n)
        xx.extend(xxx)
        yy.extend(yyy)
        cc.extend(ccc)
    xx = np.array(xx)
    yy = np.array(yy)
    cc = np.array(cc)
    stats = stats.flatten()
    plt.figure()
    plt.scatter(xx, yy, s=7*stats, color=cc)
    for i in range(9):
        plt.annotate(int(stats[i]), (xx[i] + 0.2, yy[i]), fontsize=16)
    plt.xlim([-0.5, 2.5])
    plt.ylim([-0.5, 2.5])
    plt.xticks([0, 1, 2], left_un)
    plt.yticks([0, 1, 2], left_un)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    plt.savefig(fn_out_png)
    kappa = cohen_kappa_score(left_labels, right_labels, weights='linear')
    print('kappa: ', kappa)





def draw_distribution(df, title):
    OAK_A_distribution = pd.DataFrame({'Density cut-off': df.cutoff.value_counts(),
                                       "Binned CD8 density": df.SVM.value_counts(),
                                       'MOCHA': df.spat.value_counts(),
                                       'MOCHA-BITE': df.MOCHA.value_counts(),
                                       'Manual': df.Manual.value_counts()
                                       })
    OAK_A_distribution = OAK_A_distribution.reindex(['Inflamed', 'Excluded', 'Desert'])

    OAK_A_distribution = OAK_A_distribution.T
    
    plt.subplots_adjust(left = 0.2)
    OAK_A_distribution.plot(kind='bar', stacked=True, color = colors, rot=270, fontsize = 20).legend(loc='center left',bbox_to_anchor=(1, 0.5))
    plt.savefig('/Users/lix233/gRED_DP_OAK/' + title + '_all_pred_distribution.png', dpi = 300, bbox_inches = "tight")
    
    
    
    
fn_in = '/Users/lix233/gRED_DP_OAK/pip1_4_pred_OAK.csv'

    
data = pd.read_csv(fn_in)    
    
    
manual_IP = data.Manual.to_numpy()
cutoff_pred = data.cutoff.to_numpy()
SVM_pred = data.SVM.to_numpy()
spat_pred = data.spat.to_numpy()
MOCHA_pred = data.MOCHA.to_numpy()



draw_stats(manual_IP, cutoff_pred, '/Users/lix233/gRED_DP_OAK/OAK_manual_IP_vs_cutoff_predict.png')
draw_stats(manual_IP, SVM_pred, '/Users/lix233/gRED_DP_OAK/OAK_manual_IP_vs_SVM_predict.png')
draw_stats(manual_IP, spat_pred, '/Users/lix233/gRED_DP_OAK/OAK_manual_IP_vs_spat.png')
draw_stats(manual_IP, MOCHA_pred, '/Users/lix233/gRED_DP_OAK/OAK_manual_IP_vs_MOCHA.png')


A_data = data.loc[data.ACTARM == 'atezo']
D_data = data.loc[data.ACTARM == 'doce']


draw_distribution(A_data, 'OAK_A')
draw_distribution(D_data, 'OAK_D')






fn_in = '/Users/lix233/gRED_DP_OAK/pip1_4_pred_IMP130.csv'

    
data = pd.read_csv(fn_in)    
    
    
manual_IP = data.Manual.to_numpy()
cutoff_pred = data.cutoff.to_numpy()
SVM_pred = data.SVM.to_numpy()
spat_pred = data.spat.to_numpy()
MOCHA_pred = data.MOCHA.to_numpy()



draw_stats(manual_IP, cutoff_pred.astype(str), '/Users/lix233/gRED_DP_OAK/IMP130_manual_IP_vs_cutoff_predict.png')
draw_stats(manual_IP, SVM_pred.astype(str), '/Users/lix233/gRED_DP_OAK/IMP130_manual_IP_vs_SVM_predict.png')
draw_stats(manual_IP, spat_pred, '/Users/lix233/gRED_DP_OAK/IMP130_manual_IP_vs_spat.png')
draw_stats(manual_IP, MOCHA_pred, '/Users/lix233/gRED_DP_OAK/IMP130_manual_IP_vs_MOCHA.png')


A_data = data.loc[data.ACTARM == 'MPDL3280A']
D_data = data.loc[data.ACTARM == 'PLACEBO']


draw_distribution(A_data, 'IMP130_A')
draw_distribution(D_data, 'IMP130_D')




