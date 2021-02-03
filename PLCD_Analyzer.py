'''
Created on Feb 2, 2021

@author: mina
'''
#Import modules
import numpy as np
from scipy import optimize
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import Rescaled_Csat
import vant_Hoff_Csat
from data_set_example import seq_dic

#Define constants
FONT = 'Arial'
LINE_WIDTH = 1
TITLE_SIZE = 6
LEGEND_SIZE = 6
AXES_LABEL_SIZE = 6
AXES_TICK_SIZE = 6
POINT_LABEL_SIZE = 6
ANNOTATION_SIZE = 6
MARKER_SIZE = 4
MARKER_EDGE_WIDTH = 0.4
TICK_LENGTH = 2
TICK_WIDTH = 1
TICK_PAD = 1
X_AXIS_PAD = 1
Y_AXIS_PAD = 1

#Set constants for matplotlib
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 1

#cm2inch is a helper function to convert centimeter inputs to inches for figure creation
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

#The list of variants to use in the mean-field parameterization
var_list_charged_short = ['wt', '7r', '-10r', '12d',
            '8d', '-4d', '7k12d', 'rk1', 'rk2', '7r12d']

#The list of variants to use in the mean-field final regression and plotting
var_list_charged_full = ['wt', 'k2g', '2r', '7r', '-10r', '12d',
            '8d', '4d', '-4d', '7k12d', 'rk1', 'rk2', '7r12d',
            '-2k-4r5d']

#The list of variants to use in the thermodynamic analysis
var_list_aro = ['wtAro', 'allY', 'allF', 'aro1', 'f15y', 'f2y', 'f3y']

#Two temperatures lists which can be used in the mean-field model
temp_list_short = ['4']
temp_list_full = ['4', '8', '12', '16', '20', '24']

#The penalty for having a negative rescaled csat value in the mean-field parameterization
penalty = 16

#residual_getter calculates the sum of the residuals from the regressions in the mean-field parameterizations
def residual_getter(x0, var_list, temp_list, penalty, var_list_full = 0):
    residual_class = Rescaled_Csat.rescaler(x0, var_list, temp_list, penalty, var_list_full)
    return residual_class.residual_returner()

#fit is the final parameterization of the mean-field model
fit = optimize.minimize(residual_getter, [2, 0.05], (var_list_charged_short, temp_list_full, penalty), 
    method = 'Nelder-Mead', options={'maxfev': 200, 'xatol': 0.0001, 'fatol': 0.0001})

#Initialize the figure
fig = plt.figure(figsize=cm2inch(18, 9))
gs = fig.add_gridspec(1, 2)
ax_matrix = [None] * 2
ax_matrix[0] = fig.add_subplot(gs[0, 0])
ax_matrix[1] = fig.add_subplot(gs[0, 1])
  
#Plot the mean-field V-plot onto the figure
v_plot_class = Rescaled_Csat.rescaler(fit.x, var_list_charged_short, temp_list_full, penalty, var_list_full = var_list_charged_full)
v_plot_class.ax_filler(ax_matrix[0], 'NCPR', 'Rescaled $c_{sat}$', v_plot_class.csat_list_aro_r_k_full, temp_list_full, 1)

#Create the color bar to show the temperature gradient
cax = ax_matrix[0].inset_axes([0.65, 0.95, 0.25, 0.03])
for axis in ['top','bottom','left','right']:
    cax.spines[axis].set_linewidth(10)
    
cmap = mpl.cm.binary_r
norm = mpl.colors.Normalize(vmin=4, vmax=24)
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
cb.set_label('Temperature ($^o$C)', fontname = FONT, size = LEGEND_SIZE, weight = 'bold')
cb.set_ticks(range(4, 25, 4))
cb.set_ticklabels(range(4, 25, 4))
cb.ax.tick_params(labelsize = LEGEND_SIZE, length = TICK_LENGTH, width = TICK_WIDTH,
    color = 'black', labelcolor = 'black', pad = TICK_PAD)
cb.ax.xaxis.labelpad = X_AXIS_PAD
cb.ax.yaxis.labelpad = Y_AXIS_PAD
for tick in cb.ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(AXES_TICK_SIZE)
    tick.label1.set_fontweight('bold')
    tick.label1.set_fontname(FONT)

#Create the legend for the rescaled csat data
legend_elements_full = [Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markeredgecolor='black', markersize = MARKER_SIZE, markeredgewidth=MARKER_EDGE_WIDTH, label='WT'),
                        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markeredgecolor='black', markersize = MARKER_SIZE, markeredgewidth=MARKER_EDGE_WIDTH, label='Asp variants'),
                        Line2D([0], [0], marker='o', color='w', markerfacecolor='royalblue', markeredgecolor='black', markersize = MARKER_SIZE, markeredgewidth=MARKER_EDGE_WIDTH, label='Arg variants'),
                        Line2D([0], [0], marker='o', color='w', markerfacecolor='limegreen', markeredgecolor='black', markersize = MARKER_SIZE, markeredgewidth=MARKER_EDGE_WIDTH, label='Lys variants'),
                        Line2D([0], [0], marker='D', color='w', markerfacecolor='black', markeredgecolor='black', markersize = MARKER_SIZE, markeredgewidth=MARKER_EDGE_WIDTH, label='Unused for parameterization')
                        ]

lgd = [None] * 1
for i in range(1):
    our_legend = legend_elements_full
    loc = [0.15, 0.77]
    lgd[i] = ax_matrix[i].legend(handles=our_legend, edgecolor = 'black',
        prop={'size': LEGEND_SIZE, 'family': FONT, 'weight': 'bold'},
        loc=loc)
    lgd[i].get_frame().set_linewidth(LINE_WIDTH)

ax_matrix[0].plot([min(v_plot_class.ncpr_list_left), seq_dic['7r12d']['ncpr']],
    [np.exp(v_plot_class.intercept_left + v_plot_class.slope_left * min(v_plot_class.ncpr_list_left)), 
     np.exp(v_plot_class.intercept_left + v_plot_class.slope_left * seq_dic['7r12d']['ncpr'])],
    '--', linewidth = LINE_WIDTH, color='red', label='Linear fit')
ax_matrix[0].plot([seq_dic['7r12d']['ncpr'], max(v_plot_class.ncpr_list_right)], 
    [np.exp(v_plot_class.intercept_right + v_plot_class.slope_right * seq_dic['7r12d']['ncpr']), 
     np.exp(v_plot_class.intercept_right + v_plot_class.slope_right * max(v_plot_class.ncpr_list_right))],
    '--', linewidth = LINE_WIDTH, color='red')
ax_matrix[0].annotate('$r_{left}$ = %.2f' % (v_plot_class.r_value_left), (0.05, 0.05), weight = 'bold',
    fontname = FONT, color = 'red', xycoords='axes fraction', size = ANNOTATION_SIZE)
ax_matrix[0].annotate('$r_{right}$ = %.2f' % (v_plot_class.r_value_right), (0.7, 0.05), weight = 'bold',
    fontname = FONT, color = 'red', xycoords='axes fraction', size = ANNOTATION_SIZE)

#Plot the thermodynamic analysis onto the figure
vant_hoff_aro = vant_Hoff_Csat.RT_vs_Csat(var_list_aro)
vant_hoff_aro.ax_filler(ax_matrix[1], new_check=False)

#Show the plot
plt.tight_layout()
plt.show()
