'''
Created on Aug 20, 2020

@author: mina
'''
#Import modules
import numpy as np
from scipy import stats
import matplotlib as mpl
from data_set_example import seq_dic

#temp_list is the list of all temperatures at which we have csat data
temp_list = [4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 25, 28, 30, 32, 34, 35, 36, 38, 40]

#fmt_list is a list of marker types for plotting
fmt_list = ['o', 's', 'D', '^', 'v', '>', '<', 'p']

#Define constants
R = 1.987 * 10 ** -3
REF_CONC = 0.03

FONT = 'Arial'
LINE_WIDTH = 1
TITLE_SIZE = 7
LEGEND_SIZE = 7
AXES_LABEL_SIZE = 7
AXES_TICK_SIZE = 7
POINT_LABEL_SIZE = 7
ANNOTATION_SIZE = 7
MARKER_SIZE = 4
MARKER_EDGE_WIDTH = 0.4
TICK_LENGTH = 2
TICK_WIDTH = 1
TICK_PAD = 1
X_AXIS_PAD = 1
Y_AXIS_PAD = 1

mpl.rcParams['axes.linewidth'] = LINE_WIDTH

class RT_vs_Csat:
    """ RT_vs_Csat

    Calculates thermodynamic parameters based on the dilute arm of a phase diagram

    Arguments:
    ----------
        var_list: list
            The list of variants to analyze
    """
    def __init__(self, var_list):
        #Define the dictionaries and lists that we will need
        self.rt_dic = {}
        self.csat_dic = {}
        self.temp_dic = {}
        self.fit_dic = {}
        self.var_list = var_list
        #new_var_list keep track of which variants have at least 5 data points
        self.new_var_list = []
        
        #thermo_dic stores all of the final values
        self.thermo_dic = {}
        self.thermo_dic['enthalpy'] = {}
        self.thermo_dic['entropy'] = {}

        #Iterate through the given variants
        for var in self.var_list:
            self.rt_dic[var] = []
            self.csat_dic[var] = []
            self.temp_dic[var] = []
            #Iterate through all possible temperatures
            for temp in temp_list:
                csat = seq_dic[var]['csat'][str(temp)]
                #If a csat value exists at a given temperature, it is used in the analysis 
                if csat != 0:
                    self.rt_dic[var].append(1 / (R * (temp + 273.15)))
                    self.csat_dic[var].append(np.log(csat/REF_CONC))
                    self.temp_dic[var].append((temp + 273.15))
            if len(self.csat_dic[var]) >= 5:
                self.new_var_list.append(var)
            #fit_dic is a linear regression of ln(csat) vs. 1/RT
            #slope, intercept, r_value, p_value, std_err
            self.fit_dic[var] = stats.linregress(self.rt_dic[var],
                self.csat_dic[var])
            
            #The slope and y_intercept from fit_dic correlate with the enthalpy and entropy, respecitively
            self.thermo_dic['enthalpy'][var] = self.fit_dic[var][0]
            self.thermo_dic['entropy'][var] = self.fit_dic[var][1]
            
    def tick_param_setter(self, ax, x_direction = 'in', y_direction = 'in'):
        """ tick_param_setter
    
        A helper function to improve plot appearances
    
        Arguments:
        ----------
            ax: matplotlib axes object
                The axes object on which to plot
            direction: either "in" or "out"
                Determines whether the ticks are in or out of the panel frame
        """
        ax.tick_params(axis = 'x', labelsize = AXES_TICK_SIZE, direction = x_direction,
            length = TICK_LENGTH, width = TICK_WIDTH, color = 'black',
            labelcolor = 'black', pad = TICK_PAD)
        ax.tick_params(axis = 'y', labelsize = AXES_TICK_SIZE, direction = y_direction,
            length = TICK_LENGTH, width = TICK_WIDTH, color = 'black',
            labelcolor = 'black', pad = TICK_PAD)
        ax.xaxis.labelpad = X_AXIS_PAD
        ax.yaxis.labelpad = Y_AXIS_PAD
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(AXES_TICK_SIZE)
            tick.label1.set_fontweight('bold')
            tick.label1.set_fontname(FONT)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(AXES_TICK_SIZE)
            tick.label1.set_fontweight('bold')
            tick.label1.set_fontname(FONT)
    
    def ax_filler(self, ax, new_check = True, xlim = False, ylim=False, loc='lower left'):
        """ ax_filler
    
        A helper function to plot ln(csat) vs. 1/RT
    
        Arguments:
        ----------
            ax: matplotlib axes object
                The axes object on which to plot
            new_check: boolean, optional
                If true, will only plot variants with at least 5 data points
            xlim: list, optional
                The range of the x-axis
                If omitted, will use default range
            ylim: list, optional
                The range of the y-axis
                If omitted, will use default range
            loc: see matplotlib documentation, optional
                The location for the legend
                If omitted, will use lower left of panel
        """
        xlabel = '(RT)$^{-1}$ (mol / kcal)'
        ylabel = 'ln(c$_{sat}$)'
        if xlim:
            xlim=xlim
        else:
            xlim = [1.6, 1.83]
        if ylim:
            ylim=ylim
        else:
            ylim = [-9.5, -3]
        if new_check:
            var_list = self.new_var_list
        else:
            var_list = self.var_list
        rt_dic = self.rt_dic
        csat_dic = self.csat_dic
        fit_dic = self.fit_dic
        ax.set_xlabel(xlabel, fontname = FONT, size = AXES_LABEL_SIZE, weight = 'bold')
        ax.set_ylabel(ylabel, fontname = FONT, size = AXES_LABEL_SIZE, weight = 'bold')
        self.tick_param_setter(ax)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        for i, var in enumerate(var_list):
            ax.errorbar(x = rt_dic[var], y = csat_dic[var], fmt = fmt_list[i],
                color=seq_dic[var]['color'], capsize=2, markersize=MARKER_SIZE,
                label=seq_dic[var]['name'], markeredgewidth=0,
                markeredgecolor='black')
            ax.errorbar(x = rt_dic[var],
                y = [fit_dic[var][0] * x + fit_dic[var][1] for x in rt_dic[var]],
                linewidth = LINE_WIDTH, color=seq_dic[var]['color'])
        
        ax.errorbar(x = rt_dic[var_list[0]][0], y = csat_dic[var_list[0]][0],
            linewidth = LINE_WIDTH, color='grey', label='Linear fit')
        legend = ax.legend(loc = loc, edgecolor = 'black', prop={'size': LEGEND_SIZE, 'weight': 'bold', 'family': FONT})
        legend.get_frame().set_linewidth(LINE_WIDTH)
    