'''
Created on Aug 13, 2020

@author: mina
'''
#Import modules
import numpy as np
from scipy import stats
import matplotlib.colors as colors
import matplotlib as mpl
from data_set_example import seq_dic

#Define constants
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

class rescaler:
    """ rescaler

    Calculates rescaled csat values based on mean-field model

    Arguments:
    ----------
        x0: list, shape = [aro_r_term, k_term]
            The initial lambda parameters for the mean-field model
        var_list: list
            The list of variants to parameterize the model on
        temp_list: list
            The list of temperatures at which there are csat values to analyze
        penalty: float
            The residual penalty associated with a negative rescaled csat value
            When this occurs, "penalty" will be returned as the residual value for the fitting process
            In general, this value should be relatively high so that the parameterization never causes negative csat values
        var_list_full: list, optional
            If included, will calculate final linear regression and plots using var_list_full
    """
    def __init__(self, x0, var_list, temp_list, penalty, var_list_full = 0):
        
        self.var_list = var_list
        self.temp_list = temp_list
        self.penalty = penalty
        self.var_list_full = var_list_full
        
        self.aro_aro = 1
        self.aro_r = x0[0]
        self.three_body_k = x0[1]
        
        #These dictionaries and lists allow us to build rescaled csat values with different models for plotting later
        self.ncpr_list = {}
        self.ncpr_list_left = []
        self.ncpr_list_right = []
    
        self.ncpr_list_full = {}
        self.ncpr_list_left_full = []
        self.ncpr_list_right_full = []
        
        #Raw csat lists
        self.csat_list = {}
        self.csat_list_left = []
        self.csat_list_right = []
        
        #csat lists rescaled by temperatures
        self.csat_list_temp = {}
        self.csat_list_left_temp = []
        self.csat_list_right_temp = []
    
        #csat lists rescaled using aro-aro and aro-arg interactions 
        self.csat_list_aro_r = {}
        self.csat_list_aro_r_left = []
        self.csat_list_aro_r_right = []
        
        #csat lists rescaled using aro-aro, aro-arg, and destabilizing lys interactions
        self.csat_list_aro_r_k = {}
        self.csat_list_aro_r_k_left = []
        self.csat_list_aro_r_k_right = []
        
        #rescaled csat lists using var_list_full, if applicable
        self.csat_list_aro_r_k_full = {}
        self.csat_list_aro_r_k_left_full = []
        self.csat_list_aro_r_k_right_full = []
    
        #m_dic will holds the slopes of the dilute arms
        self.m_dic = {}
        
        #res_penalty is zero unless the lambda parameters result in a negative rescaled csat value
        self.res_penalty = 0

        for temp in self.temp_list:
            
            self.ncpr_list[temp] = []
            self.ncpr_list_full[temp] = []
            self.csat_list[temp] = []
            self.csat_list_temp[temp] = []
            self.csat_list_aro_r[temp] = []
            self.csat_list_aro_r_k[temp] = []
            self.csat_list_aro_r_k_full[temp] = []
        
        if self.var_list_full == 0:
            self.var_list_used = self.var_list
        else:
            self.var_list_used = self.var_list_full
    
        #Calculate a master m-value by averaging the slopes of the dilute arms
        if len(self.temp_list) > 1:
            self.m_value_full = 0
            self.m_count_full = 0
            for var in self.var_list:
                temps = []
                csats = []
                for temp in self.temp_list:
                    csat = seq_dic[var]['csat'][temp]
                    if csat != 0:
                        temps.append(int(temp))
                        csats.append(np.log(csat))
                # slope, intercept, r_value, p_value, std_err
                fit = stats.linregress(csats, temps)
                if len(temps) > 1:
                    self.m_dic[var] = fit[0]
                    if var in self.var_list:
                        self.m_value_full += fit[0]
                        self.m_count_full += 1
            
            if self.m_value_full != 0:
                self.m_value_full /= self.m_count_full
            else:
                self.m_value_full = 1
        else:
            self.m_value_full = 1
        
        #The fit_list variables store the values for the final linear regressions of the V-shape
        self.fit_list_x_full_left = []
        self.fit_list_x_full_right = []
        self.fit_list_y_full_left = []
        self.fit_list_y_full_right = []
        
        #Iterate through the variants and fill out the previously initialized lists
        for var in self.var_list_used:
            self.fit_list_left = []
            self.fit_list_right = []        
            for temp in self.temp_list:
                if seq_dic[var]['csat'][temp] != 0:
    
                    temp_norm = np.exp((int(temp) - 4) / self.m_value_full)
                    #temp_norm = np.exp((int(temp) - 4) / m_dic[var])
                    
                    renorm_r_aro = 2 * self.aro_r * seq_dic[var]['aro'] * seq_dic[var]['r'] \
                        + self.aro_aro * seq_dic[var]['aro'] ** 2
                        
                    renorm_r_aro_k = (seq_dic[var]['aro'] ** 2) * (self.aro_aro - 3 * \
                        self.three_body_k * seq_dic[var]['k']) + (seq_dic[var]['aro'] * \
                        seq_dic[var]['r']) * (2 * self.aro_r - 6 * self.three_body_k * seq_dic[var]['k'])
                    
                    renorm_r_aro /= temp_norm
                    renorm_r_aro_k /= temp_norm
                    
                    if renorm_r_aro_k <= 0:
                        self.res_penalty = self.penalty
                    
                    if var in self.var_list_used:
                        self.ncpr_list[temp].append(seq_dic[var]['ncpr'])
                        self.csat_list[temp].append(seq_dic[var]['csat'][temp])
                        self.csat_list_temp[temp].append(seq_dic[var]['csat'][temp] / temp_norm)
                        self.csat_list_aro_r[temp].append(renorm_r_aro * seq_dic[var]['csat'][temp])
                        self.csat_list_aro_r_k[temp].append(renorm_r_aro_k * seq_dic[var]['csat'][temp])
                    self.ncpr_list_full[temp].append(seq_dic[var]['ncpr'])
                    self.csat_list_aro_r_k_full[temp].append(renorm_r_aro_k * seq_dic[var]['csat'][temp])                    
    
                    #For initial purposes, variant +7R+12D was used as the mid-point of the V-shape
                    #This can easily be changed by manipulating the following line
                    midpoint = seq_dic['7r12d']['ncpr']
                    
                    if seq_dic[var]['ncpr'] <= midpoint:
                        if var in self.var_list_used:
                            self.ncpr_list_left.append(seq_dic[var]['ncpr'])
                            self.csat_list_left.append(seq_dic[var]['csat'][temp])
                            self.csat_list_left_temp.append(seq_dic[var]['csat'][temp] / temp_norm)
                            self.csat_list_aro_r_left.append(renorm_r_aro * seq_dic[var]['csat'][temp])
                            self.csat_list_aro_r_k_left.append(renorm_r_aro_k * seq_dic[var]['csat'][temp])
                        self.ncpr_list_left_full.append(seq_dic[var]['ncpr'])
                        self.csat_list_aro_r_k_left_full.append(renorm_r_aro_k * seq_dic[var]['csat'][temp])
                        self.fit_list_left.append(renorm_r_aro_k * seq_dic[var]['csat'][temp])
            
                    if seq_dic[var]['ncpr'] >= midpoint:
                        if var in self.var_list_used:
                            self.ncpr_list_right.append(seq_dic[var]['ncpr'])
                            self.csat_list_right.append(seq_dic[var]['csat'][temp])
                            self.csat_list_right_temp.append(seq_dic[var]['csat'][temp] / temp_norm)
                            self.csat_list_aro_r_right.append(renorm_r_aro * seq_dic[var]['csat'][temp])
                            self.csat_list_aro_r_k_right.append(renorm_r_aro_k * seq_dic[var]['csat'][temp])
                        self.ncpr_list_right_full.append(seq_dic[var]['ncpr'])
                        self.csat_list_aro_r_k_right_full.append(renorm_r_aro_k * seq_dic[var]['csat'][temp])
                        self.fit_list_right.append(renorm_r_aro_k * seq_dic[var]['csat'][temp])
            
            if seq_dic[var]['ncpr'] <= midpoint:
                self.fit_list_left = np.mean(self.fit_list_left)
                self.fit_list_y_full_left.append(self.fit_list_left)
                self.fit_list_x_full_left.append(seq_dic[var]['ncpr'])
                
            if seq_dic[var]['ncpr'] >= midpoint:
                self.fit_list_right = np.mean(self.fit_list_right)
                self.fit_list_y_full_right.append(self.fit_list_right)
                self.fit_list_x_full_right.append(seq_dic[var]['ncpr'])
    
        #Perform a linear regression on either side of the V-shape
        self.slope_left, self.intercept_left, self.r_value_left, self.p_value_left, self.std_err_left = \
            stats.linregress(self.fit_list_x_full_left, np.log(self.fit_list_y_full_left))
        self.slope_right, self.intercept_right, self.r_value_right, self.p_value_right, self.std_err_right = \
            stats.linregress(self.fit_list_x_full_right, np.log(self.fit_list_y_full_right))
        
        #Calculate the total residual from the linear fits for the parameterization process 
        self.res = 0
        for i, ncpr in enumerate(self.fit_list_x_full_left):
            self.res += (np.log(self.fit_list_y_full_left[i]) - \
                (self.intercept_left + self.slope_left * ncpr)) ** 2
        for i, ncpr in enumerate(self.fit_list_x_full_right):
            self.res += (np.log(self.fit_list_y_full_right[i]) - \
                (self.intercept_right + self.slope_right * ncpr)) ** 2
    
    def tick_param_setter(self, ax, direction = 'in'):
        """ tick_param_setter
    
        A helper function to improve plot appearances
    
        Arguments:
        ----------
            ax: matplotlib axes object
                The axes object on which to plot
            direction: either "in" or "out"
                Determines whether the ticks are in or out of the panel frame
        """
        ax.tick_params(labelsize = AXES_TICK_SIZE, direction = direction,
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
                
    def ax_filler(self, ax, xlabel, ylabel, y_list, our_temp_list, full_check):
        """ ax_filler
    
        A helper function to plot rescaled csat values
    
        Arguments:
        ----------
            ax: matplotlib axes object
                The axes object on which to plot
            xlabel: string
                The label of the x-axis
            ylabel: string
                The label of the y-axis
            y_list: list
                The rescaled csat list to plot
                May be one of csat_list, csat_list_temp, csat_list_aro_r,
                csat_list_aro_r_k, or csat_list_aro_r_k_full
            our_temp_list: list
                The temperature list to use for plotting
            full_check: boolean
                Determines whether to use var_list or var_list_full
        """
        if full_check:
            var_list_used = self.var_list_full
            ncpr_list_used = self.ncpr_list_full
        else:
            var_list_used = self.var_list
            ncpr_list_used = self.ncpr_list
        
        ax.set_xlabel(xlabel, fontname = FONT, size = AXES_LABEL_SIZE, weight = 'bold')
        ax.set_ylabel(ylabel, fontname = FONT, size = AXES_LABEL_SIZE, weight = 'bold')
        ax.set(yscale="log")
        ax.minorticks_off()
        
        self.tick_param_setter(ax)
        for temp in our_temp_list:
            ml = ['o' if var in self.var_list else 'D' for var in var_list_used if seq_dic[var]['csat'][temp] != 0]
            #opacity function varies the opacity of a marker based on the given temperature
            def opacity(temp):
                if len(our_temp_list) > 1:
                    float_temp_list = [float(temp) for temp in our_temp_list]
                    return 1 + (min(float_temp_list) - float(temp)) / (max(float_temp_list) - min(float_temp_list))
                else:
                    return 1
            color_list = [colors.to_rgba(seq_dic[var]['colorFig2'], opacity(temp)) for var in var_list_used if seq_dic[var]['csat'][temp] != 0]
            for index, ncpr in enumerate(ncpr_list_used[temp]):
                ax.errorbar(x = ncpr, y = y_list[temp][index], fmt = ml[index],
                    markersize = MARKER_SIZE, markeredgecolor = 'black',
                    markeredgewidth = MARKER_EDGE_WIDTH, color = color_list[index],
                    label=temp + '$^\circ$C')
    
    #residual_returner retruns the residual calculated from above
    #If any of the rescaled csat values were negative, this function returns "res_penalty"
    def residual_returner(self):
        if self.res_penalty > 0:
            return self.res_penalty
        else:
            return self.res
