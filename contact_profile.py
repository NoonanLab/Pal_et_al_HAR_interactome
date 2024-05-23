import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import find_peaks
import numpy as np

# Read the interaction data file into a DataFrame
df = pd.read_table("CHiC_file.inter.bb")  # any CHiC file in BigInteract format with HAR and HGE DA interactions
lst_original = list(df['sourceName'])  # Extract the 'sourceName' column into a list

# Get unique HAR/DA names
new_lst = []
i = 0
while i < len(lst_original):
    name = lst_original[i]
    # Adjust names to remove suffixes
    if name[-4:] == '_sec':
        name = name[0:-4]
    elif name[-5:] == '_tert':
        name = name[0:-5]
    new_lst.append(name)
    i += 1
    if i >= len(lst_original) - 1:
        break
    # Skip to the next unique name
    while i < len(lst_original) and lst_original[i].startswith(name):
        i += 1

# Dictionary to store KDE maxima for each HAR/DA
final_dict = {}
for i in new_lst:
    element_name = i
    # Identify rows corresponding to DpnII fragments for the current HAR/DA
    is_element_1 = df['sourceName'] == i
    is_element_2 = df['sourceName'] == i + '_sec'
    is_element_3 = df['sourceName'] == i + '_tert'
    is_element = list(map(any, zip(is_element_1, is_element_2, is_element_3, is_element_4)))
    # Subset DataFrame for the current HAR/DA
    df_element = df[is_element]
    # Compute the relative distance between HAR/DA fragment midpoint and the other end midpoint (in kb)
    df_element['rel_dist (kb)'] = '0'
    a1 = (df_element['targetStart'] + df_element['targetEnd']) / 2
    b1 = (df_element['sourceStart'] + df_element['sourceEnd']) / 2
    df_element['rel_dist (kb)'] = (a1 - b1) / 1000
    
    # Create a scatter plot of interaction strength
    ax1 = df_element.plot.scatter('rel_dist (kb)', 'value', c='value', colormap='cool', sharex=False, title=i + '-ortho')
    
    # Align axes of KDE plot with scatter plot
    ax2 = ax1.twinx()
    
    minimum = min(list(df_element['rel_dist (kb)']))  # lowest score of interaction strength
    maximum = max(list(df_element['rel_dist (kb)']))  # corresponding high score
    bandwidth = 30
    
    # Adjust bandwidth based on the interaction range
    if maximum - minimum > 800:
        bandwidth = 40
    elif maximum - minimum > 500 and maximum - minimum < 800:
        bandwidth = 35
    else:
        bandwidth = 30
        
    # Plot the KDE plot with variable bandwidth
    ax_d = sns.kdeplot(df_element['rel_dist (kb)'], ax=ax2, bw=bandwidth, label="Density distribution", color='r')
    
    # Format y-axis to keep only colorbar
    ax1.set_ylabel("")
    ax1.set_xlabel("Relative distance to HAR (kb)")
    ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y))
    ax2.yaxis.set_major_formatter(ticks)
    ax2.set_yticks([], [])
    ax1.set_yticks([], [])
    
    # Set colorbar label
    cax = plt.gcf().get_axes()[1]
    plt.legend(loc='best', ncol=1, prop={'size': 8.8})
    cax.set_ylabel('CHiCAGO Score')
    
    # Save the plot (commented out in original script)
    # plt.savefig('cNSC_profiles_P1_new/' + i + '-ortho_plot.eps', format='eps')
    
    # Extract KDE curve data
    kde_curve = ax_d.lines[0]
    x = kde_curve.get_xdata()
    y = kde_curve.get_ydata()
    # Find indices for maxima of interaction strength 
    indices = find_peaks(y)
    
    # Invert the KDE curve for finding minima
    inv_data_y = y * (-1)

    non_singleton_peaks = []  # x values of non-singleton peaks

    # Filter peaks based on criteria
    if len(indices[0]) == 1:
        for i in list(df_element["value"]):
            if i >= 4.5 and x[indices[0][0]] not in non_singleton_peaks:
                non_singleton_peaks.append(x[indices[0][0]])
                
    elif len(indices[0]) == 0:
        non_singleton_peaks = []
    
    else:
        valley_indexes = find_peaks(inv_data_y)
        
        valley_indexes1 = []  # list of x values of minima
        for i in valley_indexes[0]:
            valley_indexes1.append(x[i])

        # Ensure the outer x values are included
        if not (valley_indexes1[0] < x[indices[0][0]]):
            valley_indexes1.insert(0, x[0])
        if not (valley_indexes1[-1] > x[indices[0][-1]]):
            valley_indexes1.append(x[-1])

        # Map x and y values of the scatter plot to a dictionary
        x_y_mapping_dict = dict(zip(list(df_element["rel_dist (kb)"]), df_element["value"]))

        # Make non-singleton peak list for every HAR/DA
        for i in range(len(valley_indexes1) - 1):
            small = valley_indexes1[i]
            big = valley_indexes1[i + 1]
            count = 0
            a_value_more_than_five = False
            for j in list(df_element["rel_dist (kb)"]):
                if j >= small and j <= big:
                    count += 1
                    if x_y_mapping_dict.get(j, 1) >= 4.5:
                        a_value_more_than_five = True
            if count > 1 and a_value_more_than_five:
                non_singleton_peaks.append(x[indices[0][i]])
        
        # Function to get FWHM for each peak
        def get_rel_dist(minima, peaks, xlist, ylist):
            rel_dist = []  # relative distance to HAR range for each peak
            min_index = 0
            for peak in peaks:
                while min_index >= len(minima) - 1 and not (peak < minima[min_index] and peak > minima[min_index + 1]):
                    min_index += 1
                y_left_half_max = 0.5 * (np.interp(minima[min_index], xlist, ylist) + np.interp(peak, xlist, ylist))
                y_right_half_max = 0.5 * (np.interp(peak, xlist, ylist) + np.interp(minima[min_index + 1], xlist, ylist))
                range_x = np.interp(y_right_half_max, ylist, xlist) - np.interp(y_left_half_max, ylist, xlist)
                rel_dist.append(range_x)
            return rel_dist
        
        # Function to find the closest value in a list
        def closest(lst, K):
            return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]
        
        print(get_rel_dist(valley_indexes1, non_singleton_peaks, list(x), list(y)))
