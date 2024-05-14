# iterate over all 1691 HARs/DAs with significant interactions called, and store the KDE maxima for each in a dictionary
final_dict = {}
for i in new_lst:
    HAR_name = i
    # group all fragments corresponding to a single HAR/DA together
    is_hacns457_1 =  df['sourceName']==i
    is_hacns457_2 =  df['sourceName']==i+'_sec'
    is_hacns457_3 =  df['sourceName']==i+'_tert' 
    is_hacns457_4 =  df['sourceName']==i+'_quad'
    is_hacns457 = list(map(any, zip(is_hacns457_1,is_hacns457_2,is_hacns457_3,is_hacns457_4)))
    # subset to rows containing scores for only that HAR/DA
    df_hacns457 = df[is_hacns457]
    # compute the relative distance between HAR/DA and other end (in kb)
    df_hacns457['rel_dist (kb)'] = '0'
    a1 = (df_hacns457['targetStart']+df_hacns457['targetEnd'])/2
    b1 = (df_hacns457['sourceStart']+df_hacns457['sourceEnd'])/2
    df_hacns457['rel_dist (kb)'] = (a1-b1)/1000
    
    # scatter plot of interaction strength
    ax1=df_hacns457.plot.scatter('rel_dist (kb)', 'value', c='value', colormap='cool',sharex=False,title=i+'-ortho')
    
    # align axes of KDE plot with scatter plot
    ax2 = ax1.twinx()
    
    minimum = min(list(df_hacns457['rel_dist (kb)'])) # lowest score of interaction strength
    maximum = max(list(df_hacns457['rel_dist (kb)'])) # corresponding high score
    bandwidth = 30
    
    # thresholding for bandwidth - if the interaction space of a HAR/DA is very large, take higher bandwidth, but everything is between bandwidth of 30-40
    if maximum - minimum > 800:
        bandwidth = 40
    elif maximum - minimum > 500 and maximum - minimum < 800:
        bandwidth = 35
    else:
        bandwidth = 30
        
    # plot the KDE plot with variable bandwidth
    ax_d=sns.kdeplot((df_hacns457['rel_dist (kb)']), ax=ax2,bw=bandwidth, label="Density distribution", color='r')
    
    # format y-axis to keep only colorbar
    ax1.set_ylabel("")
    ax1.set_xlabel("Relative distance to HAR (kb)")
    ticks = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y))
    ax2.yaxis.set_major_formatter(ticks)
    ax2.set_yticks([], [])
    ax1.set_yticks([], [])
    
    # get colorbar and name it correctly
    cax = plt.gcf().get_axes()[1]
    plt.legend(loc='best', ncol=1, prop={'size': 8.8})
    cax.set_ylabel('CHiCAGO Score')
    
    # save picture
    #plt.savefig('cNSC_profiles_P1_new/'+i+'-ortho_plot.eps', format='eps')
    
    # get KDE curve x- and y- coordinates
    kde_curve = ax_d.lines[0]
    # get x and y values of this curve
    x = kde_curve.get_xdata()
    y = kde_curve.get_ydata()
    # find indices for maxima of interaction strength 
    indices = find_peaks(y)
    
    # for minima, we flip the sign of the curve and run the function
    inv_data_y = y*(-1) # Tried 1/data_y but not better.

    non_singleton_peaks=[] # x values of non-singleton peaks

    # filter 1: peaks must have >=2 significant contacts in the vicinity of the peak called as the KDE (non-singleton)
    # filter 2: peaks must have a significant contacts in the vicinity of the peak with Chicago score > 5 (soft threshold 4.5)
    if (len(indices[0]) == 1):
        for i in list(df_hacns457["value"]):
            if i >= 4.5 and x[indices[0][0]] not in non_singleton_peaks:
                non_singleton_peaks.append(x[indices[0][0]])
                
    # avoid duplicate values called when just one cluster of points (one peak)       
    elif (len(indices[0]) == 0):
        non_singleton_peaks=[]
    
    # concept: take maxima and two minima closest to the maxima, scan this window to find data points, if the number > 2 then keep this maxima
    else :
        valley_indexes = find_peaks(inv_data_y)
        
        valley_indexes1 = [] # list of x values of minima
        for i in valley_indexes[0]:
            valley_indexes1.append(x[i])
            
        

        # maximas cannot be minimas so remove them: just a sanity check
        #for i in indices:
        #    if x[i] in valley_indexes:
        #        valley_indexes.remove(x[i])

        #adding the two outer x values
        if not (valley_indexes1[0] < x[indices[0][0]]):
            valley_indexes1.insert(0, x[0])
        if not (valley_indexes1[-1] > x[indices[0][-1]]):
            valley_indexes1.append(x[-1])

        # map the x and y values of scatter plot to a single dictionary
        x_y_mapping_dict = dict(zip(list(df_hacns457["rel_dist (kb)"]), df_hacns457["value"]))

        # finally put everything together and make non-singleton peak list for every HAR/DA
        for i in range(len(valley_indexes1)-1):
            small = valley_indexes1[i]
            big = valley_indexes1[i+1]
            count = 0;
            a_value_more_than_five = False
            for j in list(df_hacns457["rel_dist (kb)"]):
                if j >= small and j <= big:
                    count+=1
                    if x_y_mapping_dict.get(j, 1) >= 4.5:
                        a_value_more_than_five = True
            if count>1 and a_value_more_than_five:
                non_singleton_peaks.append(x[indices[0][i]])

                
        
        #df_hacns457['rel_dist (kb)']
        #for i in indices[0]:
        #     print(x[i])
        #for i in valley_indexes[0]:
        #    print(x[i])
        #print(non_singleton_peaks)
        #print(x_y_mapping_dict)
        #df_hacns457["value"]
        
        
        # NEW SCUM BAGGERY
        
        # all x values of minima: valley_indexes1
        
#         print(valley_indexes1)
#         print(non_singleton_peaks)
        
        def get_rel_dist(minima, peaks, xlist, ylist):
            rel_dist = [] # relative distance to HAR range for each peak
            min_index = 0
            for peak in peaks:
                while  min_index >= len(minima) - 1 and not (peak < minima[min_index] and peak > minima[min_index + 1]):
                    min_index += 1
                y_left_half_max = 0.5 * (np.interp(minima[min_index], xlist, ylist) + np.interp(peak, xlist, ylist))
                y_right_half_max = 0.5 * (np.interp(peak, xlist, ylist) + np.interp(minima[min_index + 1], xlist, ylist))
                range_x = np.interp(y_right_half_max, ylist, xlist) - np.interp(y_left_half_max, ylist, xlist)
                print(y_left_half_max, y_right_half_max, range_x)
                rel_dist.append(range_x)
            return rel_dist
        
        def closest (lst, K):
            return lst [min (range(len(lst)), key = lambda i: abs (lst[i]-K)) ]
        
        print(get_rel_dist(valley_indexes1, non_singleton_peaks, list(x), list(y)))