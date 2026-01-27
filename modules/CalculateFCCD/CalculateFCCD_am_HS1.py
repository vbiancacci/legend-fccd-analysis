import numpy as np
import matplotlib.pyplot as plt
import json
import os
from scipy import optimize


#script to determine the FCCD of a detector using the count ratio observable, ignoring TL for now

def CalculateFCCD(O_Am241_sim_list, O_Am241_err_sim_list, O_Am241_data, O_Am241_data_err, OutputFileID, FCCD_list, cuts, OutPath):

    print("working directory: ", OutPath)
    dir = OutPath+"/FCCD/am_HS1/"

    if not os.path.exists(dir+"plots"):
        os.makedirs(dir+"plots")
    
    if cuts == "False":
        cuts = False
    else:
        cuts = True


    calibration="ICPC"
    a=0
    a_err=0
    if calibration=="BEGe":
        a=62.07387874149428   
        arr=1.491153644966186 
    else: #IC
        a=59.653315556037235    
        arr= 1.5483673501187423 

    print("start...")

    #plot and fit exp decay
    xdata, ydata = np.array(FCCD_list), np.array(O_Am241_sim_list)
    print(ydata)
    yerr = O_Am241_err_sim_list 

    aguess = max(ydata)
    bguess = 1
    p_guess = [aguess,bguess]
    
    popt, pcov = optimize.curve_fit(exponential_decay, xdata, ydata, p0=p_guess, sigma = yerr, maxfev = 10**7, method ="trf") #, bounds = bounds)
    a_bad,b = popt[0],popt[1]
    xfit = np.linspace(min(xdata), max(xdata), 1000)
    
    fig, ax = plt.subplots()
    yfit = exponential_decay(xfit,a,b)
    plt.plot(xfit, yfit, "g", label = "fit: a*exp(-bx)")


    #fit exp decay of error bars:
    a_up=a+arr
    yfit_up = exponential_decay(xfit,a_up,b)
    plt.plot(xfit, yfit_up, color='grey', linestyle='dashed', linewidth=1)

    a_low=a-arr
    yfit_low = exponential_decay(xfit,a_low,b)
    plt.plot(xfit, yfit_low, color='grey', linestyle='dashed', linewidth=1)


    #calculate FCCD of data - invert eq
    FCCD_data = (1/b)*np.log(a/O_Am241_data)
    print("FCCD ", FCCD_data)

    #caluclate total error on FCCD
    FCCD_err_up = (1/b)*np.log(a_up/(O_Am241_data-O_Am241_data_err))-FCCD_data
    FCCD_err_low = FCCD_data - (1/b)*np.log(a_low/(O_Am241_data+O_Am241_data_err))
    print('total error:  + '+ str(FCCD_err_up) +" - "+str(FCCD_err_low))

    #calculate uncorrelated error on FCCD
    FCCD_uncorr_up = (1/b)*np.log(a/(O_Am241_data-O_Am241_data_err))-FCCD_data
    FCCD_uncorr_low = FCCD_data - (1/b)*np.log(a/(O_Am241_data+O_Am241_data_err))
    print('uncorrelated error:  + '+ str(FCCD_uncorr_up) +" - "+str(FCCD_uncorr_low))

    #calculate correlated error on FCCD
    FCCD_corr_up = (1/b)*np.log(a_up/(O_Am241_data))-FCCD_data
    FCCD_corr_low = FCCD_data - (1/b)*np.log(a_low/(O_Am241_data))
    print('correlated error:  + '+ str(FCCD_corr_up) +" - "+str(FCCD_corr_low))

    props = dict(boxstyle='round', alpha=0.5)
    info_str = '\n'.join((r'$a=%.3f \pm %.3f$' % (a,arr), r'$b=%.3f \pm %.3f$' % (b, np.sqrt(pcov[1][1])), r'FCCD_data=$%.2f^{+%.2f}_{-%.2f}$ mm' % (FCCD_data, FCCD_err_up, FCCD_err_low)))
    plt.text(0.62, 0.775, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

    #plot data line
    plt.hlines(O_Am241_data, 0, FCCD_list[-1], colors="orange", label = 'data')
    plt.plot(xfit, [O_Am241_data+O_Am241_data_err]*(len(xfit)), label = 'bounds', color = 'grey', linestyle = 'dashed', linewidth = '1.0')
    plt.plot(xfit, [O_Am241_data-O_Am241_data_err]*(len(xfit)), color = 'grey', linestyle = 'dashed', linewidth = '1.0')


    plt.vlines(FCCD_data, 0, O_Am241_data, colors='orange', linestyles='dashed')
    plt.vlines(FCCD_data+FCCD_err_up, 0, O_Am241_data, colors='grey', linestyles='dashed', linewidths=1)
    plt.vlines(FCCD_data-FCCD_err_low, 0, O_Am241_data, colors='grey', linestyles='dashed', linewidths=1)


    plt.ylabel(r'$O_{am\_HS1}=\frac{C_{60keV}}{C_{99keV}+C_{103keV}}$', fontsize=20)
    plt.xlabel("FCCD [mm]", fontsize=20)
    ax.tick_params(axis="both", labelsize=20)
    plt.xlim(0,FCCD_list[-1])
    plt.ylim(0,70)
    #plt.title(detector)
    plt.legend(loc="upper right", fontsize=10)
    plt.tight_layout()

    if cuts == False:
         plt.savefig(f"{dir}plots/FCCD-{OutputFileID}.png")
    else:
         plt.savefig(f"{dir}plots/FCCD-{OutputFileID}-cuts.png")


    #Save interpolated fccd for data to a json file
    FCCD_data_dict = {
        "FCCD": FCCD_data,

        "FCCD_err_up": FCCD_err_up,
        "FCCD_err_low": FCCD_err_low,
        
        "FCCD_uncorr_err_up": FCCD_uncorr_up,
        "FCCD_uncorr_err_low": FCCD_uncorr_low,
        "FCCD_corr_err_up": FCCD_corr_up,
        "FCCD_corr_err_low": FCCD_corr_low,
        
        "O_Am241_data": O_Am241_data,
        "O_Am241_data_err": O_Am241_data_err,
        
        "a": a,
        "a_err": (a_up-a_low)/2
    }

    if cuts == False:
        with open(f"{dir}/FCCD-{OutputFileID}.json", "w") as outfile:
                json.dump(FCCD_data_dict, outfile, indent=4)
    else:
        with open(f"{dir}/FCCD-{OutputFileID}-cuts.json", "w") as outfile:
            json.dump(FCCD_data_dict, outfile, indent=4)

    print("done")

def exponential_decay(x, a, b):
    f = a*np.exp(-b*x)
    return f
