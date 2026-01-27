import numpy as np
import matplotlib.pyplot as plt
import json
import os
from scipy import optimize

from utils.utils import  ChiSqCalc

#script to determine the FCCD of a detector using the count ratio observable

def CalculateFCCD(O_Am241_sim_list, O_Am241_err_sim_list, O_Am241_data, O_Am241_data_err, OutputFileID, FCCD_list, cuts, OutPath):

    print("working directory: ", OutPath)
    dir = OutPath+"/FCCD/am_HS6/"

    if not os.path.exists(dir+"plots"):
        os.makedirs(dir+"plots")
    
    if cuts == "False":
        cuts = False
    else:
        cuts = True


    print("start...")

    O_Am241_tot_err_sim_list = []
    O_Am241_corr_err_sim_list = []
    O_Am241_uncorr_err_sim_list = []

    for obaservable, error in zip(O_Am241_sim_list, O_Am241_err_sim_list):
        O_Am241_tot_err, O_Am241_corr_err, O_Am241_uncorr_err = uncertainty(obaservable, error)
        O_Am241_tot_err_sim_list.append(O_Am241_tot_err)
        O_Am241_corr_err_sim_list.append(O_Am241_corr_err)
        O_Am241_uncorr_err_sim_list.append(O_Am241_uncorr_err)


    #========= PLOTTING ===========

    #plot and fit exp decay
    xdata, ydata = np.array(FCCD_list), np.array(O_Am241_sim_list)
    y_err=O_Am241_tot_err_sim_list*ydata/100 #get absolute error, not percentage
    
    aguess = max(ydata)
    bguess = 1
    p_guess = [aguess,bguess]

    popt, pcov = optimize.curve_fit(exponential_decay, xdata, ydata, p0=p_guess, sigma = y_err, maxfev = 10**7, method ="trf")
    a,b = popt[0],popt[1]
    a_err, b_err, c_err = np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1]), np.sqrt(pcov[2][2])
    chi_sq, p_value, residuals, dof = ChiSqCalc(xdata, ydata, y_err, exponential_decay, popt)

    fig, ax = plt.subplots()
    plt.errorbar(xdata, ydata, xerr=0, yerr =y_err, label = "simulations", elinewidth = 1, fmt='x', ms = 3.0, mew = 3.0)
    xfit = np.linspace(min(xdata), max(xdata), 1000)
    yfit = exponential_decay(xfit,a,b)
    plt.plot(xfit, yfit, "g", label = "fit: a*exp(-bx)")


    #=====fit exp decay of error bars========

    #fit exp decay of total error bars:
    y_uplim = ydata+y_err
    p_guess_up = [max(y_uplim), 1]
    popt_up, pcov_up= optimize.curve_fit(exponential_decay, xdata, y_uplim, p0=p_guess_up, maxfev = 10**7, method ="trf")
    #plt.plot(xfit, yfit_up, color='grey', linestyle='dashed', linewidth=1)

    y_lowlim = ydata-y_err
    p_guess_low = [max(y_lowlim), 1]
    popt_low, pcov_low = optimize.curve_fit(exponential_decay, xdata, y_lowlim, p0=p_guess_low, maxfev = 10**7, method ="trf") 
    yfit_low = exponential_decay(xfit,*popt_low)
    #plt.plot(xfit, yfit_low, color='grey', linestyle='dashed', linewidth=1)

    a_up, b_up= popt_up[0], popt_up[1]
    a_low, b_low = popt_low[0], popt_low[1]


    #fit exp decay of correlated error bars:
    y_corr_err = O_Am241_corr_err_sim_list*ydata/100 
    y_uplim = ydata+y_corr_err
    p_guess_up = [max(y_uplim), 1]
    popt_up_corr, pcov_up_corr = optimize.curve_fit(exponential_decay, xdata, y_uplim, p0=p_guess_up, maxfev = 10**7, method ="trf") 
    yfit_up = exponential_decay(xfit,*popt_up_corr)
    #plt.plot(xfit, yfit_up, color='grey', linestyle='dashed', linewidth=1)

    y_lowlim = ydata-y_corr_err
    p_guess_low = [max(y_lowlim), 1]
    popt_low_corr, pcov_low_corr = optimize.curve_fit(exponential_decay, xdata, y_lowlim, p0=p_guess_low, maxfev = 10**7, method ="trf")
    yfit_low = exponential_decay(xfit,*popt_low_corr)
    #plt.plot(xfit, yfit_low, color='grey', linestyle='dashed', linewidth=1)

    a_up_corr, b_up_corr = popt_up_corr[0], popt_up_corr[1]
    a_low_corr, b_low_corr = popt_low_corr[0], popt_low_corr[1]


    #fit exp decay of uncorrelated error bars:
    y_uncorr_err= O_Am241_uncorr_err_sim_list*ydata/100
    y_uplim = ydata+y_uncorr_err
    p_guess_up = [max(y_uplim), 1]
    popt_up_uncorr, pcov_up_uncorr = optimize.curve_fit(exponential_decay, xdata, y_uplim, p0=p_guess_up, maxfev = 10**7, method ="trf")
    yfit_up = exponential_decay(xfit,*popt_up_uncorr)
    plt.plot(xfit, yfit_up, color='grey', linestyle='dashed', linewidth=1)

    y_lowlim = ydata-y_uncorr_err
    p_guess_low = [max(y_lowlim), 1]
    popt_low_uncorr, pcov_low_uncorr = optimize.curve_fit(exponential_decay, xdata, y_lowlim, p0=p_guess_low, maxfev = 10**7, method ="trf")
    yfit_low = exponential_decay(xfit,*popt_low_uncorr)
    plt.plot(xfit, yfit_low, color='grey', linestyle='dashed', linewidth=1)

    a_up_uncorr, b_up_uncorr = popt_up_uncorr[0], popt_up_uncorr[1]
    a_low_uncorr, b_low_uncorr = popt_low_uncorr[0], popt_low_uncorr[1]

    
    #=========Compute FCCD and all errors ===============
    
    #calculate FCCD of data - invert eq
    FCCD_data = invert_exponential(O_Am241_data,a,b) #(1/b)*np.log(a/(O_Am241_data))
    print('FCCD of data extrapolated: ')
    print(str(FCCD_data))

    #calculate total error on FCCD
    FCCD_err_up = invert_exponential((O_Am241_data-O_Am241_data_err), a_up, b_up)-FCCD_data
    FCCD_err_low = FCCD_data - invert_exponential((O_Am241_data+O_Am241_data_err), a_low, b_low)
    print('total error:  + '+ str(FCCD_err_up) +" - "+str(FCCD_err_low))

    #calculate uncorrelated error on FCCD
    FCCD_err_uncorr_up = invert_exponential((O_Am241_data-O_Am241_data_err), a_up_uncorr, b_up_uncorr)-FCCD_data    
    FCCD_err_uncorr_low = FCCD_data - invert_exponential((O_Am241_data+O_Am241_data_err), a_low_uncorr, b_low_uncorr)   
    print('uncorrelated error:  + '+ str(FCCD_err_uncorr_up) +" - "+str(FCCD_err_uncorr_low))

    #calculate correlated error on FCCD
    FCCD_err_corr_up = invert_exponential((O_Am241_data-O_Am241_data_err), a_up_corr, b_up_corr)-FCCD_data   
    FCCD_err_corr_low = FCCD_data - invert_exponential((O_Am241_data+O_Am241_data_err), a_low_corr, b_low_corr)   
    print('correlated error:  + '+ str(FCCD_err_corr_up) +" - "+str(FCCD_err_corr_low))

   
   #=============Complete Plot===================

    props = dict(boxstyle='round', alpha=0.5)
    info_str = '\n'.join((r'$a=%.3f \pm %.3f$' % (a, np.sqrt(pcov[0][0])), r'$b=%.3f \pm %.3f$' % (b, np.sqrt(pcov[1][1])), r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof), r'FCCD_data=$%.2f^{+%.2f}_{-%.2f}$ mm' % (FCCD_data, FCCD_err_up, FCCD_err_low)))
    plt.text(0.625, 0.775, info_str, transform=ax.transAxes, fontsize=9,verticalalignment='top', bbox=props) #ax.text..ax.tra

    #plot data line
    plt.hlines(O_Am241_data, 0, FCCD_list[-1], colors="orange", label = 'data')
    plt.plot(xfit, [O_Am241_data+O_Am241_data_err]*(len(xfit)), label = 'bounds', color = 'grey', linestyle = 'dashed', linewidth = '1.0')
    plt.plot(xfit, [O_Am241_data-O_Am241_data_err]*(len(xfit)), color = 'grey', linestyle = 'dashed', linewidth = '1.0')


    plt.vlines(FCCD_data, 0, O_Am241_data, colors='orange', linestyles='dashed')
    plt.vlines(FCCD_data+FCCD_err_up, 0, O_Am241_data-O_Am241_data_err, colors='grey', linestyles='dashed', linewidths=1)
    plt.vlines(FCCD_data-FCCD_err_low, 0, O_Am241_data+O_Am241_data_err, colors='grey', linestyles='dashed', linewidths=1)

    # plot syst MC/total corr error line
    #plt.vlines(FCCD_data+FCCD_err_corr_up, 0, O_Am241_data, color="red", linestyles='dashed', linewidths=1, label = "FCCD corr err")
    #plt.vlines(FCCD_data-FCCD_err_corr_low, 0, O_Am241_data, color="red", linestyles='dashed', linewidths=1)

    # #plot stat MC + stat data/ total uncorr:
    #plt.vlines(FCCD_data+FCCD_err_uncorr_up, 0, O_Am241_data-O_Am241_data_err, color="green", linestyles='dashed', linewidths=1, label = "FCCD uncorr err")
    #plt.vlines(FCCD_data-FCCD_err_uncorr_low, 0, O_Am241_data-O_Am241_data_err, color="green", linestyles='dashed', linewidths=1)


    plt.ylabel(r'$O_{am\_HS1}=\frac{C_{60keV}}{C_{99keV}+C_{103keV}}$')
    plt.xlabel("FCCD [mm]")
    plt.xlim(0,FCCD_list[-1])
    plt.ylim(0,1000)
    plt.legend(loc="upper right", fontsize=8)
    plt.show()


    if cuts == False:
         plt.savefig(f"{dir}plots/FCCD-{OutputFileID}.png")
    else:
         plt.savefig(f"{dir}plots/FCCD-{OutputFileID}-cuts.png")

    
    #Save interpolated fccd for data to a json file
    FCCD_data_dict = {
        "FCCD": FCCD_data,

        "FCCD_err_up": FCCD_err_up,
        "FCCD_err_low": FCCD_err_low,
        
        "FCCD_uncorr_err_up": FCCD_err_uncorr_up,
        "FCCD_uncorr_err_low": FCCD_err_uncorr_low,
        
        "FCCD_corr_err_up": FCCD_err_corr_up,
        "FCCD_corr_err_low": FCCD_err_corr_low,
        
        "O_Am241_data": O_Am241_data,
        "O_Am241_data_err": O_Am241_data_err,

        "a": a,
        "a_err": (a_up-a_low)/2,
        "a_corr_err": (a_up_corr-a_low_corr)/2,
        "a_uncorr_err": (a_up_uncorr-a_low_uncorr)/2,
        "b": b,
        "b_err": b_err,
        "b_corr_err": (b_up_corr-b_low_corr)/2,
        "b_uncorr_err": (b_up_uncorr-b_low_uncorr)/2
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

def invert_exponential(x,a,b):
    return (1/b)*np.log(a/x)


def uncertainty(O_Am241, O_Am241_err):

    #values from Valentina's thesis - Am source
    #all percentages
    gamma_line=1.81
    geant4=2.
    source_thickness=0.01
    source_material=0.01
    endcap_thickness=0.37
    detector_cup_thickness=0.03
    detector_cup_material=0.01

    MC_statistics = O_Am241_err/O_Am241*100

    #sum squared of all the contributions
    tot_error=np.sqrt(gamma_line**2+geant4**2+source_thickness**2+source_material**2+endcap_thickness**2+detector_cup_thickness**2+detector_cup_material**2+MC_statistics**2)

     #correlated error
    corr_error=np.sqrt(gamma_line**2+geant4**2+source_thickness**2+source_material**2+endcap_thickness**2+detector_cup_thickness**2+detector_cup_material**2)


    #uncorrelated error
    uncorr_error=MC_statistics
    
   
    return tot_error, corr_error, uncorr_error #NB: tot_error is a percentage error

