import sys
import numpy as np
import matplotlib.pyplot as plt
import json
import os
from scipy import optimize

from utils.utils import  ChiSqCalc

#script to determine the FCCD of a detector using the count ratio observable

def CalculateFCCD(O_Ba133_sim_list, O_Ba133_err_sim_list, O_Ba133_data, O_Ba133_data_err, OutputFileID, FCCD_list, cuts, OutPath):

    print("working directory: ", OutPath)
    dir = OutPath+"/FCCD/am_HS1/"

    if not os.path.exists(dir+"plots"):
        os.makedirs(dir+"plots")
    
    if cuts == "False":
        cuts = False
    else:
        cuts = True
    
    print("start...")

   
    O_Ba133_tot_err_sim_list = [] 
    O_Ba133_corr_err_sim_list = [] #corr is actually all systematic errors
    O_Ba133_uncorr_err_sim_list = [] #uncorr is actually only statistical error

    for obaservable, error in zip(O_Ba133_sim_list, O_Ba133_err_sim_list):
        O_Ba133_tot_err, O_Ba133_corr_err, O_Ba133_uncorr_err = uncertainty(obaservable, error)
        O_Ba133_tot_err_sim_list.append(O_Ba133_tot_err)
        O_Ba133_corr_err_sim_list.append(O_Ba133_corr_err)
        O_Ba133_uncorr_err_sim_list.append(O_Ba133_uncorr_err)

    #========= PLOTTING ===========

    plot_colors = {"MC":"black", "MC_fit": "black", "data": "orange", "data_err_stat": "green", "MC_err_total": "red", "FCCD": "orange", "FCCD_err_total": "blue", "FCCD_err_MCsyst": "pink", "FCCD_err_statMCstatdata": "grey"}

    #plot and fit exp decay
    print("fitting exp decay")
    xdata, ydata = np.array(FCCD_list), np.array(O_Ba133_sim_list)
    yerr = O_Ba133_tot_err_sim_list*ydata/100 #absolute total error
    aguess = max(ydata)
    bguess = 1
    cguess = min(ydata)
    p_guess = [aguess,bguess,cguess]
    popt, pcov = optimize.curve_fit(exponential_decay, xdata, ydata, p0=p_guess, sigma = yerr,absolute_sigma=False, maxfev = 10**7, method ="trf")
    a,b,c = popt[0],popt[1],popt[2]
    a_err, b_err, c_err = np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1]), np.sqrt(pcov[2][2])
    chi_sq, p_value, residuals, dof = ChiSqCalc(xdata, ydata, yerr, exponential_decay, popt)

    fig, ax = plt.subplots()
    plt.errorbar(xdata, ydata, xerr=0, yerr =yerr, label = "simulations", color= plot_colors["simulations"], elinewidth = 1, fmt='x', ms = 3.0, mew = 3.0)
    xfit = np.linspace(min(xdata), max(xdata), 1000)
    yfit = exponential_decay(xfit,*popt)
    plt.plot(xfit, yfit, color=plot_colors["MC_fit"], label = "fit: a*exp(-bx)+ c")


    #=====fit exp decay of error bars========

    # MC stat and syst
    print("fitting exp decay of total MC O_Ba133 errors")
    y_uplim = ydata+yerr
    p_guess_up = [max(y_uplim), 1, min(y_uplim)]
    popt_up, pcov_up = optimize.curve_fit(exponential_decay, xdata, y_uplim, p0=p_guess_up, maxfev = 10**7, method ="trf") 
    yfit_up = exponential_decay(xfit,*popt_up)
    plt.plot(xfit, yfit_up, color=plot_colors["MC_err_total"], linestyle='dashed', linewidth=1, label="MC err (stat/corr + syst/uncorr)")

    y_lowlim = ydata-yerr
    p_guess_low = [max(y_lowlim), 1, min(y_lowlim)]
    popt_low, pcov_low = optimize.curve_fit(exponential_decay, xdata, y_lowlim, p0=p_guess_low, maxfev = 10**7, method ="trf") 
    yfit_low = exponential_decay(xfit,*popt_low)
    plt.plot(xfit, yfit_low, color=plot_colors["MC_err_total"], linestyle='dashed', linewidth=1)

    a_up, b_up, c_up = popt_up[0], popt_up[1], popt_up[2]
    a_low, b_low, c_low = popt_low[0], popt_low[1], popt_low[2] 


    # MC corr
    print("fitting exp decay of MC stat") 
    yerr_corr = np.array(O_Ba133_corr_err_sim_list)*ydata/100
    y_uplim_corr = ydata+yerr_corr
    p_guess_up = [max(y_uplim), 1, min(y_uplim)]
    popt_up_corr, pcov_up_corr = optimize.curve_fit(exponential_decay, xdata, y_uplim_corr, p0=p_guess_up, maxfev = 10**7, method ="trf") 
    yfit_up_corr = exponential_decay(xfit,*popt_up_corr)
    # plt.plot(xfit, yfit_up_corr, color=plot_colors["MC_err_corr"], linestyle='dashed', linewidth=1, label="MC stat")

    y_lowlim_corr = ydata-yerr_corr
    p_guess_low = [max(y_lowlim), 1, min(y_lowlim)]
    popt_low_corr, pcov_low_corr = optimize.curve_fit(exponential_decay, xdata, y_lowlim_corr, p0=p_guess_low, maxfev = 10**7, method ="trf") 
    yfit_low_corr = exponential_decay(xfit,*popt_low_corr)
    # plt.plot(xfit, yfit_low_corr, color=plot_colors["MC_err_corr"],linestyle='dashed', linewidth=1)

    a_up_corr, b_up_corr, c_up_corr = popt_up_corr[0], popt_up_corr[1], popt_up_corr[2]
    a_low_corr, b_low_corr, c_low_corr = popt_low_corr[0], popt_low_corr[1], popt_low_corr[2] 

    # MC uncorr
    print("fitting exp decay of MC uncorr err")
    yerr_uncorr = O_Ba133_uncorr_err_sim_list*ydata/100
    y_uplim_uncorr = ydata+yerr_uncorr
    p_guess_up = [max(y_uplim), 1, min(y_uplim)]
    popt_up_uncorr, pcov_up_uncorr = optimize.curve_fit(exponential_decay, xdata, y_uplim_uncorr, p0=p_guess_up, maxfev = 10**7, method ="trf") 
    yfit_up_uncorr = exponential_decay(xfit,*popt_up_uncorr)
    # plt.plot(xfit, yfit_up_uncorr, color=plot_colors["MC_err_uncorr"], linestyle='dashed', linewidth=1, label="MC uncorr")

    y_lowlim_uncorr = ydata-yerr_uncorr
    p_guess_low = [max(y_lowlim), 1, min(y_lowlim)]
    popt_low_uncorr, pcov_low_uncorr = optimize.curve_fit(exponential_decay, xdata, y_lowlim_uncorr, p0=p_guess_low, maxfev = 10**7, method ="trf") 
    yfit_low_uncorr = exponential_decay(xfit,*popt_low_uncorr)
    # plt.plot(xfit, yfit_low_uncorr, color=plot_colors["MC_err_uncorr"],linestyle='dashed', linewidth=1)

    a_up_uncorr, b_up_uncorr, c_up_uncorr = popt_up_uncorr[0], popt_up_uncorr[1], popt_up_uncorr[2]
    a_low_uncorr, b_low_uncorr, c_low_uncorr = popt_low_uncorr[0], popt_low_uncorr[1], popt_low_uncorr[2] 


    #=========Compute FCCD and all errors ===============

    #calculate FCCD of data - invert eq
    FCCD_data = invert_exponential(O_Ba133_data,a,b,c)
    print('FCCD of data extrapolated: ')
    print(str(FCCD_data))

    #TOTAL ERROR = (statistical (uncorr) error on data, statistical (uncorr) and systematic (corr) on MC)
    FCCD_err_total_up = invert_exponential((O_Ba133_data-O_Ba133_data_err),a_up, b_up,c_up) - FCCD_data
    FCCD_err_total_low = FCCD_data - invert_exponential((O_Ba133_data+O_Ba133_data_err),a_low, b_low,c_low)
    print('Total error:  + '+ str(FCCD_err_total_up) +" - "+str(FCCD_err_total_low))
    

    #TOTAL CORR ERROR == syst/corr MC only
    FCCD_err_systMC_up = invert_exponential(O_Ba133_data, a_up_corr, b_up_corr, c_up_corr) -FCCD_data
    FCCD_err_systMC_low = FCCD_data - invert_exponential(O_Ba133_data, a_low_corr, b_low_corr, c_low_corr)
    print('TOTAL CORR ERROR == systematic/corr error on MC only')
    print("+ "+str(FCCD_err_systMC_up) +" - "+str(FCCD_err_systMC_low))
    FCCD_err_corr_up, FCCD_err_corr_low = FCCD_err_systMC_up, FCCD_err_systMC_low

    #TOTAL UNCORR ERROR == stat/uncorr MC + stat/uncorr data
    FCCD_err_statMCstatdata_up = invert_exponential(O_Ba133_data-O_Ba133_data_err, a_up_uncorr, b_up_uncorr, c_up_uncorr) -FCCD_data
    FCCD_err_statMCstatdata_low = FCCD_data - invert_exponential(O_Ba133_data+O_Ba133_data_err, a_low_uncorr, b_low_uncorr, c_low_uncorr)
    print('TOTAL UNCORR ERROR == stat/uncorr error on MC and stat/uncorr error on data')
    print("+ "+str(FCCD_err_statMCstatdata_up) +" - "+str(FCCD_err_statMCstatdata_low))
    FCCD_err_uncorr_up, FCCD_err_uncorr_low = FCCD_err_statMCstatdata_up, FCCD_err_statMCstatdata_low

    print("sum in quadrature of corr and uncorr error:")
    FCCD_err_corr_uncorr_quadrature_up = np.sqrt(FCCD_err_uncorr_up**2 + FCCD_err_corr_up**2)
    FCCD_err_corr_uncorr_quadrature_low = np.sqrt(FCCD_err_uncorr_low**2 + FCCD_err_corr_low**2)
    print("+ ", FCCD_err_corr_uncorr_quadrature_up, ", - ", FCCD_err_corr_uncorr_quadrature_low)

    print("linear sum of corr and uncorr error:")
    print("+ ", FCCD_err_corr_up + FCCD_err_uncorr_up,", - ", FCCD_err_corr_low + FCCD_err_uncorr_low)

    #=============Complete Plot===================

    props = dict(boxstyle='round', alpha=0.5)
    info_str = '\n'.join((r'MC fit: $y=a*e^{-bx}+c$',r'$a=%.3f \pm %.3f$' % (a, np.sqrt(pcov[0][0])), r'$b=%.3f \pm %.3f$' % (b, np.sqrt(pcov[1][1])), r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof), r'FCCD_data=$%.3f^{+%.2f}_{-%.2f}$ mm' % (FCCD_data, FCCD_err_total_up, FCCD_err_total_low)))
    plt.text(0.02, 0.98, info_str, transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props) #ax.text..ax.tra

    #plot horizontal data line and errors
    plt.hlines(O_Ba133_data, 0, FCCD_list[-1], color=plot_colors["data"], label = 'data')

    plt.plot(xfit, [O_Ba133_data+O_Ba133_data_err]*(len(xfit)), plot_colors["data_err_stat"], label = 'Data err (stat/uncorr)', linestyle = 'dashed', linewidth = '1.0')
    plt.plot(xfit, [O_Ba133_data-O_Ba133_data_err]*(len(xfit)), plot_colors["data_err_stat"], linestyle = 'dashed', linewidth = '1.0')

    #plot vertical lines
    plt.vlines(FCCD_data, 0, O_Ba133_data, color=plot_colors["FCCD"] , linestyles='dashed')

    #plot total error line
    plt.vlines(FCCD_data+FCCD_err_total_up, 0, O_Ba133_data-O_Ba133_data_err, color=plot_colors["FCCD_err_total"], linestyles='dashed', linewidths=1, label="FCCD total err")
    plt.vlines(FCCD_data-FCCD_err_total_low, 0, O_Ba133_data+O_Ba133_data_err, color=plot_colors["FCCD_err_total"], linestyles='dashed', linewidths=1)

    # plot syst MC/total corr error line
    #plt.vlines(FCCD_data+FCCD_err_corr_up, 0, O_Ba133_data, color=plot_colors["FCCD_err_MCsyst"], linestyles='dashed', linewidths=1, label = "FCCD corr err")
    #plt.vlines(FCCD_data-FCCD_err_corr_low, 0, O_Ba133_data, color=plot_colors["FCCD_err_MCsyst"], linestyles='dashed', linewidths=1)

    # #plot stat MC + stat data/ total uncorr:
    #plt.vlines(FCCD_data+FCCD_err_uncorr_up, 0, O_Ba133_data-O_Ba133_data_err, color=plot_colors["FCCD_err_statMCstatdata"], linestyles='dashed', linewidths=1, label = "FCCD uncorr err")
    #plt.vlines(FCCD_data-FCCD_err_uncorr_low, 0, O_Ba133_data-O_Ba133_data_err, color=plot_colors["FCCD_err_statMCstatdata"], linestyles='dashed', linewidths=1)

    plt.ylabel(r'$O_{Ba133} = \frac{C_{79.6keV}+C_{91keV}}{C_{356keV}}$')
    plt.xlabel("FCCD (mm)")
    plt.xlim(0,FCCD_list[-1])
    plt.ylim(0.0,1.8)
    plt.legend(loc="upper right", fontsize=8)
    plt.tight_layout()

    plt.show()


    if cuts == False:
         plt.savefig(f"{dir}plots/FCCD-{OutputFileID}.png")
    else:
         plt.savefig(f"{dir}plots/FCCD-{OutputFileID}-cuts.png")

    # Save interpolated fccd for data to a json file

    FCCD_data_dict = {
        "FCCD": FCCD_data,

        "FCC_err_up": FCCD_err_total_up,
        "FCCD_err_low": FCCD_err_total_low,

        "FCCD_err_corr_up": FCCD_err_corr_up,
        "FCCD_err_corr_low": FCCD_err_corr_low,

        "FCCD_err_uncorr_up": FCCD_err_uncorr_up,
        "FCCD_err_uncorr_low": FCCD_err_uncorr_low,

        "O_Ba133_data": O_Ba133_data,
        "O_Ba133_data_err": O_Ba133_data_err,

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



def exponential_decay(x, a, b ,c):
    return a*np.exp(-b*x) + c

def invert_exponential(x,a,b,c):
    return (1/b)*np.log(a/(x-c))


def uncertainty(O_Ba133, O_Ba133_err):
    "Returns the total, correlated and uncorrelated % errors on count ratio observable for the MC."

    #MC Systematics
    #values from Bjoern's thesis - Barium source, all percentages
    gamma_line=0.69 
    geant4=2.
    source_thickness=0.02
    source_material=0.01
    endcap_thickness=0.28
    detector_cup_thickness=0.07
    detector_cup_material=0.03

    #MC statistical
    MC_statistics = O_Ba133_err/O_Ba133*100

    #Total error: sum in quadrature of all contributions
    tot_error=np.sqrt(gamma_line**2+geant4**2+source_thickness**2+source_material**2+endcap_thickness**2+detector_cup_thickness**2+detector_cup_material**2+MC_statistics**2)

    correlated = [gamma_line, geant4, source_thickness, source_material, endcap_thickness, detector_cup_thickness, detector_cup_material]
    uncorrelated = [MC_statistics]

    #correlated error
    corr_error = np.sqrt(gamma_line**2+geant4**2+source_thickness**2+source_material**2+endcap_thickness**2+detector_cup_thickness**2+detector_cup_material**2)

    #uncorrelated error
    uncorr_error = MC_statistics
    
    return tot_error, corr_error, uncorr_error




