import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import stats
import json
from scipy.integrate import quad

from pygama.analysis import histograms as hist
from pygama.analysis import peak_fitting

from utils.utils import DefineSigma, ChiSqCalc

#Script to fit the gamma lines in the Am241 spectra for simulations

def GammaLine_Counting(energies, detector, measurement, sim_ID, energy_resolution_par, OutPath):
    
    print("sim_ID ", sim_ID)
    print("working directory: ", OutPath)
    dir = OutPath+"/PeakCounts/"+detector+"/"+measurement+"/sim/"
    print(dir)

    if not os.path.exists(dir+"plots"):
        os.makedirs(dir+"plots")

    OutputFileName = f"PeakCounts_sim-{sim_ID}"

    binwidth = 0.1 #keV

    #_________Fit 99/103 double peak____________:

    print("99/103 keV")

    #prepare histogram
    xmin_99_103, xmax_99_103 = 97, 105
    bins_peak = np.arange(xmin_99_103,xmax_99_103,binwidth)
    hist_peak, bins_peak, var_peak = hist.get_hist(energies, bins=bins_peak)
    bins_centres_peak = hist.get_bin_centers(bins_peak)

    #fit function initial guess
    R =  0.0203/0.0195
    mu_99_guess, mu_103_guess = 99., 103.
    sigma_99_guess = DefineSigma(energy_resolution_par, mu_99_guess)
    sigma_103_guess = DefineSigma(energy_resolution_par, mu_103_guess)
    a_99_guess, a_103_guess = max(hist_peak)*R, max(hist_peak)
    bkg_99_guess, s_99_guess = min(hist_peak), min(hist_peak)

    double_guess = [a_99_guess, mu_99_guess, sigma_99_guess,
                    a_103_guess, mu_103_guess, sigma_103_guess,
                    bkg_99_guess,s_99_guess]
    try:

        coeff, cov_matrix = peak_fitting.fit_hist(Am_double, hist_peak, bins_peak, var=None, guess=double_guess, poissonLL=False, integral=None, method=None, bounds=None)

        a_99, mu_99, sigma_99, a_103, mu_103, sigma_103 = coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]
        bkg_99, s_99 = coeff[6], coeff[7]
        a_99_err, mu_99_err, sigma_99_err, a_103_err, mu_103_err, sigma_103_err = np.sqrt(cov_matrix[0][0]), np.sqrt(cov_matrix[1][1]), np.sqrt(cov_matrix[2][2]), np.sqrt(cov_matrix[3][3]), np.sqrt(cov_matrix[4][4]), np.sqrt(cov_matrix[5][5])
        bkg_99_err , s_99_err  = np.sqrt(cov_matrix[6][6]), np.sqrt(cov_matrix[7][7])
        #compute chi sq of fit
        chi_sq, p_value, residuals, dof = ChiSqCalc(bins_centres_peak, hist_peak, np.sqrt(hist_peak), Am_double, coeff)
        print("r chi sq: ", chi_sq/dof)

        #Counting - integrate gaussian signal part

        C_99, C_99_err = gauss_count(a_99, mu_99, sigma_99, a_99_err, binwidth)
        C_103, C_103_err = gauss_count(a_103, mu_103, sigma_103, a_103_err, binwidth)
        C_99_103=C_99+C_103
        C_99_103_err=np.sqrt(C_99_err**2+C_103_err**2)
        print("peak count 99-103 = ", str(C_99_103)," +/- ", str(C_99_103_err))

        #plot with fit
        xfit = np.linspace(xmin_99_103, xmax_99_103, 1000)
        yfit = Am_double(xfit, *coeff)
        #yfit_step = peak_fitting.step(xfit,mu_99, sigma_99, bkg_99, s_99)
        #yfit_doubleG= double_gauss(xfit, a_99, mu_99, sigma_99, a_103, mu_103, sigma_103)

        fig, ax = plt.subplots()
        hist.plot_hist(hist_peak, bins_peak, label='Data', var=None, show_stats=False, stats_hloc=0.75, stats_vloc=0.85)
        plt.plot(xfit, yfit, label=r'Gauss$_{99kev}$ + Gauss$_{103keV}$ + Step$_{99keV}$ + Linear')
        #plt.plot(xfit, yfit_doubleG, "--", color='red', label =r'double_gauss($x,\mu_{99},\sigma_{99},a_{99},\mu_{103},\sigma_{103},a_{103}$)')
        #plt.plot(xfit, yfit_step, "--", label =r'step($x,\mu_{99},\sigma_{99},bkg,s$))')


        plt.xlim(xmin_99_103, xmax_99_103)
        plt.yscale("log")
        plt.xlabel("Energy [keV]")
        plt.ylabel("Counts / 0.1keV")
        plt.legend(loc="lower left", prop={'size': 8.5})

        props = dict(boxstyle='round', alpha=0.5)
        info_str = '\n'.join((r'$\mu_{99}=%.3g \pm %.3g$' % (mu_99, mu_99_err), r'$\mu_{103}=%.3g \pm %.3g$' % (mu_103, mu_103_err), r'$\sigma_{99}=%.3g \pm %.3g$' % (sigma_99, sigma_99_err), r'$\sigma_{103}=%.3g \pm %.3g$' % (sigma_103, sigma_103_err), r'$a_{99}=%.3g \pm %.3g$' % (a_99, a_99_err), r'$a_{103}=%.3g \pm %.3g$' % (a_103, a_103_err), r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof)))
        

    except RuntimeError:
        print("Error - curve_fit failed")

        #counting - nan values
        C_99_103, C_99_103_err = np.nan, np.nan
        print("peak count 99-103 = ", str(C_99_103)," +/- ", str(C_99_103_err))

        #plot without fit
        fig, ax = plt.subplots()
        hist.plot_hist(hist_peak, bins_peak, var=None, show_stats=False, stats_hloc=0.75, stats_vloc=0.85)
        plt.xlim(xmin_99_103, xmax_99_103)
        plt.yscale("log")
        plt.xlabel("Energy [keV]")
        plt.ylabel("Counts / 0.1keV")
        plt.legend(loc="upper left", prop={'size': 8.5})


    #Save fig
    plt.savefig(f"{dir}plots/{OutputFileName}-103keV.png")

    #___________Fit single peaks_________________:
    peak_ranges =[52,62]

    #prepare histogram

    xmin, xmax = peak_ranges[0], peak_ranges[1]
    bins_peak = np.arange(xmin,xmax,binwidth)
    hist_peak, bins_peak, var_peak = hist.get_hist(energies, bins=bins_peak)
    bins_centres_peak = hist.get_bin_centers(bins_peak)


    #fit function initial guess

    mu_59_guess, mu_57_guess, mu_53_guess = 59.5,  57.8, 53
    sigma_59_guess = DefineSigma(energy_resolution_par, mu_59_guess)
    sigma_57_guess = DefineSigma(energy_resolution_par, mu_57_guess)
    sigma_53_guess = DefineSigma(energy_resolution_par, mu_53_guess)
    a_59_guess, a_57_guess, a_53_guess = max(hist_peak), max(hist_peak)*0.8, max(hist_peak)*0.5   #1.2, *0.8
    bkg_guess, s_guess = min(hist_peak), min(hist_peak)
    gauss_step_guess = [mu_59_guess, sigma_59_guess, a_59_guess, 
                        mu_57_guess, sigma_57_guess, a_57_guess, 
                        mu_53_guess, sigma_53_guess, a_53_guess, 
                        bkg_guess, s_guess]
    bounds = ([0, 0, 0, 0, 0, 0, 0, 0, 0, -np.inf, 0], [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])

    #fit - gauss step
    try:
        coeff, cov_matrix = peak_fitting.fit_hist(Am_60, hist_peak, bins_peak, var=None, guess=gauss_step_guess, poissonLL=False, integral=None, method=None, bounds=bounds)
        mu_59, sigma_59, a_59 = coeff[0], coeff[1], coeff[2]
        mu_57, sigma_57, a_57= coeff[3], coeff[4], coeff[5]
        mu_53, sigma_53, a_53= coeff[6], coeff[7], coeff[8]
        bkg, s = coeff[9], coeff[10]
        mu_59_err, sigma_59_err, a_59_err = np.sqrt(cov_matrix[0][0]), np.sqrt(cov_matrix[1][1]), np.sqrt(cov_matrix[2][2])
        #compute chi sq of fit
        chi_sq, p_value, residuals, dof = ChiSqCalc(bins_centres_peak, hist_peak, np.sqrt(hist_peak), Am_60, coeff)
        print("r chi sq: ", chi_sq/dof)

        #Counting - integrate gaussian signal part +/- 3 sigma
        C_60, C_60_err = gauss_count(a_59, mu_59 ,sigma_59, a_59_err, binwidth)
        print("peak counts = ", str(C_60)," +/- ", str(C_60_err))

        #plot
        xfit = np.linspace(xmin, xmax, 1000)
        yfit = Am_60(xfit, *coeff)
        #yfit_step = peak_fitting.step(xfit ,mu_59, sigma_59, bkg, s)
        #yfit_gaus53 = peak_fitting.gauss(xfit ,mu_53, sigma_53, a_53)
        #yfit_gaus57 = peak_fitting.gauss(xfit ,mu_57, sigma_57, a_57)
        yfit_gaus59 = peak_fitting.gauss(xfit ,mu_59, sigma_59, a_59)

        fig, ax = plt.subplots()
        hist.plot_hist(hist_peak, bins_peak, label="Data", var=None, show_stats=False, stats_hloc=0.75, stats_vloc=0.85)
        plt.plot(xfit, yfit, label=r'Total fit: Gauss$_{59.5keV}$ + Step$_{59.5keV}$ + Gauss$_{53keV}$ + Gauss$_{57keV}$')
        #plt.plot(xfit, yfit_step, "--", label =r'step($x,\mu,\sigma,bkg,s$)')
        #plt.plot(xfit, yfit_gaus53, "--", color='blue', label =r'gaus53($x,\mu,\sigma,a)')
        #plt.plot(xfit, yfit_gaus57, "--", color='green', label =r'gaus57($x,\mu,\sigma,a)')
        #plt.plot(xfit, yfit_gaus59, "--", color='red', label =r'gauss59: gauss($x,\mu_{59},\sigma_{59},a_{59}$)')

        plt.xlim(xmin, xmax)
        #ax.set_ylim(min(hist_peak)*0.8,max(hist_peak)*1.2)
        plt.yscale("log")
        plt.xlabel("Energy [keV]")
        plt.ylabel("Counts / 0.1keV")
        plt.legend(loc="lower left", prop={'size': 8.5})

        props = dict(boxstyle='round', alpha=0.5)
        info_str = '\n'.join((r'$\mu_{59}=%.3g \pm %.3g$' % (mu_59, mu_59_err), r'$\sigma_{59}=%.3g \pm %.3g$' % (sigma_59, sigma_59_err),r'$a_{59}=%.3g \pm %.3g$' % (a_59, a_59_err), r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof)))

    except RuntimeError:
        print("Error - curve_fit failed")

        #counting - nan values
        C_60, C_60_err = np.nan, np.nan
        print("peak counts = ", str(C_60)," +/- ", str(C_60_err))

        #plot without fit
        fig, ax = plt.subplots()
        hist.plot_hist(hist_peak, bins_peak, var=None, show_stats=False, stats_hloc=0.75, stats_vloc=0.85)
        plt.xlim(xmin, xmax)
        plt.yscale("log")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Counts")
        plt.legend(loc="upper left", prop={'size': 8.5})


    #Save fig
    plt.savefig(f"{dir}plots/{OutputFileName}-59keV.png")

    #Comput count ratio O_Am241
    print("")
    if (C_60 == np.nan) or (C_99_103 == np.nan):
        O_Am241, O_Am241_err = np.nan, np.nan
    else:
        O_Am241 = C_60/C_99_103
        O_Am241_err = np.sqrt((C_60_err/C_60)**2 + (C_99_103_err/C_99_103)**2)*O_Am241
    print("O_Am241 = " , O_Am241, " +/- ", O_Am241_err)


    #Save count values to json file
    PeakCounts = {
        "C_60" : C_60,
        "C_60_err" : C_60_err,
        "C_99_103" : C_99_103,
        "C_99_103_err" : C_99_103_err,
        "O_Am241" : O_Am241,
        "O_Am241_err" : O_Am241_err,
    }

    with open(f"{dir}{OutputFileName}.json", "w") as outfile:
        json.dump(PeakCounts, outfile, indent=4)


def gauss_count(a,mu,sigma, err_a, bin_width):
    "count/integrate gaussian peak"

    integral_list = quad(peak_fitting.gauss,0,120, args=(mu,sigma,a))
    integral = integral_list[0]/bin_width
    integral_err = err_a/bin_width

    return integral, integral_err


def Am_60(x, mu_59, sigma_59, a_59, mu_57, sigma_57, a_57, mu_53, sigma_53, a_53, bkg, s):

    peak_59 = peak_fitting.gauss(x,mu_59, sigma_59, a_59)
    peak_57 = peak_fitting.gauss(x,mu_57, sigma_57, a_57)
    bump_53 = peak_fitting.gauss(x, mu_53, sigma_53, a_53)

    step = peak_fitting.step(x,mu_59,sigma_59,bkg,s)

    f = peak_59 + peak_57 + bump_53 + step

    return f


def Am_double(x,a1,mu1,sigma1,a2,mu2,sigma2,b,s,
              components=False) :

    step2 = peak_fitting.step(x,mu2,sigma2,b,s)

    gaus1 = peak_fitting.gauss(x,mu1,sigma1,a1)
    gaus2 = peak_fitting.gauss(x,mu2,sigma2,a2)
    lin=s*x

    double_f = step2 + gaus1 + gaus2 +lin

    if components:
       return double_f, gaus1, gaus2, step2,
    else:
       return double_f


