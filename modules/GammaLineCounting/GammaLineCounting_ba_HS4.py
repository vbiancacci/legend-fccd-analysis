import numpy as np
import os
import matplotlib.pyplot as plt
import json
from scipy.integrate import quad

from pygama.analysis import histograms as hist
from pygama.analysis import peak_fitting

from utils.utils import DefineSigma, ChiSqCalc

#Script to fit the gamma lines in the Ba133 spectra, for data and/or MC


def GammaLine_Counting(energies, data_sim, detector, measurement, ID, energy_filter, energy_resolution_par, cuts, OutPath):
    
    print("ID ", ID)
    print("working directory: ", OutPath)
    dir = OutPath+"/PeakCounts/"+detector+"/"+measurement
    
    if data_sim=="sim":
        dir=dir+"/sim/"
        OutputFileName = f"PeakCounts_sim-{ID}"
    else: #data
        dir=dir+"/data/"
        OutputFileName = f"PeakCounts_data-{ID}-{energy_filter}"
    
    if not os.path.exists(dir+"plots"):
        os.makedirs(dir+"plots")

    

    binwidth = 0.1 #keV

    #_________Fit 79.6/81 double peak____________:

    print("79/81 keV")

    #prepare histogram
    xmin_81, xmax_81 = 77, 84
    if data_sim=="data" and detector == "V09374A":
        xmin_81, xmax_81 = 73, 87
    bins_peak = np.arange(xmin_81,xmax_81,binwidth)
    hist_peak, bins_peak, var_peak = hist.get_hist(energies, bins=bins_peak)
    bins_centres_peak = hist.get_bin_centers(bins_peak)

    zeros = (hist_peak == 0)
    mask = ~(zeros)
    sigma = np.sqrt(hist_peak)[mask]
    hist_peak = hist_peak[mask]
    bins_centres_peak = hist.get_bin_centers(bins_peak)[mask]


    #fit function initial guess
    mu_79_guess, mu_81_guess = 79.6142, 80.9979
    sigma_79_guess = DefineSigma(energy_resolution_par, mu_79_guess)
    sigma_81_guess = DefineSigma(energy_resolution_par, mu_81_guess)
    a_guess,  bkg_guess, s_guess = max(hist_peak), min(hist_peak), min(hist_peak)
    
    double81_guess = [mu_79_guess, mu_81_guess, sigma_79_guess, sigma_81_guess, a_guess, s_guess, bkg_guess]
    try:
        coeff, cov_matrix = peak_fitting.fit_hist(Ba_double_81, hist_peak, bins_peak, var=None, guess=double81_guess, poissonLL=False, integral=None, method=None, bounds=None)
        mu_79, mu_81, sigma_79, sigma_81, a, s, bkg = coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5], coeff[6]
        mu_79_err, mu_81_err, sigma_79_err, sigma_81_err, a_err, s_err, bkg_err = np.sqrt(cov_matrix[0][0]), np.sqrt(cov_matrix[1][1]), np.sqrt(cov_matrix[2][2]), np.sqrt(cov_matrix[3][3]), np.sqrt(cov_matrix[4][4]), np.sqrt(cov_matrix[5][5]), np.sqrt(cov_matrix[6][6])

        #compute chi sq of fit
        chi_sq, p_value, residuals, dof = ChiSqCalc(bins_centres_peak, hist_peak, np.sqrt(hist_peak), Ba_double_81, coeff)
        print("r chi sq: ", chi_sq/dof)

        #Counting - integrate gaussian signal part +/- 3 sigma
        R = 2.65/32.9 #intensity ratio for Ba-133 double peak
        a_79, a_79_err, a_81, a_81_err = R*a, R*a_err, a, a_err
        C_79, C_79_err = gauss_count(a_79,mu_79,sigma_79, a_79_err, binwidth)
        print("peak count 79 = ", str(C_79)," +/- ", str(C_79_err))
        C_81, C_81_err = gauss_count(a_81,mu_81,sigma_81, a_81_err, binwidth)
        print("peak count 81 = ", str(C_81)," +/- ", str(C_81_err))

        #plot with fit
        xfit = np.linspace(xmin_81, xmax_81, 1000)
        yfit = Ba_double_81(xfit, *coeff)
        yfit_step = peak_fitting.step(xfit,mu_79, sigma_79, bkg, s)

        fig, ax = plt.subplots()
        hist.plot_hist(hist_peak, bins_peak, var=None, show_stats=False, stats_hloc=0.75, stats_vloc=0.85)
        plt.plot(xfit, yfit, label=r'Gauss$_{79kev}$ + Gauss$_{81keV}$ + Step$_{79keV}$')
        #plt.plot(xfit, yfit_step, "--", label =r'step($x,\mu_{79},\sigma_{79},bkg,s$))')

        plt.xlim(xmin_81, xmax_81)
        plt.yscale("log")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Counts / 0.1keV")
        plt.legend(loc="upper left", prop={'size': 8.5})

        props = dict(boxstyle='round', alpha=0.5)
        info_str = '\n'.join((r'$\mu_{79}=%.3g \pm %.3g$' % (mu_79, mu_79_err), r'$\mu_{81}=%.3g \pm %.3g$' % (mu_81, mu_81_err), r'$\sigma_{79}=%.3g \pm %.3g$' % (sigma_79, sigma_79_err), r'$\sigma_{81}=%.3g \pm %.3g$' % (sigma_81, sigma_81_err), r'$a=%.3g \pm %.3g$' % (a, a), r'$s=%.3g \pm %.3g$' % (s, s_err), r'$bkg=%.3g \pm %.3g$' % (bkg,bkg_err), r'$R=2.65/32.9$', r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof)))

    except RuntimeError:
        print("Error - curve_fit failed")

        #counting - nan values
        C_79, C_79_err = np.nan, np.nan
        print("peak count 79 = ", str(C_79)," +/- ", str(C_79_err))
        C_81, C_81_err = np.nan, np.nan
        print("peak count 81 = ", str(C_81)," +/- ", str(C_81_err))

        #plot without fit
        fig, ax = plt.subplots()
        hist.plot_hist(hist_peak, bins_peak, var=None, show_stats=False, stats_hloc=0.75, stats_vloc=0.85)
        plt.xlim(xmin_81, xmax_81)
        plt.yscale("log")
        plt.xlabel("Energy [keV]")
        plt.ylabel("Counts / 0.1keV")
        plt.legend(loc="upper left", prop={'size': 8.5})


    #Save fig
    if data_sim=="sim":
        plt.savefig(f"{dir}plots/{OutputFileName}-81eV.png")
    else: #data
        if cuts == False:
            plt.savefig(f"{dir}plots/{OutputFileName}-81keV.png")
        else:
            plt.savefig(f"{dir}plots/{OutputFileName}-81keV-cuts.png")

    #___________Fit single peaks_________________:
    peak_counts = []
    peak_counts_err = []
    peak_ranges = [[352, 359.5],[159,162],[221.5,225],[274,279],[300,306],[381,386.5]] #Rough by eye
    peaks = [356, 161, 223, 276, 303, 383]

    for index, i in enumerate(peak_ranges):

        #prepare histogram
        print(str(peaks[index]), " keV")
        xmin, xmax = i[0], i[1]
        bins_peak = np.arange(xmin,xmax,binwidth)
        hist_peak, bins_peak, var_peak = hist.get_hist(energies, bins=bins_peak)
        bins_centres_peak = hist.get_bin_centers(bins_peak)


        #fit function initial guess
        mu_guess, sigma_guess, a_guess, bkg_guess, s_guess = peaks[index], 1, max(hist_peak), min(hist_peak), min(hist_peak)
        gauss_step_guess = [a_guess, mu_guess, sigma_guess, bkg_guess, s_guess]
        bounds = ([0, 0, 0, -np.inf, 0], [np.inf, np.inf, np.inf, np.inf, np.inf])

        #fit - gauss step
        try:
            coeff, cov_matrix = peak_fitting.fit_hist(peak_fitting.gauss_step, hist_peak, bins_peak, var=None, guess=gauss_step_guess, poissonLL=False, integral=None, method=None, bounds=bounds)
            a, mu, sigma, bkg, s = coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]
            a_err, mu_err, sigma_err, bkg_err, s_err = np.sqrt(cov_matrix[0][0]), np.sqrt(cov_matrix[1][1]), np.sqrt(cov_matrix[2][2]), np.sqrt(cov_matrix[3][3]), np.sqrt(cov_matrix[4][4])

            #compute chi sq of fit
            chi_sq, p_value, residuals, dof = ChiSqCalc(bins_centres_peak, hist_peak, np.sqrt(hist_peak), peak_fitting.gauss_step, coeff)
            print("r chi sq: ", chi_sq/dof)

            #Counting - integrate gaussian signal part +/- 3 sigma
            C, C_err = gauss_count(a,mu,sigma, a_err, binwidth)
            print("peak counts = ", str(C)," +/- ", str(C_err))
            peak_counts.append(C)
            peak_counts_err.append(C_err)

            #plot
            xfit = np.linspace(xmin, xmax, 1000)
            yfit = peak_fitting.gauss_step(xfit, *coeff)
            yfit_step = peak_fitting.step(xfit,mu, sigma, bkg, s)

            fig, ax = plt.subplots()
            hist.plot_hist(hist_peak, bins_peak, var=None, show_stats=False, stats_hloc=0.75, stats_vloc=0.85)
            plt.plot(xfit, yfit, label=r'Total fit: Gauss$_{59.5keV}$ + Step$_{59.5keV}$')
            #plt.plot(xfit, yfit_step, "--", label =r'step($x,\mu,\sigma,bkg,s$)')

            plt.xlim(xmin, xmax)
            plt.yscale("log")
            plt.xlabel("Energy [keV]")
            plt.ylabel("Counts / 0.1keV")
            plt.legend(loc="upper left", prop={'size': 8.5})

            props = dict(boxstyle='round', alpha=0.5)
            info_str = '\n'.join((r'$a=%.3g \pm %.3g$' % (a, a_err), r'$\mu=%.3g \pm %.3g$' % (mu, mu_err), r'$\sigma=%.3g \pm %.3g$' % (sigma, sigma_err), r'$bkg=%.3g \pm %.3g$' % (bkg, bkg_err),r'$s=%.3g \pm %.3g$' % (s, s_err), r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof)))
            #plt.text(0.02, 0.8, info_str, transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)

        except RuntimeError:
            print("Error - curve_fit failed")

            #counting - nan values
            C, C_err = np.nan, np.nan
            print("peak counts = ", str(C)," +/- ", str(C_err))
            peak_counts.append(C)
            peak_counts_err.append(C_err)

            #plot without fit
            fig, ax = plt.subplots()
            hist.plot_hist(hist_peak, bins_peak, var=None, show_stats=False, stats_hloc=0.75, stats_vloc=0.85)
            plt.xlim(xmin, xmax)
            plt.yscale("log")
            plt.xlabel("Energy (keV)")
            plt.ylabel("Counts")
            plt.legend(loc="upper left", prop={'size': 8.5})


        #Save fig
        if data_sim=="sim":
            plt.savefig(f"{dir}plots/{OutputFileName}-{peaks[index]}keV.png")
        else:  #data
            if cuts == False:
                plt.savefig(f"{dir}plots/{OutputFileName}-{peaks[index]}keV.png")
            else:
                plt.savefig(f"{dir}plots/{OutputFileName}-{peaks[index]}keV-cuts.png")


    #Comput count ratio O_Ba133
    print("")
    C_356, C_356_err = peak_counts[0], peak_counts_err[0]
    if (C_356 == np.nan) or (C_79 == np.nan) or (C_81 == np.nan):
        O_Ba133, O_Ba133_err = np.nan, np.nan
    else:
        O_Ba133 = (C_79 + C_81)/C_356
        O_Ba133_err = O_Ba133*np.sqrt((C_79_err**2 + C_81_err**2)/(C_79+C_81)**2 + (C_356_err/C_356)**2)
    print("O_BA133 = " , O_Ba133, " +/- ", O_Ba133_err)


    #Save count values to json file
    PeakCounts = {
        "C_79": C_79,
        "C_79_err" : C_79_err,
        "C_81" : C_81,
        "C_81_err" : C_81_err,
        "C_356" : peak_counts[0],
        "C_356_err" : peak_counts_err[0],
        "C_161" : peak_counts[1],
        "C_161_err" : peak_counts_err[1],
        "C_223" : peak_counts[2],
        "C_223_err" : peak_counts_err[2],
        "C_276" : peak_counts[3],
        "C_276_err" : peak_counts_err[3],
        "C_303" : peak_counts[4],
        "C_303_err" : peak_counts_err[4],
        "C_383" : peak_counts[5],
        "C_383_err" : peak_counts_err[5],
        "O_Ba133" : O_Ba133,
        "O_Ba133_err" : O_Ba133_err,
    }

    if data_sim=="sim":
        with open(f"{dir}{OutputFileName}.json", "w") as outfile:
            json.dump(PeakCounts, outfile, indent=4)
    else:
        if cuts == False:
            with open(f"{dir}{OutputFileName}.json", "w") as outfile:
                json.dump(PeakCounts, outfile, indent=4)
        else:
            with open(f"{dir}{OutputFileName}-cuts.json", "w") as outfile:
                json.dump(PeakCounts, outfile, indent=4)
        


def gauss_count(a,mu,sigma, a_err, bin_width):
    "count/integrate gaussian peak"

    integral = a/bin_width
    integral_err = a_err/bin_width

    #_____3sigma_____
    integral_356_3sigma_list = quad(peak_fitting.gauss,mu-3*sigma, mu+3*sigma, args=(mu,sigma,a))
    integral = integral_356_3sigma_list[0]/bin_width
    integral_err = a_err/bin_width

    return integral, integral_err


def Ba_double_81(x, mu_79, mu_81, sigma_79, sigma_81, a, s, bkg):

    R = 2.65/32.9 #intensity ratio for Ba-133 double peak

    peak_79 = peak_fitting.gauss(x,mu_79,sigma_79,a*R)
    peak_81 = peak_fitting.gauss(x,mu_81,sigma_81,a)

    step = peak_fitting.step(x,mu_79,sigma_79,bkg,s)

    f = peak_79 + peak_81 + step

    return f
