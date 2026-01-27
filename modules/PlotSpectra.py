import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import glob
import json
import os
import sys


#Script to plot spectra of sim (best fit FCCD) and data

def PlotSpectra(energies_data, energies_sim,  detector, measurement, sim_ID_bestFCCD, peak, xmax, cuts, OutputFileID, OutPath):

    print("detector: ", detector)
    print("sim_id: ", sim_ID_bestFCCD)

    print("working directory: ", OutPath)
    dir = OutPath+"/Spectra/"+detector+"/"+measurement+"/"
    if not os.path.exists(dir):
        os.makedirs(dir)

    #initialise directories to save spectra
    if not os.path.exists(OutPath+"/Spectra/"+detector+"/"+measurement+"/"):
        os.makedirs(dir+"/Spectra/"+detector+"/"+measurement+"/")

    OutputFileName = "DataSim_"+OutputFileID

    print("start...")


    #Get peak counts C_356 for scaling
    PeakCounts_data = glob.glob(OutPath+"/PeakCounts/"+detector+"/"+measurement+"/data/PeakCounts_data*.json")[0]
    PeakCounts_sim = OutPath+"/PeakCounts/"+detector+"/"+measurement+"/sim/PeakCounts_sim_/"+sim_ID_bestFCCD+".json"

    if not os.path.exists(PeakCounts_data):
        print(f"Peak counts data file not found: {PeakCounts_data}")
        sys.exit()
    if not os.path.exists(PeakCounts_sim):
        print(f"Peak counts sim file not found: {PeakCounts_sim}")
        sys.exit()

    with open(PeakCounts_data) as json_file:
        PeakCounts = json.load(json_file)
        peak_data = PeakCounts[peak]

    with open(PeakCounts_sim) as json_file:
        PeakCounts = json.load(json_file)
        peak_sim = PeakCounts[peak]

    print("got peak counts")

    #Plot data and scaled sim
    binwidth = 0.1 #keV
    xmin=0
    bins = np.arange(xmin,xmax,binwidth)


    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex = ax0)

    counts_data, bins, bars_data = ax0.hist(energies_data, bins=bins,  label = "Data", histtype = 'step', linewidth = '0.35')
    counts_sim, bins, bars = ax0.hist(energies_sim, bins = bins, weights=(peak_data/peak_sim)*np.ones_like(energies_sim), label = "sim: FCCD "+str(FCCD)+"mm, DLF: "+str(DLF)+" (scaled)", histtype = 'step', linewidth = '0.35')

    print("basic histos complete")

    Data_sim_ratios = []
    Data_sim_ratios_err = []
    for index, bin in enumerate(bins[1:]):
        data = counts_data[index]
        sim = counts_sim[index] #This counts has already been scaled by weights
        if sim == 0:
            ratio = 0.
            error = 0.
        else:
            try:
                ratio = data/sim
                try:
                    error = np.sqrt(1/data + 1/sim)
                except:
                    error = 0.
            except:
                ratio = 0 #if sim=0 and dividing by 0
        Data_sim_ratios.append(ratio)
        Data_sim_ratios_err.append(error)

    print("errors")

    ax1.errorbar(bins[1:], Data_sim_ratios, yerr=Data_sim_ratios_err,color="green", elinewidth = 1, fmt='x', ms = 1.0, mew = 1.0)
    ax1.hlines(1, xmin, xmax, colors="gray", linestyles='dashed')


    plt.xlabel("Energy [keV]")
    ax0.set_ylabel("Counts")
    ax0.set_yscale("log")
    ax0.legend(loc = "lower left")
    ax0.set_title(detector)
    ax1.set_ylabel("data/sim")
    ax1.set_yscale("log")
    ax1.set_xlim(xmin,xmax)
    ax0.set_xlim(xmin,xmax)

    #Save fig
    if cuts == False:
        plt.savefig(f"{dir}{OutputFileName}.png")
    else:
        plt.savefig(f"{dir}{OutputFileName}-cuts.png")

    print("done")


