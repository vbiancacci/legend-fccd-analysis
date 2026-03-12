import re
import glob
import json
import os
import sys
import pandas as pd
from math import isnan

import pygama.lh5 as lh5

import utils.checks as checks

import modules.GammaLineCounting.GammaLineCounting_am_HS1_data as gl_am1_data
import modules.GammaLineCounting.GammaLineCounting_am_HS6_data as gl_am6_data
import modules.GammaLineCounting.GammaLineCounting_am_HS1_sim as gl_am1_sim
import modules.GammaLineCounting.GammaLineCounting_am_HS6_sim as gl_am6_sim
import modules.GammaLineCounting.GammaLineCounting_ba_HS4 as gl_ba

import modules.CalculateFCCD.CalculateFCCD_am_HS1 as cf_am1
import modules.CalculateFCCD.CalculateFCCD_am_HS6 as cf_am6
import modules.CalculateFCCD.CalculateFCCD_ba_HS4 as cf_ba


import modules.PlotSpectra as ps
import modules.FCCD2AV as fccd2av

from utils.utils import GetEnergyResolutionParameters

currentPath=os.path.dirname(os.path.realpath(__file__))



def GammaLineCounting(ConfigNameFile, data_or_sim, OutPath):
    with open(ConfigNameFile) as json_file: 
        config = json.load(json_file)

    detector       = config['Detector']
    campaign       = config['Campaign']
    measurement    = config['Measurement']
    energy_filter  = config['EnergyFilter']
    cuts           = config['Cuts']
    cut_parameters = config['CutParamenters']

    run, _, source_position= checks.DefineRunPosition(ConfigNameFile)
    meas_ID=f"{detector}-{campaign}-{measurement}-{run}-{source_position}"
    energy_resolution_par = GetEnergyResolutionParameters(detector, campaign, measurement, run, energy_filter)
    
    if data_or_sim=="data":
        datapath = "/global/cfs/cdirs/m2676/data/teststands/hades/prodenv/ref/v1.1.0/generated/tier/hit/"
        #ADD campaign below!
        hit_files = sorted(glob.glob(f"{datapath}/{detector}/{campaign}/{measurement}/char_data-{detector}-{measurement}-{run}*hit.lh5"))
        if len(hit_files)==0:
            print(f"Tier hit not found for this measurement: {detector}/{campaign}/{measurement}/{run}")
            sys.exit()
        else:
            sto = lh5.Store()
        hit_list = []
        if cuts == False:
            for file in hit_files:
                #get data, no cuts
                tb = sto.read_object("hit",file)[0]
                df = lh5.Table.get_dataframe(tb)
                hit_list.append(df)
            dataframe = pd.concat(hit_list, axis=0, ignore_index=True)
            energies = dataframe[energy_filter]
        else:
            for file in hit_files:
                tb = sto.read_object("hit",file)[0]
                df = lh5.Table.get_dataframe(tb)
                hit_list.append(df)
            dataframe = pd.concat(hit_list, axis=0, ignore_index=True)
            energies = dataframe.query(cut_parameters)[energy_filter]
        if measurement[:6]=="am_HS1":
            gl_am1_data.GammaLine_Counting(energies, detector, measurement, meas_ID, energy_filter, energy_resolution_par, cuts, OutPath)
        elif measurement[:6]=="am_HS6":
            gl_am6_data.GammaLine_Counting(energies, detector, measurement, meas_ID, energy_filter, energy_resolution_par, cuts, OutPath)
        else: #ba
            gl_ba.GammaLine_Counting(energies, "data", detector, measurement, meas_ID, energy_filter, energy_resolution_par, cuts, OutPath)
   
    else: #sim
        simpath= "/global/cfs/cdirs/m2676/users/biancacci/hades-sim/legend-g4simple-simulation-myfork/simulations"
        sim_files = sorted(glob.glob(f"{simpath}/{detector}/{campaign}/{measurement}/{run}_{source_position}/hdf5/AV_processed/*hdf5"))
        #print(f"{simpath}/{detector}/{campaign}/{measurement}/{run}_{source_position}/hdf5/AV_processed/*hdf5")
        if len(sim_files)==0:
            print(f"Processed MC simulations not found for this measurement: {detector}/{campaign}/{measurement}/{run}")
            sys.exit()
        for file in sim_files:
            FCCD = re.search(r'FCCD(.{3})', file)
            FCCD = FCCD.group(1)
            DLF = re.search(r'DLF(.{3})', file)
            DLF = DLF.group(1)
            print("FCCD: ", FCCD)
            print("DLF: ", DLF)
            TL_model= "notl"
            frac_FCCDbore= 0.5
            sim_ID=f"{meas_ID}-FCCD{FCCD}mm_DLF{DLF}-{TL_model}-fracFCCDbore{frac_FCCDbore}"
            dataframe =  pd.read_hdf(file, key="procdf")
            energies = dataframe['energy']
            if measurement[:6]=="am_HS1":
                gl_am1_sim.GammaLine_Counting(energies, detector, measurement, sim_ID, energy_resolution_par, OutPath)
            elif measurement[:6]=="am_HS6":
                gl_am6_sim.GammaLine_Counting(energies, detector, measurement, sim_ID, energy_resolution_par, OutPath)
            else: #ba
                gl_ba.GammaLine_Counting(energies, "sim", detector, measurement, sim_ID, None, energy_resolution_par, None, OutPath)

def BestFCCD(ConfigNameFile, option, OutPath):
    with open(ConfigNameFile) as json_file: 
        config = json.load(json_file)

    detector       = config['Detector']
    campaign       = config['Campaign']
    measurement    = config['Measurement']
    energy_filter  = config['EnergyFilter']
    cuts           = config['Cuts']
    cut_parameters = config['CutParamenters']

    dir=OutPath
    print("working directory: ", dir)

    run, _,source_position= checks.DefineRunPosition(ConfigNameFile)
    energy_resolution_par = GetEnergyResolutionParameters(detector, campaign, measurement, run, energy_filter)
    meas_ID=f"{detector}-{campaign}-{measurement}-{run}-{source_position}"
    DLF=1.0
    TL_model= "notl"
    frac_FCCDbore= 0.5
    MC_id=detector+"-"+campaign+"-"+measurement+"-"+run+"-"+source_position+"_"+TL_model+"_fracFCCDbore"+frac_FCCDbore+"_"+energy_filter
    bestFCCD_outputfile = dir+"/FCCD/"+measurement+"/FCCD_data"+MC_id+".json"
    if not os.path.exists(bestFCCD_outputfile):
        print(f"Best FCCD file not found for this measurement: {detector}/{measurement}/{run}")
        print("Run first:      python Process_Detectors.py -fccd confi.file")
        sys.exit()
    
    with open() as outfile:
        FCCD_data = json.load(outfile)
    FCCD = round(FCCD_data["FCCD"],2)
    
    sim_ID_bestFCCD=detector+"-"+campaign+"-"+measurement+"-FCCD"+str(FCCD)+"mm_DLF"+str(DLF)+"-"+TL_model+"-fracFCCDbore"+str(frac_FCCDbore)
    simpath= "/global/cfs/cdirs/m2676/users/biancacci/hades-sim/legend-g4simple-simulation-myfork/simulations"
    sim_file = glob.glob(f"{simpath}/{detector}/{campaign}/{measurement}/{run}_{source_position}/hdf5/{sim_ID_bestFCCD}.hdf5") #AV_processed/*hdf5"))
    if len(sim_file)==0:
        print(f"Processed MC simulations not found for this measurement: {detector}/{campaign}/{measurement}/{run}")
        sys.exit()
    sim_ID=f"{meas_ID}-FCCD{FCCD}mm_DLF{DLF}-{TL_model}-fracFCCDbore{frac_FCCDbore}"
    dataframe_sim =  pd.read_hdf(sim_file[0], key="procdf")
    energies_sim = dataframe_sim['energy']
    if option=="gammaline":
        if measurement[:6]=="am_HS1":
            gl_am1_sim.GammaLine_Counting(energies_sim, detector, measurement, sim_ID, energy_resolution_par, OutPath)
        elif measurement[:6]=="am_HS6":
            gl_am6_sim.GammaLine_Counting(energies_sim, detector, measurement, sim_ID, energy_resolution_par, OutPath)
        else: #ba
            gl_ba.GammaLine_Counting(energies_sim, "sim", detector, measurement, sim_ID, None, energy_resolution_par, None, OutPath)
    else:# option=="plot":
        datapath = "/global/cfs/cdirs/m2676/data/teststands/hades/prodenv/ref/v1.1.0/generated/tier/hit/"
        hit_files = sorted(glob.glob(f"{datapath}/{detector}/{campaign}/{measurement}/char_data-{detector}-{measurement}-{run}*hit.lh5"))
        if len(hit_files)==0:
            print(f"Tier hit not found for this measurement: {detector}/{campaign}/{measurement}/{run}")
            sys.exit()
        else:
            sto = lh5.Store()
        hit_list = []
        if cuts == False:
            for file in hit_files:
                #get data, no cuts
                tb = sto.read_object("hit",file)[0]
                df = lh5.Table.get_dataframe(tb)
                hit_list.append(df)
            dataframe_data = pd.concat(hit_list, axis=0, ignore_index=True)
            energies_data = dataframe_data[energy_filter]
        else: 
            for file in hit_files:
                tb = sto.read_object("hit",file)[0]
                df = lh5.Table.get_dataframe(tb)
                hit_list.append(df)
            dataframe_data = pd.concat(hit_list, axis=0, ignore_index=True)
            energies_data = dataframe_data.query(cut_parameters)[energy_filter]
        if len(hit_files)==0:
            print(f"Tier hit not found for this measurement: {detector}/{campaign}/{measurement}/{run}")
            sys.exit()
        if measurement[:2]=="am":
            peak = "C_60"
            xmax =120
        else: #ba
            peak = "C_356"
            xmax=450
        OutputFileID = f"{sim_ID}-{energy_filter}"
        ps.PlotSpectra(energies_data, energies_sim, detector, measurement, sim_ID_bestFCCD, peak, xmax,cuts, OutputFileID, OutPath)
        
  


def CalculateFCCD(ConfigNameFile, OutPath):
    with open(ConfigNameFile) as json_file: 
        config = json.load(json_file)

    detector       = config['Detector']
    campaign       = config['Campaign']
    measurement    = config['Measurement']
    energy_filter  = config['EnergyFilter']
    cuts           = config['Cuts']

    try:
        TL_model      = config['TL_model']
        frac_FCCDbore = config['frac_FCCDbore']
    except KeyError:
        TL_model      = "notl"
        frac_FCCDbore = 0.5

    run, _, source_position= checks.DefineRunPosition(ConfigNameFile)

    meas_ID=f"{detector}-{campaign}-{measurement}-{run}-{source_position}"
    OutputFileID = f"{meas_ID}-{TL_model}-fracFCCDbore{frac_FCCDbore}-{energy_filter}"

    dir=OutPath
    print("working directory: ", dir)

    #initialise directories to save
    ratio_sim_files = sorted(glob.glob(f"{dir}/PeakCounts/{detector}/{measurement}/sim/PeakCounts_sim-{meas_ID}*-{TL_model}-fracFCCDbore{frac_FCCDbore}.json")) 
    if len(ratio_sim_files)==0:
        print(f"Peak counts simulation files not found for this measurement: {detector}/{measurement}/{run}")
        print("Run first:      python Process_Detectors.py --sim config.json")
        sys.exit()
    if cuts==False:
        ratio_data_file = f"{dir}/PeakCounts/{detector}/{measurement}/data/PeakCounts_data-{meas_ID}-{energy_filter}.json"
    else:
        ratio_data_file = f"{dir}/PeakCounts/{detector}/{measurement}/data/PeakCounts_data-{meas_ID}-{energy_filter}-cuts.json"
    if not os.path.exists(ratio_data_file):
        print(f"Peak counts data file not found for this measurement: {detector}/{measurement}/{run}")
        print("Run first:      python GammaLine_Counting.py --data config.json")
        sys.exit()
  
    observable_sim_list = []
    observable_err_sim_list = []
    FCCD_list = []
    if measurement[:2]=="am":
        key="O_Am241"
        key_err="O_Am241_err"
    elif measurement[:2]=="ba":
        key="O_Ba133"
        key_err="O_Ba133_err"

    #Get count ratio for data
    with open (ratio_data_file) as json_file:
        PeakCounts = json.load(json_file)
    observable_data = PeakCounts[key]
    observable_data_err = PeakCounts [key_err]
    if isnan(observable_data) or isnan(observable_data_err):
        print("Observable equal to NaN. Please check the gamma lines fitting for the data. Run again: python3 ProcessDetector.py -d config.json")
        sys.exit()

    for rf in ratio_sim_files:
        fccd =re.search(r"FCCD([0-9.]+)", rf)
        fccd = float(fccd.group(1))
        #Get count ratio for simulations
        with open(rf) as json_file:
            PeakCounts = json.load(json_file)
        observable = PeakCounts[key]
        observable_err = PeakCounts [key_err]
        if isnan(observable) or isnan(observable_err):
            print(f"Observable equal to NaN for fccd = {fccd}.")
            continue
        else:
            observable_sim_list.append(observable)
            observable_err_sim_list.append(observable_err)
            FCCD_list.append(fccd)

    if measurement[:6]=="am_HS1":
        cf_am1.CalculateFCCD(observable_sim_list, observable_err_sim_list, observable_data, observable_data_err, OutputFileID, FCCD_list, cuts, OutPath)
    elif measurement[:6]=="am_HS6":
        cf_am6.CalculateFCCD(observable_sim_list, observable_err_sim_list, observable_data, observable_data_err, OutputFileID, FCCD_list, cuts, OutPath)
    else: #ba
        cf_ba.CalculateFCCD(observable_sim_list, observable_err_sim_list, observable_data, observable_data_err, FCCD_list, detector, campaign, measurement, run, TL_model, frac_FCCDbore, energy_filter, OutPath)
        




def CalculateAV(ConfigNameFile, OutPath):
    with open(ConfigNameFile) as json_file: 
        config = json.load(json_file)

    detector       = config['Detector']
    campaign       = config['Campaign']
    measurement    = config['Measurement']
    energy_filter  = config['EnergyFilter']
    cuts           = config['Cuts']

    try:
        TL_model      = config['TL_model']
        frac_FCCDbore = config['frac_FCCDbore']
    except KeyError:
        TL_model      = "notl"
        frac_FCCDbore = 0.5

    run, _, source_position= checks.DefineRunPosition(ConfigNameFile)

    meas_ID=f"{detector}-{campaign}-{measurement}-{run}-{source_position}"
    OutputFileID = f"{meas_ID}-{TL_model}-fracFCCDbore{frac_FCCDbore}-{energy_filter}"

    dir=OutPath
    if cuts==False:
        cuts=""
    else:
        cuts="-cuts"
    FCCDFile = sorted(glob.glob(f"{dir}/FCCD/{measurement[:6]}/FCCD-{OutputFileID}{cuts}.json"))
    if len(FCCDFile)==0:
        print(f"FCCD file not found for this measurement: {detector}/{measurement}/{run}")
        print("Run first:      python Process_Detectors.py -fccd config.json")
        sys.exit()
    fccd2av.FromFCCDToAV(FCCDFile[0], detector, measurement, cuts, OutputFileID, OutPath)
    
    
