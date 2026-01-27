import numpy as np, glob, sys, yaml
from scipy import stats

def GetEnergyResolutionParameters(detector, campaign, measurement, run, energy_filter):
    ParameterPath = "/global/cfs/cdirs/m2676/users/biancacci/hades-proc/legend-dataflow-hades/generated/par/hit"
    ParameterFileName = ParameterPath + "/" + detector + "/" + measurement + "/char_data-" + detector + "-" + measurement + "-" + run + "*-par_hit.yaml"
    ParameterFileName = ParameterPath + "/V09372A/th_HS2_lat_psa/char_data-V09372A-th_HS2_lat_psa-r001*-par_hit.yaml"
    ParameterFile = glob.glob(ParameterFileName)
    if len(ParameterFile) == 0:
        print("No energy resolution parameters found for the given measurement:", ParameterFileName)
        sys.exit()
        return
    with open(ParameterFile[0], "r") as f:
        docs = yaml.safe_load(f)
        res_parameter_a = docs["results"]["ecal"][energy_filter]["eres_linear"]["parameters"]["a"]
        res_parameter_b = docs["results"]["ecal"][energy_filter]["eres_linear"]["parameters"]["b"]
        energy_res_par = [res_parameter_a, res_parameter_b]
    return energy_res_par


def DefineSigma(energy_resolution_par, energy):
    """Define sigma of gaussian energy resolution at given energy"""
    res_parameter_a = energy_resolution_par[0]
    res_parameter_b = energy_resolution_par[1]
    FWHM = np.sqrt(res_parameter_a ** 2 + res_parameter_b * energy)
    sigma = FWHM / 2.355
    return sigma


def ChiSqCalc(xdata, ydata, yerr, fit_func, coeff):
    """calculate chi sq and p-val of a fit given the data points and fit parameters, e.g. fittype ='linear'"""
    y_obs = ydata
    y_exp = []
    y_exp = []
    for index, y_i in enumerate(y_obs):
        x_obs = xdata[index]
        y_exp_i = fit_func(x_obs, *coeff)
        y_exp.append(y_exp_i)
    else:
        chi_sq = 0.0
        residuals = []
        for i in range(len(y_obs)):
            if yerr[i] != 0:
                residual = (y_obs[i] - y_exp[i]) / yerr[i]
            else:
                residual = 0.0
            chi_sq += residual ** 2
            residuals.append(residual)
        else:
            N = len(y_exp)
            dof = N - 1
            chi_sq_red = chi_sq / dof
            p_value = 1 - stats.chi2.cdf(chi_sq, dof)
            return (
             chi_sq, p_value, residuals, dof)
