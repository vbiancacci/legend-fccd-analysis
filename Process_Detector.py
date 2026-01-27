import sys
import os
import argparse

import tools.tools as tools


CodePath=os.path.dirname(os.path.realpath(__file__))

if len(sys.argv) != 3 and sys.argv[-1]!="-h" and sys.argv[-1]!="--help":
    print('')
    print('#################################################################################')
    print('## Usage: python3 Process_Detector.py [OPERATION MODE] [CONFIG FILE] [OPTIONS]             ##')
    print('## Plese enter "python3 Process_Detectors.py --help" for more information                   ##')
    print('#################################################################################')
    print('')


def main():

    par = argparse.ArgumentParser(description="Gamma line count analysis and FCCD/AV calculation",
                                  epilog="Report bugs to <valentina.biancacci@gssi.it>",
                                  usage="python3 make.py [OPERATION: -l, -m, -p, -s, -a, -v] [CONFIG FILE] [OPTIONS: -o path]"
    )
    arg = par.add_argument
    arg("-d", "--data",  nargs=1, help="Gamma line count data", metavar="")
    arg("-s", "--sim", nargs=1, help="Gamma line count simulation", metavar="")
    arg("-f", "--fccd", nargs=1, help="Calculate FCCD", metavar="")
    arg("-v", "--av", nargs=1, help="Calculate active volume", metavar="")
    arg("-a", "--all",  nargs=1, help="Gamma line count analysis for data and simulations and calculate fccd and active volume", metavar="")
    arg("-b", "--best_fccd",  nargs=1, help="Gamma line count of best FCCD MC", metavar="")
    arg("-p", "--plot",  nargs=1, help="Plot energy spectra with best FCCD", metavar="")
    arg("-o", "--output", nargs=1, help="output path", default=[f"{CodePath}/results"], metavar="")
    
    args=vars(par.parse_args())

    if args["data"]:
        tools.GammaLineCounting(args["data"][0],"data",args["output"][0])

    if args["sim"]:
        tools.GammaLineCounting(args["sim"][0],"sim",args["output"][0])
                           
    if args["fccd"]:
        tools.CalculateFCCD(args["fccd"][0], args["output"][0])
    
    if args["av"]:
        tools.CalculateAV(args["av"][0], args["output"][0])
        
    if args["all"]:
        tools.GammaLineCounting(args["all"][0],"data",args["output"][0])
        tools.GammaLineCounting(args["all"][0],"sim",args["output"][0])
        tools.CalculateFCCD(args["all"][0], args["output"][0])
        tools.CalculateAV(args["all"][0], args["output"][0])
      
    if  args["best_fccd"]:
        tools.BestFCCD(args["best_fccd"][0],"gammaline",args["output"][0])

    if args["plot"]:
        tools.BestFCCD(args["plot"][0],"plot", args["output"][0])


if __name__ == "__main__":
   main()
