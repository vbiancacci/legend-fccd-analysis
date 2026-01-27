import sys
import os
import argparse

import modules.CompareFCCDs.CompareFCCDs_createjson as createjson
import modules.CompareFCCDs.CompareFCCDs_createplot as createplot

# Master script to automatically launch other scripts in order to automatically process multiple detectors at once

CodePath=os.path.dirname(os.path.realpath(__file__))

if len(sys.argv) != 3 and sys.argv[-1]!="-h" and sys.argv[-1]!="--help":
    print('')
    print('#################################################################################')
    print('## Usage: python3 Compare_Detectors.py [OPERATION MODE] [CONFIG FILE] [OPTIONS]             ##')
    print('## Plese enter "python3 Compare_Detectors.py --help" for more information                   ##')
    print('#################################################################################')
    print('')


def main():

    par = argparse.ArgumentParser(description="Compare FCCDs among detectors",
                                  epilog="Report bugs to <valentina.biancacci@gssi.it>",
                                  usage="python3 make.py [OPERATION: -j, -p, -a] "
                                  "[OPTIONS: -t tl_model (default 'no_tl') -f frac_FCCDbore (default 0.5)"
                                  "-e energy filter (default cuspEmax_ctc) -c cuts (default True) -o path]"
    )
    arg = par.add_argument
    arg("-t", "--tl_model", nargs=1, help="tl model", default="notl", metavar="")
    arg("-f", "--frac_FCCDbore", nargs=1, help="fraction of FCCD at the borehole", default=0.5, metavar="")
    arg("-e", "--energy_filter", nargs=1, help="energy filter", default="cuspEmax_ctc_cal", metavar="")
    arg("-c", "--cuts", nargs=1, help="apply cuts", default=True, metavar="")
    arg("-o", "--output", nargs=1, help="output path", default=f"{CodePath}/results", metavar="")
    
    args=vars(par.parse_args())

    #if args["json"]:
    if args["cuts"] == "False":
        cuts = ""
    else:
        cuts = "-cuts"

    ID=f"{args['tl_model']}-fracFCCDbore{args['frac_FCCDbore']}-{args['energy_filter']}{cuts}"
    createjson.MergeFCCDs(ID,args["output"])
    createplot.PlotFCCDs(ID,args["output"])

    print("Done!")
    print("Results stored in "+args["output"]+"/Combined/")

if __name__ == "__main__":
   main()
