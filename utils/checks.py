import sys
import json
import os
import yaml

class NoValidPosition(Exception):
    pass

class NoValidMeasurement(Exception):
    pass


currentPath=os.path.dirname(os.path.realpath(__file__))

ignored_chars = ' \n\t'


def checkMeasurement(value):
    if not os.path.isfile(value):
        raise NoValidMeasurement(f"The measurement {value} does not exist. Please check the configuration file and metadata.")


def checkInteger(value):
    try:
        value=int(value)
    except ValueError:
        print("The input events number is not a valid integer")
        raise

def checkFloat(value):
    try:
        value=float(value)
    except ValueError:
        print("The input events number is not a valid float")
        raise

def ask_and_print_runs(node):
    ans = input("Do you want to see the full list of available runs, runXXXX: [phi, x, z] ? [y/N]: ").strip().lower()
    if ans in ("y", "yes"):
        for run in node.keys():
            print(f"{run}:  {list(getattr(node, run).source_position.values())}")


def CheckPositionSource(node, user_positions,measurement,detector, campaign):
    phi_position, r_position, z_position = user_positions[0], user_positions[1], user_positions[2]
    phi_pos= list(node.group("source_position.phi_in_deg").keys())
    if phi_position not in phi_pos:
        print(f"Position ERROR: Provided phi position {phi_position} not found in the database for the given measurement {detector}/{campaign}/{measurement}.\n" 
        "Available phi positions are: " + str(phi_pos))
        ask_and_print_runs(node)
        raise NoValidPosition(f"Provided phi position {phi_position} not found in the database for the given measurement {detector}/{campaign}/{measurement}.")
    else:
        data = node.group("source_position.phi_in_deg").get(phi_position)
        r_available=[]
        z_available=[]
        matched_r = []
        for i in data.keys():
            matched_z = False
            try:
                r_pos_=data.get(i).get('source_position').get('r_in_mm')
            except AttributeError:
                r_pos_ = data.get('source_position').get('r_in_mm')
            if r_position!=r_pos_:
                matched_r.append(False)
                if r_pos_ not in r_available:
                    r_available.append(r_pos_)
                continue
            else:
                matched_r.append(True)
                try:
                    z_pos_=data.get(i).get('source_position').get('z_in_mm')
                except:
                    z_pos_ = data.get('source_position').get('z_in_mm')
                if z_position!=z_pos_:
                    matched_z = False
                    if z_pos_ not in z_available:
                        z_available.append(z_pos_)
                    continue
                else: 
                    matched_z = True
                    try:
                        run = data.get(i).get('run')
                    except:
                        run = data.get('run')
                    return run
        if not any(matched_r):
            print(f"Position ERROR: Provided r position {r_position} not found in the database for the given measurement {detector}/{campaign}/{measurement}.\n" 
            "For the provided phi position " + str(phi_position) +"\n"
            "the available r positions are: " + str(r_available))
            ask_and_print_runs(node)
            raise NoValidPosition(f"Provided r position {r_position} not found in the database for the given measurement  {detector}/{campaign}/{measurement}.")
        if matched_z==False:
            print(f"Position ERROR: Provided z position {z_position} not found in the database for the given measurement  {detector}/{campaign}/{measurement}.\n" 
            "For the provided phi position " + str(phi_position) +" and r position " + str(r_position) +"\n"
            "the available z positions are: " + str(z_available))
            ask_and_print_runs(node)
            raise NoValidPosition(f"Provided z position {z_position} not found in the database for the given measurement  {detector}/{campaign}/{measurement}.")






def DefineRunPosition(ConfigNameFile):
    with open(ConfigNameFile) as json_file: 
        config = json.load(json_file)


    detector       = config['Detector']
    campaign       = config['Campaign']
    measurement    = config ['Measurement']
    run            = config['Run']
    phi_position   = config['phiPosition']
    r_position     = config['rPosition']
    z_position     = config['zPosition']


    checkInteger(phi_position)
    checkInteger(r_position)
    checkInteger(z_position)

    MetaDataPath="/global/cfs/cdirs/m2676/users/biancacci/fccd/legend-fccd-analysis/tools/hades-metadata/hardware/configuration"
    MeasurementPath=MetaDataPath+"/"+detector+"/"+campaign+"/"+measurement+".yaml"
    checkMeasurement(MeasurementPath)
    
    if run[0]=="r":
        print(f"Run {run} not found in the metadata.  The run name format is 4-digit string with leading zeros, e.g. '0001' ")
        sys.exit()
    elif run[0].isdigit(): #the user knows the run number
        run="run"+run
        try:
            #node = getattr(node, run)
            with open(MeasurementPath, 'r') as f:
                docs = yaml.full_load(f)
                node = docs.get(run)
                phi_position = node['source_position']['phi_in_deg']
                r_position = node['source_position']['r_in_mm']
                z_position = node['source_position']['z_in_mm']
        except AttributeError:
            print(f"Run {run} not found in the metadata.  The run name format is 4-digit string with leading zeros, e.g. '0001' ")
            sys.exit()
        #phi_position = node.source_position.phi_in_deg
        #r_position = node.source_position.r_in_mm
        #z_position = node.source_position.z_in_mm

         
    else: #the user knows the source position
        user_positions = [phi_position, r_position, z_position]

        run=CheckPositionSource(node, user_positions, measurement, detector, campaign)
    run=run[:1]+run[4:]
    phi_pos=str(phi_position).split('.')[0]
    r_pos=str(r_position).split('.')[0]
    z_pos=str(z_position).split('.')[0]
    source_pos = f"{phi_pos}phi_{r_pos}r_{z_pos}z"

    return run, [phi_pos, r_pos, z_pos], source_pos

