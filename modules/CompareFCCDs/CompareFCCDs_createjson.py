
import os
import json
from pathlib import Path

#Script to compare FCCDs for each detector, for Ba_hs4, Am_HS1 
#Creates a json file of results

CodePath=os.path.dirname(os.path.realpath(__file__))


def MergeFCCDs(ID, OutPath):

    # Base directory
    base_path = Path(f"{OutPath}/AV")

    folders = {
        "am_HS1": base_path / "am_HS1",
        "ba_HS4": base_path / "ba_HS4",
    }

    combined = {}

    # Check that folders exist
    for folder_name, folder_path in folders.items():
        if not folder_path.exists():
            print(f"Folder not found: {folder_path}")

        json_files = list(folder_path.glob(f"AV-*-{ID}.json"))
        if not json_files:
            print(f"No JSON files found in {folder_path}")

        for json_file in json_files:
            with json_file.open("r", encoding="utf-8") as f:
                data = json.load(f)

            if not isinstance(data, dict) or len(data) != 1:
                raise ValueError(
                    f"Invalid format in {json_file}. "
                    "Expected exactly one top-level detector key."
                )

            detector_name, detector_data = next(iter(data.items()))

            if detector_name not in combined:
                combined[detector_name] = {}

            if folder_name in combined[detector_name]:
                raise ValueError(
                    f"Duplicate detector '{detector_name}' in folder '{folder_name}'."
                )

            combined[detector_name][folder_name] = detector_data

    # Write combined file
    if os.path.exists(OutPath+"/Combined") == False:
        os.makedirs(OutPath+"/Combined")
    output_file = Path(f"{OutPath}/Combined/combined_detectors.json")
    with output_file.open("w", encoding="utf-8") as f:
        json.dump(combined, f, indent=4)

    print(f"Combined detector JSON created: {output_file}")

    