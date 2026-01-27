import matplotlib.pyplot as plt
from pathlib import Path
import json


def PlotFCCDs(ID, OutPath):
    json_file = Path(f"{OutPath}/Combined/combined_detectors.json")
    with json_file.open("r", encoding="utf-8") as f:
        data = json.load(f)

    detectors = list(data.keys())

    marker_map = {
     'am_HS1':"o", 
     'ba_HS4':"s"}
    
    prefixes = sorted({d[:3] for d in detectors})
    cmap = plt.get_cmap("tab20")
    color_map = {prefix: cmap(i % 20) for i, prefix in enumerate(prefixes)}


    fig, ax = plt.subplots(figsize=(12, 6))

    for i, detector in enumerate(detectors):
        prefix = detector[:3]
        color = color_map[prefix]

        if "am_HS1" in data[detector] and "FCCD" in data[detector]["am_HS1"]:
            y=data[detector]["am_HS1"]["FCCD"]["Central"]
            yerr_pos=data[detector]["am_HS1"]["FCCD"]["ErrPos"]
            yerr_neg=data[detector]["am_HS1"]["FCCD"]["ErrNeg"]
            yerr = [[yerr_neg], [yerr_pos]]         
            ax.errorbar(i, y, yerr=yerr, marker=marker_map["am_HS1"],color=color, capsize=3)
        if "ba_HS4" in data[detector] and "FCCD" in data[detector]["ba_HS4"]:
            y=data[detector]["ba_HS4"]["FCCD"]["Central"]
            yerr_pos=data[detector]["ba_HS4"]["FCCD"]["ErrPos"]
            yerr_neg=data[detector]["ba_HS4"]["FCCD"]["ErrNeg"]
            yerr = [[yerr_neg], [yerr_pos]]         
            ax.errorbar(i, y, yerr=yerr, marker=marker_map["ba_HS4"],color=color, capsize=3)
        
    ax.set_xticks(range(len(detectors)))
    ax.set_xticklabels(detectors, rotation=45)
    ax.grid(linestyle="dashed", linewidth=0.5)
    ax.set_ylabel("FCCD [mm]", fontsize=20)
    ax.set_ylim([0,1.6])

    handles_marker = [plt.Line2D([0], [0], marker=m, color="k", linestyle="", label=s) for s, m in marker_map.items()]
    handles_color = [plt.Line2D([0], [0], marker="o", color=c, linestyle="", label=p) for p, c in color_map.items()]
    
    first_legend = ax.legend(handles=handles_marker, title="Source", loc="upper left")
    ax.add_artist(first_legend)
    ax.legend(handles=handles_color, title="Detector Order", loc="upper right")
    
    plt.tight_layout()
    
    plt.savefig(f"{OutPath}/Combined/FCCDs-{ID}.png")
    print(f"Combined detector plot created: {OutPath}/Combined/FCCDs-{ID}.png")
