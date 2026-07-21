import pandas as pd
import numpy as np
from pathlib import Path

from .data_structures import HistogramData
import matplotlib
import sys
if "debugpy" in sys.modules:
    matplotlib.use("Agg")
# Key is 'W' because plot_all_positions indexes colors via pos_type[0].upper(),
# which takes only the first character of the assignment string. 'WT/inactive'[0] = 'W'. - HC
colors = {'E': '#F94040', 'N':'#00C000' , "A":'#DE8BF9', 'T': '#AD07E3', 'R':'#FFDC00', 'M':'#A0FFA0', 'W': '#00FFFF'}

import matplotlib.pyplot as plt

def plot_all_positions(positions: list, hist_list: list, classifcation_dict:dict, dead_extremum, true_max, WT_value,dead_value,far_value, save_path, neutral_bin_size, prefix='', all_pos=False, is_even_bins=False):
    master_hist_data = hist_list[0]
    for i in range(len(positions)):
        name = prefix+'_pos_'+str(positions[i])
        pos_type = classifcation_dict[positions[i]]
        # Skip plotting for positions with no classification (n < 5 variants). pos_type
        # is None, which pandas converts to NaN when building the classification dict.
        # Attempting pos_type[0].upper() on a float causes a TypeError. - HC
        if all_pos:
            if not isinstance(pos_type, str):
                print(f"  Skipping plot for position {positions[i]} — no classification (too few variants)")
            else:
                make_tuning_plot_one_pos(hist_list[i], dead_extremum, WT_value,dead_value,far_value, neutral_bin_size, save_path, tle=name, true_max=true_max, position= positions[i], pos_type=pos_type[0].upper(), is_even_bins=is_even_bins)
        if i!=0:
             master_hist_data+= hist_list[i]
    all_title = prefix+'_all_positions'
    make_tuning_plot_one_pos(master_hist_data, dead_extremum, WT_value,dead_value,far_value, neutral_bin_size, save_path, tle=all_title, true_max=true_max, is_all=True, is_even_bins=is_even_bins)

    
def make_tuning_plot_one_pos(hist_data: HistogramData , dead_extremum, WT_value, dead_value, far_value, neutral_bin_size,path,     
    true_max: float = 0.0,
    position: str = None, 
    pos_type: str = None,
    tle: str = "Histogram",
    is_all: bool = False,
    is_even_bins: bool = False,
    xlabel: str = "bins",
    ylabel: str = "number of variants",
    log_y: bool = True,
    label_precision: int = 2):
        
        counts = hist_data.counts
        bin_edges = hist_data.bin_edges.copy()

        # ---- sanity check ----
        if len(bin_edges) != len(counts) + 1:
            raise ValueError("bin_edges must be one element longer than counts")
        # Moves both edges to the min/max override for plotting purposes. is_even_bins should apply to
        # both ends of the histogram, not just the dead side. - HC
        if is_even_bins:
             if dead_extremum == 'Min':
                  bin_edges[0] = dead_value
                  bin_edges[-1] = far_value
             else:
                  bin_edges[-1] = dead_value
                  bin_edges[0] = far_value
        
        # bin_widths = np.diff(bin_edges)
        # bin_centers = bin_edges[:-1] + bin_widths / 2

        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bin_widths = np.diff(bin_edges)

        
        fmt = f"{{:.{label_precision}f}}"
        bin_labels = [
            f"{fmt.format(left)}–{fmt.format(right)}"
            for left, right in zip(bin_edges[:-1], bin_edges[1:])
        ]

        # ---- create figure & axes ----
        fig, ax = plt.subplots()
        if position is not None:
            pos_color = colors[pos_type]
        else:
             pos_color= 'black'
        # ---- plot ----
        ax.bar(
            bin_centers,
            counts,
            width=bin_widths*0.85,
            align="center",
            color=pos_color,
            edgecolor='black'
        )
        
        

        WT_index =  np.digitize(WT_value, bin_edges) - 1
        WT_bin_center = bin_centers[WT_index]

        ax.axvline(dead_value, color='red')
        
        WT_bin_bottom = WT_value-(neutral_bin_size/2)
        WT_bin_top = WT_value+(neutral_bin_size/2)
        if WT_bin_top > true_max:
             WT_bin_top = true_max+(0.01*WT_bin_top)
        ax.axvspan(WT_bin_bottom, WT_bin_top, alpha=.3, color='green')

        # ---- axes formatting ----
        if not is_all:
            ax.set_ylim(0.001, 21)
            ax.set_yticks([i for i in range(5,21, 5)])
        else:
            
            ax.set_yscale('log')
            ax.set_ylim(1)
        ax.set_xticks(bin_centers)
        ax.set_xticklabels(bin_labels, rotation=45, ha="right")

        ax.set_xlabel(xlabel)
        ax.set_title(tle)

        if log_y:
            pass
            #ax.set_yscale("log")

        fig.tight_layout()
        plt.savefig(Path(path) / f'{tle}.png', dpi=600)
        
        plt.close()

        