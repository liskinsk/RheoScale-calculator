# io.py
import pandas as pd, os
from types import MappingProxyType
from typing import Optional
from pathlib import Path
from .config import RheoscaleConfig
from .plotting_histograms import plot_all_positions

def load_data(config: RheoscaleConfig) -> pd.DataFrame:

    # Otherwise load from file
    if config.input_file_name is not None:
        DMS_data = pd.read_csv(config.input_file_name)
        
        return DMS_data
    else:
        raise ValueError('Either a input file must be entered within the configeration XOR a df must be added during this call')
    
    
def validate_columns(df: pd.DataFrame, config: RheoscaleConfig):
    missing = set(config.columns.values()) - set(df.columns)

    if not missing:
        print("Column names in config match the CSV or DataFrame")
    else:
        raise ValueError(
            f"Missing columns in DataFrame: {missing}"
        )
    
def write_outputs(running_config: RheoscaleConfig, position_df: pd.DataFrame):
    out_dir = (Path(running_config.output_dir) / f'{running_config.protein_name}_Rheoscale')
    
    json_path = out_dir / f'{running_config.protein_name}_running_config.json'
    #make_config_output 
    
    running_config.to_json(json_path)

    #make data sheet output 
    data_path = out_dir / f'{running_config.protein_name}_raw_data.csv'
    position_df.to_csv(data_path)

    #just classifcations
    just_pos_and_assign = position_df[['position', 'assignment']]
    class_path = out_dir / f'{running_config.protein_name}_classifications.csv'
    just_pos_and_assign.to_csv(class_path)
    classifcation_dict =just_pos_and_assign.set_index("position")["assignment"].to_dict()
    if running_config.output_histogram_plots:
        is_hist_plots = True
    else:
        is_hist_plots= False
    if running_config.dead_extremum == "Min":
        dead_value = running_config.min_val
    else:
        dead_value = running_config.max_val

    position_list = position_df['position'].to_list()
    hist_list = position_df['histogram'].to_list()
    plot_output = out_dir / f'{running_config.protein_name}_plots'
    os.makedirs(plot_output, exist_ok=True)
    


    plot_all_positions(position_list, 
                       hist_list, 
                       classifcation_dict,
                       running_config.dead_extremum,
                       running_config._true_max,
                        running_config.WT_val,dead_value ,
                        plot_output,
                        running_config.neutral_binsize,running_config.protein_name, 
                       is_hist_plots, 
                       is_even_bins=running_config.even_bins)




