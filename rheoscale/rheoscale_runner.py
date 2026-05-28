from dataclasses import dataclass, replace
from typing import Any, Optional, NamedTuple
import pandas as pd, numpy as np

from .data_structures import HistogramData, HistogramFactory
from .config import RheoscaleConfig
from .errors import RheoscaleError
from .schemas import POSITION_COLUMNS


from .policy.check_input_cofig import check_and_update_config
from .policy.create_running_config import create_config_w_data
from .policy.rheoscores import calculate_rheoscores


from .io import load_data, validate_columns, write_outputs


class RheoscaleRunner:
    """
    High-level orchestrator for RheoScale.

    Responsibilities:
    - validate inputs
    - apply config-driven policy
    - coordinate computation
    - handle errors and IO
    """

    def __init__(self, config: RheoscaleConfig, DMS_data: Optional[pd.DataFrame] = None):
        # starting data
        self.user_config: RheoscaleConfig = config
        self.DMS_data: pd.DataFrame = self._load_data(DMS_data) 
        
        
        # added to while running
        self.running_config: RheoscaleConfig = None
        
        self.histogram_factory: HistogramFactory = None
        self.rheoscale_position_data: pd.DataFrame = None


    def _load_data(self, dms_data):
        
        if self.user_config.input_file_name is not None and dms_data is None:
            data_from_dms = load_data(self.user_config)
        elif self.user_config.input_file_name is None and dms_data is not None:
            data_from_dms =dms_data
        else:
            raise ValueError('Either a input file must be entered within the configeration XOR a df must be added during this call')   

        validate_columns(data_from_dms, self.user_config)
        
        data_from_dms = data_from_dms.dropna(how="all")

        if data_from_dms[self.user_config.columns['value']].isna().any():
            raise ValueError('The input file must have values in the functions colums')
        
        if self.user_config.log_scale:
            data_from_dms = self.transform_raw_data(data_from_dms)

        if self.user_config.error_val is not None:
            data_from_dms[self.user_config.columns['error']]  = self.user_config.error_val

        

        return data_from_dms  
    


    def run(self):
        #try:
            
            #checking or creating inputs
            update = self.test_accuracy_of_input_config()
            
            #create a running config
            self.create_analysis_config(update)
            self.histogram_factory = self.calculate_bins_and_weight()
            
            #doing rheoscale analysis
            self.rheoscale_position_data = self.rheoscale_analysis()
            
            #writing outputs
            if self.running_config.output_dir != '':
                self.output_data()

            return self.rheoscale_position_data


    def test_accuracy_of_input_config(self):
        
        return check_and_update_config(self.user_config, self.DMS_data)
        

    def transform_raw_data(self, DMS_data):
        DMS_data[self.user_config.columns['error']] = 0.434*((DMS_data[self.user_config.columns['error']])/(DMS_data[self.user_config.columns['value']]))
                                                        
        DMS_data[self.user_config.columns['value']] = np.log10(DMS_data[self.user_config.columns['value']])

        
        return DMS_data
        

    def create_analysis_config(self, update:dict):
        if self.running_config is not None: 
            raise ValueError('This should only be run once')
        
        self.running_config = replace(self.user_config,
                   number_of_positions= update['number_of_positions'],
                   WT_val= update['WT_val'],
                   WT_error= update['WT_error'],
                   min_val= update['min_val'],
                   max_val= update['max_val'],
                   error_val= update['error_val'],
                   number_of_bins= update['number_of_bins'],
                   neutral_binsize= update['neutral_binsize'],
                   _true_min  = update['_true_min'],
                   _true_max = update['_true_max']
        )
        self.running_config._validate_and_make_output()
 
    def calculate_bins_and_weight(self):
        bins_size = np.abs(self.running_config.min_val-self.running_config.max_val)/self.running_config.number_of_bins
        if self.running_config.neutral_binsize == 0.0 or self.running_config.neutral_binsize == 0:
            self.running_config = replace(self.running_config, neutral_binsize=bins_size*2)

        #add one for the final edge
        if self.running_config.dead_extremum == 'Min':
            if self.running_config._true_min > self.running_config.min_val:
                dead = np.linspace(start=self.running_config.min_val, stop=bins_size+self.running_config.min_val, num=1)
            else:
                dead = np.linspace(start=self.running_config._true_min, stop=bins_size+self.running_config.min_val, num=1)
            remaining= np.linspace(start=bins_size+self.running_config.min_val, stop=self.running_config.max_val, num=self.running_config.number_of_bins)
            dead_value = dead[0]
            #if dead_value != self.running_config.min_val:
            bin_edges = np.concatenate([dead, remaining])
            
        else: 
             dead = np.linspace(start=bins_size+self.running_config.max_val, stop=self.running_config._true_max, num=1)
             remaining= np.linspace(start=bins_size+self.running_config.min_val, stop=self.running_config.max_val, num=self.running_config.number_of_bins-1)
             bottom = np.linspace(start=remaining[0]-bins_size, stop=remaining[0], num=1)
             
             if dead != self.running_config.max_val:
                bin_edges = np.concatenate([bottom, remaining, dead])
             else: 
                bin_edges = np.concatenate([bottom, remaining])
        wt_bin = np.digitize(self.running_config.WT_val, bin_edges)-1

        if self.running_config.dead_extremum == "Min":
            
            dead_bin = 0
        else:
            
            dead_bin= np.array(range(self.running_config.number_of_bins)).max()
        
        weights = np.full(self.running_config.number_of_bins, 3)
        
        for weight in range(self.running_config.number_of_bins):
            if abs(weight - dead_bin) == 1:
                weights[weight] = 2
            if abs(weight - wt_bin) == 1:
                weights[weight] = 2
            if weight == dead_bin:
                weights[weight] = 1
            if weight == wt_bin:
                weights[weight] = 1

        hist_fact = HistogramFactory(bin_edges=bin_edges, weights=weights)

        return  hist_fact

    def rheoscale_analysis(self):
        
        
        return calculate_rheoscores(self.running_config, self.DMS_data, self.histogram_factory)
        

    def output_data(self):
        write_outputs(self.running_config, self.rheoscale_position_data)

        
