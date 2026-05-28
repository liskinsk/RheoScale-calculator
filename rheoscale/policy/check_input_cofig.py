from dataclasses import fields
import statistics, math
from typing import Optional
from ..config import RheoscaleConfig
import pandas as pd, numpy as np


def check_and_update_config(input_config: RheoscaleConfig, raw_DMS_data: pd.DataFrame):
    '''
    This is to check that the inputs in the config are possible to fit the data given
    
    This is a fail fast and fail loud method

    it needs to though each input and ensure that they make sense with the input data
     
    Simlar to the validate methods in the config its self except it now in the context of the raw data
    
    '''
    update = make_infer(input_config)

    infer_min_max(raw_DMS_data, input_config.columns, update, input_config.log_scale)
    if input_config.max_val is not None:
        if input_config.log_scale:
            update['max_val'] = np.log10(input_config.max_val)
        else:
            update['max_val']  = input_config.max_val
    if input_config.min_val is not None:
        if input_config.log_scale:
            update['min_val'] = np.log10(input_config.min_val)
        else:

            update['min_val'] = input_config.min_val

    infer_WT(raw_DMS_data, input_config, update)
    
    if input_config.WT_val is not None:
        if input_config.log_scale:
            update['WT_val'] = np.log10(input_config.WT_val)
        else:
            update['WT_val'] = input_config.WT_val
    if input_config.WT_error is not None:
        if input_config.log_scale:
            update['WT_error'] = 0.434 * (input_config.WT_error / input_config.WT_val)
        else:
            update['WT_error'] = input_config.WT_error


    infer_neutral_bin_size(update, input_config.log_scale, input_config.error_val)
    if input_config.neutral_binsize is not None:
        update['neutral_binsize']  = input_config.neutral_binsize

    infer_num_of_pos(raw_DMS_data, input_config, update)
    if input_config.number_of_positions is not None:
        update['number_of_positions']  = input_config.number_of_positions
    
    infer_error_val(raw_DMS_data, input_config, update)
    if input_config.error_val is not None:
        update['error_val']  = input_config.error_val
    
    
    infer_bins(raw_DMS_data, input_config, update)
    if input_config.number_of_bins is not None:
        update['number_of_bins']  = input_config.number_of_bins

    return update

def infer_bins(raw_DMS_data: pd.DataFrame, config:RheoscaleConfig, update:dict):
   
   dict_of_num_pos_tested = raw_DMS_data[config.columns['position']].value_counts().to_dict()
   avg_num_pos_tested = sum(dict_of_num_pos_tested.values())/len(dict_of_num_pos_tested)
   avg_count_WT_as_tested = avg_num_pos_tested+1
   
   median_pos_tested = statistics.median(dict_of_num_pos_tested.values())
   if update['error_val'] == 0.0:
        error_data_based = 20
   else:
        error_data_based = abs(update['min_val']- update['max_val'])/update['error_val']

   possible_bin_cutoffs = [avg_count_WT_as_tested,median_pos_tested,error_data_based ]

   update['number_of_bins'] = math.floor(min(possible_bin_cutoffs))

def infer_error_val(raw_DMS_data: pd.DataFrame, config:RheoscaleConfig,update:dict):
    
        # if type(value) is not float:
        #     raise
    error_list = raw_DMS_data[config.columns['error']].to_list()
    for i in range(len(error_list)):
        if error_list[i] == 0:
            error_list[i] = float(0)
        elif type(error_list[i]) is not float:
            raise TypeError(f"Error value at index {i} must be a float, got {type(error_list[i]).__name__}: {error_list[i]}")
        
    avg_error =  (sum(error_list))/(len(error_list))
    two_sigma_error = avg_error*2
    update['error_val'] = two_sigma_error




def infer_num_of_pos(raw_DMS_data: pd.DataFrame, config:RheoscaleConfig,update:dict):
    num_of_positions = len(raw_DMS_data[config.columns['position']].unique())
    update['number_of_positions'] = num_of_positions

def infer_WT(raw_DMS_data: pd.DataFrame, config: RheoscaleConfig,  update: dict):
    if config.WT_name is None:
        WT_name = 'WT'
    else: WT_name = config.WT_name
    
    if config.WT_val is None:
        just_WT = raw_DMS_data[raw_DMS_data[config.columns['position']] == WT_name]
        num_of_wt = just_WT.shape[0] #gets the number of rows
        if num_of_wt != 0:
            WT_val = (just_WT[config.columns['value']].sum())/num_of_wt
            update['WT_val'] = WT_val
            update['WT_error'] = (just_WT[config.columns['error']].sum())/num_of_wt
            #remove_values from df
            mask = raw_DMS_data[config.columns['position']] == WT_name
            raw_DMS_data.drop(raw_DMS_data[mask].index, inplace=True)
        else:
            raise ValueError(f'WT_value was not added to the config and cannot be found in the Data \nin the DATA WT values must be called "{WT_name}" in the position column')
    
    else:
        WT_val =  config.WT_val
        if WT_name in raw_DMS_data[config.columns['position']].to_list():
           mask = raw_DMS_data[config.columns['position']] == WT_name
           raw_DMS_data.drop(raw_DMS_data[mask].index, inplace=True) 

    #Check that value in 
    
    
    if update['min_val'] <= WT_val <= update['max_val']:
        pass
     
    else:
        if config.WT_val is not None:
            raise ValueError(f'the WT values provided as an input = {WT_val}\n this value is not found with in the bounds of the assay min = {update['min_val']} -- max {update['max_val']}')
        else:
            
            raise ValueError(f'the WT values calculated from the average of {num_of_wt} WT values found in the input file = {WT_val}\n this value is not found with in the bounds of the assay min = {update['min_val']} -- max {update['max_val']}')

def infer_neutral_bin_size(update: dict, is_log, error_val):
    if is_log and error_val is not None:
        WT_error = update['WT_error']
        transformed_error = np.abs((0.434*WT_error)/(10**update['WT_val']))
        update['neutral_binsize'] = transformed_error*4
    
    else:
        WT_error = update['WT_error']
        update['neutral_binsize'] = WT_error*4
    


def make_infer(input_config: RheoscaleConfig) -> dict:
    config_fields = input_config.numeric_or_none_dict()
    return config_fields

def infer_min_max(raw_DMS_data: pd.DataFrame, column_names: dict, update: dict, log_scale:bool):
  
        update['min_val'] = raw_DMS_data[column_names['value']].min()
        update['_true_min'] =raw_DMS_data[column_names['value']].min()
        update['max_val'] = raw_DMS_data[column_names['value']].max()
        update['_true_max'] =raw_DMS_data[column_names['value']].max()
