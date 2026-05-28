import pandas as pd, numpy as np
from dataclasses import dataclass, fields
from ..config import RheoscaleConfig
from typing import Union
from ..data_structures import HistogramData, HistogramFactory, RheoScores

def compute_all_rheo_scores(position:str,runnig_config: RheoscaleConfig, DMS_position_data: pd.DataFrame, hist_info: HistogramFactory) -> RheoScores:
    
    '''
    These funtions take one position and calculate the scores for every position'''
    
    hist_info = hist_info.create_hist_data(DMS_position_data[runnig_config.columns['value']])

    num_of_var = DMS_position_data.shape[0] +1

    pos_rheoscores = RheoScores(position=position, num_of_variants=num_of_var, histogram=hist_info)

    if position == '4':
        pass
    pos_rheoscores.neutral_score = compute_neutral_score(runnig_config, DMS_position_data)

    pos_rheoscores.enhancing_score = compute_enhancing_score(runnig_config, DMS_position_data)
    
    pos_rheoscores.toggle_score = compute_toggle_score(runnig_config, DMS_position_data, hist_info)

    pos_rheoscores.rheostat_score = compute_rheostat_score(runnig_config, DMS_position_data, hist_info)
    
    pos_rheoscores.weighted_rheostat_score = coumpute_weighted_rheostat_score(runnig_config, DMS_position_data, hist_info)

    pos_rheoscores.binary = is_binary(runnig_config, DMS_position_data, hist_info)
    
    pos_rheoscores.average = np.average(DMS_position_data[runnig_config.columns['value']].values)

    pos_rheoscores.st_dev = np.std(DMS_position_data[runnig_config.columns['value']].values)

    return pos_rheoscores    


def compute_neutral_score(config: RheoscaleConfig, DMS_position_data: pd.DataFrame) -> np.float64:
    num_of_variants = DMS_position_data.shape[0]
    num_of_WT_like = count_overlap_with_WT(DMS_position_data, config.columns['value'], config.columns['error'], config.WT_val, (config.neutral_binsize/2))
    score = num_of_WT_like/num_of_variants
    return np.float64(score)


def count_overlap_with_WT(df, value_column, error_column, WT_val, WT_error):
    #create an array of combined error for each measurement

    #find the abs difference of the each measurement - WT_val if this value less then the combination of the errors then they overlap 
    return (np.abs(df[value_column] - WT_val) <  WT_error).sum()
    

def compute_enhancing_score(config: RheoscaleConfig, DMS_position_data: pd.DataFrame) -> np.float64:
    num_of_variants = DMS_position_data.shape[0]
    
    num_of_WT_like = enhanced_than_WT(DMS_position_data, config.columns['value'], config.columns['error'], config.WT_val, (config.neutral_binsize/2), dead_outcome= config.dead_extremum)

    score = num_of_WT_like/num_of_variants
    return np.float64(score)
    

def enhanced_than_WT(df, value_column, error_column, WT_val, WT_error, n_sigma=1, dead_outcome= 'Min'):
    if dead_outcome == 'Min':
        #find the  difference of the each measurement - WT_val if this value greater then the combination of the errors then measuement is enhanced
        return ((df[value_column] - (WT_val+WT_error)) > 0).sum()
    elif dead_outcome == 'Max':
        #find the  difference WT_val - error  if this value greater then the combination of the errors then measuement is enhanced
        return (((WT_val-WT_error) -df[value_column]) > 0).sum()  
    else:
        raise ValueError('the dead_extemum value must be "Min" or "Max"')


def compute_toggle_score(runnig_config: RheoscaleConfig, DMS_position_data: pd.DataFrame, hist= HistogramData):
    num_of_variants = DMS_position_data.shape[0]
    bin_counts = hist.counts 
    if runnig_config.dead_extremum == 'Min':
        dead_index = 0
    elif runnig_config.dead_extremum == 'Max':
        dead_index = -1
    else:
        raise ValueError('dead_extremum must be "Min" or "Max"')
    dead_outcomes = bin_counts[dead_index]
    return dead_outcomes/num_of_variants

def compute_rheostat_score(runnig_config: RheoscaleConfig, DMS_position_data: pd.DataFrame, hist: HistogramData):
    
    num_of_bins_filled = (hist.counts > 0).sum()
   
    #find if WT bin is filled
    
    bin_idx = np.searchsorted(hist.bin_edges, runnig_config.WT_val, side="right") -1
    if bin_idx == (9):
        pass
    elif bin_idx == (8):
        pass
    if 0 <= bin_idx < len(hist.counts):
        if hist.counts[bin_idx] > 0:
            #WT bin has been filled 
            pass   
        else:
            #WT bin is not filled
            num_of_bins_filled +=1 
    else:
        raise ValueError(f"x={hist.bin_edges,runnig_config.WT_val} is outside the histogram range")
    
    rheo_score = num_of_bins_filled / runnig_config.number_of_bins

    return rheo_score
    

def coumpute_weighted_rheostat_score(runnig_config: RheoscaleConfig, DMS_position_data: pd.DataFrame, hist: HistogramData):

    
    num_of_bins_filled_weighted = ((hist.counts > 0) * hist.weights).sum()
   
    #find if WT bin is filled
    bin_idx = np.searchsorted(hist.bin_edges,runnig_config.WT_val, side="right") -1
    if 0 <= bin_idx < len(hist.counts):
        if hist.counts[bin_idx] > 0:
            #WT bin has been filled 
            pass   
        else:
            #WT bin is not filled ! FOR wieghted it is still only adding one
            num_of_bins_filled_weighted +=1 
            
    else:
        raise ValueError(f"x={hist.bin_edges,runnig_config.WT_val} is outside the histogram range")
    
    sum_w = hist.weights.sum()

    weighted_rheo_score = num_of_bins_filled_weighted / sum_w

    return weighted_rheo_score


def is_binary(runnig_config: RheoscaleConfig, DMS_position_data: pd.DataFrame, hist: HistogramData):
    
    num_of_bins_filled = (hist.counts > 0).sum()
   
    #find if WT bin is filled
    bin_idx = np.searchsorted(hist.bin_edges,runnig_config.WT_val, side="right") - 1
    if 0 <= bin_idx < len(hist.counts):
        if hist.counts[bin_idx] > 0:
            #WT bin has been filled 
            pass   
        else:
            #WT bin is not filled ! FOR wieghted it is still only adding one
            num_of_bins_filled +=1 
            
    else:
        raise ValueError(f"x={hist.bin_edges,runnig_config.WT_val} is outside the histogram range")
    
    return True if num_of_bins_filled == 2 else False
       