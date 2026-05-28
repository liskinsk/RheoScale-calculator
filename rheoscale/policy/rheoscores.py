import pandas as pd, numpy as np
from dataclasses import asdict
from  ..data_structures import RheoScores
from ..config import RheoscaleConfig
from ..schemas import POSITION_COLUMNS
from..analysis.rheo_scoring import compute_all_rheo_scores #this only takes one position at a time


def calculate_rheoscores(config: RheoscaleConfig, DMS_values: pd.DataFrame, histogram_factory) -> pd.DataFrame:
    position_rows = []
    count = 0
    for position in DMS_values[config.columns['position']].unique():
        DMS_values_one_pos = DMS_values[DMS_values[config.columns['position']] == position]
        if position == '119':
            pass
        rheo_scores: RheoScores = compute_all_rheo_scores(position, config, DMS_values_one_pos, histogram_factory)
        
        rheo_scores.assignment, rheo_scores.flag = make_assignment(config, rheo_scores)
        position_rows.append(asdict(rheo_scores))

    position_df = pd.DataFrame(position_rows, columns=POSITION_COLUMNS)


    return position_df
    
    return position_df #this needs to have position columns

def make_assignment(running_config: RheoscaleConfig, rheo_scores: RheoScores) ->str:
    if rheo_scores.assignment is not None:
        raise ValueError('assignment already made')
    
    #only things with mre than 6 muts get a score
    if rheo_scores.num_of_variants <5:
        return None, False

    if rheo_scores.neutral_score >= running_config.neutral_threshold:
        return 'neutral', False

    if rheo_scores.toggle_score >= running_config.toggle_threshold:
        return 'toggle', False

    if rheo_scores.enhancing_score >= running_config.enhancing_threshold:
        return 'enhancing', False
    
    if rheo_scores.weighted_rheostat_score >= running_config.rheostat_threshold:
        return 'rheostat', False
    
    wild_type = running_config.WT_val

    if running_config.dead_extremum == 'Min':
        dead = running_config.min_val
    elif running_config.dead_extremum == 'Max':
        dead = running_config.max_val

    wt_inactive_check =rheo_scores.toggle_score >0.35 and  rheo_scores.weighted_rheostat_score <0.35 and rheo_scores.neutral_score>0.35

    if 0.4 < rheo_scores.toggle_score < running_config.toggle_threshold and 0.4 < rheo_scores.neutral_score < running_config.neutral_threshold and  rheo_scores.weighted_rheostat_score < 0.3:
        if wt_inactive_check:
            return "WT/inactive", True 
        else:
            return "WT/inactive", False

    if running_config.error_val != 0:
        percent_allowed = running_config.error_val
    else:
        percent_allowed = np.abs(running_config.min_val-running_config.max_val)/running_config.number_of_bins



    compare_value = np.abs(np.abs(rheo_scores.average - dead) - np.abs(rheo_scores.average - wild_type))
    
    if np.abs(rheo_scores.average - dead) < np.abs(rheo_scores.average - wild_type):
        if wt_inactive_check or (compare_value<percent_allowed):
            return'adverse', True
        
        else:
            return 'adverse', False
        
    if np.abs(rheo_scores.average - dead) > np.abs(rheo_scores.average - wild_type): 
        if wt_inactive_check or (compare_value<percent_allowed):
            return'moderate', True
        else:
            return 'moderate', False
    
    


    
    
    return 'unclassified', False
