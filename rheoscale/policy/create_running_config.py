import pandas as pd
from dataclasses import replace
from typing import TypeAlias, Union


from ..config import RheoscaleConfig

Number = Union[int, float]
NumericFieldDict: TypeAlias = dict[str,None | Number]


def create_config_w_data(input_config: RheoscaleConfig, update: NumericFieldDict) -> RheoscaleConfig:
    '''
    this method goes though the input config 
    generates a value if it is not given as an input value
    then returns a new running RheoscaleConfig with all the correct values need for running the analysis
    '''
    return replace(input_config,
                   number_of_positions= update['number_of_positions'],
                   WT_val= update['WT_val'],
                   WT_error= update['WT_error'],
                   min_val= update['min_val'],
                   max_val= update['max_val'],
                   error_val= update['error_val'],
                   number_of_bins= update['number_of_bins'],
                   neutral_binsize= update['neutral_binsize']
                                      )