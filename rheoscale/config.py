from dataclasses import dataclass, field, replace, fields, asdict
from typing import Optional, Dict, Literal, Union, ClassVar,  get_origin, get_args, TypeAlias
from types import MappingProxyType
import os, pandas as pd
from pathlib import Path
#interal keys
Number = Union[int, float]
NumericFieldDict: TypeAlias = dict[str,None | Number]



ColumnKey = Literal[
    "position",
    "substitution",
    "value",
    "error"

]

class FixedKeysDict(dict):
    def __setitem__(self, key, value):
        if key not in self:
            raise KeyError(f"Key '{key}' cannot be added")
        super().__setitem__(key, value)

@dataclass (frozen=True)
class RheoscaleConfig:
    protein_name: str = field(
        metadata={"label": "Protein name"}
    )
    input_file_name: str = field(default=None, metadata={'label': 'Input File Path'})
    number_of_positions: int= field(default=None, metadata={'label': 'Number of positions'})
    log_scale: bool = field(default=False, metadata={'label': 'Convert values to log scale'})
    WT_val: float = field(default=None, metadata={'label': 'Measured value of WT'})
    WT_error: float= field(default=None, metadata={'label': 'Measured value of WT error'}) 
    WT_name: str = field(default=None, metadata={'label': 'Name of WT in position coloumn in data'})
    
    min_val: float = field(default=None, metadata={'label': 'Minimum measurable value'})
    max_val: float = field(default=None, metadata={'label': 'Maximum measurable value'}) 
    error_val: float = field(default=None, metadata={'label': 'Error override'}) 
    number_of_bins: int = field(
        default=None,
        metadata={
            "label": "Number of bins",
            "min": 2,
            "max": 20
        })
    dead_extremum: Literal["Min", "Max"]= field(default="Min", metadata={'label': 'Dead Extremum'})
    neutral_binsize: float = field(default=None, metadata={'label': 'Size of the neutral bin'}) 
    output_dir: str =field(default='Rheoscale_analysis', metadata={'label': 'Name of output dir that is added after the protein name'}) 
    
    output_histogram_plots: bool = field(default=False, metadata={'label': 'Output all the histogram plots of each position'}) 
    even_bins: bool = field(default=True, metadata={'label': 'Even bins (Even bins if True, if false dead bin will be as wide as measured values)'}) 
    _true_min: float= None
    _true_max: float= None
    '''
    potential bug is having WT val and error val at the same time'''

    # Threshold defaults
    enhancing_threshold: ClassVar[float]= 0.8
    neutral_threshold: ClassVar[float]= 0.7
    rheostat_threshold: ClassVar[float]= 0.5
    toggle_threshold: ClassVar[float]= 0.64
    
    columns: Dict[ColumnKey, str] = field(
        default_factory=lambda: {
            "position": "Position",
            "substitution": 'Substitution',
            "value": 'Value',
            "error": 'Error'
        }, metadata= {'label': 'Name of Columns'})
    



    def __post_init__(self):
        self._validate_name()
        self._validate_input_file()
        self._validate_thresholds()
        self._validate_bins()
        self._validate_num_pos()
        
        self._validate_WT()
    def _validate_name(self):
        if not isinstance(self.protein_name, str):
            raise ValueError('The input name MUST be a string!')

    def _validate_WT(self):
        if self.WT_val is None and self.WT_error is None:
            pass
        elif self.WT_error is not None and self.WT_val is not None:
            pass
        elif self.WT_val is not None and self.WT_name is None:
            pass
        elif self.WT_val is None and self.WT_name is not None:
            pass
        else:
            raise ValueError('if the WT_val AND WT_error must be entered together') 

    def _validate_input_file(self):
        if self.input_file_name is not None:
            is_path = os.path.exists(self.input_file_name)
            if is_path and self.input_file_name.endswith('.csv'):
                pass
            else:
                raise FileNotFoundError(f'the path {self.input_file_name} was not found')
    

    def _validate_bins(self):
        if self.number_of_bins is not None:
            if self.number_of_bins > 20 or self.number_of_bins <=1:
                raise ValueError(f"you entered a number_of_bins = {self.number_of_bins}, you must have at least 2 bins and no more than 20")
    

    def _validate_thresholds(self):
        for name in (
            "enhancing_threshold",
            "neutral_threshold",
            "rheostat_threshold",
            "toggle_threshold",
        ):
            value = getattr(self, name)
            if not (0.0 <= value <= 1.0):
                raise ValueError(
                    f"{name} must be between 0 and 1 inclusive, got {value}"
                )
            
    def _validate_num_pos(self):
        if self.number_of_positions is not None:
            if self.number_of_positions <1:
                raise ValueError('you must analyse at least one position')

    def _validate_and_make_output(self):
        out_dir = Path(self.output_dir)
        if self.output_dir != 'Rheoscale_analysis':
            is_path = os.path.exists(self.output_dir)
            if is_path:
                pass
            else:
                print(f'the path {self.output_dir} was not found creating dir')
                os.makedirs(out_dir)
        else:
            try: 
                print(rf'attempting to create dir: {self.output_dir}\{self.protein_name}_Rheoscale')
                os.makedirs((out_dir / f'{self.protein_name}_Rheoscale'))
                
            except OSError:
                print(rf'file: {self.output_dir}\{self.protein_name}_Rheoscale aready exists please move dir or delete or rerun with new name')

    def change_colums(self, **updates: str) -> "RheoscaleConfig":
        new_cols = dict(self.columns)
        for key, value in updates.items():
            if key not in new_cols:
                raise KeyError(f"Invalid column key: {key}")
            new_cols[key] = value
        return replace(self, columns=MappingProxyType(new_cols))
    
    def numeric_or_none_dict(self) -> NumericFieldDict:
        """
        Returns only fields whose declared type includes int or float,
        with values restricted to int | float | None.
        """
        out: dict[str, None | Number] = {}

        for f in fields(self):
            val = getattr(self, f.name)

            origin = get_origin(f.type)
            args = get_args(f.type)

            if val is None and Number in args:
                out[f.name] = None
                continue

            # Exclude bool explicitly
            if isinstance(val, bool):
                continue

            
            if f.type in (int, float):
                out[f.name] = val
            elif origin is Union and any(t in (int, float) for t in args):
                out[f.name] = val

        return FixedKeysDict(out)


    def validate(self):
        if self.requires_input_file and not self.input_file_name:
            raise ValueError(
                "Config requires an input file, but input_file_name is not set."
            )

        if self.Min_Override and self.Min_Override_Val is None:
            raise ValueError("Min_Override is True but Min_Override_Val is None")

        if self.Max_Override and self.Max_Override_Val is None:
            raise ValueError("Max_Override is True but Max_Override_Val is None")

        if self.Error_Override and self.Error_Override_Val is None:
            raise ValueError("Error_Override is True but Error_Override_Val is None")

    def to_json(self, path: str):
        import json
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=4)

    @classmethod
    def from_json(cls, path: str):
        import json
        with open(path) as f:
            return cls(**json.load(f))