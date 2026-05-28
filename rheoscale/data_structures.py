import numpy as np
from typing import NamedTuple
from dataclasses import dataclass, fields

class HistogramData(NamedTuple):
    counts: np.ndarray
    bin_edges: np.ndarray
    weights: np.ndarray = None

    def __add__(self, other):
        if not isinstance(other, HistogramData):
            return NotImplemented

        # Bin compatibility check
        if not np.array_equal(self.bin_edges, other.bin_edges):
            raise ValueError("Cannot add histograms with different bin edges")

        # Add counts
        new_counts = self.counts + other.counts

        # Handle weights
        if self.weights is None and other.weights is None:
            new_weights = None
        elif self.weights is None:
            new_weights = other.weights.copy()
        elif other.weights is None:
            new_weights = self.weights.copy()
        else:
            new_weights = self.weights + other.weights

        return HistogramData(
            counts=new_counts,
            bin_edges=self.bin_edges,
            weights=new_weights
        )     
        


class HistogramFactory:
    def __init__(self, bin_edges, weights=None):
        self.bin_edges = bin_edges
        self.weights = weights
    
    def create_hist_data(self, data):
        counts, _ = np.histogram(data, bins=self.bin_edges)
        return HistogramData(counts, self.bin_edges, self.weights)


@dataclass(slots=True)
class RheoScores:
    position: int
    num_of_variants: int
    histogram: HistogramData 
    neutral_score: float=None
    toggle_score: float=None
    rheostat_score: float=None
    weighted_rheostat_score: float=None
    enhancing_score: float=None
    binary: bool=None
    average: float = None
    st_dev: float= None
    assignment: str = None
    flag: bool = None

    def __setattr__(self, name, value):
        if name == 'position' or name == 'num_of_variants' or name == 'histogram' or name =='flag':
            pass
        elif name == 'assignment':
            if value not in ['neutral', 'toggle', 'rheostat', 'enhancing', 'adverse', 'WT/inactive','moderate', None, 'unclassified']:
                raise TypeError(f'{name} must be a name of a position type')
        elif name in {"binary"}:
            if value is not None and not isinstance(value, bool):
                raise TypeError(f"{name} must be bool")
        elif value is not None and not isinstance(value, (int, float, HistogramData)):
            raise TypeError(f"{name} must be numeric")
        
        
        super().__setattr__(name, value)

    