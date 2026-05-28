from dataclasses import dataclass
import pandas as pd 

@dataclass
class DMS_averages:
    data_frame = pd.DataFrame(columns=['positions', 'avg', 'avg_log', 'std_dev', 'std_dev_log'])
    

