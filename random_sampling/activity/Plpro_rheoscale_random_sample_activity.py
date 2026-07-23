from rheoscale.rheoscale_runner import RheoscaleRunner
from rheoscale.config import RheoscaleConfig
import pandas as pd, numpy as np
from Bio.Data import CodonTable
from collections import defaultdict

def sample_substitutions(
    df: pd.DataFrame,
    position_col: str = "Position",
    min_subs: int = 6,
    max_subs: int = 7,
    exclude_positions: list = None,
    random_state: int = None,
) -> pd.DataFrame:
    """
    Randomly downsample substitutions at each position to between min_subs and max_subs.
 
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing FRET data.
    position_col : str
        Column name identifying the position (default: "Position").
    min_subs : int
        Minimum number of substitutions to keep per position (default: 6).
    max_subs : int
        Maximum number of substitutions to keep per position (default: 7).
    exclude_positions : list, optional
        Positions to keep entirely (e.g. ["WT"]). Rows at these positions
        are passed through untouched.
    random_state : int, optional
        Seed for reproducibility. If None, sampling is different each run.
 
    Returns
    -------
    pd.DataFrame
        Downsampled dataframe, preserving original row order within each position.
    """
    if exclude_positions is None:
        exclude_positions = ["WT"]
 
    rng = np.random.default_rng(random_state)
 
    mask_exclude = df[position_col].isin(exclude_positions)
    df_excluded = df[mask_exclude]
    df_to_sample = df[~mask_exclude]
 
    sampled_parts = []
    for position, group in df_to_sample.groupby(position_col, sort=False):
        n_available = len(group)
        # Randomly pick 6 or 7 for this position
        n_keep = int(rng.integers(min_subs, max_subs + 1))
        n_keep = min(n_keep, n_available)  # never ask for more than we have
        sampled_parts.append(group.sample(n=n_keep, random_state=int(rng.integers(1e6))))
 
    df_sampled = pd.concat(sampled_parts)
 
    result = pd.concat([df_excluded, df_sampled]).sort_index().reset_index(drop=True)
 
    return result


def count_assignment_transitions(df1, df2,
                                 position_col="position",
                                 assignment_col="assignment"):
    # Keep only the needed columns
    left = df1[[position_col, assignment_col]].copy()
    right = df2[[position_col, assignment_col]].copy()

    left[position_col] = pd.to_numeric(left[position_col], errors="coerce")
    right[position_col] = pd.to_numeric(right[position_col], errors="coerce")

    left = left.dropna(subset=[position_col])
    right = right.dropna(subset=[position_col])

    left[position_col] = left[position_col].astype(int)
    right[position_col] = right[position_col].astype(int)
    # Merge on position
    merged = left.merge(
        right,
        on=position_col,
        suffixes=("_1", "_2")
    )

    counts = defaultdict(int)
    positions = defaultdict(list)

    for _, row in merged.iterrows():

        pos = int(row[position_col])

        c1 = str(row[f"{assignment_col}_1"]).strip()
        c2 = str(row[f"{assignment_col}_2"]).strip()

        if c1.lower() == c2.lower():
            transition = "Maintained"
        else:
            transition = f"{c1}->{c2}"

        counts[transition] += 1
        positions[transition].append(pos)

    results = []

    for transition in counts:
        results.append({
            "Transition": transition,
            "Count": counts[transition],
            "Positions": positions[transition]
        })

    return pd.DataFrame(results)


def make_percent(counts_dict, decimals=2):   

    total = sum(counts_dict.values())

    if total == 0:
        return {key: 0 for key in counts_dict}

    return {
        key: round(value / total, decimals)
        for key, value in counts_dict.items()
    }

def append_dict_to_df(data_dict, df=None):
    """
    Append a dictionary as a row to a DataFrame.
    
    Parameters
    ----------
    data_dict : dict
        Dictionary to append as a row.
    df : pd.DataFrame or None
        Existing DataFrame. If None, a new one is created.

    Returns
    -------
    pd.DataFrame
    """
    new_row = pd.DataFrame([data_dict])

    if df is None:
        return new_row

    return pd.concat([df, new_row], ignore_index=True).fillna(0)

def filter_repeated_values(df, column, min_count=19):
    counts = df[column].value_counts()
    
    # Keep values that meet the threshold OR are WT
    keep_values = counts[(counts >= min_count) | (counts.index == "WT")].index
    
    filtered_df = df[df[column].isin(keep_values)]
    
    print(f"Unique values remaining: {filtered_df[column].nunique()}")
    print(f"Rows remaining: {len(filtered_df)}")
    
    print("Counts of remaining values:")
    print(filtered_df[column].value_counts())
    
    return filtered_df

def main():
    plpro_data = pd.read_csv(r"C:\Lab_code\Rheoscale\random_sample_tests\activity\Plpro activity raw data.csv")


        
    plpro_data = plpro_data[plpro_data['Substitution'] != '*']
    plpro_data = plpro_data.dropna()
    df = None

    plpro_data = filter_repeated_values(plpro_data,'Position')
    ref_rheo_class =pd.read_csv(r"C:\Lab_code\Rheoscale\random_sample_tests\activity\Plpro_activity_Feb_26_classifications.csv")
    
    # values = plpro_data['Position'].astype(str)
    #plpro_data.to_csv('cut.csv')
    # ref_rheo_class = ref_rheo_class[ref_rheo_class['position'].astype(str).isin(values)].copy()
    for value in range(5,21):
        for rep in range(10):
            
            plpro_data_small = sample_substitutions(plpro_data, min_subs=value, max_subs=value)
            #plpro_data_small.to_csv('small_cut.csv')

            normalized_config = RheoscaleConfig(f'Random_{value}_rep_{rep}',min_val=0.23, max_val=1.28674797)#, number_of_bins=10)

            normalized_runner = RheoscaleRunner(normalized_config, plpro_data_small)
            rheostats =normalized_runner.run() #returns a Dataframe with positions calucations and numbers
            rheo_class = rheostats[['position', 'assignment']]
            
            if value ==20:
                pass

            result =count_assignment_transitions(ref_rheo_class, rheo_class)

            result_dict = dict(zip(result["Transition"], result["Count"]))
            
            result_dict = make_percent(result_dict)
            result_dict['Number of Substitutions'] = value
            result_dict['Rep'] = rep
            df =append_dict_to_df(result_dict, df)
    
    assignment_stats = (
    df.groupby("Number of Substitutions", as_index=False)
      .agg(["mean", "std"])
)
    
    assignment_stats.to_csv("let rheoscale choose 10 reps avg activity cut july 15 2026.csv")



main()