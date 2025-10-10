import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def print_usage():
    usage = """
    Usage: python script_name.py <input_file_name> [options]

    This script processes a CSV file containing protein data to classify positions as
    Enhancing, Neutral, Toggle, Rheostat, Moderate, Adverse, WT/inactive, or Unclassified.
    It generates histograms and bar plots for various metrics and writes results to output files.
    
    Required Arguments:
        input_file_name       The name of the input CSV file containing protein data.
    
    Options:
        --log_scale                              Apply log scale conversion (default: True; Input options: True/False)
        --min_override                           Apply minimum override (default: False; Input options: True/False)
        --min_val                                Minimum override value (Requires user value when min_override = True; Input options: any number (float data type))
        --max_override                           Apply maximum override (default: False; Input options: True/False)
        --max_val                                Maximum override value (Requires user value when max_override = True; Input options: any number (float data type))
        --error_override                         Apply error override (default: False; Input options: True/False)
        --error_val                              Error override value (Requires user value when error_override = True; Input options: any number (float data type))
        --bin_override                           Bin number override (Requires user value between 0-20)
        --dead_value                             Defines whether the value corresponding to total loss of signal (“dead”) is “Min” or “Max”(default: "Min")
        --neutral_binsize_override               Override bin size for neutral score calculation (default: (2 x standard bin size); Input options: Any number (float data type))
        --enhancing_score_threshold_override     Enhancing score threshold (default: 0.8; Input options: Any number (float data type) between 0 and 1)
        --neutral_score_threshold_override       Neutral score threshold (default: 0.7; Input options: Any number (float data type) between 0 and 1)
        --rheostat_score_threshold_override      Rheostat score threshold (default: 0.5; Input options: Any number (float data type) between 0 and 1)
        --toggle_score_threshold_override        Toggle score threshold (default: 0.64; Input options: Any number (float data type) between 0 and 1)
        --plot_prefix                            Prefix for plot filenames and folder (default: none)
    
    Example:
        python script_name.py "data.csv" --log_scale True --dead_value "Min" --neutral_binsize_override 0.5 --plot_prefix "LacI"

    Ensure the CSV file has the following columns:
        - Position: The position of the variant.
        - Substitution: The type of substitution.
        - Value: The functional value.
        - Error: The associated error value.

    Output:
        - Histogram image files in .jpg format for each position and overall.
        - A CSV file with calculated metrics and assignments for each position.
        - A pie chart image in .jpg format for the distribution of position assignments.
    """

    print(usage)

# Default values
input_file_name = ""
log_scale = True
Min_Override = False
Min_Override_Val = 0
Max_Override = False
Max_Override_Val = 0
Error_Override = False
Error_Override_Val = 0
bin_override = 0
DeadValue = "Min"
neutral_binsize_override = 0
plot_prefix = ""

# Threshold defaults
enhancing_threshold = 0.8
neutral_threshold = 0.7
rheostat_threshold = 0.5
toggle_threshold = 0.64

# Parsing command-line arguments
if len(sys.argv) < 2:
    print_usage()
    sys.exit(1)

input_file_name = sys.argv[1]
args = sys.argv[2:]

if "--log_scale" in args:
    log_scale = args[args.index("--log_scale") + 1].lower() == "true"
if "--min_override" in args:
    Min_Override = args[args.index("--min_override") + 1].lower() == "true"
if "--min_val" in args:
    Min_Override_Val = float(args[args.index("--min_val") + 1])
if "--max_override" in args:
    Max_Override = args[args.index("--max_override") + 1].lower() == "true"
if "--max_val" in args:
    Max_Override_Val = float(args[args.index("--max_val") + 1])
if "--error_override" in args:
    Error_Override = args[args.index("--error_override") + 1].lower() == "true"
if "--error_val" in args:
    Error_Override_Val = float(args[args.index("--error_val") + 1])
if "--bin_override" in args:
    bin_override = int(args[args.index("--bin_override") + 1])
if "--dead_value" in args:
    DeadValue = args[args.index("--dead_value") + 1]
if "--neutral_binsize_override" in args:
    neutral_binsize_override = float(args[args.index("--neutral_binsize_override") + 1])
if "--enhancing_score_threshold_override" in args:
    enhancing_threshold = float(args[args.index("--enhancing_score_threshold_override") + 1])
if "--neutral_score_threshold_override" in args:
    neutral_threshold = float(args[args.index("--neutral_score_threshold_override") + 1])
if "--rheostat_score_threshold_override" in args:
    rheostat_threshold = float(args[args.index("--rheostat_score_threshold_override") + 1])
if "--toggle_score_threshold_override" in args:
    toggle_threshold = float(args[args.index("--toggle_score_threshold_override") + 1])
if "--plot_prefix" in args:
    plot_prefix = args[args.index("--plot_prefix") + 1]

# Create plots directory
import os
if plot_prefix:
    plots_dir = f"{plot_prefix}_Plots"
else:
    plots_dir = "Plots"
os.makedirs(plots_dir, exist_ok=True)

# Adding User Settings to Output File
add_info = f"Input File Name: {input_file_name}\n"
add_info += f"Log Scale: {log_scale}\n"
add_info += f"Bin Number Override: {bin_override}\n"
add_info += f"Neutral Bin Size Override: {neutral_binsize_override}\n"
add_info += f"Dead Value Location: {DeadValue}\n"
add_info += f"Minimum Override: {Min_Override}\n"
if Min_Override:
    add_info += f"Minimum Override Value: {Min_Override_Val}\n"
add_info += f"Maximum Override: {Max_Override}\n"
if Max_Override:
    add_info += f"Maximum Override Value: {Max_Override_Val}\n"
add_info += f"Error Override: {Error_Override}\n"
if Error_Override:
    add_info += f"Error Override Value: {Error_Override_Val}\n"
add_info += f"Enhancing Score Threshold: {enhancing_threshold}\n"
add_info += f"Neutral Score Threshold: {neutral_threshold}\n"
add_info += f"Rheostat Score Threshold: {rheostat_threshold}\n"
add_info += f"Toggle Score Threshold: {toggle_threshold}\n"

# Log error calculation warning
if Error_Override and log_scale:
    print("Warning: LOG ERROR IS CALCULATED USING WILD TYPE VALUE")

# Reading in data file
try:
    df = pd.read_csv(input_file_name)
except FileNotFoundError:
    print(f"Error: The file '{input_file_name}' was not found.")
    sys.exit(1)

# Ensure the dataframe has the expected columns
expected_columns = ["Position", "Substitution", "Value", "Error"]
if not all(column in df.columns for column in expected_columns):
    print(f"Error: The input file must contain the following columns: {expected_columns}")
    sys.exit(1)

df.columns = expected_columns

# Convert 'Position' to string and drop rows where 'Position' or 'Value' is NaN
df["Position"] = df["Position"].astype(str)
df = df.dropna(subset=["Position", "Value"])

# Convert 'Substitution' to string and replace NaNs with space
df["Substitution"] = df["Substitution"].astype(str)
df["Substitution"].fillna(" ", inplace=True)

# Get unique positions excluding "WT"
positions_analyzed = df["Position"].unique()
positions_analyzed = positions_analyzed[positions_analyzed != "WT"]

fatal = False
if (df["Position"] == "WT").sum() != 0:
    wild_value = df.loc[df["Position"] == "WT", "Value"].values[0]
    wild_error = df.loc[df["Position"] == "WT", "Error"].values[0]
    df = df[df["Position"] != "WT"]
else:
    print("Error: NO WILD TYPE GIVEN IN DATA")
    fatal = True

if not Error_Override and df["Error"].isna().sum() > 0:
    print("Warning: ERRORS MISSING FROM ORIGINAL DATASET")
    print("Error: NO ERROR GIVEN IN DATA AND NOT OVERRIDDEN")
    fatal = True

if not fatal:
    for position in positions_analyzed:
        df = pd.concat([df, pd.DataFrame([[position, "WT", wild_value, wild_error]], columns=df.columns)], ignore_index=True)

    df["Value"] = pd.to_numeric(df["Value"])

    if Error_Override:
        df["Error"] = Error_Override_Val
    else:
        df["Error"] = pd.to_numeric(df["Error"])

    # Always calculate log values for all data
    df["log_val"] = np.log10(df["Value"])
    df["log_err"] = 0.434 * df["Error"] / df["Value"]

    # If error override is enabled, recalculate log error using wild type value
    if Error_Override:
        df["log_err"] = 0.434 * Error_Override_Val / wild_value

    count_positions = df["Position"].value_counts().reset_index()
    count_positions.columns = ["Position", "Variants"]

    if (count_positions.loc[count_positions["Variants"] < 5, "Position"] != "WT").any():
        print("Warning: The Following positions have fewer than five variants:")
        print(count_positions.loc[count_positions["Variants"] < 5, "Position"])

    # Determine min and max values based on log scale
    if log_scale:
        max_val = df["log_val"].max()
        min_val = df["log_val"].min()
        wild_value_calc = np.log10(wild_value)
        wild_error_calc = 0.434 * wild_error / wild_value if not Error_Override else 0.434 * Error_Override_Val / wild_value
    else:
        max_val = df["Value"].max()
        min_val = df["Value"].min()
        wild_value_calc = wild_value
        wild_error_calc = wild_error if not Error_Override else Error_Override_Val

    # Store global min and max for later use
    global_min = min_val
    global_max = max_val

    # Apply Min Override if enabled
    if Min_Override:
        min_val = np.log10(Min_Override_Val) if log_scale else Min_Override_Val

    # Apply Max Override if enabled
    if Max_Override:
        max_val = np.log10(Max_Override_Val) if log_scale else Max_Override_Val

    # Calculate bin number based on error and count positions
    if log_scale:
        combined_err = np.concatenate([df["log_err"].loc[df["Substitution"] != "WT"].values, df["log_err"].loc[df["Substitution"] == "WT"].values[:1]])
    else:
        combined_err = np.concatenate([df["Error"].loc[df["Substitution"] != "WT"].values, df["Error"].loc[df["Substitution"] == "WT"].values[:1]])
    
    if len(combined_err) == 0:
        print("Error: No error values found to calculate bin number.")
        sys.exit(1)

    bin_number_err = (max_val - min_val) / (np.mean(combined_err) * 2)
    bin_number_count_mean = count_positions.loc[count_positions["Variants"] >= 5, "Variants"].mean()
    bin_number_count_median = count_positions.loc[count_positions["Variants"] >= 5, "Variants"].median()

    bin_number = int(min([bin_number_err, bin_number_count_mean, bin_number_count_median]))

    # Override bin number if enabled
    if bin_override:
        bin_number = bin_override

    bin_number = min(bin_number, 20)

    # Override min and max values in the dataframe
    if log_scale:
        if Min_Override:
            df["min_overridden"] = 0
            df.loc[df["log_val"] < min_val, "min_overridden"] = 1
            df.loc[df["log_val"] < min_val, "log_val"] = min_val
            if df["min_overridden"].sum() > 0:
                add_info += f"\nThe number of values below the Minimum: {df['min_overridden'].sum()}"

        if Max_Override:
            df["max_overridden"] = 0
            df.loc[df["log_val"] > max_val, "max_overridden"] = 1
            df.loc[df["log_val"] > max_val, "log_val"] = max_val
            if df["max_overridden"].sum() > 0:
                add_info += f"\nThe number of values above the Maximum: {df['max_overridden'].sum()}"
    else:
        if Min_Override:
            df["min_overridden"] = 0
            df.loc[df["Value"] < min_val, "min_overridden"] = 1
            df.loc[df["Value"] < min_val, "Value"] = min_val
            if df["min_overridden"].sum() > 0:
                add_info += f"\nThe number of values below the Minimum: {df['min_overridden'].sum()}"

        if Max_Override:
            df["max_overridden"] = 0
            df.loc[df["Value"] > max_val, "max_overridden"] = 1
            df.loc[df["Value"] > max_val, "Value"] = max_val
            if df["max_overridden"].sum() > 0:
                add_info += f"\nThe number of values above the Maximum: {df['max_overridden'].sum()}"

    add_info += f"\nNumber of Bins Used: {bin_number}\n"

    # Determine bin size and bins
    bin_size = (max_val - min_val) / bin_number
    bins = np.linspace(min_val, max_val, bin_number + 1)

    # Calculate neutral-specific bin size - always 2x standard bin size
    if neutral_binsize_override > 0:
        neutral_bin_size = neutral_binsize_override
        if log_scale:
            print("Warning: Neutral bin size override is assumed to be in log scale units.")
        add_info += f"Neutral Bin Size (User Override): {neutral_bin_size}\n"
        # Calculate neutral bin range centered around WT
        neutral_lower = wild_value_calc - (neutral_bin_size / 2)
        neutral_upper = wild_value_calc + (neutral_bin_size / 2)
        add_info += f"Neutral Bin Range: [{neutral_lower:.6f}, {neutral_upper:.6f}]\n\n"
    elif Error_Override or not df["Error"].isna().any():
        # Use 4x WT error
        neutral_bin_size = 4 * wild_error_calc
        neutral_lower = wild_value_calc - (neutral_bin_size / 2)
        neutral_upper = wild_value_calc + (neutral_bin_size / 2)
        add_info += f"Neutral Bin Size (4× WT error): {neutral_bin_size:.6f}\n"
        add_info += f"Neutral Bin Range: [{neutral_lower:.6f}, {neutral_upper:.6f}]\n\n"
    else:
        # Use 2x standard bin size
        neutral_bin_size = 2 * bin_size
        neutral_lower = wild_value_calc - bin_size
        neutral_upper = wild_value_calc + bin_size
        add_info += f"Standard Bin Size: {bin_size:.6f}\n"
        add_info += f"Neutral Bin Size (2× standard): {neutral_bin_size:.6f}\n"
        add_info += f"Neutral Bin Range: [{neutral_lower:.6f}, {neutral_upper:.6f}]\n\n"

    # Determine dead type bin and value
    if DeadValue == "Min":
        dt_bin = 0
        dt_val = min_val
    else:
        dt_bin = len(bins) - 2
        dt_val = max_val

    # Determine wild type bin
    wt_val_for_binning = wild_value_calc
    for i in range(1, bin_number + 2):
        if wt_val_for_binning < bins[i]:
            wt_bin = i - 1
            break

    # Set default weights for bins
    weights = np.full(bin_number, 3)

    # Adjust weights for bins adjacent to wild type and dead type
    for weight in range(bin_number):
        if abs(weight - dt_bin) == 1:
            weights[weight] = 2
        if abs(weight - wt_bin) == 1:
            weights[weight] = 2
        if weight == dt_bin:
            weights[weight] = 1
        if weight == wt_bin:
            weights[weight] = 1

    bin_info = pd.DataFrame(columns=["Position", "Variants", "Enhancing", "Neutral", "Unweighted_Rheostat", 
                                     "Weighted_Rheostat", "Toggle", "Binary", "Average", "Std_Deviation", "Assignment"])
    positions_analyzed = sorted([int(pos) for pos in positions_analyzed])

    # Use the appropriate values based on log_scale setting
    value_column = "log_val" if log_scale else "Value"

    for position in positions_analyzed:
        # Create histogram with better formatting
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(df.loc[df["Position"] == str(position), value_column], bins=bins)
        ax.axvline(wild_value_calc, color="green", linewidth=3)
        ax.axvline(dt_val, color="red", linewidth=3)
        ax.set_xlabel("Log10(Value)" if log_scale else "Value", fontsize=14, fontfamily='Arial')
        ax.set_ylabel("Number of Variants", fontsize=14, fontfamily='Arial')
        ax.set_title(f"Position {position}", fontsize=16, fontfamily='Arial')
        ax.tick_params(axis='both', labelsize=12)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontfamily('Arial')
        
        # Save with prefix
        if plot_prefix:
            filename = os.path.join(plots_dir, f"{plot_prefix}_Histogram-Position{position}.jpg")
        else:
            filename = os.path.join(plots_dir, f"Histogram-Position{position}.jpg")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()

        # Get position data
        position_data = df.loc[df["Position"] == str(position)]
        position_values = position_data[value_column].values
        position_non_wt = df.loc[(df["Position"] == str(position)) & (df["Substitution"] != "WT"), value_column]
        
        # Calculate Enhancing score based on neutral bin boundaries
        if DeadValue == "Min":
            # Enhancing = variants above neutral upper limit
            enhancing_count = (position_non_wt > neutral_upper).sum()
        else:  # Max
            # Enhancing = variants below neutral lower limit
            enhancing_count = (position_non_wt < neutral_lower).sum()
        
        enhancing = enhancing_count / len(position_non_wt) if len(position_non_wt) > 0 else 0

        # Calculate Neutral score - variants within neutral bin range
        neutral_count = ((position_non_wt >= neutral_lower) & (position_non_wt <= neutral_upper)).sum()
        neutral = neutral_count / len(position_non_wt) if len(position_non_wt) > 0 else 0

        # Calculate Toggle
        bin_counts = np.histogram(position_non_wt, bins=bins)[0]
        toggle = bin_counts[dt_bin] / len(position_non_wt) if len(position_non_wt) > 0 else 0

        # Calculate Rheostat scores
        bin_counts_with_wild = np.histogram(position_values, bins=bins)[0]
        unweighted_rheostat = (bin_counts_with_wild > 0).sum() / len(bin_counts)
        weighted_rheostat = ((bin_counts_with_wild > 0) * weights).sum() / weights.sum()
        
        # Calculate Binary
        binary = (bin_counts_with_wild > 0).sum() == 2

        # Calculate statistics for assignment
        position_mean = position_non_wt.mean()
        position_std = position_non_wt.std()
        
        # For output: if log scale, keep as is; if not log scale, use original values
        if log_scale:
            output_mean = position_mean
            output_std = position_std
        else:
            output_mean = position_mean
            output_std = position_std

        # Determine Assignment
        assignment = "Unclassified"
        
        if enhancing > enhancing_threshold:
            assignment = "Enhancing"
        elif neutral > neutral_threshold:
            assignment = "Neutral"
        elif toggle > toggle_threshold:
            assignment = "Toggle"
        elif weighted_rheostat > rheostat_threshold:
            assignment = "Rheostat"
        elif (0.4 < toggle < toggle_threshold and 
              0.4 < neutral < neutral_threshold and 
              weighted_rheostat < rheostat_threshold):
            assignment = "WT/inactive"
        else:
            # Check for Moderate or Adverse
            if position_std < (global_max - global_min) / 2:
                dist_to_dt = abs(position_mean - dt_val)
                dist_to_wt = abs(position_mean - wild_value_calc)
                
                if dist_to_dt > dist_to_wt:
                    assignment = "Moderate"
                elif dist_to_dt < dist_to_wt:
                    assignment = "Adverse"

        bin_info = pd.concat([bin_info, pd.DataFrame([{
            "Position": position,
            "Variants": len(position_values),
            "Enhancing": enhancing,
            "Neutral": neutral,
            "Unweighted_Rheostat": unweighted_rheostat,
            "Weighted_Rheostat": weighted_rheostat,
            "Toggle": toggle,
            "Binary": binary,
            "Average": output_mean,
            "Std_Deviation": output_std,
            "Assignment": assignment
        }])], ignore_index=True)

    # Recombine data excluding wild type and then add a single instance of wild type
    df = pd.concat([df[df["Substitution"] != "WT"], df[df["Substitution"] == "WT"].iloc[0:1]])

    # Create "All" histogram with better formatting
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(df[value_column], bins=bins)
    ax.axvline(wild_value_calc, color="green", linewidth=3)
    ax.axvline(dt_val, color="red", linewidth=3)
    ax.set_xlabel("Log10(Value)" if log_scale else "Value", fontsize=14, fontfamily='Arial')
    ax.set_ylabel("Number of Variants", fontsize=14, fontfamily='Arial')
    ax.set_title("All", fontsize=16, fontfamily='Arial')
    ax.tick_params(axis='both', labelsize=12)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily('Arial')
    
    # Save with prefix
    if plot_prefix:
        filename = os.path.join(plots_dir, f"{plot_prefix}_Histogram-All.jpg")
    else:
        filename = os.path.join(plots_dir, "Histogram-All.jpg")
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    all_values = df[value_column].values
    all_non_wt = df.loc[df["Substitution"] != "WT", value_column]
    
    # Calculate Enhancing score for All based on neutral bin boundaries
    if DeadValue == "Min":
        # Enhancing = variants above neutral upper limit
        enhancing_count = (all_non_wt > neutral_upper).sum()
    else:
        # Enhancing = variants below neutral lower limit
        enhancing_count = (all_non_wt < neutral_lower).sum()
    
    enhancing = enhancing_count / len(all_non_wt) if len(all_non_wt) > 0 else 0

    # Calculate Neutral score for All - variants within neutral bin range
    neutral_count = ((all_non_wt >= neutral_lower) & (all_non_wt <= neutral_upper)).sum()
    neutral = neutral_count / len(all_non_wt) if len(all_non_wt) > 0 else 0

    # Calculate Toggle
    bin_counts = np.histogram(all_non_wt, bins=bins)[0]
    toggle = bin_counts[dt_bin] / len(all_non_wt) if len(all_non_wt) > 0 else 0

    # Calculate Rheostat scores
    bin_counts_with_wild = np.histogram(all_values, bins=bins)[0]
    unweighted_rheostat = (bin_counts_with_wild > 0).sum() / len(bin_counts)
    weighted_rheostat = ((bin_counts_with_wild > 0) * weights).sum() / weights.sum()
    binary = (bin_counts_with_wild > 0).sum() == 2

    # Calculate statistics for assignment
    all_mean = all_non_wt.mean()
    all_std = all_non_wt.std()
    
    # For output: if log scale, keep as is; if not log scale, use original values
    if log_scale:
        output_all_mean = all_mean
        output_all_std = all_std
    else:
        output_all_mean = all_mean
        output_all_std = all_std

    # Determine Assignment for All
    assignment = "Unclassified"
    
    if enhancing > enhancing_threshold:
        assignment = "Enhancing"
    elif neutral > neutral_threshold:
        assignment = "Neutral"
    elif toggle > toggle_threshold:
        assignment = "Toggle"
    elif weighted_rheostat > rheostat_threshold:
        assignment = "Rheostat"
    elif (0.4 < toggle < toggle_threshold and 
          0.4 < neutral < neutral_threshold and 
          weighted_rheostat < rheostat_threshold):
        assignment = "WT/inactive"
    else:
        if all_std < (global_max - global_min) / 2:
            dist_to_dt = abs(all_mean - dt_val)
            dist_to_wt = abs(all_mean - wild_value_calc)
            
            if dist_to_dt > dist_to_wt:
                assignment = "Moderate"
            elif dist_to_dt < dist_to_wt:
                assignment = "Adverse"

    bin_info = pd.concat([bin_info, pd.DataFrame([{
        "Position": "All",
        "Variants": len(all_values),
        "Enhancing": enhancing,
        "Neutral": neutral,
        "Unweighted_Rheostat": unweighted_rheostat,
        "Weighted_Rheostat": weighted_rheostat,
        "Toggle": toggle,
        "Binary": binary,
        "Average": output_all_mean,
        "Std_Deviation": output_all_std,
        "Assignment": assignment
    }])], ignore_index=True)

    # Mark positions with insufficient variants
    insufficient_mask = bin_info["Variants"] < 5
    bin_info.loc[insufficient_mask, "Variants"] = "Error: Insufficient Variants for Calculation"
    for calculation in ["Enhancing", "Neutral", "Unweighted_Rheostat", "Weighted_Rheostat", "Toggle", "Binary", "Average", "Std_Deviation", "Assignment"]:
        bin_info.loc[insufficient_mask, calculation] = np.nan

    print(bin_info)
    print(f"\nPlots saved to: {plots_dir}/")
    
    # Create output filename with prefix
    if plot_prefix:
        output_filename = f"{plot_prefix}_RheoScale-output.csv"
    else:
        output_filename = f"output-{input_file_name}"
    
    # Count assignments (excluding "All" row and insufficient variant positions)
    valid_assignments = bin_info[(bin_info["Position"] != "All") & 
                                 (bin_info["Variants"] != "Error: Insufficient Variants for Calculation")]
    assignment_counts = valid_assignments["Assignment"].value_counts().to_dict()
    
    # Define all possible position types in order
    position_types = ["Enhancing", "Neutral", "Toggle", "Rheostat", "Moderate", "Adverse", "WT/inactive", "Unclassified"]
    
    # Ensure all position types have a count (0 if not present)
    for ptype in position_types:
        if ptype not in assignment_counts:
            assignment_counts[ptype] = 0
    
    # Create assignment summary as a string for the CSV header
    assignment_summary = "\n\nAssignments,Number of positions\n"
    for ptype in position_types:
        assignment_summary += f"{ptype},{assignment_counts[ptype]}\n"
    
    # Write output file with assignment summary
    with open(output_filename, "w") as f:
        f.write(add_info)
        f.write(assignment_summary)
        f.write("\n")  # Extra line before main data
    bin_info.to_csv(output_filename, mode="a", index=False)
    print(f"Output saved to: {output_filename}")
    
    # Create pie chart with specified colors
    colors = {
        "Neutral": "#70AD47",
        "Moderate": "#C6E0B4",
        "Rheostat": "#FFFF00",
        "Adverse": "#E5CCFF",
        "Toggle": "#6600CC",
        "Enhancing": "#FF0000",
        "WT/inactive": "#D3D3D3",
        "Unclassified": "#808080"
    }
    
    # Filter out position types with 0 counts for the pie chart
    pie_data = {ptype: assignment_counts[ptype] for ptype in position_types if assignment_counts[ptype] > 0}
    
    if pie_data:  # Only create pie chart if there's data
        fig, ax = plt.subplots(figsize=(10, 8))
        
        labels = list(pie_data.keys())
        sizes = list(pie_data.values())
        pie_colors = [colors[label] for label in labels]
        
        # Create pie chart without percentages on slices
        wedges, texts = ax.pie(sizes, labels=None, colors=pie_colors, startangle=90)
        
        # Calculate percentages for legend
        total = sum(sizes)
        legend_labels = [f'{label}: {size} ({size/total*100:.1f}%)' for label, size in zip(labels, sizes)]
        
        # Create legend with percentages
        ax.legend(wedges, legend_labels, title="Position Classes", loc="center left", 
                 bbox_to_anchor=(1, 0, 0.5, 1), fontsize=11, title_fontsize=12)
        
        ax.set_title("Position Class Distribution", fontsize=16, fontfamily='Arial', fontweight='bold')
        
        # Save pie chart
        if plot_prefix:
            pie_filename = f"{plot_prefix}_position_class_PieChart.png"
        else:
            pie_filename = "position_class_PieChart.png"
        
        plt.savefig(pie_filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Pie chart saved to: {pie_filename}")
    else:
        print("No valid assignments to create pie chart.")
