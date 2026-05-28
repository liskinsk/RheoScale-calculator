# RheoScale-calculator
The Rheoscale calculator is used to analyze large datasets of protein substitutions in order to identify rheostat positions.
The 2018 version of the calculator is available as Excel(TM) files or in R. 
RheoScale 2.0 (releaseded in 2026) is available as Excel (TM) files or as a Python code.

The first worksheet in the Excel workbook contains information about how to use the calculator.  The Python code has an associated README file.
The Python code is available for download here or using PIP (Pythons installer)

For the 2018 version, please cite:

RheoScale: A tool to aggregate and quantify experimentally determined substitution outcomes for multiple variants at individual protein positions
Abby M Hodges, Aron W Fenton, Larissa L Dougherty, Andrew C Overholt, Liskin Swint-Kruse.  Hum Mutat. 2018 Dec;39(12):1814-1826. doi: 10.1002/humu.23616. 

A citation will be forthcoming for RheoScale 2.0.

Other citations associated with the 2026 version of the RheoScale calculator include
1. Identification of biochemically neutral positions in liver pyruvate kinase. 
Tyler A Martin, Tiffany Wu, Qingling Tang, Larissa L Dougherty, Daniel J Parente, Liskin Swint-Kruse, Aron W Fenton
Proteins. 2020 Oct;88(10):1340-1350. doi: 10.1002/prot.25953. Epub 2020 Jun 27.
PMID: 32449829 PMCID: PMC8990530 DOI: 10.1002/prot.25953

2. Rheostat functional outcomes occur when substitutions are introduced at nonconserved positions that diverge with speciation.
Swint-Kruse L, Martin TA, Page BM, Wu T, Gerhart PM, Dougherty LL, Tang Q, Parente DJ, Mosier BR, Bantis LE, Fenton AW.
Protein Sci. 2021 Sep;30(9):1833-1853. doi: 10.1002/pro.4136. Epub 2021 Jun 11.  PMID: 34076313 


# Python RheoScale-calculator Infomation

## Usage/Quick Start

Rheoscale2.0 has been implemented 3 different formats. 
1. As a Python Package (highly recommened; this github site)
2. As a Command line interface (CLI; this github site)
3. As an Excel workbook, which can be downloaded at https://github.com/liskinsk/RheoScale-calculator

## About the Calculations

#### Position classes:  
At least 5 amino acid variants are required for each position to be analyzed.  RheoScale first uses the value and error for the wild-type variant to determine whether a position is “neutral” or “non-neutral”.  The variant sets for non-neutral positions are then assessed with a modified histogram analysis. 

The hierarchical logic used to define position classes follows the order: Neutral>Toggle/Enhancing>Rheostat>WT-Inactive_split>Moderate/Adverse.

The descriptions of these position classes, as well as the publications describing their default score thresholds, are:
>>Neutral: At least 70% of substitutions have WT-like outcomes for the parameter measured.  The neutral score for each position, which is used to assign neutral positions, is calculated separately from all other scores, by using a neutral bin that is centered on the wild-type value and is (usually) independent of the histogram bin size.
Reference: Martin, Tyler A., Tiffany Wu, Qingling Tang, Larissa L. Dougherty, Daniel J Parente, Liskin Swint-Kruse, and Aron W. Fenton. 2020. 'Identification of biochemically neutral positions in liver pyruvate kinase', Proteins: Structure, Function, and Bioinformatics, 88: 1340-50 

>>Moderate rheostat: Substitutions have non-neutral outcomes, but the set of substitutions samples less than half of the possible range AND the values are closer to WT than to the "dead" end of the range. If positions are close to the thresholds for moderate, adverse, or WT/Inactive split, the assignment flagged, indicating that a manual inspection of the data is needed to determine the best position class.
Reference: Swint-Kruse, L., T. A. Martin, B. M. Page, T. Wu, P. M. Gerhart, L. L. Dougherty, Q. Tang, D. J. Parente, B. R. Mosier, L. E. Bantis, and A. W. Fenton. 2021. 'Rheostat functional outcomes occur when substitutions are introduced at nonconserved positions that diverge with speciation', Protein Sci, 30: 1833-53.

>>Rheostat: Simplistically, the position's set of substitutions samples at least half of the possible functional range.  When the "weighted" rheostat score is used, the sampling might be a little less than half of the range, but contributions from variants with partial loss-of-function are weighted more heavily, since they provide more confidence that a position is a rheostat position.   The weighted score has been used in almost all studies to date.  Rheostat behavior cannot be detected from the average of a position's substitution values.  
Reference: Hodges, A. M., A. W. Fenton, L. L. Dougherty, A. C. Overholt, and L. Swint-Kruse. 2018. 'RheoScale: A tool to aggregate and quantify experimentally determined substitution outcomes for multiple variants at individual protein positions', Hum Mutat, 39: 1814-26.

>>Adverse:  Like Moderate rheostat positions, substitutions have non-nuetral outcomes and the set of substitutions samples less than half of the possible range.  However, the set of values are closer to the "dead" end of the range than to WT. If positions are close to the thresholds for moderate, adverse, or WT/Inactive split, the assignment flagged, indicating that a manual inspection of the data is needed to determine the best position class.
Reference: Sreenivasan, S., Fontes, J. D., and Swint-Kruse, L. 2025. 'Dissecting the effects of single amino acid substitutions in SARS-CoV-2 Mpro', Protein Science 34, e70225.

>>Toggle: Around 2/3 of the substitutions lack detectable activity.
Reference: Miller, M., Y. Bromberg, and L. Swint-Kruse. 2017. 'Computational predictors fail to identify amino acid substitution effects at rheostat positions', Sci Rep, 7: 41329.

>>WT/Inactive split: Half of the substitutions have WT-like outcomes and the other half lack detectable activity.  This may be a hallmark of altered protein stability. These are a special case of "binary" positions (for which substitutions fall in only 2 bins).If positions are close to the thresholds for moderate, adverse, or WT/Inactive split, the assignment flagged, indicating that a manual inspection of the data is needed to determine the best position class.
Reference: Page, B. M., T. A. Martin, C. L. Wright, L. A. Fenton, M. T. Villar, Q. Tang, A. Artigues, A. Lamb, A. W. Fenton, and L. Swint-Kruse. 2022. 'Odd one out? Functional tuning of Zymomonas mobilis pyruvate kinase is narrower than its allosteric, human counterpart', Protein Sci, 31: e4336.

>>Enhancing: At least 80% of substitutions enhance the measured parameter relative to the upper limit of the neutral bin.
Reference: Sreenivasan, S., Fontes, J. D., and Swint-Kruse, L. 2025. 'Dissecting the effects of single amino acid substitutions in SARS-CoV-2 Mpro', Protein Science 34, e70225.

>>Assignment warning flags.  Assignments will be flagged when positions are close to the thresholds for moderate, adverse, or WT/Inactive split.  In these cases, manual inspection of the position's histogram is needed to determine the best position class.


#### Guidelines for histogram analyses: 
a) For linear-scale data that spans 1 order of magnitude, histogram bins should be calculated using non-transformed (linear) calculations.  If a data set covers more than two orders of magnitude, it should be converted to a log scale. Any functional  value that has already been converted to a log scale (e.g., Gibbs free energy) or doesn't cover more than two orders of magnitude (e.g., Hill number) should not be converted. A data set that spans between 1 and 2 orders of magnitude should be carefully considered to determine whether a log scale is appropriate or not.  Note that most high-throughput data (e.g., "deep mutational scanning") are already reported in log scale, and should not be converted.

b) The value that represents a completely nonfunctional protein ("dead") must be defined as either the minimum or maximum value used to analyze the data set. 

c) The default width of the neutral bin is set to 4 times the WT error (or 4 times the error override), to represent the WT value +/- 2 standard deviations. If no error is designated, the default is set to 2 times the rheostat/toggle bin size.  See Martin et al, 2020, for the statistical reasoning.  The width of the neutral bin can also be set by the user.  

d) The histogram range is the most important parameter for assigning position phenotypes. The range is determined from the minimum and maximum values of the dataset and/or known for the assay.  The min and max values of each dataset are calculated automatically.  If a data set is known to have min or max values that differ from the experimental dataset, override values can be entered.  Various examples that require override values are as follows:
-	For datasets with variants that are "better" than wild-type, it can be useful to set the relevant max/min override so that the last bin is populated by at least 5% of the total variants; this prevents a tiny set of highly-active variants from making the range artificially large and thus dominating RheoScale assignments.
-	The max/min override cells are used when the dataset being analyzed spans a smaller range of functional values than is known to be possible (e.g., when the dataset does not contain a “dead” variant.)  
-	The max/min overrides may also be used when investigators wish to designate a “dead” threshold that falls inside the measurable assay range (e.g., any variant with less than 10% activity should be considered “dead”). 
-	Examples for estimating "dead" values are described in the following publications (as well as others from the Swint-Kruse lab):  1.  Hodges, Fenton, Dougherty, Overholt, and Swint-Kruse. 2018. 'RheoScale: A tool to aggregate and quantify experimentally determined substitution outcomes for multiple variants at individual protein positions', Hum Mutat, 39: 1814-26. 2.  Sreenivasan, Fontes, and Swint-Kruse. 2025. 'Dissecting the effects of single amino acid substitutions in SARS-CoV-2 Mpro', Protein Science 34, e70225. 3.  O'Neil, Swint-Kruse, and Fenton (2024) Rheostatic contributions to protein stability can obscure a position's functional role, Protein Science 33, e5075.
-	If any data point falls outside of the max or min override value, it will be reassigned to be the override value for calculations. 

e)  An error override option  is provided for data sets that do not have an error value associated with each functional value.  Alternatively, one error value may better represent the error inherent in the experimental methodology.  This value is used in determining a recommended number of bins for analyzing the data set.  If error override is included for a data set that is converted to a log scale within the calculator, then the error value entered will be propagated using the wild type as the reference value. The formula for error propagation for log calculations in this case is 0.434*error/[WT value]. If an alternate approach to error propagation is desired (percent error, etc) then the user should amend the data set before including in the calculator.

f) The recommended number of bins is determined through a combination of the average error for the data set as well as the total number of variants at each position.  A perfect rheostat position would occupy 20 bins, but error and the number of variants available for each position constrain the number of bins that should be used.  The algorithm for the bin number recommendation c is explained in further detail in Hodges et al., 2018.  If a different number of bins is desired, that number can be entered with an override.  Empirically, 10 bins appear to work well for many datasets. Iterations with different bin numbers show that many position assignments are not very sensitive to the bin number.

g) Weighted versus unweighted rheostat scores.  Unweighted rheostat scores are calculated as the fraction of bins occupied by at least one variant. Weighted rheostat scores take into account the fact that variants that fall in bins farther from the "wild-type" and "dead" bind are more likely to be "true" intermediate outcomes than variants in bins adjacent to "wild-type" or "dead" bins. In default calculations of weighted rheostat scores, bins containing “wild-type” and “dead” values are assigned a weighting value of 1, bins adjacent to the “wild-type” and “dead” bins are assigned a weight of 2, and all other bins were weighted with a value of 3.  The weighted rheostat score is then calculated from the sum of the weighted values of each occupied bin divided by the sum of the weighted values of all bins. The final weighted score ranges from near 0.0 (no rheostat behavior) to 1.0 (perfect rheostat behavior).  The Microsoft® Excel version contains “override” cells for investigator-determined weighting factors.  In the python version, the weighting factors must be changed in the script. In general, the Swint-Kruse group has chosen to use weighted rheostat scores, and the weighted scores are used in RheoScale2.0 to make position assignments.



## Installation Instructions
### If you are new to Python
How to use this Python version of the RheoScale calculator
Python Requirements:  You need Python 3 or later and a few standard Python packages.  If you don’t have Python installed:
-	Go to https://www.python.org/downloads/
-	Download the latest Python 3 version (or higher).
-	During installation, check the box that says “Add Python to PATH”.
To confirm it’s installed, open a terminal (Command Prompt on Windows, Terminal on macOS/Linux), and run:
```bash 
python --version
```
You should see something like Python 3.10.6.

Installing required packages:  Open a terminal (Command Prompt on Windows, Terminal on macOS/Linux), and run:
```bash 
pip install pandas numpy matplotlib
```
That’s all the setup you need !!

### To install RheoScale2.0
#### To install the Python version as a package run (recommended)
```bash
python pip install rheoscale
```
This installs rheoscale, making it available for use in any Python script via import rheoscale

#### To install the Python version for CLI use OR in an editable  script
Clone this repo using git (if you need to intall git on your machine go here: https://git-scm.com/)  
then run:
```bash 
git clone https://github.com/CarterG-17/Rheoscale_2.0.git
```
This will install and create a new folder named "Rheoscale_2.0" with all the files found in this repository.   
run: 
```bash 
cd Rheoscale_2.0
```
Then you will need to install rheoscale into python using these files.  
In the newly created "Rheoscale_2.0" folder run:
```bash 
python -m pip -e .
```
Now you should be able to use the python CLI interface or make changes to the Rheoscale calculator if you wish to, using this locally installed package.

## RheoScale_2.0 Usage
### Python Package

```python 

from rheoscale.rheoscale_config import RheoscaleConfig
from rheoscale.rheoscale_runner import RheoscaleRunner
import pandas as pd
```


If you already have your data loaded as a pandas DataFrame
```python
data #<- the data you had

#create config 
config = RheoscaleConfig('Protein name')


#run script
runner = RheoscaleRunner(config, data)
position_df =runner.run() #returns a Dataframe with positions calucations and numbers

```
If you need Rheoscale to load the data from a csv
```python


#create config 
config = RheoscaleConfig('Protein name', input_file_name=data_path)


#run script
runner = RheoscaleRunner(config)
position_df =runner.run() #returns a Dataframe with positions calucations and numbers

```
Simply put there are 2 steps to running the python rheoscale

#1. Creating a Configuration (RheoscaleConfig) object that sets all the parameters for analysis

#2. Creating a Runner (RheoscaleRunner) object that ensures that all set parameters of the config make sense and then a .run() command can be run


### CLI version
To see the inputs of the CLI version run:
```bash
python -m rheoscale.cli -h
```

To use the command line interface of rheoscale you can run: 
```bash
python -m rheoscale.cli protein_name --input_file (--opitional_inputs)
```
### Excel version

To use the Excel version, download the Excel workbook from https://github.com/liskinsk/RheoScale-calculator.
To see more information please read the "How to use this calculator" worksheet in this workbook.

## Documentation

Setting up the configuration object can be the most detailed task for running rheoscale
 

**Comming Soon: Jupyter UI with widgets...**

### RheoscaleConfig()

**RheoscaleConfig** (**protein_name**, ***input_file_name***=*None*, ***number_of_positions***=*None*, ***log_scale***=*False*, ***WT_val***=*None*, ***WT_error***=*None*, ***WT_name***=*None*, ***min_val***=*None*, ***max_val***=*None*, ***error_val***=*None*, ***number_of_bins***=*None*, ***dead_extremum***=*'Min'*, ***neutral_binsize***=*None*, ***output_dir***=*None*, ***output_histogram_plots***=*False*, ***even_bins***=*True*, ***columns***=*mappingproxy({'position': 'Position', 'substitution': 'Substitution', 'value': 'Value', 'error': 'Error'}))*

#### Required Parameters:


#### Optional Parameters
**input_file_name : str (path to CSV file)** <br>
If passing in a mutational data from a CSV file provide the path to the file here
this file must contain 4 columns:  Position, Substitution, Value, Error (these do not have to be the names of the columns see **columns** parameter)

**number_of_positions : int (Defalt= None)**<br>
Use this to check that Rheoscale sees the same amount of positions as you expect

**log_scale : bool (Defalt= False)**<br>
Any data set that spans more than three orders of magnitude should be converted to a log scale. Data sets that contain negative or zero values will result in errors if converted to a log scale; substitute with a value 10x smaller than the smallest value.

**WT_val: float (Default=None)** <br>
If WT values are not found in the DataFrame or CSV file given to Rheoscale you can add them here. 

**WT_error : float (Default=None)** <br>
If WT values are given through the config you add the error here.

**WT_name : str (Default=None)**<br>
if WT values are in DataFrame or CSV but do not have the label of "WT" in the position column the alternate name can be given here (*e.g.,* if the values are named "wild-type")

**min_val : float (Default=None)**<br>
A data set may not contain the absolute minimum or maximum value associated with the experimental data for this protein. The researcher may enter a min value in the boxes below to override the min or max value found in the data set

**max_val : float (Default=None)**<br>
A data set may not contain the absolute minimum or maximum value associated with the experimental data for this protein. The researcher may enter a max value in the boxes below to override the min or max value found in the data set

**error_val : float (Default=None)**<br>
If data set does not have errors associated with each data point or if the overall data set is better represented by a predetermined error value, that value can be entered into this override box.

**number_of_bins : int (Default=None)**<br>
If the recommended bin number is too small or large to provide meaningful data, a different bin number may be entered as an override. We recommend bin numbers between 7-13 with 10 being a good starting point. Further insight is provided in the manuscript (**citation**).

**dead_extremum : Literal["Min", "Max"] (Default='Min')**<br>
Determine if the nonfunctional (dead) value associated with this protein and assay represents a minimum or a maximum extremum. 

**neutral_binsize : float (Default=None)**<br>
Override the rheoscale calculated neutral binsize value 

**output_dir : str (Default='Rheoscale_analysis')**<br>
name of the directory created by rheoscale

**output_histogram_plots : bool (Default=False)**<br>
If False (default) will only output the "all" histogram for the full dataset. If True will output the histogram of every position. 

**even_bins : bool (Default=True)**<br>
If True (default), will generating bins with uniform widths. If False, will make "dead bin" larger to reflect the size of that bin based on the data

**columns : bool (Default=dict({'position': 'Position', 'substitution': 'Substitution', 'value': 'Value', 'error': 'Error'}))**<br>
If your columns names do not match out names. Then creatr a dict that maps the that names of the your titles to each column title. Of note, the keys of this dict must always be: position, substitution, value, and error
*e.g.,*:
{"position": "Position",
 "substitution": 'Mutation',
 "value": 'Functional Value',
 "error": 'Error'}

#### Methods:

**from_json(path_to_json)**<br>
can take in a JSON file with these type of inputs

### RheoscaleRunner()
**RheoscaleRunner** (**config**, ***DMS_data***=*None*)

#### Reqired Parameters:
**config : RheoscaleConfig**<br>
a configuration from the first step

#### Optional Parameters:
**DMS_data : pd.DataFrame (Default= None)**<br>
if an input file is not given in the pandas DataFrame then a data frame can be added here


## License

This package is under a GNU AFFERO GENERAL PUBLIC LICENSE

## Contact 
For scientific questions about the dataset, please contact Dr. Liskin Swint-Kruse 
(lswint-kruse_at_kumc_dot_edu)

For techniqical Python coding questions, please contact Carter Gray (cartergray-at-ku-dot-edu)

