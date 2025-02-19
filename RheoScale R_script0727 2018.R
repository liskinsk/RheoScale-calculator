# RheoScale: A Tool to Aggregate and Quantify Experimentally-Determined Substitution Outcomes for Multiple Variants at Individual Protein Positions
# Abby M. Hodges, Aron W. Fenton, Larissa L. Dougherty, Andrew C. Overholt and Liskin Swint-Kruse
# Department of Natural, Health, and Mathematical Sciences, MidAmerica Nazarene University, Olathe, Kansas, USA
# and The Department of Biochemistry and Molecular Biology, The University of Kansas Medical Center, Kansas City, Kansas, USA
# July 27, 2018
# 
# How to use this RheoScale calculator:
#   The theory, rationale, and specifics of data analyses utilized in this calculator are explained in more detail in the associated paper and supplemental information.
# 1) Data Input File
# A csv template file is included as a supplemental file. The csv file for inputting data requires that the data be entered in sequential columns with (1) containing the amino acid position number as an integer, (2) can contain the variant amino acid for that position although this column is not used in analysis and can be left blank, (3) the functional data value and (4) the error for each data value (this column can also be left blank if an error override is entered below).  Additionally, the wild type value with appropriate error must be entered into the csv file with "WT" in the first column.
# The functional value in the third column should be the data being analyzed using this calculator.  This number can represent any type of value from any assay.  The only limitations are zero or negative values if the data is to be converted to a log scale.  In these instances, the zero and negative scores should be replaced with a numerical value one order of magnitude smaller than the smallest measurable value. If removing or replacing negative values has significant implications to the interpretation of the data, it is recommended that the user consider not converting the data to a log scale. 
# The error entered into the fourth column for each value is used to determine the optimal number of bins for the entire data set.  If these values are not included, a number must be entered into the error override line in order for the calculator to determine a recommendation for the number of bins. If the data is converted to a log scale, the errors are also converted accordingly using the formula: measured error * 0.434/functional value 
# The calculator will automatically extract the numbers of the positions varied and report them into the first column of the output file.  These values should be examined to ensure that they are consistent with the expected positions varied from the data set.  
# 
# 2) Log Scale Conversion
# For data that spans 1 or 2 orders of magnitude, histogram bins should be calculated using nontransformed (linear) calculations.  If a data set covers more than three orders of magnitude, it should be converted to a log scale by entering "TRUE" in this line. Any functional value that has already been converted to a log scale (e.g., Gibbs free energy) or doesn't cover more than two orders of magnitude (e.g., Hill number), the user should enter "FALSE" in this line. A data set that spans between 2 and 3 orders of magnitude should be carefully considered to determine whether a log scale is appropriate or not. All override values entered below should be entered in the non-log version (unless the user already converted the data values manually).
# 
# 3) Minimum/Maximum override
# The minimum and maximum values for the entire data set are necessary to determine the bins and the bin weighting for the calculator.  The min and max values are calculated automatically and used in the calculations. If a data set is known to have a different min or max value for this protein or assay, those min or max values can be provided as overrides by changing the Min or Max override from "FALSE" to "TRUE" and then providing a value in the next line.  If any data point falls outside of the maximum or minimum value, it will be reassigned to be the defined minimum/maximum value. The number of data points that exceed the maximum value are reported in the output file. 
# 
# 4) Error override
# The error override option is provided for data sets that do not have an error value associated with each functional value.  Alternatively, one error value may better represent the error inherent in the experimental methodology and can be entered in this line. If error override is desired, the "FALSE" can be changed to "TRUE" and the value in entered in the following line.  This value is used in determining a recommended number of bins for analyzing the data set. If error override is included for a data set that is converted to a log scale (see below) within the calculator, then the error value entered will be propagated using the wild type as the reference value. The formula for error propagation for log calculations in this case is 0.434*error override/wt functional value. If an alternate approach to error propagation is desired (percent error, etc) then the user should amend the data set before including in the calculator.
# 
# 5) Bin Number override
# The recommended number of bins is determined through a combination of the average error for the data set as well as the total number of variants at each position and is reported in the output file.  A perfect rheostat would occupy 20 bins, but error and representation at each position minimize the number of bins that can be used.  The algorithm for the bin number recommendation can be found in the "Bin Number Generation" section and is explained in further detail in the manuscript and the Supplemental Information.  If a different number of bins is desired, the bin override line should be changed from "FALSE" to "TRUE" and then the desired bin number entered into the line below.
# 
# 6) Nonfunctional or dead value
# The value that represents a completely nonfunctional protein ("dead") must be defined as either the minimum or maximum value used to analyze the data set. In this line, the user must enter either "MIN" or "MAX". If this value is not considered and adjusted accordingly with every calculation, the scores will be inaccurate. 

#User Defined values:
#################################################################################

#Input File Name
input_file_name <- "LacI template.csv"

#Log Scale Conversion
log_scale <- TRUE

#Minimum Override
Min_Override <- FALSE
Min_Override_Val <- 0


#Maximum Override
Max_Override <- FALSE
Max_Override_Val <- 0


#Error Override
Error_Override <- FALSE
Error_Override_Val <- 0


#Bin Number Override
bin_override <- 0


#Nonfunctional or Dead value
NonFuncProteinValue <- "Max"


#Adding User Settings to Output File
###################################################################################

add_info <- ""
add_info <- paste(add_info, paste("Input File Name: ",input_file_name))
add_info <- paste(add_info, paste("\nLog Scale: ",as.character(log_scale)))
add_info <- paste(add_info, paste("\nBin Number Override: ",as.character(bin_override)))
add_info <- paste(add_info, paste("\nNonfunctional or Dead value: ",as.character(NonFuncProteinValue)))
add_info <- paste(add_info, paste("\nMinimum Override: ",as.character(Min_Override)))
if (Min_Override) {
  add_info <- paste(add_info, paste("\nMinimum Override Value: ",as.character(Min_Override_Val)))
}
add_info <- paste(add_info, paste("\nMaximum Override: ",as.character(Max_Override)))
if (Max_Override) {
  add_info <- paste(add_info, paste("\nMaximum Override Value: ",as.character(Max_Override_Val)))
}
add_info <- paste(add_info, paste("\nError Override: ",as.character(Error_Override)))
if (Error_Override) {
  add_info <- paste(add_info, paste("\nError Override Value: ",as.character(Error_Override_Val)))
}
##################################################################################

if (Error_Override && log_scale) {
warning("LOG ERROR IS CALCULATED USING WILD TYPE VALUE")
}

#Reading in data file
##################################################################################

df <- read.csv(input_file_name)

names(df) <- c("Position","Substitution","Value","Error")

df$Position <- as.character(df$Position)

df <- df[!is.na(df$Position) & !is.na(df$Value),]

df$Substitution <- as.character(df$Substitution)
df$Substitution[is.na(df$Substitution)] <- " "

positions_analyzed <- unique(df$Position)
positions_analyzed <- positions_analyzed[positions_analyzed != "WT"]

fatal <- FALSE
if (sum(df$Position == "WT") != 0) {
  #Look for it in the dataframe
  wild_value <- df$Value[df$Position == "WT"]
  wild_error <- df$Error[df$Position == "WT"]
  df <- df[df$Position != "WT",]
} else {
  stop("NO WILD TYPE GIVEN IN DATA")
  fatal <- TRUE
}

#Check for Errors in Original Data
#################################################################################

if (!Error_Override & sum(is.na(df$Error)) > 0) {
  warning("ERRORS MISSING FROM ORIGINAL DATASET")
}

if (!Error_Override & sum(is.na(df$Error)) > 0) {
  stop("NO ERROR GIVEN IN DATA AND NOT OVERRIDDEN")
  fatal <- TRUE
}


if (!fatal) {
#Creation of Data frame 
##################################################################################

for (position in positions_analyzed) {
  df <- rbind(df, c(position, "WT", wild_value, wild_error))
}

df$Value <- as.numeric(df$Value)

if (Error_Override) {
  df$Error <- Error_Override_Val
} else {
  df$Error <- as.numeric(df$Error)
}


df$log_val <- log10(df$Value)
df$log_err <- 0.434*df$Error/df$Value

if (Error_Override) {
  df$Error <- Error_Override_Val
  df$log_err <- 0.434*Error_Override_Val/wild_value
}

count_positions <- data.frame(table(df$Position))

names(count_positions) <- c("Position","Variants")

if (nrow(count_positions[count_positions$Variants < 5 & count_positions$Positions != "WT" & count_positions$Positions != " ",])>0) {
  warning("The Following positions have fewer than five variants:")
  warning(count_positions[count_positions$Variants < 5 & count_positions$Positions != "WT" & count_positions$Positions != " ",])
}

#Min/Max Override 
###################################################################

if (log_scale) {
  max_val <- max(df$log_val)
  min_val <- min(df$log_val)
} else {
  max_val <- max(df$Value)
  min_val <- min(df$Value)
}

if (Min_Override) {
  if (log_scale) {
    min_val <- log10(Min_Override_Val)
  } else {
    min_val <- Min_Override_Val
  }
}

if (Max_Override) {
  if (log_scale) {
    max_val <- log10(Max_Override_Val)
  } else {
    max_val <- Max_Override_Val
  }
}

#Bin Number Generation
############################################################################

if (log_scale) {
  bin_number_err <- (max_val-min_val)/(mean(c(df$log_err[df$Substitution != "WT"],df$log_err[df$Substitution =="WT"][0]))*2)
} else {
  bin_number_err <- (max_val-min_val)/(mean(c(df$Error[df$Substitution != "WT"],df$Error[df$Substitution =="WT"][0]))*2)
}


bin_number_count_mean <- mean(count_positions$Variants[count_positions$Variants >= 5])
bin_number_count_median <- median(count_positions$Variants[count_positions$Variants >= 5])

bin_number <- floor(min(c(bin_number_err,bin_number_count_mean, bin_number_count_median)))

if (bin_override) {
  bin_number <- bin_override
}

bin_number <- min(bin_number, 20)

#Bin Min/Max Override
###################################################################

if(log_scale) {
  if (Min_Override) {
    df$min_overridden <- 0
    df$min_overridden[df$log_val < min_val] <- 1
    df$log_val[df$log_val < min_val] <- min_val
    if (sum(df$min_overridden) > 0) {
      add_info <- paste(add_info, paste("\n\nThe number of bins below the Minimum was: ",as.character(sum(df$min_overridden))))
    }
  }

  if (Max_Override) {
    df$max_overridden <- 0
    df$max_overridden[df$log_val > max_val] <- 1
    df$log_val[df$log_val > max_val] <- max_val
    if (sum(df$max_overridden) > 0) {
      add_info <- paste(add_info, paste("\n\nThe number of bins above the Maximum was: ",as.character(sum(df$max_overridden))))
    }
  }
} else {
  if (Min_Override) {
    df$min_overridden <- 0
    df$min_overridden[df$Value < min_val] <- 1
    df$Value[df$Value < min_val] <- min_val
    if (sum(df$min_overridden) > 0) {
      add_info <- paste(add_info, paste("\n\nThe number of bins below the Minimum was: ",as.character(sum(df$min_overridden))))
    }
  }

  if (Max_Override) {
    df$max_overridden <- 0
    df$max_overridden[df$Value > max_val] <- 1
    df$Value[df$Value > max_val] <- max_val
    if (sum(df$max_overridden) > 0) {
      add_info <- paste(add_info, paste("\n\nThe number of bins above the Maximum was: ",as.character(sum(df$max_overridden))))
    }
  }
}

add_info <- paste(add_info, paste("\nNumber of Bins Used:",as.character(bin_number)))
add_info <- paste(add_info, "\n\n")

#Bin Size calculation
###########################################################################

if (log_scale) {
  bin_size <- (max_val-min_val)/bin_number
  bins <- seq(min_val,max_val,bin_size)
} else {
  bin_size <- (max_val-min_val)/bin_number
  bins <- seq(min_val,max_val,bin_size)
}


#Determine Location of Death Type Bin and Wild Type Bin
#############################################################################

if (NonFuncProteinValue == "Min") {
  dt_bin <- 1
  dt_val <- min_val
} else {

  dt_bin <- length(bins)-1
  dt_val <- max_val
}


if (log_scale) {
  for (i in 1:bin_number+1) {
    if (log10(wild_value) < bins[i]) {
      wt_bin <- i-1
      break
    }
  }
} else {
  for (i in 1:bin_number+1) {
    if (wild_value < bins[i]) {
      wt_bin <- i-1
      break
    }
  }  
}

#Weights for Weighted Rheostat
#################################################################

#This is the default bin weight (3)
weights <- rep(3, bin_number)

for (weight in 1:bin_number) {
  if (abs(weight-dt_bin) == 1) {
    #This ist he weight of any bind which is adjacent tot he dead type
    weights[weight] <- 2
  }
  if (abs(weight-wt_bin) == 1) {
    #This is the weight of any bin which is adjacent to the wild type
    weights[weight] <- 2
  }
  if (weight == dt_bin) {
    #This is the weight of the dead type bin
    weights[weight] <- 1
  }
  
  if (weight == wt_bin) {
    #This is the weight of the wild type bin
    weights[weight] <- 1
  }

}





#Final Calculations
####################################################################

bin_info <- data.frame(Position=character(0),Variants=numeric(0),Neutral=numeric(0),Unweighted_Rheostat=numeric(0),Weighted_Rheostat=numeric(0),Toggle=numeric(0),Binary=logical(0))
positions_analyzed <- sort(as.numeric(as.character(positions_analyzed)))
for (position in positions_analyzed) {
  if (log_scale) {
    #Vertical line for Wild Type and for Dead Type
    jpeg(paste("Histogram-Position",paste((position),".jpg")))
    hist(df$log_val[df$Position==position], breaks = bins, xaxt='n', main=paste("Position ",(position)), xlab="Functional Value", ylab="Number of Variants")
    abline(v=log10(wild_value),col="green", lwd=3)
    abline(v=dt_val,col="red", lwd=3)
    axis(side=1, at=bins, labels=sprintf("%0.2f", bins))
    dev.off()
    bin_counts <- hist(df$log_val[df$Position==position & df$Substitution != "WT"], breaks = bins, plot=FALSE)$counts
    wild_bin_count <- sum(df$log_val[df$Position==position] >= log10(wild_value)-bin_size/2 & df$log_val[df$Position==position] <= log10(wild_value)+bin_size/2)
    wild_bin_count <- wild_bin_count-1
    neutral <- wild_bin_count/(sum(bin_counts))
    toggle <- bin_counts[dt_bin]/(sum(bin_counts))
    bin_counts_with_wild <- hist(df$log_val[df$Position==position], breaks = bins, plot=FALSE)$counts
    unweighted_rheostat <- sum(bin_counts_with_wild > 0)/length(bin_counts)
    weighted_rheostat <- sum((bin_counts_with_wild > 0)*weights)/sum(weights)
    binary <- sum(bin_counts_with_wild > 0) == 2
    
    bin_info <- rbind(bin_info, data.frame(Position=position,Variants=sum(bin_counts_with_wild),Neutral=neutral,Unweighted_Rheostat=unweighted_rheostat,Weighted_Rheostat=weighted_rheostat,Toggle=toggle,Binary=binary))
  } else {
    jpeg(paste("Histogram-Position",paste((position),".jpg")))
    hist(df$Value[df$Position==position], breaks = bins, xaxt='n', main=paste("Position ",(position)), xlab="Functional Value", ylab="Number of Variants")
    abline(v=wild_value,col="green", lwd=3)
    abline(v=dt_val,col="red", lwd=3)
    axis(side=1, at=bins, labels=sprintf("%0.2f", bins))
    dev.off()
    bin_counts <- hist(df$Value[df$Position==position & df$Substitution != "WT"], breaks = bins, plot=FALSE)$counts
    wild_bin_count <- sum(df$Value[df$Position==position] >= (wild_value)-bin_size/2 & df$Value[df$Position==position] <= (wild_value)+bin_size/2)
    wild_bin_count <- wild_bin_count-1
    neutral <- wild_bin_count/(sum(bin_counts))
    toggle <- bin_counts[dt_bin]/(sum(bin_counts))
    bin_counts_with_wild <- hist(df$Value[df$Position==position], breaks = bins, plot=FALSE)$counts
    unweighted_rheostat <- sum(bin_counts_with_wild > 0)/length(bin_counts)
    weighted_rheostat <- sum((bin_counts_with_wild > 0)*weights)/sum(weights)
    binary <- sum(bin_counts_with_wild > 0) == 2
    
    bin_info <- rbind(bin_info, data.frame(Position=position,Variants=sum(bin_counts_with_wild),Neutral=neutral,Unweighted_Rheostat=unweighted_rheostat,Weighted_Rheostat=weighted_rheostat,Toggle=toggle,Binary=binary))
  }
}

#All, need to include 1 wild type
df <- rbind(df[df$Substitution != "WT",],df[df$Substitution == "WT",][1,])
if (log_scale) {
  #Vertical line for Wild Type and for Dead Type
  jpeg("Histogram-All.jpg")
  hist(df$log_val, breaks = bins, xaxt='n', main="All", xlab="Functional Value", ylab="Number of Variants")
  abline(v=log10(wild_value),col="green", lwd=3)
  abline(v=dt_val,col="red", lwd=3)
  axis(side=1, at=bins, labels=sprintf("%0.2f", bins))
  dev.off()
  bin_counts <- hist(df$log_val[df$Substitution != "WT"], breaks = bins, plot=FALSE)$counts
  wild_bin_count <- sum(df$log_val >= log10(wild_value)-bin_size/2 & df$log_val <= log10(wild_value)+bin_size/2)
  wild_bin_count <- wild_bin_count-1
  neutral <- wild_bin_count/(sum(bin_counts))
  toggle <- bin_counts[dt_bin]/(sum(bin_counts))
  bin_counts_with_wild <- hist(df$log_val, breaks = bins, plot=FALSE)$counts
  unweighted_rheostat <- sum(bin_counts_with_wild > 0)/length(bin_counts)
  weighted_rheostat <- sum((bin_counts_with_wild > 0)*weights)/sum(weights)
  binary <- sum(bin_counts_with_wild > 0) == 2
  
  bin_info <- rbind(bin_info, data.frame(Position="All",Variants=sum(bin_counts_with_wild),Neutral=neutral,Unweighted_Rheostat=unweighted_rheostat,Weighted_Rheostat=weighted_rheostat,Toggle=toggle,Binary=binary))
} else {
  jpeg("Histogram-All.jpg")
  hist(df$Value, breaks = bins, xaxt='n', main="All", xlab="Functional Value", ylab="Number of Variants")
  abline(v=wild_value,col="green", lwd=3)
  abline(v=dt_val,col="red", lwd=3)
  axis(side=1, at=bins, labels=sprintf("%0.2f", bins))
  dev.off()
  bin_counts <- hist(df$Value[df$Substitution != "WT"], breaks = bins, plot=FALSE)$counts
  wild_bin_count <- sum(df$Value >= (wild_value)-bin_size/2 & df$Value <= (wild_value)+bin_size/2)
  wild_bin_count <- wild_bin_count-1
  neutral <- wild_bin_count/(sum(bin_counts))
  toggle <- bin_counts[dt_bin]/(sum(bin_counts))
  bin_counts_with_wild <- hist(df$Value, breaks = bins, plot=FALSE)$counts
  unweighted_rheostat <- sum(bin_counts_with_wild > 0)/length(bin_counts)
  weighted_rheostat <- sum((bin_counts_with_wild > 0)*weights)/sum(weights)
  binary <- sum(bin_counts_with_wild > 0) == 2
  
  bin_info <- rbind(bin_info, data.frame(Position="All",Variants=sum(bin_counts_with_wild),Neutral=neutral,Unweighted_Rheostat=unweighted_rheostat,Weighted_Rheostat=weighted_rheostat,Toggle=toggle,Binary=binary))
}
#Summary Score Plot
jpeg("SummaryScores.jpg")

barplot(t(as.matrix(bin_info[,c("Neutral","Unweighted_Rheostat","Weighted_Rheostat","Toggle")])), beside = TRUE,
        col = c("lightblue", "mistyrose", "lightcyan",
                "lavender"),
        legend = c("Neutral","Unweighted_Rheostat","Weighted_Rheostat","Toggle"),
        names.arg=bin_info$Position,
        main="Summary of All Scores")
dev.off()

#Summary Histogram
UR <- rep(0,10)
WR <- rep(0,10)
Neu <- rep(0,10)
Tog <- rep(0,10)
for (i in seq(1,10,by=1)) {
  for (value in bin_info[bin_info$Position != "All",]$Unweighted_Rheostat) {
    if (value <= i*0.1 & value > ((i-1)*0.1)) {
      UR[i] <- UR[i]+1
    }
  }
  for (value in bin_info[bin_info$Position != "All",]$Weighted_Rheostat) {
    if (value <= i*0.1 & value > ((i-1)*0.1)) {
      WR[i] <- WR[i]+1
    }
  }
  for (value in bin_info[bin_info$Position != "All",]$Neutral) {
    if (value <= i*0.1 & value > ((i-1)*0.1)) {
      Neu[i] <- Neu[i]+1
    }
  }
  for (value in bin_info[bin_info$Position != "All",]$Toggle) {
    if (value <= i*0.1 & value > ((i-1)*0.1)) {
      Tog[i] <- Tog[i]+1
    }
  }
}


for (value in bin_info[bin_info$Position != "All",]$Unweighted_Rheostat) {
  if (value == 0) {
    UR[1] <- UR[1]+1
  }
}
for (value in bin_info[bin_info$Position != "All",]$Weighted_Rheostat) {
  if (value == 0) {
    WR[1] <- WR[1]+1
  }
}
for (value in bin_info[bin_info$Position != "All",]$Neutral) {
  if (value == 0) {
    Neu[1] <- Neu[1]+1
  }
}
for (value in bin_info[bin_info$Position != "All",]$Toggle) {
  if (value == 0) {
    Tog[1] <- Tog[1]+1
  }
}

jpeg("SummaryHist.jpg")
barplot(t(as.matrix(data.frame(Neu,UR,WR,Tog))), beside = TRUE,
        col = c("lightblue", "mistyrose", "lightcyan",
                "lavender"),
        legend = c("Neutral","Unweighted_Rheostat","Weighted_Rheostat","Toggle"),
        names.arg=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1"),
        main="Histogram of All Scores")



dev.off()

#Bar plot of four scores by position
jpeg("UnweightedRheostat.jpg")
barplot(bin_info$Unweighted_Rheostat,names.arg=bin_info$Position, main="Unweighted Rheostat", xlab="Position",ylab="Value")
dev.off()
jpeg("WeightedRheostat.jpg")
barplot(bin_info$Weighted_Rheostat,names.arg=bin_info$Position, main="Weighted Rheostat", xlab="Position",ylab="Value")
dev.off()
jpeg("Toggle.jpg")
barplot(bin_info$Toggle,names.arg=bin_info$Position, main="Toggle", xlab="Position",ylab="Value")
dev.off()
jpeg("Neutral.jpg")
barplot(bin_info$Neutral,names.arg=bin_info$Position, main="Neutral", xlab="Position",ylab="Value")
dev.off()

#barplot(df$Value, grouping="Variants")
bin_info[bin_info$Variants < 5,c("Variants")] <- "Error: Insufficient Variants for Calculation"
for (calculation in c("Neutral","Unweighted_Rheostat","Weighted_Rheostat","Toggle","Binary")) {
  bin_info[bin_info$Variants == "Error: Insufficient Variants for Calculation", c(calculation)] <- NA
}

print(bin_info)
write(add_info,file=paste("output-",input_file_name))
suppressWarnings(write.table(bin_info, file=paste("output-",input_file_name),append=TRUE,sep=",",row.names=FALSE))
}