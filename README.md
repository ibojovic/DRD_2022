# Pipeline Summary

## Step 1: Data Initialization
* Loaded raw data into minfi.
* Stored data in an RGChannelSet object named RGset.
## Step 2: Fluorescence DataFrames
* Created DataFrames Red and Green to store the corresponding fluorescence data.
## Step 3: Fluorescence Table Completion
* Completed a table with Red and Green fluorescence values for a specific address.
* Identified the probe type from the manifest file and recorded its color if it was a Type I probe.
## Step 4: Preprocessing
* Constructed an object MSet.raw for further analysis.
## Step 5: Quality Control
* Generated QC plots (QCplot).
* Examined intensity of negative controls using minfi.
* Calculated detection p-values, recording the number of probes above the threshold for each sample.
## Step 6: Methylation Value Analysis
* Calculated raw beta and M values.
* Plotted density of mean methylation values, comparing WT and MUT groups.
* Assessed any observable differences between the two groups.
## Step 7: Data Normalization
* Normalized the data with a designated function.
* Compared raw and normalized data through various plots and commented on the observed changes.
* Evaluated the suitability of the normalization approach for this dataset and compared WT and MUT group distributions.
## Step 8: Principal Component Analysis
* Conducted PCA on the matrix of normalized beta values.
* Commented on the sample division based on group, sex, or batch.
## Step 9: Differential Methylation Analysis
* Identified differentially methylated probes between WT and MUT groups.
## Step 10: Statistical Testing
* Applied multiple testing corrections.
* Established a significance threshold of 0.05 and identified differentially methylated probes using nominal p-values, Bonferroni, and BH corrections.
## Step 11: Visualization of Differential Methylation
* Created a volcano plot and a Manhattan plot to visualize differential methylation results.
## Step 12: Heatmap Generation
* Generated a heatmap for the top 100 differentially methylated probes.
