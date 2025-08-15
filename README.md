# 2BDP User Manual

2BDP (Biomarker discovery process at binomial decision point) is an algorithm that utilizes machine learning 
and regression modeling to identify and validate the significance of markers related to a binomial point of interest.

Please Cite the following paper:
```
Chakraborty, Nabarun, et al. “Biomarker Discovery Process at Binomial Decision Point (2BDP): Analytical Pipeline to Construct Biomarker Panel.”
Computational and Structural Biotechnology Journal, vol. 21, 1 Jan. 2023, pp. 4729–4742, https://doi.org/10.1016/j.csbj.2023.09.025.
```

# UNDER REVISION

## Requirements
The following R libraries must be installed before running:
```
require(stringr)
require(utils)
require(openxlsx)
require(caret)
require(randomForest)
require(pROC)
require(ResourceSelection)
require(parallel)
require(ggplot2)
```
Note: It was developed on Windows OS with parallel features working on Windows. Never tested on Linux-based systems.

## Installation
1. Download R script (2bdp.R) from GitHub
2. Use source("path to script") in script that will use 2bdp algorithm

## 2BDP Operation Flow Chart

![picture alt](https://github.com/AlexLaw978/2BDP/blob/main/images/2bdpFlowChart.png)

## How To Use

1. Setup data table (df) <br />
![picture alt](https://github.com/AlexLaw978/2BDP/blob/main/images/data.png)
	- Data should be normalized and cleaned before running create2BDPClass().
	- All values need to be numeric.
 	- rownames = sample identifiers.
  	- colnames = names of factors used to build the model (Ex: genes).
2. Setup metadata table (md) <br />
![picture alt](https://github.com/AlexLaw978/2BDP/blob/main/images/metadata.png)
	- minimun of 2 columns, one for sample names and point of interest.
 	- sampleIDS = metadata column name containing all samples ("Sample").
  		- Names need to be identicial and ordered the same way as the rownames within the data table.
	- bdp = metadata column name that contains the 2 factor data ("poi1").
 	- poi = Value of interest within the bdp column that is being investigated ("T"). Factors can be named anything but only 2 can exist.
3. Source Script<br />
source(".../2bdp.R")
4. Create 2bdpObject<br />
bdpObject=create2BDPClass(data = df,metadata = md,sampleIDS = "Sample",bdp="poi1",poi="T")
5. Run 2bdp<br />
bdpObject=run2BDP(bdpObject)
6. Results are found under rsbmrValidation, kfcvValidation, bestKFOLD, bestRSMBR

```
rm(list = ls())
setwd("...")

library(openxlsx)
source("2bdp.R")


df=read.xlsx("exampleData/df.xlsx",rowNames = T)
#rowNames = T becauase sample names are first column and are required to be set as rownames
#if the sample names are not the first column then manaully set the rownames of df to that column and then remove
#ex: rownames(df)=df[,num]; df=df[,-num]. where num is the column index of samples  
md=read.xlsx("exampleData/md.xlsx")

bdpObject=create2BDPClass(data = df,metadata = md,sampleIDS = "Sample",bdp = "poi",poi = "F",
                          rfDataSize = ncol(df))

#default rfDataSize uses 50% of all factors due to the expectation that lots of factors will be used
#however this example dataset is small enough where using all factors is acceptable so the additional parameter is changed from default.


bdpObject=run2BDP(bdpObject)
```

### Running Validation Only
A formal method is in development, but its possible to manually skip to validation

1. Follow steps 1-4 in "How To Use" to generate 2bdp object
2. Create vector for all models you want tested
	- Ex: combinationVector=c("id1+id2+id3","id1+id3","..."+...)
 	- The names used should be under bdpObject[["featureNameMap"]]$ids. If orginal names contain no illegal characters and the featureNameMap was never created then use the column names from the data table to generate the formulas.
3. Assign combinationVector to biomarkerCombinations.<br />
bdpObject$biomarkerCombinations=combinationVector
4. Run the following functions as is.<br />
bdpObject=ParallelBiomarkerValidation(bdpObject)<br />
bdpObject=convertNames(bdpObject)<br />
bdpObject=trimByThreshold(bdpObject)<br />
5. Results are found under rsbmrValidation, kfcvValidation, bestKFOLD, bestRSMBR

Note: Subfunctions of run2BDP check the status of the class to allow it to restart at major checkpoints.
Run time can be extremely long depending on the values

## Options For 2BDP Class
The class is made up of several required and optional parameters which will be used to operate the algorithm.<br />
The following table describes each parameter where recommend values are placed where applicable and have been tested.

Variable | Description
------------- | -------------
*Required* data|Dataframe object where samples are rows & factors are columns.<br /> Default=NULL.
*Required* metadata|Dataframe with samples and the binomial decision point. Sample names must exist in their own column, and match the row order of the data frame.<br /> Default=NULL. 
*Required* bdp|The String within the metadata data frame column that will be used for analysis. Only 2 unique factors can exist for this.<br /> Default=NULL.
*Required* poi|The String within dbp column that is used as the point of interest.<br /> Default=NULL.
*Required* sampleIDS|Column name in the metadata data frame that contains the sample IDs.<br /> Default=NULL.
*Optional* fileBasename|A string or path to a file without the file extension is required. If provided, each major step will be saved to an excel file, otherwise, all data will only exist with the 2BDP class object in R. Note if the same file is used multiple times there is a chance for some of the save data to be lost during the overwrite process (Looking into cause).<br /> Default=NULL.
*Optional* featuresNameMap|The following regex "^[A-Za-z0-9]+$" determines if feature names are valid. If feature names are not valid, a reference map (data frame) between default ids and feature names will be made, even if you provide your own map.<br /> Default=NULL.
*Optional* ifFeatureMap| boolean to use the featureMap function. If set to false and invalid characters are present then object will not be created and you will be required to rename features in data frame before proceeding.<br /> Default=T (T or F).
*Optional* validationMethod|pick what validation method you want: "kfcv", "rsbmr", or "all". All will run both validation methods.<br /> Default="all" ("kfcv", "rsbmr", or "all").
*Optional* seed|randomization is used at several steps of this algorithm, where the seed is reset at different key steps. You can change the seed if the default seed does not work for you.<br /> Default=100.
*Optional* kfold|a number used to determine the folds with rsbmr and kfcv validation methods. Must be greater than 1.<br /> Default=10 (>1).
*Optional* threads|The amount of threads you want to use to run everything in parallel by creating the max amount of clusters to be used.<br /> Default=75% of total cores (>1). 
*Optional* totalPanels|This value determines how many times random forest should be run to achieve the requested amount of panels.<br /> Default=2000 (>1).
*Optional* topPanels|This value determines how many panels should be kept after frequency mapping for validation processing.<br /> Default=200 (>1). 
*Optional* trainSplitRatio|The ratio used to determine what samples are to be placed in the training set. Left over samples are defaulted to validation set.<br /> Default=0.7 (0-1).
*Optional* rfDataSize|Can be a float or int. Value is used to determine how many features are used when running random forest. If between 0-1, then the value will be multiplied to the total length of features. If value is >1 then that exact amount of features will be used. Recommend to keep 2x above the minimal panel size.<br /> Default=.5 (0-1, or 1>).
*Optional* amountOfFeatures|a vector of numbers for the size that panels should be. Ex: 2:10 will analyize panels between 2 - 10 features.<br /> Default=2:10 (0<x<ncol(data)). 
*Optional* subGrouping|After the frequency mapping, features in each panel are divided into subpanels. Currently 2 options are available for this: "all" or "sw". All will produce every combination of subpanels for every size provided by the amountOfFeatures vector. "sw" performs a sliding window operation over the full panels.<br /> Default=all.
*Optional* metric|Currently program only operates with the Accuracy metric. This variable is only used in randomForest generation. May be expanded to other metrics after research and testing.<br /> Default=Accuracy (Only "Accuracy").
*Optional* tuneGrid|Object used in randomForest and Caret operations. Can use your own custom object or leave null for a default use.<br /> Default=NULL.
*Optional* tuneGridParms|The value used for .mtry when creating a tuneGrid object.<br /> Default=2:20.
*Optional* trainControl|Object used in randomForest and Caret operations. Can use your own custom object or leave null for default use.<br /> Default=NULL.
*Optional* trainControlParms|a named vector (method, repeats, savePredictions) used to create the default trainControlObject for caret and randomForest.<br /> Default=c(method="repeatedcv",repeats=5,savePredictions="all").
*Optional* AUC|The AUC cutoff to identify significant panels whose value exceeds that cutoff.<br /> Default=0.8 (0-1).
*Optional* pValue|The p value cutoff to identify significant panels whose value is below that cutoff.<br /> Default=.05 (0-1).

## Developer Note
2BDP core algorithms are complete, however optimizations for program, parallel, and deployment architecture is still underway.
### Major Development Goals
- [ ] Improve/Redesign architecture to allow serial and parallel operations obtain the same results.
- [ ] Create functions to skip sections of algorithm
- [ ] Implement method to extract model objects for graphing purposes.
	- [ ] Implement method to retrieve R objects without rerunning whole program.
- [ ] Fix file data loss when overwriting existing excel file at checkpoints.
- [ ] Format program into R package.
- [ ] Linux Testing.
- [ ] Optimize for limited memory systems

### Minor Development Goals
- [ ] Dynamic or static options for seeds.
- [ ] Dynamic or static options for analysis on split data.
- [ ] Add functionality for analyzing multiple point of interests.
- [ ] Add support functions to edit constants like sheetnames.
