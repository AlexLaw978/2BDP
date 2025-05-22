# 2BDP User Manual

2BDP (Biomarker discovery process at binomial decision point) is an algorithm that utilizes machine learning 
and regression modeling to identify and validate the significance of markers related to a binomial point of interest.

Please Cite the following paper:
```
Chakraborty, Nabarun, et al. “Biomarker Discovery Process at Binomial Decision Point (2BDP): Analytical Pipeline to Construct Biomarker Panel.”
Computational and Structural Biotechnology Journal, vol. 21, 1 Jan. 2023, pp. 4729–4742, https://doi.org/10.1016/j.csbj.2023.09.025.
```

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
Note: It was developed on Windows OS with parallel features working on Windows. Parallel operations were untested on Linux-based systems.

## Installation
1. Download R script (2bdp.R) from GitHub
2. Use source("path to script") in script that will use 2bdp algorithm

## 2BDP Operation FlowChart

![picture alt](https://github.com/AlexLaw978/2BDP/blob/main/images/2bdpFlowChart.png)

## How to Use

1. create a 2BDPClass object using the create2BDPClass Function
2. pass the 2BDP object into run2BDP() function and wait to it to finish

Note: Subfunctions of run2BDP check the status of the class to allow it to restart at major checkpoints.
Run time can be extremely long depending on the values

## Options For 2BDP Class
The class is made up of several required variables which will be used to operate the 
The following table describes each variable that should be inputted. Recommend values are placed where applicable and have been tested.

Variable | Description
------------- | -------------
*Required* data|Dataframe object where samples are rows & genes for columns.<br /> Default = NULL.
*Required* metadata|Dataframe with samples and the binomial decision point. Sample names must exist in their own column, and match the row order of the data's data frame.  Default = NULL. 
*Required* bdp|The String within the metadata data frame column that will be used for analysis. Only 2 unique factors can exist for this. Default = NULL.
*Required* poi|The String within the dbp column is used as the point of interest when using the algorithm. The default value is NULL.
*Required* sampleIDS|Column name (String) in the metadata data frame that contains the sample IDs. Default = NULL.
*Optional* fileBasename|A string or path to a file without the file extension is required. If provided, each major step will be saved to an excel file, otherwise, all data will only exist with the 2BDP class object in R. Note if the same file is used multiple times there is a chance for some of the save data to be lost during the overwrite process (Looking into cause). Default = NULL.
*Optional* featuresNameMap|The following regex "^[A-Za-z0-9]+$" determines if feature names are valid. If feature names are not valid, a reference map (data frame) between default ids and feature names will be made, even if you provide your own map. Default=NULL.
*Optional* ifFeatureMap| boolean to use the featureMap function. If set to false and invalid characters are present then object will not be created and you will be required to rename features in data frame before proceeding. Default = T (T or F).
*Optional* validationMethod|pick what validation method you want: "kfcv", "rsbmr", or "all". All will run both validation methods. Default = "all" ("kfcv", "rsbmr", or "all").
*Optional* seed|randomization is used at several steps of this algorithm, where the seed is reset at different key steps. You can change the seed if the default seed does not work for you. Default = 100.
*Optional* kfold|a number used to determine the folds with rsbmr and kfcv validation methods. Must be greater than 1. Default = 10 (>1).
*Optional* threads|The amount of threads you want to use to run everything in parallel by creating the max amount of clusters to be used. Default = 75% of total cores (>1). 
*Optional* totalPanels|This value determines how many times random forest should be run to achieve the requested amount of panels. Default = 2000 (>1).
*Optional* topPanels|This value determines how many panels should be kept after frequency mapping for validation processing. Default 200 (>1). 
*Optional* trainSplitRatio|The ratio used to determine what samples are to be placed in the training set. Left over samples are defaulted to validation set. Default = 0.7 (0-1).
*Optional* rfDataSize|Can be a float or int. Value is used to determine how many features are used when running random forest. If between 0-1, then the value will be multiplied to the total length of features. If value is >1 then that exact amount of features will be used. Recommend to keep 2x above the minimal panel size. Default = .5 (0-1, or 1>).
*Optional* amountOfFeatures|a vector of numbers for the size that panels should be. Ex: 2:10 will analyize panels between 2 - 10 features. Default = 2:10 (0<x<ncol(data)). 
*Optional* subGrouping|After the frequency mapping, features in each panel are divided into subpanels. Currently 2 options are available for this: "all" or "sw". All will produce every combination of subpanels for every size provided by the amountOfFeatures vector. "sw" performs a sliding window operation over the full panels. Default = all.
*Optional* metric|Currently program only operates with the Accuracy metric. This variable is only used in randomForest generation. May be expanded to other metrics after research and testing. Default = Accuracy (Only "Accuracy").
*Optional* tuneGrid|Object used in randomForest and Caret operations. Can use your own custom object or leave null for a default use. Default=NULL.
*Optional* tuneGridParms|The value used for .mtry when creating a tuneGrid object. Default = 2:20.
*Optional* trainControl|Object used in randomForest and Caret operations. Can use your own custom object or leave null for default use. Default=NULL.
*Optional* trainControlParms|a named vector (method, repeats, savePredictions) used to create the default trainControlObject for caret and randomForest. Default =c(method="repeatedcv",repeats=5,savePredictions = "all").
*Optional* AUC|The AUC cutoff to identify significant panels whose value exceeds that cutoff. Default = 0.8 (0-1).
*Optional* pValue|The p value cutoff to identify significant panels whose value is below that cutoff. Default = .05 (0-1).

## Setup Guide
1. Setup data table <br />
![picture alt](https://github.com/AlexLaw978/2BDP/blob/main/images/data.png)
	- Data should be normalized and cleaned before running create2BDPClass().
	- All values need to be numeric.
 	- rownames = sample identifiers.
  	- colnames = names of factors used to build the model (Ex: genes).
2. Get metadata table <br />
![picture alt](https://github.com/AlexLaw978/2BDP/blob/main/images/metadata.png)
	- minimun of 2 columns, one for sample names and point of interest.
 	- sampleIDS = metadata column name containing all samples ("Sample" in this case).
  		- Names need to be identicial and ordered the same as the rownames within the data table.
	- bdp = metadata column name you want to test factors against ("poi1" for example).
 	- poi = Value of interest within the bdp column that is being investigated.
3. create 2bdpObject

## Running 2BDP

run:
```
bdpObject <- create2BDPClass(...) #fill with required data
run2BDP(bdpObject,batch=1)
```

The batch variable is used to calculate how many panels should be processed in parallel before breaking to be saved. This variable is multiplied against threads and used in parLapply() function. If value >1, then some tasks will be put on standby until a thread becomes available. 

## Developer Note
2BDP core algorithms are complete, however optimizations for program, parallel, and deployment architecture is still underway. <br />
Key development steps
- [ ] Dynamic or static options for seeds
- [ ] Dynamic or static options for analysis on split data
- [ ] Improve/Redesign parallel architecture to allow serial and parallel operations obtain the same results
- [ ] Implement method to retrieve R objects without rerunning whole program 
