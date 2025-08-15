# 2BDP User Manual

2BDP (Biomarker discovery process at binomial decision point) is an algorithm that utilizes machine learning 
and regression modeling to identify and validate the significance of markers related to a binomial point of interest.

Please Cite the following paper:
```
Chakraborty, Nabarun, et al. “Biomarker Discovery Process at Binomial Decision Point (2BDP): Analytical Pipeline to Construct Biomarker Panel.”
Computational and Structural Biotechnology Journal, vol. 21, 1 Jan. 2023, pp. 4729–4742, https://doi.org/10.1016/j.csbj.2023.09.025.
```

## Installation
1. Download most recent R library (TwoBDP_#.#.#.tar.gz) from GitHub under releases
2. Install and load library
```
install.packages("../TwoBDP_1.0.0.tar.gz",repos = NULL)
library(TwoBDP)
```

## How To Use

1. Setup data table (df) <br />
![picture alt](https://github.com/AlexLaw978/2BDP/blob/main/images/data.png)
	- Data should be normalized and cleaned before running create2BDPClass(). 
	- All values need to be numeric.
 	- rownames = sample identifiers.
  	- colnames = names of factors used to build the model (Ex: genes).
   		- Make sure names wont break formulas
2. Setup metadata table (md) <br />
![picture alt](https://github.com/AlexLaw978/2BDP/blob/main/images/metadata.png)
	- minimun of 1 column which contains the point of interest factor.
 		- Make sure factor values wont break formulas
 	- rownames = must exist within data.
  		- Names need to be identicial in order to properly merge results.
	- decisionPointColumn = metadata column name that contains the 2 factor data ("poi1").
 	- pointOfInterest = Case value within the bdp column that is being investigated. 
4. Create 2bdpObject<br />
bdp=createTwoBDPObject(data = df,metadata = md,
                   decisionPointColumn="poi",pointOfInterest="TRUE",threads = 10
                   )
5. Run 2bdp<br />
bdp=panelGeneration(batchSize = 50,bdp = bdp)
bdp=frequencyMap(bdp = bdp)
bdp=PanelExpansion(bdp = bdp)
bdp=panelValidation(bdp = bdp)

7. Results are found in bdp$validatedFeatures[[method]] or if saveDirectory was used: [method]_ValidatedFeatures.csv

```
if(F){
  install.packages("../TwoBDP_1.0.0.tar.gz",repos = NULL)
}

library(TwoBDP)

df=read.csv("TestData/df.csv")
md=read.csv("TestData/md.csv")

#setup for 2bdp
rownames(df)=df[,1];df=df[,-1]
colnames(df)=paste("V",1:ncol(df),sep="") #to make sure the names are usable
rownames(md)=md[,1]

bdp=createTwoBDPObject(data = df,metadata = md,saveDirectory="D:/Windows/Projects/CustomLibraries/TestArea",
                   decisionPointColumn="poi",pointOfInterest="TRUE",threads = 10,
                   panelGenerationParameters=TwoBDP_makePanelParameterList(totalPanels = 100)
                   )

bdp=panelGeneration(batchSize = 50,bdp = bdp) #can take a long time depending on amount of columns in data
bdp=frequencyMap(bdp = bdp)
bdp=PanelExpansion(bdp = bdp)
bdp=panelValidation(bdp = bdp)
```

### Running Validation Only
A formal method is in development, but its possible to manually skip to validation

1. Follow steps 1-4 in "How To Use" to generate 2bdp object
2. Create vector for all models you want tested
	- Ex: combinationVector=c("id1+id2+id3","id1+id3","..."+...)
 	- The names used should exist in colnames(bdp$data).
3. Add vector to 2bdp object.<br />
bdp=TwoBDP_addFeatureCombinations(bdp,combinationVector)
4. Validate panels.<br />
bdp=panelValidation(bdp = bdp)
5. Results are found in bdp$validatedFeatures[[method]] and if saveDirectory was used: [method]_ValidatedFeatures.csv

## Options For 2BDP Class
The class is made up of several required and optional parameters which will be used to operate the algorithm.<br />
The following table describes each parameter where recommend values are placed where applicable and have been tested.

Variable | Description
------------- | -------------
data | Data frame of features "columns" and samples "rows" to analyze. Rownames should be sample identifiers which are all in the rownames in metadata.
metadata | Data frame which has a 2 factor column which will be used as the point of interest.
saveDirectory | A string path where results will be saved.
decisionPointColumn | A string for the column name in metadata, which is the 2 factor for analysis.
pointOfInterest | A string for the non-control value in decisionPointColumn.
seedParms | A named list that contains the "seed" and "RNGkind". see ?RNGkind for more details and options.
threads | An int for the amount of cores use in parallel.
panelGenerationParameters | A list of required parameters to run the random forest. See TwoBDP_makePanelParameterList function for more details.
frequencyMapParameters | A list of required parameters to run the frequency map. See TwoBDP_makeFrequencyParameterList function for more details.
validationParameters | A list of required parameters to run the validation. See TwoBDP_makeValidationParameterList function for more details.


## Developer Note
2BDP core algorithms are complete, however optimizations for program, parallel, and deployment architecture is still underway.
### Major Development Goals
- [X] Improve/Redesign architecture to allow serial and parallel operations obtain the same results.
- [X] Create functions to skip sections of algorithm
- [X] Implement method to extract model objects for graphing purposes.
- [X] Fix file data loss when overwriting existing excel file at checkpoints.
- [X] Format program into R package.
- [ ] Linux Testing.

### Minor Development Goals
- [X] Dynamic or static options for seeds.
- [X] Dynamic or static options for analysis on split data.
- [ ] Improve and modelurize messages
