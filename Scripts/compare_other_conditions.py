import argparse
import pandas as pd
import numpy as np
import os
import os.path

"""
Author: Chris Magnano
2/07/2023

This script takes in raw data for both experiments, normalizes it, calculates fold changes, and uses limma in an R script to get significance values. It also creates input files for phosFate.

Outputs:
2. Files with all fold changes and p-values - 6 files, 3 10-plexes, and protein by phos data.
3. Files with only significantly changed values - 8 files, as we seperate the first 10-plex by timepoint.

"""

verbose = False

def main():
    global verbose

    #Handle command line arguments
    parser = argparse.ArgumentParser(description="Normalize and analyze mass spec data.")
    parser.add_argument("proteinData1", help="This is the file containing proteomic mass spec data for the 1st 10-plex.")
    parser.add_argument("proteinData2", help="This is the file containing proteomic mass spec data for the 2nd 10-plex.")
    parser.add_argument("proteinData3", help="This is the file containing proteomic mass spec data for the 3rd 10-plex.")
    parser.add_argument("phosphoData1", help="This is the file containing phosphoproteomic mass spec data for the 1st 10-plex.")
    parser.add_argument("phosphoData2", help="This is the file containing phosphoproteomic mass spec data for the 2nd 10-plex.")
    parser.add_argument("phosphoData3", help="This is the file containing phosphoproteomic mass spec data for the 3rd 10-plex.")
    parser.add_argument("--outName", default="output", help="This is a prefix which will be added to all output files.")
    parser.add_argument('--verbose', action='store_true', help='If true outputs additional text.' )
    parser.add_argument('--makePrizes', action='store_true', help='If true will also make merged prize lists for each time point. Expects an Uniprot column in all input data' )

    #NOTE: These arguments were included to deal with different labelings which might be more or less conservative.
    #However, DO NOT run this script multiple times with different thresholds on the same data to "see which is best"
    #unless you've talked with a statistician or are a statistician.
    parser.add_argument("--qvalThresh", default="0.1", help="This is the q value cutoff to use for significance.")
    parser.add_argument("--foldThresh", default="1.5", help="This is the fold change cutoff to use for significance.")

    #Get cmd line args
    args = parser.parse_args()
    proteinFile1=args.proteinData1
    proteinFile3=args.proteinData2
    proteinFile2=args.proteinData3
    phosphoFile1=args.phosphoData1
    phosphoFile3=args.phosphoData2
    phosphoFile2=args.phosphoData3

    outName  = args.outName
    verbose = args.verbose
    makePrizes = args.makePrizes
    qvalThresh = float(args.qvalThresh)
    foldThresh = float(args.foldThresh)

    #Load raw data. If future data has different column names, here's the place to update that.
    vPrint("Loading Protein Data")
    protDat1 = loadProteinData(proteinFile1)
    prot1ColNames = ['0minRep1 (TMT 126C PCN)','0minRep2 (TMT 127N PCN)','0minRep3 (TMT 127C PCN)','0minRep4 (TMT 128N PCN)','5minRep1 (TMT 128C PCN)','5minRep2 (TMT 129N PCN)','5minRep3 (TMT 129C PCN)','60minRep1 (TMT 130N PCN)','60minRep2 (TMT 130C PCN)','60minRep3 (TMT 131N PCN)']
    prot1IDColNames = ['Uniprot']

    protDat2 = loadProteinData(proteinFile2)
    prot2ColNames = ['WT_0Min_(TMT 126C PCN)', 'WT_0Min_(TMT 127N PCN)','_Env_0Min_(TMT 127C PCN)','_Env_0Min_(TMT 128N PCN)','WT_5Min_(TMT 128C PCN)','WT_5Min_(TMT 129N PCN)','_Env_5Min_(TMT 129C PCN)','_Env_5Min_(TMT 130N PCN)','_Genome_5Min_(TMT 130C PCN)','_Genome_5Min_(TMT 131N PCN)']
    prot2IDColNames = ['Uniprot']

    protDat3 = loadProteinData(proteinFile3)
    prot3ColNames = ['WT_0Min_(TMT 126C PCN)','WT_0Min_(TMT 127N PCN)','_Genome_0Min_(TMT 127C PCN)','_Genome_0Min_(TMT 128N PCN)','WT_60Min_(TMT 128C PCN)','WT_60Min_(TMT 129N PCN)','_Env_60Min_(TMT 129C PCN)','_Env_60Min_(TMT 130N PCN)','_Genome_60Min_(TMT 130C PCN)','_Genome_60Min_(TMT 131N PCN)']
    prot3IDColNames = ['Uniprot']

    vPrint("Loading Phospho Data")
    phosDat1 = loadPhosphoData(phosphoFile1)
    phos1ColNames = ['0minRep1 (TMT 126C PCN)','0minRep2 (TMT 127N PCN)','0minRep3 (TMT 127C PCN)','0minRep4 (TMT 128N PCN)','5minRep1 (TMT 128C PCN)','5minRep2 (TMT 129N PCN)','5minRep3 (TMT 129C PCN)','60minRep1 (TMT 130N PCN)','60minRep2 (TMT 130C PCN)','60minRep3 (TMT 131N PCN)']
    phos1IDColNames = ['Uniprot','Isoform']

    phosDat2 = loadPhosphoData(phosphoFile2)
    phos2ColNames = ['WT_0Min_(TMT 126C PCN)','WT_0Min_(TMT 127N PCN)','_Env_0Min_(TMT 127C PCN)','_Env_0Min_(TMT 128N PCN)','WT_5Min_(TMT 128C PCN)','WT_5Min_(TMT 129N PCN)','_Env_5Min_(TMT 129C PCN)','_Env_5Min_(TMT 130N PCN)','_Genome_5Min_(TMT 130C PCN)','_Genome_5Min_(TMT 131N PCN)']
    phos2IDColNames = ['Uniprot','Isoform']

    phosDat3 = loadPhosphoData(phosphoFile3)
    phos3ColNames = ['WT_0Min_(TMT 126C PCN)','WT_0Min_(TMT 127N PCN)','_Genome_0Min_(TMT 127C PCN)','_Genome_0Min_(TMT 128N PCN)','WT_60Min_(TMT 128C PCN)','WT_60Min_(TMT 129N PCN)','_Env_60Min_(TMT 129C PCN)','_Env_60Min_(TMT 130N PCN)','_Genome_60Min_(TMT 130C PCN)','_Genome_60Min_(TMT 131N PCN)']
    phos3IDColNames = ['Uniprot','Isoform']

    #Normalize data
    vPrint("Quantile Normalizing Data")
    protDat1 = normalize10Plex(protDat1, prot1ColNames, prot1IDColNames)
    protDat2 = normalize10Plex(protDat2, prot2ColNames, prot2IDColNames)
    protDat3 = normalize10Plex(protDat3, prot3ColNames, prot3IDColNames)
    phosDat1 = normalize10Plex(phosDat1, phos1ColNames, phos1IDColNames)
    phosDat2 = normalize10Plex(phosDat2, phos2ColNames, phos2IDColNames)
    phosDat3 = normalize10Plex(phosDat3, phos3ColNames, phos3IDColNames)

    #Standardize Column Names
    stdColNames(protDat1, prot1ColNames)
    stdColNames(protDat2, prot2ColNames)
    stdColNames(protDat3, prot3ColNames)
    stdColNames(phosDat1, phos1ColNames)
    stdColNames(phosDat2, phos2ColNames)
    stdColNames(phosDat3, phos3ColNames)

    experiments = []
    #Approach 1: Compare change over time within 10-plexes
    experiments.append((protDat2,"env_0min","env_5min","prot"))
    experiments.append((protDat2,"wt_0min","wt_5min","prot"))
    experiments.append((phosDat2,"env_0min","env_5min","phos"))
    experiments.append((phosDat2,"wt_0min","wt_5min","phos"))
    experiments.append((protDat3,"wt_0min","wt_60min","prot"))
    experiments.append((protDat3,"genome_0min","genome_60min","prot"))
    experiments.append((phosDat3,"wt_0min","wt_60min","phos"))
    experiments.append((phosDat3,"genome_0min","genome_60min","phos"))

    #Approach 2: Compare WT directly to other conditions)
    experiments.append((protDat2,"wt_5min","env_5min","prot"))
    experiments.append((phosDat2,"wt_5min","env_5min","phos"))
    experiments.append((protDat3,"wt_60min","genome_60min","prot"))
    experiments.append((phosDat3,"wt_60min","genome_60min","phos"))
    experiments.append((protDat2,"wt_0min","env_0min","prot"))
    experiments.append((phosDat2,"wt_0min","env_0min","phos"))
    experiments.append((protDat3,"wt_0min","genome_0min","prot"))
    experiments.append((phosDat3,"wt_0min","genome_0min","phos"))

    phosExp = []
    protExp = []
    for experiment in experiments:
        if experiment[-1] == "prot":
            protExp.append(experiment)
        elif experiment[-1] == "phos":
            phosExp.append(experiment)

    #Calculate fold changes
    vPrint("Calculating Fold Changes")
    for experiment in experiments:
        calcFoldChanges(experiment)

    #Get test results if available
    gotRes = True
    for experiment in experiments:
        gotRes = gotRes and getTestResults(outName, experiment)

    if not gotRes:
        for experiment in experiments:
            writeExpSetup(outName, experiment)
        vPrint("Run limma-test.R with same output file name argument for each experiment, then run this script again.")
        return

    vPrint("Filtering Significant Results")
    for experiment in experiments:
        filtered_data = filterRes(experiment,qvalThresh,foldThresh)

    #Save merged tables
    vPrint("Saving merged tables")
    makeMergedTables(protExp, outName+"AllProteinData.csv", qvalThresh, foldThresh)
    makeMergedTables(phosExp, outName+"AllPhosphoProteinData.csv", qvalThresh, foldThresh, True)
    return

#################################
### Main Processing Functions ###
#################################

def loadPhosphoData(phosFile):
    phosData = pd.read_csv(phosFile,sep=",")
    uniprot = phosData["Defline"].str.split("|",expand=True)[1]
    phosData["Uniprot"] = uniprot
    phosData.name = phosFile
    return phosData

def loadProteinData(protFile):
    protData = pd.read_csv(protFile, sep=",")
    uniprot = protData["Representative Protein Description"].str.split("|",expand=True)[1]
    protData["Uniprot"] = uniprot
    protData.name = protFile
    return protData

"""
Uses quantile normalization to normalize data within a 10-plex
Data is the input data
normCols is a list of columns which should be included in the normalization
idCols is a list of columns which should be used to identify duplicate entries
"""
def normalize10Plex(data, normCols, idCols):
    data = removeDuplicates(data, idCols, normCols)
    #Quantile normalization code inspired by this thread:
    #https://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe

    datToNorm = data[normCols]
    rank_mean = datToNorm.stack().groupby(datToNorm.rank(method='first').stack().astype(int)).mean()
    normedData = datToNorm.rank(method='min').stack().astype(int).map(rank_mean).unstack()

    data.drop(labels=normCols, axis="columns", inplace=True)
    data[normCols] = normedData[normCols]
    return data

"""
Takes the means of replicates in 2 conditions and calculates fold changes for each.
Calculates things as cond2/cond1, so cond1 should be the earlier condition
"""
def calcFoldChanges(experiment):
    data = experiment[0]
    cond1 = experiment[1]
    cond2 = experiment[2]
    cond1Names = getColNames(data, cond1)
    cond2Names = getColNames(data, cond2)

    if ("mean_"+cond1) not in data:
        data["mean_"+cond1] = data[cond1Names].mean(axis=1)
    if ("mean_"+cond2) not in data:
        data["mean_"+cond2] = data[cond2Names].mean(axis=1)
    data["log_"+cond1+"_"+cond2] = np.log2(data["mean_"+cond2]/data["mean_"+cond1])
    return

"""
Writes single csv files for a set of t-tests that should be performed.
"""
def writeExpSetup(outName, experiment):
    data = experiment[0]
    cond1 = experiment[1]
    cond2 = experiment[2]
    name = "exp"+experiment[3]
    cond1Names = getColNames(data, cond1)
    cond2Names = getColNames(data, cond2)
    allCols = []
    for col in cond1Names:
        if col.startswith(cond1):
            allCols.append(col)
    for col in cond2Names:
        if col.startswith(cond2):
            allCols.append(col)
    data.to_csv(outName+name+"_"+cond1+"_"+cond2+".csv",columns=allCols,index=False)
    return

"""
Gets results of t-test and fdr correction from limma
"""
def getTestResults(outName, experiment):
    data = experiment[0]
    cond1 = experiment[1]
    cond2 = experiment[2]
    name = "exp"+experiment[3]

    fname = outName+name+"_"+cond1+"_"+cond2+".csv_tmpLimma"
    if not os.path.isfile(fname):
        vPrint("Couldn't find limma result file "+fname)
        return False

    raw_data = pd.read_csv(fname)
    data["pVal_"+cond1+"_"+cond2] = raw_data["pval"].tolist()
    data["qVal_"+cond1+"_"+cond2] = raw_data["qval"].tolist()

    #Try to cleanup tmp files
    try:
        os.remove(fname)
    except OSError:
        pass

    return True

"""
Filters out significantly changed entities
"""
def filterRes(experiment,qvalThresh,foldThresh):
    data = experiment[0]
    cond1 = experiment[1]
    cond2 = experiment[2]
    name = experiment[3]

    filteredData = None
    foldThresh = np.log2(foldThresh)
    qCol = "qVal_"+cond1+"_"+cond2
    foldCol = "log_"+cond1+"_"+cond2

    filteredData = data.loc[(data[qCol] <= qvalThresh) & (np.abs(data[foldCol]) >= foldThresh)]
    vPrint("There are "+str(len(filteredData))+" significant results for "+" ".join(experiment[1:]))

    return filteredData

"""
Given lists of significant items, creates and save a prize list
"""
def createPrizeList(prot1, prot2, phos1, phos2, time, outLoc):
    timeN = "qVal"+time

    #First, we want to condense this down to a single list
    datList = [prot1[["Uniprot",timeN]],prot2[["Uniprot",timeN]],phos1[["Uniprot",timeN]],phos2[["Uniprot",timeN]]]
    sigDat = datList[0]
    for i in range(1,len(datList)):
        sigDat = sigDat.append(datList)
    sigDat = sigDat.sort_values(by=timeN, ascending=True)
    sigDat = sigDat.drop_duplicates(subset="Uniprot",keep="first")

    #Prizes are created from -log2(qVal)
    sigDat["Prize"] = -1 * np.log2(sigDat[timeN])

    #Save prize file
    sigDat.to_csv(outLoc + "prize_"+time+".csv",sep='\t',index=False, columns=["Uniprot","Prize"])
    return

"""
Makes input files for phosFate enrichment
"""
def makePhosfateInput(data,name,outName,time):
    #Get isoform as a raw number
    data = data.copy()
    data["Isoform"] = data["Isoform"].str.split(";")
    data = data.explode("Isoform")
    data = data.dropna()
    data["IsoNum"] = data["Isoform"].str.extract(".*\((.*)\).*")
    data.to_csv(outName+name+"Phosfate.csv",sep=",",columns=["Uniprot","IsoNum","log"+time],index=False,header=False)
    return

"""
Creates a single table of all protein or phospho data, determined by hasIso and the input data
"""
def makeMergedTables(experiments, fName, qvalThresh, foldThresh, hasIso=False):
    allData = []
    foldThresh = np.log2(foldThresh)

    #Rename columns for final output
    outCols = []
    idCols = ["Protein Group Name","Uniprot IDs","Representative Protein Description","Uniprot"]
    if hasIso:
        idCols = ["Protein Group","Defline","Isoform","Uniprot"]
    for experiment in experiments:
        data = experiment[0]
        cond1 = experiment[1]
        cond2 = experiment[2]
        name = experiment[3]

        colPref = cond1+" vs. "+cond2
        outColNames = ["Protein Group Name","Uniprot IDs","Representative Protein Description","Uniprot"]
        if hasIso:
            outColNames = ["Protein Group","Defline","Isoform","Uniprot"]
        outColNames += [colPref+" log2 fold change",colPref+" q-value",colPref+" significant"]
        data.rename(columns={"log_"+cond1+"_"+cond2:outColNames[-3], "qVal_"+cond1+"_"+cond2:outColNames[-2]}, inplace=True)

        data[colPref+" significant"] = (data[outColNames[-2]] <= qvalThresh) & (np.abs(data[outColNames[-3]]) >= foldThresh)
        allData.append(data[outColNames])
        outCols.append(outColNames)

    finalData = allData[0]
    for i in range(1,len(allData)):
        finalData = pd.merge(finalData, allData[i][outCols[i]],how="outer",on=idCols)
    finalData.to_csv(fName,sep=",",index=False,header=True)
    return

########################
### Helper Functions ###
########################

"""
Renames columns so they are standard across all data
"""
def stdColNames(data, colNames):

    repDict = {}

    for raw_colName in colNames:
        cName = raw_colName.lower()
        condition = ""
        if "env" in cName:
            condition = "env"
        elif "genome" in cName:
            condition = "genome"
        elif "wt" in cName:
            condition = "wt"

        time = 0
        if "60min" in cName:
            time = "60min"
        elif "5min" in cName:
            time = "5min"
        elif "0min" in cName:
            time = "0min"

        rep = 0
        if (condition, time) not in repDict:
            repDict[(condition, time)] = 1
            rep = 1
        else:
            repDict[(condition, time)] += 1
            rep = repDict[(condition, time)]
        data.rename(columns = {raw_colName: "_".join([condition, time, str(rep)])}, inplace = True)
    return

"""
Takes the larger of any identical rows and replaces them within the dataframe.
"""
def removeDuplicates(data,dupCols,normCols):

    data['tmpSum'] = data[normCols].sum(axis=1)

    dupInd = data.duplicated(subset=dupCols,keep=False)
    vPrint("Detected "+str(sum(dupInd)/2)+" duplicates")

    data = data.sort_values(by = 'tmpSum', ascending=False)
    data = data.drop_duplicates(subset=dupCols,keep="first")
    return data

"""
Gets a list of columns names for a certain timepoint
"""
def getColNames(data, condition):
    colNames = data.columns[data.columns.str.contains(pat = condition)]
    return colNames

"""
Function to print only if verbose flag is true
"""
def vPrint(pStr):
    global verbose
    if verbose:
        print(pStr)
    return

"""
Removes all non-intersecting items in dat1 and dat2
Currently not used
"""
def removeNonOverlap(dat1,dat2,idCols):
    ind1 = []
    for col in idCols:
        newInd = dat1[col].isin(dat2[col])
        if len(ind1) == 0:
            ind1 = newInd
        else:
            ind1 = ind1 & newInd

    ind2 = []
    for col in idCols:
        newInd = dat2[col].isin(dat1[col])
        if len(ind2) == 0:
            ind2 = newInd
        else:
            ind2 = ind2 & newInd

    dat2 = dat2[ind2]
    dat1 = dat1[ind1]

    return dat1,dat2
main()
