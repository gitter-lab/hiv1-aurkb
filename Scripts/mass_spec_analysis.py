import argparse
import pandas as pd
import numpy as np
import os
import os.path

"""
Author: Chris Magnano
8/17/2018

This script takes in raw data for both experiments, normalizes it, calculates fold changes, and uses limma in an R script to get significance values.

Outputs:
2. Files with all fold changes and p-values - 6 files, 3 10-plexes, and protein by phos data.
3. Files with only significantly changed values - 8 files, as we seperate the first 10-plex by timepoint.

Future improvements:
    -If we made this take in a file with lists of protein and phosphoprotein file names, and a file with column names for each, this would be fully generic.
    -Right now the script has to be run twice, and some calcualtion is re-done each time. That's pretty ugly (but hey, bibtex does it). Maybe we call R from python? Maybe call it with subprocess?
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

    #Load raw data. If future data has different column names, here's the place to that.
    vPrint("Loading Protein Data")
    protDat1 = loadProteinData(proteinFile1)
    prot1ColNames = ['0minRep1 (TMT 126C PCN)','0minRep2 (TMT 127N PCN)','0minRep3 (TMT 127C PCN)','0minRep4 (TMT 128N PCN)','5minRep1 (TMT 128C PCN)','5minRep2 (TMT 129N PCN)',
    '5minRep3 (TMT 129C PCN)','60minRep1 (TMT 130N PCN)','60minRep2 (TMT 130C PCN)','60minRep3 (TMT 131N PCN)']
    prot1IDColNames = ['Uniprot']

    protDat2 = loadProteinData(proteinFile2)
    prot2ColNames = ['WT_0Min_(TMT 126C PCN)','WT_0Min_(TMT 127N PCN)','WT_5Min_(TMT 128C PCN)','WT_5Min_(TMT 129N PCN)']
    prot2IDColNames = ['Uniprot']

    protDat3 = loadProteinData(proteinFile3)
    prot3ColNames = ['WT_0Min_(TMT 126C PCN)','WT_0Min_(TMT 127N PCN)','WT_60Min_(TMT 128C PCN)','WT_60Min_(TMT 128C PCN).1']
    prot3IDColNames = ['Uniprot']

    vPrint("Loading Phospho Data")
    phosDat1 = loadPhosphoData(phosphoFile1)
    phos1ColNames = ['0minRep1 (TMT 126C PCN)','0minRep2 (TMT 127N PCN)','0minRep3 (TMT 127C PCN)','0minRep4 (TMT 128N PCN)','5minRep1 (TMT 128C PCN)','5minRep2 (TMT 129N PCN)',
    '5minRep3 (TMT 129C PCN)','60minRep1 (TMT 130N PCN)','60minRep2 (TMT 130C PCN)','60minRep3 (TMT 131N PCN)']
    phos1IDColNames = ['Uniprot','Isoform']

    phosDat2 = loadPhosphoData(phosphoFile2)
    phos2ColNames = ['WT_0Min_(TMT 126C PCN)','WT_0Min_(TMT 127N PCN)','WT_5Min_(TMT 128C PCN)','WT_5Min_(TMT 129N PCN)']
    phos2IDColNames = ['Uniprot','Isoform']

    phosDat3 = loadPhosphoData(phosphoFile3)
    phos3ColNames = ['WT_0Min_(TMT 126C PCN)','WT_0Min_(TMT 127N PCN)','WT_60Min_(TMT 128C PCN)','WT_60Min_(TMT 128C PCN).1']
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
    stdColNames(protDat1,prot1ColNames[0:4],min5Names=prot1ColNames[4:7],min60Names=prot1ColNames[7:])
    stdColNames(protDat2,prot2ColNames[0:2],min5Names=prot2ColNames[2:])
    stdColNames(protDat3,prot3ColNames[0:2],min60Names=prot3ColNames[2:])
    stdColNames(phosDat1,phos1ColNames[0:4],min5Names=phos1ColNames[4:7],min60Names=phos1ColNames[7:])
    stdColNames(phosDat2,phos2ColNames[0:2],min5Names=phos2ColNames[2:])
    stdColNames(phosDat3,phos3ColNames[0:2],min60Names=phos3ColNames[2:])

    #Calculate fold changes
    vPrint("Calculating Fold Changes")
    calcFoldChanges(protDat1)
    calcFoldChanges(protDat2)
    calcFoldChanges(protDat3)
    calcFoldChanges(phosDat1)
    calcFoldChanges(phosDat2)
    calcFoldChanges(phosDat3)

    #Get test results if available
    gotRes = True
    gotRes = gotRes and getTestResults(protDat1,outName,"expProt1")
    gotRes = gotRes and getTestResults(protDat2,outName,"expProt2")
    gotRes = gotRes and getTestResults(protDat3,outName,"expProt3")
    gotRes = gotRes and getTestResults(phosDat1,outName,"expPhos1")
    gotRes = gotRes and getTestResults(phosDat2,outName,"expPhos2")
    gotRes = gotRes and getTestResults(phosDat3,outName,"expPhos3")

    #If we don't have results yet, write out single experiment files
    if not gotRes:
        writeExpSetup(protDat1,outName,"expProt1")
        writeExpSetup(protDat2,outName,"expProt2")
        writeExpSetup(protDat3,outName,"expProt3")
        writeExpSetup(phosDat1,outName,"expPhos1")
        writeExpSetup(phosDat2,outName,"expPhos2")
        writeExpSetup(phosDat3,outName,"expPhos3")
        vPrint("Run limma-test.R with same output file name argument for each experiment, then run this script again.")
        return

    #Save unfiltered data
    vPrint("Saving unfiltered data")
    protDat1.to_csv(outName+"unfilteredProt1.tsv",sep='\t',index=False)
    protDat2.to_csv(outName+"unfilteredProt2.tsv",sep='\t',index=False)
    protDat3.to_csv(outName+"unfilteredProt3.tsv",sep='\t',index=False)
    phosDat1.to_csv(outName+"unfilteredPhos1.tsv",sep='\t',index=False)
    phosDat2.to_csv(outName+"unfilteredPhos2.tsv",sep='\t',index=False)
    phosDat3.to_csv(outName+"unfilteredPhos3.tsv",sep='\t',index=False)

    #Save merged tables
    vPrint("Saving merged tables")
    makeMergedTables(protDat1, protDat2, protDat3, outName+"AllProteinData.csv")
    makeMergedTables(phosDat1, phosDat2, phosDat3, outName+"AllPhosphoProteinData.csv",True)

    #Create phosfate input
    vPrint("Saving phosfate input")
    makePhosfateInput(phosDat1, "phos1_05_", outName ,"05")
    makePhosfateInput(phosDat1, "phos1_060_", outName ,"060")
    makePhosfateInput(phosDat2, "phos2_05_", outName ,"05")
    makePhosfateInput(phosDat3, "phos3_060_", outName ,"060")

    #Filter out significant data
    vPrint("Filtering Significant Results")
    prot1_05, prot1_60 = filterRes(protDat1,qvalThresh,foldThresh)
    prot2_05, prot2_60 = filterRes(protDat2,qvalThresh,foldThresh)
    prot3_05, prot3_60 = filterRes(protDat3,qvalThresh,foldThresh)
    phos1_05, phos1_60 = filterRes(phosDat1,qvalThresh,foldThresh)
    phos2_05, phos2_60 = filterRes(phosDat2,qvalThresh,foldThresh)
    phos3_05, phos3_60 = filterRes(phosDat3,qvalThresh,foldThresh)

    #Save filtered data
    vPrint("Saving filtered data")
    prot1_05.to_csv(outName+'sigProt1_05.tsv',sep='\t',index=False)
    prot1_60.to_csv(outName+'sigProt1_60.tsv',sep='\t',index=False)
    prot2_05.to_csv(outName+'sigProt2_05.tsv',sep='\t',index=False)
    prot3_60.to_csv(outName+'sigProt3_60.tsv',sep='\t',index=False)
    phos1_05.to_csv(outName+'sigPhos1_05.tsv',sep='\t',index=False)
    phos1_60.to_csv(outName+'sigPhos1_60.tsv',sep='\t',index=False)
    phos2_05.to_csv(outName+'sigPhos2_05.tsv',sep='\t',index=False)
    phos3_60.to_csv(outName+'sigPhos3_60.tsv',sep='\t',index=False)

    #Create and save off prize lists
    if makePrizes:
        createPrizeList(prot1_05, prot2_05, phos1_05, phos2_05, "05", outName)
        createPrizeList(prot1_60, prot2_60, phos1_60, phos2_60, "060", outName)
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
Takes the means of replicates in the same timepoint and calculates fold changes for each.
At least 1 of min5Names min60Names should not be empty.
"""
def calcFoldChanges(data):
    min0Names = getColNames(data,0)
    min5Names = getColNames(data,5)
    min60Names = getColNames(data,60)

    data["mean0min"] = data[min0Names].mean(axis=1)
    if len(min5Names) > 0:
        data["mean5min"] = data[min5Names].mean(axis=1)
        data["log05"] = np.log2(data["mean5min"]/data["mean0min"])
    if len(min60Names) > 0:
        data["mean60min"] = data[min60Names].mean(axis=1)
        data["log060"] = np.log2(data["mean60min"]/data["mean0min"])
    return

"""
Writes single csv files for each set of t-tests that should be performed.
"""
def writeExpSetup(data,outName,name):
    min0Names = getColNames(data,0)
    min5Names = getColNames(data,5)
    min60Names = getColNames(data,60)

    if len(min5Names)>0:
        data.to_csv(outName+name+"_05.csv",columns=min0Names+min5Names,index=False)
    if len(min60Names)>0:
        data.to_csv(outName+name+"_060.csv",columns=min0Names+min60Names,index=False)
    return

"""
Gets resuls of t-test and fdr correction from limma
"""
def getTestResults(data,outName,name):
    min5Names = getColNames(data,5)
    min60Names = getColNames(data,60)

    if len(min5Names)>0:
        fname = outName+name+"_tmpLimma_05.csv"
        if not os.path.isfile(fname):
            vPrint("Couldn't find limma result file "+fname)
            return False

        dat5Min = pd.read_csv(fname)
        data["pVal05"] = dat5Min["pval"].tolist()
        data["qVal05"] = dat5Min["qval"].tolist()

        #Try to cleanup tmp files
        try:
            os.remove(fname)
            os.remove(outName+name+"_05.csv")
        except OSError:
            pass

    if len(min60Names)>0:
        fname = outName+name+"_tmpLimma_060.csv"
        if not os.path.isfile(fname):
            vPrint("Couldn't find limma result file "+fname)
            return False
        dat60Min = pd.read_csv(fname)
        data["pVal060"] = dat60Min["pval"].tolist()
        data["qVal060"] = dat60Min["qval"].tolist()

        #Try to cleanup tmp files
        try:
            os.remove(fname)
            os.remove(outName+name+"_060.csv")
        except OSError:
            pass
    return True

"""
Filters out significantly changed entities
"""
def filterRes(data,qvalThresh,foldThresh):
    filtered05 = None
    filtered060 = None
    foldThresh = np.log2(foldThresh)

    if "log05" in data:
        filtered05 = data.loc[(data.qVal05 <= qvalThresh) & (np.abs(data.log05) >= foldThresh)]
        vPrint("There are "+str(len(filtered05))+" significant results at 5 min in this file.")

    if "log060" in data:
        filtered060 = data.loc[(data.qVal060 <= qvalThresh) & (np.abs(data.log060) >= foldThresh)]
        vPrint("There are "+str(len(filtered060))+" significant results at 60 min in this file.")

    return filtered05,filtered060

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
    data["IsoNum"] = data["Isoform"].str.extract(".*\((.*)\).*")
    data.to_csv(outName+name+"Phosfate.csv",sep=",",columns=["Uniprot","IsoNum","log"+time],index=False,header=False)
    return

"""
Creates a single table of all protein or phospho data, determined by hasIso and the input data
"""
def makeMergedTables(plex1, plex2, plex3, fName, hasIso=False):
    #Rename columns for final output
    allData = [plex1, plex1, plex2, plex3]
    experimentNums = ["1","1","2","2"]
    timePoints = ["5","60","5","60"]
    outCols = []
    idCols = ["Uniprot"]
    if hasIso:
        idCols.append("Isoform")

    for i in range(len(allData)):
        exp = experimentNums[i]
        time = timePoints[i]
        colPref = "Experiment "+exp+" "+time+" min "
        outColNames = ["Uniprot"]
        if hasIso:
            outColNames.append("Isoform")
        outColNames += [colPref+"fold change",colPref+"q-value"]
        allData[i] = allData[i].rename(columns={"log0"+time:outColNames[-2], "qVal0"+time:outColNames[-1]})
        allData[i] = allData[i][outColNames]
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
Renames columns so they are standard accross all data
"""
def stdColNames(data,min0Names,min5Names=[],min60Names=[]):
    for i in range(len(min0Names)):
        data.rename(columns={min0Names[i]:"0min_"+str(i+1)},inplace=True)
    if len(min5Names)>0:
        for i in range(len(min5Names)):
            data.rename(columns={min5Names[i]:"5min_"+str(i+1)},inplace=True)
    if len(min60Names)>0:
        for i in range(len(min60Names)):
            data.rename(columns={min60Names[i]:"60min_"+str(i+1)},inplace=True)
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
Get's a list of columns names for a certain timepoint
"""
def getColNames(data,time):
    colNames = []
    hasCol = True
    i=1
    while hasCol:
        name = str(time)+"min_"+str(i)
        if name in data:
            colNames.append(name)
            i+=1
        else:
            hasCol=False
    return colNames

"""
Function to print only if verbose flag is true
"""
def vPrint(pStr):
    global verbose
    if verbose:
        print pStr
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
