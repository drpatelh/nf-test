from turtle import color
import pandas as pd 
import datapane as dp
import reportGeneration.myconstants as constants
import plotly.express as px
from reportGeneration.ReportUtil import GeneralUtils, CalculationUtils, DatapaneUtils
from reportGeneration.AbstractReport import AbstractReport
from reportGeneration.SampleReport import SampleReport
from pathlib import Path
from typing import Dict, List, Tuple
from itertools import groupby

# Threshjson, filterDf (both in QC)

class FastqReport(AbstractReport):
    fastqName: str
    sampleNames: List[str]
    sampleFiles: List[Tuple[str,Dict[str,Path]]]
    sampleQCDf: pd.DataFrame
    sampleThresholds: Dict[str, Dict[str, float]]
    def __init__(self, fastqName: str, samplesheetFn: str, resultsDir: Path, makeAdvanced:bool, demuxJson: Dict[str,float], qcDir:Path, libDir:Path):
        super().__init__(resultsDir, makeAdvanced, demuxJson, qcDir, fastqName, libDir)
        self.sampleNames = sorted(FastqReport.getAssociatedSampleNames(samplesheetFn,fastqName))
        self.sampleFiles = self.resolveSampleSpecificFilePaths()
        self.sampleThresholds = self.getSampleThresholds()
        self.fastqName = fastqName
        self.sampleQCDf = self.buildMultiSampleQCDf()
        self.splitBarcodes()
        self.validateInput()

    def resolveSampleSpecificFilePaths(self):
        sampleToFilePaths = {}
        for sample in self.sampleNames:
            filePaths = {}
            if str(self.resultsDir) == '.':
                filePaths['filterDfPath'] = self.resultsDir / f'{sample}_QC.tsv'
                filePaths['threshJsonPath'] = self.resultsDir / f'{sample}_thresholds.json'
            else:
                QCPrefix = self.resultsDir / 'QC' / sample / self.qcSubdirectory
                filePaths['filterDfPath'] = QCPrefix / f'{sample}_QC.tsv'
                filePaths['threshJsonPath'] = QCPrefix / f'{sample}_thresholds.json'
            sampleToFilePaths[sample] = filePaths
        return sampleToFilePaths

    def validateInput(self):
        GeneralUtils.makeDir(self.writeDir)
        for (_, sampleFilePaths) in self.sampleFiles.items():
            GeneralUtils.ensurePathsExist(sampleFilePaths)

    def splitBarcodes(self):
        lvl1BarcodeAliases = {}
        for sampleName in self.sampleNames:
            for name,bc in self.demuxJson["samples"][sampleName]["barcodes"].items():
                lvl1BarcodeAliases[bc['sequence']] =  name 
        libDirJSON = GeneralUtils.readJSON(self.libDir, preserveDictOrder=True)
        barcodeInfo = libDirJSON['barcodes']
        scDefiningBarcodes = [barcodeMeta for barcodeMeta in barcodeInfo if barcodeMeta.get('type', None) != 'library_index']
        groupByLevel = groupby(scDefiningBarcodes, key=lambda x: x['level'])
        expectedBarcodeCount = len(scDefiningBarcodes)
        splitByBarcode = self.sampleQCDf.Barcode.str.split('+')
        observedBarcodeCount = len(splitByBarcode[0])
        if observedBarcodeCount == expectedBarcodeCount:
            curBcIndex = 0
            barcodesToPlot = {}
            for level, group in groupByLevel:
                groupAsList = list(group)
                groupSize = len(groupAsList)
                barcodeName = groupAsList[0]['alias']
                barcodes = [''.join(x[curBcIndex:(curBcIndex+groupSize)]) for x in splitByBarcode]
                self.sampleQCDf[barcodeName] = barcodes
                barcodesToPlot[level] = barcodeName
                curBcIndex += groupSize
                if level == 1 and len(lvl1BarcodeAliases) > 0:
                    self.sampleQCDf[f"{barcodeName}_alias"] = self.sampleQCDf[barcodeName].apply(lambda x: lvl1BarcodeAliases[x])
        else: 
             raise ValueError(f"The number of sequences in '+' split cell barcodes ({observedBarcodeCount})\ndoes not match the inferred number of barcodes ({expectedBarcodeCount}).\nEnsure the specified library definition is correct for this data")
        self.barcodesToPlot = barcodesToPlot

    def build(self):
        readsPage = self.buildReadsPage()
        cellsPage = self.buildCellsPage()
        pages = [readsPage, cellsPage]
        report = dp.Report(
            blocks=pages
        )
        report.save(self.writeDir / f"{self.fastqName}.report.html") 

    def buildMultiSampleQCDf(self) -> pd.DataFrame:
        '''
        Aggreagtes QC.tsvs for each sample into a single pandas dataframe 
        '''
        accumulatedDf = pd.DataFrame()
        for sample,sampleFilePaths in self.sampleFiles.items():
            threshJson = GeneralUtils.readJSON(sampleFilePaths['threshJsonPath'])
            (df,_) = CalculationUtils.getStatsAsDataframes(sampleFilePaths['filterDfPath'], None, threshJson) 
            df['sample'] = sample
            df.sort_values(by=['UniqueReads'], inplace=True, ascending=False)
            df['Barcode'] = df.index
            index = range(len(df.index))
            # Assigning index prevents overwritting during concatenation 
            df.index = index
            accumulatedDf = pd.concat([accumulatedDf, df])
        return accumulatedDf

    def getSampleThresholds(self) -> Dict[str, int]:
        results = {}
        for sample,sampleFilePaths in self.sampleFiles.items():
            threshJson = GeneralUtils.readJSON(sampleFilePaths['threshJsonPath'])
            results[sample] = threshJson
        return results
    
    ###################### READS PAGE FUNCTIONS ##########################
    def buildMultisampleCellStatsDf(self):
        selectBlocks = []
        for sample in self.sampleNames:
            sampleSubsetDf = self.sampleQCDf[self.sampleQCDf['sample'] == sample]
            cellStatsDataFrame = SampleReport.makeCellStatsDf(sampleSubsetDf, self.sampleThresholds[sample]['UniqueReads'])
            cellStatsDataFrameStyled = cellStatsDataFrame.style.pipe(GeneralUtils.styleTable, title="Reads", hideColumnHeaders=True, boldColumn=['Metric'])
            cellStatsAsDpTable = DatapaneUtils.createTableIgnoreWarning(cellStatsDataFrameStyled, label=sample)
            selectBlocks.append(cellStatsAsDpTable)
        return dp.Select(type=dp.SelectType.TABS,blocks=selectBlocks)

    def buildDfFromDemuxSampleMetrics(self):
        totalCounts = []
        tgmtCounts = []
        lvl1Barcode = self.barcodesToPlot[1]
        for sampleName, sampleDict in self.demuxJson['samples'].items():
            readCount = sampleDict['reads'][0]
            totalCounts.append({'Sample':sampleName, 'TotalReads':readCount})
            lvl1BarcodeCounts = sampleDict['barcodes']
            if (len(lvl1BarcodeCounts) > 0):
                for lvl1Bc, metrics in lvl1BarcodeCounts.items():
                    tgmtCounts.append({'Sample': sampleName, lvl1Barcode: lvl1Bc, 'ReadCount': metrics['reads']})
        return (pd.DataFrame(totalCounts), pd.DataFrame(tgmtCounts))

    def getCellsAboveThreshold(self) -> Dict[str, int]:
        results = {}
        for sample, thresholds in self.sampleThresholds.items():
            cellsPassed = len(self.sampleQCDf[(self.sampleQCDf['UniqueReads'] >= thresholds['UniqueReads']) & (self.sampleQCDf['sample'] == sample)].index)
            results[sample] = cellsPassed
        return results

    def makeMultisampleKneePlot(self, cellsAboveDict: Dict) -> dp.Plot:
        '''
        Makes a kneeplot using @field in @data; drawing a vertical line at @threshold
        '''  
        indiciesToInclude = set(CalculationUtils.getIndicesToInclude(max(self.sampleQCDf.index)))
        self.sampleQCDf['rankOrder'] = self.sampleQCDf.index
        plottingDf = self.sampleQCDf[self.sampleQCDf['rankOrder'].isin(indiciesToInclude)]
        fig = px.line(plottingDf, x=plottingDf.index, y='UniqueReads', color='sample', log_x=True, log_y=True, template=constants.DEFAULT_FIGURE_STYLE, labels={"index": "Cell Barcodes"}, title="Per Sample: Reads Per Cell Barcode")
        colorsBeingUsed = px.colors.qualitative.D3
        # Adding vertical lines at cells above threshold for each sample (in same color)
        for i, (_, cellCount) in enumerate(sorted(cellsAboveDict.items())):
            fig.add_vline(x=cellCount, line_dash="dash", line_color=colorsBeingUsed[i])
        return dp.Plot(fig)

    def buildReadsPage(self):
        barcodeTypeStats = self.createBarcodeTypeMetricsTables()
        barcodeReadsData = self.demuxJson['reads']
        barcodeReadsPerc = FastqReport.buildDfFromJSONDict(barcodeReadsData, "Type", "list")
        barcodeReadsTotal = FastqReport.buildDfFromJSONDict(barcodeReadsData, "Type", "list", 0)
        barcodeReadsTotal = barcodeReadsTotal[['Type', 'Reads']]
        barcodeReadsTotal.rename(columns={'Type':'Status'}, inplace=True)
        barcodeReadsTotal['Percent'] = barcodeReadsPerc['Reads']
        barcodeReadsTotal = barcodeReadsTotal.style.pipe(GeneralUtils.styleTable, title="Barcode Read Status", numericCols=['Reads'])    
        barcodeReadStats = DatapaneUtils.createTableIgnoreWarning(barcodeReadsTotal)
        (countsPerSampleDf, lvl1BcCountsPerSampleDf) = self.buildDfFromDemuxSampleMetrics()
        lvl1BcCountsPerSampleDf.sort_values(by=[self.barcodesToPlot[1],'Sample'], inplace=True)
        countsPerSampleDf.sort_values(by=['Sample'], inplace=True)
        colorMapToUse = FastqReport.matchColorsToNames(px.colors.qualitative.D3, list(countsPerSampleDf['Sample'].unique()))
        readsPerSample = px.bar(countsPerSampleDf, x='TotalReads', y='Sample', color='Sample',height=500, color_discrete_map=colorMapToUse,template=constants.DEFAULT_FIGURE_STYLE, title="Reads per sample")
        lvl1BcCountsPerSampleDf.sort_values(by=['Sample'], inplace=True)
        readsPerWellBySample = px.bar(lvl1BcCountsPerSampleDf, x='ReadCount', y=self.barcodesToPlot[1], color='Sample',height=500,template=constants.DEFAULT_FIGURE_STYLE, color_discrete_map=colorMapToUse, labels={self.barcodesToPlot[1]: f'{self.barcodesToPlot[1]} Barcode ', 'ReadCount':'Total Reads'}, title=f"Reads count by {self.barcodesToPlot[1]} Barcode")
        return dp.Page(
                blocks=[
                    dp.Group(barcodeReadStats,barcodeTypeStats, readsPerSample,readsPerWellBySample, columns=2),
                ],
                title='Reads'
        )    

    def matchColorsToNames(colorPalette, names) -> Dict[str,str]:
        colorMap = {"Unknown": 'rgb(179, 188, 201)'}
        if ("Unknown" in names):
            names.remove('Unknown')
        '''
        Associated colors with categorical values @names
        for consistency between plots 
        '''
        colorsBeingUsed = px.colors.qualitative.D3
        for i,name in enumerate(sorted(names)):
            colorMap[name] = colorsBeingUsed[i]
        return colorMap


    ###################### CELLS PAGE FUNCTIONS ##########################
    def buildMultisampleOverloadingFigure(self):
        lvl2_barcode = self.barcodesToPlot[2]
        cells = self.sampleQCDf[self.sampleQCDf['pass']]
        lvl2Cnts = cells.groupby([lvl2_barcode, 'sample']).size().reset_index()
        lvl2Cnts.rename(columns={0: 'BarcodesInBead'}, inplace=True)
        cells['BarcodesInBead'] = cells.groupby([lvl2_barcode, 'sample'])[lvl2_barcode].transform('size')
        cellsPerDropBar = px.box(cells, y='UniqueReads', x='BarcodesInBead', log_y=True, color='sample', title='Unique reads per cell ',template="none", labels={"BarcodesInBead": f"Cells Per {lvl2_barcode}", "UniqueReads":"Unique Reads"}, points=False)
        DatapaneUtils.showAllTicks(cellsPerDropBar)
        barcodesPerBeadBySample = lvl2Cnts.groupby(['BarcodesInBead', 'sample']).size().reset_index()
        barcodesPerBeadBySample.rename(columns={0: 'Counts'}, inplace=True)
        cellsPerDropLine = px.bar(barcodesPerBeadBySample, x='BarcodesInBead', y='Counts', color='sample',log_y=True,  template=constants.DEFAULT_FIGURE_STYLE, title="Cells per bead", labels={"BarcodesInBead": f"Cells Per {lvl2_barcode}"})
        DatapaneUtils.showAllTicks(cellsPerDropLine)
        return dp.Group(dp.Plot(cellsPerDropLine),dp.Plot(cellsPerDropBar), columns=2)
   
    def createBarcodeTypeMetricsTables(self):
        barcodes = self.demuxJson['barcodes']
        barcodesDf = FastqReport.buildDfFromJSONDict(barcodes, "Barcode", "dict")
        #BARCODE_SHORTHAND_TO_NAME
        tableGroup = []
        allBarcodeTypes = list(barcodesDf['Barcode'].unique())
        if 'index' in allBarcodeTypes:
            allBarcodeTypes.remove('index')
        for barcodeType in allBarcodeTypes:
            fullBarcodeTypeName = constants.BARCODE_SHORTHAND_TO_NAME[barcodeType]
            subset = barcodesDf[barcodesDf['Barcode'] == barcodeType][['Match', 'Reads']]
            styledDf = subset.style.pipe(GeneralUtils.styleTable, title=fullBarcodeTypeName)                       
            table = DatapaneUtils.createTableIgnoreWarning(styledDf, label=fullBarcodeTypeName)
            tableGroup.append(table)
        return dp.Select(blocks=tableGroup)


    def createReadsPerBarcodeFigures(self):
        self.sampleQCDf.index = self.sampleQCDf.Barcode
        lvl1_Barcode = self.barcodesToPlot[1]
        lvl2_Barcode = self.barcodesToPlot[2]
        cellBeads = set(self.sampleQCDf[self.sampleQCDf['pass']][lvl2_Barcode].unique())
        beads = {}
        for bead, cnts in self.sampleQCDf.groupby(lvl2_Barcode).UniqueReads: 
            beads[bead] = {}
            beads[bead]['cell'] = bead in cellBeads
            vals = cnts.sort_values(ascending=False) 
            beads[bead]['top'] = vals[0]
            beads[bead]['total'] = sum(vals)
            beads[bead]['count'] = len(vals)
        beadStats = pd.DataFrame.from_dict(beads, orient='index')
        beadStats.sort_values(by="total", ascending=False, inplace=True)
        # Fraction of total unique reads originating from the tgmt-droplet barcode in a droplet 
        beadStats['frac'] = (beadStats.top+1) / (beadStats.total+1)
        readsPerBarcodeKnee = px.line(x=range(len(beadStats.total)), y=beadStats.total, title=f'Reads per {lvl2_Barcode} Barcode', log_x=True, log_y=True, template=constants.DEFAULT_FIGURE_STYLE, labels={"x":f"{lvl2_Barcode} Barcodes", "y":"Unique Reads"})
        readsPerBarcodeKnee.add_vline(x=beadStats.cell.sum(), line_dash="dash", line_color="green")
        readsPerBarcodeKnee.add_annotation(x=0.5,y=0, showarrow=False, text=f"{beadStats.cell.sum()} Beads", font=dict(
                        size=14,
                        color="Green"
        ))
        ReadsFromTopBcHist = px.histogram(beadStats.frac[beadStats.cell], template="none",labels={"value":"Fraction reads from top barcode"}, title=f"Mean Reads from top {lvl1_Barcode} Barcode: {beadStats.frac[beadStats.cell].mean():.1%}")
        ReadsFromTopBcHist.update_layout(showlegend=False)
        return dp.Group(blocks=[dp.Plot(readsPerBarcodeKnee), dp.Plot(ReadsFromTopBcHist)], columns=2)


    def buildCellsPage(self):    
        sampleCellsAboveThreshold = self.getCellsAboveThreshold()
        multiSampleKneePlot = self.makeMultisampleKneePlot(sampleCellsAboveThreshold)
        blocks=[dp.Text(f"## FastqName: {self.fastqName}"),multiSampleKneePlot, self.buildMultisampleOverloadingFigure()]
        if self.makeAdvanced:
            blocks.append(self.createReadsPerBarcodeFigures())
        return dp.Page(
                blocks=blocks,
                title='Cells'
    )

    ###################### STATIC HELPERS ##########################
    @staticmethod
    def getAssociatedSampleNames(sampleSheetFn: str, fastqName: str): 
        '''
        Returns a list of the values in the 'sample' column of @sampleSheetFn 
        for rows where 'fastqName' equals @fastqName
        '''
        sampleSheetDf = pd.read_csv(sampleSheetFn)
        associatedSampleMetaData = sampleSheetDf[sampleSheetDf['fastqName'] == fastqName]
        associatedSampleNames = list(associatedSampleMetaData['sample'])
        return associatedSampleNames

    @staticmethod
    def buildDfFromJSONDict(jsonDict: Dict, name: str, valueDataType: str, choiceIndex=1) -> pd.DataFrame:
        """
        Builds a dataframe from a json containing metrics  
        Used as helper to process demux metrics.json 
        """
        dictList = []
        for (keyName,valueObj) in jsonDict.items():
            if valueDataType == 'dict':
                for key,value in valueObj.items():
                    newDict={}
                    newDict[name]=keyName
                    newDict['Match']=key
                    newDict['Reads']=value[choiceIndex]
                    dictList.append(newDict)
            elif valueDataType == 'list':
                newDict = {"Reads": valueObj[choiceIndex]}
                newDict[name] = keyName
                dictList.append(newDict)
            elif valueDataType == 'int':
                newDict = {"Cells": valueObj}
                newDict[name] = keyName
                dictList.append(newDict)
        return pd.DataFrame(dictList)


