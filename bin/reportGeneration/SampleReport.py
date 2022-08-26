#!/usr/bin/env python
import numpy as np
import pandas as pd 
import datapane as dp
import reportGeneration.myconstants as constants
import plotly.express as px
import plotly.figure_factory as ff
from itertools import combinations, groupby
from reportGeneration.ReportUtil import GeneralUtils, CalculationUtils, DatapaneUtils
from reportGeneration.AbstractReport import AbstractReport
from pathlib import Path
from typing import Dict, Union, List, Tuple
from math import log10


class SampleReport(AbstractReport):
    sampleName: str
    sampleSpecificFilePaths: Dict[str,Path]
    sampleThresholds: Dict[str, float]
    unfilteredCells: pd.DataFrame
    cellStats: pd.DataFrame
    speciesList: List[str]
    

    def __init__(self, sampleName: str, resultsDir: Path, makeAdvanced:bool, demuxJson: Dict[str,float], qcDir:Path, libDir:Path):
        super().__init__(resultsDir, makeAdvanced, demuxJson, qcDir, sampleName, libDir)
        self.sampleName = sampleName
        self.sampleSpecificFilePaths = self.resolveSampleSpecificFilePaths()
        self.validateInput()
        self.sampleThresholds = GeneralUtils.readJSON(self.sampleSpecificFilePaths['threshJsonPath'])
        
    def validateInput(self):
        GeneralUtils.makeDir(self.writeDir)
        GeneralUtils.ensurePathsExist(self.sampleSpecificFilePaths)

    def build(self):
        (df, cellStats) = CalculationUtils.getStatsAsDataframes(self.sampleSpecificFilePaths['filterDfPath'], self.sampleSpecificFilePaths['cellStatsPath'], self.sampleThresholds)    
        self.unfilteredCells = df
        self.cellStats = cellStats
        self.splitBarcodes()
        self.getSpeciesList()
        dedupStats = GeneralUtils.loadStatsAsDict(self.sampleSpecificFilePaths['dedupStatsPath'])
        (summaryPage, cellAndReadStats) = self.buildSummaryPage(dedupStats)
        cellsPage = self.buildCellsPage()
        barcodesPage = self.buildBarcodesPage()
        (fragmentsPage, fragmentStats) = self.buildFragmentsPage(dedupStats)
        pages = [summaryPage, cellsPage, barcodesPage,fragmentsPage]
        if len(self.speciesList) > 1: 
            pages.append(self.buildBarnyardPage())
            #accumulatedStats.append(barnyardPage)
        (peaksPage, peakStats) = self.buildPeaksPage() 
        accumulatedStats = [cellAndReadStats, fragmentStats, peakStats]
        pages.append(peaksPage)
        statsDf = pd.concat(accumulatedStats)
        statsDf.to_csv(self.writeDir/ f"{self.sampleName}.reportStatistics.tsv", index=False, sep='\t',header=False)
        report = dp.Report(
            blocks=pages
        )

        report.save(self.writeDir / f"{self.sampleName}.report.html")
    
    ###################### SUMMARY PAGE FUNCTIONS ##########################
    @staticmethod
    def makeCellStatsDf(cellStats: pd.DataFrame, cellthres:int)-> pd.DataFrame:
        '''
        Calculates cell related statistics of ATAC run, writing them to a DataFrame
        to be displayed in the report
        '''
        cells = cellStats[cellStats['Filter'] != "LowUniqueReads"]
        stats = []
        stats.append(["Cells above threshold", len(cells.index)])
        stats.append(["Cell threshold (Unique Reads)", cellthres])
        stats.append(["Median unique reads per cell", int(round(cells.UniqueReads.median()))])
        stats.append(["Median saturation in cells", str(round(1-(cells.UniqueReads / cells.PassingReads).median(),3))])
        stats.append(["Extrapolated unique reads @ 20k total reads", CalculationUtils.extrapolate_umis_at_nreads(cells.TotalReads.median(), cells.UniqueReads.median(), 20000)])
        statsAsDictList = [GeneralUtils.makeTableDict(['Metric','Value'], valuePair) for valuePair in stats]
        statsDataFrame =  pd.DataFrame(statsAsDictList)
        statsDataFrame['Category'] = "Cells"
        return statsDataFrame 

    def makeReadsStatsDf(self,cellStats: pd.DataFrame, dedupStats: Dict, cellthres: int)-> pd.DataFrame:
        '''
        Calculates read related statistics of ATAC run, writing them to a DataFrame
        to be displayed in the report
        '''
        cellStats['pass'] = cellStats['PassingReads'] >= cellthres
        cells = cellStats[cellStats['pass']]
        if len(cells.index) < 1:
            raise ValueError(f"No cells surpass the given Unique Read Threshold: {cellthres} ")
        stats = []
        stats.append(["Total Reads", dedupStats['Input Reads']])
        stats.append(["Mean Total Reads Per Cell", int(round(cells['TotalReads'].mean()))])
        stats.append(["Aligned Passing Reads", f"{dedupStats['Passing Reads'] / dedupStats['Input Reads']:.1%}"])
        stats.append(["Reads in cells",f"{cells.PassingReads.sum() / cellStats.PassingReads.sum():.1%}"])
        stats.append(["Mito. Reads", f"{dedupStats['Filtered Mito. reads'] / dedupStats['Input Reads']:.1%}"]) 
        stats.append(["Duplicates", f"{1-dedupStats['Unique Reads'] / dedupStats['Passing Reads']:.1%}"])
        statsAsDictList = [GeneralUtils.makeTableDict(['Metric','Value'], valuePair) for valuePair in stats]
        statsDataFrame = pd.DataFrame(statsAsDictList)
        statsDataFrame['Category'] = "Reads"
        return statsDataFrame

    def buildSummaryPage(self, dedupStats: Dict) -> dp.Page:
        '''
        , , self.sampleName
        Returns df.Page containing figures for summary page
        '''
        
        kneePlot = DatapaneUtils.makeKneePlot(self.unfilteredCells, 'UniqueReads', "Unique Reads Per Cell Barcode")
        cellStatsDf = SampleReport.makeCellStatsDf(self.unfilteredCells, self.sampleThresholds['UniqueReads'])
        cellStatsTable = DatapaneUtils.createTableIgnoreWarning(cellStatsDf[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, title="Cells", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value']))
        saturationPlot = DatapaneUtils.makeScatterPlot(self.unfilteredCells, "TotalReads", "Saturation", "Filter", "Saturation", True, False, constants.QC_COLORMAP, int(round(self.sampleThresholds['UniqueReads']*0.1)))
        readStatsDf = self.makeReadsStatsDf(self.unfilteredCells, dedupStats, self.sampleThresholds['UniqueReads'])
        readStatsTable = DatapaneUtils.createTableIgnoreWarning(readStatsDf[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, title="Reads", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value']))
        pageBlocks = [kneePlot, cellStatsTable, saturationPlot, readStatsTable]
        statsDf = pd.concat([cellStatsDf, readStatsDf])
        return (dp.Page(
                dp.Text(f"## Sample: {self.sampleName}"),
                dp.Group(blocks=pageBlocks, columns=2),
                title="Summary"
        ), statsDf)

     ###################### CELL PAGE FUNCTIONS ##########################
    def makeQCCategoryCountsTable(self, dataFrame:pd.DataFrame, cellsThresh:Dict[str, float]) -> dp.Table:
        '''
        Makes a datapane table displaying each filtering criteria, its 
        threshold and the number of cells belonging falling into that criteria 
        '''
        rows = []
        for criteria in cellsThresh.keys():
            cellsInCategory = len(dataFrame[dataFrame['Filter'] == f"Low{criteria}"].index)
            rows.append([f"Low{criteria}", f"{cellsInCategory:,}", f"{cellsThresh[criteria]:,}"])
        rows.append(['Pass', f"{len(dataFrame[dataFrame['Filter'] == 'Pass'].index):,}", ''])
        statsAsDictList = [GeneralUtils.makeTableDict(['Filter','Count', 'Cut Off'], valueSet) for valueSet in rows]
        counts = pd.DataFrame(statsAsDictList)
        counts = counts.style.pipe(GeneralUtils.styleTable, title="QC Filtering Categories", hideColumnHeaders=False)
        return DatapaneUtils.createTableIgnoreWarning(counts)
        
    def buildCellsPage(self) -> dp.Page:
        '''
        , constants.QC_COLORMAP, self.sampleThresholds
        Returns df.Page containing figures related to cell statistics
        '''
        groupBlocks = []
        columnCount = 1
        if (self.unfilteredCells['FractionMito'] > 0).any():
            mitoScatter = DatapaneUtils.makeScatterPlot(self.unfilteredCells, "TotalReads", "FractionMito", "Filter", "Mitochondrial Reads", True, False, constants.QC_COLORMAP, int(round( self.sampleThresholds['UniqueReads']*0.1)))
            groupBlocks.append(mitoScatter)
            columnCount = 2
        fripScatter = DatapaneUtils.makeScatterPlot(self.unfilteredCells, "UniqueReads", "FRiP", "Filter", "Fraction of Reads in Peaks", True, False, constants.QC_COLORMAP,int(round( self.sampleThresholds['UniqueReads']*0.1)))
        groupBlocks.append(fripScatter)
        countsByFilter = self.makeQCCategoryCountsTable(self.unfilteredCells,  self.sampleThresholds)
        return dp.Page(
            dp.Text(f"## Sample: {self.sampleName}"),
            dp.Group(blocks=groupBlocks, columns=columnCount),
            countsByFilter,
            title="Cells"
        )

    ###################### BARCODE PAGE FUNCTIONS #########################
    def createReadsPerBarcodeFigures(self):
        lvl1_Barcode = self.barcodesToPlot[1]
        lvl2_Barcode = self.barcodesToPlot[2]
        cellBeads = set(self.cellStats[self.cellStats['pass']][lvl2_Barcode].unique())
        beads = {}
        for bead, cnts in self.cellStats.groupby(lvl2_Barcode).UniqueReads: 
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
        return [dp.Plot(readsPerBarcodeKnee), dp.Plot(ReadsFromTopBcHist)]



    def createLvl1Figures(self, cellStats: pd.DataFrame, barcodesToPlot):
        lvl1_Barcode = barcodesToPlot[1]
        cellStats.sort_values(by=f"{lvl1_Barcode}_alias", inplace=True)
        readsPerLvl1Box = px.box(cellStats, y=f"{lvl1_Barcode}_alias", x='UniqueReads',height=700, log_x=True, color=f"{lvl1_Barcode}_alias", template="none",  labels={'UniqueReads':"Unique Reads Per Cell"}, notched=True, boxmode="overlay", points=False, title=f"Unique Read Counts by {lvl1_Barcode} Barcode")
        readsPerLvl1Box.update_yaxes(visible=False, showticklabels=False)
        readsPerLvl1Box.update_layout(showlegend=False)
        lvl1Cnts = cellStats.groupby(f"{lvl1_Barcode}_alias").size().reset_index()
        lvl1Cnts.rename(columns={0: "Number of Cells"},inplace=True)
        lvl1Cnts.sort_values(by=f"{lvl1_Barcode}_alias", inplace=True)
        cellsPerLvl1Bar = px.bar(lvl1Cnts, y=f"{lvl1_Barcode}_alias", x="Number of Cells", height=700,template="none", color=f"{lvl1_Barcode}_alias", title=f"Cell Counts by {lvl1_Barcode} Barcode", labels={f"{lvl1_Barcode}_alias": f"{lvl1_Barcode} Barcode"})
        cellsPerLvl1Bar.update_layout(showlegend=False)
        return [dp.Plot(cellsPerLvl1Bar), dp.Plot(readsPerLvl1Box)] 


    def createBarcodeFigures(self) -> dp.Group:
        '''
        Returns a dp.Group containing: 
        1. dp.Plots created in createPerWellFigures
        2. dp.Plots created in createPerDropletBarcodeFigures (if -advanced specified on CLI)
        '''
        blocks = []
        barcodeNames = {}
        for name,bc in self.demuxJson["samples"][self.sampleName]["barcodes"].items():
            barcodeNames[bc['sequence']] =  name 
        cells = self.cellStats[self.cellStats['pass']]
        isOverloaded = len(self.barcodesToPlot) > 1
        if isOverloaded:
            blocks = self.createLvl1Figures(cells, self.barcodesToPlot)
        else: 
            return dp.Text("## Not more than one barcode by which cells are multiplexed")
        return dp.Group(blocks=blocks, columns=2)


    def buildBarcodesPage(self) -> dp.Page: 
        """
        Returns df.Page containing figures related to Bead & Barcode statistics  
        """
        barcodeFigures= self.createBarcodeFigures()
        return dp.Page(
                blocks=[
                    dp.Text(f"## Sample: {self.sampleName}"), barcodeFigures
                ],
                title="Barcodes"
    )

    ###################### FRAGMENTS PAGE FUNCTIONS ##########################
    def makeFragStatDf(self,hist: pd.Series, dedupStats: Dict) -> pd.DataFrame:
        '''
        Calculates fragments related statistics of ATAC run, writing them to a DataFrame
        to be displayed in the report
        '''
        mode = hist.argmax()
        cdf = np.cumsum(hist)
        stats = []
        stats.append(["Reads aligned as pairs", f"{dedupStats['Mapped Paired'] / dedupStats['Passing Reads']:.1%}"])
        stats.append(["Mode (in bp)", f"{mode}"])
        stats.append(["Median (in bp)", f"{np.searchsorted(cdf, cdf.iloc[-1] // 2)}"])
        stats.append(["< 147 bp", f"{cdf[147] / cdf.iloc[-1]:.1%}"])
        stats.append(["147 - 294 bp", f"{(cdf[294] - cdf[147]) / cdf.iloc[-1]:.1%}"])
        stats.append(["> 294 bp", f"{(cdf.iloc[-1] - cdf[294]) / cdf.iloc[-1]:.1%}"])
        statsAsDictList = [GeneralUtils.makeTableDict(['Metric','Value'], valuePair) for valuePair in stats]
        statsDataFrame = pd.DataFrame(statsAsDictList)
        statsDataFrame['Category'] = "Fragments"
        return statsDataFrame
    
    def makeFragSizePlot(self, fragSizeHistFn: str) -> Tuple[pd.Series, dp.Plot]:
        '''
        Processes a string containing tab delimited values representing fragment size distribution
        To produce a histogram and line plot of the distribution 
        '''
        hist = {}
        for l in open(fragSizeHistFn):
            line = l.strip().split("\t")
            hist[int(line[0])] = int(line[1])
        hist = pd.Series(hist)
        hist = hist.rolling(10, min_periods=1, win_type='gaussian').mean(std=5)
        return (hist, dp.Plot(px.line(x=hist.index[0::5], y=hist.values[0::5], template="none", labels={'x':"Fragment Size", 'y':'# of Fragments'}, title='Fragment Length Distribution')))

    def buildFragmentsPage(self, dedupStats: Dict) -> dp.Page:
        """
        Build df.Page: "Fragment" 
        containing fragment size distribution & related stats in a table
        """
        (hist, fragSizePlot) = self.makeFragSizePlot(self.sampleSpecificFilePaths['fragmentHistPath'])
        fragStatDF = self.makeFragStatDf(hist, dedupStats)
        fragStatTable = DatapaneUtils.createTableIgnoreWarning(fragStatDF[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, title="", hideColumnHeaders=True, boldColumn=['Metric']))
        return (dp.Page(
                blocks=[
                    dp.Text(f"## Sample: {self.sampleName}"),
                    dp.Group(fragStatTable,fragSizePlot,columns=2)
                ],
                title="Fragments"
        ), fragStatDF)

    ###################### BARNYARD PAGE FUNCTIONS ##########################
    def makeOverloadBarnyardPlots(self) -> dp.Group:
        '''
        Creates a dp.Group consisting of 3 figures: 
        1. backgroundPercPerWell: Box plot seperated by tgmt barcode denoting % background for each 
        2. fracMixedCellsPerWell: Bar graph showing % mixed cells for each tgmt barcode 
        3. backgroundVsCellsPerBeadBox: Box plot showing distribution of % background split by number of cells in well 
        '''
        groupBlocks = []
        cells = self.cellStats[self.cellStats.species.isin(self.speciesList + ["Ambig"])]
        lvl1_barcodes = f"{self.barcodesToPlot[1]}_alias"
        backgroundPercPerWell = px.box(x=cells.minorFrac*100, y=cells[lvl1_barcodes], points=False,height=700, title=f"% Background by {self.barcodesToPlot[1]} Barcode", boxmode="overlay", template=constants.DEFAULT_FIGURE_STYLE, color=cells[lvl1_barcodes], labels={'y':'', 'x':'Background[%]'})
        lvl1Labels = list(cells[lvl1_barcodes])
        lvl1Labels.sort(reverse=True)
        backgroundPercPerWell.update_yaxes(categoryorder='array', categoryarray=lvl1Labels)
        backgroundPercPerWell.update_layout(showlegend=False)

        cells = self.cellStats[self.cellStats.species.isin(self.speciesList + ["Mixed"])]
        lvl1Counts = cells.groupby([lvl1_barcodes, 'species']).size().reset_index()
        lvl1CountsData = self.getProportionMixed(lvl1Counts, lvl1_barcodes)
        fracMixedCellsPerWell = px.bar(lvl1CountsData, x='percentMixed', y=lvl1_barcodes, height=700, title=f"% Mixed by {self.barcodesToPlot[1]} Barcode", template=constants.DEFAULT_FIGURE_STYLE, color=lvl1_barcodes, labels={'percentMixed': 'Mixed cells [%]', lvl1_barcodes:f'{self.barcodesToPlot[1]} Barcode'})
        fracMixedCellsPerWell.update_layout(showlegend=False)
        fracMixedCellsPerWell.update_yaxes(categoryorder='array', categoryarray=lvl1Labels)
        groupBlocks.append(dp.Plot(fracMixedCellsPerWell))
        groupBlocks.append(dp.Plot(backgroundPercPerWell))
        if self.makeAdvanced:
            lvl2_barcodes = f"{self.barcodesToPlot[2]}"
            lvl2Cnts = cells[lvl2_barcodes].value_counts()
            backgroundVsCellsPerBeadBox = px.box(x=lvl2Cnts[cells[lvl2_barcodes]].array, y=cells.minorFrac*100, labels={'x':'Cells per bead', 'y':'Background[%]'}, title="Background vs cells per bead", template=constants.DEFAULT_FIGURE_STYLE, points=False)
            groupBlocks.append(backgroundVsCellsPerBeadBox)
        return dp.Group(blocks=groupBlocks,columns=2)

    def makeReadsPerSpeciesBarnyardFigures(self) -> dp.Group:
        '''
        Returns dp.Group containing the following figures
        1. Scatter plot of unique reads (per species) per cell (x-axis: species1, y-axis: species2)
        2. Scatter plot: total reads X Saturation cells colored by assigned species 
        3. Scatter plot: Total reads X FractionMitochondrial colored by assigned species 
        '''
        groupBlocks = []
        statsTable = self.makeBarnyardStatsTable()
        groupBlocks.append(statsTable)
        barnyardUniqueReadPlots = self.barnyardScaleChangeScatter("Unique", "Unique", "species", "Reads Per Species", speciesSpecificCol=True)
        groupBlocks.append(barnyardUniqueReadPlots)
        self.cellStats.loc[(self.cellStats.species == 'Ambig') | (self.cellStats.species == 'Mixed'), 'species'] = "None" 
        
        barnyardSaturationPlot = DatapaneUtils.makeScatterPlot(self.cellStats, "TotalReads", "Saturation", "species", "Saturation", True, False, constants.BARNYARD_COLORMAP, self.sampleThresholds['UniqueReads'])
        groupBlocks.append(barnyardSaturationPlot)
        barnyardMitochondrialHist = px.box(self.cellStats, y="FractionMito",  points=False, template=constants.DEFAULT_FIGURE_STYLE,color="species", color_discrete_map=constants.BARNYARD_COLORMAP, title="Mito. Reads", labels={"FractionMito":"Proportion of mitochondrial reads in cell"})
        barnyardMitochondrialHist.update_traces(selector=dict(name='None'), visible="legendonly")
        groupBlocks.append(barnyardMitochondrialHist)
        return dp.Group(blocks=groupBlocks, columns=2)

    def barnyardScaleChangeScatter(self, x:str, y:str, colorBy:str, title:str, speciesSpecificCol:bool=False,log_y:bool=True, log_x:bool=True) -> Union[dp.Select, dp.Plot]:
        '''
        Returns an interactive drop down plot with additional tab to select scale (log or linear) where:
        @data: Data being visualized 
        @discrete_color_map: Is the color map to be used to color categorical variable (specified as @colorBy)
        @speciesList: List of species in the barnyard genome
        @speciesSpecificCol: Are x and y column names species specific?  
        '''
        plotCells = self.cellStats[:self.cellStats['pass'].sum()*2]
        if len(plotCells.index) > constants.SAMPLING_NUMBER:
            plotCells = plotCells.sample(constants.SAMPLING_NUMBER)
        species1 = self.speciesList[0]
        species2 = self.speciesList[1]
        if (speciesSpecificCol):
                x = f"{species1}{x}"
                y = f"{species2}{y}"
        speciesScatterPlot = px.scatter(plotCells, x=x, y=y, color=colorBy, color_discrete_map=constants.BARNYARD_COLORMAP, log_x=False, log_y=False, template=constants.DEFAULT_FIGURE_STYLE, title=title,labels={"HumanUnique":"Human UMIs", "MouseUnique":"Mouse UMIs"})

        speciesScatterPlot.update_traces(selector=dict(name='None'), visible="legendonly")
        speciesScatterPlot.update_traces(marker=dict(size=4, opacity=0.5))
        updatemenus = [
                dict(
                    type="buttons",
                    direction="right",
                    x=0.5,
                    y=-0.4, 
                    xanchor="center",
                    yanchor="bottom",
                    buttons=list([
                        dict(
                            args=[{'yaxis.type': 'linear', 'xaxis.type':'linear'}],
                            label="Linear Scale",
                            method="relayout"
                        ),
                        dict(
                            args=[{'yaxis.type': 'log', 'xaxis.type': 'log'}],
                            label="Log Scale",
                            method="relayout"
                        )
                    ])
                ),
        ]
        speciesScatterPlot.update_layout(updatemenus=updatemenus)
        return dp.Plot(speciesScatterPlot)

    def makeBarnyardStatsTable(self)-> dp.Table:
        '''
        Creates a dataframe containing statistics to be included as a table 
        on the barnyard page
        '''
        speciesAndMixed = self.speciesList +['Mixed']
        tab = self.cellStats[(self.cellStats['pass']) & (self.cellStats['species'] != "None")]
        ncells = len(tab.index)
        stats = []
        fracCellsPerSpecies = {}
        totalFracNonMixed = 0
        for species in speciesAndMixed:
            speciesCells = (tab.species==species).sum()/ncells
            stats.append([f"{species} cells", f"{speciesCells:.1%}"])
            fracCellsPerSpecies[species] = speciesCells
            if species != "Mixed":
                totalFracNonMixed += speciesCells
        humanFrac = fracCellsPerSpecies['Human'] / totalFracNonMixed
        doublets = fracCellsPerSpecies['Mixed'] / (2*humanFrac*(1-humanFrac)) if 0 < humanFrac < 1 else np.NaN
        stats.append(['Background', f"{tab[tab.species != 'Mixed'].minorFrac.mean():.2%}"])
        stats.append(['Estimated doublets', f"{np.min([1,doublets]):.1%}"])
        statsAsDictList = [GeneralUtils.makeTableDict(['Metric','Value'], valuePair) for valuePair in stats] 
        statsDataFrame = pd.DataFrame(statsAsDictList)
        return DatapaneUtils.createTableIgnoreWarning(statsDataFrame[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, title="", hideColumnHeaders=True, boldColumn=['Metric']))
    
    def getProportionMixed(self, df: pd.DataFrame, lvl1Bc) -> pd.DataFrame:
        '''
        Returns dataframe containing proportion of 'Mixed' cells
        for each well for barnyard data where 'Mixed' cells are
        cells consisting of similar number of reads of more than 
        one species
        '''
        outputDicts = []
        df.rename(columns={0: "ReadCount"},inplace=True)
        lvl1Barcpdes = df[lvl1Bc].unique()
        for barcode in lvl1Barcpdes:
            subset = df[df[lvl1Bc] == barcode]
            mixedCount = subset[ subset['species'] == 'Mixed']['ReadCount'].sum()
            allSpeciesCount = subset['ReadCount'].sum()
            percentMixed = (mixedCount / allSpeciesCount) * 100
            outputDicts.append({lvl1Bc:barcode, 'percentMixed': percentMixed})
        return pd.DataFrame(outputDicts)

    def scoreBarn(self) -> None:
        '''
        Adds additional columns to @cellStats via mutation 
        denoting cells as belonging to one of the species in @speciesList
        '''
        self.cellStats['species'] = "None"
        speciesUniqueReadColumns = [f"{species}Unique" for species in self.speciesList]
        self.cellStats['minor'] = self.cellStats[speciesUniqueReadColumns].min(1)
        self.cellStats['minorFrac'] = self.cellStats['minor'] / self.cellStats[speciesUniqueReadColumns].sum(1)

        for species in self.speciesList: 
            otherSpeciesReadCounts = [s for s in speciesUniqueReadColumns if s != f"{species}Unique"]
            self.cellStats.loc[self.cellStats['pass'] & (self.cellStats[f"{species}Unique"] > self.cellStats[otherSpeciesReadCounts].max(1)) & (self.cellStats['minorFrac'] < 0.2), 'species'] = species

        backgroundEachSpecies = dict([(species, self.cellStats[(self.cellStats.species == species)].minorFrac.mean()) for species in self.speciesList])
        
        maximumValues = {}
        minimumValues = {}
        for species in self.speciesList:
            if (self.cellStats.species == species).any():
                uniqueReadsColName= f"{species}Unique"
                maximumValues[species] = np.percentile(self.cellStats[uniqueReadsColName][self.cellStats.species==species], 90)
                minimumValues[species] = max(10, int(maximumValues[species] * backgroundEachSpecies[species] * 5))
                self.cellStats.loc[(self.cellStats[uniqueReadsColName] >= minimumValues[species]), 'species'] = species

        for (species1, species2) in combinations(self.speciesList, 2):
            if (minimumValues.get(species1) is not None and minimumValues.get(species2) is not None):
                species1ReadCol = f"{species1}Unique"
                species2ReadCol = f"{species2}Unique"
                self.cellStats.loc[(self.cellStats[species1ReadCol] < minimumValues[species1]) & (self.cellStats[species2ReadCol] < minimumValues[species2]), 'species'] = "None"
                self.cellStats.loc[(self.cellStats[species1ReadCol] >= minimumValues[species1]), 'species'] = species1
                self.cellStats.loc[(self.cellStats[species2ReadCol] >= minimumValues[species2]), 'species'] = species2
                self.cellStats.loc[(self.cellStats[species1ReadCol] >= minimumValues[species1]) & (self.cellStats[species2ReadCol] >= minimumValues[species2]), 'species'] = "Mixed"
                self.cellStats.loc[(self.cellStats.species == 'Mixed') & (self.cellStats[species1ReadCol] > self.cellStats[species2ReadCol]) & (self.cellStats.minorFrac < min(0.25,backgroundEachSpecies[species1] * 5)), 'species'] = 'Ambig'
                self.cellStats.loc[(self.cellStats.species == 'Mixed') & (self.cellStats[species2ReadCol] > self.cellStats[species1ReadCol]) & (self.cellStats.minorFrac < min(0.25,backgroundEachSpecies[species2] * 5)), 'species'] = 'Ambig'
                self.cellStats.loc[(self.cellStats.species == species1) & (self.cellStats.minorFrac >= max(0.05, backgroundEachSpecies[species1] * 5)), 'species'] = 'Ambig'
                self.cellStats.loc[(self.cellStats.species == species2) & (self.cellStats.minorFrac >= max(0.05, backgroundEachSpecies[species2] * 5)), 'species'] = 'Ambig'

    def buildBarnyardPage(self) -> dp.Page:
        '''
        '''
        pageBlocks = [dp.Text(f"## Sample: {self.sampleName}")]
        self.scoreBarn()
        overloaded = len(self.barcodesToPlot) > 1
        presentClassificiationCategories = list(filter(lambda x: x !='None', self.cellStats['species'].unique()))
        # TODO : GENERALIZE
        if overloaded and len(presentClassificiationCategories) >= 2: 
            overloadingFiguresGroup = self.makeOverloadBarnyardPlots()
            pageBlocks.append(overloadingFiguresGroup)
        regularGroup = self.makeReadsPerSpeciesBarnyardFigures()
        pageBlocks.insert(1, regularGroup)
        return dp.Page(blocks=pageBlocks,  title="Barnyard")

    ###################### PEAKS PAGE FUNCTIONS ######################
    def buildPeaksPage(self) -> dp.Page:
        """
        Build df.Page: "Peak" 
        containing FrIP and FRiT distributions & related stats in a table
        """
        cells = self.cellStats
        isBarnyard = len(self.speciesList) > 1
        if isBarnyard:
            cells = self.cellStats[self.cellStats['species'].isin(self.speciesList)]
        peakCountStats = GeneralUtils.loadStatsAsDict(self.sampleSpecificFilePaths['peakCountsFiles']['metrics'])
        peakCellStats = pd.read_csv(self.sampleSpecificFilePaths['peakCountsFiles']['barcodeStats'], sep="\t", index_col=0)
        tssCellStats = pd.read_csv(self.sampleSpecificFilePaths['tssCountsBarcodeStats'], sep="\t", index_col=0)
        tab = cells.join(tssCellStats)
        tab.rename(columns={"Features":"UniqueTSSCounts" ,"TotalCounts": "TotalTSSCount", "FractionCounted":"FRiT"},inplace=True)
        tab = tab.join(peakCellStats)
        tab.rename(columns={"FractionCounted":"FRiP"}, inplace=True)
        tab = tab[tab['pass']]
        print(tab)
        fritDist = self.makeMultispeciesDistPlot(tab, isBarnyard, self.speciesList, "FRiT", "Reads in TSS", "Fraction of reads in TSS")
        fripDist = self.makeMultispeciesDistPlot(tab, isBarnyard, self.speciesList, "FRiP", "Reads in Peaks", "Fraction of reads in peaks")
        peakStatsDf = self.makePeakStatsDf(peakCountStats,tab)
        peakStatsTable = peakStatsDf[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, title="", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value'])
        peakStatsDPTable =  DatapaneUtils.createTableIgnoreWarning(peakStatsTable)
        peakBlocks = [dp.Text(f"## Sample: {self.sampleName}"), dp.Group(fripDist,fritDist, columns=2),peakStatsDPTable]
        return (dp.Page(
                blocks=peakBlocks,title="Peaks"
        ) ,peakStatsDf)

    def makeMultispeciesDistPlot(self, tab: pd.DataFrame, isBarnyard: bool, speciesList: List[str], field:str, title:str, xaxisTitle:str) -> dp.Plot:
        '''
        Creates distribution plots for @field categorzied by species (if more than one species exists in @speciesList) 
        '''
        plotVals = []
        groupLabels = []
        showDistLegend = False
        if isBarnyard:
            for species in sorted(speciesList):
                speciesSubset = tab[tab.species == species]
                if len(speciesSubset.index) > 0:
                    plotVals.append(list(speciesSubset[field]))
                    groupLabels.append(species)
            showDistLegend = True
        else:
            values = tab[field].replace(0, 0.01)
            plotVals.append(values)
            groupLabels.append('')
        fig = ff.create_distplot(plotVals, group_labels=groupLabels, show_hist=False, show_rug=False, curve_type='kde')
        fig.update_layout(
            title=title,
            xaxis_title=xaxisTitle,
            yaxis_title="KDE",
            showlegend=showDistLegend
        )
        return dp.Plot(fig)

    def makePeakStatsDf(self, peakCountStats: Dict, tssAndPeakDf: pd.DataFrame) -> pd.DataFrame:
        '''
        Calculates peaks related statistics of ATAC run, writing them to a DataFrame
        to be displayed in the report
        '''
        stats = []
        stats.append(["Total Peaks", int(peakCountStats.get('Unique features', 0))])
        stats.append(["Median Reads In TSS Regions", f"{tssAndPeakDf.FRiT.median():.1%}"])
        stats.append(["Median Reads in Peaks",f"{tssAndPeakDf.FRiP.median():.1%}"])
        statsAsDictList = [GeneralUtils.makeTableDict(['Metric','Value'], valuePair) for valuePair in stats]
        df = pd.DataFrame(statsAsDictList)
        df['Category'] = "Peaks"
        return df 

    def splitBarcodes(self):
        lvl1BarcodeAliases = {}
        for name,bc in self.demuxJson["samples"][self.sampleName]["barcodes"].items():
            lvl1BarcodeAliases[bc['sequence']] =  name 
        libDirJSON = GeneralUtils.readJSON(self.libDir, preserveDictOrder=True)
        barcodeInfo = libDirJSON['barcodes']
        scDefiningBarcodes = [barcodeMeta for barcodeMeta in barcodeInfo if barcodeMeta.get('type', None) != 'library_index']
        groupByLevel = groupby(scDefiningBarcodes, key=lambda x: x['level'])
        expectedBarcodeCount = len(scDefiningBarcodes)
        splitByBarcode = self.cellStats.index.str.split('+')
        observedBarcodeCount = len(splitByBarcode[0])
        if observedBarcodeCount == expectedBarcodeCount:
            curBcIndex = 0
            barcodesToPlot = {}
            for level, group in groupByLevel:
                groupAsList = list(group)
                groupSize = len(groupAsList)
                barcodeName = groupAsList[0]['alias']
                barcodes = [''.join(x[curBcIndex:(curBcIndex+groupSize)]) for x in splitByBarcode]
                self.cellStats[barcodeName] = barcodes
                barcodesToPlot[level] = barcodeName
                curBcIndex += groupSize
                if level == 1 and len(lvl1BarcodeAliases) > 0:
                    self.cellStats[f"{barcodeName}_alias"] = self.cellStats[barcodeName].apply(lambda x: lvl1BarcodeAliases[x])
        else: 
             raise ValueError(f"The number of sequences in '+' split cell barcodes ({observedBarcodeCount})\ndoes not match the inferred number of barcodes ({expectedBarcodeCount}).\nEnsure the specified library definition is correct for this data")
        self.barcodesToPlot = barcodesToPlot

    def getSpeciesList(self) -> List[str]:
        '''
        @Assumption read counts associated with all constituent species 
        are listed as <species>Unique in the passed @cellStats dataframe  
        '''
        columns = self.cellStats.columns
        self.speciesList = [x.replace("Unique", "", 1) for x in columns if x.endswith("Unique")]

    def resolveSampleSpecificFilePaths(self) -> Dict[str, Union[Path, Dict[str, Path]]] :
        """
        Returns a dictionary of the needed file paths. 
        Differs based on whether -results argument was 
        provided during invocation  
        - isn't provided: (running from within ScaleAtac workflow)
        - provided: running independently 
        """
        outputDict = {}
        if str(self.resultsDir) == '.':
            outputDict['filterDfPath'] = Path(f"{self.sampleName}_QC.tsv")
            outputDict['threshJsonPath'] = Path(f"{self.sampleName}_QC.json")
            outputDict['dedupStatsPath'] = Path(f"{self.sampleName}.dedup_stats.tsv")
            outputDict['cellStatsPath'] = Path(f"{self.sampleName}.cell_stats.tsv")
            outputDict['fragmentHistPath'] = Path(f"{self.sampleName}.fragment_hist.tsv")
            countsDir = Path(f"{self.sampleName}.counts")
            outputDict['peakCountsFiles']= {'barcodeStats': countsDir / "barcodes.stats.tsv", 'metrics': countsDir / "metrics.tsv"}
            tssCountsDir = Path(f"{self.sampleName}.tss_counts")
            outputDict['tssCountsBarcodeStats'] = tssCountsDir / "barcodes.stats.tsv"
        else:
            QCPrefix = self.resultsDir / 'QC' / self.sampleName / self.qcSubdirectory
            outputDict['filterDfPath'] = QCPrefix / f'{self.sampleName}_QC.tsv'
            outputDict['threshJsonPath'] = QCPrefix / f'{self.sampleName}_thresholds.json'
            dedupPrefix = self.resultsDir / 'align' / 'dedup'  
            outputDict['dedupStatsPath'] = dedupPrefix / f'{self.sampleName}.dedup_stats.tsv'
            outputDict['cellStatsPath'] = dedupPrefix / f'{self.sampleName}.cell_stats.tsv'
            outputDict['fragmentHistPath'] = dedupPrefix / f'{self.sampleName}.fragment_hist.tsv'
            peaksPrefix = self.resultsDir / 'peaks' 
            countsDir = peaksPrefix / f'{self.sampleName}.counts'
            outputDict['peakCountsFiles'] = {'barcodeStats': countsDir / "barcodes.stats.tsv", 'metrics': countsDir / "metrics.tsv"}
            tssCountsDir = peaksPrefix / f'{self.sampleName}.tss_counts'
            outputDict['tssCountsBarcodeStats'] = tssCountsDir / "barcodes.stats.tsv"
        return outputDict

