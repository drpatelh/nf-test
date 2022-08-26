import numpy as np 
import json
from typing import Dict, Union, Tuple, List
import pandas as pd
from os import makedirs
from os.path import exists 
from pathlib import Path
import warnings
import datapane as dp
from math import log10 
import functools 
import operator
from collections import OrderedDict
import reportGeneration.myconstants as constants
import plotly.express as px

## GENERAL PURPOSE UTILITY FUNCTIONS RELATED TO IO, INPUT PARSING, OR PANDAS MANIPULATIONS  
class GeneralUtils:
    def __init__(self):
        self=self

    @staticmethod
    def loadStatsAsDict(fn):
        '''
        Create a dictionary from a file (@fn) containing tab delimited statistics
        '''
        res = {}
        for l in open(fn):
            line = l.strip().split("\t")
            if not len(line) >= 2: continue
            stat = line[0]
            val = line[1]
            try:
                val = int(val)
            except ValueError:
                try: 
                    val = float(val)
                except ValueError:
                    pass
            res[stat] = val
        return res

    @staticmethod
    def makeTableDict(colNames, values):
        '''
        Create a simple dictionary from @pair values (for use in pd.DataFrame creation)
        Assumption: 
        - @colNames and @values are equal length where the elements at each index correspond to each other
        '''
        resultsDict = {}
        for i,col in enumerate(colNames):
            resultsDict[col] = values[i]
        return resultsDict

    @staticmethod 
    def reformatIfInt(val):
        if isinstance(val, int):
            return f"{val:,}"
        elif isinstance(val, float):
            return round(val, 2)
        else:
            return val

    @staticmethod
    def styleTable(styler, title:str, hideColumnHeaders=False, boldColumn=None, numericCols=None):
        '''
        Modifies given @pd.DataFrame.Styler 
        '''
        if (numericCols is not None):
            styler.format(GeneralUtils.reformatIfInt, subset=numericCols)
        styler.hide(axis='index')
        if hideColumnHeaders:
            styler.hide(axis='columns')
        else:
            styler.set_table_styles([{'selector':'th', 'props':[('border', 'black solid !important')]}], overwrite=False)
        if boldColumn is not None:
            styler.set_properties(
                subset=boldColumn,
                **{'font-weight': 'bold'}
            )
        if title != "":
            styler.set_caption(title)
            styler.set_table_styles([{'selector':'caption', 'props':[('border', 'black solid !important'),  ("font-weight", "bold")]}],overwrite=False)
        styler.set_properties(**{"border-color":'black', "border-style":'solid !important'})
        return styler

    @staticmethod
    def readJSON(file, preserveDictOrder:bool=False):
        """
        Reads in JSON as a python object
        """
        with open(file) as f:
            str = f.read()
            strStripped = str.rstrip()
            pairs_hook = OrderedDict if preserveDictOrder else None
            parsedJSON = json.loads(strStripped, object_pairs_hook=pairs_hook)
        return parsedJSON

    @staticmethod
    def ensurePathsExist(filePathDict:Dict[str, Union[str,Dict]]):
        '''
        Ensures all paths mentioned in the given @filePathDict exist
        '''
        for key, value in filePathDict.items():
            if type(value) is dict:
                GeneralUtils.ensurePathsExist(value)
            else: 
                if not exists(value):
                    raise FileNotFoundError(f"{key} was assumed to be located at '{str(value)}'. It is missing")

    @staticmethod
    def makeDir(dirName:Path):
        makedirs(dirName, exist_ok=True)



### FIGURE CREATION METHODS ### 
class DatapaneUtils:
    @staticmethod
    def createTableIgnoreWarning(styler, label=None) -> dp.Table:
        '''
        Wrapper around dp.Table which suppresses warning which arrises from 
        code within dp.Table calling Styler.render():
        'this method is deprecated in favour of `Styler.to_html()`' 
        '''
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            return dp.Table(styler, label=label)


    @staticmethod
    def makeScatterPlot(dataFrame: pd.DataFrame, xVar:str, yVar:str, colorBy:str, title:str, isLogX:bool, isLogY:bool, colorMap:Dict[str,str], readThresh:int) -> dp.Plot:
        '''
        returns dp.Plot rendering of an interactive plotly express scatter plot where: 
        - @dataFrame: is data being visualized
        - @xVar: x variable being plotted 
        - @yVar: y var being plotted 
        - @colorBy: categorical variable points are colored by 
        - @logX, logY: specify whether those axes are to be log scaled
        '''
        dataFrame = dataFrame[dataFrame['UniqueReads'] >= readThresh]
        cellsToPlot = dataFrame[:dataFrame['pass'].sum()*2]
        if len(cellsToPlot.index) > constants.SAMPLING_NUMBER:
            cellsToPlot = cellsToPlot.sample(constants.SAMPLING_NUMBER)
        scatterPlot = px.scatter(cellsToPlot, x=xVar, y=yVar, color=colorBy, title=title, log_x=isLogX, log_y=isLogY,color_discrete_map=colorMap, template=constants.DEFAULT_FIGURE_STYLE)
        scatterPlot.update_traces(marker=dict(size=4, opacity=0.5))
        scatterPlot.update_traces(selector=dict(name='None'), visible="legendonly")
        return dp.Plot(scatterPlot)


    @staticmethod
    def makeKneePlot(data: pd.DataFrame, field:str, title:str) -> dp.Plot:
        '''
        Makes a kneeplot using @field in @data; drawing a vertical line at @threshold
        '''    
        indices = CalculationUtils.getIndicesToInclude(len(data.index))
        vals = pd.DataFrame()
        cellbarcodeCounts = list(range(len(data.index)))
        vals['Cell Barcodes'] = [cellbarcodeCounts[i] for i in indices]
        uniqueReads = list(data.sort_values(by=field, ascending=False)[field])
        vals['Unique Reads'] = [uniqueReads[i] for i in indices]
        fig = px.line(vals, x='Cell Barcodes', y='Unique Reads',  title=title, log_x=True, log_y=True, template=constants.DEFAULT_FIGURE_STYLE)
        fig.add_vline(x=data['pass'].sum(), line_dash="dash", line_color="green")
        # Xlimit set in powers of 10 (10, maximum)
        fig.update_layout(xaxis_range=[1, log10(vals['Cell Barcodes'].max())])
        return dp.Plot(fig)

    def showAllTicks(plotlyFig):
        plotlyFig.update_layout(
            xaxis = dict(
                tickmode='linear',
                tick0 = 1,
                dtick=1
            )
        )


class CalculationUtils:
    def __init__(self):
        self=self

    @staticmethod
    def coverage_equation(x, c, n):
        return c / x - 1 + np.exp(-n / x)

    @staticmethod
    def estLibSize(reads, umis):
        '''
        given @reads total reads and @umis observer UMIS, 
        estimates library size (number of unique molecules) 
        at infinite sequencing depth
        '''
        if umis > reads: raise ValueError()
        if umis <= 0: raise ValueError()
        m = 1
        M = 100
        while CalculationUtils.coverage_equation(M * umis, umis, reads) > 0:
            M *= 10
        for _ in range(40):
            r = (m + M) / 2
            u = CalculationUtils.coverage_equation(r * umis, umis, reads)
            if (u > 0):
                m = r
            elif (u < 0):
                M = r
            else:
                break
        return (umis * (m + M) / 2)


    @staticmethod
    def getIndicesToInclude(elementCount:int, result: List[int]=[])-> List[int]:
        '''
        Where @elementCount is the total number of elements in a 
        list. Returns a list of indices in ascending order for 
        elements to be included. Used to reduce size of line 
        plots while retaining overall trend shape    
        '''
        logVal = int(round(log10(elementCount)))
        nextDown = int(logVal - 1)
        lowerLimit = int(10**nextDown)
        step = ((logVal)**3) if logVal > 0 else 1 
        lowerLimit = lowerLimit - 1 if (lowerLimit > 0) else lowerLimit
        result.append(list(range(lowerLimit, elementCount, step)))
        if (lowerLimit > 0):
            return CalculationUtils.getIndicesToInclude(lowerLimit, result)
        else:
            return sorted(functools.reduce(operator.iconcat, result, []))



    @staticmethod
    def extrapolate_umis_at_nreads(reads, umis, seqDepthInReads):
        '''
        1. estimates library size at infinite sequencing depth
        2. extrapolates the expected number of unique molecules 
        at @seqDepthInReads
        '''
        extrapolatedDiversity = CalculationUtils.estLibSize(reads, umis)
        return int(round(seqDepthInReads * (1 - ((seqDepthInReads - 1)/(seqDepthInReads))**extrapolatedDiversity)))

    @staticmethod 
    def getStatsAsDataframes(filterDfPath: str, cellStatsPath:Union[str,None], threshJson: Dict) -> Tuple[pd.DataFrame,pd.DataFrame]:
        '''
        Parses in @filterDfPath (outputted by QC filtering step) & @cellStatsPath as dataframes, adding additional columns 
        based on existing columns
        '''
        df = pd.read_csv(filterDfPath, sep='\t', index_col="Barcode")
        cellStats = None
        if (cellStatsPath is not None):
            cellStats = pd.read_csv(cellStatsPath, sep="\t", index_col=0)
            cellStats['pass'] = cellStats['UniqueReads'] >= threshJson['UniqueReads']
            cellStats['Saturation'] = np.around(1-cellStats['UniqueReads'] /cellStats['PassingReads'],decimals=3)
            cellStats['FractionMito'] = np.around(cellStats['MitoReads'] / cellStats['TotalReads'],decimals=3)
        return (df, cellStats)

