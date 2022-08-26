#!/usr/bin/env python

from reportGeneration.ReportUtil import GeneralUtils, CalculationUtils
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict
from itertools import groupby

class AbstractReport(ABC):
    resultsDir: Path
    makeAdvanced: bool
    demuxJson: Dict[str, float]
    resultsDir: Path
    writeDir: Path
    qcSubdirectory: Path
    libDir:Path
    barcodesToPlot: Dict[int, str]
    
    @abstractmethod
    def __init__(self, resultsDir: Path, makeAdvanced:bool, demuxJson: Dict[str,float], qcDir:Path, identifier: str, libDir:Path):
        self.resultsDir = resultsDir
        self.makeAdvanced = makeAdvanced
        self.demuxJson = demuxJson
        self.resultsDir = resultsDir
        self.writeDir = AbstractReport.resolveWriteDir(resultsDir, qcDir, identifier)
        self.qcSubdirectory = qcDir
        self.libDir = libDir

    @abstractmethod
    def build(self):
        '''
        Builds and saves an interactive HTML report
        '''
        pass

    @staticmethod
    def resolveWriteDir(resultsDir: Path, qcDir: Path, identifier:str) -> Path: 
        return resultsDir / 'reports' / qcDir

    @abstractmethod
    def buildCellsPage(self):
        '''
        Creates a dp.Page to be included in the HTML report labelled "Cells"
        '''
        pass
    

