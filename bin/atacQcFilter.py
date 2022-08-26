#!/usr/bin/env python
import os
import pandas as pd 
from typing import Dict, Union
import numpy as np
import json 
import argparse
import plotly.express as px
import sys 
from pathlib import Path, PosixPath


# Generates a table denoting each barcode as passed or failed 
# based on: 
# - # of UniqueReadsInPeaks -- threshold=(Median - (2*IQR)) of @topCellPercent% by UniqueReadsInPeaks)
# - Fraction of reads in peaks -- threshold=(Median - (2*IQR)) of @topCellPercent% by UniqueReadsInPeaks)
def ATACqc(outputdir:str, inputdir:str, sampleName: str, minreads: int, expectedCells: int, background:float, topCellPercent: float, hardReadLimit: int):
	minCellRatio = background * 5
	inputDirPath = Path(inputdir.rstrip("/")) if inputdir is not None else Path('.')
	makeDir(inputDirPath / 'QC' )
	makeDir(inputDirPath / 'QC' / sampleName)
	writeDir = inputDirPath / 'QC' / sampleName / outputdir if outputdir is not None else inputDirPath / 'QC' / sampleName
	makeDir(writeDir)
	all_stats = aggregateStats(inputdir, sampleName)
	if (hardReadLimit is None):
		readThreshold = int(calcFieldThreshold('UniqueReads',all_stats, expectedCells, topCellPercent, minCellRatio, minreads))
	else:
		readThreshold = hardReadLimit
	fripThreshold = calcFieldThreshold('FRiP', all_stats, expectedCells, topCellPercent, minCellRatio, minreads)
	writeFilteringResults(writeDir,sampleName, readThreshold,fripThreshold, all_stats)
	metadata = {'UniqueReads':readThreshold, "FRiP":fripThreshold}
	writeMetadataJSON(writeDir,sampleName, metadata)

# Writes the passed Dict as a JSON object for use by other scripts for further
# downstream analysis
def writeMetadataJSON(writeDir:PosixPath,sampleName: str, data: Dict[str, float]):
	json_dict = json.dumps(data)
	with open(writeDir / f'{sampleName}_thresholds.json', 'w') as outfile:
		outfile.write(json_dict)

# Writes a tsv containing information about whether filtering criteria were passed: 
# - FRiP
# - UniqueReads 
# Mutates the dataframe to contain additional columns associated with decided filters 
def writeFilteringResults(writeDir:PosixPath, sampleName: str, readThreshold:float,fripThresh:int, atacData: pd.DataFrame):
	atacData['Saturation'] = np.around(1-atacData['UniqueReads'] / atacData['PassingReads'], decimals=3)
	atacData['pass'] = atacData['UniqueReads'] >= readThreshold
	atacData['FractionMito'] = np.around(atacData['MitoReads'] / atacData['TotalReads'],decimals=3)
	# Filtering field creation
	atacData['Filter'] = 'Pass'
	uniqueReadsAboveThresh = (atacData['UniqueReads'] >= readThreshold)
	fripAboveThresh = (atacData['FRiP'] >= fripThresh)
	atacData.loc[(~uniqueReadsAboveThresh), 'Filter'] = 'LowUniqueReads'
	atacData.loc[(uniqueReadsAboveThresh & ~fripAboveThresh), 'Filter'] = "LowFRiP"
	tsvPath = writeDir  / f'{sampleName}_QC.tsv'
	atacData.sort_values(by='TotalReads', ascending=False, inplace=True)
	atacData = atacData[['TotalReads', 'PassingReads', 'UniqueReads', 'Saturation', 'MitoReads', 'FractionMito', 'PeakCount', 'ReadsInPeaks', "UniqueReadsInPeaks", "FRiP", "TSSCount", "ReadsInTSS", "FRiT", "Filter", "pass"]]
	atacData.to_csv(tsvPath, sep='\t', index_label='Barcode')
	makeKneePlots(atacData, writeDir)

# Joins statistics in cell_stats.tsv and barcodes.stats.tsv to create 
# a dataframe containing metadata about the data at the given @inputdir
# associated with the given @sampleName 
def aggregateStats(inputdir: Union[str,None], sampleName: str) -> pd.DataFrame:
	# Running w/i nextflow 
	cellStatsPath = Path(f"{sampleName}.cell_stats.tsv")
	peakStatsPath = Path(f"{sampleName}.counts", "barcodes.stats.tsv")
	tssStatsPath = Path(f"{sampleName}.tss_counts", "barcodes.stats.tsv")
	# Executing independently 
	if (inputdir is not None):
		inputdir = Path(inputdir)
		cellStatsPath = inputdir / 'align' / 'dedup'/ cellStatsPath
		peakStatsPath = inputdir /'peaks' / peakStatsPath	
		tssStatsPath = inputdir / 'peaks' / tssStatsPath
	cell_stats = pd.read_csv(cellStatsPath, sep='\t', index_col='Barcode')
	peak_stats = pd.read_csv(peakStatsPath, sep='\t', index_col='CellBarcode')
	peak_stats.rename(columns={"Features":"PeakCount", "TotalCounts": "ReadsInPeaks", "FractionCounted":"FRiP"},inplace=True)
	all_stats = pd.merge(cell_stats, peak_stats, left_index=True, right_index=True)
	tss_stats = pd.read_csv(tssStatsPath, sep='\t', index_col='CellBarcode')
	tss_stats.rename(columns={"TotalCounts":"ReadsInTSS", "FractionCounted":"FRiT", "Features":"TSSCount"}, inplace=True)
	all_stats = pd.merge(all_stats, tss_stats, left_index=True, right_index=True)
	all_stats["UniqueReadsInPeaks"] = all_stats['FRiP'] * all_stats['UniqueReads']
	return all_stats

# Calculates a field threshold value as the (100-topCellPercent) amongst cells with >= minReads 
def calcFieldThreshold(field:str, all_stats: pd.DataFrame, expectedCells: int, topCellPercent: int, minCellRatio: float, minReads:int):
	all_stats.sort_values(by=field, inplace=True, ascending=False)
	if (expectedCells == 0):
		expectedCells=(all_stats['UniqueReads'] >= minReads).sum()
	threshold = np.percentile(all_stats[field][:expectedCells], (100-topCellPercent)) * minCellRatio
	threshold = max(threshold, minReads) if field == 'UniqueReads' else threshold
	return threshold

def makeKneePlots(dataset, writeDir):
	passingCellCount = len(dataset[dataset['Filter'] == 'Pass'].index)
	vals = pd.DataFrame()
	cellbarcodeCounts = list(range(len(dataset.index)))
	vals['CellBarcodes'] = cellbarcodeCounts
	vals['FRiP'] = list(dataset.sort_values(by='FRiP', ascending=False)['FRiP'])
	uniqueReads = list(dataset.sort_values(by='UniqueReads', ascending=False)['UniqueReads'])
	vals['UniqueReads'] = uniqueReads
	uniqueReadsKneePlot = px.line(vals, x='CellBarcodes', y='UniqueReads',  log_x=True, log_y=True, template="none")
	uniqueReadsKneePlot.add_vline(x=passingCellCount, line_dash="dash", line_color="green")
	uniqueReadsKneePlot.write_image(writeDir / "uniqueReadsKneePlot.png")
	fripKneePlot = px.line(vals, x='CellBarcodes', y='FRiP',  log_x=True, template="none")
	fripKneePlot.add_vline(x=passingCellCount, line_dash="dash", line_color="green")
	fripKneePlot.write_image(writeDir / "fripKneePlot.png")
	
# If it doesn't yet exist, makes a directory with the given @dirName 
def makeDir(dirName:PosixPath):
	if not os.path.exists(dirName): 
		os.mkdir(dirName)

def main():
	parser = argparse.ArgumentParser(description="Taken in needed paramters")
	hardReadLimitNotPassed = '--hardUniqueReadsLimit' not in sys.argv
	parser.add_argument('--results', metavar='resultsDirectory', type=str, help="Directory containing completed analysis to QC", required=False, default=None)
	parser.add_argument('--sample', metavar='sampleName', type=str, help="Name of the sample to run QC on")
	parser.add_argument('--writeDir', required=False, default=None, metavar='outputDirectory',  type=str, help="Name of directory within QC/<sample> which to write output to")
	parser.add_argument('--minUniqueReads', required=hardReadLimitNotPassed, metavar='minReads', default=2,  type=int, help="(if expected cells = 0) Minimum number of unique reads a cell must to be included in threshold calculation")
	parser.add_argument('--expectedCells', required=hardReadLimitNotPassed, metavar='expectedCells', default = 0, type=int, help="Rough estimate of the cell-number. Used to help find a good cutoff threshold")
	parser.add_argument('--background', required=hardReadLimitNotPassed, metavar='background', type=float, default=0.05, help="Proportion of total reads expected to be background. (Between 0 and 1) Default=0, meaning no reads are expected to be background")
	parser.add_argument('--topPercentCells', required=hardReadLimitNotPassed, metavar='topCellPercent', type=float, default=0.10, help="Top percent of cells which will be used to calculate QC metric thresholds")
	parser.add_argument('--hardUniqueReadsLimit', required=False, metavar="minUniqueReads", type=int, default=None, help="Threshold that when specified overrides calculated threshold for minimum unique reads")
	args = parser.parse_args()
	ATACqc(args.writeDir, args.results, args.sample, args.minUniqueReads, args.expectedCells, args.background, args.topPercentCells, args.hardUniqueReadsLimit)

if __name__ == "__main__":
   main()
