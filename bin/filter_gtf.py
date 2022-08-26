#!/usr/bin/env python
import json
import pybedtools
import pybedtools.featurefuncs
import argparse


def filter_trans(x, filterCriteria):
    if 'fields' in filterCriteria.keys():
        for (key, value) in filterCriteria['fields'].items():
            if x.fields[int(key)] not in value:
                return False        
    if 'attrs' in filterCriteria.keys():
        for (key,value) in filterCriteria['attrs'].items():
            if  x.attrs.get(key, '') not in value:
                return False
    return True

# Filters @gtfFile based on criteria passed in @filterJSON. 
# outputs filtered result to @outputPath 
def main(gtfFile, filterJSON, outputPath):
    gtf = pybedtools.BedTool(gtfFile)
    with open(filterJSON) as f:
        str = f.read()
        strStripped = str.rstrip()
        filterCriteria = json.loads(strStripped)
    tss = gtf.filter(filter_trans, filterCriteria=filterCriteria)
    tss.saveas(outputPath)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GTF filtering tool. See ScaleAtac/README for usage guidance")
    parser.add_argument("-gtf", metavar="gtfFilePath", type=str, help="Path to GTF being filtered")
    parser.add_argument("-j", metavar="filterJSON", type=str, help="JSON denoting which fields to filter and by what. See README for further guidance")
    parser.add_argument("-o", metavar="outputPath", type=str, help="Path of GTF to be outputted")
    args = parser.parse_args()
    main(args.gtf, args.j, args.o)
    
