#!/usr/bin/env python

### Generate a TSS regions (BED) file from a gene annotation (GTF)
### Usage python tss_regions.py annotation.gtf
### Output: annotation.tss.bed

### TSS BED: TSS Location +/- 1kb, gene_name in name field 
import sys
import pybedtools
import pybedtools.featurefuncs
import pathlib

'''
Assumption: 
GTF has been filtered using bin/filter_gtf.py before running the pipeline 
'''
def get_tss(trans):
    pos = None
    if trans.strand == "+":
        pos = trans.start
    elif trans.strand == "-":
        pos = trans.end
    geneNameAttrs = ['gene_name', 'gene', 'gene_id']
    geneName = None
    for attr in geneNameAttrs:
        if attr in trans.attrs.keys():
            geneName = trans.attrs[attr]
            break
    if geneName is None:
        raise Exception(f"GTF: {sys.argv[1]} does not have any of the attributes associated with its features' gene names:{geneNameAttrs}\n one is needed to derive a .tss.bed")
    return pybedtools.Interval(trans.chrom, max(1,pos-1000), pos+1000, geneName, '.', trans.strand)
    
gtf = pybedtools.BedTool(sys.argv[1])
tss = gtf.each(get_tss)
mtss = tss.sort().merge(c="4,5,6", o="distinct,first,first", s=True)
mtss.saveas(pathlib.Path(sys.argv[1]).stem + ".tss.bed")
