#!/usr/bin/env python
import json
import csv
import xml.etree.ElementTree as ET
import argparse

# Fixed samplesheet.csv settings
SETTINGS = {
"CreateFastqForIndexReads":"1",
"TrimUMI":"0",
"MinimumTrimmedReadLength":"16",
"MaskShortReads":"16",
}

def load_fqBarcodes(tsv):
    """Load named fastq index sequences (e.g. S701P, ...) from a .tsv file
    Each name can refer to one or multiple sequences (when using multiple indices for one sample)
    Return dict name -> [seqs]
    """
    if not tsv: return {}
    bcs = {}
    for line in open(tsv):
        if line.startswith("#"): continue
        l = line.strip().split()
        if len(l) <= 1: continue
        bcs[l[0]] = l[1:]
    return bcs

def load_lib(jlib):
    """Load fastq generation settings from library.json"""
    res = {}
    lib = json.load(open(jlib)).get("bcl_convert", {})
    res['adapter'] = lib.get('adapter', None)
    res['OverrideCycles'] = lib.get('OverrideCycles', None)
    return res

def load_run(runInfo):
    """Load read-lengths from RunInfo.xml in Illumina RunFolder"""
    reads = []
    xml =  ET.parse(open(runInfo))
    for read in xml.getroot().findall("Run/Reads/Read"):
        reads.append((read.attrib['NumCycles'], read.attrib['IsIndexedRead'] == 'Y'))
    return reads

def load_samples(samplesCsv, bcs):
    """Load fastq_samples from samples.csv
    Returns {name: [indexSeqs], ...}
    """
    res = {}
    with open(samplesCsv) as csvfile:
        for row in csv.DictReader(csvfile):
            name = row.get("fastqName", '').strip()
            if not name:
                raise ValueError(f"Missing fastqName name for sample {row.get('sample','')}")
            for n in name:
                if not n.isalnum() or n in "-.":
                    raise ValueError(f"sample name should only contain [a-z],[A-Z],[0-9], dash (-) or dot (.)")
            seqs = row.get("index", name).strip()
            if seqs in bcs:
                seqs = bcs[seqs]
            else:
                for n in seqs:
                    if n not in 'ACTGactg;':
                        raise ValueError(f"Unknown index name / sequence: {seqs}")
                seqs = seqs.split(';')

            if name not in res:
                res[name] = seqs
            elif res[name] != seqs:
                raise ValueError(f"Mismatched index for fastqName {name}")
    return res

def print_settings(settings):
    """Samplesheet.csv settings section"""
    print("[Settings]")
    for (s,v) in settings.items(): 
        print(f"{s},{v}")


def main(samplesCsv, jlib, runInfo, fqBarcodes, settings):
    lib = load_lib(jlib)
    bcs = load_fqBarcodes(fqBarcodes)
    reads = load_run(runInfo)
    # Adapter trimming
    if (adapt := lib['adapter']):
        settings['AdapterRead1'] = adapt
        settings['AdapterRead2'] = adapt
    # If a barcode is in an index read, we need to define that read as a 'UMI' in OverrideCycles
    if (oc := lib.get('OverrideCycles', None)):
        if len(oc) != len(reads): raise ValueError("Mismatched number of reads between lib.json and runinfo.xml")
        overrideCycles = []
        for (read_type, (length, index)) in zip(oc, reads):
            overrideCycles.append(f"{read_type}{length}")
        settings['OverrideCycles'] = ';'.join(overrideCycles)

    print_settings(settings)

    print("[Data]")
    print("Sample_ID", "index", sep=',')
    for name, seqs in load_samples(samplesCsv, bcs).items():
        for i in seqs:
            # One line per index sequence (with repeated sampleName if multiple indicies)
            print(name, i, sep=',')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create bcl_convert samplesheet.csv from workflow samples.csv')
    parser.add_argument('samples', metavar='SAMPLES.csv',
                    help='CSV with samples and index sequences for scATAC workflow run')
    parser.add_argument('lib', metavar='LIBRARY.json',
                    help='Library structure definition')
    parser.add_argument('runinfo', metavar='RUNINFO.xml',
                    help='Sequencer runinfo (in runfolder)')
    parser.add_argument('--fastqIndex', 
                    help='List of fastq index sequences by name')
    args = parser.parse_args()

    main(args.samples, args.lib, args.runinfo, args.fastqIndex, SETTINGS)
