#!/usr/bin/env python
import sys
import subprocess
import shutil
import json

def subsample(in_bam, fraction, out_bam):
    print(f"Subsampling {in_bam} at {fraction:.2%} to {out_bam}")
    subprocess.run(["samtools", "view", "-b", "-s", f"{fraction:.3}", "-o", out_bam, in_bam])

def parse_read_counts(demux_metrics):
    read_counts = {}
    for l in open(demux_metrics):
        line = l.strip().split()
        if not line: continue
        if line[0] == "Demux":
            read_counts[line[1]] = int(line[2])
    return read_counts

def main(in_bam, sample, reads, demux_metrics, out_bam,):
    counts = parse_read_counts(demux_metrics)
    if not sample in counts:
        print("Missing read count for ", sample, " in ", demux_metrics, file=sys.stderr)
        sys.exit(1)
    if 0 < reads < counts[sample]:
        subsample(in_bam, reads / counts[sample], out_bam)
    else:
        print(f"Copying {in_bam} to {out_bam}")
        shutil.copy(in_bam, out_bam)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("USAGE: subsample.py INPUT.BAM SAMPLE_NAME READ_NUM DEMUX_METRICS.txt OUT.BAM", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])
