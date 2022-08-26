
#!/usr/bin/env python
# Plotly express figure style to be used for all figures 
DEFAULT_FIGURE_STYLE="none"
# Number of cells sampled for scatter plots 
SAMPLING_NUMBER = 4000
# Color mapping used for qc filter categorized scatter plots 
QC_COLORMAP = {'Pass': 'rgb(39, 139, 176)', 'LowUniqueReads': 'rgb(170, 174, 181)', 'LowReadsInPeaks':'rgb(28, 163, 89)'}
# Color mapping used for barnyard species categorization scatter plots (human and mouse hardcoded as it is the only barnyard genome used [could be changed])
BARNYARD_COLORMAP = {"None": 'rgb(179, 188, 201)', "Ambig": 'rgb(201, 147, 166)', 'Mixed': 'rgb(242, 5, 33)', 'Human': 'rgb(36, 36, 227)', 'Mouse': 'rgb(27, 99, 25)'}
# 
UNKNOWN_COLORMAP={"Unknown": 'rgb(179, 188, 201)'}

BARCODE_SHORTHAND_TO_NAME={'tgmt':'Tagmentation Barcodes', 'drop':'Droplet Barcodes', 'bead1': 'Bead 1 Barcodes', 'bead2': 'Bead 2 Barcodes', 'bead3': 'Bead 3 Barcodes'}

# Columns to display in stats tables 
DISPLAY_COLUMNS=["Metric","Value"]