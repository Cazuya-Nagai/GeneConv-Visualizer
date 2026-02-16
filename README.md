# GeneConv Visualizer

A Python-based tool to visualize **gene conversion tracts** and **concerted evolution hotspots** in mitochondrial genomes, specifically designed for analyzing duplication regions.

## Overview
This tool automates the processing of `GENECONV` output files (PI, PO, or GI fragments). It calculates stacking coverage at each nucleotide position and generates publication-ready figures that highlight the intensity of concerted evolution within specific genomic regions.

## Key Features
- **Automated Parsing**: Extracts fragment data directly from raw GENECONV text outputs.
- **Dynamic Statistics**: Automatically calculates the percentage of gene conversion coverage within the duplication region and its terminal segment.
- **Accurate Visualization**: Uses step-plot rendering to precisely represent discrete fragment coverage across the mitogenome.
- **Customizable**: Flexible P-value and fragment length filtering via command-line arguments.

## Installation
Ensure you have Python 3.x installed along with the required libraries:
```bash
pip install pandas numpy matplotlib

##requirements
pandas
numpy
matplotlib

Usage
Run the script from your terminal by providing the GENECONV output file and your duplication coordinates:

Bash
python geneconv_visualizer.py <input_file> --start <start_pos> --end <end_pos> [options]
Example:
Bash
python geneconv_visualizer.py geneconv_output.txt --start 4200 --end 6700 --out fig4_lapwing.png
Arguments:
input: Path to the raw GENECONV output file.

--start: Start coordinate of the duplication region (required).

--end: End coordinate of the duplication region (required).

--out: Name of the output image (default: geneconv_plot.png).

--p: P-value threshold for significance (default: 0.05).

--minlen: Minimum fragment length to be included (default: 500).

Citation
If you use this tool or these findings in your research, please cite our manuscript:

Nagai, K. et al. (2026). Structural Instability and Concerted Evolution in the Mitochondrial Control Region of the Grey-headed Lapwing (Vanellus cinereus) During Range Expansion. Journal of Avian Biology. (In Revision / Corresponding to Manuscript ID: JAV-03626)
