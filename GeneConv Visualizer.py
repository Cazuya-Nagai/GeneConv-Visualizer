"""
GeneConv Visualizer: Mapping Concerted Evolution in Mitogenomes
--------------------------------------------------------------
Author: Cazuya  (2026)
Description: Automatically parses GENECONV output (PI/PO/GI) and 
             visualizes gene conversion coverage with automated statistics.
"""

import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

class GeneConvVisualizer:
    def __init__(self, dup_start, dup_end, tail_len=750):
        self.dup_range = (dup_start, dup_end)
        self.tail_len = tail_len
        self.df = None
        self.counts = None

    def load_data(self, file_path):
        """Reads the GENECONV result file and extracts fragment information."""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Error: File '{file_path}' not found.")

        rows = []
        with open(file_path, 'r') as f:
            for line in f:
                # Extract only lines starting with PI, PO, or GI
                if not re.match(r'^(PI|PO|GI)', line):
                    continue
                parts = re.split(r'\s+', line.strip())
                if len(parts) < 7: continue

                try:
                    # Replace GENECONV's 'D' with 'e' for scientific notation
                    sim_p = float(parts[2].replace("D", "e").replace("d", "e"))
                    begin, end, length = int(parts[4]), int(parts[5]), int(parts[6])
                    rows.append({"SimP": sim_p, "Begin": begin, "End": end, "Length": length})
                except (ValueError, IndexError):
                    continue
        
        if not rows:
            raise ValueError("Error: No valid GENECONV fragment data found in the file.")
            
        self.df = pd.DataFrame(rows)

    def calculate_coverage(self, p_threshold=0.05, min_length=500):
        """Filters significant tracts and generates the coverage array."""
        filtered = self.df[(self.df["SimP"] < p_threshold) & (self.df["Length"] >= min_length)]
        
        if filtered.empty:
            print("Warning: No fragments met the significance criteria.")
            self.counts = np.zeros(1000) 
            return filtered

        max_pos = self.df["End"].max()
        # Allocate array with padding
        self.counts = np.zeros(max_pos + 500, dtype=int)
        for _, row in filtered.iterrows():
            self.counts[row["Begin"]:row["End"]+1] += 1
        
        return filtered

    def plot_and_save(self, output_name):
        """Generates the figure and automatically calculates statistics for the legend."""
        # Calculate statistics for the legend based on the current dataset
        total_sum = self.counts.sum()
        dup_sum = self.counts[self.dup_range[0]:self.dup_range[1]+1].sum()
        tail_start = self.dup_range[1] - self.tail_len
        tail_sum = self.counts[tail_start:self.dup_range[1]+1].sum()

        dup_pct = int(round(dup_sum / total_sum * 100)) if total_sum > 0 else 0
        tail_pct = int(round(tail_sum / total_sum * 100)) if total_sum > 0 else 0

        # Plotting
        plt.figure(figsize=(14, 4.5))
        x = np.arange(len(self.counts))
        
        # Draw coverage as a step function for scientific accuracy
        plt.plot(x, self.counts, color='blue', linewidth=1.2, 
                 drawstyle='steps-post', label='Fragment coverage')
        
        # Highlight regions
        plt.axvspan(self.dup_range[0], self.dup_range[1], color='#ffe0b2', alpha=0.6, 
                    label=f'Duplication region ({dup_pct}%)', zorder=0)
        plt.axvspan(tail_start, self.dup_range[1], color='#ffab91', alpha=0.7, 
                    label=f'Last {self.tail_len} bp of duplication ({tail_pct}%)', zorder=1)

        # Labels and Style
        plt.title("Significant Gene Conversion Fragment Coverage", fontsize=15)
        plt.xlabel("Mitochondrial DNA Position (bp)", fontsize=12)
        plt.ylabel("Stacking Coverage (Fragments)", fontsize=12)
        plt.grid(True, linestyle='-', alpha=0.3)
        plt.legend(loc='upper left', frameon=True)
        
        plt.tight_layout()
        plt.savefig(output_name, dpi=300)
        print(f"Success: Figure saved as {output_name}")

def main():
    parser = argparse.ArgumentParser(description="Visualize GENECONV coverage hotspots.")
    parser.add_argument("input", help="Path to raw GENECONV output file")
    parser.add_argument("--start", type=int, required=True, help="Duplication start coordinate")
    parser.add_argument("--end", type=int, required=True, help="Duplication end coordinate")
    parser.add_argument("--out", default="geneconv_plot.png", help="Output image name")
    parser.add_argument("--p", type=float, default=0.05, help="P-value threshold")
    parser.add_argument("--minlen", type=int, default=500, help="Min fragment length")

    args = parser.parse_args()

    try:
        viz = GeneConvVisualizer(args.start, args.end)
        viz.load_data(args.input)
        viz.calculate_coverage(p_threshold=args.p, min_length=args.minlen)
        viz.plot_and_save(args.out)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()