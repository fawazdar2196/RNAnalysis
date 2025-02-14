#!/usr/bin/env python3
"""
RNA Transcript Analysis Pipeline (Gene-Level Summary)

This pipeline processes a list of Ensembl gene IDs (from gene_list.json) and for each gene:
  - Queries Ensembl for transcript and exon metadata.
  - Retrieves transcript sequences (post-splicing, i.e., mature mRNA).
  - Computes:
      • Pre-splicing transcript length (genomic span: from first exon start to last exon end).
      • Post-splicing transcript length (length of cDNA).
      • Exon statistics: number of exons and average exon length.
      • GC content of the mature transcript.
      • RNA secondary structure free energy using ViennaRNA's RNAfold and RNAstructure's Fold.
      • Internal priming: counts of A-rich regions (6+ consecutive A's).
  - Aggregates these values across all isoforms of the gene.

The final results are saved to "transcript_analysis.csv".

Usage:
    python3 rna_transcript_pipeline.py
"""

import os
import json
import requests
import pandas as pd
import subprocess
import re
import tempfile
import time
import shutil

# -----------------------------------------------------------------------------
# Set DATAPATH so that RNAstructure finds its thermodynamic parameter files.
# Assuming that the data_tables directory is located inside the RNAstructure folder.
os.environ["DATAPATH"] = os.path.abspath("RNAstructure/data_tables")

# -----------------------------------------------------------------------------
# SET THESE VARIABLES:
#
# RNAfold is installed via Ubuntu’s vienna-rna package.
RNAFOLD_EXE = "RNAfold"

# Set FOLD_EXE to the correct path of RNAstructure's Fold executable.
# Based on our build, the executable is located in the "exe" subdirectory and named "Fold".
FOLD_EXE = "RNAstructure/exe/Fold"
# -----------------------------------------------------------------------------

def fetch_gene_data(gene_id):
    """Query Ensembl's REST API for gene info (with transcript details)."""
    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    response = requests.get(url, headers=headers)
    if response.ok:
        return response.json()
    else:
        print(f"Error fetching gene {gene_id}: {response.text}")
        return None

def fetch_transcript_sequence(transcript_id):
    """Fetch the mature (cDNA) sequence for the given transcript ID."""
    url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cdna"
    headers = {"Content-Type": "text/plain", "Accept": "text/plain"}
    response = requests.get(url, headers=headers)
    if response.ok:
        return response.text.strip()
    else:
        print(f"Error fetching transcript sequence for {transcript_id}: {response.text}")
        return ""

def compute_gc_content(seq):
    """Compute the GC content percentage of a sequence."""
    if not seq:
        return 0
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    return (gc_count / len(seq)) * 100

def analyze_structure_viennaRNA(seq):
    """
    Run RNAfold on the sequence and extract the free energy.
    If RNAfold output is non-numeric (e.g. "......"), return None.
    """
    try:
        process = subprocess.run(
            [RNAFOLD_EXE],
            input=seq,
            capture_output=True,
            text=True,
            check=True
        )
        lines = process.stdout.strip().splitlines()
        if len(lines) >= 2:
            line = lines[1]
            if '(' in line:
                match = re.search(r'\(([-\d\.]+)\)', line)
                if match:
                    energy_str = match.group(1)
                    if energy_str.replace('.', '').strip() == "":
                        print("Error: Extracted energy is non-numeric:", energy_str)
                        return None
                    try:
                        return float(energy_str)
                    except ValueError:
                        print("Error: Could not convert extracted value to float:", energy_str)
                        return None
                else:
                    print("Error: Free energy pattern not found in RNAfold output:", line)
                    return None
            else:
                print("RNAfold output does not contain free energy info:", line)
                return None
        return None
    except Exception as e:
        print("Error in ViennaRNA analysis:", e)
        return None

def analyze_structure_RNAstructure(seq):
    """
    Run RNAstructure's Fold tool on the sequence and extract the free energy.
    If the executable is not found or if there is an error reading parameters, return None.
    """
    if not shutil.which(FOLD_EXE):
        print("RNAstructure Fold executable not found:", FOLD_EXE)
        return None
    try:
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as infile:
            infile.write(seq)
            infile_name = infile.name
        
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as outfile:
            outfile_name = outfile.name

        process = subprocess.run(
            [FOLD_EXE, infile_name, outfile_name],
            capture_output=True,
            text=True
        )
        if process.returncode != 0:
            print("RNAstructure Fold error:", process.stderr)
            free_energy = None
        else:
            with open(outfile_name, "r") as f:
                first_line = f.readline().strip()
                parts = first_line.split()
                if len(parts) >= 2:
                    try:
                        free_energy = float(parts[1])
                    except ValueError:
                        free_energy = None
                else:
                    free_energy = None
        os.remove(infile_name)
        os.remove(outfile_name)
        return free_energy
    except Exception as e:
        print("Error in RNAstructure analysis:", e)
        return None

def detect_internal_priming(seq, threshold=6):
    """
    Detect internal priming sites (regions with at least 'threshold' consecutive A's).
    Returns a list of starting positions (0-indexed).
    """
    return [m.start() for m in re.finditer(r'A{' + str(threshold) + r',}', seq, flags=re.IGNORECASE)]

def process_gene(gene_id):
    """
    Process one gene: for each transcript, extract features and aggregate gene-level metrics.
    """
    gene_data = fetch_gene_data(gene_id)
    if gene_data is None:
        return None

    transcripts = gene_data.get("Transcript", [])
    transcript_metrics = []

    for transcript in transcripts:
        transcript_id = transcript.get("id")
        exon_list = transcript.get("Exon", [])
        exon_lengths = []
        exon_starts = []
        exon_ends = []
        for exon in exon_list:
            try:
                start = int(exon.get("start", 0))
                end = int(exon.get("end", 0))
                exon_starts.append(start)
                exon_ends.append(end)
                exon_lengths.append(end - start + 1)
            except Exception:
                continue
        exon_count = len(exon_lengths)
        avg_exon_length = sum(exon_lengths) / exon_count if exon_count > 0 else None

        pre_splice_length = None
        if exon_starts and exon_ends:
            pre_splice_length = max(exon_ends) - min(exon_starts) + 1

        sequence = fetch_transcript_sequence(transcript_id)
        post_splice_length = len(sequence) if sequence else None
        gc_content = compute_gc_content(sequence) if sequence else None
        viennaRNA_energy = analyze_structure_viennaRNA(sequence) if sequence else None
        RNAstructure_energy = analyze_structure_RNAstructure(sequence) if sequence else None
        priming_sites = detect_internal_priming(sequence, threshold=6) if sequence else []
        internal_priming_count = len(priming_sites)

        transcript_metrics.append({
            "transcript_id": transcript_id,
            "pre_splice_length": pre_splice_length,
            "post_splice_length": post_splice_length,
            "exon_count": exon_count,
            "avg_exon_length": avg_exon_length,
            "gc_content": gc_content,
            "viennaRNA_energy": viennaRNA_energy,
            "RNAstructure_energy": RNAstructure_energy,
            "internal_priming_count": internal_priming_count,
        })
        time.sleep(0.1)

    num_transcripts = len(transcript_metrics)
    if num_transcripts == 0:
        return None

    def avg(values):
        valid = [v for v in values if v is not None]
        return sum(valid) / len(valid) if valid else None

    gene_summary = {
        "gene_id": gene_id,
        "num_transcripts": num_transcripts,
        "avg_pre_splice_length": avg([m["pre_splice_length"] for m in transcript_metrics]),
        "avg_post_splice_length": avg([m["post_splice_length"] for m in transcript_metrics]),
        "avg_exon_count": avg([m["exon_count"] for m in transcript_metrics if m["exon_count"] is not None]),
        "avg_exon_length": avg([m["avg_exon_length"] for m in transcript_metrics if m["avg_exon_length"] is not None]),
        "avg_gc_content": avg([m["gc_content"] for m in transcript_metrics if m["gc_content"] is not None]),
        "avg_viennaRNA_energy": avg([m["viennaRNA_energy"] for m in transcript_metrics if m["viennaRNA_energy"] is not None]),
        "avg_RNAstructure_energy": avg([m["RNAstructure_energy"] for m in transcript_metrics if m["RNAstructure_energy"] is not None]),
        "avg_internal_priming_count": avg([m["internal_priming_count"] for m in transcript_metrics if m["internal_priming_count"] is not None]),
    }
    return gene_summary

def main():
    try:
        with open("gene_list.json", "r") as f:
            gene_list = json.load(f)
    except Exception as e:
        print("Error loading gene_list.json:", e)
        return

    gene_summaries = []
    for gene_id in gene_list:
        print(f"Processing gene: {gene_id}")
        summary = process_gene(gene_id)
        if summary is not None:
            gene_summaries.append(summary)

    if gene_summaries:
        df = pd.DataFrame(gene_summaries)
        df.to_csv("transcript_analysis.csv", index=False)
        print("Analysis complete. Output saved to transcript_analysis.csv")
        print(df)
    else:
        print("No gene data processed.")

if __name__ == "__main__":
    main()
