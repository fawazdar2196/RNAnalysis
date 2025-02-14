#!/usr/bin/env python3
"""
RNA Transcript Analysis Pipeline (Gene-Level Summary)

This pipeline processes a list of Ensembl gene IDs (from gene_list.json) and for each gene:
  - Queries Ensembl for transcript and exon metadata.
  - Retrieves transcript sequences (post-splicing, i.e., mature mRNA).
  - Computes:
      • Pre-splicing transcript length (from first exon start to last exon end).
      • Post-splicing transcript length (length of cDNA).
      • Exon statistics (number of exons and average exon length).
      • GC content of the mature transcript.
      • RNA secondary structure free energy using ViennaRNA's RNAfold and RNAstructure's Fold.
      • Internal priming (counts of A-rich regions, threshold = 6 consecutive A's).
  - Aggregates these values across all isoforms of the gene.

Results are saved to "transcript_analysis.csv".

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
from concurrent.futures import ThreadPoolExecutor, as_completed

# Set DATAPATH so RNAstructure finds its thermodynamic parameter files.
os.environ["DATAPATH"] = os.path.abspath("RNAstructure/data_tables")

# -----------------------------------------------------------------------------
# SET THESE VARIABLES:
RNAFOLD_EXE = "RNAfold"
# Set FOLD_EXE to the correct path of RNAstructure's Fold executable.
# In our case, it is located in the "exe" subdirectory and named "Fold".
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
    Timeout is set to 60 seconds.
    """
    try:
        proc = subprocess.run(
            [RNAFOLD_EXE],
            input=seq,
            capture_output=True,
            text=True,
            timeout=60
        )
        lines = proc.stdout.strip().splitlines()
        if len(lines) >= 2:
            line = lines[1]
            if '(' in line:
                match = re.search(r'\(([-\d\.]+)\)', line)
                if match:
                    energy_str = match.group(1)
                    if energy_str.replace('.', '').strip() == "":
                        return None
                    try:
                        return float(energy_str)
                    except ValueError:
                        return None
        return None
    except subprocess.TimeoutExpired:
        print("Error in RNAfold: Command timed out after 60 seconds")
        return None
    except Exception as e:
        print("Error in RNAfold:", e)
        return None

def analyze_structure_RNAstructure(seq):
    """
    Run RNAstructure's Fold tool on the sequence and extract the free energy.
    Uses a 60-second timeout. Reads output using utf-8 with error replacement to avoid decoding errors.
    """
    if not shutil.which(FOLD_EXE):
        print("RNAstructure Fold executable not found:", FOLD_EXE)
        return None
    try:
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as inf:
            inf.write(seq)
            inf_name = inf.name
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as outf:
            outf_name = outf.name

        proc = subprocess.run(
            [FOLD_EXE, inf_name, outf_name],
            capture_output=True,
            text=True,
            timeout=60
        )
        if proc.returncode != 0:
            print("RNAstructure Fold error:", proc.stderr)
            free_energy = None
        else:
            # Open the output file with error replacement for decoding issues
            with open(outf_name, "r", encoding="utf-8", errors="replace") as f:
                first_line = f.readline().strip()
                parts = first_line.split()
                if len(parts) >= 2:
                    try:
                        free_energy = float(parts[1])
                    except ValueError:
                        free_energy = None
                else:
                    free_energy = None
        os.remove(inf_name)
        os.remove(outf_name)
        return free_energy
    except subprocess.TimeoutExpired:
        print("Error in RNAstructure analysis: Command timed out after 60 seconds")
        return None
    except Exception as e:
        print("Error in RNAstructure analysis:", e)
        return None

def detect_internal_priming(seq, threshold=6):
    """Detect internal priming sites (regions with at least 'threshold' consecutive A's)."""
    return [m.start() for m in re.finditer(r'A{' + str(threshold) + r',}', seq, flags=re.IGNORECASE)]

def process_transcript(transcript):
    tid = transcript.get("id")
    exons = transcript.get("Exon", [])
    exon_lengths = []
    exon_starts = []
    exon_ends = []
    for exon in exons:
        try:
            s = int(exon.get("start", 0))
            e = int(exon.get("end", 0))
            exon_starts.append(s)
            exon_ends.append(e)
            exon_lengths.append(e - s + 1)
        except:
            continue
    pre_splice = max(exon_ends) - min(exon_starts) + 1 if exon_starts and exon_ends else None
    seq = fetch_transcript_sequence(tid)
    post_splice = len(seq) if seq else None
    gc = compute_gc_content(seq) if seq else None
    vRNA = analyze_structure_viennaRNA(seq) if seq else None
    rStruct = analyze_structure_RNAstructure(seq) if seq else None
    priming = detect_internal_priming(seq) if seq else []
    return {
        "transcript_id": tid,
        "pre_splice_length": pre_splice,
        "post_splice_length": post_splice,
        "exon_count": len(exon_lengths),
        "avg_exon_length": sum(exon_lengths)/len(exon_lengths) if exon_lengths else None,
        "gc_content": gc,
        "viennaRNA_energy": vRNA,
        "RNAstructure_energy": rStruct,
        "internal_priming_count": len(priming)
    }

def process_gene(gene_id):
    gene_data = fetch_gene_data(gene_id)
    if gene_data is None:
        return None
    transcripts = gene_data.get("Transcript", [])
    results = []
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = {executor.submit(process_transcript, t): t for t in transcripts}
        for future in as_completed(futures):
            res = future.result()
            if res:
                results.append(res)
    if not results:
        return None
    def avg(vals):
        valid = [v for v in vals if v is not None]
        return sum(valid)/len(valid) if valid else None
    return {
        "gene_id": gene_id,
        "num_transcripts": len(results),
        "avg_pre_splice_length": avg([r["pre_splice_length"] for r in results]),
        "avg_post_splice_length": avg([r["post_splice_length"] for r in results]),
        "avg_exon_count": avg([r["exon_count"] for r in results]),
        "avg_exon_length": avg([r["avg_exon_length"] for r in results]),
        "avg_gc_content": avg([r["gc_content"] for r in results]),
        "avg_viennaRNA_energy": avg([r["viennaRNA_energy"] for r in results]),
        "avg_RNAstructure_energy": avg([r["RNAstructure_energy"] for r in results]),
        "avg_internal_priming_count": avg([r["internal_priming_count"] for r in results])
    }

def main():
    with open("gene_list.json", "r") as f:
        gene_list = json.load(f)
    summaries = []
    for gid in gene_list:
        print(f"Processing gene: {gid}")
        summ = process_gene(gid)
        if summ:
            summaries.append(summ)
    if summaries:
        df = pd.DataFrame(summaries)
        df.to_csv("transcript_analysis.csv", index=False)
        print("Analysis complete. Output saved to transcript_analysis.csv")
        print(df)
    else:
        print("No gene data processed.")

if __name__ == "__main__":
    main()
EOF
