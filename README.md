# RNAnalysis

RNAnalysis is an RNA Transcript Analysis Pipeline that retrieves gene and transcript metadata from Ensembl, computes transcript metrics, and predicts RNA secondary structure using external tools (ViennaRNA's RNAfold and RNAstructure's Fold). The results are aggregated into a CSV file for further analysis.

## Features

- Reads gene IDs from a JSON file.
- Retrieves transcript and exon data from Ensembl's REST API.
- Fetches mature transcript (cDNA) sequences.
- Computes:
  - Pre-splicing (genomic span) and post-splicing transcript lengths.
  - Exon count and average exon length.
  - GC content of transcripts.
- Predicts RNA secondary structure free energy using:
  - **ViennaRNA (RNAfold)**
  - **RNAstructure (Fold)**
- Detects internal priming sites (6+ consecutive A’s).
- Aggregates data across all transcripts of a gene.
- Implements concurrent processing for improved performance.

## Requirements

- Python 3.6+
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) (installed via `sudo apt install vienna-rna`)
- [RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html) (built from source)
- Required Python packages: `requests`, `pandas`

## Setup

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/fawazdar2196/RNAnalysis.git
   cd RNAnalysis
