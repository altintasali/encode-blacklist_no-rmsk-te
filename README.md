# encode-blacklist_no-rmsk-te

A highly reproducible genomic data pipeline to dissect, analyze, and purge RepeatMasker features and Transposable Elements (TEs) from the ENCODE mm10 v2 genomic blacklist regions. 

The workspace leverages a master bash orchestrator for automated parallel spatial subtraction (`bedtools`) paired with an R statistical engine computing structural enrichment profiles via two-tailed Fisher's Exact Tests.

---

## 📂 Repository Architecture

```text
ENCODE_Blacklist_without_RepeatMasker-TE/
├── .gitignore
├── README.md
├── run_pipeline.sh            # Master orchestrator script
├── scripts/
│   ├── filter_blacklist.sh    # Core spatial genomics (bedtools processing)
│   └── analyze_repeats.R      # R statistical framework & data visualization
├── data/
│   ├── raw/                   # Input genomic raw source coordinates
│   │   ├── mm10-blacklist.v2.bed
│   │   └── rmsk.txt.gz
│   └── process_temp/          # Clean intermediate alignment BED tracking
│       ├── mm10_all-rmsk.bed
│       ├── mm10_only-TEs.bed
│       ├── mm10_all-rmsk_trapped_in_blacklist.bed
│       └── mm10_TEs_trapped_in_blacklist.bed
├── results/
│   ├── cleaned_blacklist/     # Final, rescued blacklists (Outputs)
│   │   ├── mm10-blacklist.v2_no-TEs.bed
│   │   └── mm10-blacklist.v2_no-all-rmsk.bed
│   └── enrichment/            # Analytical summaries and diagnostic plots
│       ├── enrichment_results.tsv.gz     # Gzipped statistical calculations
│       ├── mm10_TE_pool_horizontal.png   # Comparative composition bars
│       └── repeat_blacklist_forest_plot.png # Log-scaled forest plot
└── logs/                      # Execution tracking and environment diagnostics
    ├── pipeline_execution.log # Full stdout/stderr master execution log
    └── r_analysis.log         # Dedicated R environment metrics log
```

## 📊 Pipeline Workflow

1. **Orchestration**: `run_pipeline.sh` manages environment steps, guarantees directory setups, checks exits, and safely pipes global stream logs to files.
2. **Genomic Processing (`scripts/filter_blacklist.sh`)**:
   - Compiles raw UCSC RepeatMasker streams into standardized 8-column BED sets.
   - Isolates pure Transposable Elements (TE) from background features.
   - Leverages `bedtools intersect` to map items trapped inside the blacklists.
   - Generates two rescued versions via `bedtools subtract` saved directly to `results/cleaned_blacklist/`: one removing core TEs, and one dropping all repetitive elements completely.
3. **Statistical Diagnostics (`scripts/analyze_repeats.R`)**:
   - Tests structural enrichment by executing iterative two-tailed Fisher’s Exact Tests per feature class.
   - Controls Type I errors across multi-testing rows using Benjamini-Hochberg ($p$-adjustment).
   - Generates production-ready structural profiles without trailing junk graphics (`Rplots.pdf` generation is suppressed).

## 📦 Requirements

### Command Line Binaries
* `bash` (>= 4.0)
* `wget`
* `bedtools` (>= v2.30.0)
* `awk`

### R Packages
* `R` (>= 4.1)
* `data.table` 
* `ggplot2`
* `scales`
* `magrittr`
* `pacman` (Handles automated missing-dependency deployment)

---

## 🚀 Execution

To trigger the full spatial subtraction, run statistical profiling, export files, and render data visualization from scratch, execute the top-level script from the repository root:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

## 📥 Data Outputs

All metrics and cleaned outputs land directly inside `results/`:
* **`results/cleaned_blacklist/`**: Contains your custom-filtered, ready-to-use background files with repetitive structures subtracted out.
* **`results/enrichment/enrichment_results.tsv.gz`**: Tab-separated, gzipped spreadsheet containing detailed statistical summaries (Odds Ratios, Confidence Intervals, and adjusted $p$-values) sorted by enrichment intensity.
* **`results/enrichment/mm10_TE_pool_horizontal.png`**: Side-by-side faceted composition plot showing repeat class ratios globally vs. inside the blacklist.
* **`results/enrichment/repeat_blacklist_forest_plot.png`**: A log-scaled, color-coded forest plot highlighting significant enrichments with explicit text annotations (`OR = ...`).

## 👤 Author
* **Ali Altintas** * Date: July 2026
