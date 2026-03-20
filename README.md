# CdiffFragR

R tools for *Clostridioides difficile* fluorescent PCR ribotyping. Matches query chromatograms (`.fsa` format) generated on capillary electrophoresis sequencers to a database of known *C. difficile* ribotypes.

Developed and maintained by the [Walk Lab](https://www.walklab.org) at Montana State University, Bozeman. Compatible with Mac OS X and Windows 10.

---

## Installation

### 1. R
Download and install R (≥ 4.0): https://cran.r-project.org

### 2. Rtools *(Windows only)*
Download Rtools (≥ 4.0): https://cran.r-project.org/bin/windows/Rtools/

### 3. RStudio
Download RStudio Desktop (≥ 1.4): https://www.rstudio.com/products/rstudio/download/

### 4. CdiffFragR
Click the **Code** button on this page and select **Download ZIP**. Extract the contents to a directory on your computer.

---

## Usage

1. Place `.fsa` files to analyze in the `Files_to_analyze` directory.
2. Open `CdiffFragR.Rproj` in RStudio.
3. Open `Call_FSA.R` via **File > Open File…**
4. Highlight all text and click **Run**.

After processing, each `.fsa` file is moved to `Files_analyzed/`. Results are written to a timestamped directory (`Results_YYYY.MM.DD-Hour/`) containing:

- `chrom_*.jpeg` — chromatogram summaries
- `hits_*.jpeg` — overlay plots comparing each query to its closest database match
- `SUMMARY.csv` — results table

---

## Interpreting Results

Open `SUMMARY.csv` and review the `Dist_1` column, which reports the Bray-Curtis dissimilarity between the query and the closest database match (0 = identical, 1 = completely dissimilar).

| Dist_1 | Confidence |
|---|---|
| < 0.10 | Good match |
| 0.10 – 0.20 | Questionable match |
| > 0.20 | Poor match |

**Good match** — Solid match with a very low chance of a false positive.

**Questionable match** — May or may not be a match. Open the hit image and check whether peaks align (peaks within 1–2 bp should be considered the same).

**Poor match** — Not a match to the database. Check the chromatogram: either the ribotype is not yet in the database, or the sample should be re-analyzed.

### New Ribotypes
If you believe you have identified a new ribotype, re-analyze to confirm, then contact the Walk Lab — we welcome additions to the database.

---

## Manual

See [`cdifffragr_12.06.pdf`](cdifffragr_12.06.pdf) for the full user manual.

---

## License

[GPL-3.0](LICENSE)
