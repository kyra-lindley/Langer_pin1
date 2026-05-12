# Langer_pin1

Code supporting the analysis in:

> Bowman CL, Daniel CJ, Carlson EJ, Shah VM, Farrell AS, Kresse KM, Wang X, **Lindley KA**, et al. **The Prolyl Isomerase PIN1 Affects Fibroblast Differentiation States and Cross-talk in Pancreatic Cancer.** *Cancer Res.* 2025;85(24):4899–4917. https://doi.org/10.1158/0008-5472.CAN-24-3437

This repository contains shell scripts, R analysis code, and ImageJ macros used to process and analyze bulk RNA-seq data from cancer-associated fibroblast (CAF) and pancreatic stellate cell (PSC) models of PIN1 knockdown in pancreatic ductal adenocarcinoma (PDAC).

---

## Repository Structure

| File | Language | Description |
|------|----------|-------------|
| `trim.sh` | Bash | Adapter trimming with Trimmomatic (paired-end) |
| `align.sh` | Bash | SLURM array job for STAR alignment to human genome |
| `feature_counts.sh` | Bash | Read counting with featureCounts |
| `deseq_init.R` | R | DESeq2 object initialization, rlog normalization, and results extraction |
| `sig_deg_graph.R` | R | Visualization of significant differentially expressed genes |
| `genelist_heatmap` | R | Heatmap generation from curated gene lists |
| `CHP.ijm` | ImageJ Macro | Image analysis macro |
| `IHC_pMET_batch_threshold_analysis.ijm` | ImageJ Macro | Batch IHC thresholding and quantification for pMET staining |

---

## Workflow Overview

Bulk RNA-seq data were processed following this general order:

```
Raw FASTQ
    │
    ▼
trim.sh          # Trimmomatic: adapter removal, quality filtering
    │
    ▼
align.sh         # STAR: alignment to human reference genome (GRCh38)
    │
    ▼
feature_counts.sh  # featureCounts: gene-level read quantification
    │
    ▼
deseq_init.R     # DESeq2: normalization, DE analysis, rlog transformation
    │
    ▼
sig_deg_graph.R / genelist_heatmap  # Visualization
```

---

## Dependencies

**Shell scripts** (run on OHSU ARC/Exacloud HPC via SLURM):
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [STAR](https://github.com/alexdobin/STAR)
- [Subread/featureCounts](https://subread.sourceforge.net/)

**R packages:**
- `DESeq2`
- `dplyr`
- `readxl`
- `AnnotationDbi`
- `org.Hs.eg.db`

**Image analysis:**
- [ImageJ/Fiji](https://fiji.sc/)

---

## Usage Notes

- `align.sh` is configured as a SLURM array job. Update `--array` to match your sample count and set paths for `GENOME_DIR`, `FASTQ_DIR`, and `OUT_DIR` before running.
- `trim.sh` expects paired-end reads named `*_R1.fastq.gz` / `*_R2.fastq.gz`. Update `ADAPTERS` and directory paths as needed.
- `deseq_init.R` requires a counts matrix and a metadata `.xlsx` file. Update the `contrast` vector and output paths to match your comparison of interest. The script saves DESeq2 (`dds`), results (`res`), and rlog (`rld`) objects as `.rds` files for downstream use.
- Sample IDs in `deseq_init.R` follow the naming convention used in the paper (e.g., `CAFPINT1` = CAF, PIN1 knockdown, TGF-β treated, replicate 1).

---

## Contact

Kyra Lindley — kyra.a.lindley@gmail.com  
Langer Lab, CEDAR, Oregon Health and Science University
