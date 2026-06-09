# AlphaMicrobiome

AlphaMicrobiome is a Python‑based toolkit for microbiome data analysis, designed to deliver efficient analytical solutions. It integrates a variety of microbiome analysis algorithms and supports end‑to‑end analysis workflows, ranging from raw data preprocessing to advanced statistical analysis.

---

# Features

- **End‑to‑end workflows** – Support for raw sequencing data (FASTQ, FASTA), OTU/ASV tables, and metadata integration.
- **Alpha diversity estimation** – Shannon, Simpson, Chao1, ACE, Faith’s PD, and more.
- **Beta diversity analysis** – UniFrac (weighted/unweighted), Bray‑Curtis, Jaccard.
- **Taxonomic composition** – Phylum‑to‑genus level summarization with interactive bar plots.
- **Statistical comparisons** – PERMANOVA, ANOSIM, DESeq2‑like differential abundance, LEfSe wrapper.
- **Visualization** – PCA, PCoA, NMDS, heatmaps, rarefaction curves.
- **Python‑first design** – Seamless integration with `pandas`, `scikit‑bio`, `numpy`, and `matplotlib`.

---

# Installation

Download and install from GitHub:
```bash
git clone https://github.com/zhaoy2020/AlphaMicrobiome.git
cd AlphaMicrobiome
pip install -e .
```

Install directly from PyPi:
```bash
pip install AlphaMicrobiome
```

# Tutorial
Detailed tutorials can be found in the notebooks.

# Contributing
We welcome contributions! Please follow these steps:
1. Fork the repository and clone it locally.
2. Create a new branch: git checkout -b feature/your-feature-name.
3. Install development dependencies: pip install -e .[dev].
4. Write tests under tests/ and ensure all tests pass: pytest tests/.
5. Update documentation and examples as needed.
6. Push to your fork and open a Pull Request.