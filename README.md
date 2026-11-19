<div align="center">

![visitors](https://visitor-badge.laobi.icu/badge?page_id=zhaoy2020/AlphaMicrobiome)
[![GitHub Repo stars](https://img.shields.io/github/stars/zhaoy2020/AlphaMicrobiome?style=social)](https://github.com/zhaoy2020/AlphaMicrobiome/stargazers)
[![GitHub Code License](https://img.shields.io/github/license/zhaoy2020/AlphaMicrobiome)](LICENSE)
[![GitHub pull request](https://img.shields.io/badge/PRs-welcome-blue)](https://github.com/zhaoy2020/AlphaMicrobiome/pulls)

</div>

# AlphaMicrobiome

AlphaMicrobiome is an integrated Python and R toolkit for microbiome data analysis, designed to provide efficient and reproducible analytical solutions. It supports workflows from microbiome data preprocessing and diversity analysis to statistical testing, ecological analysis, machine learning, and visualization.

---

# Features

- **End‑to‑end workflows** – Support for raw sequencing data (FASTQ, FASTA), OTU/ASV tables, and metadata integration.
- **Alpha diversity estimation** – Shannon, Simpson, Chao1, ACE, Faith’s PD, and more.
- **Beta diversity analysis** – UniFrac (weighted/unweighted), Bray‑Curtis, Jaccard.
- **Taxonomic composition** – Phylum‑to‑genus level summarization with interactive bar plots.
- **Statistical comparisons** – PERMANOVA, ANOSIM, DESeq2‑like differential abundance, LEfSe wrapper.
- **Visualization** – PCA, PCoA, NMDS, heatmaps, rarefaction curves.
- **ython and R integration** – Combined use of Python and R for microbiome and ecological analyses.

---

# Study-specific code
The repository includes the Python and R scripts used for the analyses and figure generation in our potato soil microbiome study. These scripts cover major analyses including diversity, community composition, environmental associations, community assembly, microbial networks, functional prediction, and figure generation.

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
Detailed tutorials and examples can be found in the notebooks and analysis scripts.

# Contributing
We welcome contributions! Please follow these steps:
1. Fork the repository and clone it locally.
2. Create a new branch: git checkout -b feature/your-feature-name.
3. Install development dependencies: pip install -e .[dev].
4. Write tests under tests/ and ensure all tests pass: pytest tests/.
5. Update documentation and examples as needed.
6. Push to your fork and open a Pull Request.