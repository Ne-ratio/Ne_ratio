# Ne Ratio Calculator

A computational framework for estimating effective population size (Nₑ) ratios across genomic compartments (autosomes, X chromosome, Y chromosome, mitochondrial DNA) from VCF files.

## Description

This tool implements the phylogenetic framework described in Lu et al. for robust estimation of relative effective population size ratios. The method uses interspecies divergence to calibrate for mutation rate variation, enabling comparable Nₑ ratio estimates across species.

**Key Features:**
- Calculate nucleotide diversity (π) for genomic compartments
- Estimate sequence divergence (D) using outgroup species
- Compute Nₑ ratios with mutation rate calibration
- Support for sliding window analysis
- Bootstrap confidence intervals
- Memory-efficient VCF processing

## Installation

```bash
git clone https://github.com/yourusername/ne_ratio_calculator.git
cd ne_ratio_calculator
pip install -e .