# NE Ratio Calculator

A phylogenetically informed framework for estimating effective population size (*Nₑ*) ratios across genomic compartments.
## Features
- Calculate nucleotide diversity (π) from VCF files- Estimate sequence divergence (D) using outgroups- Compute *Nₑ* ratios for autosomes, X, Y, and mitochondrial DNA- Phylogenetic calibration for mutation rate variation- Support for multiple species and populations
## Installation

```bash
pip install ne-ratio-calculator

git clone https://github.com/Ne-ratio/Ne_ratio.git
cd Ne_retio.git
pip install -e 

##Quick Start
python
from ne_calculator import NECalculator
# Initialize calculator
calculator = NECalculator()
# Calculate diversity from VCF
diversity_df = calculator.calculate_diversity("data/species.vcf.gz")
# Calculate divergence using outgroup
divergence_df = calculator.calculate_divergence("data/alignment.fasta", outgroup="macaque")
# Compute Nₑ ratios
ne_ratios = calculator.compute_ne_ratios(diversity_df, divergence_df)print(ne_ratios)

##Command Line Usage
bash
# Process VCF and calculate diversity
python scripts/process_vcf.py --vcf input.vcf.gz --out diversity.csv
# Run complete Nₑ analysis
python scripts/run_ne_analysis.py --vcf data.vcf --outgroup macaque --out results/

# Citation
Lu, L., Dai, W., Pan, Y., & Yan, Z. (in preparation). A Divergence-Calibrated Framework for Robust Cross-Species Comparison of Effective Population Size Ratios. Manuscript in preparation for Methods in Ecology and Evolution.

(The manuscript and associated tool are available at: https://github.com/Ne-ratio/Ne_ratio)

# Contact
For questions, please contact Zheng Yan (yanz@lzu.edu.cn).

