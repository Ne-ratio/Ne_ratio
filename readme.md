# NE Ratio Calculator

A phylogenetically informed framework for estimating effective population size (*Nₑ*) ratios across genomic compartments.
## Features
- Calculate nucleotide diversity (π) from VCF files- Estimate sequence divergence (D) using outgroups- Compute *Nₑ* ratios for autosomes, X, Y, and mitochondrial DNA- Phylogenetic calibration for mutation rate variation- Support for multiple species and populations
- All samples (e.g., human, bonobo, chimp) need to be mapped to the outgroup genome in order to calculate divergence.
## Installation and dependencies

There is no system installation required. Just download this entire repository using the green "Code" button at the top of this page, or with the bash command `git clone https://github.com/Ne-ratio/Ne_ratio.`. If you prefer to move scripts to the directory where you run them, bear in mind that many of these scripts require the script `genomics.py` to be present in the same directory (or on your PYTHONPATH). The easiest is just to leave them in the main directory and specify the full path to the script when running it from elsewhere.

The only dependency is [numpy](https://numpy.org/).

Most of the scripts now run in python 3. Some are still written for python 2, but those will be updated soon.

___

## Processing VCF files

Most of my scripts use a processed `.vcf` format that I call `.geno`. This looks something like this:

```
#CHROM      POS      ind1      ind2      ind3
scaffold1   1        A/A       A/G       G|A
scaffold1   1        N/N       T/T       T|C
```

Missing data is denoted as `N`, and phased and unphased genotypes are shown conventionally with `|` and `/`.

The script `parseVCF.py` in the [`VCF_processing`](https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing) directory, will convert vcf to this format. It has various options for filtering based on read depth, genotype quality or any other flag in the `FORMAT` column of the vcf.

#### Example command

```bash
python VCF_processing/parseVCF.py -i input.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 max=50 -o output.geno.gz
```

You can read more about this script in the [`VCF_processing`](https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing) directory.

---

## Filtering genotype files prior to further analysis

If you vcf file was not already filtered, or you would like to filter further. The script `filtergenotypes.py` has many options for filtering. Some examples include:
* Number of individuals with non-missing `N/N` genotypes at a site (`--minCalls`)
* Number of alleles observed at a site across all individuals (`--minAlleles` and `--maxAlleles`)
* Minor allele count (`--minVarCount`)
* Distance between sites for thinning (`--thinDist`)

It requires the script `genomics.py` to be present in the same directory, or in your Python path.

#### Example command
```bash
python filterGenotypes.py --threads 4 -i input.geno.gz -o output.geno.gz --minAlleles 2 --minCalls 10 --thinDist 1000
```

#### Notes

`python filterGenotypes.py -h` will print a full list of command options.

By default, this script assumes that the `.geno` input file is is encoded as diploid genotypes with a phase operator (`/` or `|`):
```
scaffold1  1        A/A       G/G       G|A
``` 
You can specify a different input formats using the `-if`, but this is not recommended.

You can also specify various putput formats using `-of`.

| Output format | Description | Example |
| :-----------: | ----------- | -------- |
| `phased` (default) | Alleles separates by a phase operator. This doesn't mean the phase is known, just that it is indicated | `A/A    G/G    G\|A` |
| `diplo`     | For diploids only. Genotypes are single bases denoting the diploid genotype, using ambiguity codes for heterozygotes | `A       G       R` |
| `alleles`  | as above but without the phase operator | `AA    GG    GA` |
| `randomAllele` | Randomly pick one allele per individual | `A    G    A` |
| `coded` | Coded numerically as in the VCF | `0/0    1/1    1\|0` |
| `bases` | Separate the alleles for each individual into different columns and also give different headers for each | `A    A    G    G    G    A` |
| `counts` | Count of the minor allele in each individual | `0    2    1` |

___

## Diversity and divergence analyses in sliding windows

The script `popgenWindows.py` computes some standard population genomic statistics in sliding windows:  *pi*,  *F*<sub>ST</sub> and *D<sub>XY</sub>*. It requires the script `genomics.py` to be present in the same directory, or in your Python path.

#### Example command
```bash
python popgenWindows.py -w 50000 -m 5000 -g input.geno.gz -o output.csv.gz -f phased -T 5 -p popA A1,A2,A3,A4 -p popB B1,B2,B3,B4,B6,B6 -p popC -p popD --popsFile pops.txt 
```

#### Notes

`python popgenWindows.py -h` Will print a full list of command arguments.

* Input is a `.geno` file as shown above. This can be gzipped (`.geno.gz`).
Output is a `.csv`. If you add `.gz` it will be gzipped.

* Genotype encoding is indicated by the `-f` flag. `-f phased` is normally used, see the table above. Other options are `-f haplo` for haploid data (although `phased` will also interpret haploid data correctly), `-f diplo` and `f pairs` which is like the `alleles` output in the table above, but assumes diplid data.

* There are three options for defining the windows to be analysed, using the `--windType` argument.

| Wndow Type | Description |
| :-----------: | ----------- | 
| `coordinate`     | Windows will cover a fixed range in the genome, which is defined as the window size. If there is missing data, this can lead to variable numbers of sites used for each window. |
| `sites`        | Each window will have the same number of sites. If there is missing data, this can lead to different absolute sizes for windows in terms of genome coordinates. |
| `predefined`          | This will analyse predefined windows provided using the `--windCoords` flag. |

* You can either include sample names after the population name, separated by commas, or provide only the population name, along with a populations file, with the flag `--popsFile `, which has two columns: the first gives sample names and teh second gives population name:

```
C1  popC
C2  popC
C3  popC
C4  popC
D1  popD
D2  popD
D3  popD
D4  popD
```

* The most common source of errors here involve the `-m` (`--minSites`) flag. If you are useing coordinate windows and have any sites with missing data, then `-m` must be set to a value smaller than the window size. If you have reduced representation data such as RADseq, you will need a much lower `-m` value (more like 1% of the window size or even less).

* If some samples are haploid and others are diploid, you can use one of the diploid formats, but indicate that certain samples are haploid by listing them after the `--haploid` flag. The script will force them to have haploid genotyps, and any apparently heterozygous genotype will be converted to `N`.

* The script can run on multiple cores (`-T` flag). Try different numbers, as using too many can slow the script down (due to the difficulty in sorting the outputs coming from the different cores).



# Citation
Lu, L., Dai, W., Pan, Y., & Yan, Z. (in preparation). A Divergence-Calibrated Framework for Robust Cross-Species Comparison of Effective Population Size Ratios. Manuscript in preparation for Methods in Ecology and Evolution.

(The manuscript and associated tool are available at: https://github.com/Ne-ratio/Ne_ratio)

# Contact
For questions, please contact Zheng Yan (yanz@lzu.edu.cn).



