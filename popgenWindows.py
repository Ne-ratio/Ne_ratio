# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 09:59:39 2025

@author: Liangyu Lu
"""



#!/usr/bin/env python3
"""
Calculate nucleotide diversity (π) in non-overlapping windows from VCF files.

This is a replacement for the original popgenWindows.py from genomics_general,
optimized for calculating π across different genomic compartments.
"""

import argparse
import gzip
import math
import sys
import numpy as np
from collections import defaultdict
import pandas as pd


def parse_vcf(vcf_file, samples=None, ploidy=2):
    """
    Parse VCF file and extract genotype information.
    
    Args:
        vcf_file (str): Path to VCF file
        samples (list): List of samples to include (None = all)
        ploidy (int): Ploidy of the organism
    
    Returns:
        dict: Dictionary with chromosome -> positions -> genotypes
    """
    genotypes = defaultdict(dict)
    samples_list = []
    
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                header = line.strip().split('\t')
                all_samples = header[9:]
                if samples is None:
                    samples_list = all_samples
                else:
                    samples_list = [s for s in all_samples if s in samples]
                sample_indices = [header.index(s) for s in samples_list]
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(',')
            
            # Skip multi-allelic sites for simplicity
            if len(alts) > 1:
                continue
            
            # Parse genotypes
            gt_list = []
            for idx in sample_indices:
                gt_field = fields[idx]
                gt_parts = gt_field.split(':')[0]  # Take GT part only
                
                if '/' in gt_parts:
                    alleles = gt_parts.split('/')
                elif '|' in gt_parts:
                    alleles = gt_parts.split('|')
                else:
                    alleles = [gt_parts]
                
                # Convert to numerical representation
                num_gt = []
                for allele in alleles:
                    if allele == '.' or allele == './.':
                        num_gt.append(-1)  # Missing data
                    else:
                        try:
                            num_gt.append(int(allele))
                        except ValueError:
                            num_gt.append(-1)
                
                gt_list.append(num_gt)
            
            genotypes[chrom][pos] = {
                'ref': ref,
                'alt': alts[0] if alts[0] != '.' else '',
                'genotypes': gt_list
            }
    
    return genotypes, samples_list


def calculate_pi_window(genotypes, ploidy=2):
    """
    Calculate nucleotide diversity (π) for a window of sites.
    
    Formula: π = (n/(n-1)) * ΣΣ_{i<j} π_ij
    where π_ij is the number of pairwise differences between sequences i and j
    """
    if not genotypes:
        return 0.0, 0
    
    n_samples = len(genotypes[0]['genotypes'])
    if n_samples < 2:
        return 0.0, 0
    
    total_pairs = (n_samples * (n_samples - 1)) / 2
    total_sites = 0
    total_differences = 0
    
    for site_data in genotypes:
        gts = site_data['genotypes']
        
        # Count differences for this site
        site_differences = 0
        valid_pairs = 0
        
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                gt_i = gts[i]
                gt_j = gts[j]
                
                # Skip if any genotype is missing
                if -1 in gt_i or -1 in gt_j:
                    continue
                
                # Count differences (simple biallelic case)
                diff_count = 0
                for a in gt_i:
                    for b in gt_j:
                        if a != b:
                            diff_count += 1
                
                site_differences += diff_count
                valid_pairs += 1
        
        if valid_pairs > 0:
            total_differences += site_differences
            total_sites += 1
    
    if total_sites == 0:
        return 0.0, 0
    
    # Average differences per site
    pi = total_differences / (total_sites * total_pairs) if total_pairs > 0 else 0.0
    
    return pi, total_sites


def process_chromosome_windows(genotypes, window_size=50000, min_sites=1):
    """
    Process a chromosome in non-overlapping windows.
    """
    windows = []
    
    if not genotypes:
        return windows
    
    chrom = list(genotypes.keys())[0] if len(genotypes) == 1 else 'unknown'
    positions = sorted(genotypes.keys())
    
    if not positions:
        return windows
    
    min_pos = positions[0]
    max_pos = positions[-1]
    
    for window_start in range(min_pos, max_pos, window_size):
        window_end = window_start + window_size
        window_genotypes = []
        
        for pos in positions:
            if window_start <= pos < window_end:
                window_genotypes.append(genotypes[pos])
        
        pi, n_sites = calculate_pi_window(window_genotypes)
        
        if n_sites >= min_sites:
            windows.append({
                'chrom': chrom,
                'start': window_start,
                'end': window_end,
                'pi': pi,
                'sites': n_sites
            })
    
    return windows


def main():
    parser = argparse.ArgumentParser(description='Calculate nucleotide diversity in non-overlapping windows')
    parser.add_argument('--vcf', required=True, help='Input VCF file (can be gzipped)')
    parser.add_argument('--windType', choices=['coord', 'sites'], default='coord',
                       help='Window type: coordinate-based or site-based')
    parser.add_argument('--windSize', type=int, default=50000,
                       help='Window size in bp (for coord) or sites (for sites)')
    parser.add_argument('--stepSize', type=int,
                       help='Step size for sliding windows (default: windSize for non-overlapping)')
    parser.add_argument('--minSites', type=int, default=1,
                       help='Minimum number of sites in window to output')
    parser.add_argument('--outFile', required=True,
                       help='Output file name')
    parser.add_argument('--samples', nargs='+',
                       help='Specific samples to include (default: all)')
    parser.add_argument('--ploidy', type=int, default=2,
                       help='Ploidy of samples (default: 2)')
    parser.add_argument('--chromosomes', nargs='+',
                       help='Specific chromosomes to process (default: all)')
    
    args = parser.parse_args()
    
    if args.stepSize is None:
        args.stepSize = args.windSize
    
    print(f"Loading VCF: {args.vcf}", file=sys.stderr)
    genotypes, samples = parse_vcf(args.vcf, args.samples, args.ploidy)
    print(f"Loaded {len(samples)} samples", file=sys.stderr)
    
    all_windows = []
    
    for chrom in genotypes:
        if args.chromosomes and chrom not in args.chromosomes:
            continue
        
        print(f"Processing chromosome {chrom}...", file=sys.stderr)
        chrom_windows = process_chromosome_windows(
            {chrom: genotypes[chrom]}, 
            window_size=args.windSize,
            min_sites=args.minSites
        )
        all_windows.extend(chrom_windows)
    
    # Create DataFrame and save
    df = pd.DataFrame(all_windows)
    
    # Reorder columns
    if not df.empty:
        df = df[['chrom', 'start', 'end', 'pi', 'sites']]
        df.to_csv(args.outFile, index=False, sep='\t')
        print(f"Results saved to {args.outFile}", file=sys.stderr)
        print(f"Total windows: {len(df)}", file=sys.stderr)
        print(f"Mean π: {df['pi'].mean():.6f}", file=sys.stderr)
    else:
        print("No windows passed filters", file=sys.stderr)
        # Create empty output file with header
        empty_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'pi', 'sites'])
        empty_df.to_csv(args.outFile, index=False, sep='\t')


if __name__ == '__main__':
    main()