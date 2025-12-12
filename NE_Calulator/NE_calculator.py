
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 15:34:57 2025

@author: Liangyu Lu
"""


import pandas as pd
import numpy as np
from itertools import product
from scipy.stats import mannwhitneyu



df = pd.read_csv("C:/Users/finni/Desktop/sliding_windows.csv")
agg = df.groupby(['species', 'chr_type']).agg(
        pi_mean=('pi_value', 'mean'),
        pi_sd=  ('pi_value', 'std'),
        div_mean=('div_value', 'mean'),
        div_sd=  ('div_value', 'std'))


nest = {}
for (sp, chrom), row in agg.iterrows():
    for x in ['pi_mean', 'pi_sd', 'div_mean', 'div_sd']:
        nest.setdefault(x, {}).setdefault(chrom, {})[sp] = row[x]

# Ne ratio and SE（Delta method）
def ratio_se(chrom, sp):
    p_q  = nest['pi_mean'][chrom][sp]
    s_q  = nest['pi_sd'][chrom][sp]
    d_q  = nest['div_mean'][chrom][sp]
    ds_q = nest['div_sd'][chrom][sp]
    p_a  = nest['pi_mean']['Autosome'][sp]
    s_a  = nest['pi_sd']['Autosome'][sp]
    d_a  = nest['div_mean']['Autosome'][sp]
    ds_a = nest['div_sd']['Autosome'][sp]

    R = (p_q / p_a) * (d_a / d_q)            # Ne ratio
    rel_var = (s_q / p_q)**2 + (s_a / p_a)**2 + (ds_a / d_a)**2 + (ds_q / d_q)**2
    SE = R * np.sqrt(rel_var)
    return R, SE

# Changes based on species name and gene name
species = ['Bonobo', 'Chimpanzee', 'Human']
chroms  = ['X', 'Y', 'MT']

res = pd.DataFrame(index=species,
                   columns=pd.MultiIndex.from_product([chroms, ['Ne_ratio', 'SE']]))
for chrom in chroms:
    for sp in species:
        r, se = ratio_se(chrom, sp)
        res.loc[sp, (chrom, 'Ne_ratio')] = r
        res.loc[sp, (chrom, 'SE')]       = se



def effective_se(x, n_block=20):
    x = np.asarray(x)
    block_size = max(1, len(x) // n_block)
    block_means = [x[i:i+block_size].mean() for i in range(0, len(x), block_size)]
    return np.std(block_means, ddof=1) / np.sqrt(len(block_means))
se_agg = (df.groupby(['species', 'chr_type'])
            .agg(pi_mean=('pi_value', 'mean'),
                 pi_se=  ('pi_value', effective_se),
                 div_mean=('div_value', 'mean'),
                 div_se=  ('div_value', effective_se)))
def ratio_se_block(chrom, sp):
    p_q  = se_agg.loc[(sp, chrom), 'pi_mean']
    se_pq= se_agg.loc[(sp, chrom), 'pi_se']
    d_q  = se_agg.loc[(sp, chrom), 'div_mean']
    se_dq= se_agg.loc[(sp, chrom), 'div_se']
    p_a  = se_agg.loc[(sp, 'Autosome'), 'pi_mean']
    se_pa= se_agg.loc[(sp, 'Autosome'), 'pi_se']
    d_a  = se_agg.loc[(sp, 'Autosome'), 'div_mean']
    se_da= se_agg.loc[(sp, 'Autosome'), 'div_se']

    R = (p_q / p_a) * (d_a / d_q)
    rel_var = (se_pq/p_q)**2 + (se_pa/p_a)**2 + (se_da/d_a)**2 + (se_dq/d_q)**2
    return R, R * np.sqrt(rel_var)
block_res = pd.DataFrame(index=species,
                         columns=pd.MultiIndex.from_product([chroms, ['Ne_ratio', 'SE']]))
for chrom in chroms:
    for sp in species:
        r, se = ratio_se_block(chrom, sp)
        block_res.loc[sp, (chrom, 'Ne_ratio')] = r
        block_res.loc[sp, (chrom, 'SE')]       = se

print(block_res.round(4))

def calc_window_ne_ratio(df_sub):
    auto = df_sub[df_sub['chr_type'] == 'Autosome']
    pi_a  = auto['pi_value'].mean()
    div_a = auto['div_value'].mean()
    df_sub = df_sub.assign(
        Ne_ratio=lambda x: (x['pi_value'] / pi_a) * (div_a / x['div_value'])
    )
    return df_sub


focus_chr = ['X', 'Y', 'MT']
df_test = df[df['chr_type'].isin(focus_chr + ['Autosome'])].copy()


window_ne_list = []
for sp in species:
    sp_df = calc_window_ne_ratio(df_test[df_test['species'] == sp])
    for chrom in focus_chr:
        tmp = sp_df[sp_df['chr_type'] == chrom]['Ne_ratio'].dropna()
        if len(tmp) == 0:
            continue
        window_ne_list.append({'species': sp, 'chr_type': chrom, 'Ne_ratio': tmp.values})

window_ne_df = pd.DataFrame(window_ne_list)

#Wilcoxon test
from scipy.stats import mannwhitneyu
import itertools

def calc_window_ne_ratio(df_sub):
    auto = df_sub[df_sub['chr_type'] == 'Autosome']
    pi_a, div_a = auto['pi_value'].mean(), auto['div_value'].mean()
    df_sub = df_sub.assign(
        Ne_ratio=lambda x: (x['pi_value'] / pi_a) * (div_a / x['div_value'])
    )
    return df_sub

focus_chr = ['X', 'Y', 'MT']
df_test = df[df['chr_type'].isin(focus_chr + ['Autosome'])].copy()

window_ne_list = []
for sp in species:
    sp_df = calc_window_ne_ratio(df_test[df_test['species'] == sp])
    for chrom in focus_chr:
        vals = sp_df.loc[sp_df['chr_type'] == chrom, 'Ne_ratio'].dropna().values
        if vals.size:
            window_ne_list.append({'species': sp, 'chr_type': chrom, 'Ne_ratio': vals})

window_ne_df = pd.DataFrame(window_ne_list)

test_res = []
for chrom in focus_chr:
    sub = window_ne_df[window_ne_df['chr_type'] == chrom]
    for sp1, sp2 in itertools.combinations(sub['species'].unique(), 2):
        x = sub[sub['species'] == sp1]['Ne_ratio'].iloc[0]
        y = sub[sub['species'] == sp2]['Ne_ratio'].iloc[0]
        stat, p = mannwhitneyu(x, y, alternative='two-sided')
        test_res.append({
            'chr_type': chrom,
            'species1': sp1,
            'species2': sp2,
            'n1': len(x),
            'n2': len(y),
            'statistic': stat,
            'p_value': p
        })

test_df = pd.DataFrame(test_res)
print('\nWilcoxon rank-sum test for window-level Ne ratio (no BH correction)')
print(test_df.round(6))

out_csv = r'C:/Users/finni/Desktop/Ne_ratio_wilcoxon_bonobo.csv'
test_df.to_csv(out_csv, index=False, float_format='%.6f')
print(f'\n {out_csv}')


block_res.to_csv(r'C:/Users/finni/Desktop/Ne_ratio_with_SE_bonobo.csv', float_format='%.4f')
long_for_r = (
    block_res
    .stack(level=0)         
    .reset_index()
    .rename(columns={'level_0': 'species',
                      'level_1': 'chr_type',
                      'Ne_ratio': 'ne_ratio',
                      'SE': 'se'})
    .assign(chr_type=lambda x: pd.Categorical(x['chr_type'],
                                              categories=['X','Y','MT'],
                                              ordered=True))
    .sort_values(['species', 'chr_type'])
    .reset_index(drop=True)
)

long_for_r.to_csv(r'C:/Users/finni/Desktop/Ne_ratio_long_for_R_bonobo.csv',
                  index=False,
                  float_format='%.6f')