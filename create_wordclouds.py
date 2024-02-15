"""
This script creates wordclouds from the results of the pubmed_trends_analysis tool
"""

from collections import defaultdict
import re

import numpy as np
import matplotlib.pyplot as plt
from wordcloud import WordCloud

# Result file created by dig_analysis.sh
RESULTS = "results/omics_trend_2024-02-09_15-21_dig_analysis.tsv"

with open(RESULTS, "r", encoding='utf-8') as f:
    words, keywords = np.genfromtxt(f, usecols=[2, 4], dtype='str', delimiter="\t").transpose()

# Keywords of interest (given in omics_keywords.txt)
targets = ['proteomics', 'metabolomics', 'lipidomics', 'glycomics']

# Function create wordclound in loop
def word_cloud(dict_):
    wc = WordCloud(max_words=90,
                   background_color="white",
                   contour_width=0,
                   colormap='Dark2').generate_from_frequencies(dict_)
    plt.axis('off')
    return wc

# Dict of words to replace in wordclouds
to_replace = {
            'ms': 'MS',
            'gc-ms': 'GC-MS',
            'lc-ms/ms': 'LC-MS/MS', 
            'lc-ms': 'LC-MS',
            'maldi-ms': 'MALDI-MS',
            'maldi': 'MALDI',
            'maldi-tof ms': 'MALDI-TOF MS',
            'maldi-tog-ms': 'MALDI-TOF MS',
            'nmr': 'NMR',
            'silac': 'SILAC',
            'itraq': 'iTraq',
            'tmt': 'TMT',
            'covid-19': 'COVID-19', 
            'sars-cov-2': 'SARS-CoV-2',
            'igg': 'IgG',
            'o-glycosylation': 'O-glycosylation',
            'n-glycosylation': 'N-glycosylation',
            'n-glycan': 'N-glycan',
            'n-glycome': 'N-glycome',
            'n-linked glycans': 'N-linked glycans',
            'n-linked glycosylation': 'N-linked glycosylation',
            'o-glycan': 'O-glycan',
            'rna-seq': 'RNA-Seq',
            'nafld': 'NAFLD',
            'hplc': 'HPLC',
            'immunoglubulin g': 'immunoglobulin G',
            'omic': 'omics',
            'multi-omic': 'multi-omics'
          }

# Count words frequencies and create wordclouds
for target in targets:
    # Get words related to the target keyword (but remove the keyword itself)
    pat = re.compile(f'{target}|{target[:-1]}')
    list_words = [pat.sub('', word.lower())
                  for word_string in words[keywords == target] for word in word_string.split(';')]
    word, count = np.unique(list_words, return_counts=True)
    res = dict(zip(word, count))\

    key_to_remove = ['']
    items_to_add = defaultdict(int)
    for w, c in res.items():

        if  w+'s' in res:
            res[w] += res[w+'s']
            key_to_remove.append(w+'s')

        # Change some words/accronyms
        if w in to_replace:
            items_to_add[to_replace[w]] = res[w]
            key_to_remove.append(w)

    # Remove lower case accronyms and words plural form
    for k in np.unique(key_to_remove):
        res.pop(k)

    # Update dict with replaced terms
    res = res | items_to_add

    if 'lc-ms/m' in res:
        res['LC-MS/MS'] += res['lc-ms/m']
        res.pop('lc-ms/m')

    wc = word_cloud(res)
    wordcloud_svg = wc.to_svg(embed_font=True)
    with open(f'{target}_wc.svg', 'w') as f:
        f.write(wordcloud_svg)
