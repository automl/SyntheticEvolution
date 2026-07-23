from pathlib import Path
import logging
from typing import Dict, List, Any, Optional, Sequence, Tuple
import sys
import json
import numpy as np

p: Dict[str, List[str]] = {
    'CU': ['GU', 'AU', 'CG'],
    'CA': ['UA', 'CG'],
    'GA': ['UA', 'GC', 'GU'],
    'CC': ['CG', 'GC'],
    'AA': ['AU', 'UA'],
    'UU': ['AU', 'UA'],
    'GG': ['GC', 'CG', 'GU', 'UG'],
    'GC': ['AU', 'CG', 'GC'],
    'CG': ['GC', 'AU'],
    'AU': ['UA', 'GC'],
    'UA': ['AU', 'CG'],
    'GU': ['GC', 'AU', 'UG'],
    'UG': ['GC', 'AU', 'GU']
}

WC_PAIRS: List[str] = ['AU', 'UA', 'GC', 'CG']

STEM_KEEP_PROB = 1 - 0.1

for base1 in "ACGU":
    for base2 in "ACGU":
        options = p.get(base1 + base2, WC_PAIRS)
        bin_arr = np.array([[pair[0] == base1, pair[1] == base2] for pair in options])
        cov = np.cov(bin_arr, rowvar=False, ddof=0)
        print(f"{base1 + base2}: var1: {cov[0][0]:.2f}, var2: {cov[1][1]:.2f}, cov: {cov[0][1]:.2f}")


# RNA_SEQ =   "AAACAGAUCACCCGCUGAGCGGGUUAUCUGUU"
# STRUCTURE = "()()()()()()()()()()()()()()()()"

# AU GC - GU