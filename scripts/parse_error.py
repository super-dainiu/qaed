import re
import pandas as pd

type = 'hess'
with open(f'{type}rand_err.txt') as f:
    logger_output = f.read()

data = []
# using aed, n = 64, uni_err = 9.18e-15, schur_err = 6.42e-15, eigvec_err = 3.60e-13
# parse AED
aed_pattern = re.compile(r'using aed, n = (\d+), uni_err = (\d+\.\d+e-\d+), schur_err = (\d+\.\d+e-\d+), eigvec_err = (\d+\.\d+e(\+|-)\d+)', re.DOTALL)
for aed_match in aed_pattern.finditer(logger_output):
    data.append({
        'deflation strategy': 'AED',
        'matrix size': int(aed_match.group(1)),
        'unitary error': float(aed_match.group(2)),
        'schur error': float(aed_match.group(3)),
        'eigvec error': float(aed_match.group(4)),
    })

# using iqr, n = 64, uni_err = 9.18e-15, schur_err = 6.42e-15, eigvec_err = 3.60e-13

# parse IQR
iqr_pattern = re.compile(r'using iqr, n = (\d+), uni_err = (\d+\.\d+e-\d+), schur_err = (\d+\.\d+e-\d+), eigvec_err = (\d+\.\d+e(\+|-)\d+)', re.DOTALL)
for iqr_match in iqr_pattern.finditer(logger_output):
    data.append({
        'deflation strategy': 'w/o AED',
        'matrix size': int(iqr_match.group(1)),
        'unitary error': float(iqr_match.group(2)),
        'schur error': float(iqr_match.group(3)),
        'eigvec error': float(iqr_match.group(4)),
    })

# sort deflation strategy and matrix size
df = pd.DataFrame(data).sort_values(['matrix size', 'deflation strategy']).reset_index(drop=True)
df.to_latex(f'{type}rand_err.tex', index=False, float_format='%.3e', na_rep='--')
df.to_csv(f'{type}rand_err.csv', index=False, float_format='%.3e', na_rep='--')