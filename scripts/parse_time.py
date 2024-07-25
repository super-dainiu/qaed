import re
import pandas as pd


type = 'skew'
with open(f'{type}rand.txt') as f:
    logger_output = f.read()

# Using aed ...
# Total steps of AED: 125
# Time for H matrix multiplications: 1.172720 s
# Time for computing Q matrix: 1.228778 s
# Time for Householder computations: 0.685299 s
# Time for AED: 0.030266 s, AED deflation: 6
# Total time: 3.739278 s, total deflation: 64
# Time without constructing Q: 2.510528 s
# Elapsed time is 3.777767 seconds.

# parse AED
data = []
aed_pattern = re.compile(r'Using aed ...\nTotal steps of AED: (\d+).*?Time for H matrix multiplications: (\d+\.\d+) s.*?Time for computing Q matrix: (\d+\.\d+) s.*?Time for Householder computations: (\d+\.\d+) s.*?Time for AED: (\d+\.\d+) s, AED deflation: (\d+).*?Total time: (\d+\.\d+) s, total deflation: (\d+).*?Time without constructing Q: (\d+\.\d+) s', re.DOTALL)
for aed_match in aed_pattern.finditer(logger_output):
    data.append({
        'deflation strategy': 'AED',
        'matrix size': int(aed_match.group(8)),
        'total QR sweeps': int(aed_match.group(1)),
        'total time': float(aed_match.group(7)),
        'time to construct Q': float(aed_match.group(3)),
        'time for AED': float(aed_match.group(5)),
    })
    

# Using iqr ...
# Total steps of IQR: 110
# Time for H matrix multiplications: 0.905635 s
# Time for computing Q matrix: 0.466697 s
# Time for Householder computations: 0.530911 s
# Total time: 2.417120 s, total deflation: 64
# Time without constructing Q: 1.950463 s
# Elapsed time is 2.418864 seconds.

# parse IQR
iqr_pattern = re.compile(r'Using iqr ...\nTotal steps of IQR: (\d+).*?Time for H matrix multiplications: (\d+\.\d+) s.*?Time for computing Q matrix: (\d+\.\d+) s.*?Time for Householder computations: (\d+\.\d+) s.*?Total time: (\d+\.\d+) s, total deflation: (\d+).*?Time without constructing Q: (\d+\.\d+) s', re.DOTALL)
for iqr_match in iqr_pattern.finditer(logger_output):
    data.append({
        'deflation strategy': 'w/o AED',
        'matrix size': int(iqr_match.group(6)),
        'total QR sweeps': int(iqr_match.group(1)),
        'total time': float(iqr_match.group(5)),
        'time to construct Q': float(iqr_match.group(3)),
        'time for AED': None,
    })

# sort deflation strategy and matrix size
df = pd.DataFrame(data).sort_values(['matrix size', 'deflation strategy']).reset_index(drop=True)
df.to_latex(f'{type}rand.tex', index=False, float_format='%.3e', na_rep='N/A')
df.to_csv(f'{type}rand.csv', index=False, float_format='%.3e', na_rep='N/A')