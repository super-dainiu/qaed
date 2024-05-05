import re
import pandas as pd


type = 'hess'
with open(f'{type}rand.txt') as f:
    logger_output = f.read()

# Using aed ...
# Total step of AED 190
# Time for H matrix multiplications: 1.895366 s
# Time for computing Q matrix: 0.804389 s
# Time for Householder computations: 0.895097 s
# Time for AED = 2.227258 s, aed deflation = 44
# Total time = 6.221653 s, total deflation = 64
# Time without constructing Q = 5.417292 s

# parse AED
data = []
aed_pattern = re.compile(r'Using aed ...\nTotal step of AED (\d+).*?Time for H matrix multiplications: (\d+\.\d+) s.*?Time for computing Q matrix: (\d+\.\d+) s.*?Time for Householder computations: (\d+\.\d+) s.*?Time for AED = (\d+\.\d+) s, aed deflation = (\d+).*?Total time = (\d+\.\d+) s, total deflation = (\d+).*?Time without constructing Q = (\d+\.\d+) s', re.DOTALL)
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
# Total steps of IQR: 200
# Time for H matrix multiplications: 1.936159 s
# Time for computing Q matrix: 0.811810 s
# Time for Householder computations: 0.918215 s
# Total time: 4.111669 s, total deflation: 64
# Time without constructing Q: 3.299899 s

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