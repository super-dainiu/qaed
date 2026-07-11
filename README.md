# qaed

Quaternion eigenvalue decomposition via the implicit QR algorithm with
aggressive early deflation (AED), in MATLAB and quaternion-native C++.

- `matlab/` — the algorithms: `hessq` (Hessenberg reduction), `iqrq`
  (implicit QR), `aedq` (QR + AED), `eigvec` (eigenvectors), `eigq`
  (driver), `ordschurq`/`swapq`/`sylvesterc*` (eigenvalue reordering), and
  the skew-Hermitian tridiagonal specializations `skew_iqrq`/`skew_aedq`.
  Each function transparently dispatches to the C++ core when the mex file
  is built; `qaed_accel(false)` forces the pure-MATLAB reference path.
- `cpp/` — the C++ core (all computations in quaternion arithmetic), a
  standalone benchmark `qaed_bench`, and the MEX gateway. See `cpp/README.md`.
- `tests/` — the experiment drivers (`test_time_*`, `test_relerr_*`,
  `test_cpp_vs_matlab`, `find_best_alpha`).
- `scripts/` — slurm scripts for building and benchmarking, and plotting
  helpers.

## Setup

QTFM (Quaternion Toolbox for MATLAB) is not vendored; install it once with

```sh
scripts/setup_qtfm.sh     # downloads QTFM 3.4 into qtfm/ and applies patches/
```

Build the C++ core:

```sh
cd cpp && make            # standalone benchmark (needs BLAS/LAPACK)
sbatch scripts/build_mex_and_check.sh   # MEX + MATLAB cross-check (slurm)
```
