# qaed C++ core

C++ port of the hot path of qaed (`hessq → iqrq / aedq → eigvec`), a faithful
line-by-line translation of the MATLAB reference. The MATLAB front-end is
unchanged; use the `*_cpp.m` drop-in wrappers in the repo root to dispatch to
this core.

## Files

- `quat.hpp` — quaternion scalar; closed-form `schur1` (standardization
  `u' q u = w + |v| i`) and the scalar Sylvester solver `sylvesterc`.
- `qmat.hpp` — column-major quaternion matrix; quaternion GEMM via 4 complex
  `zgemm` on the split `A = A1 + A2 j`; Householder kernels
  (qtfm-compatible `householder_vector` with `v = e1`); `schur2` (2×2 ordered
  quaternion Schur via `zgeev` on the 4×4 complex adjoint + quaternion
  Householder), `swapq`, `ordschurq`, `shift2`.
- `qaed.hpp` — `hessq`, `iqrq`, `aedq` (+AED window step), `sylvesterc_tri`,
  `eigvec`, `eigq`. Same control flow, deflation criteria and shift selection
  as the .m files.
- `bench.cpp` — standalone benchmark / self-test driver (no MATLAB needed).
- `qaed_mex.cpp` — MEX gateway used by `hessq_cpp / iqrq_cpp / aedq_cpp /
  eigq_cpp`.

Not yet ported: the skew/tridiagonal specializations (`skew_iqrq`,
`skew_aedq`). `eigq_cpp` handles those inputs through the general path
(correct, just not specialized).

## Build & run (standalone)

```sh
module load GCCcore/13.3.0
make            # links FlexiBLAS (-l:libflexiblas.so.3)
./qaed_bench --selftest
./qaed_bench --alg aed --n 512            # random matrix benchmark
./qaed_bench --alg aed --load H.qmat --save-out out   # cross-check with MATLAB
```

`.qmat` is a raw binary format (int64 m, n; then w,x,y,z column-major doubles);
`scripts/save_qmat.m` / `scripts/load_qmat.m` convert to/from qtfm quaternions.

## Build (MEX) and MATLAB comparison

```sh
sbatch scripts/build_mex_and_check.sh   # builds qaed_mex, runs test_cpp_vs_matlab.m
sbatch scripts/cpp_bench.sh             # standalone benchmark on a compute node
```

Inside MATLAB the build command is:

```matlab
cd cpp
mex -R2018a -DQAED_ILP64 CXXFLAGS='$CXXFLAGS -O3 -std=c++17' qaed_mex.cpp -lmwlapack -lmwblas
```

`-DQAED_ILP64` is required with MATLAB's `libmwlapack/libmwblas` (64-bit BLAS
integers). OpenMP is intentionally disabled in the MEX build (MATLAB ships
libiomp5; mixing it with libgomp is fragile) — heavy GEMMs still run on
MATLAB's multithreaded MKL.

MATLAB bundles an old `libstdc++` that lacks the `GLIBCXX` version required by
mex files built with the system g++; launch MATLAB with
`export LD_PRELOAD=/usr/lib64/libstdc++.so.6`
(as `scripts/build_mex_and_check.sh` does) or the mex file will fail to load.

## Notes / conventions

- Standardized eigenvalues use the upper-half-plane representative
  `w + |v| i` (`Im >= 0`), matching what the deflation tests and
  `sylvesterc` assume.
- `schur2` orders the two standardized eigenvalues by decreasing modulus,
  matching MATLAB's `ordschur` + `groupByValue` behavior.
- All algorithm code uses 1-based inclusive index ranges to mirror the
  MATLAB sources; compare side by side when auditing.
