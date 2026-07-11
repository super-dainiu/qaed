# qaed C++ core

C++ port of the hot path of qaed (`hessq → iqrq / aedq → eigvec`), a faithful
line-by-line translation of the MATLAB reference. The MATLAB front-end is
unchanged; each function in `matlab/` dispatches here transparently when the
mex file is built (`qaed_accel(false)` forces the pure-MATLAB path).

All computation is native quaternion arithmetic ("quaternion BLAS"): the
scalar type is `Quat` (4 doubles) with direct quaternion multiplication, and
matrix products use the hand-written `qgemm` kernel (column-major gaxpy with
hoisted sign-permuted vectors, OpenMP over columns). Nothing is mapped to
complex or real BLAS. The only complex computation is the 4x4 complex adjoint
eigenproblem inside `schur2` (2x2 quaternion blocks), which is exactly what
the MATLAB reference does via `schur(adjoint(A))`. A complex-split GEMM
variant exists behind `-DQAED_COMPLEX_GEMM` for benchmark comparison only.

## Files

- `quat.hpp` — quaternion scalar; closed-form `schur1` (standardization
  `u' q u = w + |v| i`) and the scalar Sylvester solver `sylvesterc`.
- `qmat.hpp` — column-major quaternion matrix; native quaternion GEMM
  (`qgemm`); Householder kernels (qtfm-compatible `householder_vector` with
  `v = e1`); `schur2` (2×2 ordered quaternion Schur via `zgeev` on the 4×4
  complex adjoint, mirroring the MATLAB reference, + quaternion Householder),
  `swapq`, `ordschurq`, `shift2`.
- `qaed.hpp` — `hessq`, `iqrq`, `aedq` (+AED window step), `sylvesterc_tri`,
  `eigvec`, `eigq`. Same control flow, deflation criteria and shift selection
  as the .m files.
- `bench.cpp` — standalone benchmark / self-test driver (no MATLAB needed).
- `qaed_mex.cpp` — MEX gateway used by the dispatch in `matlab/`.

## Build & run (standalone)

```sh
module load GCCcore/13.3.0
make            # links FlexiBLAS (-l:libflexiblas.so.3)
./qaed_bench --selftest
./qaed_bench --alg aed --n 512            # random matrix benchmark
```

Cross-checking against MATLAB goes through the MEX gateway
(`test_cpp_vs_matlab.m`), not through files.

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
