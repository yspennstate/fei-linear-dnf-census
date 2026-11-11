# The Fourier entropy of linear monotone DNFs on few variables

Paper, programs and data for the census of linear monotone DNFs with at most eight relevant variables.

A monotone DNF is linear when any two of its terms share at most one variable. Up to permutation of variables there are
552,251 such functions with exactly eight relevant variables (18,785,639,553 labelled ones). Over all of them

* the largest value of H/I lies in [3.88910541981418, 3.88910541981709] and is attained only by the OR of four variables
  with the 4-cycle DNF on four others;
* the largest value of min_i E_i/Inf_i lies in [3.30470132244542, 3.30470132245649] and is attained only by two disjoint
  4-cycle DNFs.

The paper also proves a load certificate for a single term, a bias bound for H/I, and shows that affine transversal
designs defeat all of these local criteria while their ratios stay small (H/I = 4.1345862211 on TD(3,7)).

## Contents

* `main.tex`, `main.pdf` - the paper.
* `census/` - the census and its checks.
  * `census.cpp` enumerates linear hypergraphs with edges of size at least two by exact canonical forms and writes
    `census<n>.families` (one family per line, terms as bitmasks) and `census<n>.json` (counts and maximisers).
  * `labeled_check.cpp` compares the pair-recursion count of labelled families with the orbit-stabilizer sum over the
    census (`census<n>.labeled_check.json`).
  * `certify.cpp` and `log_intervals.py` give outward-rounded intervals for both ratios (`census<n>.certified.json`,
    `log_intervals.json`).
  * `structured.cpp` computes the transversal designs (`TD_<w>_<Q>.json`).
  * `analyze.py`, `proof_checks.py`, `final_audit.py` recompute the maximisers with exact arithmetic and check the
    numerical claims used in the proofs (`analysis.json`, `exact_values.json`, `proof_checks.json`, `final_audit.json`).
  * `census8.families.gz` is the eight-variable list, compressed (21 MB unpacked, sha256
    c9b67ca7be231bf1f5016ff464e750a69a2a7ba1d80b210083ce5251a483593f).
* `check/` - checks written separately from the census code.
  * `indep_census_check.c` reads a census file, checks that every family is a linear antichain covering all n variables,
    recomputes canonical forms by a different method, the automorphism counts and both ratios from integer Walsh
    transforms (`out_n<n>.json`).
  * `td_check.py` recomputes Table 2 of the paper (`td_check_out.txt`).

## Reproducing

```
cd census
g++ -O2 -std=c++17 census.cpp -o census && ./census 8 census8
g++ -O2 -std=c++17 labeled_check.cpp -o labeled_check && ./labeled_check 8 census8
g++ -O2 -std=c++17 certify.cpp -o certify && ./certify 8 census8
cd ../check
gcc -O2 indep_census_check.c -o indep -lm
gunzip -k ../census/census8.families.gz && ./indep 8 ../census/census8.families
python3 td_check.py 3 7
```

The eight-variable census takes a few hours on one core; n <= 7 takes minutes. `td_check.py` needs numpy and about
1 GB of memory for TD(5,5) (25 variables).

## Related papers

* [The one-coordinate step in the Fourier entropy-influence conjecture](https://github.com/yspennstate/fourier-entropy-influence-one-coordinate)
* [Good coordinates in the Fourier entropy-influence conjecture](https://github.com/yspennstate/fei-good-coordinates)
* [A variance form of Han's bound on the Fourier entropy](https://github.com/yspennstate/fei-variance-han-bound)
* [Four equivalent forms of the Fourier entropy-influence conjecture](https://github.com/yspennstate/fei-equivalent-forms)
* [A filtration bound for the erasure entropies of Boolean functions](https://github.com/yspennstate/fei-filtration-witness)
* [Erasure entropies of sparse linear DNFs](https://github.com/yspennstate/fei-sparse-linear-dnf)
* [Hidden coordinates and selection rules in the Fourier entropy-influence conjecture](https://github.com/yspennstate/fei-hidden-coordinates)
* [A forced coordinate with nearly maximal erasure entropy](https://github.com/yspennstate/fei-forced-coordinate)
* [The erasure entropy of a forced coordinate under every bias](https://github.com/yspennstate/fei-forced-coordinate-bias)

## License

The text and data are released under CC BY 4.0 and the code under the MIT license (see [LICENSE](LICENSE)).
