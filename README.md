# Success probability of the four-prime isogeny search

This package reproduces the experiment and the figure supporting Proposition 1
(four-prime success density). It measures how often a random supersingular
curve E has an isogeny E → E^(p) with degree in the four-prime family, and
compares that frequency with the lower bound proved in the paper.

## What is computed

Let p ≡ 3 (mod 4) be prime. Fix λ > 0, put B = λ·p^{1/12}, and split (B/2, B]
into four quarters I_1, …, I_4. The family of degrees is

    D = { q1 q2 q3 q4 : q_i prime, q_i ∈ I_i }.

Let R(E) be the number of isogenies E → E^(p) with degree in D. The search in
the paper succeeds on E exactly when R(E) > 0.

The code works on the quaternion side, with no elliptic curves and no modular
polynomials:

* **Random curves.** A curve is represented by its endomorphism ring, a maximal
  order O in the quaternion algebra ramified at p and ∞. Each sample is a fresh
  random walk of length 3·log₂p + 20 in the 2-isogeny graph, starting from the
  order of j = 1728. One step replaces O by the right order of the ideal
  Oβ + 2O for a random β ∈ O with 2 | nrd(β).
* **Counting.** For squarefree d with 4d < p, the isogenies E → E^(p) of degree d
  correspond to the elements α ∈ O, up to sign, with trd(α) = 0 and
  nrd(α) = p·d. These are the vectors of value d of the ternary quadratic form
  Q = nrd/p on {α ∈ O : trd α = 0, p | nrd α}. The code enumerates all values
  of Q up to max D and counts those in D.
  - The enumeration uses only integer arithmetic. Every coordinate range comes
    from an integer square root, so no vector is missed.
  - Membership in D is tested for all values at once with a remainder tree:
    first d | ∏(primes in (B/2, B]), then the few survivors are factored.
* **Lower bound.** Proposition 1 with explicit constants:
  Pr[R(E) > 0] ≥ S₀²/(S₀ + K₀).
  - S₀ is the first moment bounded below by Tatuzawa's theorem, dropping the
    one possible exceptional degree.
  - K₀ is the collision count bounded above by [OWX26, Thm 5.5], i.e.
    288·Σ_{d,e∈D} Σ_{t²<4de} τ(4de − t²).
  - This double sum is a fixed integer, too large to sum directly. It is
    estimated by uniform random sampling, with relative standard error about 2%.

We use λ ∈ {8, 10}. With λ = 1, the family of the paper, successes are too rare
to measure; the proof applies verbatim to every fixed λ.

## Requirements

Python 3.9+ and PARI/GP through one of:

* `pip install cypari2` (Linux, macOS);
* `pip install cypari` (also works on native Windows).

The figure also needs `pip install matplotlib`.

## Quick check (a few seconds)

    python3 density.py validate

This samples 1000 random curves for p = 10007 and counts all isogenies E → E^(p)
of squarefree degree in [101, 1999]. It compares their average with the exact
value given by class numbers [OWX26, Prop. 3.6]. The difference should be within
about ±2 standard errors. Expected output:

    sampled mean of R           = 1982.7480 +- 0.7404
    exact first moment          = 1981.9692
    difference                  = +1.05 standard errors

To check the sampling used for K₀ against an exact sum on a small family, run
`python3 density.py selftest`.

## Reproducing the data

Each run uses p = the first prime ≡ 3 (mod 4) above 2^bits. Samples are
independent, so split them over cores with different seeds. On 8 cores
(macOS/Linux shell):

    for b in 60 80 100 120; do
      for s in 1 2 3 4 5 6 7 8; do
        python3 density.py run --bits $b --lambdas 8 10 --samples 375 --seed $s \
          --out runs/run_${b}_seed$s.csv > log_${b}_$s.txt &
      done; wait
    done
    python3 density.py merge runs/*.csv --out results.csv

Create the `runs` folder first (`mkdir runs`). Approximate wall-clock time with
8 parallel jobs on a laptop:

| bits | time per sample, one job | 3000 samples on 8 cores |
|---|---|---|
| 60 | ~2.5 s | ~15 min |
| 80 | ~4 s | ~25 min |
| 100 | ~6 s | ~40 min |
| 120 | ~4–6 s | ~25–40 min |

Nearly all of the time goes to λ = 10, about 10⁶ vectors per sample. The cost
per sample hardly depends on p: only the walk gets longer. The success
probability decays only polylogarithmically, so larger p is also feasible.

On Windows PowerShell, start the jobs with:

    $procs = 1..8 | ForEach-Object { Start-Process python -NoNewWindow -PassThru `
      -ArgumentList "density.py","run","--bits","120","--lambdas","8","10","--samples","375","--seed","$_" `
      -RedirectStandardOutput "log120_$_.txt" }
    $procs | Wait-Process

Then compute the lower bound (about a minute) and draw the figure:

    python3 density.py bound --bits 60 70 80 90 100 110 120 --lambdas 8 10 --out lower_bound.csv
    python3 make_figure.py --results results.csv --bound lower_bound.csv

## Data used in the paper (`data/`)

* `runs/` — the run files: two seeds at 61 bits (3000 samples in total), and
  pooled seeds 1–8 at 81, 101 and 121 bits (3000 samples each).
* `results.csv` — the output of `merge`, i.e. Pr[R(E) > 0] and its standard
  error for each (p, λ).
* `lower_bound.csv` — the output of `bound`.

Running `python3 make_figure.py` with no arguments redraws the figure from these
files.

| bits | λ = 8 | λ = 10 | proven lower bound |
|---|---|---|---|
| 61 | 0.205 ± 0.007 | 0.596 ± 0.009 | 2.2·10⁻¹¹ |
| 81 | 0.110 ± 0.006 | 0.314 ± 0.008 | 1.0·10⁻¹¹ |
| 101 | 0.061 ± 0.004 | 0.180 ± 0.007 | 5.5·10⁻¹² |
| 121 | 0.035 ± 0.003 | 0.128 ± 0.006 | 3.3·10⁻¹² |

## Files

| file | contents |
|---|---|
| `quat.gp` | quaternion arithmetic, the 2-isogeny walk on maximal orders, the ternary form, exact enumeration, family membership |
| `density.py` | commands `validate`, `run`, `merge`, `bound`, `selftest` |
| `make_figure.py` | the figure (log₁₀ Pr[R(E) > 0] against the bit length of p) |
| `data/` | the data behind the figure in the paper |

## Notes

* **Sampling distribution.** The walk samples curves with weight proportional
  to 1/#Aut(E), the stationary distribution of the 2-isogeny walk. This differs
  from the uniform distribution only at j = 0 and j = 1728, an O(1/p) effect.
* **Same results everywhere.** Results depend only on `--seed` and are the same
  under cypari2 and cypari.
