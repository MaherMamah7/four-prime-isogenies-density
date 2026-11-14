#!/usr/bin/env python3
"""
Success probability of the four-prime search, measured on the quaternion side.

For a prime p = 3 (mod 4) and lambda > 0, put B = lambda * p^(1/12), split
(B/2, B] into four quarters I_1..I_4, and let
    D = { q1 q2 q3 q4 : q_i prime, q_i in I_i }.
For a supersingular curve E, R(E) = sum_{d in D} r_E(d) is the number of
isogenies E -> E^(p) with degree in D, and the search succeeds on E iff
R(E) > 0.  This program estimates Pr_E[R(E) > 0] and evaluates the lower bound
of Proposition 1 with explicit constants.

Commands
  validate   check the counting against exact class-number sums (small p)
  run        sample random curves and count those with R(E) > 0
  merge      pool several run files (same p, different seeds)
  bound      the lower bound of Proposition 1 with explicit constants

Examples
  python3 density.py validate
  python3 density.py run --bits 60 --lambdas 8 10 --samples 375 --seed 1
  python3 density.py merge runs/*.csv --out results.csv
  python3 density.py bound --bits 60 80 100 120 --lambdas 8 10 --out lower_bound.csv

Requires PARI/GP through cypari2 (Linux, macOS) or cypari (also Windows).
"""
import argparse
import csv
import math
import os
import random
import sys
import time

try:
    import cypari2
    pari = cypari2.Pari()
except ImportError:
    from cypari import pari
try:
    pari.allocatemem(4 * 10**9, silent=True)
except TypeError:
    pari.allocatemem(4 * 10**9)
_GP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "quat.gp").replace("\\", "/")
pari(f'read("{_GP}")')


# ------------------------------------------------------------------ setup
def prime_3mod4(bits):
    """Smallest prime p = 3 (mod 4) with p >= 2^bits."""
    q = 2 ** bits
    while not (pari.isprime(q) and q % 4 == 3):
        q += 1
    return q


class Family:
    """Q_i = primes in (a_i, a_{i+1}], with a_i = floor(B/2 + (i-1) B/8)."""

    def __init__(self, p, lam):
        self.p, self.lam = p, lam
        B = lam * p ** (1 / 12)
        self.a = [math.floor(B / 2 + i * B / 8) for i in range(5)]
        self.Q = [int(pari(f'primepi({self.a[i + 1]}) - primepi({self.a[i]})')) for i in range(4)]
        self.numD = math.prod(self.Q)
        if self.numD == 0:
            raise ValueError(f"lambda={lam}: an interval contains no prime")
        # largest degree: the product of the largest prime in each interval
        self.X = math.prod(int(pari(f'precprime({self.a[i + 1]})')) for i in range(4))
        if 4 * self.X >= p:
            raise ValueError("need 4d < p for every d in D")

    def gp(self):
        return f'famprep({self.a})'


def walk_length(p):
    """Each sample is a fresh walk from the order of j = 1728.  The simple
    random walk on the 3-regular 2-isogeny graph mixes after about 3 log2 p
    steps (cutoff, Lubetzky-Peres); we take 3 log2 p + 20."""
    return 3 * p.bit_length() + 20


class Sampler:
    def __init__(self, p, seed, t):
        self.t = t
        pari(f'setrand({seed}); qinit({p}); O0 = order0(); D0 = abs(matdet(ordgram(O0)));')

    def next(self):
        pari(f'Om = walk(O0, {self.t});')

    def check(self):
        return bool(pari('ordcheck(Om, D0)'))


# --------------------------------------------------------------- commands
def cmd_validate(a):
    """Mean of R over random curves vs the exact first moment, for all
    squarefree d in [dmin, dmax] at a small prime.  Under the walk's
    stationary distribution (weight 12/((p-1) #Aut(E)/2)) the mean of R
    counted up to sign equals 12/(p-1) * sum_d sum_E r_E(d)."""
    p = a.p
    assert p % 4 == 3 and pari.isprime(p)
    D = [d for d in range(a.dmin, min(a.dmax, (p - 1) // 4) + 1) if pari.issquarefree(d)]
    exact = 12 * sum(int(pari(f'Aexact({p}, {d})')) for d in D) / (p - 1)
    pari(f'DV = {D};')
    s = Sampler(p, a.seed, walk_length(p))
    tot = tot2 = 0
    t0 = time.time()
    for i in range(a.samples):
        s.next()
        if i % 100 == 0:
            assert s.check()
        r = int(pari('Rlist(Om, DV)'))
        tot += r; tot2 += r * r
    n = a.samples
    mean = tot / n
    se = math.sqrt(max(tot2 / n - mean ** 2, 0) / (n - 1))
    print(f"p = {p}, {len(D)} squarefree degrees in [{D[0]}, {D[-1]}], {n} samples ({time.time() - t0:.1f}s)")
    print(f"  sampled mean of R           = {mean:.4f} +- {se:.4f}")
    print(f"  exact first moment          = {exact:.4f}")
    print(f"  difference                  = {(mean - exact) / se:+.2f} standard errors")


def cmd_run(a):
    p = prime_3mod4(a.bits)
    fams = [Family(p, lam) for lam in a.lambdas]
    X = max(f.X for f in fams)
    pari('FL = [' + ', '.join(f.gp() for f in fams) + '];')
    s = Sampler(p, a.seed, walk_length(p))
    hits = [0] * len(fams); sumR = [0] * len(fams)
    t0 = time.time()
    for i in range(a.samples):
        s.next()
        R = [int(r) for r in pari(f'Rfam(Om, FL, {X})')]
        for k, r in enumerate(R):
            hits[k] += r > 0; sumR[k] += r
        if (i + 1) % max(1, a.samples // 10) == 0:
            print(f"  {i + 1}/{a.samples} samples, {time.time() - t0:.0f}s", flush=True)
    out = a.out or f"run_{p.bit_length()}bits_seed{a.seed}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["p", "bits", "lambda", "numD", "samples", "hits", "sumR", "seed"])
        for f, h, t in zip(fams, hits, sumR):
            w.writerow([p, p.bit_length(), f.lam, f.numD, a.samples, h, t, a.seed])
    n = a.samples
    print(f"p = {p} ({p.bit_length()} bits), {n} samples")
    for f, h in zip(fams, hits):
        pr = h / n
        print(f"  lambda = {f.lam:g}: |D| = {f.numD}, Pr[R>0] = {pr:.5f} +- {math.sqrt(pr * (1 - pr) / n):.5f}")
    print(f"wrote {out}")


def cmd_merge(a):
    agg = {}
    for fn in a.files:
        for r in csv.DictReader(open(fn)):
            k = (int(r["p"]), float(r["lambda"]))
            g = agg.setdefault(k, [0, 0])
            g[0] += int(r["samples"]); g[1] += int(r["hits"])
    rows = []
    print("bits  lambda  samples   Pr[R>0]   SE")
    for (p, lam), (n, h) in sorted(agg.items()):
        pr = h / n; se = math.sqrt(pr * (1 - pr) / n)
        rows.append([p, p.bit_length(), lam, n, h, f"{pr:.6g}", f"{se:.4g}"])
        print(f"{p.bit_length():4} {lam:7g} {n:8}   {pr:.5f}   {se:.5f}")
    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["p", "bits", "lambda", "samples", "hits", "pr", "se"])
            w.writerows(rows)
        print(f"wrote {a.out}")


# Tatuzawa: L(1, chi) > 0.655 e^{-1} / log|D| for all but at most one primitive
# real character; with |D| <= 4pd < p^2 and the factor 1/2 for imprimitive
# chi_{-4pd} absorbed, every A(d) >= TAT * sqrt(pd) / log p.
TAT = 0.655 * math.exp(-2) / (2 * math.pi)


def lower_bound(p, F, mc, rng):
    """Proposition 1 with explicit constants, per curve:
         S0 = 12/(p-1) * TAT sqrt(p)/log p * (sum_D sqrt(d) - max_D sqrt(d))
         K0 = 12/(p-1) * 288 * sum_{(d,e) in D x D} T(de)       [OWX26, Thm 5.5]
         Pr[R>0] >= S0^2 / (S0 + K0)                              (Cauchy-Schwarz)
    The double sum is a fixed integer, estimated by uniform sampling of
    (d, e, t) with |t| < 2 sqrt(de); returns (S0, K0, relative SE of K0, bound)."""
    sq = math.prod(float(pari(f'sqrtsum({F.a[i]}, {F.a[i + 1]})')) for i in range(4))
    S0 = 12 * TAT * math.sqrt(p) * (sq - math.sqrt(F.X)) / math.log(p) / (p - 1)
    primes = [[int(q) for q in pari(f'primes([{F.a[i] + 1}, {F.a[i + 1]}])')] for i in range(4)]
    vals = []
    for _ in range(mc):
        m = math.prod(rng.choice(Q) for Q in primes) * math.prod(rng.choice(Q) for Q in primes)
        tm = math.isqrt(4 * m - 1)
        t = rng.randint(-tm, tm)
        vals.append((2 * tm + 1) * int(pari(f'numdiv({4 * m - t * t})')))
    mean = sum(vals) / mc
    rel = math.sqrt(sum((v - mean) ** 2 for v in vals) / (mc - 1) / mc) / mean
    K0 = 12 * 288 * F.numD ** 2 * mean / (p - 1)
    return S0, K0, rel, S0 ** 2 / (S0 + K0)


def cmd_bound(a):
    rng = random.Random(a.seed)
    rows = []
    print("bits  lambda        |D|   S0          K0                lower bound on Pr[R>0]")
    for bits in a.bits:
        p = prime_3mod4(bits)
        for lam in a.lambdas:
            F = Family(p, lam)
            S0, K0, rel, L = lower_bound(p, F, a.mc, rng)
            rows.append([p, p.bit_length(), lam, F.numD, f"{S0:.6g}", f"{K0:.6g}", f"{rel:.3g}", f"{L:.6g}"])
            print(f"{p.bit_length():4} {lam:7g} {F.numD:10}   {S0:.3e}   {K0:.3e} (+-{rel:.1%})   {L:.3e}", flush=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["p", "bits", "lambda", "numD", "S0", "K0", "K0_rel_se", "lower_bound"])
        w.writerows(rows)
    print(f"wrote {a.out}")


def cmd_selftest(a):
    """Checks the sampled K0 against the exact double sum on a tiny family."""
    Q = [[5, 7], [11, 13], [17, 19], [23, 29]]
    D = sorted(w * x * y * z for w in Q[0] for x in Q[1] for y in Q[2] for z in Q[3])
    exact = sum(int(pari(f'Texact({d * e})')) for d in D for e in D)
    rng = random.Random(1)
    vals = []
    for _ in range(a.mc):
        m = math.prod(rng.choice(q) for q in Q) * math.prod(rng.choice(q) for q in Q)
        tm = math.isqrt(4 * m - 1); t = rng.randint(-tm, tm)
        vals.append((2 * tm + 1) * int(pari(f'numdiv({4 * m - t * t})')))
    mean = sum(vals) / a.mc
    se = math.sqrt(sum((v - mean) ** 2 for v in vals) / (a.mc - 1) / a.mc)
    est = len(D) ** 2 * mean
    print(f"sum T(de) over a 16-element family: exact {exact}, sampled {est:.6g} "
          f"+- {len(D) ** 2 * se:.3g} ({(est - exact) / (len(D) ** 2 * se):+.2f} SE)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="check the counting against class numbers")
    v.add_argument("--p", type=int, default=10007)
    v.add_argument("--dmin", type=int, default=101)
    v.add_argument("--dmax", type=int, default=1999)
    v.add_argument("--samples", type=int, default=1000)
    v.add_argument("--seed", type=int, default=1)
    r = sub.add_parser("run", help="sample curves and count successes")
    r.add_argument("--bits", type=int, required=True, help="p is the first prime = 3 mod 4 above 2^bits")
    r.add_argument("--lambdas", type=float, nargs="+", default=[8, 10])
    r.add_argument("--samples", type=int, default=375)
    r.add_argument("--seed", type=int, default=1, help="use a different seed for every job")
    r.add_argument("--out", default="")
    m = sub.add_parser("merge", help="pool run files")
    m.add_argument("files", nargs="+")
    m.add_argument("--out", default="")
    b = sub.add_parser("bound", help="lower bound of Proposition 1, explicit constants")
    b.add_argument("--bits", type=int, nargs="+", default=[60, 70, 80, 90, 100, 110, 120])
    b.add_argument("--lambdas", type=float, nargs="+", default=[8, 10])
    b.add_argument("--mc", type=int, default=10000, help="samples for the divisor sum in K0")
    b.add_argument("--seed", type=int, default=1)
    b.add_argument("--out", default="lower_bound.csv")
    t = sub.add_parser("selftest", help="check the sampled divisor sum against the exact one")
    t.add_argument("--mc", type=int, default=20000)
    a = ap.parse_args()
    {"validate": cmd_validate, "run": cmd_run, "merge": cmd_merge,
     "bound": cmd_bound, "selftest": cmd_selftest}[a.cmd](a)


if __name__ == "__main__":
    main()
