\\ quat.gp -- quaternion-side computation of R(E) for random supersingular curves.
\\ Loaded by density.py; can also be used directly in gp:  read("quat.gp")
\\
\\ Algebra: B = (-1,-p), p == 3 (mod 4), i^2 = -1, j^2 = -p, k = ij.
\\ A quaternion is a column vector [a,b,c,d]~ = a + b i + c j + d k, and a
\\ lattice (order or ideal) is a 4x4 matrix whose columns form a Z-basis.
\\
\\ Deuring correspondence: for squarefree d with 4d < p, the isogenies
\\ E -> E^(p) of degree d (up to automorphisms) correspond to pairs +-alpha
\\ in O = End(E) with trd(alpha) = 0 and nrd(alpha) = p d.  So r_E(d) is the
\\ number of vectors (up to sign) of value d of the ternary form Q = nrd/p on
\\ the lattice {alpha in O : trd(alpha) = 0, p | nrd(alpha)}.

\\ ---------------------------------------------------------------- arithmetic
qinit(pp) = { P = pp; }
qmul(x, y) =
{
  [x[1]*y[1] - x[2]*y[2] - P*x[3]*y[3] - P*x[4]*y[4],
   x[1]*y[2] + x[2]*y[1] + P*(x[3]*y[4] - x[4]*y[3]),
   x[1]*y[3] + x[3]*y[1] - x[2]*y[4] + x[4]*y[2],
   x[1]*y[4] + x[4]*y[1] + x[2]*y[3] - x[3]*y[2]]~;
}
qconj(x) = [x[1], -x[2], -x[3], -x[4]]~;
qnrd(x)  = x[1]^2 + x[2]^2 + P*x[3]^2 + P*x[4]^2;
qtrd(x)  = 2*x[1];

\\ Z-span of the columns of a rational matrix, as a basis.
lat(M) = { my(d = denominator(M)); mathnf(d*M)/d; }

\\ ------------------------------------------------------- orders and the walk
\\ Maximal order Z<1, i, (i+j)/2, (1+k)/2>  (End of the curve j = 1728).
order0() = [1, 0, 0, 1/2; 0, 1, 1/2, 0; 0, 0, 1/2, 0; 0, 0, 0, 1/2];
ordgram(Om) = matrix(4, 4, m, n, qtrd(qmul(Om[,m], qconj(Om[,n]))));
ordreduce(Om) = Om * qflllgram(ordgram(Om));

\\ One step of the random 2-isogeny walk: O -> O_R(O beta + 2 O) for a
\\ uniformly random beta in O/2O with 2 | nrd(beta), beta != 0 (each of the
\\ three left ideals of norm 2 is hit by exactly three such beta).
\\ For a left ideal I of norm N, O_R(I) = conj(I) I / N.
step2(Om) =
{
  my(c, beta, I, G);
  until (c != [0,0,0,0]~ && qnrd(beta) % 2 == 0,
    c = vectorv(4, m, random(2)); beta = Om*c);
  I = lat(matconcat([matrix(4, 4, r, m, qmul(Om[,m], beta)[r]), 2*Om]));
  G = matrix(4, 16, r, t, my(m = (t-1)\4 + 1, n = (t-1)%4 + 1);
                          qmul(qconj(I[,m]), I[,n])[r]);
  ordreduce(lat(G) / 2);
}
walk(Om, t) = { for (s = 1, t, Om = step2(Om)); Om; }

\\ Sanity check: O is a ring containing 1 with the discriminant of order0().
ordcheck(Om, D0) =
{
  my(Oi = matsolve(Om, matid(4)), z);
  if (abs(matdet(ordgram(Om))) != D0, error("wrong discriminant"));
  for (m = 1, 4, for (n = 1, 4, z = Oi * qmul(Om[,m], Om[,n]);
    for (r = 1, 4, if (denominator(z[r]) != 1, error("not closed under *")))));
  z = Oi * [1,0,0,0]~;
  for (r = 1, 4, if (denominator(z[r]) != 1, error("1 not in order")));
  1;
}

\\ ----------------------------------------------------------- ternary form
\\ Integral Gram matrix G (= 2 x bilinear form) of Q = nrd/p on
\\ {alpha in O : trd(alpha) = 0, p | nrd(alpha)}, so that Q(v) = v~*G*v/2.
qform(Om) =
{
  my(K, W, G, Gi, Kp, V, GQ, U);
  K = matkerint(Mat(Om[1,]));                  \\ trace-zero combinations
  W = Om * K;
  G = matrix(3, 3, m, n, (qnrd(W[,m] + W[,n]) - qnrd(W[,m]) - qnrd(W[,n])) / 2);
  Gi = denominator(G) * G;
  if (matrank(Gi * Mod(1, P)) != 1, error("trace-zero form mod p not of rank 1"));
  Kp = lift(matker(Gi * Mod(1, P)));          \\ where p | nrd
  V = lat(matconcat([Kp, P*matid(3)]));
  GQ = 2 * (V~ * G * V / P);
  U = qflllgram(GQ);
  U~ * GQ * U;
}

\\ All values Q(v) <= X, one for each pair +-v of nonzero vectors, for a
\\ positive definite integral G.  Pure integer arithmetic: with
\\   g11 * (v~Gv) = (g11 x1 + g12 x2 + g13 x3)^2 + A x2^2 + 2 Bc x2 x3 + C x3^2
\\ every coordinate range comes from sqrtint, so no vector can be missed.
qvalues(G, X) =
{
  my(N = 2*X, g11 = G[1,1], g12 = G[1,2], g13 = G[1,3], A, Bc, C, res = List(),
     D2, s, lo2, hi2, q23, rem, r, t, lo1, hi1);
  A = g11*G[2,2] - g12^2; Bc = g11*G[2,3] - g12*g13; C = g11*G[3,3] - g13^2;
  for (x3 = 0, sqrtint((N * A) \ matdet(G)),
    D2 = (Bc*x3)^2 - A*(C*x3^2 - g11*N);
    if (D2 < 0, next);
    s = sqrtint(D2);
    lo2 = ceil((-Bc*x3 - s - 1) / A); hi2 = floor((-Bc*x3 + s + 1) / A);
    if (x3 == 0, lo2 = max(lo2, 0));
    for (x2 = lo2, hi2,
      q23 = A*x2^2 + 2*Bc*x2*x3 + C*x3^2;
      rem = g11*N - q23;
      if (rem < 0, next);
      r = sqrtint(rem); t = g12*x2 + g13*x3;
      lo1 = ceil((-t - r) / g11); hi1 = floor((-t + r) / g11);
      if (x3 == 0 && x2 == 0, lo1 = max(lo1, 1));
      for (x1 = lo1, hi1, listput(res, ((g11*x1 + t)^2 + q23) / (2*g11)))));
  Vec(res);
}

\\ -------------------------------------------------------- four-prime family
\\ A family is F = [a, SMALL, BIG] with integer edges a = [a1,...,a5]:
\\ Q_i = primes in (a_i, a_{i+1}],  D = {q1 q2 q3 q4 : q_i in Q_i}.
\\ SMALL = product of the primes <= min(a1, 10^5) (cheap filter),
\\ BIG   = product of all primes in (a1, a5]  (d in D implies d | BIG).
famprep(a) = [a, vecprod(primes([2, min(a[1], 10^5)])), vecprod(primes([a[1] + 1, a[5]]))];

\\ Exact membership test (used only on the few candidates that divide BIG).
inFam(d, F) =
{
  my(a = F[1], f, q);
  if (d <= a[1]^4 || d > a[5]^4 || gcd(d, F[2]) != 1, return(0));
  f = factor(d);
  if (#f[,1] != 4 || vecmax(f[,2]) > 1, return(0));
  q = f[,1];
  for (i = 1, 4, if (q[i] <= a[i] || q[i] > a[i+1], return(0)));
  1;
}

\\ Remainder tree: [B mod S[1], ..., B mod S[n]] in quasi-linear time.
remtree(B, S) =
{
  my(T = List([S]), L, m, r);
  if (#S == 0, return([]));
  while (#T[#T] > 1,
    L = T[#T]; m = (#L + 1) \ 2;
    listput(T, vector(m, k, if (2*k <= #L, L[2*k-1] * L[2*k], L[2*k-1]))));
  r = [B % T[#T][1]];
  forstep (l = #T - 1, 1, -1, L = T[l]; r = vector(#L, k, r[(k+1)\2] % L[k]));
  r;
}

\\ Number of entries of Qv that lie in D.
famcount(Qv, F) =
{
  my(a = F[1], S, r, c = 0);
  S = select(d -> d > a[1]^4 && d <= a[5]^4 && gcd(d, F[2]) == 1, Qv);
  r = remtree(F[3], S);
  for (k = 1, #S, if (r[k] == 0 && inFam(S[k], F), c++));
  c;
}

\\ R(E) for the order O and each family in Fl (X = largest degree over Fl).
Rfam(Om, Fl, X) = { my(Qv = qvalues(qform(Om), X)); vector(#Fl, f, famcount(Qv, Fl[f])); }

\\ R(E) for an explicit sorted list of degrees (used by the validation).
Rlist(Om, Dv) = { my(Qv = qvalues(qform(Om), Dv[#Dv])); sum(n = 1, #Qv, setsearch(Dv, Qv[n]) > 0); }

\\ ----------------------------------------------------------------- moments
\\ sum_E r_E(d) exactly, by class numbers [OWX26, Prop. 3.6] (small p only).
Aexact(p, d) =
{
  if (d % 4 == 1, (qfbclassno(-4*p*d) + qfbclassno(-p*d)) / 2,
                  qfbclassno(-4*p*d) / 2);
}
\\ sum of sqrt(q) over primes q in (a, b]
sqrtsum(a, b) = { my(s = 0.); forprime(q = a + 1, b, s += sqrt(q)); s; }
\\ T(m) = sum_{t^2 < 4m} tau(4m - t^2)  (exact; small m only)
Texact(m) = { my(s = numdiv(4*m)); for (t = 1, sqrtint(4*m - 1), s += 2*numdiv(4*m - t^2)); s; }
