# Introduction #

I need to use some cryptic name for code variables and functions for the sake of duplicate naming.
In general, when the variable have indexes, its corresponding variable in the code will be a matrix of [ex](http://www.ginac.de/reference/ex_8h-source.html) with the name of prefix **mat`_`**. Example would help:
> F<sub>nmlk</sub>  <=>  `mat_F[n][m][l][k]`, and often have a function to calculate its value `ex fn_F(int n, int m, int n, int k)`



# Tables of Variables and equivalent #
| **Analytics Variable/Function** | **Code Variable** | **Calculate Function** | **Note** | **Header/Implement File** |
|:--------------------------------|:------------------|:-----------------------|:---------|:--------------------------|
| m<sup>2</sup><sub>k</sub><br>q<sub>k</sub> <table><thead><th> <code>mat_msquare[k]</code><br><code>mat_q[k]</code> </th><th>  </th><th>   Internal mass in square<br>   Momenta in parallel and orthogonal space </th><th> terms.h </th></thead><tbody>
<tr><td> a<sub>lk</sub> to d<sub>lk</sub></td><td> <code>mat_a[l][k]</code><br><code>mat_d[l][k]</code> </td><td> <code>fn_a(int l, int k)</code><br><code>fn_d(int l, int k)</code> </td><td> Level 1 integral terms </td><td> lev1.h </td></tr>
<tr><td> AC<sub>lk</sub><br>alpha<sub>lk</sub><br>A<sub>lk</sub><br>B<sub>lk</sub><br>C<sub>lk</sub><br>D<sub>lk</sub><br>f<sub>lk</sub><br>f<sup>-</sup><sub>lk</sub></td><td> <code>mat_AC[l][k]</code><br><code>mat_alpha[l][k]</code><br><code>mat_A[l][k]</code><br><code>mat_B[l][k]</code><br><code>mat_X[l][k]</code><br><code>mat_D[l][k]</code><br><code>mat_f[l][k]</code><br><code>mat_fminus[l][k]</code> </td><td> <code>fn_alpha(int l, int k)</code> </td><td> Level 2 integral terms </td><td> lev2.h </td></tr>
<tr><td> beta<sub>mlk</sub><br>phi<sub>mlk</sub><br>g<sub>mlk</sub><br>g<sup>-</sup><sub>mlk</sub><br>Q<sub>mlk</sub><br>P<sub>mlk</sub><br>E<sub>mlk</sub><br>F<sub>nmlk</sub> </td><td> <code>mat_beta[m][l][k]</code><br><code>mat_phi[m][l][k]</code><br><code>mat_g[m][l][k]</code><br><code>mat_gminus[m][l][k]</code><br><code>mat_Q[m][l][k]</code><br><code>mat_P[m][l][k]</code><br><code>mat_E[m][l][k]</code><br><code>mat_F[n][m][l][k]</code></td><td> fn_beta(int m, int l, int k)<br>fn_phi(int m, int l, int k)<br>fn_g(int m, int l, int k)<br>fn_gminus(int m, int l, int k)<br>fn_Q(int m, int l, int k)<br>fn_P(int m, int l, int k)<br>fn_E(int m, int l, int k)<br>fn_F(int n, int m, int l, int k) </td><td> Level 3 integral terms </td><td> lev3.h</td></tr>