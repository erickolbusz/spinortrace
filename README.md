# spinortrace

Computes (scalar) traces of gamma matrices including 4-vector and metric tensor contractions.

## Quick Usage

The main function is `fulltrace` which takes the arguments:
- <u>terms (required):</u> an iterable (usually a list) of the indices and momenta to trace. 
	- A $\gamma^+$ is denoted by `"+"`, a $\gamma^-$ by `"-"`, and *any other* string is interpreted as a label for a perp index, unless it begins with a `?` to show that it takes all possible values.
	- A Feynman-slashed momentum will automatically be evaluated over all possible indices.
- <u>symmetries (optional):</u> a list of tuples `(index1, index2)` that enforce that the answer is symmetric in each pair of indices
- <u>Gterms (optional):</u> a list of tuples `(vector, (index1, index2))` representing gluon propagators (e.g. $G_{\beta\gamma}(l)$ ) contracted with the trace. 
- <u>guppers (optional):</u> a list of tuples `(index1, index2)` representing metric tensors with two upstairs indices contracting with the trace.
- <u>glowers (optional):</u> same as `guppers` but with two downstairs indices.

The program loops over all possibilities for the indices and computes each term, then prints the sum in LaTeX code.

## 4-Vectors

Vectors are defined as instances of the `lcvec` class. Vector components should first be defined as sympy symbols (see the definition of `p1`, etc.) then the vector is defined from its $+,-,\perp$ components as e.g. `p2v = lcvec(0,p2,p2perp)`. The following vectors are included by default:
- $p_2$ = `p2v` = $\langle 0, p_2^-, p_{2\perp} \rangle$
- $p_1$ = `p1v` = $\langle p_1^+, 0, p_{1\perp} \rangle$
- $p_1'$ = `p1pv` = $\langle p_1^{'+}, 0, p_{1\perp} \rangle$
- $l$ = `lv` = $\langle \frac{l_\perp^2}{2l^-}, l^-, l_\perp \rangle$
- $k$ = `kv` = $\langle k^+, 0, k_\perp \rangle$
- $k_1$ = `k1v` = $\langle k_1^+, 0, k_{1\perp} \rangle$
- $k_2$ = `k2v` = $\langle k_2^+, 0, k_{2\perp} \rangle$

with $l^-$ being substituted with $yp_2^-$ when printing the final expressions.

## Gluon Propagators

All gluon propagators are in the $A^- = 0$ gauge with $n^+ = 1$. The nonzero terms are:
```math
G_{--}(l) = \frac{l_\perp^2}{(l^-)^2}\ ,\ G_{-\perp a}(l) = \frac{l_{\perp a}}{l^-}\ ,\ G_{\perp a\perp b}(l) = -g_{\perp a\perp b}\ .
```
To properly include a gluon term, the associated indices in `terms` **must** start with a `?`, otherwise they will be considered to be perp indices and only $G_{\perp\perp}$ will enter.

## $g_{\mu\nu}$ Contraction

If the trace includes a (non-perp) $\gamma^\mu$ and $\gamma^\nu$ that are contracted with a $g_{\mu\nu}$, use `MULABEL` and `NULABEL` in `terms` and the contraction will automatically be calculated by anticommuting the two until they are consecutive terms.

## Quirks/Bugs

- Currently wildcard indices are only useful for gluon propagators. They are not fully supported for use with `guppers` or `glowers` (though the code is technically there already). This shouldn't matter because currently the only non-perp metric tensors being contracted is $g_{\mu\nu}$.

---

## Examples

Basic usage:

`fulltrace(['+','-'])` or `fulltrace('+-')`
```math
\textrm{tr}\,\left[\gamma^+\gamma^-\right] = 4
```

With a perp metric tensor:

`fulltrace(['+','-','a','b'],glowers=[('a','b')])` or `fulltrace('+-ab',glowers=[('ab')])`
```math
\textrm{tr}\,\left[\gamma^+\gamma^-\gamma^{\perp a}\gamma^{\perp b}\right]g_{\perp a\perp b} = 8
```

With 4-vectors and mu-nu contraction:

`fulltrace(['+',MULABEL,kv,lv,NULABEL,'-'])`
```math
\textrm{tr}\,\left[\gamma^+\gamma^\mu k\!\!\!\!/ l\!\!\!\!\!\:/ \gamma^\nu \gamma^-\right]g_{\mu\nu} = -16(l_{\perp}\cdot k_{\perp}) + 16k^+p_2^-y
```

With 4-vectors, mu-nu contraction, and a gluon propagator, diagram 12b

`fulltrace(['-','?\\beta',lv-p1v,NULABEL,p2v+kv,'-','+','-',p2v+kv,MULABEL,lv-p1pv,'?\\gamma'], Gterms=[ (lv,('?\\gamma','?\\beta')) ])`
```math
\textrm{tr}\,\left[\gamma^-\gamma^{\beta}(l\!\!\!\!\!\:/-p_1\!\!\!\!\!\!\!/)\gamma^\nu(p_2\!\!\!\!\!\!\!/+k\!\!\!\!/)\gamma^-\gamma^+\gamma^-(p_2\!\!\!\!\!\!\!/+k\!\!\!\!/)\gamma^\mu(l\!\!\!\!\!\:/-p_1'\!\!\!\!\!\!\!/\gamma^\gamma) \right]G_{\gamma\beta}(l) g_{\mu\nu} = -64(p_2^-)^2y^2(k_{\perp} + p_{2\perp})^2 + 128(p_2^-)^2y(k_{\perp} + p_{2\perp})(l_{\perp} - p_{1\perp}) - 64(p_2^-)^2(l_{\perp} - p_{1\perp})^2
```