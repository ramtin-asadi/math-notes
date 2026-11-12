# Numerical Methods — Part 1

A mathematical formula and a reliable numerical procedure are different objects. These notes track the approximation error, the sensitivity of the underlying problem, and the behavior of the numerical method separately.

Sources:

- Tobin Driscoll and Richard Braun, [Fundamentals of Numerical Computation, online Python edition](https://fncbook.com/python). Exercise numbers below refer to the online sections linked individually.
- Laurent Demanet, [MIT 18.330, Introduction to Numerical Analysis, Spring 2012](https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/), notes and assignments.
- Boyd and Vandenberghe, [Convex Optimization](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf), Chapter 9, for descent methods.

Small computations are shown by hand. Selected source exercises identify theoretical subparts when a larger source problem also requests programming. No software project is needed to follow the notes.

[Contents](README.md) · [Part 2](07-numerical-methods-part-2.md)

Contents:

- [1. Error, floating point, and conditioning](#1-error-floating-point-and-conditioning)
- [2. Numerical linear algebra](#2-numerical-linear-algebra)
- [3. Root finding](#3-root-finding)
- [4. Numerical optimization](#4-numerical-optimization)
- [5. Differentiation and integration](#5-differentiation-and-integration)
- [6. Numerical ordinary differential equations](#6-numerical-ordinary-differential-equations)

## 1. Error, floating point, and conditioning

### 1.1 Theory

If $`x`$ is exact and $`\widehat x`$ is an approximation, absolute error is $`\vert \widehat x-x\vert`$ and relative error is that quantity divided by $`\vert x\vert`$, when $`x\ne0`$. A relative error at a zero answer is undefined; choose an absolute scale instead.

A floating-point number stores a sign, a finite significand, and an exponent. Binary64 has 53 binary digits of significand precision, including the leading bit for normal numbers. The spacing immediately above one is $`2^{-52}`$, often called machine epsilon. Unit roundoff for rounding to nearest is $`u=2^{-53}`$. Naming the convention avoids a factor-of-two ambiguity.

For an ordinary operation whose exact result remains in the normal representable range, the usual model is $`\mathop{\mathrm{fl}} (a\circ b)=(a\circ b)(1+\delta)`$ with $`\vert \delta\vert \le u`$. Overflow, underflow, subnormal results, division by zero, and invalid operations require separate treatment.

Truncation error comes from replacing an exact mathematical operation with a finite approximation, such as a finite Taylor sum or a difference quotient. Roundoff error comes from finite representation and arithmetic. Refining an approximation can decrease truncation error while increasing the influence of roundoff.

Conditioning belongs to the mathematical input-output problem. Stability belongs to an algorithm used to solve it. Forward error compares the computed output with the true output. Backward error asks for a nearby input problem that the computed output solves exactly. A backward-stable method can have a large forward error on an ill-conditioned problem.

### 1.2 Derivations and proofs

#### Scalar relative condition number

For differentiable $`f`$ and nonzero $`x,f(x)`$, a small relative input perturbation $`x\varepsilon`$ gives $`f(x+x\varepsilon)-f(x)=f'(x)x\varepsilon+o(\varepsilon)`$. Hence

```math
\kappa_f(x)=\left\vert \frac{xf'(x)}{f(x)}\right\vert .
```

The absolute condition number is $`\vert f'(x)\vert`$. A large relative number can arise because the output is close to zero even when the absolute sensitivity is moderate.

#### Cancellation and reformulation

For small $`x`$, evaluating $`\sqrt{1+x}-1`$ subtracts two numbers near one to produce a much smaller answer. Rationalize:

```math
\sqrt{1+x}-1=\frac{x}{\sqrt{1+x}+1}.
```

The exact problem is well-conditioned near zero in relative terms, since the output behaves like $`x/2`$. The direct algorithm can first round $`1+x`$ to one and return zero. The rationalized form retains the input in the numerator and divides by a quantity near two.

Subtraction itself may be exactly representable; the danger is amplification of errors already present in two nearly equal operands. A cancellation diagnosis should explain where those operand errors came from.

#### Unit roundoff and accumulated summation error

In binary64 arithmetic with rounding to nearest, the gap from one to the next larger representable number is $`2^{-52}`$. The unit roundoff used in standard relative-error bounds is half that gap, $`u=2^{-53}`$. For an elementary operation whose exact result remains in the normal range, a common model is $`\mathop{\mathrm{fl}} (a\circ b)=(a\circ b)(1+\delta)`$ with $`\vert \delta\vert \le u`$. Overflow, underflow, and exceptional inputs require separate treatment.

Repeated rounding errors do not simply add as signed relative errors. A useful product bound is $`\prod_{i=1}^n(1+\delta_i)=1+\theta_n`$ with $`\vert \theta_n\vert \le\gamma_n=nu/(1-nu)`$ when $`nu\lt 1`$. Expanding the products or bounding them by a geometric series gives this estimate.

Sequential summation then satisfies

```math
\left\vert \mathop{\mathrm{fl}} \left(\sum_{i=1}^n x_i\right)-\sum_{i=1}^n x_i\right\vert
\le\gamma_{n-1}\sum_{i=1}^n\vert x_i\vert .
```

The relative error bound divides by the magnitude of the exact sum. Its condition factor is $`\sum\vert x_i\vert /\vert \sum x_i\vert`$ when the sum is nonzero. Cancellation can make this factor large. If the exact sum is zero, relative error against that sum is undefined and an absolute scale must be used.

Pairwise summation arranges additions in a balanced tree. Each input passes through only about $`\lceil\log_2n\rceil`$ rounding stages, producing a corresponding smaller depth-based bound. This changes the algorithmic error without changing the mathematical sum or its conditioning.

#### Example NM-E17: stable quadratic roots

Solve $`x^2+10^8x+1=0`$ accurately.

Solution. The direct formula computes one root as $`(-10^8+\sqrt{10^{16}-4})/2`$. Its numerator subtracts nearly equal numbers, and the small root can lose many digits. Instead set

```math
q=-\frac12\left(b+\mathop{\mathrm{sign}} (b)\sqrt{b^2-4ac}\right),
\qquad x_1=q/a,\qquad x_2=c/q.
```

Here $`q`$ uses addition of same-sign magnitudes, and the second root uses the product relation $`x_1x_2=c/a`$. The roots are approximately $`-10^8`$ and $`-10^{-8}`$.

The formula assumes a nonzero $`q`$ and a suitable real discriminant for the real-root problem. Degenerate cases such as zero coefficients must be handled from their actual equation. The correction is an algebraic reformulation, not a change to the polynomial being solved.

### 1.3 Worked examples

#### Example NM-E01: a lost unit

In binary64 round-to-nearest, adding one to $`2^{53}`$ lies exactly between adjacent representable numbers. Ties-to-even rounds it back to $`2^{53}`$, so the computed expression $`(2^{53}+1)-2^{53}`$ can equal zero rather than one. The first addition loses the unit; the final subtraction cannot recover information no longer stored.

#### Example NM-E02: relative sensitivity near a zero

For $`f(x)=\log x`$, $`\kappa_f(x)=1/\vert \log x\vert`$. At $`x`$ near one the output is near zero, making relative sensitivity large. The absolute condition number $`1/x`$ is near one. An absolute error can remain small while its fraction of a tiny exact logarithm is large.

### 1.4 Selected exercises

#### Problem NM-01: four elementary condition numbers

Source: [FNC, Exercise 1.2.1(a–d)](https://fncbook.com/python/conditioning/).

Problem. Find the relative condition numbers of a power, logarithm, cosine, and exponential, specifying where the relative definitions apply.

Solution. Apply $`\vert xf'/f\vert`$ on domains where it is defined. For $`f=x^p`$ on $`x\gt 0`$, it is $`\vert p\vert`$. For $`\log x`$, with $`x\gt 0`$ and $`x\ne1`$, it is $`1/\vert \log x\vert`$. For $`\cos x`$ away from its zeros, it is $`\vert x\tan x\vert`$. For $`e^x`$ at nonzero input, it is $`\vert x\vert`$. At zero input the relative-input definition needs a limiting convention; the corresponding absolute condition numbers remain well-defined. These are local sensitivities, not finite-perturbation error guarantees of arbitrary size.

#### Problem NM-02: composed functions

Source: [FNC, Exercise 1.2.2(a–c)](https://fncbook.com/python/conditioning/).

Problem. Find the relative condition numbers of $`\sqrt{x+5}`$, $`\cos(2\pi x)`$, and $`e^{-x^2}`$ and identify excluded inputs or outputs.

Solution. Direct differentiation gives, respectively,

```math
\kappa_{\sqrt{x+5}}=\left\vert \frac{x}{2(x+5)}\right\vert ,\quad x\gt -5,
```

```math
\kappa_{\cos(2\pi x)}=\vert 2\pi x\tan(2\pi x)\vert ,\qquad
\kappa_{e^{-x^2}}=2x^2.
```

The cosine expression excludes zero outputs; the square-root expression excludes its boundary. Applying the chain rule for condition numbers gives the same expressions wherever the intermediate relative quantities are defined. At an intermediate zero, compute the direct derivative or take a justified limit instead of multiplying undefined factors.

#### Problem NM-03: the chain rule for conditioning

Source: [FNC, Exercise 1.2.4](https://fncbook.com/python/conditioning/).

Problem. Derive the relative condition-number chain rule for a composition at points where all required ratios are defined.

Solution. For $`h=f\circ g`$, with all denominators nonzero,

```math
\kappa_h(x)=\left\vert \frac{xf'(g(x))g'(x)}{f(g(x))}\right\vert
=\left\vert \frac{g(x)f'(g(x))}{f(g(x))}\right\vert
\left\vert \frac{xg'(x)}{g(x)}\right\vert .
```

This is $`\kappa_f(g(x))\kappa_g(x)`$. Differentiability and nonzero intermediate input and output values justify the division. The result describes the composite mathematical map, not every numerical implementation of the composition.

#### Problem NM-04: a cancellation-heavy exponential series

Source: [MIT 18.330, Homework 1, Problem 4](https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/d8531883a96517b39f1b6640fa983d65_MIT18_330S12_hw1.pdf), explanation of the finite-precision failure.

Problem. Explain why direct finite-precision summation of the Taylor series for $`e^{-25}`$ can be inaccurate even after including many terms.

Solution. The exact series $`\sum_{n\ge0}(-25)^n/n!`$ equals $`e^{-25}\approx1.39\times10^{-11}`$. Its terms initially grow in magnitude, reaching sizes of order $`e^{25}/\sqrt{2\pi25}`$ near $`n=25`$, before shrinking. Positive and negative terms of order billions must cancel to leave an answer of order $`10^{-11}`$. Tiny relative errors in the large terms can overwhelm the result.

Adding more terms reduces truncation error but cannot recover digits lost through that cancellation. Computing $`e^{25}`$ with positive terms and then taking its reciprocal is a more stable conceptual route; a practical exponential implementation also uses range reduction and carefully controlled approximations.

### 1.5 Additional solved problems

#### Problem NM-05: a stable quadratic root formula

Constructed problem. For $`x^2-bx+1=0`$ with large $`b\gt 2`$, compute the small root without subtracting nearly equal numbers.

Solution. The direct small root $`(b-\sqrt{b^2-4})/2`$ suffers cancellation. The large root is $`r_+=(b+\sqrt{b^2-4})/2`$, and the product of the roots is one, so $`r_-=1/r_+`$. The latter avoids the dangerous subtraction. Separate overflow protection may still be needed when forming $`b^2`$ for extremely large $`b`$; an algebraically improved formula is not automatically safe over the entire machine range.

#### Problem NM-06: perturbing a repeated root

Constructed problem. Compare roots of $`(x-1)^2=0`$ and $`(x-1)^2-\varepsilon=0`$ for $`\varepsilon\gt 0`$.

Solution. The perturbed roots are $`1\pm\sqrt\varepsilon`$. A coefficient perturbation of order $`\varepsilon`$ produces root changes of order $`\sqrt\varepsilon`$, whose ratio to $`\varepsilon`$ diverges. The repeated-root problem is ill-conditioned. A stable algorithm cannot remove this mathematical sensitivity.

## 2. Numerical linear algebra

### 2.1 Theory

To solve $`Ax=b`$, factor the matrix and solve triangular systems rather than explicitly forming an inverse. LU elimination writes $`PA=LU`$ with a permutation $`P`$ for pivoting. Forward substitution solves $`Ly=Pb`$ and backward substitution solves $`Ux=y`$. A zero or tiny pivot can require a row swap; nonsingularity of the full matrix does not guarantee every unpivoted leading pivot is safe.

For symmetric PD $`A`$, Cholesky writes $`A=LL^\top`$ with positive diagonal entries. Positive definiteness guarantees the square roots used in the factorization remain positive in exact arithmetic. A semidefinite matrix can require pivoting or a rank-revealing factorization.

For full-column-rank least squares, thin QR writes $`A=QR`$ with $`Q^\top Q=I`$ and square invertible upper-triangular $`R`$. The minimizer solves $`R\widehat x=Q^\top b`$. Orthogonal transformations preserve Euclidean norms and avoid explicitly squaring the condition number through $`A^\top A`$.

Householder reflection $`H=I-2vv^\top/(v^\top v)`$ is orthogonal and can zero a column's lower entries in one transformation. Successive reflections produce QR. Modified Gram-Schmidt is another construction, but strong near-dependence can require reorthogonalization.

### 2.2 Derivations and proofs

#### QR least squares without normal equations

Complete the thin $`Q`$ to a square orthogonal matrix $`(Q,Q_\perp)`$. Then

```math
\Vert Ax-b\Vert ^2=\Vert Rx-Q^\top b\Vert ^2+\Vert Q_\perp^\top b\Vert ^2.
```

The second term is fixed. The first becomes zero at the unique triangular-system solution. This is the same geometric projection as the normal equations, expressed in numerically safer coordinates.

#### Residual and forward error

If $`A`$ is invertible, define $`r=b-A\widehat x`$. Then $`A(x-\widehat x)=r`$, so $`\Vert x-\widehat x\Vert \le\Vert A^{-1}\Vert \Vert r\Vert`$. Since $`\Vert b\Vert \le\Vert A\Vert \Vert x\Vert`$,

```math
\frac{\Vert x-\widehat x\Vert }{\Vert x\Vert }
\le\kappa(A)\frac{\Vert r\Vert }{\Vert b\Vert }.
```

The residual also means $`\widehat x`$ exactly solves $`A\widehat x=b-r`$, identifying a backward perturbation in the right-hand side. A small relative residual alone does not guarantee a small forward error.

#### Why normal equations square a condition number

For full-column-rank $`A`$, the singular values of $`A^\top A`$ are the squares of those of $`A`$. Hence $`\kappa_2(A^\top A)=\kappa_2(A)^2`$. The algebraic normal equations remain correct, but forming the Gram matrix can discard small distinctions between nearly dependent columns.

A thin QR factorization $`A=QR`$ gives $`\Vert Ax-b\Vert ^2=\Vert Rx-Q^\top b\Vert ^2+\Vert (I-QQ^\top)b\Vert ^2`$. Solving the triangular system uses the geometry of $`A`$ directly. Householder QR is backward stable under the standard arithmetic model, whereas explicitly forming normal equations can magnify the effect of rounding through the squared condition number.

Cholesky is appropriate for a supplied SPD matrix, but applying it to a rounded Gram matrix does not recover information already lost in forming that matrix. The Cholesky recurrence is

```math
L_{jj}=\sqrt{a_{jj}-\sum_{k\lt j}L_{jk}^2},\qquad
L_{ij}=\frac{a_{ij}-\sum_{k\lt j}L_{ik}L_{jk}}{L_{jj}}.
```

Positive Schur complements of an SPD matrix make the exact radicands positive. A computed nonpositive radicand can reflect loss of definiteness through rounding or a matrix that was never SPD.

#### Example NM-E18: an independent pair hidden by its Gram matrix

Take

```math
A=\begin{pmatrix}1&1\\[0pt] \epsilon&0\\[0pt] 0&\epsilon\end{pmatrix},
\qquad b=\begin{pmatrix}0\\[0pt] \epsilon\\[0pt] -\epsilon\end{pmatrix},
\qquad \epsilon=10^{-8}.
```

Solution. The exact solution is $`x=(1,-1)^\top`$. The columns are independent. Their Gram matrix is

```math
A^\top A=
\begin{pmatrix}1+\epsilon^2&1\\[0pt] 1&1+\epsilon^2\end{pmatrix},
```

with eigenvalues $`2+\epsilon^2`$ and $`\epsilon^2`$. Thus $`\kappa_2(A)\approx\sqrt2\,10^8`$, while the Gram condition number is about $`2\cdot10^{16}`$.

In binary64, adding $`10^{-16}`$ to one rounds to one. Directly formed diagonal entries can therefore equal the off-diagonal entries, making the rounded Gram matrix singular. The original rectangular matrix still stores the small entries explicitly.

QR avoids this specific loss, though it cannot remove the original condition number. A backward error of machine scale can still create a forward error amplified by about $`10^8`$. Stability of the solver and conditioning of the supplied data remain separate issues.

### 2.3 Worked examples

#### Example NM-E03: Cholesky by hand

For

```math
A=\left(\begin{smallmatrix}4&2\\[0pt] 2&3\end{smallmatrix}\right)
```

, write

```math
L=\left(\begin{smallmatrix}l_{11}&0\\[0pt] l_{21}&l_{22}\end{smallmatrix}\right)
```

. Matching entries gives $`l_{11}=2`$, $`l_{21}=1`$, and $`l_{22}=\sqrt2`$. All roots are positive. To solve $`Ax=(6,5)^\top`$, forward substitution gives $`y=(3,\sqrt2)`$ and backward substitution gives $`x=(1,1)`$.

#### Example NM-E04: why pivoting is useful

For

```math
A=\left(\begin{smallmatrix}\varepsilon&1\\[0pt] 1&1\end{smallmatrix}\right)
```

 with small nonzero $`\varepsilon`$, unpivoted elimination uses multiplier $`1/\varepsilon`$, producing large intermediate entries. Swapping rows first uses multiplier $`\varepsilon`$, keeping them moderate. The original matrix remains well-conditioned as $`\varepsilon\to0`$; the large intermediates are a feature of the unpivoted algorithm.

### 2.4 Selected exercises

#### Problem NM-07: orthogonal-matrix properties

Source: [FNC, Exercise 3.3.2](https://fncbook.com/python/qr/).

Problem. For square orthogonal matrices, prove the inverse, transpose, product, norm-preservation, and condition-number identities and show invariance of $`\Vert AQ\Vert _2`$.

Solution. For square $`Q`$ with $`Q^\top Q=I`$, invertibility gives $`Q^{-1}=Q^\top`$, and hence $`QQ^\top=I`$, proving $`Q^\top`$ orthogonal. For any $`x`$, $`\Vert Qx\Vert _2^2=x^\top Q^\top Qx=\Vert x\Vert _2^2`$, so $`\Vert Q\Vert _2=\Vert Q^{-1}\Vert _2=1`$ and $`\kappa_2(Q)=1`$. The substitution $`y=Qx`$ maps the unit sphere onto itself, so $`\Vert AQ\Vert _2=\sup_{\Vert y\Vert =1}\Vert Ay\Vert _2=\Vert A\Vert _2`$. Finally $`(QU)^\top QU=U^\top U=I`$ for orthogonal $`U`$, so their product is orthogonal.

#### Problem NM-08: invertibility of the thin triangular factor

Source: [FNC, Exercise 3.3.4](https://fncbook.com/python/qr/).

Problem. For a thin QR factorization of a full-column-rank matrix, prove that its square triangular factor is invertible.

Solution. Suppose full-column-rank $`A=QR`$ had singular square $`R`$. There would be $`z\ne0`$ with $`Rz=0`$, giving $`Az=QRz=0`$, contrary to column independence. Hence $`R`$ is invertible. Conversely, invertible $`R`$ and orthonormal $`Q`$ make $`A`$ injective, since $`QRz=0`$ implies $`Rz=Q^\top QRz=0`$.

#### Problem NM-09: a pseudoinverse from QR

Source: [FNC, Exercise 3.3.6](https://fncbook.com/python/qr/).

Problem. Derive the pseudoinverse of a full-column-rank matrix from its thin QR factorization and verify the Moore-Penrose identities.

Solution. For thin full-rank QR, every least-squares fit solves $`Rx=Q^\top b`$, so $`x=R^{-1}Q^\top b`$. This defines $`A^+=R^{-1}Q^\top`$. To check directly, $`A^+A=I`$ and $`AA^+=QQ^\top`$, a symmetric projector. The two remaining Moore-Penrose identities follow by multiplying these products by $`A`$ and $`A^+`$ respectively.

#### Problem NM-10: a triangular inverse with large norm

Source: [FNC, Exercise 2.8.6(a–c)](https://fncbook.com/python/condition-number/).

Problem. Let $`A_n=I-2N`$ where $`N`$ has ones on the first superdiagonal. Find its inverse and infinity-norm condition number for $`n\gt 1`$.

Solution. Let $`A_n=I-2N`$, where $`N`$ has ones on its first superdiagonal and $`N^n=0`$. The finite geometric identity gives

```math
A_n^{-1}=I+2N+4N^2+\cdots+2^{n-1}N^{n-1}.
```

Thus $`(A_n^{-1})_{ij}=2^{j-i}`$ for $`j\ge i`$ and zero below. In particular

```math
A_2^{-1}=\left(\begin{smallmatrix}1&2\\[0pt] 0&1\end{smallmatrix}\right)
```

 and

```math
A_3^{-1}=\left(\begin{smallmatrix}1&2&4\\[0pt] 0&1&2\\[0pt] 0&0&1\end{smallmatrix}\right)
```

. The largest inverse row sum is $`2^n-1`$, while $`\Vert A_n\Vert _\infty=3`$ for $`n\gt 1`$. Hence $`\kappa_\infty(A_n)=3(2^n-1)`$: a triangular matrix can be exponentially ill-conditioned.

### 2.5 Additional solved problems

#### Problem NM-11: conditioning of normal equations

Constructed problem. If full-column-rank $`A`$ has singular values one and $`10^{-4}`$, compare its condition number with that of $`A^\top A`$.

Solution. The first is $`10^4`$. The Gram matrix has eigenvalues one and $`10^{-8}`$, so its condition number is $`10^8`$. Forming the Gram matrix can make finite-precision loss more severe. Solving by QR avoids that explicit squaring, although the original problem still contains a weak direction with singular value $`10^{-4}`$.

#### Problem NM-12: iterative refinement as an error equation

Constructed problem. Given an approximate $`\widehat x`$, explain one refinement step.

Solution. Compute $`r=b-A\widehat x`$, solve $`A\delta=r`$, and set $`\widehat x_{\mathrm{new}}=\widehat x+\delta`$. In exact arithmetic this reaches the exact solution in one step because $`\delta=x-\widehat x`$. In finite precision the correction solve and residual computation are approximate. Higher-precision residuals and a sufficiently accurate factorization can improve the result, but severe conditioning can limit the attainable accuracy.

## 3. Root finding

### 3.1 Theory

Bisection applies to a continuous real function with opposite signs at the endpoints. The intermediate value theorem guarantees a root. Halve the interval and retain a half with a sign change. It gives a reliable location bound without derivatives, although it does not exploit a particularly smooth root.

Fixed-point iteration rewrites a root equation as $`x=g(x)`$ and iterates $`x_{k+1}=g(x_k)`$. A closed invariant interval and a uniform derivative bound $`\vert g'\vert \le q\lt 1`$ give a contraction and convergence. A derivative bound without invariance or completeness is insufficient for a global existence claim.

Newton uses the tangent line: $`x_{k+1}=x_k-f(x_k)/f'(x_k)`$. It is locally quadratic near a sufficiently smooth simple root, but a poor starting point, a small derivative, or a multiple root can change its behavior. Bracketing safeguards combine reliability with fast local steps.

### 3.2 Derivations and proofs

#### Bisection error

Starting with width $`b-a`$, after $`n`$ halvings the retained interval has width $`(b-a)2^{-n}`$. Its midpoint $`m_n`$ lies at distance at most half that width from every point of the interval, including its retained root:

```math
\vert m_n-r\vert \le\frac{b-a}{2^{n+1}}.
```

This bound assumes exact sign decisions and arithmetic. In finite precision, an interval can become too narrow for its computed midpoint to differ from an endpoint.

#### Local Newton convergence

Assume $`f(r)=0`$, $`\vert f'\vert \ge m\gt 0`$, and $`\vert f''\vert \le M`$ in a neighborhood of $`r`$. Taylor expansion from $`x_k`$ to $`r`$ gives

```math
0=f(x_k)+f'(x_k)(r-x_k)+\frac12f''(\xi_k)(r-x_k)^2.
```

Rearrange after the Newton update:

```math
\vert x_{k+1}-r\vert \le\frac{M}{2m}\vert x_k-r\vert ^2.
```

Choose the initial error small enough that the bound remains inside the neighborhood and the factor times the initial error is below one. Induction then keeps the iterates there and proves convergence. The derivative lower bound rules out a multiple root and prevents division by zero.

#### A local Newton theorem with a neighborhood check

Let $`r`$ be a simple root, and suppose $`f`$ is twice continuously differentiable on a closed interval around it, with $`\vert f'(x)\vert \ge m\gt 0`$ and $`\vert f''(x)\vert \le M`$. Taylor-expand $`f(r)`$ about the current iterate $`x`$:

```math
0=f(x)+f'(x)(r-x)+\frac12f''(\xi)(r-x)^2.
```

For the Newton update $`N(x)=x-f(x)/f'(x)`$, this gives

```math
\vert N(x)-r\vert \le\frac M{2m}\vert x-r\vert ^2.
```

Choose a radius $`\delta`$ inside the regularity interval such that $`C\delta\lt 1`$, where $`C=M/(2m)`$. If $`\vert x_0-r\vert \le\delta`$, the estimate keeps the next iterate in the same interval and reduces the error by at least factor $`C\delta`$. Induction proves convergence and the quadratic error estimate.

The neighborhood check closes the proof: it prevents an estimate requiring local derivative bounds from being applied after an iterate has left their domain. A formula for the asymptotic error constant alone does not supply that invariance.

#### Multiple roots reduce the order

If $`f(x)=(x-r)^m g(x)`$ with $`g(r)\ne0`$, then near the root $`f/f'=(x-r)/m+O((x-r)^2)`$. Ordinary Newton therefore has

```math
x_{n+1}-r=\frac{m-1}{m}(x_n-r)+O((x_n-r)^2).
```

For multiplicity greater than one, the leading error is linear. If the multiplicity is known, the modified step $`x-mf(x)/f'(x)`$ cancels that leading term and can restore quadratic convergence under the corresponding smoothness conditions.

#### Example NM-E19: a triple root

For $`f(x)=(x-1)^3`$, Newton's method gives $`x_{n+1}=(2x_n+1)/3`$ whenever $`x_n\ne1`$. Thus $`e_{n+1}=(2/3)e_n`$ exactly. Starting at four gives $`e_n=3(2/3)^n`$.

The step size tends to zero, and the residual $`e_n^3`$ tends to zero faster than the location error. A residual tolerance can therefore look impressive even while the root error is substantially larger. The derivative also tends to zero, violating the simple-root lower bound used above.

Multiplying the Newton correction by three lands at one in a single step for this pure power. A general multiple-root function has additional higher-order terms, so it need not terminate in one step even with the modification.

### 3.3 Worked examples

#### Example NM-E05: Newton for a square root

For $`f(x)=x^2-2`$, the update is $`x_{k+1}=(x_k+2/x_k)/2`$. Starting from $`3/2`$ gives $`17/12`$, then $`577/408`$. The exact error identity is $`x_{k+1}-\sqrt2=(x_k-\sqrt2)^2/(2x_k)`$, displaying quadratic convergence and showing that positive iterates stay at or above the positive root after the first step.

#### Example NM-E06: a bisection count

For $`f(x)=x^3-2`$ on $`[1,2]`$, require midpoint error at most $`10^{-3}`$. The bound is $`2^{-(n+1)}`$. Since $`2^{-10}\lt 10^{-3}`$, nine halvings suffice for this midpoint convention. Specifying whether $`n`$ counts midpoint evaluations or completed halvings avoids an indexing ambiguity.

### 3.4 Selected exercises

#### Problem NM-13: local fixed-point convergence

Source: [FNC, Exercise 4.2.1(a–c)](https://fncbook.com/python/fixed-point/).

Problem. Check local convergence of the fixed-point maps $`1+x-x^2/9`$ at three, $`\pi+\sin x/4`$ at $`\pi`$, and $`x+1-\tan(x/4)`$ at $`\pi`$.

Solution. For $`g=1+x-x^2/9`$ at $`p=3`$, $`g(3)=3`$ and $`g'(3)=1/3`$. For $`g=\pi+\sin x/4`$ at $`p=\pi`$, the derivative is $`-1/4`$. For $`g=x+1-\tan(x/4)`$ at $`p=\pi`$, it is $`1-\sec^2(\pi/4)/4=1/2`$. All three absolute derivatives are below one. Continuity of the derivatives gives a small invariant contraction neighborhood around each fixed point, so sufficiently close starts converge. The negative derivative in the second case makes small errors alternate signs.

#### Problem NM-14: root extraction by a poor fixed-point map

Source: [FNC, Exercise 4.2.5](https://fncbook.com/python/fixed-point/).

Problem. Examine whether the iteration $`x_{k+1}=a/x_k^{m-1}`$ locally computes the positive $`m`$th root of $`a\gt 0`$ for integer $`m\gt 1`$.

Solution. For $`g(x)=a/x^{m-1}`$ and positive root $`p=a^{1/m}`$, $`g'(p)=-(m-1)`$. For $`m\gt 2`$ the fixed point is repelling, so nearby nonexact starts generally move away. At $`m=2`$ the derivative is $`-1`$, and the derivative test alone is inconclusive. But $`g(g(x))=x`$, producing a two-cycle unless $`x=p`$. Thus none of the integer values $`m\gt 1`$ yields local convergence from all sufficiently close starts, apart from starting at the exact root itself.

#### Problem NM-15: choosing a constant derivative scale

Source: [FNC, Exercise 4.2.7](https://fncbook.com/python/fixed-point/).

Problem. For $`g(x)=x-f(x)/c`$ near a simple root $`p`$ with $`f'(p)\gt 0`$, determine the constant scales $`c`$ giving a local contraction.

Solution. For $`g(x)=x-f(x)/c`$, a root $`p`$ is a fixed point and $`g'(p)=1-f'(p)/c`$. With $`f'(p)\gt 0`$, the condition $`\vert 1-f'(p)/c\vert \lt 1`$ is equivalent to $`c\gt f'(p)/2`$. Derivative continuity then supplies a local contraction neighborhood. A negative $`c`$ fails this criterion. Choosing $`c=f'(p)`$ makes the derivative zero, but the true derivative at the unknown root is not usually available in advance.

#### Problem NM-16: division-free reciprocal Newton iteration

Source: [FNC, Exercise 4.3.5](https://fncbook.com/python/newton/).

Problem. Derive a division-free Newton iteration for the reciprocal of $`b\ne0`$ and prove its error recursion and convergence condition.

Solution. For $`f(x)=1/x-b`$, $`f'=-1/x^2`$. Newton simplifies to $`x_{k+1}=2x_k-bx_k^2=x_k(2-bx_k)`$, using only multiplication and subtraction. If $`e_k=1-bx_k`$, then $`e_{k+1}=e_k^2`$. Thus $`\vert e_0\vert \lt 1`$ gives quadratic convergence to the reciprocal. The iteration formula can be evaluated without division, though a suitable initial scale is still required.

#### Problem NM-17: Newton can cycle

Source: [FNC, Exercise 4.3.6](https://fncbook.com/python/newton/).

Problem. Analyze Newton's method for $`f(x)=\mathop{\mathrm{sign}} (x)\sqrt{\vert x\vert }`$ from a nonzero starting value.

Solution. For $`f(x)=\mathop{\mathrm{sign}} (x)\sqrt{\vert x\vert }`$ and $`x\ne0`$, $`f'(x)=1/(2\sqrt{\vert x\vert })`$. Hence $`f/f'=2x`$ and Newton gives $`x_{k+1}=-x_k`$. Every nonzero start alternates forever between two points. The root is zero, where the function is not differentiable with a finite derivative. The smooth simple-root theorem does not apply.

### 3.5 Additional solved problems

#### Problem NM-18: small steps without convergence

Source: [FNC, Exercise 4.3.8](https://fncbook.com/python/newton/).

Problem. Give a sequence whose consecutive steps tend to zero but which fails to converge, and explain the implication for an iterative stopping rule.

Solution. Let $`x_k=\sum_{j=1}^k1/j`$. Then $`x_k-x_{k-1}=1/k\to0`$, but the harmonic sums diverge. For example, each block from $`2^{m-1}+1`$ through $`2^m`$ contributes at least $`1/2`$, so the total is unbounded. A small iterate change does not alone certify proximity to a limit or a root.

## 4. Numerical optimization

### 4.1 Theory

Gradient descent uses $`x_{k+1}=x_k-\alpha\nabla f(x_k)`$. A positive step size does not automatically give descent; it must be small enough relative to curvature or chosen by a line search. If the gradient is $`L`$-Lipschitz, Taylor's integral bound gives

```math
f(x-\alpha\nabla f(x))
\le f(x)-\alpha(1-\alpha L/2)\Vert \nabla f(x)\Vert ^2.
```

Thus $`0\lt \alpha\lt 2/L`$ ensures a strict decrease whenever the gradient is nonzero. For general nonconvex functions, decreasing values need not imply finding a global optimum.

Projected gradient uses $`x_{k+1}=\Pi_C(x_k-\alpha\nabla f(x_k))`$ for a closed convex feasible set $`C`$. Its fixed-point condition is equivalent to $`\nabla f(x^*)^\top(y-x^*)\ge0`$ for every $`y\in C`$. This is the first-order optimality condition for a differentiable convex objective.

### 4.2 Derivations and proofs

#### Gradient descent on a PD quadratic

For $`f(x)=\tfrac12x^\top Hx-b^\top x`$ with $`H`$ PD, let $`x^*=H^{-1}b`$. The error satisfies $`e_{k+1}=(I-\alpha H)e_k`$. In an orthonormal eigenbasis, each error coordinate is multiplied by $`1-\alpha\lambda_i`$.

Convergence for all starts therefore holds exactly when $`\vert 1-\alpha\lambda_i\vert \lt 1`$ for every eigenvalue, or $`0\lt \alpha\lt 2/\lambda_{\max}`$. The optimal fixed step for the worst Euclidean contraction factor is $`\alpha^*=2/(\lambda_{\min}+\lambda_{\max})`$, giving factor $`(\kappa-1)/(\kappa+1)`$. Ill-conditioning slows descent because no single step can move equally efficiently in both steep and flat directions.

#### Projection characterization

If $`p=\Pi_C(z)`$, differentiating the squared distance along the segment from $`p`$ to any $`y\in C`$ gives $`(z-p)^\top(y-p)\le0`$. Conversely this inequality makes $`\Vert z-y\Vert ^2=\Vert z-p\Vert ^2+\Vert y-p\Vert ^2-2(z-p)^\top(y-p)\ge\Vert z-p\Vert ^2`$. Apply it to $`z=x^*-\alpha\nabla f(x^*)`$ and $`p=x^*`$ to obtain the fixed-point optimality condition.

#### The descent lemma and a provable step

Suppose $`\nabla f`$ is Lipschitz with constant $`L`$ on the relevant convex domain. Integrating the gradient along the segment from $`x`$ to $`x+h`$ gives

```math
f(x+h)\le f(x)+\nabla f(x)^\top h+\frac L2\Vert h\Vert ^2.
```

For $`h=-\eta\nabla f(x)`$,

```math
f(x-\eta\nabla f(x))
\le f(x)-\eta(1-L\eta/2)\Vert \nabla f(x)\Vert ^2.
```

Thus every nonstationary step decreases the objective for $`0\lt \eta\lt 2/L`$. With $`\eta=1/L`$, the decrease is at least $`\Vert \nabla f(x)\Vert ^2/(2L)`$.

If $`f`$ is also $`\mu`$-strongly convex, minimizing its quadratic lower bound gives $`f(x)-f^*\le\Vert \nabla f(x)\Vert ^2/(2\mu)`$. Combining yields $`f(x_{k+1})-f^*\le(1-\mu/L)(f(x_k)-f^*)`$. This is a global linear-rate bound under both assumptions.

#### Projection supplies a constrained fixed-point condition

For a nonempty closed convex set $`C`$, $`p=\Pi_C(z)`$ is characterized by $`(z-p)^\top(y-p)\le0`$ for all $`y\in C`$. Apply the characterization to two projections and add the resulting inequalities to obtain

```math
\Vert \Pi_C(z)-\Pi_C(w)\Vert ^2
\le(\Pi_C(z)-\Pi_C(w))^\top(z-w).
```

Cauchy–Schwarz proves nonexpansiveness. Also $`x=\Pi_C(x-\eta\nabla f(x))`$ is equivalent to $`\nabla f(x)^\top(y-x)\ge0`$ for all feasible $`y`$. For convex $`f`$, this is precisely the global first-order optimality condition.

Projected gradient therefore checks a constrained stationarity condition. The gradient itself need not vanish at a boundary optimum.

#### Example NM-E20: projecting onto a simplex

Project $`z=(1.2,0.3,-0.2)`$ onto $`C=\lbrace x:x_i\ge0,\ \sum_i x_i=1\rbrace`$.

Solution. The KKT conditions for minimizing $`\tfrac12\Vert x-z\Vert ^2`$ give $`x_i=\max(z_i-\tau,0)`$ with a scalar $`\tau`$ chosen to make the coordinates sum to one.

Assume the first two coordinates remain positive and the third is zero. Then $`1.2-\tau+0.3-\tau=1`$ gives $`\tau=0.25`$. The proposed projection is $`(0.95,0.05,0)`$. Both active positive coordinates exceed zero, and the excluded coordinate has $`z_3-\tau=-0.45\le0`$, so the assumed active set is consistent.

Simply clipping negative entries and renormalizing would give $`(0.8,0.2,0)`$, a different vector. Euclidean projection subtracts a common threshold from the positive coordinates; it does not generally rescale them proportionally.

### 4.3 Worked examples

#### Example NM-E07: one fast and one slow direction

Let $`H=\mathop{\mathrm{diag}} (1,9)`$ and choose $`\alpha=0.2`$. Error multipliers are $`0.8`$ and $`-0.8`$, so both coordinates converge at rate $`0.8`$, with the second alternating. At $`\alpha=0.25`$, the second multiplier is $`-1.25`$ and generic iterates diverge. The optimal constant step is $`2/(1+9)=0.2`$.

#### Example NM-E08: projecting onto a half-line

Minimize $`(x+1)^2/2`$ over $`x\ge0`$. With step one, the unconstrained update is $`x-(x+1)=-1`$ for any start. Projection clips it to zero, the constrained optimum. At zero the gradient is one rather than zero; the constraint blocks the downhill negative direction.

### 4.4 Additional solved problems

#### Problem NM-19: exact line search on a quadratic

Constructed problem. For gradient $`g=Hx-b\ne0`$, minimize $`f(x-\alpha g)`$ over $`\alpha`$.

Solution. Expand the quadratic: $`f(x-\alpha g)=f(x)-\alpha g^\top g+\tfrac12\alpha^2g^\top Hg`$. The coefficient $`g^\top Hg`$ is positive, so the unique minimizer is $`\alpha=(g^\top g)/(g^\top Hg)`$. The next gradient is orthogonal to the old gradient because $`g^\top[g-\alpha Hg]=0`$. This does not mean consecutive parameter errors are orthogonal.

#### Problem NM-20: a box projection

Constructed problem. Project $`z=(-2,1/2,4)`$ onto $`[0,1]^3`$.

Solution. The squared distance is the sum of three coordinate squared distances, and the box constraints separate by coordinate. Each coordinate is clipped to the interval, giving $`(0,1/2,1)`$. The minimum squared distance is $`4+0+9=13`$. For a simplex or an ellipsoid, coordinatewise clipping generally fails because the constraints couple coordinates.

#### Problem NM-21: a simplex projection

Constructed problem. Project $`(2,1)`$ onto $`\lbrace x_1,x_2\ge0:x_1+x_2=1\rbrace`$.

Solution. Substitute $`x_2=1-x_1`$ with $`0\le x_1\le1`$. The distance objective is $`(x_1-2)^2+(-x_1)^2=2x_1^2-4x_1+4`$. Its unconstrained minimum is at one, which is feasible. The projection is $`(1,0)`$. This happens to be an endpoint; a general simplex projection finds a common threshold and clips, rather than clipping to a box and then normalizing arbitrarily.

#### Problem NM-22: a stopping bound using strong convexity

Constructed problem. If $`f`$ is $`m`$-strongly convex and differentiable, bound its optimality gap using the gradient.

Solution. Strong convexity gives $`f(y)\ge f(x)+\nabla f(x)^\top(y-x)+(m/2)\Vert y-x\Vert ^2`$. Minimize the right-hand side over all $`y`$; completing the square gives lower bound $`f(x)-\Vert \nabla f(x)\Vert ^2/(2m)`$ on $`f^*`$. Hence $`f(x)-f^*\le\Vert \nabla f(x)\Vert ^2/(2m)`$. The estimate is for unconstrained minimization; constrained stationary points need a projected-gradient or KKT residual instead.

## 5. Differentiation and integration

### 5.1 Theory

The forward difference $`[f(x+h)-f(x)]/h`$ approximates $`f'(x)`$ with first-order truncation error for a sufficiently smooth function. The centered difference $`[f(x+h)-f(x-h)]/(2h)`$ cancels even Taylor terms and has second-order error.

The step cannot usually be made arbitrarily small in floating-point arithmetic. If function evaluations have absolute errors of order $`u`$, division by $`h`$ can amplify those errors to order $`u/h`$. A centered derivative's total error often has the form $`C_1h^2+C_2u/h`$, suggesting a balanced scale $`h\asymp u^{1/3}`$ after accounting for the function's size and derivatives.

Quadrature approximates an integral by a weighted sum of function values. Composite trapezoidal integration is second-order for functions with bounded continuous second derivative. Composite Simpson is fourth-order under appropriate fourth-derivative smoothness and an even number of equal subintervals. Endpoint singularities can reduce these nominal orders.

Monte Carlo integration represents an integral as an expectation and averages random evaluations. For iid finite-variance draws, its standard deviation decreases as $`N^{-1/2}`$. This rate depends on variance rather than a spatial grid dimension, although high dimension can make the variance large or the integrand hard to sample.

### 5.2 Derivations and proofs

#### Difference truncation errors

Taylor expansion gives

```math
\frac{f(x+h)-f(x)}h
=f'(x)+\frac h2f''(x)+O(h^2),
```

```math
\frac{f(x+h)-f(x-h)}{2h}
=f'(x)+\frac{h^2}{6}f'''(x)+O(h^4).
```

The second expansion with its $`O(h^4)`$ remainder requires sufficient fifth-order smoothness. A weaker third-derivative bound still yields a second-order error bound without asserting all higher coefficients.

If each of two centered-function evaluations has absolute error at most $`\eta`$, their difference divided by $`2h`$ has perturbation at most $`\eta/h`$. Combining with a third-derivative bound $`M_3`$ gives $`M_3h^2/6+\eta/h`$. Differentiating this bound gives $`h^3=3\eta/M_3`$ when both terms are positive and the local smoothness region permits that step.

#### Trapezoidal error

On a single interval of length $`h`$, expand about its midpoint. The exact integral and endpoint trapezoid share the constant and odd contributions, while their quadratic contributions differ by $`h^3f''/12`$. With $`\vert f''\vert \le M_2`$, the local error magnitude is at most $`M_2h^3/12`$. Summing $`(b-a)/h`$ equal cells gives

```math
\vert T_h-I\vert \le\frac{(b-a)M_2h^2}{12}.
```

For sufficient higher smoothness, the leading signed term is $`T_h-I=(h^2/12)[f'(b)-f'(a)]+O(h^4)`$. This expansion supports extrapolation and endpoint correction.

#### Monte Carlo variance and a control variate

For $`I=E[f(U)]`$ and iid draws, $`\widehat I=N^{-1}\sum_if(U_i)`$ is unbiased with variance $`\mathop{\mathrm{Var}} (f(U))/N`$. If $`g`$ has known mean $`\mu_g`$, use

```math
\widehat I_c=\frac1N\sum_i[f(U_i)-c(g(U_i)-\mu_g)].
```

It remains unbiased for fixed $`c`$. Its per-draw variance is $`\mathop{\mathrm{Var}} (f)-2c\mathop{\mathrm{Cov}} (f,g)+c^2\mathop{\mathrm{Var}} (g)`$. Minimize this quadratic to obtain $`c^*=\mathop{\mathrm{Cov}} (f,g)/\mathop{\mathrm{Var}} (g)`$ when the denominator is positive. The minimized variance is $`\mathop{\mathrm{Var}} (f)(1-\rho_{fg}^2)`$. Estimating $`c`$ on the same sample requires an additional finite-sample analysis, though regular plug-in choices can retain first-order efficiency.

#### Quadrature error comes from interpolation error

On one interval of width $`h`$, the trapezoidal rule integrates the straight-line interpolant through the endpoint values. For a twice continuously differentiable function, the interpolation error and integration give a local trapezoidal excess $`h^3f''(\xi)/12`$. Summing equal panels yields

```math
T_h-\int_a^bf(x)dx=\frac{(b-a)h^2}{12}f''(\xi)
```

for some $`\xi\in(a,b)`$. In particular, the rule overestimates a convex function's integral. The sign is a useful check.

Composite Simpson integration on an even number of equal panels has error $`I-S_h=-(b-a)h^4f^{(4)}(\xi)/180`$ under a continuous fourth derivative. The stronger smoothness condition is essential. Applying a nominal fourth-order formula to an endpoint singularity does not create four derivatives there.

For centered differentiation, truncation is typically $`C_1h^2`$ while subtraction and division can produce arithmetic error of order $`C_2u/h`$. Balancing gives $`h`$ on the scale of $`u^{1/3}`$ after including the problem-dependent constants. This is a scale argument, not a universal numerical step for every function.

#### Example NM-E21: removing an integrable endpoint singularity

Evaluate $`\int_0^1x^{-1/2}dx`$ through a change of variables.

Solution. The integral is improper but finite. Set $`x=t^2`$ for $`t\in[0,1]`$. For positive $`t`$, $`dx=2t\,dt`$ and $`x^{-1/2}=1/t`$, so the transformed integrand is the constant two. Its continuous extension at zero is also two. Thus the integral is two.

An endpoint rule applied directly to $`x^{-1/2}`$ attempts to evaluate an infinite value at zero. The transformed rule has no singularity and integrates exactly even with one trapezoidal panel.

For $`\int_0^1\sqrt x\,dx`$, the same substitution produces $`\int_0^12t^2dt=2/3`$. The original derivative is singular at zero; the transformed integrand is a polynomial. A substitution can improve smoothness before discretization rather than asking a high-order rule to compensate for missing regularity.

The Jacobian and orientation of the substitution are parts of the integral. Replacing the function values without transforming $`dx`$ would solve a different problem.

### 5.3 Worked examples

#### Example NM-E09: a centered derivative of an exponential

For $`f(x)=e^x`$ at zero, the centered derivative is $`(e^h-e^{-h})/(2h)=\sinh(h)/h=1+h^2/6+O(h^4)`$. Halving $`h`$ reduces the leading truncation error by approximately four before roundoff becomes dominant.

#### Example NM-E10: integrating a square

For $`\int_0^1x^2dx=1/3`$, the composite trapezoid on $`n`$ equal cells satisfies $`T_n=1/3+1/(6n^2)`$. This follows by substituting $`\sum_{j=1}^{n-1}j^2=(n-1)n(2n-1)/6`$ into the quadrature sum. At $`n=2`$, the result is $`3/8`$ and error $`1/24`$.

#### Example NM-E11: controlling Monte Carlo variance

For $`U\sim U(0,1)`$, take $`f(U)=U^2`$ and $`g(U)=U`$, whose mean is one half. The covariance is $`EU^3-EU^2EU=1/4-1/6=1/12`$, equal to $`\mathop{\mathrm{Var}} (U)`$. Thus $`c^*=1`$. The raw per-draw variance is $`1/5-1/9=4/45`$, while the controlled variance is $`4/45-1/12=1/180`$, a factor-sixteen reduction.

### 5.4 Selected exercises

#### Problem NM-23: a derivative from a quadratic interpolant

Source: [FNC, Exercise 5.4.1(a,b)](https://fncbook.com/python/finitediffs/).

Problem. Construct the quadratic interpolant through values at $`-h,0,h`$ and differentiate it to obtain a centered derivative approximation.

Solution. With $`f_-=f(-h)`$, $`f_0=f(0)`$, $`f_+=f(h)`$, write

```math
Q(x)=f_0+\frac{f_+-f_-}{2h}x
+\frac{f_--2f_0+f_+}{2h^2}x^2.
```

Substitution at zero gives $`f_0`$. At $`h`$, the linear and quadratic contributions combine to $`f_+-f_0`$; at $`-h`$, they combine to $`f_--f_0`$. Thus it interpolates all three data. Differentiating at zero gives $`Q'(0)=(f_+-f_-)/(2h)`$, the centered difference. The approximation arises by differentiating the interpolating polynomial.

#### Problem NM-24: consistency conditions for difference weights

Source: [FNC, Exercise 5.4.7(a,b)](https://fncbook.com/python/finitediffs/).

Problem. Establish centered-difference consistency and derive the constant- and linear-moment conditions for a general finite first-derivative stencil.

Solution. For the centered formula, the numerator and denominator vanish at $`h=0`$. Differentiating with respect to $`h`$ gives limit $`[f'(x+h)+f'(x-h)]/2\to f'(x)`$.

For a general formula $`h^{-1}\sum_ka_kf(x+kh)`$, Taylor expansion begins with $`f(x)\sum_ka_k/h+f'(x)\sum_kka_k+O(h)`$. To converge to the derivative for every smooth function, the necessary conditions are $`\sum_ka_k=0`$ and $`\sum_kka_k=1`$. They are also sufficient for first-order consistency with fixed finite weights and enough smoothness. The first removes constant functions; the second differentiates a linear function correctly.

#### Problem NM-25: a fourth-order centered stencil

Source: [FNC, Exercise 5.5.6(a,b)](https://fncbook.com/python/fd-converge/).

Problem. Derive the fourth-order centered five-point first-derivative weights and its leading truncation term.

Solution. On offsets $`-2,-1,0,1,2`$, require $`\sum k^ja_k=0`$ for $`j=0,2,3,4`$ and one for $`j=1`$. These are five linear equations for five weights. The antisymmetric choice $`a_{-k}=-a_k`$, $`a_0=0`$ automatically satisfies the even-power conditions. The remaining equations are $`2a_1+4a_2=1`$ and $`2a_1+16a_2=0`$, yielding $`a_1=2/3`$, $`a_2=-1/12`$. Hence

```math
f'(x)\approx
\frac{f(x-2h)-8f(x-h)+8f(x+h)-f(x+2h)}{12h}.
```

The fifth moment of the weights is $`-4`$, so its leading error is $`-h^4f^{(5)}(x)/30`$ when the expansion is justified.

#### Problem NM-26: endpoint-corrected trapezoidal integration

Source: [FNC, Exercise 5.6.2](https://fncbook.com/python/integration/).

Problem. Cancel the leading composite trapezoidal error with endpoint derivatives, then approximate those derivatives to obtain an endpoint-corrected Gregory rule.

Solution. Subtract $`h^2[f'(b)-f'(a)]/12`$ from $`T_h`$ to remove its leading error. Approximate the endpoint derivatives by second-order one-sided formulas:

```math
f'(a)\approx\frac{-3f_0+4f_1-f_2}{2h},\qquad
f'(b)\approx\frac{3f_n-4f_{n-1}+f_{n-2}}{2h}.
```

The resulting Gregory formula is

```math
G_h=T_h-\frac h{24}
[3(f_n+f_0)-4(f_{n-1}+f_1)+(f_{n-2}+f_2)].
```

Each endpoint derivative error is $`O(h^2)`$, and multiplying it by $`h^2`$ leaves an $`O(h^4)`$ effect. Combined with the trapezoid expansion's next term, this gives fourth-order accuracy under the required smoothness.

#### Problem NM-27: deriving Simpson weights

Source: [FNC, Exercise 5.6.4(a–d)](https://fncbook.com/python/integration/).

Problem. Integrate a quadratic interpolant over a pair of panels to derive Simpson's local and composite weights.

Solution. The quadratic interpolating values $`\alpha,\beta,\gamma`$ at $`-h,0,h`$ is $`p(x)=\beta+(\gamma-\alpha)x/(2h)+(\alpha-2\beta+\gamma)x^2/(2h^2)`$. Substitution verifies its three values. Integrating over $`[-h,h]`$ removes the odd term and gives

```math
\int_{-h}^hp(x)dx
=2h\beta+\frac h3(\alpha-2\beta+\gamma)
=\frac h3(\alpha+4\beta+\gamma).
```

Translate each consecutive pair of mesh intervals to this centered coordinate and sum. Shared endpoints receive weight two, odd-index interior midpoints receive weight four, and outer endpoints receive weight one. Thus for even $`n`$,

```math
S_h=\frac h3[f_0+4f_1+2f_2+\cdots+4f_{n-1}+f_n].
```

### 5.5 Additional solved problems

#### Problem NM-28: Richardson extrapolation

Constructed problem. If $`A(h)=L+ch^2+dh^4+O(h^6)`$, eliminate the $`h^2`$ term.

Solution. $`A(h/2)=L+ch^2/4+dh^4/16+O(h^6)`$. Therefore $`[4A(h/2)-A(h)]/3=L-dh^4/4+O(h^6)`$, a fourth-order approximation. The improvement depends on the stated error expansion; it need not occur for a nonsmooth problem whose leading error is not a fixed multiple of $`h^2`$.

#### Problem NM-29: antithetic integration

Constructed problem. For $`\int_0^1u^2du`$, use paired average $`Y=[U^2+(1-U)^2]/2`$.

Solution. $`Y=U^2-U+1/2`$ has mean $`1/3`$ and variance $`1/180`$. Averaging $`N`$ independent pairs gives variance $`1/(180N)`$. At the same $`2N`$ function-evaluation count, ordinary Monte Carlo has variance $`(4/45)/(2N)=2/(45N)`$, eight times larger. Count computational effort consistently when comparing variance-reduction procedures.

#### Problem NM-30: an importance-sampling condition

Constructed problem. For $`I=\int f(x)p(x)dx`$, derive an estimator using draws from density $`q`$.

Solution. Where $`q\gt 0`$ whenever $`f p\ne0`$, write $`I=E_q[f(X)p(X)/q(X)]`$. The iid sample average of these weights is unbiased when integrable. Its variance is $`N^{-1}[\int f(x)^2p(x)^2/q(x)dx-I^2]`$ when finite. A proposal with very small $`q`$ in an important region can produce enormous or infinite variance even if it visits that region occasionally.

#### Problem NM-31: quadrature with a nonsmooth endpoint

Constructed problem. Explain why a fourth-order rule need not achieve fourth order for $`\sqrt x`$ on $`[0,1]`$.

Solution. The derivatives required by the usual Simpson error bound are unbounded at zero. A numerical result can still converge, but the proof giving $`O(h^4)`$ does not apply. The substitution $`x=t^2`$ changes $`\int_0^1\sqrt x\,dx`$ into $`\int_0^12t^2dt`$, a smooth polynomial integral. A variable transformation can improve the mathematical approximation problem before changing the quadrature rule.

## 6. Numerical ordinary differential equations

### 6.1 Theory

For $`y'=f(t,y)`$, $`y(0)=y_0`$, explicit Euler uses $`y_{n+1}=y_n+hf(t_n,y_n)`$. Implicit Euler uses $`y_{n+1}=y_n+hf(t_{n+1},y_{n+1})`$, which requires solving an equation at every step.

The local one-step defect compares an exact solution advanced one step with the numerical formula started from its exact current value. For Euler it is $`O(h^2)`$. Some texts divide this defect by $`h`$ and call that the local truncation error, giving $`O(h)`$. The global error at a fixed time after many steps is first-order under the usual Lipschitz and smoothness assumptions.

Runge-Kutta methods combine several derivative evaluations. The explicit midpoint and Heun methods are second-order; merely using two evaluations does not guarantee that order. Coefficients must satisfy the appropriate Taylor-matching conditions.

Absolute stability concerns the numerical solution of $`y'=\lambda y`$ for $`\mathop{\mathrm{Re}} \lambda\lt 0`$. A method can be high-order on smooth solutions and still unstable for a step too large relative to rapidly decaying components. Stiffness makes explicit stability limits much smaller than the time scale on which the main solution changes.

### 6.2 Derivations and proofs

#### Euler's global error

Let $`f`$ be Lipschitz in $`y`$ with constant $`L`$, and let the exact one-step defect be at most $`Ch^2`$. If $`e_n=y(t_n)-y_n`$, subtraction gives

```math
\vert e_{n+1}\vert \le(1+Lh)\vert e_n\vert +Ch^2.
```

With exact initial data, iterate the bound:

```math
\vert e_n\vert \le Ch^2\sum_{j=0}^{n-1}(1+Lh)^j
\le\frac{Ch}{L}(e^{LT}-1),\quad nh\le T,
```

for $`L\gt 0`$. For $`L=0`$, the sum gives $`CTh`$. A second-order local defect has accumulated over order $`1/h`$ steps, yielding first-order global error.

#### Stability of Euler methods

For $`y'=\lambda y`$, explicit Euler gives $`y_n=(1+h\lambda)^ny_0`$. Stability requires $`\vert 1+h\lambda\vert \le1`$; strict inequality gives decay of that mode. For real $`\lambda\lt 0`$, this means $`0\lt h\le2/\vert \lambda\vert`$, with a nondecaying oscillation at the upper endpoint.

Implicit Euler gives amplification $`1/(1-h\lambda)`$. If $`\mathop{\mathrm{Re}} \lambda\le0`$, then $`\vert 1-h\lambda\vert \ge1`$, so it is stable throughout the left half-plane. This is A-stability. It also damps very negative real modes to zero as $`h\vert \lambda\vert \to\infty`$.

#### A second-order step

Heun uses predictor $`\widetilde y=y_n+hf(t_n,y_n)`$ and correction $`y_{n+1}=y_n+\tfrac h2[f(t_n,y_n)+f(t_n+h,\widetilde y)]`$. Expand the second evaluation:

```math
f(t+h,y+hf)=f+h(f_t+f_yf)+O(h^2).
```

Thus the update is $`y+hf+\tfrac12h^2(f_t+f_yf)+O(h^3)`$, matching the exact Taylor expansion through second order. The local defect is third-order and the global order is two under stability and smoothness assumptions.

#### Local defect accumulates into global error

Suppose a one-step method has local error at most $`Ch^{p+1}`$ when started from the exact solution and a perturbation amplification bound $`1+Lh`$. With global error $`e_n`$,

```math
\vert e_{n+1}\vert \le(1+Lh)\vert e_n\vert +Ch^{p+1}.
```

Starting from exact initial data and summing the geometric series gives, for $`L\gt 0`$,

```math
\vert e_n\vert \le
\frac{(1+Lh)^n-1}{L}Ch^p
\le\frac{e^{LT}-1}{L}Ch^p,\qquad nh\le T.
```

For $`L=0`$, the bound becomes $`CT h^p`$. One factor of $`h`$ is lost because about $`T/h`$ local defects accumulate. This explains why an Euler local defect of order $`h^2`$ produces first-order global error.

The bound also separates consistency from stability. Small local defects can still become large if the numerical propagation amplifies them strongly. A stable method's amplification estimate must be suitable for the equation and the chosen step.

#### Why Heun has second-order accuracy

Let $`k_1=f(t,y)`$ and $`k_2=f(t+h,y+hk_1)`$. A smooth expansion gives $`k_2=f+h(f_t+f_yf)+O(h^2)`$. The update

```math
y_{\mathrm{new}}=y+\frac h2(k_1+k_2)
=y+hf+\frac{h^2}{2}(f_t+f_yf)+O(h^3)
```

matches the exact solution's Taylor expansion through order two because $`y''=f_t+f_yf`$. Thus the local defect is third order and, under a stability estimate, the global method is second order. Counting two evaluations alone would not prove this matching.

#### Example NM-E22: stiffness despite a smooth exact path

Consider $`y'=-1000(y-\cos t)-\sin t`$ with $`y(0)=1`$. The exact solution is $`y(t)=\cos t`$, as direct differentiation verifies.

A perturbation $`e`$ from this path obeys $`e'=-1000e`$. Explicit Euler multiplies that perturbation by $`1-1000h`$, so stability requires $`h\le0.002`$ on the negative real test axis. With $`h=0.01`$ the multiplier is minus nine.

The desired path changes on an order-one time scale, but tiny truncation or rounding perturbations can excite a much faster mode. This is a basic stiffness mechanism. Backward Euler damps that mode for every positive step, although its step still must resolve the cosine path accurately.

A stable implicit method removes this particular stability restriction; it does not make arbitrary large steps accurate. The linear solve used at each implicit step also has its own error, which must remain controlled relative to the time discretization.

### 6.3 Worked examples

#### Example NM-E12: stable and unstable decay

For $`y'=-10y`$ and $`y_0=1`$, explicit Euler with $`h=0.1`$ has amplification zero. With $`h=0.25`$ it has amplification $`-1.5`$ and grows in magnitude while the exact solution decays. Implicit Euler at $`h=0.25`$ has amplification $`1/3.5`$, which decays but may have significant discretization error.

#### Example NM-E13: a midpoint step

For $`y'=y`$, $`y_0=1`$, and $`h=0.2`$, explicit midpoint evaluates the slope at predicted midpoint $`1+0.1=1.1`$, giving $`y_1=1+0.2(1.1)=1.22`$. Euler gives $`1.2`$, while $`e^{0.2}\approx1.22140`$. The second-order method is more accurate for this step, with no claim that it is unconditionally stable.

### 6.4 Selected exercises

#### Problem NM-32: four two-step Euler calculations

Source: [FNC, Exercise 6.2.1(a–d)](https://fncbook.com/python/euler/).

Problem. Take two explicit Euler steps for each initial-value problem: (a) $`u'=-2tu`$, $`u(0)=2`$, $`h=0.1`$; (b) $`u'=u+t`$, $`u(0)=2`$, $`h=0.2`$; (c) $`u'=(1-u)/t`$, $`u(1)=6`$, $`h=0.25`$; (d) $`u'=2u(1-u)`$, $`u(0)=1/2`$, $`h=0.25`$. Compare each terminal approximation with its exact solution.

Solution. In each case update from the old time and old state.

(a) For $`u'=-2tu`$, $`u(0)=2`$, $`h=0.1`$, the first step stays at two and the second is $`2+0.1[-2(0.1)(2)]=1.96`$. The exact value at $`0.2`$ is $`2e^{-0.04}\approx1.921579`$, so numerical-minus-exact error is about $`0.038421`$.

(b) For $`u'=u+t`$, $`u(0)=2`$, $`h=0.2`$, the steps are $`2.4`$ and $`2.4+0.2(2.4+0.2)=2.92`$. Exact $`u(0.4)=-1.4+3e^{0.4}\approx3.075474`$, so the error is about $`-0.155474`$.

(c) For $`u'=(1-u)/t`$, $`u(1)=6`$, $`h=0.25`$, the first value is $`4.75`$. The next is $`4.75+0.25(1-4.75)/1.25=4`$. Exact $`u(1.5)=1+5/1.5=13/3`$, so the error is $`-1/3`$.

(d) For $`u'=2u(1-u)`$, $`u(0)=1/2`$, $`h=0.25`$, the first value is $`0.625`$. The second is $`0.625+0.25[2(0.625)(0.375)]=0.7421875`$. Exact $`u(0.5)=1/(1+e^{-1})\approx0.731059`$, so the error is about $`0.011129`$.

#### Problem NM-33: two evaluations without second-order accuracy

Source: [FNC, Exercise 6.2.3(a,b)](https://fncbook.com/python/euler/).

Problem. For the method predicting $`v=y+hf(t,y)`$ and updating $`y^+=y+hf(t+h,v)`$, establish consistency and determine its local and global orders.

Solution. The source method predicts $`v=y+hf(t,y)`$, then uses $`y^+=y+hf(t+h,v)`$. Its one-step function is $`\phi(t,y,h)=f(t+h,y+hf(t,y))`$. As $`h\to0`$, $`\phi\to f(t,y)`$, establishing consistency.

Taylor expansion gives $`y^+=y+hf+h^2(f_t+f_yf)+O(h^3)`$. The exact coefficient of the second-order term is one half, not one. Thus the local defect is generally $`O(h^2)`$ and the global method is only first-order. The average of the two slopes in Heun is what supplies the missing factor.

#### Problem NM-34: explicit Euler's unstable test equation

Source: [FNC, Exercise 6.2.4(a,b)](https://fncbook.com/python/euler/).

Problem. Solve explicit Euler's recurrence on $`u'=ku`$ and give a negative $`k`$ and fixed step for which the exact solution decays but the numerical values diverge.

Solution. The recurrence $`u_{i+1}=(1+hk)u_i`$, $`u_0=1`$, gives $`u_i=(1+hk)^i`$ at time $`ih`$. Choose $`k=-1`$ and $`h=3`$. The exact solution $`e^{-t}`$ is bounded and decays, but the numerical solution is $`(-2)^i`$, whose magnitude diverges. The fixed-step long-time instability is distinct from convergence as $`h\to0`$ over a fixed finite interval.

### 6.5 Additional solved problems

#### Problem NM-35: accumulated roundoff in Euler

Source: [FNC, Exercise 6.2.6](https://fncbook.com/python/euler/).

Problem. Include a uniformly bounded roundoff perturbation in each step of an Euler-type method and derive the resulting global truncation-plus-roundoff bound.

Solution. Let the per-step defect be bounded by $`Ch^{p+1}+u`$, with roundoff perturbation magnitude at most $`u`$. The error recurrence becomes $`\vert e_{n+1}\vert \le(1+Lh)\vert e_n\vert +Ch^{p+1}+u`$. With exact initial value,

```math
\vert e_n\vert \le\frac{e^{LT}-1}{L}(Ch^p+u/h)
```

for $`L\gt 0`$, and the limiting $`T(Ch^p+u/h)`$ form for $`L=0`$. This worst-case bound allows all roundoff errors to align. Smaller $`h`$ reduces the truncation term but enlarges the accumulated-roundoff term.

#### Problem NM-36: error versus stability

Constructed problem. Is an implicit Euler solution accurate merely because it is stable at a large step?

Solution. For $`y'=-y`$ with $`h=10`$, implicit Euler gives one-step multiplier $`1/11\approx0.0909`$, while the exact multiplier is $`e^{-10}\approx0.0000454`$. Both decay and the numerical method is stable, but the relative discrepancy in this mode is large. Stability controls amplification; accuracy also requires sufficiently small discretization error for the intended quantity.

#### Problem NM-37: a stiff pair

Constructed problem. Consider $`\dot x=-x`$, $`\dot y=-1000y`$.

Solution. Explicit Euler must use $`h\lt 0.002`$ for strict stability of the fast component, even if only the slowly varying $`x`$ is of primary interest after the rapid transient. An implicit method can remain stable at a larger step, but should still resolve any transient whose accuracy matters. Stiffness describes the mismatch between fast numerical stability restrictions and slower solution features.

#### Problem NM-38: order estimation from step halving

Constructed problem. If global errors at $`h,h/2,h/4`$ are approximately $`0.04,0.01,0.0025`$, estimate the order.

Solution. The error ratio is four under each halving. For a leading form $`Ch^p`$, the ratio is $`2^p`$, so $`p\approx2`$. Without an exact reference, differences of successive approximations can estimate the ratio if they are in the asymptotic regime and share the same leading coefficient. A short table alone does not prove the order for every problem.

[Contents](README.md) · [Part 2](07-numerical-methods-part-2.md)
