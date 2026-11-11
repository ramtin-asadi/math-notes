# Numerical Methods

A mathematical formula and a reliable numerical procedure are different objects. These notes track the approximation error, the sensitivity of the underlying problem, and the behavior of the numerical method separately.

Sources:

- Tobin Driscoll and Richard Braun, [Fundamentals of Numerical Computation, online Python edition](https://fncbook.com/python). Exercise numbers below refer to the online sections linked individually.
- Laurent Demanet, [MIT 18.330, Introduction to Numerical Analysis, Spring 2012](https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/), notes and assignments.
- Boyd and Vandenberghe, [Convex Optimization](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf), Chapter 9, for descent methods.

Small computations are shown by hand. Selected source exercises identify theoretical subparts when a larger source problem also requests programming. No software project is needed to follow the notes.

Contents:

1. [Error, floating point, and conditioning](#1-error-floating-point-and-conditioning)
2. [Numerical linear algebra](#2-numerical-linear-algebra)
3. [Root finding](#3-root-finding)
4. [Numerical optimization](#4-numerical-optimization)
5. [Differentiation and integration](#5-differentiation-and-integration)
6. [Numerical ordinary differential equations](#6-numerical-ordinary-differential-equations)
7. [Numerical partial differential equations](#7-numerical-partial-differential-equations)
8. [Further worked problems](#8-further-worked-problems)

## 1. Error, floating point, and conditioning

### 1.1 Theory

If $x$ is exact and $\widehat x$ is an approximation, absolute error is $|\widehat x-x|$ and relative error is that quantity divided by $|x|$, when $x\ne0$. A relative error at a zero answer is undefined; choose an absolute scale instead.

A floating-point number stores a sign, a finite significand, and an exponent. Binary64 has 53 binary digits of significand precision, including the leading bit for normal numbers. The spacing immediately above one is $2^{-52}$, often called machine epsilon. Unit roundoff for rounding to nearest is $u=2^{-53}$. Naming the convention avoids a factor-of-two ambiguity.

For an ordinary operation whose exact result remains in the normal representable range, the usual model is $\operatorname{fl}(a\circ b)=(a\circ b)(1+\delta)$ with $|\delta|\le u$. Overflow, underflow, subnormal results, division by zero, and invalid operations require separate treatment.

Truncation error comes from replacing an exact mathematical operation with a finite approximation, such as a finite Taylor sum or a difference quotient. Roundoff error comes from finite representation and arithmetic. Refining an approximation can decrease truncation error while increasing the influence of roundoff.

Conditioning belongs to the mathematical input-output problem. Stability belongs to an algorithm used to solve it. Forward error compares the computed output with the true output. Backward error asks for a nearby input problem that the computed output solves exactly. A backward-stable method can have a large forward error on an ill-conditioned problem.

### 1.2 Derivations and proofs

#### Scalar relative condition number

For differentiable $f$ and nonzero $x,f(x)$, a small relative input perturbation $x\varepsilon$ gives $f(x+x\varepsilon)-f(x)=f'(x)x\varepsilon+o(\varepsilon)$. Hence

$$
\kappa_f(x)=\left|\frac{xf'(x)}{f(x)}\right|.
$$

The absolute condition number is $|f'(x)|$. A large relative number can arise because the output is close to zero even when the absolute sensitivity is moderate.

#### Cancellation and reformulation

For small $x$, evaluating $\sqrt{1+x}-1$ subtracts two numbers near one to produce a much smaller answer. Rationalize:

$$
\sqrt{1+x}-1=\frac{x}{\sqrt{1+x}+1}.
$$

The exact problem is well-conditioned near zero in relative terms, since the output behaves like $x/2$. The direct algorithm can first round $1+x$ to one and return zero. The rationalized form retains the input in the numerator and divides by a quantity near two.

Subtraction itself may be exactly representable; the danger is amplification of errors already present in two nearly equal operands. A cancellation diagnosis should explain where those operand errors came from.

#### Unit roundoff and accumulated summation error

In binary64 arithmetic with rounding to nearest, the gap from one to the next larger representable number is $2^{-52}$. The unit roundoff used in standard relative-error bounds is half that gap, $u=2^{-53}$. For an elementary operation whose exact result remains in the normal range, a common model is $\operatorname{fl}(a\circ b)=(a\circ b)(1+\delta)$ with $|\delta|\le u$. Overflow, underflow, and exceptional inputs require separate treatment.

Repeated rounding errors do not simply add as signed relative errors. A useful product bound is $\prod_{i=1}^n(1+\delta_i)=1+\theta_n$ with $|\theta_n|\le\gamma_n=nu/(1-nu)$ when $nu<1$. Expanding the products or bounding them by a geometric series gives this estimate.

Sequential summation then satisfies

$$
\left|\operatorname{fl}\left(\sum_{i=1}^n x_i\right)-\sum_{i=1}^n x_i\right|
\le\gamma_{n-1}\sum_{i=1}^n|x_i|.
$$

The relative error bound divides by the magnitude of the exact sum. Its condition factor is $\sum|x_i|/|\sum x_i|$ when the sum is nonzero. Cancellation can make this factor large. If the exact sum is zero, relative error against that sum is undefined and an absolute scale must be used.

Pairwise summation arranges additions in a balanced tree. Each input passes through only about $\lceil\log_2n\rceil$ rounding stages, producing a corresponding smaller depth-based bound. This changes the algorithmic error without changing the mathematical sum or its conditioning.

#### Example NM-E17: stable quadratic roots

Solve $x^2+10^8x+1=0$ accurately.

Solution. The direct formula computes one root as $(-10^8+\sqrt{10^{16}-4})/2$. Its numerator subtracts nearly equal numbers, and the small root can lose many digits. Instead set

$$
q=-\frac12\left(b+\operatorname{sign}(b)\sqrt{b^2-4ac}\right),
\qquad x_1=q/a,\qquad x_2=c/q.
$$

Here $q$ uses addition of same-sign magnitudes, and the second root uses the product relation $x_1x_2=c/a$. The roots are approximately $-10^8$ and $-10^{-8}$.

The formula assumes a nonzero $q$ and a suitable real discriminant for the real-root problem. Degenerate cases such as zero coefficients must be handled from their actual equation. The correction is an algebraic reformulation, not a change to the polynomial being solved.

### 1.3 Worked examples

#### Example NM-E01: a lost unit

In binary64 round-to-nearest, adding one to $2^{53}$ lies exactly between adjacent representable numbers. Ties-to-even rounds it back to $2^{53}$, so the computed expression $(2^{53}+1)-2^{53}$ can equal zero rather than one. The first addition loses the unit; the final subtraction cannot recover information no longer stored.

#### Example NM-E02: relative sensitivity near a zero

For $f(x)=\log x$, $\kappa_f(x)=1/|\log x|$. At $x$ near one the output is near zero, making relative sensitivity large. The absolute condition number $1/x$ is near one. An absolute error can remain small while its fraction of a tiny exact logarithm is large.

### 1.4 Selected exercises

#### Problem NM-01: four elementary condition numbers

Source: [FNC, Exercise 1.2.1(a–d)](https://fncbook.com/python/conditioning/).

Problem. Find the relative condition numbers of a power, logarithm, cosine, and exponential, specifying where the relative definitions apply.

Solution. Apply $|xf'/f|$ on domains where it is defined. For $f=x^p$ on $x>0$, it is $|p|$. For $\log x$, with $x>0$ and $x\ne1$, it is $1/|\log x|$. For $\cos x$ away from its zeros, it is $|x\tan x|$. For $e^x$ at nonzero input, it is $|x|$. At zero input the relative-input definition needs a limiting convention; the corresponding absolute condition numbers remain well-defined. These are local sensitivities, not finite-perturbation error guarantees of arbitrary size.

#### Problem NM-02: composed functions

Source: [FNC, Exercise 1.2.2(a–c)](https://fncbook.com/python/conditioning/).

Problem. Find the relative condition numbers of $\sqrt{x+5}$, $\cos(2\pi x)$, and $e^{-x^2}$ and identify excluded inputs or outputs.

Solution. Direct differentiation gives, respectively,

$$
\kappa_{\sqrt{x+5}}=\left|\frac{x}{2(x+5)}\right|,\quad x>-5,
$$

$$
\kappa_{\cos(2\pi x)}=|2\pi x\tan(2\pi x)|,\qquad
\kappa_{e^{-x^2}}=2x^2.
$$

The cosine expression excludes zero outputs; the square-root expression excludes its boundary. Applying the chain rule for condition numbers gives the same expressions wherever the intermediate relative quantities are defined. At an intermediate zero, compute the direct derivative or take a justified limit instead of multiplying undefined factors.

#### Problem NM-03: the chain rule for conditioning

Source: [FNC, Exercise 1.2.4](https://fncbook.com/python/conditioning/).

Problem. Derive the relative condition-number chain rule for a composition at points where all required ratios are defined.

Solution. For $h=f\circ g$, with all denominators nonzero,

$$
\kappa_h(x)=\left|\frac{xf'(g(x))g'(x)}{f(g(x))}\right|
=\left|\frac{g(x)f'(g(x))}{f(g(x))}\right|
\left|\frac{xg'(x)}{g(x)}\right|.
$$

This is $\kappa_f(g(x))\kappa_g(x)$. Differentiability and nonzero intermediate input and output values justify the division. The result describes the composite mathematical map, not every numerical implementation of the composition.

#### Problem NM-04: a cancellation-heavy exponential series

Source: [MIT 18.330, Homework 1, Problem 4](https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/d8531883a96517b39f1b6640fa983d65_MIT18_330S12_hw1.pdf), explanation of the finite-precision failure.

Problem. Explain why direct finite-precision summation of the Taylor series for $e^{-25}$ can be inaccurate even after including many terms.

Solution. The exact series $\sum_{n\ge0}(-25)^n/n!$ equals $e^{-25}\approx1.39\times10^{-11}$. Its terms initially grow in magnitude, reaching sizes of order $e^{25}/\sqrt{2\pi25}$ near $n=25$, before shrinking. Positive and negative terms of order billions must cancel to leave an answer of order $10^{-11}$. Tiny relative errors in the large terms can overwhelm the result.

Adding more terms reduces truncation error but cannot recover digits lost through that cancellation. Computing $e^{25}$ with positive terms and then taking its reciprocal is a more stable conceptual route; a practical exponential implementation also uses range reduction and carefully controlled approximations.

### 1.5 Additional solved problems

#### Problem NM-05: a stable quadratic root formula

Constructed problem. For $x^2-bx+1=0$ with large $b>2$, compute the small root without subtracting nearly equal numbers.

Solution. The direct small root $(b-\sqrt{b^2-4})/2$ suffers cancellation. The large root is $r_+=(b+\sqrt{b^2-4})/2$, and the product of the roots is one, so $r_-=1/r_+$. The latter avoids the dangerous subtraction. Separate overflow protection may still be needed when forming $b^2$ for extremely large $b$; an algebraically improved formula is not automatically safe over the entire machine range.

#### Problem NM-06: perturbing a repeated root

Constructed problem. Compare roots of $(x-1)^2=0$ and $(x-1)^2-\varepsilon=0$ for $\varepsilon>0$.

Solution. The perturbed roots are $1\pm\sqrt\varepsilon$. A coefficient perturbation of order $\varepsilon$ produces root changes of order $\sqrt\varepsilon$, whose ratio to $\varepsilon$ diverges. The repeated-root problem is ill-conditioned. A stable algorithm cannot remove this mathematical sensitivity.

## 2. Numerical linear algebra

### 2.1 Theory

To solve $Ax=b$, factor the matrix and solve triangular systems rather than explicitly forming an inverse. LU elimination writes $PA=LU$ with a permutation $P$ for pivoting. Forward substitution solves $Ly=Pb$ and backward substitution solves $Ux=y$. A zero or tiny pivot can require a row swap; nonsingularity of the full matrix does not guarantee every unpivoted leading pivot is safe.

For symmetric PD $A$, Cholesky writes $A=LL^\top$ with positive diagonal entries. Positive definiteness guarantees the square roots used in the factorization remain positive in exact arithmetic. A semidefinite matrix can require pivoting or a rank-revealing factorization.

For full-column-rank least squares, thin QR writes $A=QR$ with $Q^\top Q=I$ and square invertible upper-triangular $R$. The minimizer solves $R\widehat x=Q^\top b$. Orthogonal transformations preserve Euclidean norms and avoid explicitly squaring the condition number through $A^\top A$.

Householder reflection $H=I-2vv^\top/(v^\top v)$ is orthogonal and can zero a column's lower entries in one transformation. Successive reflections produce QR. Modified Gram-Schmidt is another construction, but strong near-dependence can require reorthogonalization.

### 2.2 Derivations and proofs

#### QR least squares without normal equations

Complete the thin $Q$ to a square orthogonal matrix $(Q,Q_\perp)$. Then

$$
\|Ax-b\|^2=\|Rx-Q^\top b\|^2+\|Q_\perp^\top b\|^2.
$$

The second term is fixed. The first becomes zero at the unique triangular-system solution. This is the same geometric projection as the normal equations, expressed in numerically safer coordinates.

#### Residual and forward error

If $A$ is invertible, define $r=b-A\widehat x$. Then $A(x-\widehat x)=r$, so $\|x-\widehat x\|\le\|A^{-1}\|\|r\|$. Since $\|b\|\le\|A\|\|x\|$,

$$
\frac{\|x-\widehat x\|}{\|x\|}
\le\kappa(A)\frac{\|r\|}{\|b\|}.
$$

The residual also means $\widehat x$ exactly solves $A\widehat x=b-r$, identifying a backward perturbation in the right-hand side. A small relative residual alone does not guarantee a small forward error.

#### Why normal equations square a condition number

For full-column-rank $A$, the singular values of $A^\top A$ are the squares of those of $A$. Hence $\kappa_2(A^\top A)=\kappa_2(A)^2$. The algebraic normal equations remain correct, but forming the Gram matrix can discard small distinctions between nearly dependent columns.

A thin QR factorization $A=QR$ gives $\|Ax-b\|^2=\|Rx-Q^\top b\|^2+\|(I-QQ^\top)b\|^2$. Solving the triangular system uses the geometry of $A$ directly. Householder QR is backward stable under the standard arithmetic model, whereas explicitly forming normal equations can magnify the effect of rounding through the squared condition number.

Cholesky is appropriate for a supplied SPD matrix, but applying it to a rounded Gram matrix does not recover information already lost in forming that matrix. The Cholesky recurrence is

$$
L_{jj}=\sqrt{a_{jj}-\sum_{k<j}L_{jk}^2},\qquad
L_{ij}=\frac{a_{ij}-\sum_{k<j}L_{ik}L_{jk}}{L_{jj}}.
$$

Positive Schur complements of an SPD matrix make the exact radicands positive. A computed nonpositive radicand can reflect loss of definiteness through rounding or a matrix that was never SPD.

#### Example NM-E18: an independent pair hidden by its Gram matrix

Take

$$
A=\begin{pmatrix}1&1\\\epsilon&0\\0&\epsilon\end{pmatrix},
\qquad b=\begin{pmatrix}0\\\epsilon\\-\epsilon\end{pmatrix},
\qquad \epsilon=10^{-8}.
$$

Solution. The exact solution is $x=(1,-1)^\top$. The columns are independent. Their Gram matrix is

$$
A^\top A=
\begin{pmatrix}1+\epsilon^2&1\\1&1+\epsilon^2\end{pmatrix},
$$

with eigenvalues $2+\epsilon^2$ and $\epsilon^2$. Thus $\kappa_2(A)\approx\sqrt2\,10^8$, while the Gram condition number is about $2\cdot10^{16}$.

In binary64, adding $10^{-16}$ to one rounds to one. Directly formed diagonal entries can therefore equal the off-diagonal entries, making the rounded Gram matrix singular. The original rectangular matrix still stores the small entries explicitly.

QR avoids this specific loss, though it cannot remove the original condition number. A backward error of machine scale can still create a forward error amplified by about $10^8$. Stability of the solver and conditioning of the supplied data remain separate issues.

### 2.3 Worked examples

#### Example NM-E03: Cholesky by hand

For $A=\left(\begin{smallmatrix}4&2\\2&3\end{smallmatrix}\right)$, write $L=\left(\begin{smallmatrix}l_{11}&0\\l_{21}&l_{22}\end{smallmatrix}\right)$. Matching entries gives $l_{11}=2$, $l_{21}=1$, and $l_{22}=\sqrt2$. All roots are positive. To solve $Ax=(6,5)^\top$, forward substitution gives $y=(3,\sqrt2)$ and backward substitution gives $x=(1,1)$.

#### Example NM-E04: why pivoting is useful

For $A=\left(\begin{smallmatrix}\varepsilon&1\\1&1\end{smallmatrix}\right)$ with small nonzero $\varepsilon$, unpivoted elimination uses multiplier $1/\varepsilon$, producing large intermediate entries. Swapping rows first uses multiplier $\varepsilon$, keeping them moderate. The original matrix remains well-conditioned as $\varepsilon\to0$; the large intermediates are a feature of the unpivoted algorithm.

### 2.4 Selected exercises

#### Problem NM-07: orthogonal-matrix properties

Source: [FNC, Exercise 3.3.2](https://fncbook.com/python/qr/).

Problem. For square orthogonal matrices, prove the inverse, transpose, product, norm-preservation, and condition-number identities and show invariance of $\|AQ\|_2$.

Solution. For square $Q$ with $Q^\top Q=I$, invertibility gives $Q^{-1}=Q^\top$, and hence $QQ^\top=I$, proving $Q^\top$ orthogonal. For any $x$, $\|Qx\|_2^2=x^\top Q^\top Qx=\|x\|_2^2$, so $\|Q\|_2=\|Q^{-1}\|_2=1$ and $\kappa_2(Q)=1$. The substitution $y=Qx$ maps the unit sphere onto itself, so $\|AQ\|_2=\sup_{\|y\|=1}\|Ay\|_2=\|A\|_2$. Finally $(QU)^\top QU=U^\top U=I$ for orthogonal $U$, so their product is orthogonal.

#### Problem NM-08: invertibility of the thin triangular factor

Source: [FNC, Exercise 3.3.4](https://fncbook.com/python/qr/).

Problem. For a thin QR factorization of a full-column-rank matrix, prove that its square triangular factor is invertible.

Solution. Suppose full-column-rank $A=QR$ had singular square $R$. There would be $z\ne0$ with $Rz=0$, giving $Az=QRz=0$, contrary to column independence. Hence $R$ is invertible. Conversely, invertible $R$ and orthonormal $Q$ make $A$ injective, since $QRz=0$ implies $Rz=Q^\top QRz=0$.

#### Problem NM-09: a pseudoinverse from QR

Source: [FNC, Exercise 3.3.6](https://fncbook.com/python/qr/).

Problem. Derive the pseudoinverse of a full-column-rank matrix from its thin QR factorization and verify the Moore-Penrose identities.

Solution. For thin full-rank QR, every least-squares fit solves $Rx=Q^\top b$, so $x=R^{-1}Q^\top b$. This defines $A^+=R^{-1}Q^\top$. To check directly, $A^+A=I$ and $AA^+=QQ^\top$, a symmetric projector. The two remaining Moore-Penrose identities follow by multiplying these products by $A$ and $A^+$ respectively.

#### Problem NM-10: a triangular inverse with large norm

Source: [FNC, Exercise 2.8.6(a–c)](https://fncbook.com/python/condition-number/).

Problem. Let $A_n=I-2N$ where $N$ has ones on the first superdiagonal. Find its inverse and infinity-norm condition number for $n>1$.

Solution. Let $A_n=I-2N$, where $N$ has ones on its first superdiagonal and $N^n=0$. The finite geometric identity gives

$$
A_n^{-1}=I+2N+4N^2+\cdots+2^{n-1}N^{n-1}.
$$

Thus $(A_n^{-1})_{ij}=2^{j-i}$ for $j\ge i$ and zero below. In particular $A_2^{-1}=\left(\begin{smallmatrix}1&2\\0&1\end{smallmatrix}\right)$ and $A_3^{-1}=\left(\begin{smallmatrix}1&2&4\\0&1&2\\0&0&1\end{smallmatrix}\right)$. The largest inverse row sum is $2^n-1$, while $\|A_n\|_\infty=3$ for $n>1$. Hence $\kappa_\infty(A_n)=3(2^n-1)$: a triangular matrix can be exponentially ill-conditioned.

### 2.5 Additional solved problems

#### Problem NM-11: conditioning of normal equations

Constructed problem. If full-column-rank $A$ has singular values one and $10^{-4}$, compare its condition number with that of $A^\top A$.

Solution. The first is $10^4$. The Gram matrix has eigenvalues one and $10^{-8}$, so its condition number is $10^8$. Forming the Gram matrix can make finite-precision loss more severe. Solving by QR avoids that explicit squaring, although the original problem still contains a weak direction with singular value $10^{-4}$.

#### Problem NM-12: iterative refinement as an error equation

Constructed problem. Given an approximate $\widehat x$, explain one refinement step.

Solution. Compute $r=b-A\widehat x$, solve $A\delta=r$, and set $\widehat x_{\mathrm{new}}=\widehat x+\delta$. In exact arithmetic this reaches the exact solution in one step because $\delta=x-\widehat x$. In finite precision the correction solve and residual computation are approximate. Higher-precision residuals and a sufficiently accurate factorization can improve the result, but severe conditioning can limit the attainable accuracy.

## 3. Root finding

### 3.1 Theory

Bisection applies to a continuous real function with opposite signs at the endpoints. The intermediate value theorem guarantees a root. Halve the interval and retain a half with a sign change. It gives a reliable location bound without derivatives, although it does not exploit a particularly smooth root.

Fixed-point iteration rewrites a root equation as $x=g(x)$ and iterates $x_{k+1}=g(x_k)$. A closed invariant interval and a uniform derivative bound $|g'|\le q<1$ give a contraction and convergence. A derivative bound without invariance or completeness is insufficient for a global existence claim.

Newton uses the tangent line: $x_{k+1}=x_k-f(x_k)/f'(x_k)$. It is locally quadratic near a sufficiently smooth simple root, but a poor starting point, a small derivative, or a multiple root can change its behavior. Bracketing safeguards combine reliability with fast local steps.

### 3.2 Derivations and proofs

#### Bisection error

Starting with width $b-a$, after $n$ halvings the retained interval has width $(b-a)2^{-n}$. Its midpoint $m_n$ lies at distance at most half that width from every point of the interval, including its retained root:

$$
|m_n-r|\le\frac{b-a}{2^{n+1}}.
$$

This bound assumes exact sign decisions and arithmetic. In finite precision, an interval can become too narrow for its computed midpoint to differ from an endpoint.

#### Local Newton convergence

Assume $f(r)=0$, $|f'|\ge m>0$, and $|f''|\le M$ in a neighborhood of $r$. Taylor expansion from $x_k$ to $r$ gives

$$
0=f(x_k)+f'(x_k)(r-x_k)+\frac12f''(\xi_k)(r-x_k)^2.
$$

Rearrange after the Newton update:

$$
|x_{k+1}-r|\le\frac{M}{2m}|x_k-r|^2.
$$

Choose the initial error small enough that the bound remains inside the neighborhood and the factor times the initial error is below one. Induction then keeps the iterates there and proves convergence. The derivative lower bound rules out a multiple root and prevents division by zero.

#### A local Newton theorem with a neighborhood check

Let $r$ be a simple root, and suppose $f$ is twice continuously differentiable on a closed interval around it, with $|f'(x)|\ge m>0$ and $|f''(x)|\le M$. Taylor-expand $f(r)$ about the current iterate $x$:

$$
0=f(x)+f'(x)(r-x)+\frac12f''(\xi)(r-x)^2.
$$

For the Newton update $N(x)=x-f(x)/f'(x)$, this gives

$$
|N(x)-r|\le\frac M{2m}|x-r|^2.
$$

Choose a radius $\delta$ inside the regularity interval such that $C\delta<1$, where $C=M/(2m)$. If $|x_0-r|\le\delta$, the estimate keeps the next iterate in the same interval and reduces the error by at least factor $C\delta$. Induction proves convergence and the quadratic error estimate.

The neighborhood check closes the proof: it prevents an estimate requiring local derivative bounds from being applied after an iterate has left their domain. A formula for the asymptotic error constant alone does not supply that invariance.

#### Multiple roots reduce the order

If $f(x)=(x-r)^m g(x)$ with $g(r)\ne0$, then near the root $f/f'=(x-r)/m+O((x-r)^2)$. Ordinary Newton therefore has

$$
x_{n+1}-r=\frac{m-1}{m}(x_n-r)+O((x_n-r)^2).
$$

For multiplicity greater than one, the leading error is linear. If the multiplicity is known, the modified step $x-mf(x)/f'(x)$ cancels that leading term and can restore quadratic convergence under the corresponding smoothness conditions.

#### Example NM-E19: a triple root

For $f(x)=(x-1)^3$, Newton's method gives $x_{n+1}=(2x_n+1)/3$ whenever $x_n\ne1$. Thus $e_{n+1}=(2/3)e_n$ exactly. Starting at four gives $e_n=3(2/3)^n$.

The step size tends to zero, and the residual $e_n^3$ tends to zero faster than the location error. A residual tolerance can therefore look impressive even while the root error is substantially larger. The derivative also tends to zero, violating the simple-root lower bound used above.

Multiplying the Newton correction by three lands at one in a single step for this pure power. A general multiple-root function has additional higher-order terms, so it need not terminate in one step even with the modification.

### 3.3 Worked examples

#### Example NM-E05: Newton for a square root

For $f(x)=x^2-2$, the update is $x_{k+1}=(x_k+2/x_k)/2$. Starting from $3/2$ gives $17/12$, then $577/408$. The exact error identity is $x_{k+1}-\sqrt2=(x_k-\sqrt2)^2/(2x_k)$, displaying quadratic convergence and showing that positive iterates stay at or above the positive root after the first step.

#### Example NM-E06: a bisection count

For $f(x)=x^3-2$ on $[1,2]$, require midpoint error at most $10^{-3}$. The bound is $2^{-(n+1)}$. Since $2^{-10}<10^{-3}$, nine halvings suffice for this midpoint convention. Specifying whether $n$ counts midpoint evaluations or completed halvings avoids an indexing ambiguity.

### 3.4 Selected exercises

#### Problem NM-13: local fixed-point convergence

Source: [FNC, Exercise 4.2.1(a–c)](https://fncbook.com/python/fixed-point/).

Problem. Check local convergence of the fixed-point maps $1+x-x^2/9$ at three, $\pi+\sin x/4$ at $\pi$, and $x+1-\tan(x/4)$ at $\pi$.

Solution. For $g=1+x-x^2/9$ at $p=3$, $g(3)=3$ and $g'(3)=1/3$. For $g=\pi+\sin x/4$ at $p=\pi$, the derivative is $-1/4$. For $g=x+1-\tan(x/4)$ at $p=\pi$, it is $1-\sec^2(\pi/4)/4=1/2$. All three absolute derivatives are below one. Continuity of the derivatives gives a small invariant contraction neighborhood around each fixed point, so sufficiently close starts converge. The negative derivative in the second case makes small errors alternate signs.

#### Problem NM-14: root extraction by a poor fixed-point map

Source: [FNC, Exercise 4.2.5](https://fncbook.com/python/fixed-point/).

Problem. Examine whether the iteration $x_{k+1}=a/x_k^{m-1}$ locally computes the positive $m$th root of $a>0$ for integer $m>1$.

Solution. For $g(x)=a/x^{m-1}$ and positive root $p=a^{1/m}$, $g'(p)=-(m-1)$. For $m>2$ the fixed point is repelling, so nearby nonexact starts generally move away. At $m=2$ the derivative is $-1$, and the derivative test alone is inconclusive. But $g(g(x))=x$, producing a two-cycle unless $x=p$. Thus none of the integer values $m>1$ yields local convergence from all sufficiently close starts, apart from starting at the exact root itself.

#### Problem NM-15: choosing a constant derivative scale

Source: [FNC, Exercise 4.2.7](https://fncbook.com/python/fixed-point/).

Problem. For $g(x)=x-f(x)/c$ near a simple root $p$ with $f'(p)>0$, determine the constant scales $c$ giving a local contraction.

Solution. For $g(x)=x-f(x)/c$, a root $p$ is a fixed point and $g'(p)=1-f'(p)/c$. With $f'(p)>0$, the condition $|1-f'(p)/c|<1$ is equivalent to $c>f'(p)/2$. Derivative continuity then supplies a local contraction neighborhood. A negative $c$ fails this criterion. Choosing $c=f'(p)$ makes the derivative zero, but the true derivative at the unknown root is not usually available in advance.

#### Problem NM-16: division-free reciprocal Newton iteration

Source: [FNC, Exercise 4.3.5](https://fncbook.com/python/newton/).

Problem. Derive a division-free Newton iteration for the reciprocal of $b\ne0$ and prove its error recursion and convergence condition.

Solution. For $f(x)=1/x-b$, $f'=-1/x^2$. Newton simplifies to $x_{k+1}=2x_k-bx_k^2=x_k(2-bx_k)$, using only multiplication and subtraction. If $e_k=1-bx_k$, then $e_{k+1}=e_k^2$. Thus $|e_0|<1$ gives quadratic convergence to the reciprocal. The iteration formula can be evaluated without division, though a suitable initial scale is still required.

#### Problem NM-17: Newton can cycle

Source: [FNC, Exercise 4.3.6](https://fncbook.com/python/newton/).

Problem. Analyze Newton's method for $f(x)=\operatorname{sign}(x)\sqrt{|x|}$ from a nonzero starting value.

Solution. For $f(x)=\operatorname{sign}(x)\sqrt{|x|}$ and $x\ne0$, $f'(x)=1/(2\sqrt{|x|})$. Hence $f/f'=2x$ and Newton gives $x_{k+1}=-x_k$. Every nonzero start alternates forever between two points. The root is zero, where the function is not differentiable with a finite derivative. The smooth simple-root theorem does not apply.

### 3.5 Additional solved problems

#### Problem NM-18: small steps without convergence

Source: [FNC, Exercise 4.3.8](https://fncbook.com/python/newton/).

Problem. Give a sequence whose consecutive steps tend to zero but which fails to converge, and explain the implication for an iterative stopping rule.

Solution. Let $x_k=\sum_{j=1}^k1/j$. Then $x_k-x_{k-1}=1/k\to0$, but the harmonic sums diverge. For example, each block from $2^{m-1}+1$ through $2^m$ contributes at least $1/2$, so the total is unbounded. A small iterate change does not alone certify proximity to a limit or a root.

## 4. Numerical optimization

### 4.1 Theory

Gradient descent uses $x_{k+1}=x_k-\alpha\nabla f(x_k)$. A positive step size does not automatically give descent; it must be small enough relative to curvature or chosen by a line search. If the gradient is $L$-Lipschitz, Taylor's integral bound gives

$$
f(x-\alpha\nabla f(x))
\le f(x)-\alpha(1-\alpha L/2)\|\nabla f(x)\|^2.
$$

Thus $0<\alpha<2/L$ ensures a strict decrease whenever the gradient is nonzero. For general nonconvex functions, decreasing values need not imply finding a global optimum.

Projected gradient uses $x_{k+1}=\Pi_C(x_k-\alpha\nabla f(x_k))$ for a closed convex feasible set $C$. Its fixed-point condition is equivalent to $\nabla f(x^*)^\top(y-x^*)\ge0$ for every $y\in C$. This is the first-order optimality condition for a differentiable convex objective.

### 4.2 Derivations and proofs

#### Gradient descent on a PD quadratic

For $f(x)=\tfrac12x^\top Hx-b^\top x$ with $H$ PD, let $x^*=H^{-1}b$. The error satisfies $e_{k+1}=(I-\alpha H)e_k$. In an orthonormal eigenbasis, each error coordinate is multiplied by $1-\alpha\lambda_i$.

Convergence for all starts therefore holds exactly when $|1-\alpha\lambda_i|<1$ for every eigenvalue, or $0<\alpha<2/\lambda_{\max}$. The optimal fixed step for the worst Euclidean contraction factor is $\alpha^*=2/(\lambda_{\min}+\lambda_{\max})$, giving factor $(\kappa-1)/(\kappa+1)$. Ill-conditioning slows descent because no single step can move equally efficiently in both steep and flat directions.

#### Projection characterization

If $p=\Pi_C(z)$, differentiating the squared distance along the segment from $p$ to any $y\in C$ gives $(z-p)^\top(y-p)\le0$. Conversely this inequality makes $\|z-y\|^2=\|z-p\|^2+\|y-p\|^2-2(z-p)^\top(y-p)\ge\|z-p\|^2$. Apply it to $z=x^*-\alpha\nabla f(x^*)$ and $p=x^*$ to obtain the fixed-point optimality condition.

#### The descent lemma and a provable step

Suppose $\nabla f$ is Lipschitz with constant $L$ on the relevant convex domain. Integrating the gradient along the segment from $x$ to $x+h$ gives

$$
f(x+h)\le f(x)+\nabla f(x)^\top h+\frac L2\|h\|^2.
$$

For $h=-\eta\nabla f(x)$,

$$
f(x-\eta\nabla f(x))
\le f(x)-\eta(1-L\eta/2)\|\nabla f(x)\|^2.
$$

Thus every nonstationary step decreases the objective for $0<\eta<2/L$. With $\eta=1/L$, the decrease is at least $\|\nabla f(x)\|^2/(2L)$.

If $f$ is also $\mu$-strongly convex, minimizing its quadratic lower bound gives $f(x)-f^*\le\|\nabla f(x)\|^2/(2\mu)$. Combining yields $f(x_{k+1})-f^*\le(1-\mu/L)(f(x_k)-f^*)$. This is a global linear-rate bound under both assumptions.

#### Projection supplies a constrained fixed-point condition

For a nonempty closed convex set $C$, $p=\Pi_C(z)$ is characterized by $(z-p)^\top(y-p)\le0$ for all $y\in C$. Apply the characterization to two projections and add the resulting inequalities to obtain

$$
\|\Pi_C(z)-\Pi_C(w)\|^2
\le(\Pi_C(z)-\Pi_C(w))^\top(z-w).
$$

Cauchy–Schwarz proves nonexpansiveness. Also $x=\Pi_C(x-\eta\nabla f(x))$ is equivalent to $\nabla f(x)^\top(y-x)\ge0$ for all feasible $y$. For convex $f$, this is precisely the global first-order optimality condition.

Projected gradient therefore checks a constrained stationarity condition. The gradient itself need not vanish at a boundary optimum.

#### Example NM-E20: projecting onto a simplex

Project $z=(1.2,0.3,-0.2)$ onto $C=\{x:x_i\ge0,\ \sum_i x_i=1\}$.

Solution. The KKT conditions for minimizing $\tfrac12\|x-z\|^2$ give $x_i=\max(z_i-\tau,0)$ with a scalar $\tau$ chosen to make the coordinates sum to one.

Assume the first two coordinates remain positive and the third is zero. Then $1.2-\tau+0.3-\tau=1$ gives $\tau=0.25$. The proposed projection is $(0.95,0.05,0)$. Both active positive coordinates exceed zero, and the excluded coordinate has $z_3-\tau=-0.45\le0$, so the assumed active set is consistent.

Simply clipping negative entries and renormalizing would give $(0.8,0.2,0)$, a different vector. Euclidean projection subtracts a common threshold from the positive coordinates; it does not generally rescale them proportionally.

### 4.3 Worked examples

#### Example NM-E07: one fast and one slow direction

Let $H=\operatorname{diag}(1,9)$ and choose $\alpha=0.2$. Error multipliers are $0.8$ and $-0.8$, so both coordinates converge at rate $0.8$, with the second alternating. At $\alpha=0.25$, the second multiplier is $-1.25$ and generic iterates diverge. The optimal constant step is $2/(1+9)=0.2$.

#### Example NM-E08: projecting onto a half-line

Minimize $(x+1)^2/2$ over $x\ge0$. With step one, the unconstrained update is $x-(x+1)=-1$ for any start. Projection clips it to zero, the constrained optimum. At zero the gradient is one rather than zero; the constraint blocks the downhill negative direction.

### 4.4 Additional solved problems

#### Problem NM-19: exact line search on a quadratic

Constructed problem. For gradient $g=Hx-b\ne0$, minimize $f(x-\alpha g)$ over $\alpha$.

Solution. Expand the quadratic: $f(x-\alpha g)=f(x)-\alpha g^\top g+\tfrac12\alpha^2g^\top Hg$. The coefficient $g^\top Hg$ is positive, so the unique minimizer is $\alpha=(g^\top g)/(g^\top Hg)$. The next gradient is orthogonal to the old gradient because $g^\top[g-\alpha Hg]=0$. This does not mean consecutive parameter errors are orthogonal.

#### Problem NM-20: a box projection

Constructed problem. Project $z=(-2,1/2,4)$ onto $[0,1]^3$.

Solution. The squared distance is the sum of three coordinate squared distances, and the box constraints separate by coordinate. Each coordinate is clipped to the interval, giving $(0,1/2,1)$. The minimum squared distance is $4+0+9=13$. For a simplex or an ellipsoid, coordinatewise clipping generally fails because the constraints couple coordinates.

#### Problem NM-21: a simplex projection

Constructed problem. Project $(2,1)$ onto $\{x_1,x_2\ge0:x_1+x_2=1\}$.

Solution. Substitute $x_2=1-x_1$ with $0\le x_1\le1$. The distance objective is $(x_1-2)^2+(-x_1)^2=2x_1^2-4x_1+4$. Its unconstrained minimum is at one, which is feasible. The projection is $(1,0)$. This happens to be an endpoint; a general simplex projection finds a common threshold and clips, rather than clipping to a box and then normalizing arbitrarily.

#### Problem NM-22: a stopping bound using strong convexity

Constructed problem. If $f$ is $m$-strongly convex and differentiable, bound its optimality gap using the gradient.

Solution. Strong convexity gives $f(y)\ge f(x)+\nabla f(x)^\top(y-x)+(m/2)\|y-x\|^2$. Minimize the right-hand side over all $y$; completing the square gives lower bound $f(x)-\|\nabla f(x)\|^2/(2m)$ on $f^*$. Hence $f(x)-f^*\le\|\nabla f(x)\|^2/(2m)$. The estimate is for unconstrained minimization; constrained stationary points need a projected-gradient or KKT residual instead.

## 5. Differentiation and integration

### 5.1 Theory

The forward difference $[f(x+h)-f(x)]/h$ approximates $f'(x)$ with first-order truncation error for a sufficiently smooth function. The centered difference $[f(x+h)-f(x-h)]/(2h)$ cancels even Taylor terms and has second-order error.

The step cannot usually be made arbitrarily small in floating-point arithmetic. If function evaluations have absolute errors of order $u$, division by $h$ can amplify those errors to order $u/h$. A centered derivative's total error often has the form $C_1h^2+C_2u/h$, suggesting a balanced scale $h\asymp u^{1/3}$ after accounting for the function's size and derivatives.

Quadrature approximates an integral by a weighted sum of function values. Composite trapezoidal integration is second-order for functions with bounded continuous second derivative. Composite Simpson is fourth-order under appropriate fourth-derivative smoothness and an even number of equal subintervals. Endpoint singularities can reduce these nominal orders.

Monte Carlo integration represents an integral as an expectation and averages random evaluations. For iid finite-variance draws, its standard deviation decreases as $N^{-1/2}$. This rate depends on variance rather than a spatial grid dimension, although high dimension can make the variance large or the integrand hard to sample.

### 5.2 Derivations and proofs

#### Difference truncation errors

Taylor expansion gives

$$
\frac{f(x+h)-f(x)}h
=f'(x)+\frac h2f''(x)+O(h^2),
$$

$$
\frac{f(x+h)-f(x-h)}{2h}
=f'(x)+\frac{h^2}{6}f'''(x)+O(h^4).
$$

The second expansion with its $O(h^4)$ remainder requires sufficient fifth-order smoothness. A weaker third-derivative bound still yields a second-order error bound without asserting all higher coefficients.

If each of two centered-function evaluations has absolute error at most $\eta$, their difference divided by $2h$ has perturbation at most $\eta/h$. Combining with a third-derivative bound $M_3$ gives $M_3h^2/6+\eta/h$. Differentiating this bound gives $h^3=3\eta/M_3$ when both terms are positive and the local smoothness region permits that step.

#### Trapezoidal error

On a single interval of length $h$, expand about its midpoint. The exact integral and endpoint trapezoid share the constant and odd contributions, while their quadratic contributions differ by $h^3f''/12$. With $|f''|\le M_2$, the local error magnitude is at most $M_2h^3/12$. Summing $(b-a)/h$ equal cells gives

$$
|T_h-I|\le\frac{(b-a)M_2h^2}{12}.
$$

For sufficient higher smoothness, the leading signed term is $T_h-I=(h^2/12)[f'(b)-f'(a)]+O(h^4)$. This expansion supports extrapolation and endpoint correction.

#### Monte Carlo variance and a control variate

For $I=E[f(U)]$ and iid draws, $\widehat I=N^{-1}\sum_if(U_i)$ is unbiased with variance $\operatorname{Var}(f(U))/N$. If $g$ has known mean $\mu_g$, use

$$
\widehat I_c=\frac1N\sum_i[f(U_i)-c(g(U_i)-\mu_g)].
$$

It remains unbiased for fixed $c$. Its per-draw variance is $\operatorname{Var}(f)-2c\operatorname{Cov}(f,g)+c^2\operatorname{Var}(g)$. Minimize this quadratic to obtain $c^*=\operatorname{Cov}(f,g)/\operatorname{Var}(g)$ when the denominator is positive. The minimized variance is $\operatorname{Var}(f)(1-\rho_{fg}^2)$. Estimating $c$ on the same sample requires an additional finite-sample analysis, though regular plug-in choices can retain first-order efficiency.

#### Quadrature error comes from interpolation error

On one interval of width $h$, the trapezoidal rule integrates the straight-line interpolant through the endpoint values. For a twice continuously differentiable function, the interpolation error and integration give a local trapezoidal excess $h^3f''(\xi)/12$. Summing equal panels yields

$$
T_h-\int_a^bf(x)dx=\frac{(b-a)h^2}{12}f''(\xi)
$$

for some $\xi\in(a,b)$. In particular, the rule overestimates a convex function's integral. The sign is a useful check.

Composite Simpson integration on an even number of equal panels has error $I-S_h=-(b-a)h^4f^{(4)}(\xi)/180$ under a continuous fourth derivative. The stronger smoothness condition is essential. Applying a nominal fourth-order formula to an endpoint singularity does not create four derivatives there.

For centered differentiation, truncation is typically $C_1h^2$ while subtraction and division can produce arithmetic error of order $C_2u/h$. Balancing gives $h$ on the scale of $u^{1/3}$ after including the problem-dependent constants. This is a scale argument, not a universal numerical step for every function.

#### Example NM-E21: removing an integrable endpoint singularity

Evaluate $\int_0^1x^{-1/2}dx$ through a change of variables.

Solution. The integral is improper but finite. Set $x=t^2$ for $t\in[0,1]$. For positive $t$, $dx=2t\,dt$ and $x^{-1/2}=1/t$, so the transformed integrand is the constant two. Its continuous extension at zero is also two. Thus the integral is two.

An endpoint rule applied directly to $x^{-1/2}$ attempts to evaluate an infinite value at zero. The transformed rule has no singularity and integrates exactly even with one trapezoidal panel.

For $\int_0^1\sqrt x\,dx$, the same substitution produces $\int_0^12t^2dt=2/3$. The original derivative is singular at zero; the transformed integrand is a polynomial. A substitution can improve smoothness before discretization rather than asking a high-order rule to compensate for missing regularity.

The Jacobian and orientation of the substitution are parts of the integral. Replacing the function values without transforming $dx$ would solve a different problem.

### 5.3 Worked examples

#### Example NM-E09: a centered derivative of an exponential

For $f(x)=e^x$ at zero, the centered derivative is $(e^h-e^{-h})/(2h)=\sinh(h)/h=1+h^2/6+O(h^4)$. Halving $h$ reduces the leading truncation error by approximately four before roundoff becomes dominant.

#### Example NM-E10: integrating a square

For $\int_0^1x^2dx=1/3$, the composite trapezoid on $n$ equal cells satisfies $T_n=1/3+1/(6n^2)$. This follows by substituting $\sum_{j=1}^{n-1}j^2=(n-1)n(2n-1)/6$ into the quadrature sum. At $n=2$, the result is $3/8$ and error $1/24$.

#### Example NM-E11: controlling Monte Carlo variance

For $U\sim U(0,1)$, take $f(U)=U^2$ and $g(U)=U$, whose mean is one half. The covariance is $EU^3-EU^2EU=1/4-1/6=1/12$, equal to $\operatorname{Var}(U)$. Thus $c^*=1$. The raw per-draw variance is $1/5-1/9=4/45$, while the controlled variance is $4/45-1/12=1/180$, a factor-sixteen reduction.

### 5.4 Selected exercises

#### Problem NM-23: a derivative from a quadratic interpolant

Source: [FNC, Exercise 5.4.1(a,b)](https://fncbook.com/python/finitediffs/).

Problem. Construct the quadratic interpolant through values at $-h,0,h$ and differentiate it to obtain a centered derivative approximation.

Solution. With $f_-=f(-h)$, $f_0=f(0)$, $f_+=f(h)$, write

$$
Q(x)=f_0+\frac{f_+-f_-}{2h}x
+\frac{f_--2f_0+f_+}{2h^2}x^2.
$$

Substitution at zero gives $f_0$. At $h$, the linear and quadratic contributions combine to $f_+-f_0$; at $-h$, they combine to $f_--f_0$. Thus it interpolates all three data. Differentiating at zero gives $Q'(0)=(f_+-f_-)/(2h)$, the centered difference. The approximation arises by differentiating the interpolating polynomial.

#### Problem NM-24: consistency conditions for difference weights

Source: [FNC, Exercise 5.4.7(a,b)](https://fncbook.com/python/finitediffs/).

Problem. Establish centered-difference consistency and derive the constant- and linear-moment conditions for a general finite first-derivative stencil.

Solution. For the centered formula, the numerator and denominator vanish at $h=0$. Differentiating with respect to $h$ gives limit $[f'(x+h)+f'(x-h)]/2\to f'(x)$.

For a general formula $h^{-1}\sum_ka_kf(x+kh)$, Taylor expansion begins with $f(x)\sum_ka_k/h+f'(x)\sum_kka_k+O(h)$. To converge to the derivative for every smooth function, the necessary conditions are $\sum_ka_k=0$ and $\sum_kka_k=1$. They are also sufficient for first-order consistency with fixed finite weights and enough smoothness. The first removes constant functions; the second differentiates a linear function correctly.

#### Problem NM-25: a fourth-order centered stencil

Source: [FNC, Exercise 5.5.6(a,b)](https://fncbook.com/python/fd-converge/).

Problem. Derive the fourth-order centered five-point first-derivative weights and its leading truncation term.

Solution. On offsets $-2,-1,0,1,2$, require $\sum k^ja_k=0$ for $j=0,2,3,4$ and one for $j=1$. These are five linear equations for five weights. The antisymmetric choice $a_{-k}=-a_k$, $a_0=0$ automatically satisfies the even-power conditions. The remaining equations are $2a_1+4a_2=1$ and $2a_1+16a_2=0$, yielding $a_1=2/3$, $a_2=-1/12$. Hence

$$
f'(x)\approx
\frac{f(x-2h)-8f(x-h)+8f(x+h)-f(x+2h)}{12h}.
$$

The fifth moment of the weights is $-4$, so its leading error is $-h^4f^{(5)}(x)/30$ when the expansion is justified.

#### Problem NM-26: endpoint-corrected trapezoidal integration

Source: [FNC, Exercise 5.6.2](https://fncbook.com/python/integration/).

Problem. Cancel the leading composite trapezoidal error with endpoint derivatives, then approximate those derivatives to obtain an endpoint-corrected Gregory rule.

Solution. Subtract $h^2[f'(b)-f'(a)]/12$ from $T_h$ to remove its leading error. Approximate the endpoint derivatives by second-order one-sided formulas:

$$
f'(a)\approx\frac{-3f_0+4f_1-f_2}{2h},\qquad
f'(b)\approx\frac{3f_n-4f_{n-1}+f_{n-2}}{2h}.
$$

The resulting Gregory formula is

$$
G_h=T_h-\frac h{24}
[3(f_n+f_0)-4(f_{n-1}+f_1)+(f_{n-2}+f_2)].
$$

Each endpoint derivative error is $O(h^2)$, and multiplying it by $h^2$ leaves an $O(h^4)$ effect. Combined with the trapezoid expansion's next term, this gives fourth-order accuracy under the required smoothness.

#### Problem NM-27: deriving Simpson weights

Source: [FNC, Exercise 5.6.4(a–d)](https://fncbook.com/python/integration/).

Problem. Integrate a quadratic interpolant over a pair of panels to derive Simpson's local and composite weights.

Solution. The quadratic interpolating values $\alpha,\beta,\gamma$ at $-h,0,h$ is $p(x)=\beta+(\gamma-\alpha)x/(2h)+(\alpha-2\beta+\gamma)x^2/(2h^2)$. Substitution verifies its three values. Integrating over $[-h,h]$ removes the odd term and gives

$$
\int_{-h}^hp(x)dx
=2h\beta+\frac h3(\alpha-2\beta+\gamma)
=\frac h3(\alpha+4\beta+\gamma).
$$

Translate each consecutive pair of mesh intervals to this centered coordinate and sum. Shared endpoints receive weight two, odd-index interior midpoints receive weight four, and outer endpoints receive weight one. Thus for even $n$,

$$
S_h=\frac h3[f_0+4f_1+2f_2+\cdots+4f_{n-1}+f_n].
$$

### 5.5 Additional solved problems

#### Problem NM-28: Richardson extrapolation

Constructed problem. If $A(h)=L+ch^2+dh^4+O(h^6)$, eliminate the $h^2$ term.

Solution. $A(h/2)=L+ch^2/4+dh^4/16+O(h^6)$. Therefore $[4A(h/2)-A(h)]/3=L-dh^4/4+O(h^6)$, a fourth-order approximation. The improvement depends on the stated error expansion; it need not occur for a nonsmooth problem whose leading error is not a fixed multiple of $h^2$.

#### Problem NM-29: antithetic integration

Constructed problem. For $\int_0^1u^2du$, use paired average $Y=[U^2+(1-U)^2]/2$.

Solution. $Y=U^2-U+1/2$ has mean $1/3$ and variance $1/180$. Averaging $N$ independent pairs gives variance $1/(180N)$. At the same $2N$ function-evaluation count, ordinary Monte Carlo has variance $(4/45)/(2N)=2/(45N)$, eight times larger. Count computational effort consistently when comparing variance-reduction procedures.

#### Problem NM-30: an importance-sampling condition

Constructed problem. For $I=\int f(x)p(x)dx$, derive an estimator using draws from density $q$.

Solution. Where $q>0$ whenever $f p\ne0$, write $I=E_q[f(X)p(X)/q(X)]$. The iid sample average of these weights is unbiased when integrable. Its variance is $N^{-1}[\int f(x)^2p(x)^2/q(x)dx-I^2]$ when finite. A proposal with very small $q$ in an important region can produce enormous or infinite variance even if it visits that region occasionally.

#### Problem NM-31: quadrature with a nonsmooth endpoint

Constructed problem. Explain why a fourth-order rule need not achieve fourth order for $\sqrt x$ on $[0,1]$.

Solution. The derivatives required by the usual Simpson error bound are unbounded at zero. A numerical result can still converge, but the proof giving $O(h^4)$ does not apply. The substitution $x=t^2$ changes $\int_0^1\sqrt x\,dx$ into $\int_0^12t^2dt$, a smooth polynomial integral. A variable transformation can improve the mathematical approximation problem before changing the quadrature rule.

## 6. Numerical ordinary differential equations

### 6.1 Theory

For $y'=f(t,y)$, $y(0)=y_0$, explicit Euler uses $y_{n+1}=y_n+hf(t_n,y_n)$. Implicit Euler uses $y_{n+1}=y_n+hf(t_{n+1},y_{n+1})$, which requires solving an equation at every step.

The local one-step defect compares an exact solution advanced one step with the numerical formula started from its exact current value. For Euler it is $O(h^2)$. Some texts divide this defect by $h$ and call that the local truncation error, giving $O(h)$. The global error at a fixed time after many steps is first-order under the usual Lipschitz and smoothness assumptions.

Runge-Kutta methods combine several derivative evaluations. The explicit midpoint and Heun methods are second-order; merely using two evaluations does not guarantee that order. Coefficients must satisfy the appropriate Taylor-matching conditions.

Absolute stability concerns the numerical solution of $y'=\lambda y$ for $\operatorname{Re}\lambda<0$. A method can be high-order on smooth solutions and still unstable for a step too large relative to rapidly decaying components. Stiffness makes explicit stability limits much smaller than the time scale on which the main solution changes.

### 6.2 Derivations and proofs

#### Euler's global error

Let $f$ be Lipschitz in $y$ with constant $L$, and let the exact one-step defect be at most $Ch^2$. If $e_n=y(t_n)-y_n$, subtraction gives

$$
|e_{n+1}|\le(1+Lh)|e_n|+Ch^2.
$$

With exact initial data, iterate the bound:

$$
|e_n|\le Ch^2\sum_{j=0}^{n-1}(1+Lh)^j
\le\frac{Ch}{L}(e^{LT}-1),\quad nh\le T,
$$

for $L>0$. For $L=0$, the sum gives $CTh$. A second-order local defect has accumulated over order $1/h$ steps, yielding first-order global error.

#### Stability of Euler methods

For $y'=\lambda y$, explicit Euler gives $y_n=(1+h\lambda)^ny_0$. Stability requires $|1+h\lambda|\le1$; strict inequality gives decay of that mode. For real $\lambda<0$, this means $0<h\le2/|\lambda|$, with a nondecaying oscillation at the upper endpoint.

Implicit Euler gives amplification $1/(1-h\lambda)$. If $\operatorname{Re}\lambda\le0$, then $|1-h\lambda|\ge1$, so it is stable throughout the left half-plane. This is A-stability. It also damps very negative real modes to zero as $h|\lambda|\to\infty$.

#### A second-order step

Heun uses predictor $\widetilde y=y_n+hf(t_n,y_n)$ and correction $y_{n+1}=y_n+\tfrac h2[f(t_n,y_n)+f(t_n+h,\widetilde y)]$. Expand the second evaluation:

$$
f(t+h,y+hf)=f+h(f_t+f_yf)+O(h^2).
$$

Thus the update is $y+hf+\tfrac12h^2(f_t+f_yf)+O(h^3)$, matching the exact Taylor expansion through second order. The local defect is third-order and the global order is two under stability and smoothness assumptions.

#### Local defect accumulates into global error

Suppose a one-step method has local error at most $Ch^{p+1}$ when started from the exact solution and a perturbation amplification bound $1+Lh$. With global error $e_n$,

$$
|e_{n+1}|\le(1+Lh)|e_n|+Ch^{p+1}.
$$

Starting from exact initial data and summing the geometric series gives, for $L>0$,

$$
|e_n|\le
\frac{(1+Lh)^n-1}{L}Ch^p
\le\frac{e^{LT}-1}{L}Ch^p,\qquad nh\le T.
$$

For $L=0$, the bound becomes $CT h^p$. One factor of $h$ is lost because about $T/h$ local defects accumulate. This explains why an Euler local defect of order $h^2$ produces first-order global error.

The bound also separates consistency from stability. Small local defects can still become large if the numerical propagation amplifies them strongly. A stable method's amplification estimate must be suitable for the equation and the chosen step.

#### Why Heun has second-order accuracy

Let $k_1=f(t,y)$ and $k_2=f(t+h,y+hk_1)$. A smooth expansion gives $k_2=f+h(f_t+f_yf)+O(h^2)$. The update

$$
y_{\mathrm{new}}=y+\frac h2(k_1+k_2)
=y+hf+\frac{h^2}{2}(f_t+f_yf)+O(h^3)
$$

matches the exact solution's Taylor expansion through order two because $y''=f_t+f_yf$. Thus the local defect is third order and, under a stability estimate, the global method is second order. Counting two evaluations alone would not prove this matching.

#### Example NM-E22: stiffness despite a smooth exact path

Consider $y'=-1000(y-\cos t)-\sin t$ with $y(0)=1$. The exact solution is $y(t)=\cos t$, as direct differentiation verifies.

A perturbation $e$ from this path obeys $e'=-1000e$. Explicit Euler multiplies that perturbation by $1-1000h$, so stability requires $h\le0.002$ on the negative real test axis. With $h=0.01$ the multiplier is minus nine.

The desired path changes on an order-one time scale, but tiny truncation or rounding perturbations can excite a much faster mode. This is a basic stiffness mechanism. Backward Euler damps that mode for every positive step, although its step still must resolve the cosine path accurately.

A stable implicit method removes this particular stability restriction; it does not make arbitrary large steps accurate. The linear solve used at each implicit step also has its own error, which must remain controlled relative to the time discretization.

### 6.3 Worked examples

#### Example NM-E12: stable and unstable decay

For $y'=-10y$ and $y_0=1$, explicit Euler with $h=0.1$ has amplification zero. With $h=0.25$ it has amplification $-1.5$ and grows in magnitude while the exact solution decays. Implicit Euler at $h=0.25$ has amplification $1/3.5$, which decays but may have significant discretization error.

#### Example NM-E13: a midpoint step

For $y'=y$, $y_0=1$, and $h=0.2$, explicit midpoint evaluates the slope at predicted midpoint $1+0.1=1.1$, giving $y_1=1+0.2(1.1)=1.22$. Euler gives $1.2$, while $e^{0.2}\approx1.22140$. The second-order method is more accurate for this step, with no claim that it is unconditionally stable.

### 6.4 Selected exercises

#### Problem NM-32: four two-step Euler calculations

Source: [FNC, Exercise 6.2.1(a–d)](https://fncbook.com/python/euler/).

Problem. Take two explicit Euler steps for each initial-value problem: (a) $u'=-2tu$, $u(0)=2$, $h=0.1$; (b) $u'=u+t$, $u(0)=2$, $h=0.2$; (c) $u'=(1-u)/t$, $u(1)=6$, $h=0.25$; (d) $u'=2u(1-u)$, $u(0)=1/2$, $h=0.25$. Compare each terminal approximation with its exact solution.

Solution. In each case update from the old time and old state.

(a) For $u'=-2tu$, $u(0)=2$, $h=0.1$, the first step stays at two and the second is $2+0.1[-2(0.1)(2)]=1.96$. The exact value at $0.2$ is $2e^{-0.04}\approx1.921579$, so numerical-minus-exact error is about $0.038421$.

(b) For $u'=u+t$, $u(0)=2$, $h=0.2$, the steps are $2.4$ and $2.4+0.2(2.4+0.2)=2.92$. Exact $u(0.4)=-1.4+3e^{0.4}\approx3.075474$, so the error is about $-0.155474$.

(c) For $u'=(1-u)/t$, $u(1)=6$, $h=0.25$, the first value is $4.75$. The next is $4.75+0.25(1-4.75)/1.25=4$. Exact $u(1.5)=1+5/1.5=13/3$, so the error is $-1/3$.

(d) For $u'=2u(1-u)$, $u(0)=1/2$, $h=0.25$, the first value is $0.625$. The second is $0.625+0.25[2(0.625)(0.375)]=0.7421875$. Exact $u(0.5)=1/(1+e^{-1})\approx0.731059$, so the error is about $0.011129$.

#### Problem NM-33: two evaluations without second-order accuracy

Source: [FNC, Exercise 6.2.3(a,b)](https://fncbook.com/python/euler/).

Problem. For the method predicting $v=y+hf(t,y)$ and updating $y^+=y+hf(t+h,v)$, establish consistency and determine its local and global orders.

Solution. The source method predicts $v=y+hf(t,y)$, then uses $y^+=y+hf(t+h,v)$. Its one-step function is $\phi(t,y,h)=f(t+h,y+hf(t,y))$. As $h\to0$, $\phi\to f(t,y)$, establishing consistency.

Taylor expansion gives $y^+=y+hf+h^2(f_t+f_yf)+O(h^3)$. The exact coefficient of the second-order term is one half, not one. Thus the local defect is generally $O(h^2)$ and the global method is only first-order. The average of the two slopes in Heun is what supplies the missing factor.

#### Problem NM-34: explicit Euler's unstable test equation

Source: [FNC, Exercise 6.2.4(a,b)](https://fncbook.com/python/euler/).

Problem. Solve explicit Euler's recurrence on $u'=ku$ and give a negative $k$ and fixed step for which the exact solution decays but the numerical values diverge.

Solution. The recurrence $u_{i+1}=(1+hk)u_i$, $u_0=1$, gives $u_i=(1+hk)^i$ at time $ih$. Choose $k=-1$ and $h=3$. The exact solution $e^{-t}$ is bounded and decays, but the numerical solution is $(-2)^i$, whose magnitude diverges. The fixed-step long-time instability is distinct from convergence as $h\to0$ over a fixed finite interval.

### 6.5 Additional solved problems

#### Problem NM-35: accumulated roundoff in Euler

Source: [FNC, Exercise 6.2.6](https://fncbook.com/python/euler/).

Problem. Include a uniformly bounded roundoff perturbation in each step of an Euler-type method and derive the resulting global truncation-plus-roundoff bound.

Solution. Let the per-step defect be bounded by $Ch^{p+1}+u$, with roundoff perturbation magnitude at most $u$. The error recurrence becomes $|e_{n+1}|\le(1+Lh)|e_n|+Ch^{p+1}+u$. With exact initial value,

$$
|e_n|\le\frac{e^{LT}-1}{L}(Ch^p+u/h)
$$

for $L>0$, and the limiting $T(Ch^p+u/h)$ form for $L=0$. This worst-case bound allows all roundoff errors to align. Smaller $h$ reduces the truncation term but enlarges the accumulated-roundoff term.

#### Problem NM-36: error versus stability

Constructed problem. Is an implicit Euler solution accurate merely because it is stable at a large step?

Solution. For $y'=-y$ with $h=10$, implicit Euler gives one-step multiplier $1/11\approx0.0909$, while the exact multiplier is $e^{-10}\approx0.0000454$. Both decay and the numerical method is stable, but the relative discrepancy in this mode is large. Stability controls amplification; accuracy also requires sufficiently small discretization error for the intended quantity.

#### Problem NM-37: a stiff pair

Constructed problem. Consider $\dot x=-x$, $\dot y=-1000y$.

Solution. Explicit Euler must use $h<0.002$ for strict stability of the fast component, even if only the slowly varying $x$ is of primary interest after the rapid transient. An implicit method can remain stable at a larger step, but should still resolve any transient whose accuracy matters. Stiffness describes the mismatch between fast numerical stability restrictions and slower solution features.

#### Problem NM-38: order estimation from step halving

Constructed problem. If global errors at $h,h/2,h/4$ are approximately $0.04,0.01,0.0025$, estimate the order.

Solution. The error ratio is four under each halving. For a leading form $Ch^p$, the ratio is $2^p$, so $p\approx2$. Without an exact reference, differences of successive approximations can estimate the ratio if they are in the asymptotic regime and share the same leading coefficient. A short table alone does not prove the order for every problem.

## 7. Numerical partial differential equations

### 7.1 Theory

For the heat equation $u_t=\nu u_{xx}$, use space grid $x_j=jh$ and time grid $t_n=nk$. Centered space differences and a forward time difference give

$$
U_j^{n+1}=U_j^n+\lambda(U_{j+1}^n-2U_j^n+U_{j-1}^n),
\qquad\lambda=\nu k/h^2.
$$

Initial values populate the first time level. Boundary conditions supply values or relations at the spatial edges at every later time. An interior formula alone does not define a complete numerical problem.

Consistency means substituting a smooth exact solution into the discrete equation leaves a defect tending to zero. Stability bounds how errors or perturbations are amplified by repeated discrete steps. Convergence means the numerical solution approaches the exact one as the grid is refined. For suitable linear well-posed initial-value problems and consistent linear schemes, the Lax equivalence theorem connects stability and convergence; its functional setting and boundary treatment are part of the theorem.

For the heat scheme above, a sufficient maximum-norm stability condition is $0\le\lambda\le1/2$. At that condition the coefficients $\lambda,1-2\lambda,\lambda$ are nonnegative and sum to one. Each interior update is a weighted average of nearby values.

### 7.2 Derivations and proofs

#### Fourier-mode stability

On an infinite or periodic spatial grid, insert $U_j^n=G^ne^{ij\theta}$. The amplification factor is

$$
G=1+\lambda(e^{i\theta}-2+e^{-i\theta})
=1-4\lambda\sin^2(\theta/2).
$$

For all modes to have $|G|\le1$, require $0\le\lambda\le1/2$. If $\lambda>1/2$, the highest-frequency mode $\theta=\pi$ has $G=1-4\lambda<-1$ and grows with alternating sign.

#### Consistency and a convergence bound

For a sufficiently smooth exact solution, the time difference is $u_t+(k/2)u_{tt}+O(k^2)$ and the centered space difference is $u_{xx}+(h^2/12)u_{xxxx}+O(h^4)$. The normalized equation defect is $O(k+h^2)$.

With exact compatible initial and Dirichlet boundary values and $\lambda\le1/2$, the error update is a sup-norm contraction plus a per-step forcing of size at most $Ck(k+h^2)$. Summing over at most $T/k$ steps gives $\|e^n\|_\infty\le CT(k+h^2)$, proving convergence on a fixed horizon. Boundary approximation errors must be added if the boundary conditions are themselves discretized.

#### Backward time in an option PDE

For a Black-Scholes terminal-value problem, set $\tau=T-t$. Then $v_\tau=\tfrac12\sigma^2S^2v_{SS}+rSv_S-rv$. This is forward evolution in time-to-maturity with payoff as the initial condition.

On a log-price grid $x=\log S$, the equation becomes $v_\tau=av_{xx}+bv_x-rv$, with $a=\sigma^2/2$ and $b=r-\sigma^2/2$. Put $v=e^{\alpha x+\gamma\tau}u$, where $\alpha=-b/(2a)$ and $\gamma=-b^2/(4a)-r$. Differentiating and substituting cancels the first-derivative and zeroth-order terms, leaving $u_\tau=au_{xx}$. The heat restriction is then $ak/h^2\le1/2$.

A direct central-difference scheme with convection has additional sign and stability issues. The pure heat restriction should not be quoted as a universal sufficient rule for every Black-Scholes discretization. Boundary truncation at finite log prices also creates an error separate from the interior grid error.

#### Fourier modes and a full convergence estimate

For the explicit heat scheme, substitute a spatial mode $u_j^n=g^ne^{ij\theta}$. The amplification factor is

$$
g(\theta)=1+r(e^{i\theta}-2+e^{-i\theta})
=1-4r\sin^2(\theta/2).
$$

Requiring $|g(\theta)|\le1$ for every frequency gives $0\le r\le1/2$. The same restriction makes the real-space stencil a convex average, proving maximum-norm stability for compatible boundary values without Fourier analysis.

For a smooth exact solution sampled on the grid, define the divided local residual $\tau_j^n=D_tu-\nu D_{xx}u$. Taylor expansion gives $|\tau_j^n|\le C(k+h^2)$. The one-step error equation has forcing $k\tau_j^n$. Under the maximum-norm stability bound and exact boundary data,

$$
\|e^{n+1}\|_\infty
\le\|e^n\|_\infty+kC(k+h^2).
$$

Summing through $nk\le T$ gives $\|e^n\|_\infty\le\|e^0\|_\infty+TC(k+h^2)$. This is a convergence proof using the actual local residual and stable propagation. Boundary approximation errors must be included if the supplied boundary values are not exact.

On a fixed spatial mesh, diffusion is a stiff ODE system: the most negative eigenvalue has magnitude on the scale $\nu/h^2$. Refining space consequently forces an explicit time step on the scale $h^2$. This is the same stiffness issue seen in the ODE section.

#### Crank–Nicolson stability and high-frequency damping

For $y'=\lambda y$, Crank–Nicolson has factor $R(z)=(1+z/2)/(1-z/2)$, $z=h\lambda$. It is stable throughout the left half-plane. Along the negative real axis, however, $R(z)\to-1$ as $z\to-\infty$. Extremely fast modes can alternate with little decay under a large step.

Backward Euler has factor $1/(1-z)\to0$ on that axis, so it damps such modes much more strongly. Both are stable, but their treatment of unresolved high frequencies differs. Nonsmooth initial data, such as a payoff corner after reversing pricing time, can expose this difference.

#### Example NM-E23: a strongly decaying mesh mode

For a heat-system eigenvalue with $z=-20$, backward Euler gives factor $1/21$. Crank–Nicolson gives $(1-10)/(1+10)=-9/11$. Both magnitudes are below one, but after two steps the factors are $1/441$ and $81/121$.

The exact factor over two steps is $e^{-40}$, far smaller than either. Crank–Nicolson's second-order accuracy is an asymptotic statement as the step becomes small; it does not imply accurate damping of a mode with $|z|=20$. The example explains why stability, order, and damping should be inspected separately.

### 7.3 Worked examples

#### Example NM-E14: one explicit diffusion step

Let $\nu=1$, $h=1/4$, and $k=1/32$, so $\lambda=1/2$. With zero endpoint values and interior values $(1,2,1)$, the next interior values are $(1,1,1)$: the center becomes the average of its two neighbors, and each side becomes the average of the center and its zero boundary. All values remain within the previous range.

#### Example NM-E15: a forbidden refinement

If $h$ is halved while $k$ is held fixed, $\lambda$ grows by four. Starting from $\lambda=1/4$, this produces $\lambda=1$, which is unstable. To preserve a fixed explicit heat CFL number, reducing spatial spacing by two requires reducing the time step by four.

#### Example NM-E16: an exact decaying mode

For $u_t=\nu u_{xx}$ on $[0,1]$ with zero endpoints, $u(x,t)=e^{-\nu\pi^2t}\sin(\pi x)$ satisfies the PDE and boundary conditions. It provides an exact reference for checking the scheme's error without estimating a reference from another numerical grid.

### 7.4 Selected exercises

#### Problem NM-39: verifying a Neumann heat solution

Source: [FNC, Exercise 11.1.1](https://fncbook.com/python/blackscholes/).

Problem. Verify that $u(x,t)=e^{-4\pi^2t}\cos(2\pi x)$ solves the heat equation on $[0,1]$ with zero Neumann boundary conditions and its cosine initial value.

Solution. For $u=e^{-4\pi^2t}\cos(2\pi x)$, $u_t=-4\pi^2u$ and $u_{xx}=-4\pi^2u$, so $u_t=u_{xx}$. The spatial derivative is $-2\pi e^{-4\pi^2t}\sin(2\pi x)$, zero at both $x=0$ and $x=1$. At time zero the stated cosine initial condition is recovered. All three parts of the initial-boundary-value problem are satisfied.

#### Problem NM-40: a heat kernel calculation

Source: [FNC, Exercise 11.1.2](https://fncbook.com/python/blackscholes/).

Problem. Verify the heat equation for $t^{-1/2}e^{-x^2/(4t)}$ at positive times and explain its singular initial behavior.

Solution. Let $u=t^{-1/2}e^{-x^2/(4t)}$ for $t>0$. Logarithmic differentiation gives $u_t/u=-1/(2t)+x^2/(4t^2)$. Also $u_x/u=-x/(2t)$ and $u_{xx}/u=-1/(2t)+x^2/(4t^2)$. Thus $u_t=u_{xx}$. The time-zero limit is singular; it is not a classical initial function at the origin. Multiplication by $1/\sqrt{4\pi}$ gives the normalized Gaussian heat kernel.

#### Problem NM-41: a mixed-boundary eigenproblem

Source: [MIT 18.330, Homework 6, Problem 1(a)](https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/28d3c9b92505167a23832f5e9c1d11d1_MIT18_330S12_hw6.pdf), eigenvalues and eigenfunctions.

Problem. Find all eigenvalues and eigenfunctions of $-v''=\lambda v$ on $(0,1)$ with $v'(0)=0$ and $v(1)=0$, excluding trivial solutions.

Solution. Solve $-v''=\lambda v$ with $v'(0)=0$, $v(1)=0$. For $\lambda=\omega^2>0$, the general solution is $A\cos\omega x+B\sin\omega x$. The first boundary sets $B=0$ and the second requires $\cos\omega=0$. Thus $\omega_n=(n+1/2)\pi$, $\lambda_n=(n+1/2)^2\pi^2$, and $v_n(x)=\cos[(n+1/2)\pi x]$, $n=0,1,\ldots$. At $\lambda=0$, the affine solution is forced to zero. For negative $\lambda$, hyperbolic solutions and the two boundary conditions likewise give only zero. The first three nonzero modes have zero, one, and two interior sign changes respectively.

#### Problem NM-42: a Neumann compatibility condition

Source: [MIT 18.330, Homework 6, Problem 2(a,b)](https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/28d3c9b92505167a23832f5e9c1d11d1_MIT18_330S12_hw6.pdf).

Problem. For $-[(1+x^2)u']'=f$ on $(0,1)$ with zero Neumann derivatives, derive the compatibility condition, test $f=1$, and establish sufficiency and nonuniqueness for continuous compatible forcing.

Solution. Integrate $-[(1+x^2)u']'=f$ over $[0,1]$. The zero Neumann derivatives make the integrated left side zero, so $\int_0^1f(x)dx=0$ is necessary. Therefore $f=1$ is impossible.

For continuous $f$ satisfying that condition, integrate once:

$$
u'(x)=-\frac{\int_0^xf(s)ds}{1+x^2}.
$$

It satisfies the left derivative condition and, by zero total integral, the right one. Integrating once more constructs $u$ up to an arbitrary additive constant. This proves sufficiency in the stated regularity class and explains nonuniqueness without a level normalization.

### 7.5 Additional solved problems

#### Problem NM-43: implicit diffusion

Constructed problem. Derive backward Euler with centered differences for the heat equation.

Solution. Use $(U_j^{n+1}-U_j^n)/k=\nu(U_{j+1}^{n+1}-2U_j^{n+1}+U_{j-1}^{n+1})/h^2$. Rearranging gives

$$
-\lambda U_{j-1}^{n+1}+(1+2\lambda)U_j^{n+1}
-\lambda U_{j+1}^{n+1}=U_j^n.
$$

At each step solve a tridiagonal system, moving known new-time boundary values to the right. Fourier amplification is $1/[1+4\lambda\sin^2(\theta/2)]$, bounded by one for all $\lambda\ge0$. Unconditional stability does not remove its first-order time error.

#### Problem NM-44: Crank-Nicolson

Constructed problem. Average explicit and implicit spatial terms.

Solution. The scheme is $(U^{n+1}-U^n)/k=(\nu/2)D_{xx}(U^{n+1}+U^n)$. Its amplification is

$$
G=\frac{1-2\lambda\sin^2(\theta/2)}
{1+2\lambda\sin^2(\theta/2)}.
$$

It has magnitude at most one for every nonnegative $\lambda$ and is second-order in time for smooth solutions. For large $\lambda$ and high frequencies, $G$ approaches minus one, so oscillatory components can decay slowly. Stability and strong damping are separate properties.

#### Problem NM-45: the role of a boundary error

Constructed problem. Interior differences are second-order, but a Neumann boundary derivative is approximated by a first-order one-sided difference. State what must be analyzed.

Solution. The full discrete problem includes its boundary equations. An $O(h)$ boundary defect can reduce the global accuracy unless the inverse discrete operator and the defect's localization provide a sharper compensating estimate. A uniform bound on the discrete solution operator in the norm being used connects defects to solution errors. Quoting only the interior stencil's $O(h^2)$ consistency is insufficient.

#### Problem NM-46: a maximum principle

Constructed problem. Prove range preservation for explicit heat steps with zero boundary values and $0\le\lambda\le1/2$.

Solution. Each interior value is a convex combination of three old values with coefficients $\lambda,1-2\lambda,\lambda$. It lies between their minimum and maximum. Including the fixed zero boundaries, the new global maximum cannot exceed the old maximum or zero, and the new minimum cannot be below the old minimum or zero. Induction proves range preservation at every time step. Negative weights remove this argument.

#### Problem NM-47: a convergence rate under coupled refinement

Constructed problem. If the explicit heat error is $O(k+h^2)$ and stability is maintained with $k=ch^2$, find the combined rate.

Solution. Substitute to obtain $O(ch^2+h^2)=O(h^2)$ at a fixed final time. Halving $h$ requires roughly four times as many time steps and twice as many one-dimensional spatial points, making work grow by roughly eight for a simple explicit implementation. A mathematical accuracy order and a computational cost order describe different tradeoffs.

#### Problem NM-48: an error budget

Constructed problem. A numerical option value uses a finite spatial domain, a mesh, and an iterative linear solve. Identify separate errors.

Solution. Domain truncation approximates behavior outside the chosen price interval through artificial boundaries. Discretization replaces the PDE derivatives and time evolution by finite steps. Algebraic error comes from solving the resulting finite system only approximately. Roundoff contributes further arithmetic error. Tightening the linear-solver tolerance below the dominant domain or mesh error gives little improvement in the final price; each component needs a scale appropriate to the desired accuracy.

<!-- extended-practice-start -->

## 8. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 8.1 Cancellation checked against high precision

Equivalent real-number formulas need not behave equally in floating-point arithmetic. Subtracting two nearby numbers can discard meaningful digits. A high-precision reference helps separate this arithmetic effect from approximation error in the mathematical formula.

#### Problem NM-49: small square-root difference at 10 to minus 4

Evaluate $f(x)=\sqrt{1+x}-1$ at $x=10^{-4}$ using the direct formula and the rationalized formula. Compare relative errors with a high-precision evaluation.

Solution. Multiply numerator and denominator by $\sqrt{1+x}+1$:

$$
\sqrt{1+x}-1=\frac x{\sqrt{1+x}+1}.
$$

This is an exact identity for the positive input. The direct form subtracts quantities close to one. The rationalized denominator is close to two and contains no subtraction of nearby values.

Using binary64 arithmetic for both formulas and a 90-digit decimal calculation as reference gives:

| Quantity | Value |
|---|---:|
| High-precision reference, rounded for display | $4.99987500624960929e-05$ |
| Direct result | $4.99987500623966241e-05$ |
| Rationalized result | $4.99987500624960997e-05$ |
| Direct relative error | $1.989448e-12$ |
| Rationalized relative error | $1.128499e-16$ |

The direct formula first rounds $1+x$ and its square root near one. The subtraction then reveals the absolute error in those rounded quantities relative to a much smaller true result. At sufficiently small $x$, the rounded intermediate quantity can equal one exactly, producing zero.

The problem itself is well conditioned near zero when measured relatively for nonzero $x$. Indeed,

$$
\left|\frac{xf'(x)}{f(x)}\right|
=\frac{\sqrt{1+x}+1}{2\sqrt{1+x}}\longrightarrow1.
$$

Thus the lost accuracy comes from the chosen expression, rather than large sensitivity of the mathematical function to relative input error. The stable form still rounds and is not exact, but its arithmetic avoids the damaging cancellation. For inputs near the domain boundary $x=-1$, the conditioning analysis would be different.

#### Problem NM-50: small square-root difference at 10 to minus 8

Evaluate $f(x)=\sqrt{1+x}-1$ at $x=10^{-8}$ using the direct formula and the rationalized formula. Compare relative errors with a high-precision evaluation.

Solution. Multiply numerator and denominator by $\sqrt{1+x}+1$:

$$
\sqrt{1+x}-1=\frac x{\sqrt{1+x}+1}.
$$

This is an exact identity for the positive input. The direct form subtracts quantities close to one. The rationalized denominator is close to two and contains no subtraction of nearby values.

Using binary64 arithmetic for both formulas and a 90-digit decimal calculation as reference gives:

| Quantity | Value |
|---|---:|
| High-precision reference, rounded for display | $4.99999998749999989e-09$ |
| Direct result | $4.99999996961264515e-09$ |
| Rationalized result | $4.99999998749999989e-09$ |
| Direct relative error | $3.577471e-09$ |
| Rationalized relative error | $3.373340e-17$ |

The direct formula first rounds $1+x$ and its square root near one. The subtraction then reveals the absolute error in those rounded quantities relative to a much smaller true result. At sufficiently small $x$, the rounded intermediate quantity can equal one exactly, producing zero.

The problem itself is well conditioned near zero when measured relatively for nonzero $x$. Indeed,

$$
\left|\frac{xf'(x)}{f(x)}\right|
=\frac{\sqrt{1+x}+1}{2\sqrt{1+x}}\longrightarrow1.
$$

Thus the lost accuracy comes from the chosen expression, rather than large sensitivity of the mathematical function to relative input error. The stable form still rounds and is not exact, but its arithmetic avoids the damaging cancellation. For inputs near the domain boundary $x=-1$, the conditioning analysis would be different.

#### Problem NM-51: small square-root difference at 10 to minus 12

Evaluate $f(x)=\sqrt{1+x}-1$ at $x=10^{-12}$ using the direct formula and the rationalized formula. Compare relative errors with a high-precision evaluation.

Solution. Multiply numerator and denominator by $\sqrt{1+x}+1$:

$$
\sqrt{1+x}-1=\frac x{\sqrt{1+x}+1}.
$$

This is an exact identity for the positive input. The direct form subtracts quantities close to one. The rationalized denominator is close to two and contains no subtraction of nearby values.

Using binary64 arithmetic for both formulas and a 90-digit decimal calculation as reference gives:

| Quantity | Value |
|---|---:|
| High-precision reference, rounded for display | $4.99999999999874984e-13$ |
| Direct result | $5.00044450291170506e-13$ |
| Rationalized result | $4.99999999999874984e-13$ |
| Direct relative error | $8.890058e-05$ |
| Rationalized relative error | $3.222232e-17$ |

The direct formula first rounds $1+x$ and its square root near one. The subtraction then reveals the absolute error in those rounded quantities relative to a much smaller true result. At sufficiently small $x$, the rounded intermediate quantity can equal one exactly, producing zero.

The problem itself is well conditioned near zero when measured relatively for nonzero $x$. Indeed,

$$
\left|\frac{xf'(x)}{f(x)}\right|
=\frac{\sqrt{1+x}+1}{2\sqrt{1+x}}\longrightarrow1.
$$

Thus the lost accuracy comes from the chosen expression, rather than large sensitivity of the mathematical function to relative input error. The stable form still rounds and is not exact, but its arithmetic avoids the damaging cancellation. For inputs near the domain boundary $x=-1$, the conditioning analysis would be different.

#### Problem NM-52: small square-root difference at 10 to minus 16

Evaluate $f(x)=\sqrt{1+x}-1$ at $x=10^{-16}$ using the direct formula and the rationalized formula. Compare relative errors with a high-precision evaluation.

Solution. Multiply numerator and denominator by $\sqrt{1+x}+1$:

$$
\sqrt{1+x}-1=\frac x{\sqrt{1+x}+1}.
$$

This is an exact identity for the positive input. The direct form subtracts quantities close to one. The rationalized denominator is close to two and contains no subtraction of nearby values.

Using binary64 arithmetic for both formulas and a 90-digit decimal calculation as reference gives:

| Quantity | Value |
|---|---:|
| High-precision reference, rounded for display | $4.99999999999999990e-17$ |
| Direct result | $0.00000000000000000e+00$ |
| Rationalized result | $4.99999999999999990e-17$ |
| Direct relative error | $1.000000e+00$ |
| Rationalized relative error | $4.097787e-18$ |

The direct formula first rounds $1+x$ and its square root near one. The subtraction then reveals the absolute error in those rounded quantities relative to a much smaller true result. At sufficiently small $x$, the rounded intermediate quantity can equal one exactly, producing zero.

The problem itself is well conditioned near zero when measured relatively for nonzero $x$. Indeed,

$$
\left|\frac{xf'(x)}{f(x)}\right|
=\frac{\sqrt{1+x}+1}{2\sqrt{1+x}}\longrightarrow1.
$$

Thus the lost accuracy comes from the chosen expression, rather than large sensitivity of the mathematical function to relative input error. The stable form still rounds and is not exact, but its arithmetic avoids the damaging cancellation. For inputs near the domain boundary $x=-1$, the conditioning analysis would be different.

### 8.2 Elimination, pivots, and residual checks

An LU factorization records elimination in a lower-triangular multiplier matrix and an upper-triangular reduced system. Solving by forward and back substitution gives an exact algebraic answer. A small pivot can nevertheless make the unpivoted arithmetic poor, even when the final linear system is not very ill conditioned.

#### Problem NM-53: two-by-two LU solve, case 1

Solve $Ax=b$ for $A=\begin{pmatrix}2&1\\4&3\end{pmatrix}$ and $b=\begin{pmatrix}3\\7\end{pmatrix}$ using an unpivoted LU factorization. Check the residual and discuss the first pivot.

Solution. The first pivot is $a_{11}=2\ne0$. Eliminate the lower-left entry using multiplier $\ell=a_{21}/a_{11}=2$. This gives

$$
L=\begin{pmatrix}1&0\\2&1\end{pmatrix},\qquad
U=\begin{pmatrix}2&1\\0&1\end{pmatrix}.
$$

Multiplication verifies $LU=A$. Both diagonal entries of $U$ are nonzero, so the solve can proceed. Forward substitution in $Lz=b$ gives $z_1=b_1$ and $z_2=b_2-\ell b_1$, hence $z=\begin{pmatrix}3\\1\end{pmatrix}$.

Back substitution in $Ux=z$ yields

$$
x_2=z_2/u_{22}=1,\qquad
x_1=(z_1-u_{12}x_2)/u_{11}=1.
$$

Substituting the vector $x=\begin{pmatrix}1\\1\end{pmatrix}$ into the original equations gives $Ax=b$ exactly. With exact fractions, the residual is zero; floating-point substitution should instead report a small norm scaled by $\|A\|\|x\|+\|b\|$.

For comparison, the infinity-norm condition number is $\kappa_\infty(A)=21\approx 21.000000$. This measures sensitivity of the supplied system. The multiplier magnitude $|\ell|=2$ describes this particular elimination step. A tiny first pivot can create a large multiplier even when the condition number is moderate.

Partial pivoting compares the magnitudes in the current column and swaps rows before elimination when appropriate. It avoids unnecessarily large multipliers in these two-row examples. Pivoting does not remove intrinsic ill conditioning, and a small residual alone does not guarantee a small forward error for an ill-conditioned system.

#### Problem NM-54: two-by-two LU solve, case 2

Solve $Ax=b$ for $A=\begin{pmatrix}4&-1\\2&3\end{pmatrix}$ and $b=\begin{pmatrix}1\\5\end{pmatrix}$ using an unpivoted LU factorization. Check the residual and discuss the first pivot.

Solution. The first pivot is $a_{11}=4\ne0$. Eliminate the lower-left entry using multiplier $\ell=a_{21}/a_{11}=\frac{1}{2}$. This gives

$$
L=\begin{pmatrix}1&0\\\frac{1}{2}&1\end{pmatrix},\qquad
U=\begin{pmatrix}4&-1\\0&\frac{7}{2}\end{pmatrix}.
$$

Multiplication verifies $LU=A$. Both diagonal entries of $U$ are nonzero, so the solve can proceed. Forward substitution in $Lz=b$ gives $z_1=b_1$ and $z_2=b_2-\ell b_1$, hence $z=\begin{pmatrix}1\\\frac{9}{2}\end{pmatrix}$.

Back substitution in $Ux=z$ yields

$$
x_2=z_2/u_{22}=\frac{9}{7},\qquad
x_1=(z_1-u_{12}x_2)/u_{11}=\frac{4}{7}.
$$

Substituting the vector $x=\begin{pmatrix}\frac{4}{7}\\\frac{9}{7}\end{pmatrix}$ into the original equations gives $Ax=b$ exactly. With exact fractions, the residual is zero; floating-point substitution should instead report a small norm scaled by $\|A\|\|x\|+\|b\|$.

For comparison, the infinity-norm condition number is $\kappa_\infty(A)=\frac{15}{7}\approx 2.142857$. This measures sensitivity of the supplied system. The multiplier magnitude $|\ell|=\frac{1}{2}$ describes this particular elimination step. A tiny first pivot can create a large multiplier even when the condition number is moderate.

Partial pivoting compares the magnitudes in the current column and swaps rows before elimination when appropriate. It avoids unnecessarily large multipliers in these two-row examples. Pivoting does not remove intrinsic ill conditioning, and a small residual alone does not guarantee a small forward error for an ill-conditioned system.

#### Problem NM-55: two-by-two LU solve, case 3

Solve $Ax=b$ for $A=\begin{pmatrix}1&2\\3&7\end{pmatrix}$ and $b=\begin{pmatrix}4\\11\end{pmatrix}$ using an unpivoted LU factorization. Check the residual and discuss the first pivot.

Solution. The first pivot is $a_{11}=1\ne0$. Eliminate the lower-left entry using multiplier $\ell=a_{21}/a_{11}=3$. This gives

$$
L=\begin{pmatrix}1&0\\3&1\end{pmatrix},\qquad
U=\begin{pmatrix}1&2\\0&1\end{pmatrix}.
$$

Multiplication verifies $LU=A$. Both diagonal entries of $U$ are nonzero, so the solve can proceed. Forward substitution in $Lz=b$ gives $z_1=b_1$ and $z_2=b_2-\ell b_1$, hence $z=\begin{pmatrix}4\\-1\end{pmatrix}$.

Back substitution in $Ux=z$ yields

$$
x_2=z_2/u_{22}=-1,\qquad
x_1=(z_1-u_{12}x_2)/u_{11}=6.
$$

Substituting the vector $x=\begin{pmatrix}6\\-1\end{pmatrix}$ into the original equations gives $Ax=b$ exactly. With exact fractions, the residual is zero; floating-point substitution should instead report a small norm scaled by $\|A\|\|x\|+\|b\|$.

For comparison, the infinity-norm condition number is $\kappa_\infty(A)=90\approx 90.000000$. This measures sensitivity of the supplied system. The multiplier magnitude $|\ell|=3$ describes this particular elimination step. A tiny first pivot can create a large multiplier even when the condition number is moderate.

Partial pivoting compares the magnitudes in the current column and swaps rows before elimination when appropriate. It avoids unnecessarily large multipliers in these two-row examples. Pivoting does not remove intrinsic ill conditioning, and a small residual alone does not guarantee a small forward error for an ill-conditioned system.

#### Problem NM-56: two-by-two LU solve, case 4

Solve $Ax=b$ for $A=\begin{pmatrix}\frac{1}{10000}&1\\1&1\end{pmatrix}$ and $b=\begin{pmatrix}1\\2\end{pmatrix}$ using an unpivoted LU factorization. Check the residual and discuss the first pivot.

Solution. The first pivot is $a_{11}=\frac{1}{10000}\ne0$. Eliminate the lower-left entry using multiplier $\ell=a_{21}/a_{11}=10000$. This gives

$$
L=\begin{pmatrix}1&0\\10000&1\end{pmatrix},\qquad
U=\begin{pmatrix}\frac{1}{10000}&1\\0&-9999\end{pmatrix}.
$$

Multiplication verifies $LU=A$. Both diagonal entries of $U$ are nonzero, so the solve can proceed. Forward substitution in $Lz=b$ gives $z_1=b_1$ and $z_2=b_2-\ell b_1$, hence $z=\begin{pmatrix}1\\-9998\end{pmatrix}$.

Back substitution in $Ux=z$ yields

$$
x_2=z_2/u_{22}=\frac{9998}{9999},\qquad
x_1=(z_1-u_{12}x_2)/u_{11}=\frac{10000}{9999}.
$$

Substituting the vector $x=\begin{pmatrix}\frac{10000}{9999}\\\frac{9998}{9999}\end{pmatrix}$ into the original equations gives $Ax=b$ exactly. With exact fractions, the residual is zero; floating-point substitution should instead report a small norm scaled by $\|A\|\|x\|+\|b\|$.

For comparison, the infinity-norm condition number is $\kappa_\infty(A)=\frac{40000}{9999}\approx 4.000400$. This measures sensitivity of the supplied system. The multiplier magnitude $|\ell|=10000$ describes this particular elimination step. A tiny first pivot can create a large multiplier even when the condition number is moderate.

Partial pivoting compares the magnitudes in the current column and swaps rows before elimination when appropriate. It avoids unnecessarily large multipliers in these two-row examples. Pivoting does not remove intrinsic ill conditioning, and a small residual alone does not guarantee a small forward error for an ill-conditioned system.

### 8.3 Bisection guarantees and Newton speed

Bisection reduces a known sign-changing bracket at a predictable rate. Newton's method can be much faster near a simple root, but its speed needs a separate convergence argument. A positive square root gives an exact error identity that makes the comparison transparent.

#### Problem NM-57: computing the square root of 2

Find the positive root of $f(x)=x^2-2$. Give a bisection iteration count ensuring midpoint error at most $10^{-6}$ and show Newton iterates starting from $x_0=2$.

Solution. The function is continuous and strictly increasing for positive $x$. Its values at $1$ and $2$ have opposite signs, so this interval contains exactly one positive root.

After $n$ bisection reductions, the bracket width is $(2-1)/2^n$. The midpoint error is at most half that width. Thus it is enough to choose

$$
n\ge\log_2\frac{1}{2\cdot10^{-6}},
\qquad n=19.
$$

This count depends on the initial bracket and requested absolute error, not on the derivative at the root.

Newton's update simplifies to $x_{k+1}=(x_k+2/x_k)/2$. Starting from the upper endpoint gives:

| Iteration | Newton value | Absolute error |
|---|---:|---:|
| 0 | $2.000000000000$ | $5.857864e-01$ |
| 1 | $1.500000000000$ | $8.578644e-02$ |
| 2 | $1.416666666667$ | $2.453104e-03$ |
| 3 | $1.414215686275$ | $2.123901e-06$ |
| 4 | $1.414213562375$ | $1.594724e-12$ |

Let $r=\sqrt{2}$. Direct subtraction, using $r^2=2$, gives

$$
x_{k+1}-r
=\frac{(x_k-r)^2}{2x_k}.
$$

Positive iterates remain on or above the root after one step. An iterate above the root decreases unless already equal to it. The identity then proves convergence and quadratic error reduction near the root, since the denominator approaches $2r>0$.

A general Newton iteration does not inherit these properties automatically. A derivative close to zero or a poor initial point can cause a large step, a cycle, or departure from the domain. A safeguarded method can retain a sign-changing bracket and accept a Newton step only when it remains suitable.

#### Problem NM-58: computing the square root of 3

Find the positive root of $f(x)=x^2-3$. Give a bisection iteration count ensuring midpoint error at most $10^{-6}$ and show Newton iterates starting from $x_0=2$.

Solution. The function is continuous and strictly increasing for positive $x$. Its values at $1$ and $2$ have opposite signs, so this interval contains exactly one positive root.

After $n$ bisection reductions, the bracket width is $(2-1)/2^n$. The midpoint error is at most half that width. Thus it is enough to choose

$$
n\ge\log_2\frac{1}{2\cdot10^{-6}},
\qquad n=19.
$$

This count depends on the initial bracket and requested absolute error, not on the derivative at the root.

Newton's update simplifies to $x_{k+1}=(x_k+3/x_k)/2$. Starting from the upper endpoint gives:

| Iteration | Newton value | Absolute error |
|---|---:|---:|
| 0 | $2.000000000000$ | $2.679492e-01$ |
| 1 | $1.750000000000$ | $1.794919e-02$ |
| 2 | $1.732142857143$ | $9.204957e-05$ |
| 3 | $1.732050810015$ | $2.445850e-09$ |
| 4 | $1.732050807569$ | $0.000000e+00$ |

Let $r=\sqrt{3}$. Direct subtraction, using $r^2=3$, gives

$$
x_{k+1}-r
=\frac{(x_k-r)^2}{2x_k}.
$$

Positive iterates remain on or above the root after one step. An iterate above the root decreases unless already equal to it. The identity then proves convergence and quadratic error reduction near the root, since the denominator approaches $2r>0$.

A general Newton iteration does not inherit these properties automatically. A derivative close to zero or a poor initial point can cause a large step, a cycle, or departure from the domain. A safeguarded method can retain a sign-changing bracket and accept a Newton step only when it remains suitable.

#### Problem NM-59: computing the square root of 5

Find the positive root of $f(x)=x^2-5$. Give a bisection iteration count ensuring midpoint error at most $10^{-6}$ and show Newton iterates starting from $x_0=3$.

Solution. The function is continuous and strictly increasing for positive $x$. Its values at $2$ and $3$ have opposite signs, so this interval contains exactly one positive root.

After $n$ bisection reductions, the bracket width is $(3-2)/2^n$. The midpoint error is at most half that width. Thus it is enough to choose

$$
n\ge\log_2\frac{1}{2\cdot10^{-6}},
\qquad n=19.
$$

This count depends on the initial bracket and requested absolute error, not on the derivative at the root.

Newton's update simplifies to $x_{k+1}=(x_k+5/x_k)/2$. Starting from the upper endpoint gives:

| Iteration | Newton value | Absolute error |
|---|---:|---:|
| 0 | $3.000000000000$ | $7.639320e-01$ |
| 1 | $2.333333333333$ | $9.726536e-02$ |
| 2 | $2.238095238095$ | $2.027261e-03$ |
| 3 | $2.236068895643$ | $9.181436e-07$ |
| 4 | $2.236067977500$ | $1.882938e-13$ |

Let $r=\sqrt{5}$. Direct subtraction, using $r^2=5$, gives

$$
x_{k+1}-r
=\frac{(x_k-r)^2}{2x_k}.
$$

Positive iterates remain on or above the root after one step. An iterate above the root decreases unless already equal to it. The identity then proves convergence and quadratic error reduction near the root, since the denominator approaches $2r>0$.

A general Newton iteration does not inherit these properties automatically. A derivative close to zero or a poor initial point can cause a large step, a cycle, or departure from the domain. A safeguarded method can retain a sign-changing bracket and accept a Newton step only when it remains suitable.

#### Problem NM-60: computing the square root of 10

Find the positive root of $f(x)=x^2-10$. Give a bisection iteration count ensuring midpoint error at most $10^{-6}$ and show Newton iterates starting from $x_0=4$.

Solution. The function is continuous and strictly increasing for positive $x$. Its values at $3$ and $4$ have opposite signs, so this interval contains exactly one positive root.

After $n$ bisection reductions, the bracket width is $(4-3)/2^n$. The midpoint error is at most half that width. Thus it is enough to choose

$$
n\ge\log_2\frac{1}{2\cdot10^{-6}},
\qquad n=19.
$$

This count depends on the initial bracket and requested absolute error, not on the derivative at the root.

Newton's update simplifies to $x_{k+1}=(x_k+10/x_k)/2$. Starting from the upper endpoint gives:

| Iteration | Newton value | Absolute error |
|---|---:|---:|
| 0 | $4.000000000000$ | $8.377223e-01$ |
| 1 | $3.250000000000$ | $8.772234e-02$ |
| 2 | $3.163461538462$ | $1.183878e-03$ |
| 3 | $3.162277881693$ | $2.215244e-07$ |
| 4 | $3.162277660168$ | $7.549517e-15$ |

Let $r=\sqrt{10}$. Direct subtraction, using $r^2=10$, gives

$$
x_{k+1}-r
=\frac{(x_k-r)^2}{2x_k}.
$$

Positive iterates remain on or above the root after one step. An iterate above the root decreases unless already equal to it. The identity then proves convergence and quadratic error reduction near the root, since the denominator approaches $2r>0$.

A general Newton iteration does not inherit these properties automatically. A derivative close to zero or a poor initial point can cause a large step, a cycle, or departure from the domain. A safeguarded method can retain a sign-changing bracket and accept a Newton step only when it remains suitable.

### 8.4 Gradient-descent stability along eigenvectors

On a positive-definite quadratic, a fixed gradient step acts independently on each eigenvector direction. This reduces convergence to scalar multiplier conditions. A step that improves one direction can still make the high-curvature direction oscillate or diverge.

#### Problem NM-61: fixed-step quadratic descent, case 1

Minimize $f(x)=\tfrac12(1x_1^2+4x_2^2)$ by gradient descent from $x_0=(1,1)^\top$ with step $\eta=\frac{1}{4}$. Derive the exact iterates and classify convergence.

Solution. The gradient is $(1x_1,4x_2)^\top$. The update is diagonal:

$$
x_{k+1}=
\begin{pmatrix}1-\eta1&0\\0&1-\eta4\end{pmatrix}x_k.
$$

For the supplied step, the two multipliers are $m_1=\frac{3}{4}$ and $m_2=0$. Since both initial coordinates equal one,

$$
x_k=\begin{pmatrix}(\frac{3}{4})^k\\(0)^k\end{pmatrix}.
$$

Both multipliers have magnitude below one, so every starting vector converges to zero.

The exact objective along the iterates is

$$
f(x_k)=\frac12\left[1(\frac{3}{4})^{2k}
+4(0)^{2k}\right].
$$

For convergence from every starting vector, both multiplier magnitudes must be less than one. Since the largest eigenvalue is $4$, the condition is $0<\eta<2/4$. Equality at the upper endpoint gives a persistent oscillation in that eigenvector direction.

For this two-eigenvalue quadratic, the step minimizing the worst absolute multiplier is $\eta_*=2/(\lambda_{\min}+\lambda_{\max})$, with contraction factor $(\kappa-1)/(\kappa+1)$. A high condition number makes that best factor close to one. Preconditioning changes the relative curvature before applying descent; reducing the step alone does not remove the slow weak-curvature direction.

This exact calculation is stronger than observing a few decreasing objective values. It identifies the behavior for every iteration and explains the strict step-size restriction.

#### Problem NM-62: fixed-step quadratic descent, case 2

Minimize $f(x)=\tfrac12(2x_1^2+5x_2^2)$ by gradient descent from $x_0=(1,1)^\top$ with step $\eta=\frac{1}{5}$. Derive the exact iterates and classify convergence.

Solution. The gradient is $(2x_1,5x_2)^\top$. The update is diagonal:

$$
x_{k+1}=
\begin{pmatrix}1-\eta2&0\\0&1-\eta5\end{pmatrix}x_k.
$$

For the supplied step, the two multipliers are $m_1=\frac{3}{5}$ and $m_2=0$. Since both initial coordinates equal one,

$$
x_k=\begin{pmatrix}(\frac{3}{5})^k\\(0)^k\end{pmatrix}.
$$

Both multipliers have magnitude below one, so every starting vector converges to zero.

The exact objective along the iterates is

$$
f(x_k)=\frac12\left[2(\frac{3}{5})^{2k}
+5(0)^{2k}\right].
$$

For convergence from every starting vector, both multiplier magnitudes must be less than one. Since the largest eigenvalue is $5$, the condition is $0<\eta<2/5$. Equality at the upper endpoint gives a persistent oscillation in that eigenvector direction.

For this two-eigenvalue quadratic, the step minimizing the worst absolute multiplier is $\eta_*=2/(\lambda_{\min}+\lambda_{\max})$, with contraction factor $(\kappa-1)/(\kappa+1)$. A high condition number makes that best factor close to one. Preconditioning changes the relative curvature before applying descent; reducing the step alone does not remove the slow weak-curvature direction.

This exact calculation is stronger than observing a few decreasing objective values. It identifies the behavior for every iteration and explains the strict step-size restriction.

#### Problem NM-63: fixed-step quadratic descent, case 3

Minimize $f(x)=\tfrac12(1x_1^2+3x_2^2)$ by gradient descent from $x_0=(1,1)^\top$ with step $\eta=\frac{2}{3}$. Derive the exact iterates and classify convergence.

Solution. The gradient is $(1x_1,3x_2)^\top$. The update is diagonal:

$$
x_{k+1}=
\begin{pmatrix}1-\eta1&0\\0&1-\eta3\end{pmatrix}x_k.
$$

For the supplied step, the two multipliers are $m_1=\frac{1}{3}$ and $m_2=-1$. Since both initial coordinates equal one,

$$
x_k=\begin{pmatrix}(\frac{1}{3})^k\\(-1)^k\end{pmatrix}.
$$

The second multiplier equals minus one. With a nonzero second initial coordinate, that coordinate alternates without decaying, so the iterates do not converge to the minimum.

The exact objective along the iterates is

$$
f(x_k)=\frac12\left[1(\frac{1}{3})^{2k}
+3(-1)^{2k}\right].
$$

For convergence from every starting vector, both multiplier magnitudes must be less than one. Since the largest eigenvalue is $3$, the condition is $0<\eta<2/3$. Equality at the upper endpoint gives a persistent oscillation in that eigenvector direction.

For this two-eigenvalue quadratic, the step minimizing the worst absolute multiplier is $\eta_*=2/(\lambda_{\min}+\lambda_{\max})$, with contraction factor $(\kappa-1)/(\kappa+1)$. A high condition number makes that best factor close to one. Preconditioning changes the relative curvature before applying descent; reducing the step alone does not remove the slow weak-curvature direction.

This exact calculation is stronger than observing a few decreasing objective values. It identifies the behavior for every iteration and explains the strict step-size restriction.

#### Problem NM-64: fixed-step quadratic descent, case 4

Minimize $f(x)=\tfrac12(1x_1^2+4x_2^2)$ by gradient descent from $x_0=(1,1)^\top$ with step $\eta=\frac{3}{5}$. Derive the exact iterates and classify convergence.

Solution. The gradient is $(1x_1,4x_2)^\top$. The update is diagonal:

$$
x_{k+1}=
\begin{pmatrix}1-\eta1&0\\0&1-\eta4\end{pmatrix}x_k.
$$

For the supplied step, the two multipliers are $m_1=\frac{2}{5}$ and $m_2=-\frac{7}{5}$. Since both initial coordinates equal one,

$$
x_k=\begin{pmatrix}(\frac{2}{5})^k\\(-\frac{7}{5})^k\end{pmatrix}.
$$

The second multiplier has magnitude above one. From the supplied nonzero second coordinate, its magnitude grows geometrically and the method diverges.

The exact objective along the iterates is

$$
f(x_k)=\frac12\left[1(\frac{2}{5})^{2k}
+4(-\frac{7}{5})^{2k}\right].
$$

For convergence from every starting vector, both multiplier magnitudes must be less than one. Since the largest eigenvalue is $4$, the condition is $0<\eta<2/4$. Equality at the upper endpoint gives a persistent oscillation in that eigenvector direction.

For this two-eigenvalue quadratic, the step minimizing the worst absolute multiplier is $\eta_*=2/(\lambda_{\min}+\lambda_{\max})$, with contraction factor $(\kappa-1)/(\kappa+1)$. A high condition number makes that best factor close to one. Preconditioning changes the relative curvature before applying descent; reducing the step alone does not remove the slow weak-curvature direction.

This exact calculation is stronger than observing a few decreasing objective values. It identifies the behavior for every iteration and explains the strict step-size restriction.

### 8.5 Finite differences and removable truncation error

A polynomial permits an exact expansion of finite-difference error. Forward and centered differences use different cancellations. Richardson extrapolation removes a known leading power of the step size; it does not remove floating-point error or an incorrectly assumed error expansion.

#### Problem NM-65: differentiating a cubic, case 1

For $f(x)=x^3$ at $x=2$ with $h=\frac{1}{10}$, compute forward and centered derivative approximations, their errors, and a Richardson extrapolation.

Solution. Expand $(x+h)^3=x^3+3x^2h+3xh^2+h^3$. The forward difference is therefore

$$
D_+(h)=\frac{f(x+h)-f(x)}h
=3x^2+3xh+h^2=\frac{1261}{100}.
$$

The exact derivative is $f'(x)=3x^2=12$, so the signed forward error is $3xh+h^2=\frac{61}{100}$. At a special point such as zero, the first-order coefficient can vanish; an observed higher rate there does not make the method second order for general functions and points.

Subtract the expansions at $x+h$ and $x-h$. Even powers of $h$ cancel before division:

$$
D_0(h)=\frac{f(x+h)-f(x-h)}{2h}
=3x^2+h^2=\frac{1201}{100}.
$$

Its signed error is $h^2=\frac{1}{100}$. At half the step, $D_0(h/2)=3x^2+h^2/4=\frac{4801}{400}$. Thus

$$
\frac{4D_0(h/2)-D_0(h)}3
=3x^2=12.
$$

The extrapolated derivative is exact for this cubic in exact arithmetic because there are no remaining higher-order terms. For a general sufficiently smooth function, the same combination cancels an $h^2$ term and typically leaves an $O(h^4)$ error.

Every difference divides a subtraction error by $h$. Consequently, taking a smaller step eventually increases floating-point error even while reducing truncation error. The algebraic cancellation in Richardson extrapolation also combines rounded quantities and can magnify their noise. Its accuracy depends on being in a range where the assumed smooth truncation expansion still dominates arithmetic error.

#### Problem NM-66: differentiating a cubic, case 2

For $f(x)=x^3$ at $x=-1$ with $h=\frac{1}{20}$, compute forward and centered derivative approximations, their errors, and a Richardson extrapolation.

Solution. Expand $(x+h)^3=x^3+3x^2h+3xh^2+h^3$. The forward difference is therefore

$$
D_+(h)=\frac{f(x+h)-f(x)}h
=3x^2+3xh+h^2=\frac{1141}{400}.
$$

The exact derivative is $f'(x)=3x^2=3$, so the signed forward error is $3xh+h^2=-\frac{59}{400}$. At a special point such as zero, the first-order coefficient can vanish; an observed higher rate there does not make the method second order for general functions and points.

Subtract the expansions at $x+h$ and $x-h$. Even powers of $h$ cancel before division:

$$
D_0(h)=\frac{f(x+h)-f(x-h)}{2h}
=3x^2+h^2=\frac{1201}{400}.
$$

Its signed error is $h^2=\frac{1}{400}$. At half the step, $D_0(h/2)=3x^2+h^2/4=\frac{4801}{1600}$. Thus

$$
\frac{4D_0(h/2)-D_0(h)}3
=3x^2=3.
$$

The extrapolated derivative is exact for this cubic in exact arithmetic because there are no remaining higher-order terms. For a general sufficiently smooth function, the same combination cancels an $h^2$ term and typically leaves an $O(h^4)$ error.

Every difference divides a subtraction error by $h$. Consequently, taking a smaller step eventually increases floating-point error even while reducing truncation error. The algebraic cancellation in Richardson extrapolation also combines rounded quantities and can magnify their noise. Its accuracy depends on being in a range where the assumed smooth truncation expansion still dominates arithmetic error.

#### Problem NM-67: differentiating a cubic, case 3

For $f(x)=x^3$ at $x=0$ with $h=\frac{1}{100}$, compute forward and centered derivative approximations, their errors, and a Richardson extrapolation.

Solution. Expand $(x+h)^3=x^3+3x^2h+3xh^2+h^3$. The forward difference is therefore

$$
D_+(h)=\frac{f(x+h)-f(x)}h
=3x^2+3xh+h^2=\frac{1}{10000}.
$$

The exact derivative is $f'(x)=3x^2=0$, so the signed forward error is $3xh+h^2=\frac{1}{10000}$. At a special point such as zero, the first-order coefficient can vanish; an observed higher rate there does not make the method second order for general functions and points.

Subtract the expansions at $x+h$ and $x-h$. Even powers of $h$ cancel before division:

$$
D_0(h)=\frac{f(x+h)-f(x-h)}{2h}
=3x^2+h^2=\frac{1}{10000}.
$$

Its signed error is $h^2=\frac{1}{10000}$. At half the step, $D_0(h/2)=3x^2+h^2/4=\frac{1}{40000}$. Thus

$$
\frac{4D_0(h/2)-D_0(h)}3
=3x^2=0.
$$

The extrapolated derivative is exact for this cubic in exact arithmetic because there are no remaining higher-order terms. For a general sufficiently smooth function, the same combination cancels an $h^2$ term and typically leaves an $O(h^4)$ error.

Every difference divides a subtraction error by $h$. Consequently, taking a smaller step eventually increases floating-point error even while reducing truncation error. The algebraic cancellation in Richardson extrapolation also combines rounded quantities and can magnify their noise. Its accuracy depends on being in a range where the assumed smooth truncation expansion still dominates arithmetic error.

#### Problem NM-68: differentiating a cubic, case 4

For $f(x)=x^3$ at $x=3$ with $h=\frac{1}{5}$, compute forward and centered derivative approximations, their errors, and a Richardson extrapolation.

Solution. Expand $(x+h)^3=x^3+3x^2h+3xh^2+h^3$. The forward difference is therefore

$$
D_+(h)=\frac{f(x+h)-f(x)}h
=3x^2+3xh+h^2=\frac{721}{25}.
$$

The exact derivative is $f'(x)=3x^2=27$, so the signed forward error is $3xh+h^2=\frac{46}{25}$. At a special point such as zero, the first-order coefficient can vanish; an observed higher rate there does not make the method second order for general functions and points.

Subtract the expansions at $x+h$ and $x-h$. Even powers of $h$ cancel before division:

$$
D_0(h)=\frac{f(x+h)-f(x-h)}{2h}
=3x^2+h^2=\frac{676}{25}.
$$

Its signed error is $h^2=\frac{1}{25}$. At half the step, $D_0(h/2)=3x^2+h^2/4=\frac{2701}{100}$. Thus

$$
\frac{4D_0(h/2)-D_0(h)}3
=3x^2=27.
$$

The extrapolated derivative is exact for this cubic in exact arithmetic because there are no remaining higher-order terms. For a general sufficiently smooth function, the same combination cancels an $h^2$ term and typically leaves an $O(h^4)$ error.

Every difference divides a subtraction error by $h$. Consequently, taking a smaller step eventually increases floating-point error even while reducing truncation error. The algebraic cancellation in Richardson extrapolation also combines rounded quantities and can magnify their noise. Its accuracy depends on being in a range where the assumed smooth truncation expansion still dominates arithmetic error.

### 8.6 Time-stepping amplification

For $y'=\lambda y$, a one-step method multiplies the current value by an amplification factor. Comparing that factor with $e^{h\lambda}$ separates accuracy from stability. Negative real $\lambda$ gives a decaying exact solution, so a growing numerical factor is immediately suspect.

#### Problem NM-69: three methods on a decaying equation, case 1

Solve $y'=-2y$, $y(0)=1$, using step $h=\frac{1}{4}$ for $N=4$ steps. Compare explicit Euler, backward Euler, and explicit trapezoidal (Heun) values and stability.

Solution. The exact terminal value is $e^{\lambda Nh}=e^{-2}\approx 1.35335283e-01$. Put $z=h\lambda=-\frac{1}{2}$.

Explicit Euler gives $y_{n+1}=(1+z)y_n$. Backward Euler solves $y_{n+1}=y_n+zy_{n+1}$, so its factor is $1/(1-z)$. Heun predicts $\tilde y=y_n+zy_n$ and averages the two slopes, giving factor $1+z+z^2/2$.

| Method | One-step factor | Value after $N$ steps |
|---|---:|---:|
| Explicit Euler | $\frac{1}{2}$ | $6.25000000e-02$ |
| Backward Euler | $\frac{2}{3}$ | $1.97530864e-01$ |
| Heun | $\frac{5}{8}$ | $1.52587891e-01$ |

On the negative real axis, explicit Euler is stable when $-2\le z\le0$. Heun has the same real-axis interval even though its local accuracy is higher. Backward Euler has magnitude below one for every negative real $z$.

Stability does not by itself guarantee an accurate value at this fixed step. A factor of zero can erase the numerical solution in one step even though the exact solution is positive. A factor of minus one gives a nondecaying sign alternation. A large stable implicit step may suppress or distort a transient that matters to the application.

Under refinement with fixed final time, explicit and backward Euler converge with first-order global error, while Heun has second-order error for this smooth problem. Stiffness concerns the practical step restriction needed to control rapidly decaying modes, which can be much smaller than the time scale on which the desired solution changes.

#### Problem NM-70: three methods on a decaying equation, case 2

Solve $y'=-10y$, $y(0)=1$, using step $h=\frac{1}{10}$ for $N=5$ steps. Compare explicit Euler, backward Euler, and explicit trapezoidal (Heun) values and stability.

Solution. The exact terminal value is $e^{\lambda Nh}=e^{-5}\approx 6.73794700e-03$. Put $z=h\lambda=-1$.

Explicit Euler gives $y_{n+1}=(1+z)y_n$. Backward Euler solves $y_{n+1}=y_n+zy_{n+1}$, so its factor is $1/(1-z)$. Heun predicts $\tilde y=y_n+zy_n$ and averages the two slopes, giving factor $1+z+z^2/2$.

| Method | One-step factor | Value after $N$ steps |
|---|---:|---:|
| Explicit Euler | $0$ | $0.00000000e+00$ |
| Backward Euler | $\frac{1}{2}$ | $3.12500000e-02$ |
| Heun | $\frac{1}{2}$ | $3.12500000e-02$ |

On the negative real axis, explicit Euler is stable when $-2\le z\le0$. Heun has the same real-axis interval even though its local accuracy is higher. Backward Euler has magnitude below one for every negative real $z$.

Stability does not by itself guarantee an accurate value at this fixed step. A factor of zero can erase the numerical solution in one step even though the exact solution is positive. A factor of minus one gives a nondecaying sign alternation. A large stable implicit step may suppress or distort a transient that matters to the application.

Under refinement with fixed final time, explicit and backward Euler converge with first-order global error, while Heun has second-order error for this smooth problem. Stiffness concerns the practical step restriction needed to control rapidly decaying modes, which can be much smaller than the time scale on which the desired solution changes.

#### Problem NM-71: three methods on a decaying equation, case 3

Solve $y'=-10y$, $y(0)=1$, using step $h=\frac{1}{4}$ for $N=4$ steps. Compare explicit Euler, backward Euler, and explicit trapezoidal (Heun) values and stability.

Solution. The exact terminal value is $e^{\lambda Nh}=e^{-10}\approx 4.53999298e-05$. Put $z=h\lambda=-\frac{5}{2}$.

Explicit Euler gives $y_{n+1}=(1+z)y_n$. Backward Euler solves $y_{n+1}=y_n+zy_{n+1}$, so its factor is $1/(1-z)$. Heun predicts $\tilde y=y_n+zy_n$ and averages the two slopes, giving factor $1+z+z^2/2$.

| Method | One-step factor | Value after $N$ steps |
|---|---:|---:|
| Explicit Euler | $-\frac{3}{2}$ | $5.06250000e+00$ |
| Backward Euler | $\frac{2}{7}$ | $6.66389005e-03$ |
| Heun | $\frac{13}{8}$ | $6.97290039e+00$ |

On the negative real axis, explicit Euler is stable when $-2\le z\le0$. Heun has the same real-axis interval even though its local accuracy is higher. Backward Euler has magnitude below one for every negative real $z$.

Stability does not by itself guarantee an accurate value at this fixed step. A factor of zero can erase the numerical solution in one step even though the exact solution is positive. A factor of minus one gives a nondecaying sign alternation. A large stable implicit step may suppress or distort a transient that matters to the application.

Under refinement with fixed final time, explicit and backward Euler converge with first-order global error, while Heun has second-order error for this smooth problem. Stiffness concerns the practical step restriction needed to control rapidly decaying modes, which can be much smaller than the time scale on which the desired solution changes.

#### Problem NM-72: three methods on a decaying equation, case 4

Solve $y'=-100y$, $y(0)=1$, using step $h=\frac{1}{10}$ for $N=10$ steps. Compare explicit Euler, backward Euler, and explicit trapezoidal (Heun) values and stability.

Solution. The exact terminal value is $e^{\lambda Nh}=e^{-100}\approx 3.72007598e-44$. Put $z=h\lambda=-10$.

Explicit Euler gives $y_{n+1}=(1+z)y_n$. Backward Euler solves $y_{n+1}=y_n+zy_{n+1}$, so its factor is $1/(1-z)$. Heun predicts $\tilde y=y_n+zy_n$ and averages the two slopes, giving factor $1+z+z^2/2$.

| Method | One-step factor | Value after $N$ steps |
|---|---:|---:|
| Explicit Euler | $-9$ | $3.48678440e+09$ |
| Backward Euler | $\frac{1}{11}$ | $3.85543289e-11$ |
| Heun | $41$ | $1.34226593e+16$ |

On the negative real axis, explicit Euler is stable when $-2\le z\le0$. Heun has the same real-axis interval even though its local accuracy is higher. Backward Euler has magnitude below one for every negative real $z$.

Stability does not by itself guarantee an accurate value at this fixed step. A factor of zero can erase the numerical solution in one step even though the exact solution is positive. A factor of minus one gives a nondecaying sign alternation. A large stable implicit step may suppress or distort a transient that matters to the application.

Under refinement with fixed final time, explicit and backward Euler converge with first-order global error, while Heun has second-order error for this smooth problem. Stiffness concerns the practical step restriction needed to control rapidly decaying modes, which can be much smaller than the time scale on which the desired solution changes.

### 8.7 A discrete heat mode and the mesh restriction

The highest-frequency mesh modes are often the first to become unstable in an explicit diffusion scheme. A three-point interior grid is small enough to calculate an eigenmode exactly and still display the stability mechanism.

#### Problem NM-73: explicit heat update with mesh ratio \frac{1}{4}

For zero-boundary heat diffusion, take three interior points and $r=\nu k/h^2=\frac{1}{4}$. Start from $u^0=(1/\sqrt2,-1,1/\sqrt2)^\top$. Find three explicit time steps and analyze the mode.

Solution. The interior second-difference matrix without the $h^{-2}$ factor is

$$
L=\begin{pmatrix}-2&1&0\\1&-2&1\\0&1&-2\end{pmatrix}.
$$

For the supplied alternating vector $v$, direct multiplication gives $Lv=(-2-\sqrt2)v$. The explicit update $u^{n+1}=(I+rL)u^n$ therefore preserves this direction and changes only its magnitude:

$$
u^n=g^nv,\qquad
g=1-r(2+\sqrt2)\approx 0.146447.
$$

After three steps,

$$
u^3\approx
\begin{pmatrix}0.002221\\-0.003141\\0.002221\end{pmatrix}.
$$

This mode decays in magnitude precisely when $|g|<1$; a negative factor makes its sign alternate between time levels. At $|g|>1$, this numerical high-frequency component grows even though diffusion should dissipate it.

On this particular three-interior-point grid, the most negative eigenvalue gives stability limit $r\le2/(2+\sqrt2)\approx0.585786$. On increasingly fine grids, the most negative scaled eigenvalue approaches minus four, so the familiar uniform-grid limit approaches $r\le1/2$.

The stronger maximum-principle condition $0\le r\le1/2$ follows directly from the stencil weights $r,1-2r,r$ being nonnegative and summing to one. It is a convenient grid-independent sufficient condition and ensures that one update cannot create a new interior extremum beyond the supplied values and boundaries.

Backward Euler solves $(I-rL)u^{n+1}=u^n$ and has factor $1/[1+r(2+\sqrt2)]$ on this mode, below one for every positive $r$. Its unconditional linear stability still leaves temporal accuracy to be controlled by the step size.

#### Problem NM-74: explicit heat update with mesh ratio \frac{1}{2}

For zero-boundary heat diffusion, take three interior points and $r=\nu k/h^2=\frac{1}{2}$. Start from $u^0=(1/\sqrt2,-1,1/\sqrt2)^\top$. Find three explicit time steps and analyze the mode.

Solution. The interior second-difference matrix without the $h^{-2}$ factor is

$$
L=\begin{pmatrix}-2&1&0\\1&-2&1\\0&1&-2\end{pmatrix}.
$$

For the supplied alternating vector $v$, direct multiplication gives $Lv=(-2-\sqrt2)v$. The explicit update $u^{n+1}=(I+rL)u^n$ therefore preserves this direction and changes only its magnitude:

$$
u^n=g^nv,\qquad
g=1-r(2+\sqrt2)\approx -0.707107.
$$

After three steps,

$$
u^3\approx
\begin{pmatrix}-0.250000\\0.353553\\-0.250000\end{pmatrix}.
$$

This mode decays in magnitude precisely when $|g|<1$; a negative factor makes its sign alternate between time levels. At $|g|>1$, this numerical high-frequency component grows even though diffusion should dissipate it.

On this particular three-interior-point grid, the most negative eigenvalue gives stability limit $r\le2/(2+\sqrt2)\approx0.585786$. On increasingly fine grids, the most negative scaled eigenvalue approaches minus four, so the familiar uniform-grid limit approaches $r\le1/2$.

The stronger maximum-principle condition $0\le r\le1/2$ follows directly from the stencil weights $r,1-2r,r$ being nonnegative and summing to one. It is a convenient grid-independent sufficient condition and ensures that one update cannot create a new interior extremum beyond the supplied values and boundaries.

Backward Euler solves $(I-rL)u^{n+1}=u^n$ and has factor $1/[1+r(2+\sqrt2)]$ on this mode, below one for every positive $r$. Its unconditional linear stability still leaves temporal accuracy to be controlled by the step size.

#### Problem NM-75: explicit heat update with mesh ratio \frac{3}{5}

For zero-boundary heat diffusion, take three interior points and $r=\nu k/h^2=\frac{3}{5}$. Start from $u^0=(1/\sqrt2,-1,1/\sqrt2)^\top$. Find three explicit time steps and analyze the mode.

Solution. The interior second-difference matrix without the $h^{-2}$ factor is

$$
L=\begin{pmatrix}-2&1&0\\1&-2&1\\0&1&-2\end{pmatrix}.
$$

For the supplied alternating vector $v$, direct multiplication gives $Lv=(-2-\sqrt2)v$. The explicit update $u^{n+1}=(I+rL)u^n$ therefore preserves this direction and changes only its magnitude:

$$
u^n=g^nv,\qquad
g=1-r(2+\sqrt2)\approx -1.048528.
$$

After three steps,

$$
u^3\approx
\begin{pmatrix}-0.815127\\1.152764\\-0.815127\end{pmatrix}.
$$

This mode decays in magnitude precisely when $|g|<1$; a negative factor makes its sign alternate between time levels. At $|g|>1$, this numerical high-frequency component grows even though diffusion should dissipate it.

On this particular three-interior-point grid, the most negative eigenvalue gives stability limit $r\le2/(2+\sqrt2)\approx0.585786$. On increasingly fine grids, the most negative scaled eigenvalue approaches minus four, so the familiar uniform-grid limit approaches $r\le1/2$.

The stronger maximum-principle condition $0\le r\le1/2$ follows directly from the stencil weights $r,1-2r,r$ being nonnegative and summing to one. It is a convenient grid-independent sufficient condition and ensures that one update cannot create a new interior extremum beyond the supplied values and boundaries.

Backward Euler solves $(I-rL)u^{n+1}=u^n$ and has factor $1/[1+r(2+\sqrt2)]$ on this mode, below one for every positive $r$. Its unconditional linear stability still leaves temporal accuracy to be controlled by the step size.

#### Problem NM-76: explicit heat update with mesh ratio 1

For zero-boundary heat diffusion, take three interior points and $r=\nu k/h^2=1$. Start from $u^0=(1/\sqrt2,-1,1/\sqrt2)^\top$. Find three explicit time steps and analyze the mode.

Solution. The interior second-difference matrix without the $h^{-2}$ factor is

$$
L=\begin{pmatrix}-2&1&0\\1&-2&1\\0&1&-2\end{pmatrix}.
$$

For the supplied alternating vector $v$, direct multiplication gives $Lv=(-2-\sqrt2)v$. The explicit update $u^{n+1}=(I+rL)u^n$ therefore preserves this direction and changes only its magnitude:

$$
u^n=g^nv,\qquad
g=1-r(2+\sqrt2)\approx -2.414214.
$$

After three steps,

$$
u^3\approx
\begin{pmatrix}-9.949747\\14.071068\\-9.949747\end{pmatrix}.
$$

This mode decays in magnitude precisely when $|g|<1$; a negative factor makes its sign alternate between time levels. At $|g|>1$, this numerical high-frequency component grows even though diffusion should dissipate it.

On this particular three-interior-point grid, the most negative eigenvalue gives stability limit $r\le2/(2+\sqrt2)\approx0.585786$. On increasingly fine grids, the most negative scaled eigenvalue approaches minus four, so the familiar uniform-grid limit approaches $r\le1/2$.

The stronger maximum-principle condition $0\le r\le1/2$ follows directly from the stencil weights $r,1-2r,r$ being nonnegative and summing to one. It is a convenient grid-independent sufficient condition and ensures that one update cannot create a new interior extremum beyond the supplied values and boundaries.

Backward Euler solves $(I-rL)u^{n+1}=u^n$ and has factor $1/[1+r(2+\sqrt2)]$ on this mode, below one for every positive $r$. Its unconditional linear stability still leaves temporal accuracy to be controlled by the step size.

### 8.8 A polynomial control variate with known moments

Monte Carlo uncertainty can be reduced by subtracting a correlated quantity with known expectation. An optimal coefficient follows by minimizing a quadratic variance expression. A fair comparison records how many random inputs and function evaluations are used.

#### Problem NM-77: integrating a uniform power 2

Estimate $I=\int_0^1u^{2}\,du$ using independent uniform draws. Use $U-1/2$ as a control variate. Derive the optimal coefficient and exact variance reduction.

Solution. For $f(U)=U^p$ with $p=2$, uniform moments give

$$
I=E[U^p]=\frac1{p+1}=\frac{1}{3},\qquad
\operatorname{Var}(U^p)=\frac1{2p+1}-\frac1{(p+1)^2}=\frac{4}{45}.
$$

Define one controlled observation $Y_c=U^p-c(U-1/2)$. Its expectation remains $I$ for every fixed $c$ because the control has expectation zero.

The required covariance and control variance are

$$
\operatorname{Cov}(U^p,U)
=\frac1{p+2}-\frac1{2(p+1)}=\frac{1}{12},\qquad
\operatorname{Var}(U)=\frac1{12}.
$$

Expanding the controlled variance gives $\operatorname{Var}(Y_c)=\operatorname{Var}(U^p)-2c\operatorname{Cov}(U^p,U)+c^2/12$. Its unique minimum occurs at

$$
c_*=\frac{\operatorname{Cov}(U^p,U)}{\operatorname{Var}(U)}=1.
$$

The resulting one-draw variance is $\frac{1}{180}$, so a mean of $N$ independent controlled observations has variance $\frac{1}{180}/N$. The plain Monte Carlo variance is $\frac{4}{45}/N$. Their ratio is $16\approx 16.000000$.

Both estimators use the same number of independent uniform draws and one evaluation of $U^p$ per draw. The control only adds inexpensive arithmetic here. The variance calculation assumes the coefficient is fixed at its known optimal value. Estimating a coefficient on a pilot sample, or estimating it jointly on the main sample, requires accounting for that procedure. A known expectation for the control is essential; subtracting a correlated variable with an unknown mean can introduce bias.

#### Problem NM-78: integrating a uniform power 3

Estimate $I=\int_0^1u^{3}\,du$ using independent uniform draws. Use $U-1/2$ as a control variate. Derive the optimal coefficient and exact variance reduction.

Solution. For $f(U)=U^p$ with $p=3$, uniform moments give

$$
I=E[U^p]=\frac1{p+1}=\frac{1}{4},\qquad
\operatorname{Var}(U^p)=\frac1{2p+1}-\frac1{(p+1)^2}=\frac{9}{112}.
$$

Define one controlled observation $Y_c=U^p-c(U-1/2)$. Its expectation remains $I$ for every fixed $c$ because the control has expectation zero.

The required covariance and control variance are

$$
\operatorname{Cov}(U^p,U)
=\frac1{p+2}-\frac1{2(p+1)}=\frac{3}{40},\qquad
\operatorname{Var}(U)=\frac1{12}.
$$

Expanding the controlled variance gives $\operatorname{Var}(Y_c)=\operatorname{Var}(U^p)-2c\operatorname{Cov}(U^p,U)+c^2/12$. Its unique minimum occurs at

$$
c_*=\frac{\operatorname{Cov}(U^p,U)}{\operatorname{Var}(U)}=\frac{9}{10}.
$$

The resulting one-draw variance is $\frac{9}{700}$, so a mean of $N$ independent controlled observations has variance $\frac{9}{700}/N$. The plain Monte Carlo variance is $\frac{9}{112}/N$. Their ratio is $\frac{25}{4}\approx 6.250000$.

Both estimators use the same number of independent uniform draws and one evaluation of $U^p$ per draw. The control only adds inexpensive arithmetic here. The variance calculation assumes the coefficient is fixed at its known optimal value. Estimating a coefficient on a pilot sample, or estimating it jointly on the main sample, requires accounting for that procedure. A known expectation for the control is essential; subtracting a correlated variable with an unknown mean can introduce bias.

#### Problem NM-79: integrating a uniform power 4

Estimate $I=\int_0^1u^{4}\,du$ using independent uniform draws. Use $U-1/2$ as a control variate. Derive the optimal coefficient and exact variance reduction.

Solution. For $f(U)=U^p$ with $p=4$, uniform moments give

$$
I=E[U^p]=\frac1{p+1}=\frac{1}{5},\qquad
\operatorname{Var}(U^p)=\frac1{2p+1}-\frac1{(p+1)^2}=\frac{16}{225}.
$$

Define one controlled observation $Y_c=U^p-c(U-1/2)$. Its expectation remains $I$ for every fixed $c$ because the control has expectation zero.

The required covariance and control variance are

$$
\operatorname{Cov}(U^p,U)
=\frac1{p+2}-\frac1{2(p+1)}=\frac{1}{15},\qquad
\operatorname{Var}(U)=\frac1{12}.
$$

Expanding the controlled variance gives $\operatorname{Var}(Y_c)=\operatorname{Var}(U^p)-2c\operatorname{Cov}(U^p,U)+c^2/12$. Its unique minimum occurs at

$$
c_*=\frac{\operatorname{Cov}(U^p,U)}{\operatorname{Var}(U)}=\frac{4}{5}.
$$

The resulting one-draw variance is $\frac{4}{225}$, so a mean of $N$ independent controlled observations has variance $\frac{4}{225}/N$. The plain Monte Carlo variance is $\frac{16}{225}/N$. Their ratio is $4\approx 4.000000$.

Both estimators use the same number of independent uniform draws and one evaluation of $U^p$ per draw. The control only adds inexpensive arithmetic here. The variance calculation assumes the coefficient is fixed at its known optimal value. Estimating a coefficient on a pilot sample, or estimating it jointly on the main sample, requires accounting for that procedure. A known expectation for the control is essential; subtracting a correlated variable with an unknown mean can introduce bias.

#### Problem NM-80: integrating a uniform power 5

Estimate $I=\int_0^1u^{5}\,du$ using independent uniform draws. Use $U-1/2$ as a control variate. Derive the optimal coefficient and exact variance reduction.

Solution. For $f(U)=U^p$ with $p=5$, uniform moments give

$$
I=E[U^p]=\frac1{p+1}=\frac{1}{6},\qquad
\operatorname{Var}(U^p)=\frac1{2p+1}-\frac1{(p+1)^2}=\frac{25}{396}.
$$

Define one controlled observation $Y_c=U^p-c(U-1/2)$. Its expectation remains $I$ for every fixed $c$ because the control has expectation zero.

The required covariance and control variance are

$$
\operatorname{Cov}(U^p,U)
=\frac1{p+2}-\frac1{2(p+1)}=\frac{5}{84},\qquad
\operatorname{Var}(U)=\frac1{12}.
$$

Expanding the controlled variance gives $\operatorname{Var}(Y_c)=\operatorname{Var}(U^p)-2c\operatorname{Cov}(U^p,U)+c^2/12$. Its unique minimum occurs at

$$
c_*=\frac{\operatorname{Cov}(U^p,U)}{\operatorname{Var}(U)}=\frac{5}{7}.
$$

The resulting one-draw variance is $\frac{100}{4851}$, so a mean of $N$ independent controlled observations has variance $\frac{100}{4851}/N$. The plain Monte Carlo variance is $\frac{25}{396}/N$. Their ratio is $\frac{49}{16}\approx 3.062500$.

Both estimators use the same number of independent uniform draws and one evaluation of $U^p$ per draw. The control only adds inexpensive arithmetic here. The variance calculation assumes the coefficient is fixed at its known optimal value. Estimating a coefficient on a pilot sample, or estimating it jointly on the main sample, requires accounting for that procedure. A known expectation for the control is essential; subtracting a correlated variable with an unknown mean can introduce bias.

<!-- further-distinct-problems-NM -->

### 8.9 Further accuracy and convergence problems

These additional problems are constructed for these notes.

#### Problem NM-81: residual correction as an iteration

Solve $Ax=b$ for $A=\begin{pmatrix}4&1\\1&3\end{pmatrix}$ and $b=(1,2)^\top$ using $x_{k+1}=x_k+B(b-Ax_k)$ with $B=\operatorname{diag}(1/4,1/3)$ and $x_0=0$. Compute two steps and prove convergence.

Solution. The exact inverse gives $x^*=(1/11,7/11)^\top$. The first residual is $b$, so $x_1=(1/4,2/3)^\top$. The recurrence can also be written

$$
x_{1,k+1}=\frac{1-x_{2,k}}4,\qquad
x_{2,k+1}=\frac{2-x_{1,k}}3.
$$

Thus $x_2=(1/12,7/12)^\top$. Each update uses both old coordinates, rather than substituting an already updated coordinate into the second equation.

Let $e_k=x_k-x^*$. Since $Ax^*=b$,

$$
e_{k+1}=(I-BA)e_k,\qquad
I-BA=\begin{pmatrix}0&-1/4\\-1/3&0\end{pmatrix}.
$$

The infinity norm of this matrix is $1/3<1$, so $\|e_k\|_\infty\le(1/3)^k\|e_0\|_\infty$. Its spectral radius is $1/\sqrt{12}$, consistent with convergence.

The same algebra describes a residual-correction step using an approximate inverse or approximate linear solve. If the correction used the exact inverse and the residual were exact, one step would eliminate the error. In floating-point iterative refinement, residual accuracy and correction-solve accuracy both matter. Recomputing a residual more accurately can recover information lost in the initial solve, but convergence requires control of the resulting error-propagation map.

#### Problem NM-82: a Householder reflection

Construct an orthogonal reflection taking $x=(3,4,0)^\top$ to a multiple of the first coordinate vector. Apply it to $y=(1,0,2)^\top$ and explain its first QR elimination step.

Solution. The norm of $x$ is five. Choose $v=x+5e_1=(8,4,0)^\top$ and set $H=I-2vv^\top/(v^\top v)$. Since $v^\top v=80$,

$$
H=\begin{pmatrix}-3/5&-4/5&0\\-4/5&3/5&0\\0&0&1\end{pmatrix}.
$$

Multiplication gives $Hx=(-5,0,0)^\top$ and $Hy=(-3/5,-4/5,2)^\top$.

The matrix is symmetric. Let $P=vv^\top/(v^\top v)$, an orthogonal rank-one projection. Then $H^2=(I-2P)^2=I$, so $H^\top H=I$. Lengths and Euclidean condition scales are preserved by this transformation.

For the matrix with columns $x,y$, the transformed first column has zeros below its leading entry. A second reflection acting only on the last two rows can eliminate the remaining bottom entry of the second column. The magnitude of the second diagonal factor is $\sqrt{(-4/5)^2+2^2}=2\sqrt{29}/5$. The product of the reflections supplies the orthogonal factor in QR.

The sign choice $v=x+\operatorname{sign}(x_1)\|x\|e_1$ avoids subtracting nearly equal first-coordinate magnitudes when the vector already points close to an axis. Either mathematical reflection can be valid, but their computed vectors can have different cancellation behavior.

#### Problem NM-83: Newton iteration for two equations

Solve $x^2+y^2=1$ and $x-y=0$ by Newton's method from $(1,0)$. Compute three updates, explain convergence after the first update, and examine the start $(0,0)$.

Solution. The vector function and Jacobian are

$$
F(x,y)=\binom{x^2+y^2-1}{x-y},\qquad
J(x,y)=\begin{pmatrix}2x&2y\\1&-1\end{pmatrix}.
$$

Solve $J h=-F$ rather than forming a Jacobian inverse explicitly. At $(1,0)$, the equations give $h_1=0,h_2=1$, producing $(1,1)$. At that point, $h_1=h_2=-1/4$, producing $(3/4,3/4)$. The next correction is $(-1/24,-1/24)$, producing $(17/24,17/24)$.

The linear equation $x-y=0$ is satisfied exactly after the first step. Every subsequent step stays on that line. Writing $x=y=z>0$ reduces Newton to the scalar equation $2z^2-1=0$:

$$
z_{\mathrm{new}}=\frac12\left(z+\frac1{2z}\right).
$$

The positive square-root argument proves convergence to $z=1/\sqrt2$ with quadratic local error reduction. The other system solution is the negative pair and is not reached from these positive iterates.

At zero the Jacobian's first row vanishes, while the first required Newton right side equals one. The linear system has no solution, so the Newton step is undefined there. The existence of a system solution does not make every starting point admissible or every Jacobian invertible.

#### Problem NM-84: trapezoidal refinement and extrapolation

Approximate $\int_0^1e^x dx$ with composite trapezoids at $h=1/2$ and $h=1/4$. Combine them to cancel the leading error and relate the result to Simpson's rule.

Solution. The two trapezoidal values are

$$
T_{1/2}=\frac14(1+2e^{1/2}+e)\approx1.753931092,
$$

$$
T_{1/4}=\frac18(1+2e^{1/4}+2e^{1/2}+2e^{3/4}+e)
\approx1.727221905.
$$

The exact integral is $e-1\approx1.718281828$. Both trapezoids overestimate because the exponential is convex.

If $T_h=I+ch^2+O(h^4)$, then $T_{h/2}=I+ch^2/4+O(h^4)$. Hence

$$
R=\frac{4T_{1/4}-T_{1/2}}3
\approx1.718318842.
$$

Substitution of the weights gives $R=(1/12)[1+4e^{1/4}+2e^{1/2}+4e^{3/4}+e]$, exactly composite Simpson with four panels.

The resulting error is about $3.70\cdot10^{-5}$. Since the fourth derivative is at most $e$ on the interval, the Simpson bound is $e(1/4)^4/180\approx5.90\cdot10^{-5}$, covering the actual error.

The cancellation used the even-power trapezoidal expansion for a smooth function. Observing two approximations close to one another does not by itself establish that expansion. Endpoint singularities, arithmetic error, or an unresolved oscillation could invalidate the extrapolated error expectation.

#### Problem NM-85: importance sampling repairs infinite variance

Estimate $I=\int_0^1x^{-1/2}dx$ by Monte Carlo. Determine the variance under uniform sampling and under density $q(x)=1/(2\sqrt x)$, and give a way to sample from the latter.

Solution. Under $U\sim\operatorname{Uniform}(0,1)$, the observation $U^{-1/2}$ has mean two. Its second moment is $\int_0^1x^{-1}dx=\infty$. The usual finite-variance Monte Carlo standard error and square-root-normal CLT therefore do not apply.

For any density positive where the integrand is nonzero, the importance-weighted observation is $f(X)/q(X)$ and has expectation $\int f(x)dx$. Its second moment is $\int f(x)^2/q(x)dx$, which must be checked separately.

The proposed density integrates to one and has CDF $Q(x)=\sqrt x$. If $V$ is uniform, set $X=V^2$. Then

$$
\frac{X^{-1/2}}{1/(2\sqrt X)}=2
$$

almost surely, producing zero variance.

This density is proportional to the positive integrand, the theoretical zero-variance choice. In a general integration problem its normalizing constant is the unknown integral, so the exact ideal density may not be available. Here the elementary integral is known and the example exposes the variance mechanism directly.

If a mixture density $q_\epsilon=(1-\epsilon)q+\epsilon$ is used with $0<\epsilon<1$, then $f/q_\epsilon\le2/(1-\epsilon)$. The estimator remains unbiased and now has a bounded observation, hence finite variance. Matching the sampling density to a singularity can be more useful than merely increasing the number of uniform draws.

#### Problem NM-86: antithetic exponential integration at equal cost

For $I=\int_0^1e^u du$, compare a plain Monte Carlo estimator using $2N$ independent function evaluations with $N$ antithetic pairs $(e^U+e^{1-U})/2$. Derive exact variances.

Solution. Put $m=e-1$. One plain observation has variance

$$
V=\frac{e^2-1}{2}-(e-1)^2.
$$

The reflected input is also uniform, so both members of an antithetic pair have this same mean and variance. Their product is identically $e$, giving covariance $C=e-(e-1)^2<0$.

One pair average has variance $(V+C)/2$. Averaging $N$ independent pairs therefore has variance $(V+C)/(2N)$. The plain estimator using $2N$ independent observations has variance $V/(2N)$. The equal-evaluation variance-reduction factor is

$$
\frac V{V+C}>1.
$$

Both estimators are unbiased. The antithetic calculation uses one random draw to generate two function inputs, while the plain calculation uses two. The comparison has fixed the number of function evaluations, which is often the main cost; fixing random draws alone would answer a different cost question.

Negative covariance follows here from the exact product and can also be understood from monotonicity: a large $U$ raises the first function value and lowers the reflected one. Antithetic pairing does not guarantee improvement for an arbitrary nonmonotone integrand. Its covariance must be analyzed for the actual transformation.

#### Problem NM-87: a manufactured heat mode and refinement order

For $u_t=u_{xx}$ on $[0,1]$ with zero boundaries and initial value $\sin(\pi x)$, derive the exact grid amplitude of the explicit scheme. Under $k=h^2/4$, determine the leading decay-rate error and the resulting spatial refinement order.

Solution. The sampled sine is an eigenvector of the second-difference matrix. Its discrete eigenvalue is

$$
\lambda_h=-\frac4{h^2}\sin^2(\pi h/2)
=-\pi^2+\frac{\pi^4h^2}{12}+O(h^4).
$$

Each time step multiplies this vector by $g=1+k\lambda_h$. After $N=T/k$ integer steps, its amplitude is $g^N$, while the exact amplitude is $e^{-\pi^2T}$.

The numerical exponential decay rate is

$$
\frac{\log g}{k}
=\lambda_h-\frac k2\lambda_h^2+O(k^2)
=-\pi^2+\pi^4\left(\frac{h^2}{12}-\frac k2\right)
+O(h^4+kh^2+k^2).
$$

With $k=h^2/4$, the leading rate error is $-\pi^4h^2/24$. The scheme damps this mode slightly too rapidly. At fixed $T$, the amplitude error is order $h^2$, so halving $h$ asymptotically reduces it by factor four.

Choosing $k=h^2/6$ cancels this leading term for the smooth heat mode. Such cancellation can produce a higher observed coupled-refinement order in this special setting. It does not make the method second order in independently varied time step, nor does it establish the same behavior for a different PDE or nonsmooth initial data.

Here $k/h^2=1/4$ satisfies the stability restriction, so the error expansion is describing a convergent scheme rather than an unstable refinement path.

#### Problem NM-88: a discrete boundary problem and its inverse bound

Discretize $-u''=2$ on $(0,1)$ with $u(0)=u(1)=0$ using centered differences. Solve for two interior points and derive a maximum-norm stability bound for the general zero-boundary discrete inverse.

Solution. With two interior points, $h=1/3$ and the system is

$$
9\begin{pmatrix}2&-1\\-1&2\end{pmatrix}
\binom{u_1}{u_2}=\binom22.
$$

Symmetry or subtraction gives $u_1=u_2$, and each equals $2/9$. The exact solution is $u(x)=x(1-x)$, which has the same sampled values. Centered second differences are exact for quadratics, so this agreement holds on every uniform mesh.

For a general grid forcing $f$, the discrete maximum principle implies that nonnegative forcing produces nonnegative solution. Thus the inverse matrix has nonnegative entries. Its row sums equal the solution for forcing identically one.

That solution is exactly $w_j=x_j(1-x_j)/2$, again because it is quadratic. Therefore

$$
\|A_h^{-1}\|_\infty=\max_j w_j\le\frac18.
$$

If a smooth exact solution has divided spatial residual bounded by $Ch^2$ and exact boundary samples, the error solves $A_he=\tau$. The inverse bound gives $\|e\|_\infty\le Ch^2/8$. This proves convergence from consistency plus a uniform stability bound for the boundary-value problem.

The inverse norm is bounded as the mesh is refined even though individual matrix entries scale like $h^{-2}$. A matrix's entry size alone does not determine its inverse's action on the relevant forcing norm.
