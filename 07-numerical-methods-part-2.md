# Numerical Methods — Part 2

A mathematical formula and a reliable numerical procedure are different objects. These notes track the approximation error, the sensitivity of the underlying problem, and the behavior of the numerical method separately.

Sources:

- Tobin Driscoll and Richard Braun, [Fundamentals of Numerical Computation, online Python edition](https://fncbook.com/python). Exercise numbers below refer to the online sections linked individually.
- Laurent Demanet, [MIT 18.330, Introduction to Numerical Analysis, Spring 2012](https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/), notes and assignments.
- Boyd and Vandenberghe, [Convex Optimization](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf), Chapter 9, for descent methods.

Small computations are shown by hand. Selected source exercises identify theoretical subparts when a larger source problem also requests programming. No software project is needed to follow the notes.

[Contents](README.md) · [Part 1](07-numerical-methods-part-1.md)

Contents:

- [7. Numerical partial differential equations](#7-numerical-partial-differential-equations)
- [8. Further worked problems](#8-further-worked-problems)

## 7. Numerical partial differential equations

### 7.1 Theory

For the heat equation $`u_t=\nu u_{xx}`$, use space grid $`x_j=jh`$ and time grid $`t_n=nk`$. Centered space differences and a forward time difference give

```math
U_j^{n+1}=U_j^n+\lambda(U_{j+1}^n-2U_j^n+U_{j-1}^n),
\qquad\lambda=\nu k/h^2.
```

Initial values populate the first time level. Boundary conditions supply values or relations at the spatial edges at every later time. An interior formula alone does not define a complete numerical problem.

Consistency means substituting a smooth exact solution into the discrete equation leaves a defect tending to zero. Stability bounds how errors or perturbations are amplified by repeated discrete steps. Convergence means the numerical solution approaches the exact one as the grid is refined. For suitable linear well-posed initial-value problems and consistent linear schemes, the Lax equivalence theorem connects stability and convergence; its functional setting and boundary treatment are part of the theorem.

For the heat scheme above, a sufficient maximum-norm stability condition is $`0\le\lambda\le1/2`$. At that condition the coefficients $`\lambda,1-2\lambda,\lambda`$ are nonnegative and sum to one. Each interior update is a weighted average of nearby values.

### 7.2 Derivations and proofs

#### Fourier-mode stability

On an infinite or periodic spatial grid, insert $`U_j^n=G^ne^{ij\theta}`$. The amplification factor is

```math
G=1+\lambda(e^{i\theta}-2+e^{-i\theta})
=1-4\lambda\sin^2(\theta/2).
```

For all modes to have $`\vert G\vert \le1`$, require $`0\le\lambda\le1/2`$. If $`\lambda\gt 1/2`$, the highest-frequency mode $`\theta=\pi`$ has $`G=1-4\lambda\lt -1`$ and grows with alternating sign.

#### Consistency and a convergence bound

For a sufficiently smooth exact solution, the time difference is $`u_t+(k/2)u_{tt}+O(k^2)`$ and the centered space difference is $`u_{xx}+(h^2/12)u_{xxxx}+O(h^4)`$. The normalized equation defect is $`O(k+h^2)`$.

With exact compatible initial and Dirichlet boundary values and $`\lambda\le1/2`$, the error update is a sup-norm contraction plus a per-step forcing of size at most $`Ck(k+h^2)`$. Summing over at most $`T/k`$ steps gives $`\Vert e^n\Vert _\infty\le CT(k+h^2)`$, proving convergence on a fixed horizon. Boundary approximation errors must be added if the boundary conditions are themselves discretized.

#### Backward time in an option PDE

For a Black-Scholes terminal-value problem, set $`\tau=T-t`$. Then $`v_\tau=\tfrac12\sigma^2S^2v_{SS}+rSv_S-rv`$. This is forward evolution in time-to-maturity with payoff as the initial condition.

On a log-price grid $`x=\log S`$, the equation becomes $`v_\tau=av_{xx}+bv_x-rv`$, with $`a=\sigma^2/2`$ and $`b=r-\sigma^2/2`$. Put $`v=e^{\alpha x+\gamma\tau}u`$, where $`\alpha=-b/(2a)`$ and $`\gamma=-b^2/(4a)-r`$. Differentiating and substituting cancels the first-derivative and zeroth-order terms, leaving $`u_\tau=au_{xx}`$. The heat restriction is then $`ak/h^2\le1/2`$.

A direct central-difference scheme with convection has additional sign and stability issues. The pure heat restriction should not be quoted as a universal sufficient rule for every Black-Scholes discretization. Boundary truncation at finite log prices also creates an error separate from the interior grid error.

#### Fourier modes and a full convergence estimate

For the explicit heat scheme, substitute a spatial mode $`u_j^n=g^ne^{ij\theta}`$. The amplification factor is

```math
g(\theta)=1+r(e^{i\theta}-2+e^{-i\theta})
=1-4r\sin^2(\theta/2).
```

Requiring $`\vert g(\theta)\vert \le1`$ for every frequency gives $`0\le r\le1/2`$. The same restriction makes the real-space stencil a convex average, proving maximum-norm stability for compatible boundary values without Fourier analysis.

For a smooth exact solution sampled on the grid, define the divided local residual $`\tau_j^n=D_tu-\nu D_{xx}u`$. Taylor expansion gives $`\vert \tau_j^n\vert \le C(k+h^2)`$. The one-step error equation has forcing $`k\tau_j^n`$. Under the maximum-norm stability bound and exact boundary data,

```math
\Vert e^{n+1}\Vert _\infty
\le\Vert e^n\Vert _\infty+kC(k+h^2).
```

Summing through $`nk\le T`$ gives $`\Vert e^n\Vert _\infty\le\Vert e^0\Vert _\infty+TC(k+h^2)`$. This is a convergence proof using the actual local residual and stable propagation. Boundary approximation errors must be included if the supplied boundary values are not exact.

On a fixed spatial mesh, diffusion is a stiff ODE system: the most negative eigenvalue has magnitude on the scale $`\nu/h^2`$. Refining space consequently forces an explicit time step on the scale $`h^2`$. This is the same stiffness issue seen in the ODE section.

#### Crank–Nicolson stability and high-frequency damping

For $`y'=\lambda y`$, Crank–Nicolson has factor $`R(z)=(1+z/2)/(1-z/2)`$, $`z=h\lambda`$. It is stable throughout the left half-plane. Along the negative real axis, however, $`R(z)\to-1`$ as $`z\to-\infty`$. Extremely fast modes can alternate with little decay under a large step.

Backward Euler has factor $`1/(1-z)\to0`$ on that axis, so it damps such modes much more strongly. Both are stable, but their treatment of unresolved high frequencies differs. Nonsmooth initial data, such as a payoff corner after reversing pricing time, can expose this difference.

#### Example NM-E23: a strongly decaying mesh mode

For a heat-system eigenvalue with $`z=-20`$, backward Euler gives factor $`1/21`$. Crank–Nicolson gives $`(1-10)/(1+10)=-9/11`$. Both magnitudes are below one, but after two steps the factors are $`1/441`$ and $`81/121`$.

The exact factor over two steps is $`e^{-40}`$, far smaller than either. Crank–Nicolson's second-order accuracy is an asymptotic statement as the step becomes small; it does not imply accurate damping of a mode with $`\vert z\vert =20`$. The example explains why stability, order, and damping should be inspected separately.

### 7.3 Worked examples

#### Example NM-E14: one explicit diffusion step

Let $`\nu=1`$, $`h=1/4`$, and $`k=1/32`$, so $`\lambda=1/2`$. With zero endpoint values and interior values $`(1,2,1)`$, the next interior values are $`(1,1,1)`$: the center becomes the average of its two neighbors, and each side becomes the average of the center and its zero boundary. All values remain within the previous range.

#### Example NM-E15: a forbidden refinement

If $`h`$ is halved while $`k`$ is held fixed, $`\lambda`$ grows by four. Starting from $`\lambda=1/4`$, this produces $`\lambda=1`$, which is unstable. To preserve a fixed explicit heat CFL number, reducing spatial spacing by two requires reducing the time step by four.

#### Example NM-E16: an exact decaying mode

For $`u_t=\nu u_{xx}`$ on $`[0,1]`$ with zero endpoints, $`u(x,t)=e^{-\nu\pi^2t}\sin(\pi x)`$ satisfies the PDE and boundary conditions. It provides an exact reference for checking the scheme's error without estimating a reference from another numerical grid.

### 7.4 Selected exercises

#### Problem NM-39: verifying a Neumann heat solution

Source: [FNC, Exercise 11.1.1](https://fncbook.com/python/blackscholes/).

Problem. Verify that $`u(x,t)=e^{-4\pi^2t}\cos(2\pi x)`$ solves the heat equation on $`[0,1]`$ with zero Neumann boundary conditions and its cosine initial value.

Solution. For $`u=e^{-4\pi^2t}\cos(2\pi x)`$, $`u_t=-4\pi^2u`$ and $`u_{xx}=-4\pi^2u`$, so $`u_t=u_{xx}`$. The spatial derivative is $`-2\pi e^{-4\pi^2t}\sin(2\pi x)`$, zero at both $`x=0`$ and $`x=1`$. At time zero the stated cosine initial condition is recovered. All three parts of the initial-boundary-value problem are satisfied.

#### Problem NM-40: a heat kernel calculation

Source: [FNC, Exercise 11.1.2](https://fncbook.com/python/blackscholes/).

Problem. Verify the heat equation for $`t^{-1/2}e^{-x^2/(4t)}`$ at positive times and explain its singular initial behavior.

Solution. Let $`u=t^{-1/2}e^{-x^2/(4t)}`$ for $`t\gt 0`$. Logarithmic differentiation gives $`u_t/u=-1/(2t)+x^2/(4t^2)`$. Also $`u_x/u=-x/(2t)`$ and $`u_{xx}/u=-1/(2t)+x^2/(4t^2)`$. Thus $`u_t=u_{xx}`$. The time-zero limit is singular; it is not a classical initial function at the origin. Multiplication by $`1/\sqrt{4\pi}`$ gives the normalized Gaussian heat kernel.

#### Problem NM-41: a mixed-boundary eigenproblem

Source: [MIT 18.330, Homework 6, Problem 1(a)](https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/28d3c9b92505167a23832f5e9c1d11d1_MIT18_330S12_hw6.pdf), eigenvalues and eigenfunctions.

Problem. Find all eigenvalues and eigenfunctions of $`-v''=\lambda v`$ on $`(0,1)`$ with $`v'(0)=0`$ and $`v(1)=0`$, excluding trivial solutions.

Solution. Solve $`-v''=\lambda v`$ with $`v'(0)=0`$, $`v(1)=0`$. For $`\lambda=\omega^2\gt 0`$, the general solution is $`A\cos\omega x+B\sin\omega x`$. The first boundary sets $`B=0`$ and the second requires $`\cos\omega=0`$. Thus $`\omega_n=(n+1/2)\pi`$, $`\lambda_n=(n+1/2)^2\pi^2`$, and $`v_n(x)=\cos[(n+1/2)\pi x]`$, $`n=0,1,\ldots`$. At $`\lambda=0`$, the affine solution is forced to zero. For negative $`\lambda`$, hyperbolic solutions and the two boundary conditions likewise give only zero. The first three nonzero modes have zero, one, and two interior sign changes respectively.

#### Problem NM-42: a Neumann compatibility condition

Source: [MIT 18.330, Homework 6, Problem 2(a,b)](https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/28d3c9b92505167a23832f5e9c1d11d1_MIT18_330S12_hw6.pdf).

Problem. For $`-[(1+x^2)u']'=f`$ on $`(0,1)`$ with zero Neumann derivatives, derive the compatibility condition, test $`f=1`$, and establish sufficiency and nonuniqueness for continuous compatible forcing.

Solution. Integrate $`-[(1+x^2)u']'=f`$ over $`[0,1]`$. The zero Neumann derivatives make the integrated left side zero, so $`\int_0^1f(x)dx=0`$ is necessary. Therefore $`f=1`$ is impossible.

For continuous $`f`$ satisfying that condition, integrate once:

```math
u'(x)=-\frac{\int_0^xf(s)ds}{1+x^2}.
```

It satisfies the left derivative condition and, by zero total integral, the right one. Integrating once more constructs $`u`$ up to an arbitrary additive constant. This proves sufficiency in the stated regularity class and explains nonuniqueness without a level normalization.

### 7.5 Additional solved problems

#### Problem NM-43: implicit diffusion

Constructed problem. Derive backward Euler with centered differences for the heat equation.

Solution. Use $`(U_j^{n+1}-U_j^n)/k=\nu(U_{j+1}^{n+1}-2U_j^{n+1}+U_{j-1}^{n+1})/h^2`$. Rearranging gives

```math
-\lambda U_{j-1}^{n+1}+(1+2\lambda)U_j^{n+1}
-\lambda U_{j+1}^{n+1}=U_j^n.
```

At each step solve a tridiagonal system, moving known new-time boundary values to the right. Fourier amplification is $`1/[1+4\lambda\sin^2(\theta/2)]`$, bounded by one for all $`\lambda\ge0`$. Unconditional stability does not remove its first-order time error.

#### Problem NM-44: Crank-Nicolson

Constructed problem. Average explicit and implicit spatial terms.

Solution. The scheme is $`(U^{n+1}-U^n)/k=(\nu/2)D_{xx}(U^{n+1}+U^n)`$. Its amplification is

```math
G=\frac{1-2\lambda\sin^2(\theta/2)}
{1+2\lambda\sin^2(\theta/2)}.
```

It has magnitude at most one for every nonnegative $`\lambda`$ and is second-order in time for smooth solutions. For large $`\lambda`$ and high frequencies, $`G`$ approaches minus one, so oscillatory components can decay slowly. Stability and strong damping are separate properties.

#### Problem NM-45: the role of a boundary error

Constructed problem. Interior differences are second-order, but a Neumann boundary derivative is approximated by a first-order one-sided difference. State what must be analyzed.

Solution. The full discrete problem includes its boundary equations. An $`O(h)`$ boundary defect can reduce the global accuracy unless the inverse discrete operator and the defect's localization provide a sharper compensating estimate. A uniform bound on the discrete solution operator in the norm being used connects defects to solution errors. Quoting only the interior stencil's $`O(h^2)`$ consistency is insufficient.

#### Problem NM-46: a maximum principle

Constructed problem. Prove range preservation for explicit heat steps with zero boundary values and $`0\le\lambda\le1/2`$.

Solution. Each interior value is a convex combination of three old values with coefficients $`\lambda,1-2\lambda,\lambda`$. It lies between their minimum and maximum. Including the fixed zero boundaries, the new global maximum cannot exceed the old maximum or zero, and the new minimum cannot be below the old minimum or zero. Induction proves range preservation at every time step. Negative weights remove this argument.

#### Problem NM-47: a convergence rate under coupled refinement

Constructed problem. If the explicit heat error is $`O(k+h^2)`$ and stability is maintained with $`k=ch^2`$, find the combined rate.

Solution. Substitute to obtain $`O(ch^2+h^2)=O(h^2)`$ at a fixed final time. Halving $`h`$ requires roughly four times as many time steps and twice as many one-dimensional spatial points, making work grow by roughly eight for a simple explicit implementation. A mathematical accuracy order and a computational cost order describe different tradeoffs.

#### Problem NM-48: an error budget

Constructed problem. A numerical option value uses a finite spatial domain, a mesh, and an iterative linear solve. Identify separate errors.

Solution. Domain truncation approximates behavior outside the chosen price interval through artificial boundaries. Discretization replaces the PDE derivatives and time evolution by finite steps. Algebraic error comes from solving the resulting finite system only approximately. Roundoff contributes further arithmetic error. Tightening the linear-solver tolerance below the dominant domain or mesh error gives little improvement in the final price; each component needs a scale appropriate to the desired accuracy.

<!-- extended-practice-start -->

## 8. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 8.1 Cancellation checked against high precision

Equivalent real-number formulas need not behave equally in floating-point arithmetic. Subtracting two nearby numbers can discard meaningful digits. A high-precision reference helps separate this arithmetic effect from approximation error in the mathematical formula.

#### Problem NM-49: small square-root difference at 10 to minus 4

Evaluate $`f(x)=\sqrt{1+x}-1`$ at $`x=10^{-4}`$ using the direct formula and the rationalized formula. Compare relative errors with a high-precision evaluation.

Solution. Multiply numerator and denominator by $`\sqrt{1+x}+1`$:

```math
\sqrt{1+x}-1=\frac x{\sqrt{1+x}+1}.
```

This is an exact identity for the positive input. The direct form subtracts quantities close to one. The rationalized denominator is close to two and contains no subtraction of nearby values.

Using binary64 arithmetic for both formulas and a 90-digit decimal calculation as reference gives:

| Quantity | Value |
|---|---:|
| High-precision reference, rounded for display | $`4.99987500624960929e-05`$ |
| Direct result | $`4.99987500623966241e-05`$ |
| Rationalized result | $`4.99987500624960997e-05`$ |
| Direct relative error | $`1.989448e-12`$ |
| Rationalized relative error | $`1.128499e-16`$ |

The direct formula first rounds $`1+x`$ and its square root near one. The subtraction then reveals the absolute error in those rounded quantities relative to a much smaller true result. At sufficiently small $`x`$, the rounded intermediate quantity can equal one exactly, producing zero.

The problem itself is well conditioned near zero when measured relatively for nonzero $`x`$. Indeed,

```math
\left\vert \frac{xf'(x)}{f(x)}\right\vert
=\frac{\sqrt{1+x}+1}{2\sqrt{1+x}}\longrightarrow1.
```

Thus the lost accuracy comes from the chosen expression, rather than large sensitivity of the mathematical function to relative input error. The stable form still rounds and is not exact, but its arithmetic avoids the damaging cancellation. For inputs near the domain boundary $`x=-1`$, the conditioning analysis would be different.

#### Problem NM-50: small square-root difference at 10 to minus 8

Evaluate $`f(x)=\sqrt{1+x}-1`$ at $`x=10^{-8}`$ using the direct formula and the rationalized formula. Compare relative errors with a high-precision evaluation.

Solution. Multiply numerator and denominator by $`\sqrt{1+x}+1`$:

```math
\sqrt{1+x}-1=\frac x{\sqrt{1+x}+1}.
```

This is an exact identity for the positive input. The direct form subtracts quantities close to one. The rationalized denominator is close to two and contains no subtraction of nearby values.

Using binary64 arithmetic for both formulas and a 90-digit decimal calculation as reference gives:

| Quantity | Value |
|---|---:|
| High-precision reference, rounded for display | $`4.99999998749999989e-09`$ |
| Direct result | $`4.99999996961264515e-09`$ |
| Rationalized result | $`4.99999998749999989e-09`$ |
| Direct relative error | $`3.577471e-09`$ |
| Rationalized relative error | $`3.373340e-17`$ |

The direct formula first rounds $`1+x`$ and its square root near one. The subtraction then reveals the absolute error in those rounded quantities relative to a much smaller true result. At sufficiently small $`x`$, the rounded intermediate quantity can equal one exactly, producing zero.

The problem itself is well conditioned near zero when measured relatively for nonzero $`x`$. Indeed,

```math
\left\vert \frac{xf'(x)}{f(x)}\right\vert
=\frac{\sqrt{1+x}+1}{2\sqrt{1+x}}\longrightarrow1.
```

Thus the lost accuracy comes from the chosen expression, rather than large sensitivity of the mathematical function to relative input error. The stable form still rounds and is not exact, but its arithmetic avoids the damaging cancellation. For inputs near the domain boundary $`x=-1`$, the conditioning analysis would be different.

#### Problem NM-51: small square-root difference at 10 to minus 12

Evaluate $`f(x)=\sqrt{1+x}-1`$ at $`x=10^{-12}`$ using the direct formula and the rationalized formula. Compare relative errors with a high-precision evaluation.

Solution. Multiply numerator and denominator by $`\sqrt{1+x}+1`$:

```math
\sqrt{1+x}-1=\frac x{\sqrt{1+x}+1}.
```

This is an exact identity for the positive input. The direct form subtracts quantities close to one. The rationalized denominator is close to two and contains no subtraction of nearby values.

Using binary64 arithmetic for both formulas and a 90-digit decimal calculation as reference gives:

| Quantity | Value |
|---|---:|
| High-precision reference, rounded for display | $`4.99999999999874984e-13`$ |
| Direct result | $`5.00044450291170506e-13`$ |
| Rationalized result | $`4.99999999999874984e-13`$ |
| Direct relative error | $`8.890058e-05`$ |
| Rationalized relative error | $`3.222232e-17`$ |

The direct formula first rounds $`1+x`$ and its square root near one. The subtraction then reveals the absolute error in those rounded quantities relative to a much smaller true result. At sufficiently small $`x`$, the rounded intermediate quantity can equal one exactly, producing zero.

The problem itself is well conditioned near zero when measured relatively for nonzero $`x`$. Indeed,

```math
\left\vert \frac{xf'(x)}{f(x)}\right\vert
=\frac{\sqrt{1+x}+1}{2\sqrt{1+x}}\longrightarrow1.
```

Thus the lost accuracy comes from the chosen expression, rather than large sensitivity of the mathematical function to relative input error. The stable form still rounds and is not exact, but its arithmetic avoids the damaging cancellation. For inputs near the domain boundary $`x=-1`$, the conditioning analysis would be different.

#### Problem NM-52: small square-root difference at 10 to minus 16

Evaluate $`f(x)=\sqrt{1+x}-1`$ at $`x=10^{-16}`$ using the direct formula and the rationalized formula. Compare relative errors with a high-precision evaluation.

Solution. Multiply numerator and denominator by $`\sqrt{1+x}+1`$:

```math
\sqrt{1+x}-1=\frac x{\sqrt{1+x}+1}.
```

This is an exact identity for the positive input. The direct form subtracts quantities close to one. The rationalized denominator is close to two and contains no subtraction of nearby values.

Using binary64 arithmetic for both formulas and a 90-digit decimal calculation as reference gives:

| Quantity | Value |
|---|---:|
| High-precision reference, rounded for display | $`4.99999999999999990e-17`$ |
| Direct result | $`0.00000000000000000e+00`$ |
| Rationalized result | $`4.99999999999999990e-17`$ |
| Direct relative error | $`1.000000e+00`$ |
| Rationalized relative error | $`4.097787e-18`$ |

The direct formula first rounds $`1+x`$ and its square root near one. The subtraction then reveals the absolute error in those rounded quantities relative to a much smaller true result. At sufficiently small $`x`$, the rounded intermediate quantity can equal one exactly, producing zero.

The problem itself is well conditioned near zero when measured relatively for nonzero $`x`$. Indeed,

```math
\left\vert \frac{xf'(x)}{f(x)}\right\vert
=\frac{\sqrt{1+x}+1}{2\sqrt{1+x}}\longrightarrow1.
```

Thus the lost accuracy comes from the chosen expression, rather than large sensitivity of the mathematical function to relative input error. The stable form still rounds and is not exact, but its arithmetic avoids the damaging cancellation. For inputs near the domain boundary $`x=-1`$, the conditioning analysis would be different.

### 8.2 Elimination, pivots, and residual checks

An LU factorization records elimination in a lower-triangular multiplier matrix and an upper-triangular reduced system. Solving by forward and back substitution gives an exact algebraic answer. A small pivot can nevertheless make the unpivoted arithmetic poor, even when the final linear system is not very ill conditioned.

#### Problem NM-53: two-by-two LU solve, case 1

Solve $`Ax=b`$ for

```math
A=\begin{pmatrix}2&1\\[0pt] 4&3\end{pmatrix}
```

 and

```math
b=\begin{pmatrix}3\\[0pt] 7\end{pmatrix}
```

 using an unpivoted LU factorization. Check the residual and discuss the first pivot.

Solution. The first pivot is $`a_{11}=2\ne0`$. Eliminate the lower-left entry using multiplier $`\ell=a_{21}/a_{11}=2`$. This gives

```math
L=\begin{pmatrix}1&0\\[0pt] 2&1\end{pmatrix},\qquad
U=\begin{pmatrix}2&1\\[0pt] 0&1\end{pmatrix}.
```

Multiplication verifies $`LU=A`$. Both diagonal entries of $`U`$ are nonzero, so the solve can proceed. Forward substitution in $`Lz=b`$ gives $`z_1=b_1`$ and $`z_2=b_2-\ell b_1`$, hence

```math
z=\begin{pmatrix}3\\[0pt] 1\end{pmatrix}
```

.

Back substitution in $`Ux=z`$ yields

```math
x_2=z_2/u_{22}=1,\qquad
x_1=(z_1-u_{12}x_2)/u_{11}=1.
```

Substituting the vector

```math
x=\begin{pmatrix}1\\[0pt] 1\end{pmatrix}
```

 into the original equations gives $`Ax=b`$ exactly. With exact fractions, the residual is zero; floating-point substitution should instead report a small norm scaled by $`\Vert A\Vert \Vert x\Vert +\Vert b\Vert`$.

For comparison, the infinity-norm condition number is $`\kappa_\infty(A)=21\approx 21.000000`$. This measures sensitivity of the supplied system. The multiplier magnitude $`\vert \ell\vert =2`$ describes this particular elimination step. A tiny first pivot can create a large multiplier even when the condition number is moderate.

Partial pivoting compares the magnitudes in the current column and swaps rows before elimination when appropriate. It avoids unnecessarily large multipliers in these two-row examples. Pivoting does not remove intrinsic ill conditioning, and a small residual alone does not guarantee a small forward error for an ill-conditioned system.

#### Problem NM-54: two-by-two LU solve, case 2

Solve $`Ax=b`$ for

```math
A=\begin{pmatrix}4&-1\\[0pt] 2&3\end{pmatrix}
```

 and

```math
b=\begin{pmatrix}1\\[0pt] 5\end{pmatrix}
```

 using an unpivoted LU factorization. Check the residual and discuss the first pivot.

Solution. The first pivot is $`a_{11}=4\ne0`$. Eliminate the lower-left entry using multiplier $`\ell=a_{21}/a_{11}=\frac{1}{2}`$. This gives

```math
L=\begin{pmatrix}1&0\\[0pt] \frac{1}{2}&1\end{pmatrix},\qquad
U=\begin{pmatrix}4&-1\\[0pt] 0&\frac{7}{2}\end{pmatrix}.
```

Multiplication verifies $`LU=A`$. Both diagonal entries of $`U`$ are nonzero, so the solve can proceed. Forward substitution in $`Lz=b`$ gives $`z_1=b_1`$ and $`z_2=b_2-\ell b_1`$, hence

```math
z=\begin{pmatrix}1\\[0pt] \frac{9}{2}\end{pmatrix}
```

.

Back substitution in $`Ux=z`$ yields

```math
x_2=z_2/u_{22}=\frac{9}{7},\qquad
x_1=(z_1-u_{12}x_2)/u_{11}=\frac{4}{7}.
```

Substituting the vector

```math
x=\begin{pmatrix}\frac{4}{7}\\[0pt] \frac{9}{7}\end{pmatrix}
```

 into the original equations gives $`Ax=b`$ exactly. With exact fractions, the residual is zero; floating-point substitution should instead report a small norm scaled by $`\Vert A\Vert \Vert x\Vert +\Vert b\Vert`$.

For comparison, the infinity-norm condition number is $`\kappa_\infty(A)=\frac{15}{7}\approx 2.142857`$. This measures sensitivity of the supplied system. The multiplier magnitude $`\vert \ell\vert =\frac{1}{2}`$ describes this particular elimination step. A tiny first pivot can create a large multiplier even when the condition number is moderate.

Partial pivoting compares the magnitudes in the current column and swaps rows before elimination when appropriate. It avoids unnecessarily large multipliers in these two-row examples. Pivoting does not remove intrinsic ill conditioning, and a small residual alone does not guarantee a small forward error for an ill-conditioned system.

#### Problem NM-55: two-by-two LU solve, case 3

Solve $`Ax=b`$ for

```math
A=\begin{pmatrix}1&2\\[0pt] 3&7\end{pmatrix}
```

 and

```math
b=\begin{pmatrix}4\\[0pt] 11\end{pmatrix}
```

 using an unpivoted LU factorization. Check the residual and discuss the first pivot.

Solution. The first pivot is $`a_{11}=1\ne0`$. Eliminate the lower-left entry using multiplier $`\ell=a_{21}/a_{11}=3`$. This gives

```math
L=\begin{pmatrix}1&0\\[0pt] 3&1\end{pmatrix},\qquad
U=\begin{pmatrix}1&2\\[0pt] 0&1\end{pmatrix}.
```

Multiplication verifies $`LU=A`$. Both diagonal entries of $`U`$ are nonzero, so the solve can proceed. Forward substitution in $`Lz=b`$ gives $`z_1=b_1`$ and $`z_2=b_2-\ell b_1`$, hence

```math
z=\begin{pmatrix}4\\[0pt] -1\end{pmatrix}
```

.

Back substitution in $`Ux=z`$ yields

```math
x_2=z_2/u_{22}=-1,\qquad
x_1=(z_1-u_{12}x_2)/u_{11}=6.
```

Substituting the vector

```math
x=\begin{pmatrix}6\\[0pt] -1\end{pmatrix}
```

 into the original equations gives $`Ax=b`$ exactly. With exact fractions, the residual is zero; floating-point substitution should instead report a small norm scaled by $`\Vert A\Vert \Vert x\Vert +\Vert b\Vert`$.

For comparison, the infinity-norm condition number is $`\kappa_\infty(A)=90\approx 90.000000`$. This measures sensitivity of the supplied system. The multiplier magnitude $`\vert \ell\vert =3`$ describes this particular elimination step. A tiny first pivot can create a large multiplier even when the condition number is moderate.

Partial pivoting compares the magnitudes in the current column and swaps rows before elimination when appropriate. It avoids unnecessarily large multipliers in these two-row examples. Pivoting does not remove intrinsic ill conditioning, and a small residual alone does not guarantee a small forward error for an ill-conditioned system.

#### Problem NM-56: two-by-two LU solve, case 4

Solve $`Ax=b`$ for

```math
A=\begin{pmatrix}\frac{1}{10000}&1\\[0pt] 1&1\end{pmatrix}
```

 and

```math
b=\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

 using an unpivoted LU factorization. Check the residual and discuss the first pivot.

Solution. The first pivot is $`a_{11}=\frac{1}{10000}\ne0`$. Eliminate the lower-left entry using multiplier $`\ell=a_{21}/a_{11}=10000`$. This gives

```math
L=\begin{pmatrix}1&0\\[0pt] 10000&1\end{pmatrix},\qquad
U=\begin{pmatrix}\frac{1}{10000}&1\\[0pt] 0&-9999\end{pmatrix}.
```

Multiplication verifies $`LU=A`$. Both diagonal entries of $`U`$ are nonzero, so the solve can proceed. Forward substitution in $`Lz=b`$ gives $`z_1=b_1`$ and $`z_2=b_2-\ell b_1`$, hence

```math
z=\begin{pmatrix}1\\[0pt] -9998\end{pmatrix}
```

.

Back substitution in $`Ux=z`$ yields

```math
x_2=z_2/u_{22}=\frac{9998}{9999},\qquad
x_1=(z_1-u_{12}x_2)/u_{11}=\frac{10000}{9999}.
```

Substituting the vector

```math
x=\begin{pmatrix}\frac{10000}{9999}\\[0pt] \frac{9998}{9999}\end{pmatrix}
```

 into the original equations gives $`Ax=b`$ exactly. With exact fractions, the residual is zero; floating-point substitution should instead report a small norm scaled by $`\Vert A\Vert \Vert x\Vert +\Vert b\Vert`$.

For comparison, the infinity-norm condition number is $`\kappa_\infty(A)=\frac{40000}{9999}\approx 4.000400`$. This measures sensitivity of the supplied system. The multiplier magnitude $`\vert \ell\vert =10000`$ describes this particular elimination step. A tiny first pivot can create a large multiplier even when the condition number is moderate.

Partial pivoting compares the magnitudes in the current column and swaps rows before elimination when appropriate. It avoids unnecessarily large multipliers in these two-row examples. Pivoting does not remove intrinsic ill conditioning, and a small residual alone does not guarantee a small forward error for an ill-conditioned system.

### 8.3 Bisection guarantees and Newton speed

Bisection reduces a known sign-changing bracket at a predictable rate. Newton's method can be much faster near a simple root, but its speed needs a separate convergence argument. A positive square root gives an exact error identity that makes the comparison transparent.

#### Problem NM-57: computing the square root of 2

Find the positive root of $`f(x)=x^2-2`$. Give a bisection iteration count ensuring midpoint error at most $`10^{-6}`$ and show Newton iterates starting from $`x_0=2`$.

Solution. The function is continuous and strictly increasing for positive $`x`$. Its values at $`1`$ and $`2`$ have opposite signs, so this interval contains exactly one positive root.

After $`n`$ bisection reductions, the bracket width is $`(2-1)/2^n`$. The midpoint error is at most half that width. Thus it is enough to choose

```math
n\ge\log_2\frac{1}{2\cdot10^{-6}},
\qquad n=19.
```

This count depends on the initial bracket and requested absolute error, not on the derivative at the root.

Newton's update simplifies to $`x_{k+1}=(x_k+2/x_k)/2`$. Starting from the upper endpoint gives:

| Iteration | Newton value | Absolute error |
|---|---:|---:|
| 0 | $`2.000000000000`$ | $`5.857864e-01`$ |
| 1 | $`1.500000000000`$ | $`8.578644e-02`$ |
| 2 | $`1.416666666667`$ | $`2.453104e-03`$ |
| 3 | $`1.414215686275`$ | $`2.123901e-06`$ |
| 4 | $`1.414213562375`$ | $`1.594724e-12`$ |

Let $`r=\sqrt{2}`$. Direct subtraction, using $`r^2=2`$, gives

```math
x_{k+1}-r
=\frac{(x_k-r)^2}{2x_k}.
```

Positive iterates remain on or above the root after one step. An iterate above the root decreases unless already equal to it. The identity then proves convergence and quadratic error reduction near the root, since the denominator approaches $`2r\gt 0`$.

A general Newton iteration does not inherit these properties automatically. A derivative close to zero or a poor initial point can cause a large step, a cycle, or departure from the domain. A safeguarded method can retain a sign-changing bracket and accept a Newton step only when it remains suitable.

#### Problem NM-58: computing the square root of 3

Find the positive root of $`f(x)=x^2-3`$. Give a bisection iteration count ensuring midpoint error at most $`10^{-6}`$ and show Newton iterates starting from $`x_0=2`$.

Solution. The function is continuous and strictly increasing for positive $`x`$. Its values at $`1`$ and $`2`$ have opposite signs, so this interval contains exactly one positive root.

After $`n`$ bisection reductions, the bracket width is $`(2-1)/2^n`$. The midpoint error is at most half that width. Thus it is enough to choose

```math
n\ge\log_2\frac{1}{2\cdot10^{-6}},
\qquad n=19.
```

This count depends on the initial bracket and requested absolute error, not on the derivative at the root.

Newton's update simplifies to $`x_{k+1}=(x_k+3/x_k)/2`$. Starting from the upper endpoint gives:

| Iteration | Newton value | Absolute error |
|---|---:|---:|
| 0 | $`2.000000000000`$ | $`2.679492e-01`$ |
| 1 | $`1.750000000000`$ | $`1.794919e-02`$ |
| 2 | $`1.732142857143`$ | $`9.204957e-05`$ |
| 3 | $`1.732050810015`$ | $`2.445850e-09`$ |
| 4 | $`1.732050807569`$ | $`0.000000e+00`$ |

Let $`r=\sqrt{3}`$. Direct subtraction, using $`r^2=3`$, gives

```math
x_{k+1}-r
=\frac{(x_k-r)^2}{2x_k}.
```

Positive iterates remain on or above the root after one step. An iterate above the root decreases unless already equal to it. The identity then proves convergence and quadratic error reduction near the root, since the denominator approaches $`2r\gt 0`$.

A general Newton iteration does not inherit these properties automatically. A derivative close to zero or a poor initial point can cause a large step, a cycle, or departure from the domain. A safeguarded method can retain a sign-changing bracket and accept a Newton step only when it remains suitable.

#### Problem NM-59: computing the square root of 5

Find the positive root of $`f(x)=x^2-5`$. Give a bisection iteration count ensuring midpoint error at most $`10^{-6}`$ and show Newton iterates starting from $`x_0=3`$.

Solution. The function is continuous and strictly increasing for positive $`x`$. Its values at $`2`$ and $`3`$ have opposite signs, so this interval contains exactly one positive root.

After $`n`$ bisection reductions, the bracket width is $`(3-2)/2^n`$. The midpoint error is at most half that width. Thus it is enough to choose

```math
n\ge\log_2\frac{1}{2\cdot10^{-6}},
\qquad n=19.
```

This count depends on the initial bracket and requested absolute error, not on the derivative at the root.

Newton's update simplifies to $`x_{k+1}=(x_k+5/x_k)/2`$. Starting from the upper endpoint gives:

| Iteration | Newton value | Absolute error |
|---|---:|---:|
| 0 | $`3.000000000000`$ | $`7.639320e-01`$ |
| 1 | $`2.333333333333`$ | $`9.726536e-02`$ |
| 2 | $`2.238095238095`$ | $`2.027261e-03`$ |
| 3 | $`2.236068895643`$ | $`9.181436e-07`$ |
| 4 | $`2.236067977500`$ | $`1.882938e-13`$ |

Let $`r=\sqrt{5}`$. Direct subtraction, using $`r^2=5`$, gives

```math
x_{k+1}-r
=\frac{(x_k-r)^2}{2x_k}.
```

Positive iterates remain on or above the root after one step. An iterate above the root decreases unless already equal to it. The identity then proves convergence and quadratic error reduction near the root, since the denominator approaches $`2r\gt 0`$.

A general Newton iteration does not inherit these properties automatically. A derivative close to zero or a poor initial point can cause a large step, a cycle, or departure from the domain. A safeguarded method can retain a sign-changing bracket and accept a Newton step only when it remains suitable.

#### Problem NM-60: computing the square root of 10

Find the positive root of $`f(x)=x^2-10`$. Give a bisection iteration count ensuring midpoint error at most $`10^{-6}`$ and show Newton iterates starting from $`x_0=4`$.

Solution. The function is continuous and strictly increasing for positive $`x`$. Its values at $`3`$ and $`4`$ have opposite signs, so this interval contains exactly one positive root.

After $`n`$ bisection reductions, the bracket width is $`(4-3)/2^n`$. The midpoint error is at most half that width. Thus it is enough to choose

```math
n\ge\log_2\frac{1}{2\cdot10^{-6}},
\qquad n=19.
```

This count depends on the initial bracket and requested absolute error, not on the derivative at the root.

Newton's update simplifies to $`x_{k+1}=(x_k+10/x_k)/2`$. Starting from the upper endpoint gives:

| Iteration | Newton value | Absolute error |
|---|---:|---:|
| 0 | $`4.000000000000`$ | $`8.377223e-01`$ |
| 1 | $`3.250000000000`$ | $`8.772234e-02`$ |
| 2 | $`3.163461538462`$ | $`1.183878e-03`$ |
| 3 | $`3.162277881693`$ | $`2.215244e-07`$ |
| 4 | $`3.162277660168`$ | $`7.549517e-15`$ |

Let $`r=\sqrt{10}`$. Direct subtraction, using $`r^2=10`$, gives

```math
x_{k+1}-r
=\frac{(x_k-r)^2}{2x_k}.
```

Positive iterates remain on or above the root after one step. An iterate above the root decreases unless already equal to it. The identity then proves convergence and quadratic error reduction near the root, since the denominator approaches $`2r\gt 0`$.

A general Newton iteration does not inherit these properties automatically. A derivative close to zero or a poor initial point can cause a large step, a cycle, or departure from the domain. A safeguarded method can retain a sign-changing bracket and accept a Newton step only when it remains suitable.

### 8.4 Gradient-descent stability along eigenvectors

On a positive-definite quadratic, a fixed gradient step acts independently on each eigenvector direction. This reduces convergence to scalar multiplier conditions. A step that improves one direction can still make the high-curvature direction oscillate or diverge.

#### Problem NM-61: fixed-step quadratic descent, case 1

Minimize $`f(x)=\tfrac12(1x_1^2+4x_2^2)`$ by gradient descent from $`x_0=(1,1)^\top`$ with step $`\eta=\frac{1}{4}`$. Derive the exact iterates and classify convergence.

Solution. The gradient is $`(1x_1,4x_2)^\top`$. The update is diagonal:

```math
x_{k+1}=
\begin{pmatrix}1-\eta1&0\\[0pt] 0&1-\eta4\end{pmatrix}x_k.
```

For the supplied step, the two multipliers are $`m_1=\frac{3}{4}`$ and $`m_2=0`$. Since both initial coordinates equal one,

```math
x_k=\begin{pmatrix}(\frac{3}{4})^k\\[0pt] (0)^k\end{pmatrix}.
```

Both multipliers have magnitude below one, so every starting vector converges to zero.

The exact objective along the iterates is

```math
f(x_k)=\frac12\left[1(\frac{3}{4})^{2k}
+4(0)^{2k}\right].
```

For convergence from every starting vector, both multiplier magnitudes must be less than one. Since the largest eigenvalue is $`4`$, the condition is $`0\lt \eta\lt 2/4`$. Equality at the upper endpoint gives a persistent oscillation in that eigenvector direction.

For this two-eigenvalue quadratic, the step minimizing the worst absolute multiplier is $`\eta_*=2/(\lambda_{\min}+\lambda_{\max})`$, with contraction factor $`(\kappa-1)/(\kappa+1)`$. A high condition number makes that best factor close to one. Preconditioning changes the relative curvature before applying descent; reducing the step alone does not remove the slow weak-curvature direction.

This exact calculation is stronger than observing a few decreasing objective values. It identifies the behavior for every iteration and explains the strict step-size restriction.

#### Problem NM-62: fixed-step quadratic descent, case 2

Minimize $`f(x)=\tfrac12(2x_1^2+5x_2^2)`$ by gradient descent from $`x_0=(1,1)^\top`$ with step $`\eta=\frac{1}{5}`$. Derive the exact iterates and classify convergence.

Solution. The gradient is $`(2x_1,5x_2)^\top`$. The update is diagonal:

```math
x_{k+1}=
\begin{pmatrix}1-\eta2&0\\[0pt] 0&1-\eta5\end{pmatrix}x_k.
```

For the supplied step, the two multipliers are $`m_1=\frac{3}{5}`$ and $`m_2=0`$. Since both initial coordinates equal one,

```math
x_k=\begin{pmatrix}(\frac{3}{5})^k\\[0pt] (0)^k\end{pmatrix}.
```

Both multipliers have magnitude below one, so every starting vector converges to zero.

The exact objective along the iterates is

```math
f(x_k)=\frac12\left[2(\frac{3}{5})^{2k}
+5(0)^{2k}\right].
```

For convergence from every starting vector, both multiplier magnitudes must be less than one. Since the largest eigenvalue is $`5`$, the condition is $`0\lt \eta\lt 2/5`$. Equality at the upper endpoint gives a persistent oscillation in that eigenvector direction.

For this two-eigenvalue quadratic, the step minimizing the worst absolute multiplier is $`\eta_*=2/(\lambda_{\min}+\lambda_{\max})`$, with contraction factor $`(\kappa-1)/(\kappa+1)`$. A high condition number makes that best factor close to one. Preconditioning changes the relative curvature before applying descent; reducing the step alone does not remove the slow weak-curvature direction.

This exact calculation is stronger than observing a few decreasing objective values. It identifies the behavior for every iteration and explains the strict step-size restriction.

#### Problem NM-63: fixed-step quadratic descent, case 3

Minimize $`f(x)=\tfrac12(1x_1^2+3x_2^2)`$ by gradient descent from $`x_0=(1,1)^\top`$ with step $`\eta=\frac{2}{3}`$. Derive the exact iterates and classify convergence.

Solution. The gradient is $`(1x_1,3x_2)^\top`$. The update is diagonal:

```math
x_{k+1}=
\begin{pmatrix}1-\eta1&0\\[0pt] 0&1-\eta3\end{pmatrix}x_k.
```

For the supplied step, the two multipliers are $`m_1=\frac{1}{3}`$ and $`m_2=-1`$. Since both initial coordinates equal one,

```math
x_k=\begin{pmatrix}(\frac{1}{3})^k\\[0pt] (-1)^k\end{pmatrix}.
```

The second multiplier equals minus one. With a nonzero second initial coordinate, that coordinate alternates without decaying, so the iterates do not converge to the minimum.

The exact objective along the iterates is

```math
f(x_k)=\frac12\left[1(\frac{1}{3})^{2k}
+3(-1)^{2k}\right].
```

For convergence from every starting vector, both multiplier magnitudes must be less than one. Since the largest eigenvalue is $`3`$, the condition is $`0\lt \eta\lt 2/3`$. Equality at the upper endpoint gives a persistent oscillation in that eigenvector direction.

For this two-eigenvalue quadratic, the step minimizing the worst absolute multiplier is $`\eta_*=2/(\lambda_{\min}+\lambda_{\max})`$, with contraction factor $`(\kappa-1)/(\kappa+1)`$. A high condition number makes that best factor close to one. Preconditioning changes the relative curvature before applying descent; reducing the step alone does not remove the slow weak-curvature direction.

This exact calculation is stronger than observing a few decreasing objective values. It identifies the behavior for every iteration and explains the strict step-size restriction.

#### Problem NM-64: fixed-step quadratic descent, case 4

Minimize $`f(x)=\tfrac12(1x_1^2+4x_2^2)`$ by gradient descent from $`x_0=(1,1)^\top`$ with step $`\eta=\frac{3}{5}`$. Derive the exact iterates and classify convergence.

Solution. The gradient is $`(1x_1,4x_2)^\top`$. The update is diagonal:

```math
x_{k+1}=
\begin{pmatrix}1-\eta1&0\\[0pt] 0&1-\eta4\end{pmatrix}x_k.
```

For the supplied step, the two multipliers are $`m_1=\frac{2}{5}`$ and $`m_2=-\frac{7}{5}`$. Since both initial coordinates equal one,

```math
x_k=\begin{pmatrix}(\frac{2}{5})^k\\[0pt] (-\frac{7}{5})^k\end{pmatrix}.
```

The second multiplier has magnitude above one. From the supplied nonzero second coordinate, its magnitude grows geometrically and the method diverges.

The exact objective along the iterates is

```math
f(x_k)=\frac12\left[1(\frac{2}{5})^{2k}
+4(-\frac{7}{5})^{2k}\right].
```

For convergence from every starting vector, both multiplier magnitudes must be less than one. Since the largest eigenvalue is $`4`$, the condition is $`0\lt \eta\lt 2/4`$. Equality at the upper endpoint gives a persistent oscillation in that eigenvector direction.

For this two-eigenvalue quadratic, the step minimizing the worst absolute multiplier is $`\eta_*=2/(\lambda_{\min}+\lambda_{\max})`$, with contraction factor $`(\kappa-1)/(\kappa+1)`$. A high condition number makes that best factor close to one. Preconditioning changes the relative curvature before applying descent; reducing the step alone does not remove the slow weak-curvature direction.

This exact calculation is stronger than observing a few decreasing objective values. It identifies the behavior for every iteration and explains the strict step-size restriction.

### 8.5 Finite differences and removable truncation error

A polynomial permits an exact expansion of finite-difference error. Forward and centered differences use different cancellations. Richardson extrapolation removes a known leading power of the step size; it does not remove floating-point error or an incorrectly assumed error expansion.

#### Problem NM-65: differentiating a cubic, case 1

For $`f(x)=x^3`$ at $`x=2`$ with $`h=\frac{1}{10}`$, compute forward and centered derivative approximations, their errors, and a Richardson extrapolation.

Solution. Expand $`(x+h)^3=x^3+3x^2h+3xh^2+h^3`$. The forward difference is therefore

```math
D_+(h)=\frac{f(x+h)-f(x)}h
=3x^2+3xh+h^2=\frac{1261}{100}.
```

The exact derivative is $`f'(x)=3x^2=12`$, so the signed forward error is $`3xh+h^2=\frac{61}{100}`$. At a special point such as zero, the first-order coefficient can vanish; an observed higher rate there does not make the method second order for general functions and points.

Subtract the expansions at $`x+h`$ and $`x-h`$. Even powers of $`h`$ cancel before division:

```math
D_0(h)=\frac{f(x+h)-f(x-h)}{2h}
=3x^2+h^2=\frac{1201}{100}.
```

Its signed error is $`h^2=\frac{1}{100}`$. At half the step, $`D_0(h/2)=3x^2+h^2/4=\frac{4801}{400}`$. Thus

```math
\frac{4D_0(h/2)-D_0(h)}3
=3x^2=12.
```

The extrapolated derivative is exact for this cubic in exact arithmetic because there are no remaining higher-order terms. For a general sufficiently smooth function, the same combination cancels an $`h^2`$ term and typically leaves an $`O(h^4)`$ error.

Every difference divides a subtraction error by $`h`$. Consequently, taking a smaller step eventually increases floating-point error even while reducing truncation error. The algebraic cancellation in Richardson extrapolation also combines rounded quantities and can magnify their noise. Its accuracy depends on being in a range where the assumed smooth truncation expansion still dominates arithmetic error.

#### Problem NM-66: differentiating a cubic, case 2

For $`f(x)=x^3`$ at $`x=-1`$ with $`h=\frac{1}{20}`$, compute forward and centered derivative approximations, their errors, and a Richardson extrapolation.

Solution. Expand $`(x+h)^3=x^3+3x^2h+3xh^2+h^3`$. The forward difference is therefore

```math
D_+(h)=\frac{f(x+h)-f(x)}h
=3x^2+3xh+h^2=\frac{1141}{400}.
```

The exact derivative is $`f'(x)=3x^2=3`$, so the signed forward error is $`3xh+h^2=-\frac{59}{400}`$. At a special point such as zero, the first-order coefficient can vanish; an observed higher rate there does not make the method second order for general functions and points.

Subtract the expansions at $`x+h`$ and $`x-h`$. Even powers of $`h`$ cancel before division:

```math
D_0(h)=\frac{f(x+h)-f(x-h)}{2h}
=3x^2+h^2=\frac{1201}{400}.
```

Its signed error is $`h^2=\frac{1}{400}`$. At half the step, $`D_0(h/2)=3x^2+h^2/4=\frac{4801}{1600}`$. Thus

```math
\frac{4D_0(h/2)-D_0(h)}3
=3x^2=3.
```

The extrapolated derivative is exact for this cubic in exact arithmetic because there are no remaining higher-order terms. For a general sufficiently smooth function, the same combination cancels an $`h^2`$ term and typically leaves an $`O(h^4)`$ error.

Every difference divides a subtraction error by $`h`$. Consequently, taking a smaller step eventually increases floating-point error even while reducing truncation error. The algebraic cancellation in Richardson extrapolation also combines rounded quantities and can magnify their noise. Its accuracy depends on being in a range where the assumed smooth truncation expansion still dominates arithmetic error.

#### Problem NM-67: differentiating a cubic, case 3

For $`f(x)=x^3`$ at $`x=0`$ with $`h=\frac{1}{100}`$, compute forward and centered derivative approximations, their errors, and a Richardson extrapolation.

Solution. Expand $`(x+h)^3=x^3+3x^2h+3xh^2+h^3`$. The forward difference is therefore

```math
D_+(h)=\frac{f(x+h)-f(x)}h
=3x^2+3xh+h^2=\frac{1}{10000}.
```

The exact derivative is $`f'(x)=3x^2=0`$, so the signed forward error is $`3xh+h^2=\frac{1}{10000}`$. At a special point such as zero, the first-order coefficient can vanish; an observed higher rate there does not make the method second order for general functions and points.

Subtract the expansions at $`x+h`$ and $`x-h`$. Even powers of $`h`$ cancel before division:

```math
D_0(h)=\frac{f(x+h)-f(x-h)}{2h}
=3x^2+h^2=\frac{1}{10000}.
```

Its signed error is $`h^2=\frac{1}{10000}`$. At half the step, $`D_0(h/2)=3x^2+h^2/4=\frac{1}{40000}`$. Thus

```math
\frac{4D_0(h/2)-D_0(h)}3
=3x^2=0.
```

The extrapolated derivative is exact for this cubic in exact arithmetic because there are no remaining higher-order terms. For a general sufficiently smooth function, the same combination cancels an $`h^2`$ term and typically leaves an $`O(h^4)`$ error.

Every difference divides a subtraction error by $`h`$. Consequently, taking a smaller step eventually increases floating-point error even while reducing truncation error. The algebraic cancellation in Richardson extrapolation also combines rounded quantities and can magnify their noise. Its accuracy depends on being in a range where the assumed smooth truncation expansion still dominates arithmetic error.

#### Problem NM-68: differentiating a cubic, case 4

For $`f(x)=x^3`$ at $`x=3`$ with $`h=\frac{1}{5}`$, compute forward and centered derivative approximations, their errors, and a Richardson extrapolation.

Solution. Expand $`(x+h)^3=x^3+3x^2h+3xh^2+h^3`$. The forward difference is therefore

```math
D_+(h)=\frac{f(x+h)-f(x)}h
=3x^2+3xh+h^2=\frac{721}{25}.
```

The exact derivative is $`f'(x)=3x^2=27`$, so the signed forward error is $`3xh+h^2=\frac{46}{25}`$. At a special point such as zero, the first-order coefficient can vanish; an observed higher rate there does not make the method second order for general functions and points.

Subtract the expansions at $`x+h`$ and $`x-h`$. Even powers of $`h`$ cancel before division:

```math
D_0(h)=\frac{f(x+h)-f(x-h)}{2h}
=3x^2+h^2=\frac{676}{25}.
```

Its signed error is $`h^2=\frac{1}{25}`$. At half the step, $`D_0(h/2)=3x^2+h^2/4=\frac{2701}{100}`$. Thus

```math
\frac{4D_0(h/2)-D_0(h)}3
=3x^2=27.
```

The extrapolated derivative is exact for this cubic in exact arithmetic because there are no remaining higher-order terms. For a general sufficiently smooth function, the same combination cancels an $`h^2`$ term and typically leaves an $`O(h^4)`$ error.

Every difference divides a subtraction error by $`h`$. Consequently, taking a smaller step eventually increases floating-point error even while reducing truncation error. The algebraic cancellation in Richardson extrapolation also combines rounded quantities and can magnify their noise. Its accuracy depends on being in a range where the assumed smooth truncation expansion still dominates arithmetic error.

### 8.6 Time-stepping amplification

For $`y'=\lambda y`$, a one-step method multiplies the current value by an amplification factor. Comparing that factor with $`e^{h\lambda}`$ separates accuracy from stability. Negative real $`\lambda`$ gives a decaying exact solution, so a growing numerical factor is immediately suspect.

#### Problem NM-69: three methods on a decaying equation, case 1

Solve $`y'=-2y`$, $`y(0)=1`$, using step $`h=\frac{1}{4}`$ for $`N=4`$ steps. Compare explicit Euler, backward Euler, and explicit trapezoidal (Heun) values and stability.

Solution. The exact terminal value is $`e^{\lambda Nh}=e^{-2}\approx 1.35335283e-01`$. Put $`z=h\lambda=-\frac{1}{2}`$.

Explicit Euler gives $`y_{n+1}=(1+z)y_n`$. Backward Euler solves $`y_{n+1}=y_n+zy_{n+1}`$, so its factor is $`1/(1-z)`$. Heun predicts $`\tilde y=y_n+zy_n`$ and averages the two slopes, giving factor $`1+z+z^2/2`$.

| Method | One-step factor | Value after $`N`$ steps |
|---|---:|---:|
| Explicit Euler | $`\frac{1}{2}`$ | $`6.25000000e-02`$ |
| Backward Euler | $`\frac{2}{3}`$ | $`1.97530864e-01`$ |
| Heun | $`\frac{5}{8}`$ | $`1.52587891e-01`$ |

On the negative real axis, explicit Euler is stable when $`-2\le z\le0`$. Heun has the same real-axis interval even though its local accuracy is higher. Backward Euler has magnitude below one for every negative real $`z`$.

Stability does not by itself guarantee an accurate value at this fixed step. A factor of zero can erase the numerical solution in one step even though the exact solution is positive. A factor of minus one gives a nondecaying sign alternation. A large stable implicit step may suppress or distort a transient that matters to the application.

Under refinement with fixed final time, explicit and backward Euler converge with first-order global error, while Heun has second-order error for this smooth problem. Stiffness concerns the practical step restriction needed to control rapidly decaying modes, which can be much smaller than the time scale on which the desired solution changes.

#### Problem NM-70: three methods on a decaying equation, case 2

Solve $`y'=-10y`$, $`y(0)=1`$, using step $`h=\frac{1}{10}`$ for $`N=5`$ steps. Compare explicit Euler, backward Euler, and explicit trapezoidal (Heun) values and stability.

Solution. The exact terminal value is $`e^{\lambda Nh}=e^{-5}\approx 6.73794700e-03`$. Put $`z=h\lambda=-1`$.

Explicit Euler gives $`y_{n+1}=(1+z)y_n`$. Backward Euler solves $`y_{n+1}=y_n+zy_{n+1}`$, so its factor is $`1/(1-z)`$. Heun predicts $`\tilde y=y_n+zy_n`$ and averages the two slopes, giving factor $`1+z+z^2/2`$.

| Method | One-step factor | Value after $`N`$ steps |
|---|---:|---:|
| Explicit Euler | $`0`$ | $`0.00000000e+00`$ |
| Backward Euler | $`\frac{1}{2}`$ | $`3.12500000e-02`$ |
| Heun | $`\frac{1}{2}`$ | $`3.12500000e-02`$ |

On the negative real axis, explicit Euler is stable when $`-2\le z\le0`$. Heun has the same real-axis interval even though its local accuracy is higher. Backward Euler has magnitude below one for every negative real $`z`$.

Stability does not by itself guarantee an accurate value at this fixed step. A factor of zero can erase the numerical solution in one step even though the exact solution is positive. A factor of minus one gives a nondecaying sign alternation. A large stable implicit step may suppress or distort a transient that matters to the application.

Under refinement with fixed final time, explicit and backward Euler converge with first-order global error, while Heun has second-order error for this smooth problem. Stiffness concerns the practical step restriction needed to control rapidly decaying modes, which can be much smaller than the time scale on which the desired solution changes.

#### Problem NM-71: three methods on a decaying equation, case 3

Solve $`y'=-10y`$, $`y(0)=1`$, using step $`h=\frac{1}{4}`$ for $`N=4`$ steps. Compare explicit Euler, backward Euler, and explicit trapezoidal (Heun) values and stability.

Solution. The exact terminal value is $`e^{\lambda Nh}=e^{-10}\approx 4.53999298e-05`$. Put $`z=h\lambda=-\frac{5}{2}`$.

Explicit Euler gives $`y_{n+1}=(1+z)y_n`$. Backward Euler solves $`y_{n+1}=y_n+zy_{n+1}`$, so its factor is $`1/(1-z)`$. Heun predicts $`\tilde y=y_n+zy_n`$ and averages the two slopes, giving factor $`1+z+z^2/2`$.

| Method | One-step factor | Value after $`N`$ steps |
|---|---:|---:|
| Explicit Euler | $`-\frac{3}{2}`$ | $`5.06250000e+00`$ |
| Backward Euler | $`\frac{2}{7}`$ | $`6.66389005e-03`$ |
| Heun | $`\frac{13}{8}`$ | $`6.97290039e+00`$ |

On the negative real axis, explicit Euler is stable when $`-2\le z\le0`$. Heun has the same real-axis interval even though its local accuracy is higher. Backward Euler has magnitude below one for every negative real $`z`$.

Stability does not by itself guarantee an accurate value at this fixed step. A factor of zero can erase the numerical solution in one step even though the exact solution is positive. A factor of minus one gives a nondecaying sign alternation. A large stable implicit step may suppress or distort a transient that matters to the application.

Under refinement with fixed final time, explicit and backward Euler converge with first-order global error, while Heun has second-order error for this smooth problem. Stiffness concerns the practical step restriction needed to control rapidly decaying modes, which can be much smaller than the time scale on which the desired solution changes.

#### Problem NM-72: three methods on a decaying equation, case 4

Solve $`y'=-100y`$, $`y(0)=1`$, using step $`h=\frac{1}{10}`$ for $`N=10`$ steps. Compare explicit Euler, backward Euler, and explicit trapezoidal (Heun) values and stability.

Solution. The exact terminal value is $`e^{\lambda Nh}=e^{-100}\approx 3.72007598e-44`$. Put $`z=h\lambda=-10`$.

Explicit Euler gives $`y_{n+1}=(1+z)y_n`$. Backward Euler solves $`y_{n+1}=y_n+zy_{n+1}`$, so its factor is $`1/(1-z)`$. Heun predicts $`\tilde y=y_n+zy_n`$ and averages the two slopes, giving factor $`1+z+z^2/2`$.

| Method | One-step factor | Value after $`N`$ steps |
|---|---:|---:|
| Explicit Euler | $`-9`$ | $`3.48678440e+09`$ |
| Backward Euler | $`\frac{1}{11}`$ | $`3.85543289e-11`$ |
| Heun | $`41`$ | $`1.34226593e+16`$ |

On the negative real axis, explicit Euler is stable when $`-2\le z\le0`$. Heun has the same real-axis interval even though its local accuracy is higher. Backward Euler has magnitude below one for every negative real $`z`$.

Stability does not by itself guarantee an accurate value at this fixed step. A factor of zero can erase the numerical solution in one step even though the exact solution is positive. A factor of minus one gives a nondecaying sign alternation. A large stable implicit step may suppress or distort a transient that matters to the application.

Under refinement with fixed final time, explicit and backward Euler converge with first-order global error, while Heun has second-order error for this smooth problem. Stiffness concerns the practical step restriction needed to control rapidly decaying modes, which can be much smaller than the time scale on which the desired solution changes.

### 8.7 A discrete heat mode and the mesh restriction

The highest-frequency mesh modes are often the first to become unstable in an explicit diffusion scheme. A three-point interior grid is small enough to calculate an eigenmode exactly and still display the stability mechanism.

#### Problem NM-73: explicit heat update with mesh ratio \frac{1}{4}

For zero-boundary heat diffusion, take three interior points and $`r=\nu k/h^2=\frac{1}{4}`$. Start from $`u^0=(1/\sqrt2,-1,1/\sqrt2)^\top`$. Find three explicit time steps and analyze the mode.

Solution. The interior second-difference matrix without the $`h^{-2}`$ factor is

```math
L=\begin{pmatrix}-2&1&0\\[0pt] 1&-2&1\\[0pt] 0&1&-2\end{pmatrix}.
```

For the supplied alternating vector $`v`$, direct multiplication gives $`Lv=(-2-\sqrt2)v`$. The explicit update $`u^{n+1}=(I+rL)u^n`$ therefore preserves this direction and changes only its magnitude:

```math
u^n=g^nv,\qquad
g=1-r(2+\sqrt2)\approx 0.146447.
```

After three steps,

```math
u^3\approx
\begin{pmatrix}0.002221\\[0pt] -0.003141\\[0pt] 0.002221\end{pmatrix}.
```

This mode decays in magnitude precisely when $`\vert g\vert \lt 1`$; a negative factor makes its sign alternate between time levels. At $`\vert g\vert \gt 1`$, this numerical high-frequency component grows even though diffusion should dissipate it.

On this particular three-interior-point grid, the most negative eigenvalue gives stability limit $`r\le2/(2+\sqrt2)\approx0.585786`$. On increasingly fine grids, the most negative scaled eigenvalue approaches minus four, so the familiar uniform-grid limit approaches $`r\le1/2`$.

The stronger maximum-principle condition $`0\le r\le1/2`$ follows directly from the stencil weights $`r,1-2r,r`$ being nonnegative and summing to one. It is a convenient grid-independent sufficient condition and ensures that one update cannot create a new interior extremum beyond the supplied values and boundaries.

Backward Euler solves $`(I-rL)u^{n+1}=u^n`$ and has factor $`1/[1+r(2+\sqrt2)]`$ on this mode, below one for every positive $`r`$. Its unconditional linear stability still leaves temporal accuracy to be controlled by the step size.

#### Problem NM-74: explicit heat update with mesh ratio \frac{1}{2}

For zero-boundary heat diffusion, take three interior points and $`r=\nu k/h^2=\frac{1}{2}`$. Start from $`u^0=(1/\sqrt2,-1,1/\sqrt2)^\top`$. Find three explicit time steps and analyze the mode.

Solution. The interior second-difference matrix without the $`h^{-2}`$ factor is

```math
L=\begin{pmatrix}-2&1&0\\[0pt] 1&-2&1\\[0pt] 0&1&-2\end{pmatrix}.
```

For the supplied alternating vector $`v`$, direct multiplication gives $`Lv=(-2-\sqrt2)v`$. The explicit update $`u^{n+1}=(I+rL)u^n`$ therefore preserves this direction and changes only its magnitude:

```math
u^n=g^nv,\qquad
g=1-r(2+\sqrt2)\approx -0.707107.
```

After three steps,

```math
u^3\approx
\begin{pmatrix}-0.250000\\[0pt] 0.353553\\[0pt] -0.250000\end{pmatrix}.
```

This mode decays in magnitude precisely when $`\vert g\vert \lt 1`$; a negative factor makes its sign alternate between time levels. At $`\vert g\vert \gt 1`$, this numerical high-frequency component grows even though diffusion should dissipate it.

On this particular three-interior-point grid, the most negative eigenvalue gives stability limit $`r\le2/(2+\sqrt2)\approx0.585786`$. On increasingly fine grids, the most negative scaled eigenvalue approaches minus four, so the familiar uniform-grid limit approaches $`r\le1/2`$.

The stronger maximum-principle condition $`0\le r\le1/2`$ follows directly from the stencil weights $`r,1-2r,r`$ being nonnegative and summing to one. It is a convenient grid-independent sufficient condition and ensures that one update cannot create a new interior extremum beyond the supplied values and boundaries.

Backward Euler solves $`(I-rL)u^{n+1}=u^n`$ and has factor $`1/[1+r(2+\sqrt2)]`$ on this mode, below one for every positive $`r`$. Its unconditional linear stability still leaves temporal accuracy to be controlled by the step size.

#### Problem NM-75: explicit heat update with mesh ratio \frac{3}{5}

For zero-boundary heat diffusion, take three interior points and $`r=\nu k/h^2=\frac{3}{5}`$. Start from $`u^0=(1/\sqrt2,-1,1/\sqrt2)^\top`$. Find three explicit time steps and analyze the mode.

Solution. The interior second-difference matrix without the $`h^{-2}`$ factor is

```math
L=\begin{pmatrix}-2&1&0\\[0pt] 1&-2&1\\[0pt] 0&1&-2\end{pmatrix}.
```

For the supplied alternating vector $`v`$, direct multiplication gives $`Lv=(-2-\sqrt2)v`$. The explicit update $`u^{n+1}=(I+rL)u^n`$ therefore preserves this direction and changes only its magnitude:

```math
u^n=g^nv,\qquad
g=1-r(2+\sqrt2)\approx -1.048528.
```

After three steps,

```math
u^3\approx
\begin{pmatrix}-0.815127\\[0pt] 1.152764\\[0pt] -0.815127\end{pmatrix}.
```

This mode decays in magnitude precisely when $`\vert g\vert \lt 1`$; a negative factor makes its sign alternate between time levels. At $`\vert g\vert \gt 1`$, this numerical high-frequency component grows even though diffusion should dissipate it.

On this particular three-interior-point grid, the most negative eigenvalue gives stability limit $`r\le2/(2+\sqrt2)\approx0.585786`$. On increasingly fine grids, the most negative scaled eigenvalue approaches minus four, so the familiar uniform-grid limit approaches $`r\le1/2`$.

The stronger maximum-principle condition $`0\le r\le1/2`$ follows directly from the stencil weights $`r,1-2r,r`$ being nonnegative and summing to one. It is a convenient grid-independent sufficient condition and ensures that one update cannot create a new interior extremum beyond the supplied values and boundaries.

Backward Euler solves $`(I-rL)u^{n+1}=u^n`$ and has factor $`1/[1+r(2+\sqrt2)]`$ on this mode, below one for every positive $`r`$. Its unconditional linear stability still leaves temporal accuracy to be controlled by the step size.

#### Problem NM-76: explicit heat update with mesh ratio 1

For zero-boundary heat diffusion, take three interior points and $`r=\nu k/h^2=1`$. Start from $`u^0=(1/\sqrt2,-1,1/\sqrt2)^\top`$. Find three explicit time steps and analyze the mode.

Solution. The interior second-difference matrix without the $`h^{-2}`$ factor is

```math
L=\begin{pmatrix}-2&1&0\\[0pt] 1&-2&1\\[0pt] 0&1&-2\end{pmatrix}.
```

For the supplied alternating vector $`v`$, direct multiplication gives $`Lv=(-2-\sqrt2)v`$. The explicit update $`u^{n+1}=(I+rL)u^n`$ therefore preserves this direction and changes only its magnitude:

```math
u^n=g^nv,\qquad
g=1-r(2+\sqrt2)\approx -2.414214.
```

After three steps,

```math
u^3\approx
\begin{pmatrix}-9.949747\\[0pt] 14.071068\\[0pt] -9.949747\end{pmatrix}.
```

This mode decays in magnitude precisely when $`\vert g\vert \lt 1`$; a negative factor makes its sign alternate between time levels. At $`\vert g\vert \gt 1`$, this numerical high-frequency component grows even though diffusion should dissipate it.

On this particular three-interior-point grid, the most negative eigenvalue gives stability limit $`r\le2/(2+\sqrt2)\approx0.585786`$. On increasingly fine grids, the most negative scaled eigenvalue approaches minus four, so the familiar uniform-grid limit approaches $`r\le1/2`$.

The stronger maximum-principle condition $`0\le r\le1/2`$ follows directly from the stencil weights $`r,1-2r,r`$ being nonnegative and summing to one. It is a convenient grid-independent sufficient condition and ensures that one update cannot create a new interior extremum beyond the supplied values and boundaries.

Backward Euler solves $`(I-rL)u^{n+1}=u^n`$ and has factor $`1/[1+r(2+\sqrt2)]`$ on this mode, below one for every positive $`r`$. Its unconditional linear stability still leaves temporal accuracy to be controlled by the step size.

### 8.8 A polynomial control variate with known moments

Monte Carlo uncertainty can be reduced by subtracting a correlated quantity with known expectation. An optimal coefficient follows by minimizing a quadratic variance expression. A fair comparison records how many random inputs and function evaluations are used.

#### Problem NM-77: integrating a uniform power 2

Estimate $`I=\int_0^1u^{2}\,du`$ using independent uniform draws. Use $`U-1/2`$ as a control variate. Derive the optimal coefficient and exact variance reduction.

Solution. For $`f(U)=U^p`$ with $`p=2`$, uniform moments give

```math
I=E[U^p]=\frac1{p+1}=\frac{1}{3},\qquad
\mathop{\mathrm{Var}} (U^p)=\frac1{2p+1}-\frac1{(p+1)^2}=\frac{4}{45}.
```

Define one controlled observation $`Y_c=U^p-c(U-1/2)`$. Its expectation remains $`I`$ for every fixed $`c`$ because the control has expectation zero.

The required covariance and control variance are

```math
\mathop{\mathrm{Cov}} (U^p,U)
=\frac1{p+2}-\frac1{2(p+1)}=\frac{1}{12},\qquad
\mathop{\mathrm{Var}} (U)=\frac1{12}.
```

Expanding the controlled variance gives $`\mathop{\mathrm{Var}} (Y_c)=\mathop{\mathrm{Var}} (U^p)-2c\mathop{\mathrm{Cov}} (U^p,U)+c^2/12`$. Its unique minimum occurs at

```math
c_*=\frac{\mathop{\mathrm{Cov}} (U^p,U)}{\mathop{\mathrm{Var}} (U)}=1.
```

The resulting one-draw variance is $`\frac{1}{180}`$, so a mean of $`N`$ independent controlled observations has variance $`\frac{1}{180}/N`$. The plain Monte Carlo variance is $`\frac{4}{45}/N`$. Their ratio is $`16\approx 16.000000`$.

Both estimators use the same number of independent uniform draws and one evaluation of $`U^p`$ per draw. The control only adds inexpensive arithmetic here. The variance calculation assumes the coefficient is fixed at its known optimal value. Estimating a coefficient on a pilot sample, or estimating it jointly on the main sample, requires accounting for that procedure. A known expectation for the control is essential; subtracting a correlated variable with an unknown mean can introduce bias.

#### Problem NM-78: integrating a uniform power 3

Estimate $`I=\int_0^1u^{3}\,du`$ using independent uniform draws. Use $`U-1/2`$ as a control variate. Derive the optimal coefficient and exact variance reduction.

Solution. For $`f(U)=U^p`$ with $`p=3`$, uniform moments give

```math
I=E[U^p]=\frac1{p+1}=\frac{1}{4},\qquad
\mathop{\mathrm{Var}} (U^p)=\frac1{2p+1}-\frac1{(p+1)^2}=\frac{9}{112}.
```

Define one controlled observation $`Y_c=U^p-c(U-1/2)`$. Its expectation remains $`I`$ for every fixed $`c`$ because the control has expectation zero.

The required covariance and control variance are

```math
\mathop{\mathrm{Cov}} (U^p,U)
=\frac1{p+2}-\frac1{2(p+1)}=\frac{3}{40},\qquad
\mathop{\mathrm{Var}} (U)=\frac1{12}.
```

Expanding the controlled variance gives $`\mathop{\mathrm{Var}} (Y_c)=\mathop{\mathrm{Var}} (U^p)-2c\mathop{\mathrm{Cov}} (U^p,U)+c^2/12`$. Its unique minimum occurs at

```math
c_*=\frac{\mathop{\mathrm{Cov}} (U^p,U)}{\mathop{\mathrm{Var}} (U)}=\frac{9}{10}.
```

The resulting one-draw variance is $`\frac{9}{700}`$, so a mean of $`N`$ independent controlled observations has variance $`\frac{9}{700}/N`$. The plain Monte Carlo variance is $`\frac{9}{112}/N`$. Their ratio is $`\frac{25}{4}\approx 6.250000`$.

Both estimators use the same number of independent uniform draws and one evaluation of $`U^p`$ per draw. The control only adds inexpensive arithmetic here. The variance calculation assumes the coefficient is fixed at its known optimal value. Estimating a coefficient on a pilot sample, or estimating it jointly on the main sample, requires accounting for that procedure. A known expectation for the control is essential; subtracting a correlated variable with an unknown mean can introduce bias.

#### Problem NM-79: integrating a uniform power 4

Estimate $`I=\int_0^1u^{4}\,du`$ using independent uniform draws. Use $`U-1/2`$ as a control variate. Derive the optimal coefficient and exact variance reduction.

Solution. For $`f(U)=U^p`$ with $`p=4`$, uniform moments give

```math
I=E[U^p]=\frac1{p+1}=\frac{1}{5},\qquad
\mathop{\mathrm{Var}} (U^p)=\frac1{2p+1}-\frac1{(p+1)^2}=\frac{16}{225}.
```

Define one controlled observation $`Y_c=U^p-c(U-1/2)`$. Its expectation remains $`I`$ for every fixed $`c`$ because the control has expectation zero.

The required covariance and control variance are

```math
\mathop{\mathrm{Cov}} (U^p,U)
=\frac1{p+2}-\frac1{2(p+1)}=\frac{1}{15},\qquad
\mathop{\mathrm{Var}} (U)=\frac1{12}.
```

Expanding the controlled variance gives $`\mathop{\mathrm{Var}} (Y_c)=\mathop{\mathrm{Var}} (U^p)-2c\mathop{\mathrm{Cov}} (U^p,U)+c^2/12`$. Its unique minimum occurs at

```math
c_*=\frac{\mathop{\mathrm{Cov}} (U^p,U)}{\mathop{\mathrm{Var}} (U)}=\frac{4}{5}.
```

The resulting one-draw variance is $`\frac{4}{225}`$, so a mean of $`N`$ independent controlled observations has variance $`\frac{4}{225}/N`$. The plain Monte Carlo variance is $`\frac{16}{225}/N`$. Their ratio is $`4\approx 4.000000`$.

Both estimators use the same number of independent uniform draws and one evaluation of $`U^p`$ per draw. The control only adds inexpensive arithmetic here. The variance calculation assumes the coefficient is fixed at its known optimal value. Estimating a coefficient on a pilot sample, or estimating it jointly on the main sample, requires accounting for that procedure. A known expectation for the control is essential; subtracting a correlated variable with an unknown mean can introduce bias.

#### Problem NM-80: integrating a uniform power 5

Estimate $`I=\int_0^1u^{5}\,du`$ using independent uniform draws. Use $`U-1/2`$ as a control variate. Derive the optimal coefficient and exact variance reduction.

Solution. For $`f(U)=U^p`$ with $`p=5`$, uniform moments give

```math
I=E[U^p]=\frac1{p+1}=\frac{1}{6},\qquad
\mathop{\mathrm{Var}} (U^p)=\frac1{2p+1}-\frac1{(p+1)^2}=\frac{25}{396}.
```

Define one controlled observation $`Y_c=U^p-c(U-1/2)`$. Its expectation remains $`I`$ for every fixed $`c`$ because the control has expectation zero.

The required covariance and control variance are

```math
\mathop{\mathrm{Cov}} (U^p,U)
=\frac1{p+2}-\frac1{2(p+1)}=\frac{5}{84},\qquad
\mathop{\mathrm{Var}} (U)=\frac1{12}.
```

Expanding the controlled variance gives $`\mathop{\mathrm{Var}} (Y_c)=\mathop{\mathrm{Var}} (U^p)-2c\mathop{\mathrm{Cov}} (U^p,U)+c^2/12`$. Its unique minimum occurs at

```math
c_*=\frac{\mathop{\mathrm{Cov}} (U^p,U)}{\mathop{\mathrm{Var}} (U)}=\frac{5}{7}.
```

The resulting one-draw variance is $`\frac{100}{4851}`$, so a mean of $`N`$ independent controlled observations has variance $`\frac{100}{4851}/N`$. The plain Monte Carlo variance is $`\frac{25}{396}/N`$. Their ratio is $`\frac{49}{16}\approx 3.062500`$.

Both estimators use the same number of independent uniform draws and one evaluation of $`U^p`$ per draw. The control only adds inexpensive arithmetic here. The variance calculation assumes the coefficient is fixed at its known optimal value. Estimating a coefficient on a pilot sample, or estimating it jointly on the main sample, requires accounting for that procedure. A known expectation for the control is essential; subtracting a correlated variable with an unknown mean can introduce bias.

<!-- further-distinct-problems-NM -->

### 8.9 Further accuracy and convergence problems

These additional problems are constructed for these notes.

#### Problem NM-81: residual correction as an iteration

Solve $`Ax=b`$ for

```math
A=\begin{pmatrix}4&1\\[0pt] 1&3\end{pmatrix}
```

 and $`b=(1,2)^\top`$ using $`x_{k+1}=x_k+B(b-Ax_k)`$ with $`B=\mathop{\mathrm{diag}} (1/4,1/3)`$ and $`x_0=0`$. Compute two steps and prove convergence.

Solution. The exact inverse gives $`x^*=(1/11,7/11)^\top`$. The first residual is $`b`$, so $`x_1=(1/4,2/3)^\top`$. The recurrence can also be written

```math
x_{1,k+1}=\frac{1-x_{2,k}}4,\qquad
x_{2,k+1}=\frac{2-x_{1,k}}3.
```

Thus $`x_2=(1/12,7/12)^\top`$. Each update uses both old coordinates, rather than substituting an already updated coordinate into the second equation.

Let $`e_k=x_k-x^*`$. Since $`Ax^*=b`$,

```math
e_{k+1}=(I-BA)e_k,\qquad
I-BA=\begin{pmatrix}0&-1/4\\[0pt] -1/3&0\end{pmatrix}.
```

The infinity norm of this matrix is $`1/3\lt 1`$, so $`\Vert e_k\Vert _\infty\le(1/3)^k\Vert e_0\Vert _\infty`$. Its spectral radius is $`1/\sqrt{12}`$, consistent with convergence.

The same algebra describes a residual-correction step using an approximate inverse or approximate linear solve. If the correction used the exact inverse and the residual were exact, one step would eliminate the error. In floating-point iterative refinement, residual accuracy and correction-solve accuracy both matter. Recomputing a residual more accurately can recover information lost in the initial solve, but convergence requires control of the resulting error-propagation map.

#### Problem NM-82: a Householder reflection

Construct an orthogonal reflection taking $`x=(3,4,0)^\top`$ to a multiple of the first coordinate vector. Apply it to $`y=(1,0,2)^\top`$ and explain its first QR elimination step.

Solution. The norm of $`x`$ is five. Choose $`v=x+5e_1=(8,4,0)^\top`$ and set $`H=I-2vv^\top/(v^\top v)`$. Since $`v^\top v=80`$,

```math
H=\begin{pmatrix}-3/5&-4/5&0\\[0pt] -4/5&3/5&0\\[0pt] 0&0&1\end{pmatrix}.
```

Multiplication gives $`Hx=(-5,0,0)^\top`$ and $`Hy=(-3/5,-4/5,2)^\top`$.

The matrix is symmetric. Let $`P=vv^\top/(v^\top v)`$, an orthogonal rank-one projection. Then $`H^2=(I-2P)^2=I`$, so $`H^\top H=I`$. Lengths and Euclidean condition scales are preserved by this transformation.

For the matrix with columns $`x,y`$, the transformed first column has zeros below its leading entry. A second reflection acting only on the last two rows can eliminate the remaining bottom entry of the second column. The magnitude of the second diagonal factor is $`\sqrt{(-4/5)^2+2^2}=2\sqrt{29}/5`$. The product of the reflections supplies the orthogonal factor in QR.

The sign choice $`v=x+\mathop{\mathrm{sign}} (x_1)\Vert x\Vert e_1`$ avoids subtracting nearly equal first-coordinate magnitudes when the vector already points close to an axis. Either mathematical reflection can be valid, but their computed vectors can have different cancellation behavior.

#### Problem NM-83: Newton iteration for two equations

Solve $`x^2+y^2=1`$ and $`x-y=0`$ by Newton's method from $`(1,0)`$. Compute three updates, explain convergence after the first update, and examine the start $`(0,0)`$.

Solution. The vector function and Jacobian are

```math
F(x,y)=\binom{x^2+y^2-1}{x-y},\qquad
J(x,y)=\begin{pmatrix}2x&2y\\[0pt] 1&-1\end{pmatrix}.
```

Solve $`J h=-F`$ rather than forming a Jacobian inverse explicitly. At $`(1,0)`$, the equations give $`h_1=0,h_2=1`$, producing $`(1,1)`$. At that point, $`h_1=h_2=-1/4`$, producing $`(3/4,3/4)`$. The next correction is $`(-1/24,-1/24)`$, producing $`(17/24,17/24)`$.

The linear equation $`x-y=0`$ is satisfied exactly after the first step. Every subsequent step stays on that line. Writing $`x=y=z\gt 0`$ reduces Newton to the scalar equation $`2z^2-1=0`$:

```math
z_{\mathrm{new}}=\frac12\left(z+\frac1{2z}\right).
```

The positive square-root argument proves convergence to $`z=1/\sqrt2`$ with quadratic local error reduction. The other system solution is the negative pair and is not reached from these positive iterates.

At zero the Jacobian's first row vanishes, while the first required Newton right side equals one. The linear system has no solution, so the Newton step is undefined there. The existence of a system solution does not make every starting point admissible or every Jacobian invertible.

#### Problem NM-84: trapezoidal refinement and extrapolation

Approximate $`\int_0^1e^x dx`$ with composite trapezoids at $`h=1/2`$ and $`h=1/4`$. Combine them to cancel the leading error and relate the result to Simpson's rule.

Solution. The two trapezoidal values are

```math
T_{1/2}=\frac14(1+2e^{1/2}+e)\approx1.753931092,
```

```math
T_{1/4}=\frac18(1+2e^{1/4}+2e^{1/2}+2e^{3/4}+e)
\approx1.727221905.
```

The exact integral is $`e-1\approx1.718281828`$. Both trapezoids overestimate because the exponential is convex.

If $`T_h=I+ch^2+O(h^4)`$, then $`T_{h/2}=I+ch^2/4+O(h^4)`$. Hence

```math
R=\frac{4T_{1/4}-T_{1/2}}3
\approx1.718318842.
```

Substitution of the weights gives $`R=(1/12)[1+4e^{1/4}+2e^{1/2}+4e^{3/4}+e]`$, exactly composite Simpson with four panels.

The resulting error is about $`3.70\cdot10^{-5}`$. Since the fourth derivative is at most $`e`$ on the interval, the Simpson bound is $`e(1/4)^4/180\approx5.90\cdot10^{-5}`$, covering the actual error.

The cancellation used the even-power trapezoidal expansion for a smooth function. Observing two approximations close to one another does not by itself establish that expansion. Endpoint singularities, arithmetic error, or an unresolved oscillation could invalidate the extrapolated error expectation.

#### Problem NM-85: importance sampling repairs infinite variance

Estimate $`I=\int_0^1x^{-1/2}dx`$ by Monte Carlo. Determine the variance under uniform sampling and under density $`q(x)=1/(2\sqrt x)`$, and give a way to sample from the latter.

Solution. Under $`U\sim\mathop{\mathrm{Uniform}} (0,1)`$, the observation $`U^{-1/2}`$ has mean two. Its second moment is $`\int_0^1x^{-1}dx=\infty`$. The usual finite-variance Monte Carlo standard error and square-root-normal CLT therefore do not apply.

For any density positive where the integrand is nonzero, the importance-weighted observation is $`f(X)/q(X)`$ and has expectation $`\int f(x)dx`$. Its second moment is $`\int f(x)^2/q(x)dx`$, which must be checked separately.

The proposed density integrates to one and has CDF $`Q(x)=\sqrt x`$. If $`V`$ is uniform, set $`X=V^2`$. Then

```math
\frac{X^{-1/2}}{1/(2\sqrt X)}=2
```

almost surely, producing zero variance.

This density is proportional to the positive integrand, the theoretical zero-variance choice. In a general integration problem its normalizing constant is the unknown integral, so the exact ideal density may not be available. Here the elementary integral is known and the example exposes the variance mechanism directly.

If a mixture density $`q_\epsilon=(1-\epsilon)q+\epsilon`$ is used with $`0\lt \epsilon\lt 1`$, then $`f/q_\epsilon\le2/(1-\epsilon)`$. The estimator remains unbiased and now has a bounded observation, hence finite variance. Matching the sampling density to a singularity can be more useful than merely increasing the number of uniform draws.

#### Problem NM-86: antithetic exponential integration at equal cost

For $`I=\int_0^1e^u du`$, compare a plain Monte Carlo estimator using $`2N`$ independent function evaluations with $`N`$ antithetic pairs $`(e^U+e^{1-U})/2`$. Derive exact variances.

Solution. Put $`m=e-1`$. One plain observation has variance

```math
V=\frac{e^2-1}{2}-(e-1)^2.
```

The reflected input is also uniform, so both members of an antithetic pair have this same mean and variance. Their product is identically $`e`$, giving covariance $`C=e-(e-1)^2\lt 0`$.

One pair average has variance $`(V+C)/2`$. Averaging $`N`$ independent pairs therefore has variance $`(V+C)/(2N)`$. The plain estimator using $`2N`$ independent observations has variance $`V/(2N)`$. The equal-evaluation variance-reduction factor is

```math
\frac V{V+C}\gt 1.
```

Both estimators are unbiased. The antithetic calculation uses one random draw to generate two function inputs, while the plain calculation uses two. The comparison has fixed the number of function evaluations, which is often the main cost; fixing random draws alone would answer a different cost question.

Negative covariance follows here from the exact product and can also be understood from monotonicity: a large $`U`$ raises the first function value and lowers the reflected one. Antithetic pairing does not guarantee improvement for an arbitrary nonmonotone integrand. Its covariance must be analyzed for the actual transformation.

#### Problem NM-87: a manufactured heat mode and refinement order

For $`u_t=u_{xx}`$ on $`[0,1]`$ with zero boundaries and initial value $`\sin(\pi x)`$, derive the exact grid amplitude of the explicit scheme. Under $`k=h^2/4`$, determine the leading decay-rate error and the resulting spatial refinement order.

Solution. The sampled sine is an eigenvector of the second-difference matrix. Its discrete eigenvalue is

```math
\lambda_h=-\frac4{h^2}\sin^2(\pi h/2)
=-\pi^2+\frac{\pi^4h^2}{12}+O(h^4).
```

Each time step multiplies this vector by $`g=1+k\lambda_h`$. After $`N=T/k`$ integer steps, its amplitude is $`g^N`$, while the exact amplitude is $`e^{-\pi^2T}`$.

The numerical exponential decay rate is

```math
\frac{\log g}{k}
=\lambda_h-\frac k2\lambda_h^2+O(k^2)
=-\pi^2+\pi^4\left(\frac{h^2}{12}-\frac k2\right)
+O(h^4+kh^2+k^2).
```

With $`k=h^2/4`$, the leading rate error is $`-\pi^4h^2/24`$. The scheme damps this mode slightly too rapidly. At fixed $`T`$, the amplitude error is order $`h^2`$, so halving $`h`$ asymptotically reduces it by factor four.

Choosing $`k=h^2/6`$ cancels this leading term for the smooth heat mode. Such cancellation can produce a higher observed coupled-refinement order in this special setting. It does not make the method second order in independently varied time step, nor does it establish the same behavior for a different PDE or nonsmooth initial data.

Here $`k/h^2=1/4`$ satisfies the stability restriction, so the error expansion is describing a convergent scheme rather than an unstable refinement path.

#### Problem NM-88: a discrete boundary problem and its inverse bound

Discretize $`-u''=2`$ on $`(0,1)`$ with $`u(0)=u(1)=0`$ using centered differences. Solve for two interior points and derive a maximum-norm stability bound for the general zero-boundary discrete inverse.

Solution. With two interior points, $`h=1/3`$ and the system is

```math
9\begin{pmatrix}2&-1\\[0pt] -1&2\end{pmatrix}
\binom{u_1}{u_2}=\binom22.
```

Symmetry or subtraction gives $`u_1=u_2`$, and each equals $`2/9`$. The exact solution is $`u(x)=x(1-x)`$, which has the same sampled values. Centered second differences are exact for quadratics, so this agreement holds on every uniform mesh.

For a general grid forcing $`f`$, the discrete maximum principle implies that nonnegative forcing produces nonnegative solution. Thus the inverse matrix has nonnegative entries. Its row sums equal the solution for forcing identically one.

That solution is exactly $`w_j=x_j(1-x_j)/2`$, again because it is quadratic. Therefore

```math
\Vert A_h^{-1}\Vert _\infty=\max_j w_j\le\frac18.
```

If a smooth exact solution has divided spatial residual bounded by $`Ch^2`$ and exact boundary samples, the error solves $`A_he=\tau`$. The inverse bound gives $`\Vert e\Vert _\infty\le Ch^2/8`$. This proves convergence from consistency plus a uniform stability bound for the boundary-value problem.

The inverse norm is bounded as the mesh is refined even though individual matrix entries scale like $`h^{-2}`$. A matrix's entry size alone does not determine its inverse's action on the relevant forcing norm.

[Contents](README.md) · [Part 1](07-numerical-methods-part-1.md)
