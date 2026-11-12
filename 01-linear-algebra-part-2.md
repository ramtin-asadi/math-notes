# Linear Algebra — Part 2

These notes use real finite-dimensional spaces unless a different field is stated. Vectors are columns, $`A^\top`$ is the transpose, $`I`$ is the identity, and $`\Vert x\Vert _2=(x^\top x)^{1/2}`$. For complex vectors, replace transposes by conjugate transposes and use the Hermitian inner product.

Sources:

- Gilbert Strang, [MIT 18.06SC, Fall 2011](https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/), especially Units II and III.
- Sheldon Axler, [Linear Algebra Done Right, fourth edition](https://linear.axler.net/LADR4e.pdf), Chapters 2, 3, 6, and 7. Exercise numbers below refer to this edition.

The explanations and solutions are written for these notes. Source problems are briefly paraphrased and linked; examples and additional problems are constructed here.

[Contents](README.md) · [Part 1](01-linear-algebra-part-1.md)

Contents:

- [6. Matrix differentiation](#6-matrix-differentiation)
- [7. Further worked problems](#7-further-worked-problems)

## 6. Matrix differentiation

### Further problems on sensitivity and regularization

#### Problem LA-53: submultiplicativity

Source: [Axler, Exercises 7F, 5](https://linear.axler.net/LADR4e.pdf).

Problem. Prove submultiplicativity $`\Vert ST\Vert \le\Vert S\Vert \Vert T\Vert`$ for compatible induced operator norms.

Solution. For every unit vector $`x`$, $`\Vert STx\Vert \le\Vert S\Vert \Vert Tx\Vert \le\Vert S\Vert \Vert T\Vert`$. Taking the supremum over $`x`$ gives $`\Vert ST\Vert \le\Vert S\Vert \Vert T\Vert`$. If $`Tx=0`$, the first inequality still holds. The proof uses induced norms with compatible input, intermediate, and output vector norms; arbitrary matrix-entry norms need not obey the same formula without an additional constant.

#### Problem LA-54: the norm of a Gram matrix

Source: [Axler, Exercises 7F, 19](https://linear.axler.net/LADR4e.pdf).

Problem. Prove $`\Vert A^\top A\Vert _2=\Vert A\Vert _2^2`$ and explain the role of the Gram matrix's positivity.

Solution. The eigenvalues of $`A^\top A`$ are $`\sigma_i^2`$, all nonnegative, so its operator norm is $`\sigma_{\max}^2=\Vert A\Vert _2^2`$. Alternatively, the Rayleigh quotient gives $`\Vert A^\top A\Vert _2=\max_{\Vert x\Vert =1}x^\top A^\top Ax=\max_{\Vert x\Vert =1}\Vert Ax\Vert ^2`$. The second equality relies on symmetry and PSD of the Gram matrix.

#### Problem LA-55: the smallest possible condition number

Source: [Axler, Exercises 7F, 24](https://linear.axler.net/LADR4e.pdf).

Problem. Characterize the invertible matrices satisfying $`\Vert A^{-1}\Vert _2=\Vert A\Vert _2^{-1}`$, and relate the result to condition number one.

Solution. For invertible $`A`$, $`\Vert A^{-1}\Vert _2=\Vert A\Vert _2^{-1}`$ exactly when $`\sigma_{\max}/\sigma_{\min}=1`$. Since every singular value lies between those two, they must all equal a common positive number $`s=\Vert A\Vert _2`$. The SVD then gives $`A=sUV^\top`$, so $`A/\Vert A\Vert _2`$ is orthogonal. Conversely a positive scalar times an orthogonal matrix has all singular values equal and condition number one.

#### Problem LA-56: product order and operator norm

Source: [Axler, Exercises 7F, 6](https://linear.axler.net/LADR4e.pdf).

Problem. Determine whether $`\Vert ST\Vert _2=\Vert TS\Vert _2`$ must hold for two square matrices. Supply a proof or a counterexample.

Solution. The assertion $`\Vert ST\Vert =\Vert TS\Vert`$ is false. Set

```math
S=\begin{pmatrix}1&0\\[0pt] 0&0\end{pmatrix},\qquad
T=\begin{pmatrix}0&1\\[0pt] 0&0\end{pmatrix}.
```

Then $`ST=T`$, with norm one, but $`TS=0`$, with norm zero. Matrix products act in order: $`TS`$ first removes the second input coordinate, leaving nothing for $`T`$ to transfer into the first coordinate.

#### Problem LA-57: a nearly singular pair of equations

Constructed problem. Solve $`x+y=2`$, $`x+(1+\varepsilon)y=2+\varepsilon`$ for $`\varepsilon\gt 0`$, and then perturb the second right-hand side by $`\delta`$.

Solution. Subtraction gives $`\varepsilon y=\varepsilon`$, so the exact solution is $`(1,1)`$. After perturbation, $`\varepsilon\widetilde y=\varepsilon+\delta`$, so $`\widetilde y=1+\delta/\varepsilon`$ and $`\widetilde x=1-\delta/\varepsilon`$. The error norm is $`\sqrt2\vert \delta\vert /\varepsilon`$. The equations become nearly duplicates as $`\varepsilon\downarrow0`$, and the subtraction needed to identify $`y`$ divides a small difference by another small number.

#### Problem LA-58: a residual-based error bound

Constructed problem. If $`A`$ is PD with $`\lambda_{\min}(A)=m\gt 0`$, bound the error of an approximate solution using its residual.

Solution. With $`r=b-A\widetilde x`$ and exact solution $`x`$, $`A(x-\widetilde x)=r`$. The inverse norm is $`1/m`$, so $`\Vert x-\widetilde x\Vert _2\le\Vert r\Vert _2/m`$. Equality is possible when $`r`$ lies along an eigenvector for $`m`$. To guarantee solution error at most $`\eta`$, the residual should therefore be at most $`m\eta`$. A stopping rule based only on a fixed residual tolerance can fail to control error uniformly across differently conditioned systems.

#### Problem LA-59: ridge bias in each direction

Constructed problem. Suppose $`b=Ax_0`$ exactly and compute $`x_\lambda-x_0`$.

Solution. Expand $`x_0=\sum_i c_iv_i`$ in a full right singular basis. For each positive singular direction,

```math
(x_\lambda)_i=\frac{\sigma_i^2}{\sigma_i^2+\lambda}c_i,
\qquad
(x_\lambda-x_0)_i=-\frac{\lambda}{\sigma_i^2+\lambda}c_i.
```

For a zero singular direction, $`x_\lambda`$ has coefficient zero and the error is $`-c_i`$. Therefore exact data do not make ridge unbiased relative to an arbitrary true coefficient vector. As $`\lambda\downarrow0`$, the solution approaches $`A^+b`$, which recovers only the component of $`x_0`$ perpendicular to the kernel.

#### Problem LA-60: shrinkage toward a general positive target

Constructed problem. Let $`S`$ be PSD, $`T`$ PD, and $`S_\alpha=(1-\alpha)S+\alpha T`$ for $`0\lt \alpha\le1`$. Prove invertibility and give a lower eigenvalue bound.

Solution. For a unit vector $`x`$,

```math
x^\top S_\alpha x=(1-\alpha)x^\top Sx+\alpha x^\top Tx
\ge\alpha\lambda_{\min}(T)\gt 0.
```

Thus $`S_\alpha`$ is PD and $`\lambda_{\min}(S_\alpha)\ge\alpha\lambda_{\min}(T)`$. Its largest eigenvalue is at most $`(1-\alpha)\lambda_{\max}(S)+\alpha\lambda_{\max}(T)`$, giving a condition-number upper bound by dividing by the earlier lower bound. Unlike an identity target, a general $`T`$ need not share eigenvectors with $`S`$.

### 6.1 Theory

For a differentiable scalar function $`f:\mathbb R^n\to\mathbb R`$, the gradient is defined by

```math
df=\nabla f(x)^\top dx.
```

The Hessian is the matrix of second partial derivatives. For twice continuously differentiable $`f`$, it is symmetric. A positive-definite Hessian everywhere on a convex domain makes $`f`$ strictly convex; PSD gives convexity.

Useful identities are

```math
\nabla_x(a^\top x)=a,\qquad
\nabla_x(x^\top Ax)=(A+A^\top)x,\qquad
\nabla_x\Vert Ax-b\Vert ^2=2A^\top(Ax-b).
```

For a matrix variable, the Frobenius inner product is $`\langle G,H\rangle_F=\mathop{\mathrm{tr}} (G^\top H)`$. Define the matrix gradient by $`df=\mathop{\mathrm{tr}} ((\nabla_Xf)^\top dX)`$. This convention avoids ambiguous numerator-layout and denominator-layout derivative tables.

### 6.2 Derivations and proofs

#### Dimensions before derivatives

For $`A\in\mathbb R^{m\times n}`$ and $`x\in\mathbb R^n`$, the residual $`Ax-b`$ has $`m`$ coordinates. Its gradient with respect to $`x`$ must have $`n`$ coordinates. The expression $`A^\top(Ax-b)`$ has the right size; $`A(Ax-b)`$ may not even be defined. Checking dimensions catches many derivative mistakes before any algebra.

For $`X\in\mathbb R^{n\times k}`$, a gradient with respect to $`X`$ has the same shape as $`X`$. The differential convention $`\mathop{\mathrm{tr}} (G^\top dX)`$ makes this visible. A trace is a scalar, and cyclic rearrangement is allowed only when the resulting products are defined. One cannot arbitrarily swap factors: $`\mathop{\mathrm{tr}} (ABC)=\mathop{\mathrm{tr}} (BCA)`$, but usually $`\mathop{\mathrm{tr}} (ABC)\ne\mathop{\mathrm{tr}} (ACB)`$.

#### A quadratic minimum along feasible directions

Positive definiteness of $`H`$ on all of $`\mathbb R^n`$ is sufficient for the constrained quadratic, but can be stronger than needed. Let the columns of $`Z`$ form a basis of $`\ker C^\top`$, and let $`x_p`$ be one feasible point. Every feasible $`x`$ is $`x_p+Zz`$. Substitution turns the problem into

```math
\frac12 z^\top(Z^\top HZ)z+
z^\top Z^\top(Hx_p-g)+\text{constant}.
```

If $`Z^\top HZ`$ is PD, this reduced problem has a unique minimizer even if $`H`$ is singular or indefinite in infeasible directions. Its stationarity equation is

```math
Z^\top HZ\,z=Z^\top(g-Hx_p).
```

If the reduced Hessian is merely PSD, use the range condition from Problem LA-17. A linear term outside its range makes the objective unbounded below along an allowed flat direction. This approach handles singular covariance examples without an unjustified matrix inverse.

#### Differentiating an inverse and a log determinant

For an invertible matrix $`X`$, differentiate $`XX^{-1}=I`$:

```math
(dX)X^{-1}+X\,d(X^{-1})=0,
\qquad d(X^{-1})=-X^{-1}(dX)X^{-1}.
```

For PD $`X`$, the determinant expansion $`\det(I+\varepsilon H)=1+\varepsilon\mathop{\mathrm{tr}} H+O(\varepsilon^2)`$ gives

```math
d\log\det X=\mathop{\mathrm{tr}} (X^{-1}dX).
```

To justify that expansion, use the permutation definition of a determinant. The linear terms come from choosing a single diagonal $`\varepsilon H_{ii}`$ and all other diagonal ones; any nonidentity permutation requires at least two off-diagonal factors and contributes at order at least two. The gradient of $`\log\det X`$ is $`X^{-\top}`$, which equals $`X^{-1}`$ on symmetric PD matrices.

#### Differentiating quadratic forms

Expand in a direction $`h`$:

```math
(x+h)^\top A(x+h)-x^\top Ax
=h^\top Ax+x^\top Ah+h^\top Ah.
```

The linear term is $`[(A+A^\top)x]^\top h`$ and the remainder is $`O(\Vert h\Vert ^2)`$. Thus the gradient is $`(A+A^\top)x`$ and the Hessian is $`A+A^\top`$. For $`f(x)=\tfrac12x^\top Hx-g^\top x+c`$ with symmetric PD $`H`$,

```math
\nabla f=Hx-g,\quad x^*=H^{-1}g,\quad
f(x)-f(x^*)=\tfrac12(x-x^*)^\top H(x-x^*).
```

The last identity proves global optimality and uniqueness.

#### Equality-constrained quadratic

Minimize $`\tfrac12x^\top Hx-g^\top x`$ subject to $`C^\top x=d`$, where $`H`$ is PD and $`C`$ has independent columns. With multiplier $`\nu`$,

```math
L=\tfrac12x^\top Hx-g^\top x+\nu^\top(C^\top x-d).
```

Stationarity gives $`Hx-g+C\nu=0`$, so $`x=H^{-1}(g-C\nu)`$. Substitution yields

```math
\nu=(C^\top H^{-1}C)^{-1}(C^\top H^{-1}g-d).
```

The inverse exists because $`z^\top C^\top H^{-1}Cz\gt 0`$ for nonzero $`z`$. For a feasible displacement $`h`$ with $`C^\top h=0`$, the linear objective change vanishes at this point and the remaining change is $`\tfrac12h^\top Hh\gt 0`$. This proves uniqueness on the feasible set.

#### Matrix-valued unknowns and the Frobenius gradient

For a matrix variable $`X`$, a differential is a scalar linear function of a perturbation $`H`$. The Frobenius gradient $`G`$ is defined by $`Df(X)[H]=\mathop{\mathrm{tr}} (G^\top H)`$. This convention specifies both the dimensions and transpose placement.

Let $`A\in\mathbb R^{m\times p}`$, $`X\in\mathbb R^{p\times q}`$, $`B\in\mathbb R^{q\times n}`$, and $`C\in\mathbb R^{m\times n}`$. For

```math
f(X)=\frac12\Vert AXB-C\Vert _F^2,
```

put $`R=AXB-C`$. Expanding $`f(X+hH)`$ to first order gives $`Df(X)[H]=\mathop{\mathrm{tr}} (R^\top AHB)`$. Cyclic movement inside the trace yields

```math
\nabla_Xf=A^\top(AXB-C)B^\top.
```

The gradient has size $`p\times q`$, matching $`X`$. Its Hessian acts on a perturbation as $`H\mapsto A^\top AHBB^\top`$. The second directional derivative is $`\Vert AHB\Vert _F^2\ge0`$, proving convexity without vectorizing the matrix.

The objective is strictly convex when the map $`H\mapsto AHB`$ has no nonzero kernel. Full column rank of $`A`$ and full row rank of $`B`$ ensure this. Rank deficiency can leave whole perturbation directions invisible to the objective.

For the separate expression $`f(X)=\mathop{\mathrm{tr}} (X^\top AXB)`$ with square compatible $`A,B`$, differentiating both copies of $`X`$ gives

```math
\nabla_Xf=AXB+A^\top XB^\top.
```

Only when both fixed matrices are symmetric does this simplify to $`2AXB`$. Treating a quadratic-looking expression as though its coefficients were symmetric is a common transpose error.

#### Example LA-E23: a rectangular-objective derivative

Take $`A=\mathop{\mathrm{diag}} (1,2)`$, $`B=\mathop{\mathrm{diag}} (1,3)`$, and let $`C`$ be the two-by-two all-ones matrix. Minimize $`\tfrac12\Vert AXB-C\Vert _F^2`$.

Solution. Both fixed matrices are invertible, so exact fit is possible and unique:

```math
X^*=A^{-1}CB^{-1}
=\begin{pmatrix}1&1/3\\[0pt] 1/2&1/6\end{pmatrix}.
```

Multiplying $`AX^*B`$ gives $`C`$. The residual and the Frobenius gradient are both zero. For any perturbation $`H\ne0`$, invertibility gives $`AHB\ne0`$ and the loss increases by $`\tfrac12\Vert AHB\Vert _F^2`$. This checks uniqueness through the objective itself.

For symmetric positive-definite $`X`$, another useful differential is $`D\log\det X[H]=\mathop{\mathrm{tr}} (X^{-1}H)`$. In a symmetric direction its second derivative is $`-\mathop{\mathrm{tr}} (X^{-1}HX^{-1}H)`$. Writing $`K=X^{-1/2}HX^{-1/2}`$ makes this $`-\mathop{\mathrm{tr}} (K^2)\le0`$, proving concavity of log determinant on the positive-definite cone. The domain is essential: an arbitrary real matrix need not have a positive determinant or a real logarithm.

### 6.3 Worked examples

#### Example LA-E16: differentiating a nonsymmetric expression

Let

```math
A=\left(\begin{smallmatrix}1&3\\[0pt] 0&2\end{smallmatrix}\right)
```

.

Solution. Since $`x^\top Ax=x_1^2+3x_1x_2+2x_2^2`$, the gradient is $`(2x_1+3x_2,3x_1+4x_2)^\top`$. It equals $`(A+A^\top)x`$, not $`2Ax`$. The Hessian has determinant $`8-9=-1`$, so the quadratic is indefinite.

#### Example LA-E17: a budget-constrained minimum

Minimize $`w^\top\Sigma w`$ subject to $`\mathbf1^\top w=1`$, with $`\Sigma`$ PD.

Solution. The Lagrangian $`L=w^\top\Sigma w-\eta(\mathbf1^\top w-1)`$ gives $`2\Sigma w=\eta\mathbf1`$. Enforcing the budget gives

```math
w^*=\frac{\Sigma^{-1}\mathbf1}{\mathbf1^\top\Sigma^{-1}\mathbf1}.
```

The denominator is positive because $`\Sigma^{-1}`$ is PD. Every other feasible $`w=w^*+h`$ has $`\mathbf1^\top h=0`$, and the objective difference is $`h^\top\Sigma h\gt 0`$ unless $`h=0`$. With $`\Sigma=\mathop{\mathrm{diag}} (1,4)`$, the weights are $`(4/5,1/5)`$ and the minimum variance is $`4/5`$.

#### Example LA-E18: a matrix gradient

For $`F(X)=\tfrac12\Vert AX-B\Vert _F^2`$, find $`\nabla_XF`$.

Solution. Write $`R=AX-B`$. Then

```math
dF=\mathop{\mathrm{tr}} (R^\top dR)
=\mathop{\mathrm{tr}} (R^\top A\,dX)
=\mathop{\mathrm{tr}} ([A^\top R]^\top dX).
```

Hence $`\nabla_XF=A^\top(AX-B)`$. Each column is an ordinary least-squares gradient, consistent with the sum-of-column-errors interpretation.

### 6.4 Additional solved problems

These problems are constructed for these notes.

#### Problem LA-17: a singular quadratic objective

For symmetric PSD $`H`$, characterize the minimizers of $`f(x)=\tfrac12x^\top Hx-g^\top x`$.

Solution. Decompose $`g=g_R+g_N`$ into $`\mathop{\mathrm{range}} H`$ and $`\ker H`$. If $`g_N\ne0`$, set $`x=tg_N`$. The quadratic term is zero and $`f(x)=-t\Vert g_N\Vert ^2\to-\infty`$, so no minimum exists. If $`g_N=0`$, then $`g`$ lies in the range and the equation $`Hx=g`$ is solvable. Its solutions are $`H^+g+z`$, $`z\in\ker H`$. Expanding around any such solution gives a nonnegative difference $`\tfrac12h^\top Hh`$. All these solutions are minimizers; uniqueness occurs exactly when $`H`$ is PD.

#### Problem LA-18: a singular covariance with a budget constraint

Set $`\Sigma=\mathop{\mathrm{diag}} (1,0)`$ and minimize $`w^\top\Sigma w`$ under $`w_1+w_2=1`$.

Solution. Substitution gives the objective $`w_1^2`$, uniquely minimized at $`(0,1)`$ with value zero. The inverse-based formula cannot be used because $`\Sigma`$ is singular. Replacing the inverse by $`\Sigma^+`$ in that formula would give $`(1,0)`$, which is wrong. A null direction satisfying the budget can dominate every positive-variance direction. In general solve the constrained problem using its null space or its KKT system, and check curvature along feasible directions.

#### Problem LA-19: an active nonnegativity constraint

Minimize $`\tfrac12[(x_1+1)^2+(x_2-2)^2]`$ subject to $`x_1,x_2\ge0`$.

Solution. The unconstrained minimizer $`(-1,2)`$ is infeasible. Use constraints $`-x_i\le0`$ and multipliers $`\mu_i\ge0`$. Stationarity reads $`x_1+1-\mu_1=0`$ and $`x_2-2-\mu_2=0`$. At $`(0,2)`$ choose $`(\mu_1,\mu_2)=(1,0)`$. Feasibility, stationarity, and $`\mu_ix_i=0`$ all hold. Strict convexity gives the unique global minimum, whose value is $`1/2`$.

#### Problem LA-20: ridge with an unpenalized intercept

Minimize $`\Vert y-a\mathbf1-Xb\Vert ^2+\lambda\Vert b\Vert ^2`$, with $`\lambda\gt 0`$.

Solution. Differentiating in $`a`$ gives $`a=\bar y-\bar x^\top b`$, where $`\bar x^\top=\mathbf1^\top X/n`$. Substituting centers both sides:

```math
\min_b\Vert Hy-HXb\Vert ^2+\lambda\Vert b\Vert ^2,\qquad
b^*=(X^\top HX+\lambda I)^{-1}X^\top Hy.
```

The intercept follows from the earlier equation. Since $`X^\top HX`$ is PSD, adding $`\lambda I`$ makes the slope solution unique. Penalizing the intercept as well would define a different problem and would no longer force the residuals to sum to zero.

### Further problems on matrix differentiation

#### Problem LA-61: weighted ridge regression

Constructed problem. Minimize $`(y-Xb)^\top W(y-Xb)+\lambda b^\top Rb`$, with $`W`$ PD, $`R`$ PSD, and $`\lambda\gt 0`$.

Solution. Differentiate the two quadratic terms:

```math
\nabla f(b)=-2X^\top W(y-Xb)+2\lambda Rb.
```

Stationarity gives $`(X^\top WX+\lambda R)b=X^\top Wy`$. For any $`z`$, the coefficient matrix has quadratic form $`\Vert W^{1/2}Xz\Vert ^2+\lambda z^\top Rz`$. Its kernel is $`\ker X\cap\ker R`$. Thus the solution is unique exactly when this intersection is zero. Merely calling the penalty ridge is insufficient to guarantee uniqueness if $`R`$ leaves some unidentified directions unpenalized.

#### Problem LA-62: sensitivity of a quadratic optimum

Constructed problem. Let $`x^*(g)=H^{-1}g`$ for fixed PD $`H`$. Differentiate the optimum and the minimized value with respect to $`g`$.

Solution. A perturbation $`dg`$ changes the optimum by $`dx^*=H^{-1}dg`$. The minimized value is $`v(g)=c-\tfrac12g^\top H^{-1}g`$, so $`\nabla_gv=-H^{-1}g=-x^*`$. Its Hessian is $`-H^{-1}`$, hence the value is concave in $`g`$. The decision is most sensitive along the smallest-eigenvalue directions of $`H`$, agreeing with the conditioning calculation.

#### Problem LA-63: optimizing a trace expression

Constructed problem. For symmetric PD $`H`$ and a given $`G\in\mathbb R^{n\times k}`$, minimize $`F(X)=\tfrac12\mathop{\mathrm{tr}} (X^\top HX)-\mathop{\mathrm{tr}} (G^\top X)`$.

Solution. Write $`X`$ and $`G`$ in columns. The trace expression is $`\sum_{j=1}^k[\tfrac12x_j^\top Hx_j-g_j^\top x_j]`$. Each column has the unique minimum $`x_j=H^{-1}g_j`$. Therefore $`X^*=H^{-1}G`$. Equivalently the matrix gradient is $`HX-G`$, and

```math
F(X)-F(X^*)=\tfrac12\mathop{\mathrm{tr}} ((X-X^*)^\top H(X-X^*))\ge0,
```

with equality only at $`X=X^*`$.

#### Problem LA-64: a quadratic with two equality constraints

Constructed problem. Minimize $`x_1^2+x_2^2+x_3^2`$ subject to $`x_1+x_2+x_3=1`$ and $`x_1-x_3=0`$.

Solution. The second equality sets $`x_3=x_1=t`$, and the first gives $`x_2=1-2t`$. Substitute:

```math
f(t)=2t^2+(1-2t)^2=6t^2-4t+1.
```

Its derivative is $`12t-4`$, so $`t=1/3`$ and $`x=(1/3,1/3,1/3)`$. The second derivative is $`12\gt 0`$, proving uniqueness on the feasible line. The minimum is $`1/3`$. Eliminating constraints first has reduced a three-variable constrained problem to an ordinary scalar quadratic.

#### Problem LA-65: curvature only on a feasible line

Constructed problem. Minimize $`x_1^2-x_2^2`$ subject to $`x_2=0`$, then subject to $`x_1=0`$.

Solution. Under the first constraint the objective becomes $`x_1^2`$, so zero is the unique minimizer even though the full Hessian is indefinite. Under the second constraint it becomes $`-x_2^2`$, which tends to $`-\infty`$ as $`\vert x_2\vert \to\infty`$. The two feasible tangent spaces select opposite curvature directions. A constrained second-order test must inspect feasible directions, not merely classify the unconstrained Hessian.

#### Problem LA-66: the derivative of an inverse quadratic

Constructed problem. For fixed $`a\ne0`$ and PD $`X`$, differentiate $`f(X)=a^\top X^{-1}a`$.

Solution. The inverse differential gives

```math
df=-a^\top X^{-1}(dX)X^{-1}a
=-\mathop{\mathrm{tr}} (X^{-1}aa^\top X^{-1}dX).
```

Thus $`\nabla_Xf=-X^{-1}aa^\top X^{-1}`$ on symmetric matrices. Along a PSD direction $`H`$, $`df[H]=-(X^{-1}a)^\top H(X^{-1}a)\le0`$. Increasing the matrix in the PSD order can only decrease this inverse quadratic, as long as the path stays PD.

#### Problem LA-67: a log-determinant optimum

Constructed problem. For PD $`S`$, minimize $`F(X)=\mathop{\mathrm{tr}} (SX)-\log\det X`$ over PD $`X`$.

Solution. Stationarity gives $`S-X^{-1}=0`$, so $`X^*=S^{-1}`$. To verify global optimality without assuming it, set $`Y=S^{1/2}XS^{1/2}`$. Then $`\mathop{\mathrm{tr}} (SX)=\mathop{\mathrm{tr}} Y`$ and $`\log\det X=\log\det Y-\log\det S`$. Hence, up to a constant,

```math
F(X)=\sum_i(\lambda_i(Y)-\log\lambda_i(Y)).
```

For a positive scalar $`t`$, $`t-\log t`$ has unique minimum one at $`t=1`$. Thus all eigenvalues of $`Y`$ must equal one, giving $`Y=I`$ and $`X=S^{-1}`$. This proves existence and uniqueness.

#### Problem LA-68: differentiating a normalized quadratic solution

Constructed problem. Let $`w(\Sigma)=\Sigma^{-1}\mathbf1/(\mathbf1^\top\Sigma^{-1}\mathbf1)`$ for PD $`\Sigma`$. Find its first-order change under a symmetric perturbation $`E`$.

Solution. Write $`a=\Sigma^{-1}\mathbf1`$ and $`s=\mathbf1^\top a\gt 0`$. The inverse formula gives $`da=-\Sigma^{-1}Ea`$ and $`ds=-a^\top Ea`$. Applying the scalar quotient rule componentwise,

```math
dw=-\frac{\Sigma^{-1}Ea}{s}
+\frac{a(a^\top Ea)}{s^2}.
```

Multiplying by $`\mathbf1^\top`$ gives zero: the first term sums to $`-a^\top Ea/s`$, and the second to $`a^\top Ea/s`$. Thus the first-order perturbation respects the budget constraint. If $`\Sigma`$ is nearly singular, the inverse factors can make this change large even for small $`E`$.

<!-- extended-practice-start -->

## 7. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 7.1 Weighted projection onto a line

For an ordinary Euclidean projection, all coordinate errors receive the same weight. With a positive diagonal matrix $`W`$, the loss is $`(y-\hat y)^\top W(y-\hat y)`$. A large diagonal entry makes the corresponding error more expensive. The fitted vector still lies in the chosen subspace, but the residual is orthogonal using the $`W`$ inner product. These problems compute both objects and check that orthogonality.

#### Problem LA-69: weighted line fit, case 1

Let

```math
x=\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

,

```math
y=\begin{pmatrix}4\\[0pt] 1\end{pmatrix}
```

, and

```math
W=\begin{pmatrix}3&0\\[0pt] 0&1\end{pmatrix}
```

. Find the scalar $`a_*`$ minimizing $`g(a)=(y-ax)^\top W(y-ax)`$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

```math
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
```

Here $`x^\top Wx=7\gt 0`$ and $`x^\top Wy=14`$. The derivative is $`g'(a)=-2x^\top Wy+2a x^\top Wx`$, so

```math
a_*=\frac{x^\top Wy}{x^\top Wx}=2,\qquad
\hat y=a_*x=\begin{pmatrix}2\\[0pt] 4\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}2\\[0pt] -3\end{pmatrix}.
```

The second derivative is $`2x^\top Wx=14\gt 0`$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $`x`$.

For an arbitrary input vector $`z`$, the same derivation gives $`P_Wz=x(x^\top Wz)/(x^\top Wx)`$. Therefore

```math
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{3}{7}&\frac{2}{7}\\[0pt] \frac{6}{7}&\frac{4}{7}\end{pmatrix}.
```

Multiplying this matrix by $`y`$ recovers the fitted vector above. Also $`P_Wx=x`$, so applying the projection twice changes nothing: $`P_W^2=P_W`$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $`P_W^\top W=WP_W`$.

The residual check is $`x^\top Wr=x^\top Wy-a_*x^\top Wx=0`$. The minimum loss is

```math
r^\top Wr=21.
```

Finally, writing $`a=a_*+h`$ and expanding the loss shows $`g(a)-g(a_*)=7h^2`$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-70: weighted line fit, case 2

Let

```math
x=\begin{pmatrix}2\\[0pt] -1\end{pmatrix}
```

,

```math
y=\begin{pmatrix}1\\[0pt] 3\end{pmatrix}
```

, and

```math
W=\begin{pmatrix}1&0\\[0pt] 0&4\end{pmatrix}
```

. Find the scalar $`a_*`$ minimizing $`g(a)=(y-ax)^\top W(y-ax)`$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

```math
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
```

Here $`x^\top Wx=8\gt 0`$ and $`x^\top Wy=-10`$. The derivative is $`g'(a)=-2x^\top Wy+2a x^\top Wx`$, so

```math
a_*=\frac{x^\top Wy}{x^\top Wx}=-\frac{5}{4},\qquad
\hat y=a_*x=\begin{pmatrix}-\frac{5}{2}\\[0pt] \frac{5}{4}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{7}{2}\\[0pt] \frac{7}{4}\end{pmatrix}.
```

The second derivative is $`2x^\top Wx=16\gt 0`$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $`x`$.

For an arbitrary input vector $`z`$, the same derivation gives $`P_Wz=x(x^\top Wz)/(x^\top Wx)`$. Therefore

```math
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{1}{2}&-1\\[0pt] -\frac{1}{4}&\frac{1}{2}\end{pmatrix}.
```

Multiplying this matrix by $`y`$ recovers the fitted vector above. Also $`P_Wx=x`$, so applying the projection twice changes nothing: $`P_W^2=P_W`$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $`P_W^\top W=WP_W`$.

The residual check is $`x^\top Wr=x^\top Wy-a_*x^\top Wx=0`$. The minimum loss is

```math
r^\top Wr=\frac{49}{2}.
```

Finally, writing $`a=a_*+h`$ and expanding the loss shows $`g(a)-g(a_*)=8h^2`$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-71: weighted line fit, case 3

Let

```math
x=\begin{pmatrix}1\\[0pt] 1\end{pmatrix}
```

,

```math
y=\begin{pmatrix}-1\\[0pt] 4\end{pmatrix}
```

, and

```math
W=\begin{pmatrix}2&0\\[0pt] 0&5\end{pmatrix}
```

. Find the scalar $`a_*`$ minimizing $`g(a)=(y-ax)^\top W(y-ax)`$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

```math
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
```

Here $`x^\top Wx=7\gt 0`$ and $`x^\top Wy=18`$. The derivative is $`g'(a)=-2x^\top Wy+2a x^\top Wx`$, so

```math
a_*=\frac{x^\top Wy}{x^\top Wx}=\frac{18}{7},\qquad
\hat y=a_*x=\begin{pmatrix}\frac{18}{7}\\[0pt] \frac{18}{7}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}-\frac{25}{7}\\[0pt] \frac{10}{7}\end{pmatrix}.
```

The second derivative is $`2x^\top Wx=14\gt 0`$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $`x`$.

For an arbitrary input vector $`z`$, the same derivation gives $`P_Wz=x(x^\top Wz)/(x^\top Wx)`$. Therefore

```math
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{2}{7}&\frac{5}{7}\\[0pt] \frac{2}{7}&\frac{5}{7}\end{pmatrix}.
```

Multiplying this matrix by $`y`$ recovers the fitted vector above. Also $`P_Wx=x`$, so applying the projection twice changes nothing: $`P_W^2=P_W`$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $`P_W^\top W=WP_W`$.

The residual check is $`x^\top Wr=x^\top Wy-a_*x^\top Wx=0`$. The minimum loss is

```math
r^\top Wr=\frac{250}{7}.
```

Finally, writing $`a=a_*+h`$ and expanding the loss shows $`g(a)-g(a_*)=7h^2`$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-72: weighted line fit, case 4

Let

```math
x=\begin{pmatrix}3\\[0pt] 1\end{pmatrix}
```

,

```math
y=\begin{pmatrix}2\\[0pt] -2\end{pmatrix}
```

, and

```math
W=\begin{pmatrix}1&0\\[0pt] 0&2\end{pmatrix}
```

. Find the scalar $`a_*`$ minimizing $`g(a)=(y-ax)^\top W(y-ax)`$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

```math
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
```

Here $`x^\top Wx=11\gt 0`$ and $`x^\top Wy=2`$. The derivative is $`g'(a)=-2x^\top Wy+2a x^\top Wx`$, so

```math
a_*=\frac{x^\top Wy}{x^\top Wx}=\frac{2}{11},\qquad
\hat y=a_*x=\begin{pmatrix}\frac{6}{11}\\[0pt] \frac{2}{11}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{16}{11}\\[0pt] -\frac{24}{11}\end{pmatrix}.
```

The second derivative is $`2x^\top Wx=22\gt 0`$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $`x`$.

For an arbitrary input vector $`z`$, the same derivation gives $`P_Wz=x(x^\top Wz)/(x^\top Wx)`$. Therefore

```math
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{9}{11}&\frac{6}{11}\\[0pt] \frac{3}{11}&\frac{2}{11}\end{pmatrix}.
```

Multiplying this matrix by $`y`$ recovers the fitted vector above. Also $`P_Wx=x`$, so applying the projection twice changes nothing: $`P_W^2=P_W`$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $`P_W^\top W=WP_W`$.

The residual check is $`x^\top Wr=x^\top Wy-a_*x^\top Wx=0`$. The minimum loss is

```math
r^\top Wr=\frac{128}{11}.
```

Finally, writing $`a=a_*+h`$ and expanding the loss shows $`g(a)-g(a_*)=11h^2`$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-73: weighted line fit, case 5

Let

```math
x=\begin{pmatrix}1\\[0pt] -2\end{pmatrix}
```

,

```math
y=\begin{pmatrix}3\\[0pt] 2\end{pmatrix}
```

, and

```math
W=\begin{pmatrix}4&0\\[0pt] 0&1\end{pmatrix}
```

. Find the scalar $`a_*`$ minimizing $`g(a)=(y-ax)^\top W(y-ax)`$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

```math
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
```

Here $`x^\top Wx=8\gt 0`$ and $`x^\top Wy=8`$. The derivative is $`g'(a)=-2x^\top Wy+2a x^\top Wx`$, so

```math
a_*=\frac{x^\top Wy}{x^\top Wx}=1,\qquad
\hat y=a_*x=\begin{pmatrix}1\\[0pt] -2\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}2\\[0pt] 4\end{pmatrix}.
```

The second derivative is $`2x^\top Wx=16\gt 0`$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $`x`$.

For an arbitrary input vector $`z`$, the same derivation gives $`P_Wz=x(x^\top Wz)/(x^\top Wx)`$. Therefore

```math
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{1}{2}&-\frac{1}{4}\\[0pt] -1&\frac{1}{2}\end{pmatrix}.
```

Multiplying this matrix by $`y`$ recovers the fitted vector above. Also $`P_Wx=x`$, so applying the projection twice changes nothing: $`P_W^2=P_W`$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $`P_W^\top W=WP_W`$.

The residual check is $`x^\top Wr=x^\top Wy-a_*x^\top Wx=0`$. The minimum loss is

```math
r^\top Wr=32.
```

Finally, writing $`a=a_*+h`$ and expanding the loss shows $`g(a)-g(a_*)=8h^2`$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-74: weighted line fit, case 6

Let

```math
x=\begin{pmatrix}2\\[0pt] 3\end{pmatrix}
```

,

```math
y=\begin{pmatrix}-1\\[0pt] 5\end{pmatrix}
```

, and

```math
W=\begin{pmatrix}3&0\\[0pt] 0&2\end{pmatrix}
```

. Find the scalar $`a_*`$ minimizing $`g(a)=(y-ax)^\top W(y-ax)`$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

```math
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
```

Here $`x^\top Wx=30\gt 0`$ and $`x^\top Wy=24`$. The derivative is $`g'(a)=-2x^\top Wy+2a x^\top Wx`$, so

```math
a_*=\frac{x^\top Wy}{x^\top Wx}=\frac{4}{5},\qquad
\hat y=a_*x=\begin{pmatrix}\frac{8}{5}\\[0pt] \frac{12}{5}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}-\frac{13}{5}\\[0pt] \frac{13}{5}\end{pmatrix}.
```

The second derivative is $`2x^\top Wx=60\gt 0`$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $`x`$.

For an arbitrary input vector $`z`$, the same derivation gives $`P_Wz=x(x^\top Wz)/(x^\top Wx)`$. Therefore

```math
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{2}{5}&\frac{2}{5}\\[0pt] \frac{3}{5}&\frac{3}{5}\end{pmatrix}.
```

Multiplying this matrix by $`y`$ recovers the fitted vector above. Also $`P_Wx=x`$, so applying the projection twice changes nothing: $`P_W^2=P_W`$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $`P_W^\top W=WP_W`$.

The residual check is $`x^\top Wr=x^\top Wy-a_*x^\top Wx=0`$. The minimum loss is

```math
r^\top Wr=\frac{169}{5}.
```

Finally, writing $`a=a_*+h`$ and expanding the loss shows $`g(a)-g(a_*)=30h^2`$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

### 7.2 Least squares with a visible residual

At the design points $`-1,0,1`$, the intercept column and slope column are orthogonal. This makes the normal equations diagonal without changing the meaning of least squares. The middle observation is especially useful: a straight line can fit all three observations exactly only if the second difference $`y_0-2y_1+y_2`$ vanishes.

#### Problem LA-75: three-observation regression, case 1

Fit $`y_i=a+bt_i+r_i`$ at $`t=(-1,0,1)`$ to

```math
y=\begin{pmatrix}1\\[0pt] 2\\[0pt] 5\end{pmatrix}
```

 by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

```math
X=\begin{pmatrix}1&-1\\[0pt] 1&0\\[0pt] 1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\[0pt] 0&2\end{pmatrix}.
```

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $`2X^\top X`$. The normal equations consequently have one solution:

```math
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
```

Substitution gives $`a=\frac{8}{3}`$ and $`b=2`$. Evaluating the fitted line at each of the original design points gives

```math
\hat y=\begin{pmatrix}\frac{2}{3}\\[0pt] \frac{8}{3}\\[0pt] \frac{14}{3}\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}\frac{1}{3}\\[0pt] -\frac{2}{3}\\[0pt] \frac{1}{3}\end{pmatrix}.
```

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $`-r_0+r_2=0`$. These are the two entries of $`X^\top r=0`$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $`d=y_0-2y_1+y_2=2`$. Direct substitution into the fitted values gives

```math
r=\frac d6\begin{pmatrix}1\\[0pt] -2\\[0pt] 1\end{pmatrix},\qquad
\mathop{\mathrm{SSE}} =\frac{d^2}6=\frac{2}{3}.
```

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $`d`$.

For any proposed change $`h`$ in the two coefficients, orthogonality gives $`\Vert y-X(\hat\beta+h)\Vert _2^2=\Vert r\Vert _2^2+\Vert Xh\Vert _2^2`$. The second term is strictly positive for nonzero $`h`$ because $`X`$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-76: three-observation regression, case 2

Fit $`y_i=a+bt_i+r_i`$ at $`t=(-1,0,1)`$ to

```math
y=\begin{pmatrix}4\\[0pt] 1\\[0pt] 2\end{pmatrix}
```

 by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

```math
X=\begin{pmatrix}1&-1\\[0pt] 1&0\\[0pt] 1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\[0pt] 0&2\end{pmatrix}.
```

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $`2X^\top X`$. The normal equations consequently have one solution:

```math
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
```

Substitution gives $`a=\frac{7}{3}`$ and $`b=-1`$. Evaluating the fitted line at each of the original design points gives

```math
\hat y=\begin{pmatrix}\frac{10}{3}\\[0pt] \frac{7}{3}\\[0pt] \frac{4}{3}\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}\frac{2}{3}\\[0pt] -\frac{4}{3}\\[0pt] \frac{2}{3}\end{pmatrix}.
```

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $`-r_0+r_2=0`$. These are the two entries of $`X^\top r=0`$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $`d=y_0-2y_1+y_2=4`$. Direct substitution into the fitted values gives

```math
r=\frac d6\begin{pmatrix}1\\[0pt] -2\\[0pt] 1\end{pmatrix},\qquad
\mathop{\mathrm{SSE}} =\frac{d^2}6=\frac{8}{3}.
```

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $`d`$.

For any proposed change $`h`$ in the two coefficients, orthogonality gives $`\Vert y-X(\hat\beta+h)\Vert _2^2=\Vert r\Vert _2^2+\Vert Xh\Vert _2^2`$. The second term is strictly positive for nonzero $`h`$ because $`X`$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-77: three-observation regression, case 3

Fit $`y_i=a+bt_i+r_i`$ at $`t=(-1,0,1)`$ to

```math
y=\begin{pmatrix}-2\\[0pt] 0\\[0pt] 4\end{pmatrix}
```

 by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

```math
X=\begin{pmatrix}1&-1\\[0pt] 1&0\\[0pt] 1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\[0pt] 0&2\end{pmatrix}.
```

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $`2X^\top X`$. The normal equations consequently have one solution:

```math
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
```

Substitution gives $`a=\frac{2}{3}`$ and $`b=3`$. Evaluating the fitted line at each of the original design points gives

```math
\hat y=\begin{pmatrix}-\frac{7}{3}\\[0pt] \frac{2}{3}\\[0pt] \frac{11}{3}\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}\frac{1}{3}\\[0pt] -\frac{2}{3}\\[0pt] \frac{1}{3}\end{pmatrix}.
```

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $`-r_0+r_2=0`$. These are the two entries of $`X^\top r=0`$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $`d=y_0-2y_1+y_2=2`$. Direct substitution into the fitted values gives

```math
r=\frac d6\begin{pmatrix}1\\[0pt] -2\\[0pt] 1\end{pmatrix},\qquad
\mathop{\mathrm{SSE}} =\frac{d^2}6=\frac{2}{3}.
```

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $`d`$.

For any proposed change $`h`$ in the two coefficients, orthogonality gives $`\Vert y-X(\hat\beta+h)\Vert _2^2=\Vert r\Vert _2^2+\Vert Xh\Vert _2^2`$. The second term is strictly positive for nonzero $`h`$ because $`X`$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-78: three-observation regression, case 4

Fit $`y_i=a+bt_i+r_i`$ at $`t=(-1,0,1)`$ to

```math
y=\begin{pmatrix}3\\[0pt] 3\\[0pt] 3\end{pmatrix}
```

 by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

```math
X=\begin{pmatrix}1&-1\\[0pt] 1&0\\[0pt] 1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\[0pt] 0&2\end{pmatrix}.
```

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $`2X^\top X`$. The normal equations consequently have one solution:

```math
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
```

Substitution gives $`a=3`$ and $`b=0`$. Evaluating the fitted line at each of the original design points gives

```math
\hat y=\begin{pmatrix}3\\[0pt] 3\\[0pt] 3\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}0\\[0pt] 0\\[0pt] 0\end{pmatrix}.
```

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $`-r_0+r_2=0`$. These are the two entries of $`X^\top r=0`$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $`d=y_0-2y_1+y_2=0`$. Direct substitution into the fitted values gives

```math
r=\frac d6\begin{pmatrix}1\\[0pt] -2\\[0pt] 1\end{pmatrix},\qquad
\mathop{\mathrm{SSE}} =\frac{d^2}6=0.
```

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $`d`$.

For any proposed change $`h`$ in the two coefficients, orthogonality gives $`\Vert y-X(\hat\beta+h)\Vert _2^2=\Vert r\Vert _2^2+\Vert Xh\Vert _2^2`$. The second term is strictly positive for nonzero $`h`$ because $`X`$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-79: three-observation regression, case 5

Fit $`y_i=a+bt_i+r_i`$ at $`t=(-1,0,1)`$ to

```math
y=\begin{pmatrix}2\\[0pt] -1\\[0pt] 2\end{pmatrix}
```

 by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

```math
X=\begin{pmatrix}1&-1\\[0pt] 1&0\\[0pt] 1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\[0pt] 0&2\end{pmatrix}.
```

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $`2X^\top X`$. The normal equations consequently have one solution:

```math
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
```

Substitution gives $`a=1`$ and $`b=0`$. Evaluating the fitted line at each of the original design points gives

```math
\hat y=\begin{pmatrix}1\\[0pt] 1\\[0pt] 1\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}1\\[0pt] -2\\[0pt] 1\end{pmatrix}.
```

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $`-r_0+r_2=0`$. These are the two entries of $`X^\top r=0`$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $`d=y_0-2y_1+y_2=6`$. Direct substitution into the fitted values gives

```math
r=\frac d6\begin{pmatrix}1\\[0pt] -2\\[0pt] 1\end{pmatrix},\qquad
\mathop{\mathrm{SSE}} =\frac{d^2}6=6.
```

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $`d`$.

For any proposed change $`h`$ in the two coefficients, orthogonality gives $`\Vert y-X(\hat\beta+h)\Vert _2^2=\Vert r\Vert _2^2+\Vert Xh\Vert _2^2`$. The second term is strictly positive for nonzero $`h`$ because $`X`$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-80: three-observation regression, case 6

Fit $`y_i=a+bt_i+r_i`$ at $`t=(-1,0,1)`$ to

```math
y=\begin{pmatrix}5\\[0pt] 0\\[0pt] -3\end{pmatrix}
```

 by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

```math
X=\begin{pmatrix}1&-1\\[0pt] 1&0\\[0pt] 1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\[0pt] 0&2\end{pmatrix}.
```

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $`2X^\top X`$. The normal equations consequently have one solution:

```math
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
```

Substitution gives $`a=\frac{2}{3}`$ and $`b=-4`$. Evaluating the fitted line at each of the original design points gives

```math
\hat y=\begin{pmatrix}\frac{14}{3}\\[0pt] \frac{2}{3}\\[0pt] -\frac{10}{3}\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}\frac{1}{3}\\[0pt] -\frac{2}{3}\\[0pt] \frac{1}{3}\end{pmatrix}.
```

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $`-r_0+r_2=0`$. These are the two entries of $`X^\top r=0`$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $`d=y_0-2y_1+y_2=2`$. Direct substitution into the fitted values gives

```math
r=\frac d6\begin{pmatrix}1\\[0pt] -2\\[0pt] 1\end{pmatrix},\qquad
\mathop{\mathrm{SSE}} =\frac{d^2}6=\frac{2}{3}.
```

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $`d`$.

For any proposed change $`h`$ in the two coefficients, orthogonality gives $`\Vert y-X(\hat\beta+h)\Vert _2^2=\Vert r\Vert _2^2+\Vert Xh\Vert _2^2`$. The second term is strictly positive for nonzero $`h`$ because $`X`$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

### 7.3 Common and opposite directions in a quadratic form

For a symmetric matrix with equal diagonal entries, the directions $`(1,1)`$ and $`(1,-1)`$ separate common movement from opposite movement. A covariance matrix must be nonnegative in both directions. A matrix can have positive diagonal entries and still fail this requirement.

#### Problem LA-81: equal-diagonal quadratic form, case 1

Diagonalize

```math
A=\begin{pmatrix}5&3\\[0pt] 3&5\end{pmatrix}
```

. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

```math
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] 1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] -1\end{pmatrix}.
```

Their inner product is zero and each has norm one. Multiplication gives $`Aq_+=8q_+`$ and $`Aq_-=2q_-`$. Hence, with $`Q=(q_+,q_-)`$,

```math
A=Q\begin{pmatrix}8&0\\[0pt] 0&2\end{pmatrix}Q^\top.
```

Every vector has a unique expansion $`x=u q_++w q_-`$, where $`u=(x_1+x_2)/\sqrt2`$ and $`w=(x_1-x_2)/\sqrt2`$. The quadratic form becomes

```math
x^\top Ax=8u^2+(2)w^2.
```

This expression proves that $`A`$ is positive definite. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

Both eigenvalues are positive, so the inverse exists and is

```math
\begin{pmatrix}\frac{5}{16}&-\frac{3}{16}\\[0pt] -\frac{3}{16}&\frac{5}{16}\end{pmatrix}
```

.

As independent checks, the trace is $`10=8+(2)`$ and the determinant is $`16=(8)(2)`$. For nonzero $`x`$, the Rayleigh quotient is a weighted average of the two eigenvalues:

```math
\frac{x^\top Ax}{x^\top x}
=\frac{8u^2+(2)w^2}{u^2+w^2}.
```

It lies between $`2`$ and $`8`$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-82: equal-diagonal quadratic form, case 2

Diagonalize

```math
A=\begin{pmatrix}2&2\\[0pt] 2&2\end{pmatrix}
```

. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

```math
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] 1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] -1\end{pmatrix}.
```

Their inner product is zero and each has norm one. Multiplication gives $`Aq_+=4q_+`$ and $`Aq_-=0q_-`$. Hence, with $`Q=(q_+,q_-)`$,

```math
A=Q\begin{pmatrix}4&0\\[0pt] 0&0\end{pmatrix}Q^\top.
```

Every vector has a unique expansion $`x=u q_++w q_-`$, where $`u=(x_1+x_2)/\sqrt2`$ and $`w=(x_1-x_2)/\sqrt2`$. The quadratic form becomes

```math
x^\top Ax=4u^2+(0)w^2.
```

This expression proves that $`A`$ is positive semidefinite and singular. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

The zero eigenvalue gives a nonzero direction with zero quadratic cost. The matrix has no inverse; a covariance interpretation is possible, with an exact linear relation between the two centered coordinates.

As independent checks, the trace is $`4=4+(0)`$ and the determinant is $`0=(4)(0)`$. For nonzero $`x`$, the Rayleigh quotient is a weighted average of the two eigenvalues:

```math
\frac{x^\top Ax}{x^\top x}
=\frac{4u^2+(0)w^2}{u^2+w^2}.
```

It lies between $`0`$ and $`4`$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-83: equal-diagonal quadratic form, case 3

Diagonalize

```math
A=\begin{pmatrix}2&3\\[0pt] 3&2\end{pmatrix}
```

. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

```math
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] 1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] -1\end{pmatrix}.
```

Their inner product is zero and each has norm one. Multiplication gives $`Aq_+=5q_+`$ and $`Aq_-=-1q_-`$. Hence, with $`Q=(q_+,q_-)`$,

```math
A=Q\begin{pmatrix}5&0\\[0pt] 0&-1\end{pmatrix}Q^\top.
```

Every vector has a unique expansion $`x=u q_++w q_-`$, where $`u=(x_1+x_2)/\sqrt2`$ and $`w=(x_1-x_2)/\sqrt2`$. The quadratic form becomes

```math
x^\top Ax=5u^2+(-1)w^2.
```

This expression proves that $`A`$ is indefinite. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

The negative eigenvalue gives a direction with negative quadratic cost. Such a matrix cannot be a covariance matrix, since a variance cannot be negative.

As independent checks, the trace is $`4=5+(-1)`$ and the determinant is $`-5=(5)(-1)`$. For nonzero $`x`$, the Rayleigh quotient is a weighted average of the two eigenvalues:

```math
\frac{x^\top Ax}{x^\top x}
=\frac{5u^2+(-1)w^2}{u^2+w^2}.
```

It lies between $`-1`$ and $`5`$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-84: equal-diagonal quadratic form, case 4

Diagonalize

```math
A=\begin{pmatrix}4&-3\\[0pt] -3&4\end{pmatrix}
```

. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

```math
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] 1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] -1\end{pmatrix}.
```

Their inner product is zero and each has norm one. Multiplication gives $`Aq_+=1q_+`$ and $`Aq_-=7q_-`$. Hence, with $`Q=(q_+,q_-)`$,

```math
A=Q\begin{pmatrix}1&0\\[0pt] 0&7\end{pmatrix}Q^\top.
```

Every vector has a unique expansion $`x=u q_++w q_-`$, where $`u=(x_1+x_2)/\sqrt2`$ and $`w=(x_1-x_2)/\sqrt2`$. The quadratic form becomes

```math
x^\top Ax=1u^2+(7)w^2.
```

This expression proves that $`A`$ is positive definite. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

Both eigenvalues are positive, so the inverse exists and is

```math
\begin{pmatrix}\frac{4}{7}&\frac{3}{7}\\[0pt] \frac{3}{7}&\frac{4}{7}\end{pmatrix}
```

.

As independent checks, the trace is $`8=1+(7)`$ and the determinant is $`7=(1)(7)`$. For nonzero $`x`$, the Rayleigh quotient is a weighted average of the two eigenvalues:

```math
\frac{x^\top Ax}{x^\top x}
=\frac{1u^2+(7)w^2}{u^2+w^2}.
```

It lies between $`1`$ and $`7`$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-85: equal-diagonal quadratic form, case 5

Diagonalize

```math
A=\begin{pmatrix}3&0\\[0pt] 0&3\end{pmatrix}
```

. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

```math
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] 1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] -1\end{pmatrix}.
```

Their inner product is zero and each has norm one. Multiplication gives $`Aq_+=3q_+`$ and $`Aq_-=3q_-`$. Hence, with $`Q=(q_+,q_-)`$,

```math
A=Q\begin{pmatrix}3&0\\[0pt] 0&3\end{pmatrix}Q^\top.
```

Every vector has a unique expansion $`x=u q_++w q_-`$, where $`u=(x_1+x_2)/\sqrt2`$ and $`w=(x_1-x_2)/\sqrt2`$. The quadratic form becomes

```math
x^\top Ax=3u^2+(3)w^2.
```

This expression proves that $`A`$ is positive definite. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

Both eigenvalues are positive, so the inverse exists and is

```math
\begin{pmatrix}\frac{1}{3}&0\\[0pt] 0&\frac{1}{3}\end{pmatrix}
```

.

As independent checks, the trace is $`6=3+(3)`$ and the determinant is $`9=(3)(3)`$. For nonzero $`x`$, the Rayleigh quotient is a weighted average of the two eigenvalues:

```math
\frac{x^\top Ax}{x^\top x}
=\frac{3u^2+(3)w^2}{u^2+w^2}.
```

It lies between $`3`$ and $`3`$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-86: equal-diagonal quadratic form, case 6

Diagonalize

```math
A=\begin{pmatrix}1&-1\\[0pt] -1&1\end{pmatrix}
```

. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

```math
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] 1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\[0pt] -1\end{pmatrix}.
```

Their inner product is zero and each has norm one. Multiplication gives $`Aq_+=0q_+`$ and $`Aq_-=2q_-`$. Hence, with $`Q=(q_+,q_-)`$,

```math
A=Q\begin{pmatrix}0&0\\[0pt] 0&2\end{pmatrix}Q^\top.
```

Every vector has a unique expansion $`x=u q_++w q_-`$, where $`u=(x_1+x_2)/\sqrt2`$ and $`w=(x_1-x_2)/\sqrt2`$. The quadratic form becomes

```math
x^\top Ax=0u^2+(2)w^2.
```

This expression proves that $`A`$ is positive semidefinite and singular. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

The zero eigenvalue gives a nonzero direction with zero quadratic cost. The matrix has no inverse; a covariance interpretation is possible, with an exact linear relation between the two centered coordinates.

As independent checks, the trace is $`2=0+(2)`$ and the determinant is $`0=(0)(2)`$. For nonzero $`x`$, the Rayleigh quotient is a weighted average of the two eigenvalues:

```math
\frac{x^\top Ax}{x^\top x}
=\frac{0u^2+(2)w^2}{u^2+w^2}.
```

It lies between $`0`$ and $`2`$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

### 7.4 Rank-one least squares and minimum norm

A rank-one matrix can produce only multiples of one vector. Fitting data first projects the target onto that output line. Recovering coefficients then leaves an entire null-space direction undetermined. The pseudoinverse chooses the coefficient vector perpendicular to that null space.

#### Problem LA-87: rank-one pseudoinverse, case 1

Let

```math
a=\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

,

```math
b=\begin{pmatrix}1\\[0pt] 1\end{pmatrix}
```

, $`A=ab^\top`$, and

```math
y=\begin{pmatrix}3\\[0pt] 1\end{pmatrix}
```

. Find every least-squares solution of $`Ax\approx y`$ and identify the one with smallest Euclidean norm.

Solution. Since $`Ax=a(b^\top x)`$, the scalar $`b^\top x`$ determines the entire output. The column space of

```math
A=\begin{pmatrix}1&1\\[0pt] 2&2\end{pmatrix}
```

 is the line spanned by $`a`$. Projecting $`y`$ onto this line gives

```math
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}1\\[0pt] 2\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}2\\[0pt] -1\end{pmatrix}.
```

Here $`a^\top a=5`$, $`b^\top b=2`$, and $`a^\top y=5`$. The check $`a^\top r=0`$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $`b^\top x=1`$.

The vector parallel to $`b`$ satisfying that equation is

```math
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{1}{2}\\[0pt] \frac{1}{2}\end{pmatrix}.
```

Consequently the pseudoinverse is

```math
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{1}{10}&\frac{1}{5}\\[0pt] \frac{1}{10}&\frac{1}{5}\end{pmatrix}.
```

The null space consists of vectors perpendicular to $`b`$. In two dimensions it is spanned by

```math
z=\begin{pmatrix}-1\\[0pt] 1\end{pmatrix}
```

. All least-squares solutions are thus $`x=x_++tz`$ for real $`t`$. They give exactly the same fitted vector because $`Az=0`$.

The minimum-norm claim requires one further check. Since $`x_+`$ is parallel to $`b`$ while $`z`$ is perpendicular to $`b`$, their inner product is zero. Hence

```math
\Vert x_++tz\Vert _2^2=\Vert x_+\Vert _2^2+t^2\Vert z\Vert _2^2.
```

The minimum occurs uniquely at $`t=0`$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-88: rank-one pseudoinverse, case 2

Let

```math
a=\begin{pmatrix}2\\[0pt] -1\end{pmatrix}
```

,

```math
b=\begin{pmatrix}1\\[0pt] -2\end{pmatrix}
```

, $`A=ab^\top`$, and

```math
y=\begin{pmatrix}1\\[0pt] 4\end{pmatrix}
```

. Find every least-squares solution of $`Ax\approx y`$ and identify the one with smallest Euclidean norm.

Solution. Since $`Ax=a(b^\top x)`$, the scalar $`b^\top x`$ determines the entire output. The column space of

```math
A=\begin{pmatrix}2&-4\\[0pt] -1&2\end{pmatrix}
```

 is the line spanned by $`a`$. Projecting $`y`$ onto this line gives

```math
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}-\frac{4}{5}\\[0pt] \frac{2}{5}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{9}{5}\\[0pt] \frac{18}{5}\end{pmatrix}.
```

Here $`a^\top a=5`$, $`b^\top b=5`$, and $`a^\top y=-2`$. The check $`a^\top r=0`$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $`b^\top x=-\frac{2}{5}`$.

The vector parallel to $`b`$ satisfying that equation is

```math
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}-\frac{2}{25}\\[0pt] \frac{4}{25}\end{pmatrix}.
```

Consequently the pseudoinverse is

```math
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{2}{25}&-\frac{1}{25}\\[0pt] -\frac{4}{25}&\frac{2}{25}\end{pmatrix}.
```

The null space consists of vectors perpendicular to $`b`$. In two dimensions it is spanned by

```math
z=\begin{pmatrix}2\\[0pt] 1\end{pmatrix}
```

. All least-squares solutions are thus $`x=x_++tz`$ for real $`t`$. They give exactly the same fitted vector because $`Az=0`$.

The minimum-norm claim requires one further check. Since $`x_+`$ is parallel to $`b`$ while $`z`$ is perpendicular to $`b`$, their inner product is zero. Hence

```math
\Vert x_++tz\Vert _2^2=\Vert x_+\Vert _2^2+t^2\Vert z\Vert _2^2.
```

The minimum occurs uniquely at $`t=0`$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-89: rank-one pseudoinverse, case 3

Let

```math
a=\begin{pmatrix}1\\[0pt] 3\end{pmatrix}
```

,

```math
b=\begin{pmatrix}2\\[0pt] 1\end{pmatrix}
```

, $`A=ab^\top`$, and

```math
y=\begin{pmatrix}-1\\[0pt] 2\end{pmatrix}
```

. Find every least-squares solution of $`Ax\approx y`$ and identify the one with smallest Euclidean norm.

Solution. Since $`Ax=a(b^\top x)`$, the scalar $`b^\top x`$ determines the entire output. The column space of

```math
A=\begin{pmatrix}2&1\\[0pt] 6&3\end{pmatrix}
```

 is the line spanned by $`a`$. Projecting $`y`$ onto this line gives

```math
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}\frac{1}{2}\\[0pt] \frac{3}{2}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}-\frac{3}{2}\\[0pt] \frac{1}{2}\end{pmatrix}.
```

Here $`a^\top a=10`$, $`b^\top b=5`$, and $`a^\top y=5`$. The check $`a^\top r=0`$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $`b^\top x=\frac{1}{2}`$.

The vector parallel to $`b`$ satisfying that equation is

```math
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{1}{5}\\[0pt] \frac{1}{10}\end{pmatrix}.
```

Consequently the pseudoinverse is

```math
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{1}{25}&\frac{3}{25}\\[0pt] \frac{1}{50}&\frac{3}{50}\end{pmatrix}.
```

The null space consists of vectors perpendicular to $`b`$. In two dimensions it is spanned by

```math
z=\begin{pmatrix}-1\\[0pt] 2\end{pmatrix}
```

. All least-squares solutions are thus $`x=x_++tz`$ for real $`t`$. They give exactly the same fitted vector because $`Az=0`$.

The minimum-norm claim requires one further check. Since $`x_+`$ is parallel to $`b`$ while $`z`$ is perpendicular to $`b`$, their inner product is zero. Hence

```math
\Vert x_++tz\Vert _2^2=\Vert x_+\Vert _2^2+t^2\Vert z\Vert _2^2.
```

The minimum occurs uniquely at $`t=0`$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-90: rank-one pseudoinverse, case 4

Let

```math
a=\begin{pmatrix}3\\[0pt] 1\end{pmatrix}
```

,

```math
b=\begin{pmatrix}1\\[0pt] -1\end{pmatrix}
```

, $`A=ab^\top`$, and

```math
y=\begin{pmatrix}2\\[0pt] 5\end{pmatrix}
```

. Find every least-squares solution of $`Ax\approx y`$ and identify the one with smallest Euclidean norm.

Solution. Since $`Ax=a(b^\top x)`$, the scalar $`b^\top x`$ determines the entire output. The column space of

```math
A=\begin{pmatrix}3&-3\\[0pt] 1&-1\end{pmatrix}
```

 is the line spanned by $`a`$. Projecting $`y`$ onto this line gives

```math
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}\frac{33}{10}\\[0pt] \frac{11}{10}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}-\frac{13}{10}\\[0pt] \frac{39}{10}\end{pmatrix}.
```

Here $`a^\top a=10`$, $`b^\top b=2`$, and $`a^\top y=11`$. The check $`a^\top r=0`$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $`b^\top x=\frac{11}{10}`$.

The vector parallel to $`b`$ satisfying that equation is

```math
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{11}{20}\\[0pt] -\frac{11}{20}\end{pmatrix}.
```

Consequently the pseudoinverse is

```math
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{3}{20}&\frac{1}{20}\\[0pt] -\frac{3}{20}&-\frac{1}{20}\end{pmatrix}.
```

The null space consists of vectors perpendicular to $`b`$. In two dimensions it is spanned by

```math
z=\begin{pmatrix}1\\[0pt] 1\end{pmatrix}
```

. All least-squares solutions are thus $`x=x_++tz`$ for real $`t`$. They give exactly the same fitted vector because $`Az=0`$.

The minimum-norm claim requires one further check. Since $`x_+`$ is parallel to $`b`$ while $`z`$ is perpendicular to $`b`$, their inner product is zero. Hence

```math
\Vert x_++tz\Vert _2^2=\Vert x_+\Vert _2^2+t^2\Vert z\Vert _2^2.
```

The minimum occurs uniquely at $`t=0`$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-91: rank-one pseudoinverse, case 5

Let

```math
a=\begin{pmatrix}1\\[0pt] -1\end{pmatrix}
```

,

```math
b=\begin{pmatrix}3\\[0pt] 2\end{pmatrix}
```

, $`A=ab^\top`$, and

```math
y=\begin{pmatrix}4\\[0pt] 1\end{pmatrix}
```

. Find every least-squares solution of $`Ax\approx y`$ and identify the one with smallest Euclidean norm.

Solution. Since $`Ax=a(b^\top x)`$, the scalar $`b^\top x`$ determines the entire output. The column space of

```math
A=\begin{pmatrix}3&2\\[0pt] -3&-2\end{pmatrix}
```

 is the line spanned by $`a`$. Projecting $`y`$ onto this line gives

```math
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}\frac{3}{2}\\[0pt] -\frac{3}{2}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{5}{2}\\[0pt] \frac{5}{2}\end{pmatrix}.
```

Here $`a^\top a=2`$, $`b^\top b=13`$, and $`a^\top y=3`$. The check $`a^\top r=0`$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $`b^\top x=\frac{3}{2}`$.

The vector parallel to $`b`$ satisfying that equation is

```math
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{9}{26}\\[0pt] \frac{3}{13}\end{pmatrix}.
```

Consequently the pseudoinverse is

```math
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{3}{26}&-\frac{3}{26}\\[0pt] \frac{1}{13}&-\frac{1}{13}\end{pmatrix}.
```

The null space consists of vectors perpendicular to $`b`$. In two dimensions it is spanned by

```math
z=\begin{pmatrix}-2\\[0pt] 3\end{pmatrix}
```

. All least-squares solutions are thus $`x=x_++tz`$ for real $`t`$. They give exactly the same fitted vector because $`Az=0`$.

The minimum-norm claim requires one further check. Since $`x_+`$ is parallel to $`b`$ while $`z`$ is perpendicular to $`b`$, their inner product is zero. Hence

```math
\Vert x_++tz\Vert _2^2=\Vert x_+\Vert _2^2+t^2\Vert z\Vert _2^2.
```

The minimum occurs uniquely at $`t=0`$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-92: rank-one pseudoinverse, case 6

Let

```math
a=\begin{pmatrix}2\\[0pt] 3\end{pmatrix}
```

,

```math
b=\begin{pmatrix}2\\[0pt] -3\end{pmatrix}
```

, $`A=ab^\top`$, and

```math
y=\begin{pmatrix}1\\[0pt] -2\end{pmatrix}
```

. Find every least-squares solution of $`Ax\approx y`$ and identify the one with smallest Euclidean norm.

Solution. Since $`Ax=a(b^\top x)`$, the scalar $`b^\top x`$ determines the entire output. The column space of

```math
A=\begin{pmatrix}4&-6\\[0pt] 6&-9\end{pmatrix}
```

 is the line spanned by $`a`$. Projecting $`y`$ onto this line gives

```math
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}-\frac{8}{13}\\[0pt] -\frac{12}{13}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{21}{13}\\[0pt] -\frac{14}{13}\end{pmatrix}.
```

Here $`a^\top a=13`$, $`b^\top b=13`$, and $`a^\top y=-4`$. The check $`a^\top r=0`$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $`b^\top x=-\frac{4}{13}`$.

The vector parallel to $`b`$ satisfying that equation is

```math
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}-\frac{8}{169}\\[0pt] \frac{12}{169}\end{pmatrix}.
```

Consequently the pseudoinverse is

```math
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{4}{169}&\frac{6}{169}\\[0pt] -\frac{6}{169}&-\frac{9}{169}\end{pmatrix}.
```

The null space consists of vectors perpendicular to $`b`$. In two dimensions it is spanned by

```math
z=\begin{pmatrix}3\\[0pt] 2\end{pmatrix}
```

. All least-squares solutions are thus $`x=x_++tz`$ for real $`t`$. They give exactly the same fitted vector because $`Az=0`$.

The minimum-norm claim requires one further check. Since $`x_+`$ is parallel to $`b`$ while $`z`$ is perpendicular to $`b`$, their inner product is zero. Hence

```math
\Vert x_++tz\Vert _2^2=\Vert x_+\Vert _2^2+t^2\Vert z\Vert _2^2.
```

The minimum occurs uniquely at $`t=0`$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

### 7.5 Curvature along a constraint

A quadratic form need not be positive definite on the whole plane to have a unique minimum on a line. Under $`x_1+x_2=s`$, the only permitted changes are multiples of $`(1,-1)`$. The relevant curvature is the quadratic form evaluated along that direction.

#### Problem LA-93: constrained diagonal quadratic, case 1

Minimize $`f(x)=1x_1^2+(4)x_2^2`$ subject to $`x_1+x_2=1`$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $`x_2=1-x_1`$. The objective on the feasible line is

```math
g(x_1)=1x_1^2+(4)(1-x_1)^2.
```

Its derivative is $`g'(x_1)=2(5)x_1-2(4)(1)`$ and its second derivative is $`2(5)\gt 0`$. Therefore

```math
x_1^*=\frac{h_2s}{h_1+h_2}=\frac{4}{5},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=\frac{1}{5}.
```

These coordinates sum to $`1`$, as required. The objective value is $`\frac{4}{5}`$.

For a separate derivation, use the convention $`L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)`$. Stationarity gives $`2h_1x_1=\lambda`$ and $`2h_2x_2=\lambda`$. Together with the equality constraint, these equations yield the same coordinates and $`\lambda=\frac{8}{5}`$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $`x^*+(t,-t)`$. Expanding the objective and using stationarity cancels the linear terms:

```math
f(x^*+(t,-t))-f(x^*)=(1+4)t^2=5t^2.
```

This difference is positive whenever $`t\ne0`$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-94: constrained diagonal quadratic, case 2

Minimize $`f(x)=2x_1^2+(3)x_2^2`$ subject to $`x_1+x_2=2`$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $`x_2=2-x_1`$. The objective on the feasible line is

```math
g(x_1)=2x_1^2+(3)(2-x_1)^2.
```

Its derivative is $`g'(x_1)=2(5)x_1-2(3)(2)`$ and its second derivative is $`2(5)\gt 0`$. Therefore

```math
x_1^*=\frac{h_2s}{h_1+h_2}=\frac{6}{5},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=\frac{4}{5}.
```

These coordinates sum to $`2`$, as required. The objective value is $`\frac{24}{5}`$.

For a separate derivation, use the convention $`L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)`$. Stationarity gives $`2h_1x_1=\lambda`$ and $`2h_2x_2=\lambda`$. Together with the equality constraint, these equations yield the same coordinates and $`\lambda=\frac{24}{5}`$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $`x^*+(t,-t)`$. Expanding the objective and using stationarity cancels the linear terms:

```math
f(x^*+(t,-t))-f(x^*)=(2+3)t^2=5t^2.
```

This difference is positive whenever $`t\ne0`$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-95: constrained diagonal quadratic, case 3

Minimize $`f(x)=0x_1^2+(5)x_2^2`$ subject to $`x_1+x_2=1`$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $`x_2=1-x_1`$. The objective on the feasible line is

```math
g(x_1)=0x_1^2+(5)(1-x_1)^2.
```

Its derivative is $`g'(x_1)=2(5)x_1-2(5)(1)`$ and its second derivative is $`2(5)\gt 0`$. Therefore

```math
x_1^*=\frac{h_2s}{h_1+h_2}=1,\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=0.
```

These coordinates sum to $`1`$, as required. The objective value is $`0`$.

For a separate derivation, use the convention $`L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)`$. Stationarity gives $`2h_1x_1=\lambda`$ and $`2h_2x_2=\lambda`$. Together with the equality constraint, these equations yield the same coordinates and $`\lambda=0`$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $`x^*+(t,-t)`$. Expanding the objective and using stationarity cancels the linear terms:

```math
f(x^*+(t,-t))-f(x^*)=(0+5)t^2=5t^2.
```

This difference is positive whenever $`t\ne0`$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-96: constrained diagonal quadratic, case 4

Minimize $`f(x)=-1x_1^2+(3)x_2^2`$ subject to $`x_1+x_2=1`$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $`x_2=1-x_1`$. The objective on the feasible line is

```math
g(x_1)=-1x_1^2+(3)(1-x_1)^2.
```

Its derivative is $`g'(x_1)=2(2)x_1-2(3)(1)`$ and its second derivative is $`2(2)\gt 0`$. Therefore

```math
x_1^*=\frac{h_2s}{h_1+h_2}=\frac{3}{2},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=-\frac{1}{2}.
```

These coordinates sum to $`1`$, as required. The objective value is $`-\frac{3}{2}`$.

For a separate derivation, use the convention $`L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)`$. Stationarity gives $`2h_1x_1=\lambda`$ and $`2h_2x_2=\lambda`$. Together with the equality constraint, these equations yield the same coordinates and $`\lambda=-3`$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $`x^*+(t,-t)`$. Expanding the objective and using stationarity cancels the linear terms:

```math
f(x^*+(t,-t))-f(x^*)=(-1+3)t^2=2t^2.
```

This difference is positive whenever $`t\ne0`$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-97: constrained diagonal quadratic, case 5

Minimize $`f(x)=4x_1^2+(1)x_2^2`$ subject to $`x_1+x_2=-2`$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $`x_2=-2-x_1`$. The objective on the feasible line is

```math
g(x_1)=4x_1^2+(1)(-2-x_1)^2.
```

Its derivative is $`g'(x_1)=2(5)x_1-2(1)(-2)`$ and its second derivative is $`2(5)\gt 0`$. Therefore

```math
x_1^*=\frac{h_2s}{h_1+h_2}=-\frac{2}{5},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=-\frac{8}{5}.
```

These coordinates sum to $`-2`$, as required. The objective value is $`\frac{16}{5}`$.

For a separate derivation, use the convention $`L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)`$. Stationarity gives $`2h_1x_1=\lambda`$ and $`2h_2x_2=\lambda`$. Together with the equality constraint, these equations yield the same coordinates and $`\lambda=-\frac{16}{5}`$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $`x^*+(t,-t)`$. Expanding the objective and using stationarity cancels the linear terms:

```math
f(x^*+(t,-t))-f(x^*)=(4+1)t^2=5t^2.
```

This difference is positive whenever $`t\ne0`$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-98: constrained diagonal quadratic, case 6

Minimize $`f(x)=3x_1^2+(3)x_2^2`$ subject to $`x_1+x_2=1`$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $`x_2=1-x_1`$. The objective on the feasible line is

```math
g(x_1)=3x_1^2+(3)(1-x_1)^2.
```

Its derivative is $`g'(x_1)=2(6)x_1-2(3)(1)`$ and its second derivative is $`2(6)\gt 0`$. Therefore

```math
x_1^*=\frac{h_2s}{h_1+h_2}=\frac{1}{2},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=\frac{1}{2}.
```

These coordinates sum to $`1`$, as required. The objective value is $`\frac{3}{2}`$.

For a separate derivation, use the convention $`L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)`$. Stationarity gives $`2h_1x_1=\lambda`$ and $`2h_2x_2=\lambda`$. Together with the equality constraint, these equations yield the same coordinates and $`\lambda=3`$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $`x^*+(t,-t)`$. Expanding the objective and using stationarity cancels the linear terms:

```math
f(x^*+(t,-t))-f(x^*)=(3+3)t^2=6t^2.
```

This difference is positive whenever $`t\ne0`$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

### 7.6 A small residual in a weak direction

A residual measures how well a computed vector satisfies the supplied equations. A forward error measures its distance from the true solution. If the matrix almost removes a direction, a large forward error in that direction can produce a small residual. These cases keep the distinction visible using exact arithmetic.

#### Problem LA-99: weak-direction sensitivity, case 1

Let $`A=\mathop{\mathrm{diag}} (1,\frac{1}{10})`$ and $`b=(1,\frac{1}{10})^\top`$. Replace the second component of $`b`$ by $`\frac{11}{100}`$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $`x_1=1`$ and $`\frac{1}{10}x_2=\frac{1}{10}`$, so $`x=(1,1)^\top`$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

```math
\tilde x=\begin{pmatrix}1\\[0pt] \frac{11}{10}\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\[0pt] \frac{1}{10}\end{pmatrix}.
```

The input perturbation has norm $`\vert \delta\vert =\frac{1}{100}`$, while the solution perturbation has norm $`\vert \delta\vert /\varepsilon=\frac{1}{10}`$. Division by $`\varepsilon`$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $`r=b-A\tilde x`$,

```math
r=\begin{pmatrix}0\\[0pt] -\frac{1}{100}\end{pmatrix},\qquad
\Vert r\Vert _2=\frac{1}{100}.
```

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $`A`$ are $`1`$ and $`\varepsilon=\frac{1}{10}`$. Therefore

```math
\Vert A\Vert _2=1,\qquad \Vert A^{-1}\Vert _2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=10.
```

The general absolute bound $`\Vert \Delta x\Vert _2\le\Vert A^{-1}\Vert _2\Vert \Delta b\Vert _2`$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-100: weak-direction sensitivity, case 2

Let $`A=\mathop{\mathrm{diag}} (1,\frac{1}{100})`$ and $`b=(1,\frac{1}{100})^\top`$. Replace the second component of $`b`$ by $`\frac{11}{1000}`$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $`x_1=1`$ and $`\frac{1}{100}x_2=\frac{1}{100}`$, so $`x=(1,1)^\top`$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

```math
\tilde x=\begin{pmatrix}1\\[0pt] \frac{11}{10}\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\[0pt] \frac{1}{10}\end{pmatrix}.
```

The input perturbation has norm $`\vert \delta\vert =\frac{1}{1000}`$, while the solution perturbation has norm $`\vert \delta\vert /\varepsilon=\frac{1}{10}`$. Division by $`\varepsilon`$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $`r=b-A\tilde x`$,

```math
r=\begin{pmatrix}0\\[0pt] -\frac{1}{1000}\end{pmatrix},\qquad
\Vert r\Vert _2=\frac{1}{1000}.
```

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $`A`$ are $`1`$ and $`\varepsilon=\frac{1}{100}`$. Therefore

```math
\Vert A\Vert _2=1,\qquad \Vert A^{-1}\Vert _2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=100.
```

The general absolute bound $`\Vert \Delta x\Vert _2\le\Vert A^{-1}\Vert _2\Vert \Delta b\Vert _2`$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-101: weak-direction sensitivity, case 3

Let $`A=\mathop{\mathrm{diag}} (1,\frac{1}{1000})`$ and $`b=(1,\frac{1}{1000})^\top`$. Replace the second component of $`b`$ by $`\frac{1}{500}`$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $`x_1=1`$ and $`\frac{1}{1000}x_2=\frac{1}{1000}`$, so $`x=(1,1)^\top`$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

```math
\tilde x=\begin{pmatrix}1\\[0pt] 2\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\[0pt] 1\end{pmatrix}.
```

The input perturbation has norm $`\vert \delta\vert =\frac{1}{1000}`$, while the solution perturbation has norm $`\vert \delta\vert /\varepsilon=1`$. Division by $`\varepsilon`$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $`r=b-A\tilde x`$,

```math
r=\begin{pmatrix}0\\[0pt] -\frac{1}{1000}\end{pmatrix},\qquad
\Vert r\Vert _2=\frac{1}{1000}.
```

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $`A`$ are $`1`$ and $`\varepsilon=\frac{1}{1000}`$. Therefore

```math
\Vert A\Vert _2=1,\qquad \Vert A^{-1}\Vert _2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=1000.
```

The general absolute bound $`\Vert \Delta x\Vert _2\le\Vert A^{-1}\Vert _2\Vert \Delta b\Vert _2`$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-102: weak-direction sensitivity, case 4

Let $`A=\mathop{\mathrm{diag}} (1,\frac{1}{10000})`$ and $`b=(1,\frac{1}{10000})^\top`$. Replace the second component of $`b`$ by $`\frac{101}{1000000}`$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $`x_1=1`$ and $`\frac{1}{10000}x_2=\frac{1}{10000}`$, so $`x=(1,1)^\top`$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

```math
\tilde x=\begin{pmatrix}1\\[0pt] \frac{101}{100}\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\[0pt] \frac{1}{100}\end{pmatrix}.
```

The input perturbation has norm $`\vert \delta\vert =\frac{1}{1000000}`$, while the solution perturbation has norm $`\vert \delta\vert /\varepsilon=\frac{1}{100}`$. Division by $`\varepsilon`$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $`r=b-A\tilde x`$,

```math
r=\begin{pmatrix}0\\[0pt] -\frac{1}{1000000}\end{pmatrix},\qquad
\Vert r\Vert _2=\frac{1}{1000000}.
```

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $`A`$ are $`1`$ and $`\varepsilon=\frac{1}{10000}`$. Therefore

```math
\Vert A\Vert _2=1,\qquad \Vert A^{-1}\Vert _2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=10000.
```

The general absolute bound $`\Vert \Delta x\Vert _2\le\Vert A^{-1}\Vert _2\Vert \Delta b\Vert _2`$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-103: weak-direction sensitivity, case 5

Let $`A=\mathop{\mathrm{diag}} (1,\frac{1}{1000000})`$ and $`b=(1,\frac{1}{1000000})^\top`$. Replace the second component of $`b`$ by $`\frac{1}{500000}`$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $`x_1=1`$ and $`\frac{1}{1000000}x_2=\frac{1}{1000000}`$, so $`x=(1,1)^\top`$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

```math
\tilde x=\begin{pmatrix}1\\[0pt] 2\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\[0pt] 1\end{pmatrix}.
```

The input perturbation has norm $`\vert \delta\vert =\frac{1}{1000000}`$, while the solution perturbation has norm $`\vert \delta\vert /\varepsilon=1`$. Division by $`\varepsilon`$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $`r=b-A\tilde x`$,

```math
r=\begin{pmatrix}0\\[0pt] -\frac{1}{1000000}\end{pmatrix},\qquad
\Vert r\Vert _2=\frac{1}{1000000}.
```

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $`A`$ are $`1`$ and $`\varepsilon=\frac{1}{1000000}`$. Therefore

```math
\Vert A\Vert _2=1,\qquad \Vert A^{-1}\Vert _2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=1000000.
```

The general absolute bound $`\Vert \Delta x\Vert _2\le\Vert A^{-1}\Vert _2\Vert \Delta b\Vert _2`$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-104: weak-direction sensitivity, case 6

Let $`A=\mathop{\mathrm{diag}} (1,\frac{1}{1000000})`$ and $`b=(1,\frac{1}{1000000})^\top`$. Replace the second component of $`b`$ by $`-\frac{9}{1000000}`$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $`x_1=1`$ and $`\frac{1}{1000000}x_2=\frac{1}{1000000}`$, so $`x=(1,1)^\top`$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

```math
\tilde x=\begin{pmatrix}1\\[0pt] -9\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\[0pt] -10\end{pmatrix}.
```

The input perturbation has norm $`\vert \delta\vert =\frac{1}{100000}`$, while the solution perturbation has norm $`\vert \delta\vert /\varepsilon=10`$. Division by $`\varepsilon`$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $`r=b-A\tilde x`$,

```math
r=\begin{pmatrix}0\\[0pt] \frac{1}{100000}\end{pmatrix},\qquad
\Vert r\Vert _2=\frac{1}{100000}.
```

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $`A`$ are $`1`$ and $`\varepsilon=\frac{1}{1000000}`$. Therefore

```math
\Vert A\Vert _2=1,\qquad \Vert A^{-1}\Vert _2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=1000000.
```

The general absolute bound $`\Vert \Delta x\Vert _2\le\Vert A^{-1}\Vert _2\Vert \Delta b\Vert _2`$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

[Contents](README.md) · [Part 1](01-linear-algebra-part-1.md)
