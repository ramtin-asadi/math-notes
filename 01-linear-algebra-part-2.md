# Linear Algebra — Part 2

These notes use real finite-dimensional spaces unless a different field is stated. Vectors are columns, $`A^\top`$ is the transpose, $`I`$ is the identity, and $`\Vert x\Vert _2=(x^\top x)^{1/2}`$. For complex vectors, replace transposes by conjugate transposes and use the Hermitian inner product.

Sources:

- Gilbert Strang, [MIT 18.06SC, Fall 2011](https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/), especially Units II and III.
- Sheldon Axler, [Linear Algebra Done Right, fourth edition](https://linear.axler.net/LADR4e.pdf), Chapters 2, 3, 6, and 7. Exercise numbers below refer to this edition.

The explanations and solutions are written for these notes. Source problems are briefly paraphrased and linked; examples and additional problems are constructed here.

[Contents](README.md) · [Part 1](01-linear-algebra-part-1.md) · [Part 3](01-linear-algebra-part-3.md)

Contents:

- [4. SVD and low-rank approximation](#4-svd-and-low-rank-approximation)
- [5. Conditioning and regularization](#5-conditioning-and-regularization)
- [6. Matrix differentiation](#6-matrix-differentiation)

## 4. SVD and low-rank approximation

### Further problems on eigenvalues and quadratic forms

#### Problem LA-37: an all-ones matrix

Source: [Axler, Exercises 7C, 3](https://linear.axler.net/LADR4e.pdf).

Problem. Prove positivity of the $`n\times n`$ all-ones matrix and determine its rank, null space, and eigenvalues.

Solution. For $`J=\mathbf1\mathbf1^\top`$,

```math
x^\top Jx=(\mathbf1^\top x)^2\ge0.
```

It is symmetric and PSD. The vector $`\mathbf1`$ has eigenvalue $`n`$, since $`J\mathbf1=n\mathbf1`$. Every vector whose coordinates sum to zero has eigenvalue zero. Thus the rank is one, except in the vacuous zero-dimensional setting, and the nullity is $`n-1`$. Positive entries throughout the matrix have not prevented a large null space.

#### Problem LA-38: adding a positive semidefinite matrix

Source: [Axler, Exercises 7C, 7](https://linear.axler.net/LADR4e.pdf).

Problem. Show that adding a positive-semidefinite operator to a positive-definite operator preserves positive definiteness and invertibility.

Solution. Let $`S`$ be PD and $`T`$ PSD. For nonzero $`x`$, $`x^\top(S+T)x=x^\top Sx+x^\top Tx\gt 0+0`$. Hence $`S+T`$ is PD and invertible. A quantitative form is $`\lambda_{\min}(S+T)\ge\lambda_{\min}(S)`$, obtained by minimizing the Rayleigh quotient over unit vectors. Commutativity is unnecessary for this lower bound.

#### Problem LA-39: changing variables in a positive form

Source: [Axler, Exercises 7C, 9](https://linear.axler.net/LADR4e.pdf).

Problem. For PSD $`T`$ and any compatible linear map $`S`$, prove positivity of $`S^\top TS`$. When $`T`$ is PD, determine the condition for the new form to be PD.

Solution. If $`T`$ is PSD, then $`(S^\top TS)^\top=S^\top TS`$ and

```math
x^\top S^\top TSx=(Sx)^\top T(Sx)\ge0.
```

Thus $`S^\top TS`$ is PSD. If $`T`$ is PD, this new matrix is PD precisely when $`S`$ is injective: any nonzero kernel vector of $`S`$ makes the form zero, while injectivity ensures $`Sx\ne0`$ for nonzero $`x`$. This is the algebra behind both weighted least squares and covariance transformations.

#### Problem LA-40: products of positive operators

Source: [Axler, Exercises 7C, 18](https://linear.axler.net/LADR4e.pdf).

Problem. For self-adjoint PSD operators $`S,T`$ on the same space, prove that $`ST`$ is positive exactly when $`S`$ and $`T`$ commute.

Solution. If $`ST`$ is positive, it is self-adjoint, so $`ST=(ST)^\top=TS`$. Conversely suppose PSD matrices $`S,T`$ commute. Each eigenspace of $`S`$ is invariant under $`T`$, because $`S(Tv)=T(Sv)=\lambda Tv`$. Diagonalize the symmetric restriction of $`T`$ within each such eigenspace. This produces a common orthonormal eigenbasis in which both matrices are diagonal with nonnegative entries. Their product has products of those entries on the diagonal, all nonnegative. Thus $`ST`$ is PSD. Separate positivity does not suffice unless the product is also symmetric.

#### Problem LA-41: a three-dimensional diagonalization

Constructed problem. Diagonalize

```math
A=\begin{pmatrix}2&1&0\\[0pt] 1&2&0\\[0pt] 0&0&4\end{pmatrix}
```

and minimize its quadratic form on the unit sphere.

Solution. Choose $`q_1=(1,1,0)/\sqrt2`$, $`q_2=(1,-1,0)/\sqrt2`$, $`q_3=(0,0,1)`$. Their eigenvalues are $`3,1,4`$. For $`x=a q_1+b q_2+c q_3`$, the unit constraint is $`a^2+b^2+c^2=1`$ and the form is $`3a^2+b^2+4c^2`$. Its minimum is one, achieved only when $`a=c=0`$ and $`b=\pm1`$. The maximum is four at $`\pm q_3`$.

#### Problem LA-42: a rank-one covariance update

Constructed problem. Show that $`S+vv^\top`$ is PSD if $`S`$ is PSD, and determine when a null direction disappears.

Solution. For every $`x`$,

```math
x^\top(S+vv^\top)x=x^\top Sx+(v^\top x)^2\ge0.
```

Since both terms are nonnegative, this expression is zero exactly when $`x\in\ker S`$ and $`v^\top x=0`$. Hence $`\ker(S+vv^\top)=\ker S\cap v^\perp`$. A one-dimensional null space generated by $`z`$ disappears if $`v^\top z\ne0`$. A null space of dimension at least two cannot be removed completely by one rank-one update, because at least one independent null direction remains perpendicular to $`v`$.

#### Problem LA-43: a two-by-two positivity criterion

Constructed problem. Characterize PSD matrices

```math
A=\left(\begin{smallmatrix}a&c\\[0pt] c&b\end{smallmatrix}\right)
```

.

Solution. Necessity gives $`a,b\ge0`$ by testing the coordinate vectors. If $`a\gt 0`$, complete the square:

```math
ax^2+2cxy+by^2=a(x+cy/a)^2+(b-c^2/a)y^2.
```

Nonnegativity for every $`(x,y)`$ requires and is implied by $`ab-c^2\ge0`$. If $`a=0`$, nonnegativity forces $`c=0`$: otherwise fix $`y=1`$ and send $`x`$ in the direction making $`2cx+b`$ negative. Thus the full PSD criterion is $`a\ge0`$, $`b\ge0`$, and $`ab-c^2\ge0`$. PD requires $`a\gt 0`$ and $`ab-c^2\gt 0`$.

#### Problem LA-44: a diagonalizable matrix without an orthogonal eigenbasis

Constructed problem. Study

```math
A=\left(\begin{smallmatrix}1&1\\[0pt] 0&2\end{smallmatrix}\right)
```

.

Solution. The eigenvector directions are $`v_1=(1,0)`$ for eigenvalue one and $`v_2=(1,1)`$ for eigenvalue two. They are independent, so $`A=V\mathop{\mathrm{diag}} (1,2)V^{-1}`$. Their inner product is one, and scaling either vector cannot make it zero. Thus this diagonalizable matrix has no orthogonal eigenbasis. Indeed

```math
A^\top A=\left(\begin{smallmatrix}1&1\\[0pt] 1&5\end{smallmatrix}\right)
```

 differs from

```math
AA^\top=\left(\begin{smallmatrix}2&2\\[0pt] 2&4\end{smallmatrix}\right)
```

, so it is not normal.

### 4.1 Theory

Every real $`m\times n`$ matrix admits a singular value decomposition

```math
A=U\Sigma V^\top.
```

Here $`U`$ and $`V`$ are orthogonal, and the rectangular diagonal matrix $`\Sigma`$ has entries $`\sigma_1\ge\cdots\ge\sigma_r\gt 0`$ followed by zeros. The rank is $`r`$. Singular values are the square roots of the eigenvalues of $`A^\top A`$.

The SVD separates a map into input directions, stretches, and output directions:

```math
Av_i=\sigma_i u_i,\qquad
A=\sum_{i=1}^r\sigma_i u_iv_i^\top.
```

The Moore-Penrose pseudoinverse reverses the nonzero stretches:

```math
A^+=\sum_{i=1}^r\sigma_i^{-1}v_iu_i^\top.
```

Then $`AA^+`$ projects onto $`\mathop{\mathrm{range}} A`$, and $`A^+A`$ projects onto $`(\ker A)^\perp`$. The vector $`A^+b`$ has the smallest Euclidean norm among all least-squares solutions.

### 4.2 Derivations and proofs

#### Singular vectors belong to two different spaces

For a rectangular matrix, right singular vectors $`v_i`$ belong to the input space $`\mathbb R^n`$, while left singular vectors $`u_i`$ belong to the output space $`\mathbb R^m`$. The map takes $`v_i`$ to the scaled output direction $`\sigma_i u_i`$. These vectors should not be identified merely because both lists are orthonormal.

The matrix $`A^\top A`$ acts on the input space and identifies which input directions are stretched strongly. The matrix $`AA^\top`$ acts on the output space and identifies the matching output directions. Both have the same positive eigenvalues $`\sigma_i^2`$, but they may have different numbers of zero eigenvalues when $`m\ne n`$.

For $`b`$, first compute its output coordinates $`u_i^\top b`$. Divide each attainable coordinate by its singular value, then place it along $`v_i`$. Coordinates in the orthogonal complement of the output range cannot be produced and are discarded. This sequence is exactly the formula for $`A^+b`$.

#### All four Moore-Penrose identities

The SVD immediately gives

```math
AA^+A=A,\qquad A^+AA^+=A^+,\qquad
(AA^+)^\top=AA^+,\qquad(A^+A)^\top=A^+A.
```

For the first two identities, a positive singular value is multiplied by its reciprocal and then itself again. A zero singular value remains zero. For the symmetry identities, $`AA^+=\sum_{i=1}^r u_iu_i^\top`$ and $`A^+A=\sum_{i=1}^r v_iv_i^\top`$ are sums of symmetric rank-one projectors.

To see the smallest-norm property explicitly, let $`x_0=A^+b`$. Every other least-squares solution is $`x_0+z`$ for $`z\in\ker A`$. The SVD puts $`x_0`$ in the span of the positive right singular directions, which is perpendicular to $`\ker A`$. Therefore $`\Vert x_0+z\Vert ^2=\Vert x_0\Vert ^2+\Vert z\Vert ^2`$. The unique minimum norm occurs at $`z=0`$.

#### Reconstruction error in PCA

For a centered observation vector $`x`$, reconstruction using a unit direction $`w`$ is $`ww^\top x`$. Orthogonal decomposition gives

```math
\Vert x-ww^\top x\Vert ^2=\Vert x\Vert ^2-(w^\top x)^2.
```

Average this identity over the sample. The first term is fixed, so minimizing reconstruction error is equivalent to maximizing projected variance. For $`k`$ orthonormal columns of $`W`$, the same calculation gives residual sum of squares $`\Vert X\Vert _F^2-\mathop{\mathrm{tr}} (W^\top X^\top XW)`$. The spectral argument selects the leading $`k`$ eigenvectors.

Centering is part of this statement. Without centering, the leading singular direction can describe the mean offset instead of variation around the mean. If a leading eigenvalue is repeated, the best principal subspace may be unique at a larger dimension even though an individual first principal direction is not.

#### Constructing the SVD

Diagonalize $`A^\top A`$ using the spectral theorem. It is PSD because $`x^\top A^\top Ax=\Vert Ax\Vert ^2`$. For each positive eigenvalue $`\sigma_i^2`$, define $`u_i=Av_i/\sigma_i`$. Then

```math
u_i^\top u_j=\frac{v_i^\top A^\top Av_j}{\sigma_i\sigma_j}=\delta_{ij}.
```

For a zero eigenvalue, $`\Vert Av_i\Vert ^2=0`$ implies $`Av_i=0`$. Extend the positive-direction $`u_i`$ to an orthonormal basis of the output space. Applying both sides to the input basis proves the decomposition.

#### Best low-rank approximation

The truncated matrix $`A_k=\sum_{i=1}^k\sigma_i u_iv_i^\top`$ has rank at most $`k`$ and error $`\Vert A-A_k\Vert _2=\sigma_{k+1}`$. For any matrix $`B`$ of rank at most $`k`$, its restriction to $`\mathop{\mathrm{span}} (v_1,\ldots,v_{k+1})`$ has a nonzero kernel vector. Normalize one to $`z`$. Then $`Bz=0`$ and

```math
\Vert A-B\Vert _2\ge\Vert (A-B)z\Vert =\Vert Az\Vert \ge\sigma_{k+1}.
```

This proves spectral-norm optimality. In Frobenius norm, orthogonal projection of $`A`$ onto a fixed output subspace $`S`$ is $`P_SA`$, with residual squared $`\Vert A\Vert _F^2-\mathop{\mathrm{tr}} (P_SAA^\top)`$. Among spaces of dimension at most $`k`$, the trace is largest for the first $`k`$ left singular vectors: in an eigenbasis it is $`\sum_i\sigma_i^2p_i`$, where $`0\le p_i\le1`$ and $`\sum_i p_i\le k`$. Allocating those weights to the largest eigenvalues gives

```math
\min_{\mathop{\mathrm{rank}} B\le k}\Vert A-B\Vert _F^2=\sum_{i\gt k}\sigma_i^2.
```

#### PCA as maximum variance

For centered data with covariance $`S`$, the variance along a unit direction $`w`$ is $`w^\top Sw`$. Maximizing it subject to $`w^\top w=1`$ gives $`Sw=\lambda w`$ by a Lagrange multiplier. The Rayleigh bound identifies the largest eigenvalue as the maximum. Subsequent directions maximize variance under orthogonality to earlier ones. For centered $`X=U\Sigma V^\top`$, the principal directions are the columns of $`V`$ and the variance explained by direction $`i`$ is $`\sigma_i^2/(n-1)`$.

#### Centering, covariance rank, and reconstruction

Let $`X`$ be an $`n\times p`$ data matrix with observations in rows. The centering matrix is $`H=I-\mathbf1\mathbf1^\top/n`$. It is symmetric and idempotent, kills the constant vector, and has rank $`n-1`$. The centered matrix $`Z=HX`$ has zero column means.

The sample covariance is $`S=Z^\top Z/(n-1)`$ for $`n\gt 1`$. Since a Gram matrix has the same rank as its underlying matrix,

```math
\mathop{\mathrm{rank}} (S)=\mathop{\mathrm{rank}} (Z)
\le\min(p,n-1).
```

Thus more features than independent centered observations necessarily produce a singular covariance. Even with many observations, exact column relations can lower the rank further.

For a unit direction $`q`$, the rank-one reconstruction of all centered rows is $`Zqq^\top`$. Expanding its squared error gives

```math
\Vert Z-Zqq^\top\Vert _F^2
=\mathop{\mathrm{tr}} (Z^\top Z)-q^\top Z^\top Zq.
```

Minimizing reconstruction error is therefore equivalent to maximizing projected variance. With $`k`$ orthonormal directions collected in $`Q_k`$, the captured sum is $`\mathop{\mathrm{tr}} (Q_k^\top Z^\top ZQ_k)`$.

To justify choosing the top eigenvectors, let $`v_i`$ be an orthonormal eigenbasis with descending eigenvalues $`\lambda_i`$. Put $`w_i=\Vert Q_k^\top v_i\Vert ^2`$. These numbers satisfy $`0\le w_i\le1`$ and $`\sum_iw_i=k`$. The captured variance is $`\sum_i\lambda_iw_i`$, which cannot exceed $`\sum_{i=1}^k\lambda_i`$: moving any weight from a lower eigenvalue to an available higher one cannot decrease the sum. Equality is attained by the leading eigenvectors.

If the cutoff eigenvalue is tied, different subspaces within the tied eigenspace can be equally good. Individual eigenvector signs never change the reconstruction. PCA identifies directions or subspaces, rather than a unique signed list in every dataset.

#### Example LA-E21: covariance PCA and rescaled PCA

Use the four observations $`(2,0),(-2,0),(0,1),(0,-1)`$.

Solution. The mean is zero and the scatter matrix is $`Z^\top Z=\mathop{\mathrm{diag}} (8,2)`$. The sample covariance eigenvalues are $`8/3`$ and $`2/3`$. The first coordinate is the leading principal direction and captures $`8/(8+2)=80\%`$ of the variance.

Keeping only that coordinate reconstructs the first two observations exactly and maps the other two to zero. Their squared errors are one each, totaling two, equal to the discarded scatter eigenvalue.

If each coordinate is first divided by its sample standard deviation, the covariance becomes identity. Both standardized directions now have equal variance; there is no unique leading direction. Standardizing changed the metric in which reconstruction errors are measured. A PCA result therefore needs its centering and scaling convention stated along with the computed eigenvectors.

### 4.3 Worked examples

#### Example LA-E10: an explicit rectangular SVD

Let

```math
A=\left(\begin{smallmatrix}3&0\\[0pt] 0&2\\[0pt] 0&0\end{smallmatrix}\right)
```

.

Solution. The coordinate directions already give the SVD, with singular values $`3,2`$. The best rank-one approximation keeps only the entry $`3`$. Its spectral and Frobenius errors are both $`2`$. The pseudoinverse is

```math
\left(\begin{smallmatrix}1/3&0&0\\[0pt] 0&1/2&0\end{smallmatrix}\right)
```

, so $`A^+(6,6,7)^\top=(2,3)^\top`$ and the residual $`(0,0,7)^\top`$ cannot be fitted.

#### Example LA-E11: rank-deficient least squares

Minimize $`\Vert Ax-b\Vert ^2`$ for

```math
A=\left(\begin{smallmatrix}1&1\\[0pt] 2&2\end{smallmatrix}\right)
```

 and $`b=(1,0)^\top`$.

Solution. Only $`s=x_1+x_2`$ affects the fit. Minimize $`(s-1)^2+4s^2`$, whose derivative is $`10s-2`$. Hence $`s=1/5`$. All minimizers are $`(t,1/5-t)`$; the smallest norm occurs at equal coordinates, $`x^+=(1/10,1/10)`$. The fitted vector is $`(1/5,2/5)`$, and the residual $`(4/5,-2/5)`$ is perpendicular to $`(1,2)`$.

#### Example LA-E12: two principal directions

Let

```math
S=\left(\begin{smallmatrix}5&4\\[0pt] 4&5\end{smallmatrix}\right)
```

.

Solution. The unit eigenvectors are $`(1,1)/\sqrt2`$ and $`(1,-1)/\sqrt2`$, with eigenvalues $`9`$ and $`1`$. The first component explains $`9/(9+1)=90\%`$ of total variance. A one-dimensional reconstruction keeps the common movement of the coordinates and removes their difference. If the variables have different units, standardizing them changes the covariance matrix and can change these directions.

### 4.4 Selected exercises

#### Problem LA-11: pseudoinverse of a surjective map

Source: [Axler, Exercises 6C, 21(a,b,c)](https://linear.axler.net/LADR4e.pdf). Let $`T(a,b,c)=(a+b+c,2b+3c)`$.

Solution. Its matrix is

```math
A=\left(\begin{smallmatrix}1&1&1\\[0pt] 0&2&3\end{smallmatrix}\right)
```

. Independent rows make it surjective. Therefore

```math
A^+=A^\top(AA^\top)^{-1},\quad
AA^\top=\begin{pmatrix}3&5\\[0pt] 5&13\end{pmatrix},\quad
(AA^\top)^{-1}=\frac1{14}\begin{pmatrix}13&-5\\[0pt] -5&3\end{pmatrix}.
```

Multiplying yields

```math
A^+\binom xy=\frac1{14}
\begin{pmatrix}13x-5y\\[0pt] 3x+y\\[0pt] -2x+4y\end{pmatrix}.
```

Direct multiplication gives $`AA^+=I_2`$, the projection onto the full range. A kernel basis is $`n=(1,-3,2)^\top`$, with $`\Vert n\Vert ^2=14`$. Another multiplication gives

```math
A^+A=I_3-\frac{nn^\top}{14}.
```

This is the projection onto the kernel's orthogonal complement, establishing both requested identities and the minimum-norm property.

#### Problem LA-12: singular values without eigenvalue calculations

Source: [Axler, Exercises 7E, 5](https://linear.axler.net/LADR4e.pdf). For $`T(x,y)=(-4y,x)`$, find the singular values.

Solution. Its matrix satisfies $`A^\top A=\mathop{\mathrm{diag}} (1,16)`$, so the singular values in decreasing order are $`4,1`$. One SVD uses $`v_1=e_2,v_2=e_1`$ and $`u_1=-e_1,u_2=e_2`$. Then $`Av_1=4u_1`$ and $`Av_2=u_2`$. The eigenvalues of $`A`$ itself are imaginary, which does not affect the real nonnegative singular values.

#### Problem LA-13: singular values of a self-adjoint operator

Source: [Axler, Exercises 7E, 7](https://linear.axler.net/LADR4e.pdf).

Problem. Relate the singular values of a self-adjoint operator to its eigenvalues. Extend the argument to a complex normal operator.

Solution. For a real symmetric matrix, $`A=Q\Lambda Q^\top`$ gives $`A^\top A=Q\Lambda^2Q^\top`$. The singular values are therefore $`\vert \lambda_i\vert`$, sorted in decreasing order. For a complex normal operator, the spectral theorem gives $`A=Q\Lambda Q^*`$ and $`A^*A=Q\vert \Lambda\vert ^2Q^*`$, proving the same conclusion. Self-adjointness or normality is needed: a nonzero nilpotent matrix can have only zero eigenvalues and still have a positive singular value.

## 5. Conditioning and regularization

### Further problems on SVD and approximation

#### Problem LA-45: transposition preserves positive singular values

Source: [Axler, Exercises 7E, 9](https://linear.axler.net/LADR4e.pdf).

Problem. Prove that a matrix and its transpose have the same positive singular values, including the rectangular case.

Solution. Transposing $`A=\sum_{i=1}^r\sigma_i u_iv_i^\top`$ gives $`A^\top=\sum_{i=1}^r\sigma_i v_iu_i^\top`$. Both orthonormal lists simply exchange roles, leaving the positive singular values unchanged. If the matrix is rectangular, one of $`A^\top A`$ and $`AA^\top`$ has additional zeros, which is why the claim is stated for positive singular values.

#### Problem LA-46: singular values of an inverse

Source: [Axler, Exercises 7E, 10](https://linear.axler.net/LADR4e.pdf).

Problem. Express the singular values of an invertible matrix's inverse in terms of those of the original matrix, with the correct ordering.

Solution. If $`A`$ is square and invertible, all singular values are positive. Its inverse is $`A^{-1}=V\Sigma^{-1}U^\top`$. Thus the singular values are the reciprocals, reordered as $`1/\sigma_n\ge\cdots\ge1/\sigma_1`$. In particular, the inverse amplifies most strongly the output direction associated with the smallest original singular value.

#### Problem LA-47: a lower singular-value bound

Source: [Axler, Exercises 7E, 14](https://linear.axler.net/LADR4e.pdf).

Problem. Prove $`\Vert Ax\Vert _2\ge\sigma_{\min}\Vert x\Vert _2`$, including zero singular directions when the map is not injective.

Solution. Expand $`x=\sum_i c_iv_i`$ in the full right singular basis, including zero singular directions. Then

```math
\Vert Ax\Vert ^2=\sum_i\sigma_i^2c_i^2
\ge\sigma_{\min}^2\sum_i c_i^2
=\sigma_{\min}^2\Vert x\Vert ^2.
```

Taking square roots proves the bound. When $`\sigma_{\min}\gt 0`$, the map is injective and cannot collapse any nonzero direction. When $`\sigma_{\min}=0`$, the inequality still holds but provides no positive lower bound.

#### Problem LA-48: transposing a pseudoinverse

Source: [Axler, Exercises 7E, 16](https://linear.axler.net/LADR4e.pdf).

Problem. Prove $`(A^\top)^+=(A^+)^\top`$ for a matrix that may be rank deficient.

Solution. From the SVD, $`A^+=\sum_i\sigma_i^{-1}v_iu_i^\top`$. Its transpose is $`\sum_i\sigma_i^{-1}u_iv_i^\top`$. Applying the pseudoinverse construction to $`A^\top=\sum_i\sigma_i v_iu_i^\top`$ produces that same expression. Hence $`(A^\top)^+=(A^+)^\top`$. The argument includes rank-deficient matrices because the sums only use positive singular values.

#### Problem LA-49: a non-diagonal rank-one approximation

Constructed problem. Find the best rank-one approximation to

```math
A=\left(\begin{smallmatrix}2&1\\[0pt] 1&2\end{smallmatrix}\right)
```

.

Solution. Its symmetric eigendecomposition has eigenvalues $`3,1`$ with unit directions $`q_1=(1,1)/\sqrt2`$, $`q_2=(1,-1)/\sqrt2`$. Both eigenvalues are positive, so this is also an SVD. Keep the first term:

```math
A_1=3q_1q_1^\top
=\frac32\begin{pmatrix}1&1\\[0pt] 1&1\end{pmatrix}.
```

The residual is

```math
q_2q_2^\top=\tfrac12\left(\begin{smallmatrix}1&-1\\[0pt] -1&1\end{smallmatrix}\right)
```

, with spectral and Frobenius norm one. The best-approximation theorem proves that no rank-one matrix can reduce either error below one.

#### Problem LA-50: a rank-one pseudoinverse

Constructed problem. For nonzero $`a\in\mathbb R^m`$ and $`b\in\mathbb R^n`$, find $`(ab^\top)^+`$.

Solution. Normalize the vectors: $`u=a/\Vert a\Vert`$, $`v=b/\Vert b\Vert`$. Then $`ab^\top=(\Vert a\Vert \Vert b\Vert )uv^\top`$ is a one-term SVD. Reversing its only nonzero singular value gives

```math
(ab^\top)^+=\frac{ba^\top}{\Vert a\Vert ^2\Vert b\Vert ^2}.
```

Multiplication gives $`AA^+=aa^\top/\Vert a\Vert ^2`$ and $`A^+A=bb^\top/\Vert b\Vert ^2`$, the appropriate line projectors. If either original vector is zero, the matrix and its pseudoinverse are zero, and the displayed division formula must not be used.

#### Problem LA-51: PCA from actual observations

Constructed problem. The four centered observations are $`(2,1),(-2,-1),(1,2),(-1,-2)`$. Find both principal directions and explained-variance shares.

Solution. The column sums vanish. Summing outer products gives

```math
X^\top X=\begin{pmatrix}10&8\\[0pt] 8&10\end{pmatrix},\qquad
S=\frac13\begin{pmatrix}10&8\\[0pt] 8&10\end{pmatrix}.
```

The eigenvectors are $`(1,1)/\sqrt2`$ and $`(1,-1)/\sqrt2`$, with covariance eigenvalues $`6`$ and $`2/3`$. Total variance is $`20/3`$, so the shares are $`9/10`$ and $`1/10`$. The leading scores are $`3/\sqrt2,-3/\sqrt2,3/\sqrt2,-3/\sqrt2`$. Reconstructing each observation from its leading score produces $`(3/2,3/2)`$ or its negative. The residual sum of squares is the discarded eigenvalue of $`X^\top X`$, namely two.

#### Problem LA-52: nonuniqueness of a best rank-one approximation

Constructed problem. Find two distinct best rank-one approximations to $`I_2`$.

Solution. Both $`\mathop{\mathrm{diag}} (1,0)`$ and $`\mathop{\mathrm{diag}} (0,1)`$ have rank one and leave a residual with spectral norm one and Frobenius norm one. The two singular values of $`I_2`$ are tied at one, so either direction can be retained. More generally $`ww^\top`$ is best for every unit $`w`$. Its residual is the orthogonal projector onto the perpendicular line. A spectral gap between the last retained and first discarded singular values removes this particular ambiguity in the optimal subspace.

### 5.1 Theory

The induced matrix norm is $`\Vert A\Vert =\sup_{\Vert x\Vert =1}\Vert Ax\Vert`$. For the Euclidean norm it equals the largest singular value. For invertible $`A`$,

```math
\kappa_2(A)=\Vert A\Vert _2\Vert A^{-1}\Vert _2
=\frac{\sigma_{\max}}{\sigma_{\min}}\ge1.
```

Conditioning describes sensitivity of the mathematical problem to changes in its inputs. Stability describes how an algorithm handles rounding. A stable algorithm can still produce a sensitive answer when the problem is ill-conditioned.

Ridge regularization minimizes $`\Vert Ax-b\Vert ^2+\lambda\Vert x\Vert ^2`$ for $`\lambda\gt 0`$. It adds curvature in directions where the data give little or no information. It also changes the objective being solved; reducing sensitivity generally introduces shrinkage.

### 5.2 Derivations and proofs

#### Absolute and relative sensitivity

An absolute perturbation bound compares $`\Vert \delta x\Vert`$ with the size of the input change. A relative bound divides by the original sizes. These answer different questions when the true output is small. Relative error is undefined at a zero exact output, so an absolute error or a scale chosen independently of that output is needed there.

For fixed invertible $`A`$, $`x=A^{-1}b`$ and $`\delta x=A^{-1}\delta b`$ exactly. Thus the absolute condition number with respect to $`b`$ in the chosen norm is $`\Vert A^{-1}\Vert`$. At a particular nonzero $`b`$, the relative condition number is $`\Vert A^{-1}\Vert \Vert b\Vert /\Vert A^{-1}b\Vert`$, which is at most $`\kappa(A)`$. The usual matrix condition number is a bound over right-hand-side directions; it need not equal the sensitivity at every particular $`b`$.

#### Normal equations square the condition number

For a full-column-rank matrix, the positive eigenvalues of $`A^\top A`$ are $`\sigma_i^2`$. Hence

```math
\kappa_2(A^\top A)=\frac{\sigma_{\max}^2}{\sigma_{\min}^2}
=\kappa_2(A)^2.
```

The normal equations are exact algebra, but explicitly forming $`A^\top A`$ in finite precision can lose information in small singular directions. QR works with orthogonal transformations whose Euclidean condition number is one. It avoids this squaring at the level of the matrix being solved, although the underlying least-squares problem can still be sensitive.

#### Choosing a regularization scale

The ridge transition occurs around $`\sigma_i^2\approx\lambda`$. If $`\sigma_i^2\gg\lambda`$, then $`\sigma_i/(\sigma_i^2+\lambda)`$ is close to $`1/\sigma_i`$. If $`\sigma_i^2\ll\lambda`$, then the factor is approximately $`\sigma_i/\lambda`$, which strongly suppresses that direction.

This also means that units affect the penalty. Rescaling one column of $`A`$ changes the numerical size of its coefficient without changing the fitted values available to the unpenalized model. A penalty $`\lambda\Vert x\Vert ^2`$ is applied to coefficients in the chosen units, so rescaling columns can change the regularized fitted result. Standardization or a deliberately chosen matrix penalty $`x^\top Rx`$ makes that choice explicit.

#### Perturbing a linear system

Suppose $`Ax=b`$ and $`(A+E)(x+\delta x)=b+e`$. Subtracting gives

```math
(A+E)\delta x=e-Ex.
```

If $`\Vert A^{-1}E\Vert \lt 1`$, the series $`\sum_{j\ge0}(-A^{-1}E)^j`$ converges and inverts $`I+A^{-1}E`$. Thus

```math
\Vert \delta x\Vert \le
\frac{\Vert A^{-1}\Vert }{1-\Vert A^{-1}E\Vert }
(\Vert e\Vert +\Vert E\Vert \Vert x\Vert ).
```

For $`b\ne0`$, $`x\ne0`$, and $`\kappa(A)\Vert E\Vert /\Vert A\Vert \lt 1`$, this implies

```math
\frac{\Vert \delta x\Vert }{\Vert x\Vert }
\le
\frac{\kappa(A)}
{1-\kappa(A)\Vert E\Vert /\Vert A\Vert }
\left(\frac{\Vert e\Vert }{\Vert b\Vert }+\frac{\Vert E\Vert }{\Vert A\Vert }\right).
```

We used $`\Vert b\Vert \le\Vert A\Vert \Vert x\Vert`$ to express the right-hand-side perturbation relatively. This is an upper bound; particular perturbation directions may have much smaller effects.

#### Ridge in singular directions

Differentiating the ridge objective gives $`(A^\top A+\lambda I)x=A^\top b`$. Its Hessian is $`2(A^\top A+\lambda I)`$ and is PD, so the unique solution is

```math
x_\lambda=\sum_{i=1}^r
\frac{\sigma_i}{\sigma_i^2+\lambda}(u_i^\top b)v_i.
```

Unregularized inversion uses $`1/\sigma_i`$; ridge replaces it by the bounded factor $`\sigma_i/(\sigma_i^2+\lambda)`$. The Hessian eigenvalues are $`2(\sigma_i^2+\lambda)`$, including $`2\lambda`$ in null directions. Consequently its condition number is $`(\sigma_{\max}^2+\lambda)/(\sigma_{\min}^2+\lambda)`$, with $`\sigma_{\min}=0`$ if the column rank is deficient.

#### General quadratic regularization

Consider $`\min_x\Vert Ax-b\Vert ^2+\lambda\Vert Lx\Vert ^2`$ with $`\lambda\gt 0`$. The matrix $`L`$ can penalize selected coordinates or differences between neighboring coordinates. Its normal equation is

```math
(A^\top A+\lambda L^\top L)x=A^\top b.
```

For any direction $`h`$, the associated quadratic form is $`\Vert Ah\Vert ^2+\lambda\Vert Lh\Vert ^2`$. It vanishes exactly when $`h`$ belongs to both kernels. Consequently the objective has a unique minimizer precisely when

```math
\ker A\cap\ker L=\lbrace 0\rbrace .
```

Penalizing some directions does not cure an unidentified direction left unpenalized by $`L`$. The standard ridge choice $`L=I`$ always gives a positive-definite system, but a penalty excluding an intercept requires checking the joint-kernel condition.

When $`L=I`$ and $`A=U\Sigma V^\top`$, ridge uses filter coefficients $`\sigma_i/(\sigma_i^2+\lambda)`$ on the left-singular data coordinates. The pseudoinverse instead uses $`1/\sigma_i`$ on positive singular values. Ridge suppresses amplification in small-singular-value directions by accepting some residual error.

#### Bias and variance are separate effects

Suppose $`b=Ax_0+\epsilon`$ with $`E[\epsilon]=0`$ and covariance $`\sigma_\epsilon^2I`$. The ridge estimator has expectation

```math
E[\hat x_\lambda]
=(A^\top A+\lambda I)^{-1}A^\top Ax_0.
```

In a right-singular direction, the true coefficient is multiplied by $`\sigma_i^2/(\sigma_i^2+\lambda)`$. This produces shrinkage bias. The sampling variance in that direction is

```math
\sigma_\epsilon^2\frac{\sigma_i^2}{(\sigma_i^2+\lambda)^2}.
```

The mean squared error is squared bias plus variance, summed over orthogonal directions. A null-space component of $`x_0`$ is not identified by the data; ridge sets it to zero and retains its squared contribution as bias. Increasing the penalty reduces variance but need not reduce total error for every possible true vector.

#### Example LA-E22: penalizing differences instead of levels

Let $`A=(1,1)`$ and $`L=(1,-1)`$. Minimize $`(x_1+x_2-b)^2+\lambda(x_1-x_2)^2`$.

Solution. The data identify the sum, while the penalty controls the difference. Their kernels are the opposite lines and intersect only at zero. The unique solution is $`x_1=x_2=b/2`$ for every positive $`\lambda`$, with objective zero. The normal matrix is

```math
\begin{pmatrix}1+\lambda&1-\lambda\\[0pt] 1-\lambda&1+\lambda\end{pmatrix},
```

whose eigenvalues are two in the common direction and $`2\lambda`$ in the difference direction. If $`L`$ were instead equal to $`A`$, both terms would ignore the difference direction and uniqueness would fail. The direction penalized matters as much as the size of the penalty.

### 5.3 Worked examples

#### Example LA-E13: a sensitive direction

For $`A=\mathop{\mathrm{diag}} (1,10^{-6})`$ and $`b=(1,10^{-6})`$, the solution is $`(1,1)`$. Change only the second coordinate of $`b`$ by $`10^{-6}`$.

Solution. The new solution is $`(1,2)`$. An input perturbation of norm $`10^{-6}`$ has produced an output perturbation of norm one. The condition number is $`10^6`$. With $`\lambda=10^{-4}`$, the second ridge coordinate is $`10^{-12}/(10^{-12}+10^{-4})`$, about $`10^{-8}`$. It is insensitive but strongly shrunk relative to the unregularized value.

#### Example LA-E14: covariance shrinkage

Let $`S=Q\mathop{\mathrm{diag}} (4,0)Q^\top`$. Add $`\lambda I`$.

Solution. The eigenvectors stay the same and the eigenvalues become $`4+\lambda,\lambda`$. For $`\lambda=1/10`$, the condition number is $`4.1/0.1=41`$. Adding a scalar multiple of the identity makes a PSD covariance matrix invertible. Shrinkage toward a different target need not preserve the original eigenvectors.

#### Example LA-E15: a small residual with a large error

Use $`A=\mathop{\mathrm{diag}} (1,\varepsilon)`$, $`b=(1,0)`$, and an approximate solution $`\widetilde x=(1,1)`$.

Solution. The exact solution is $`(1,0)`$. The error has norm one, but the residual $`b-A\widetilde x=(0,-\varepsilon)`$ has norm $`\varepsilon`$. The relation $`x-\widetilde x=A^{-1}r`$ explains this. Reporting a tiny residual without the conditioning can overstate accuracy.

### 5.4 Selected exercises

#### Problem LA-14: reverse triangle inequality for operators

Source: [Axler, Exercises 7F, 1](https://linear.axler.net/LADR4e.pdf).

Problem. Prove the reverse triangle inequality $`\lvert\Vert S\Vert -\Vert T\Vert \rvert\le\Vert S-T\Vert`$ for induced operator norms.

Solution. For unit $`x`$, $`\Vert Sx\Vert \le\Vert (S-T)x\Vert +\Vert Tx\Vert \le\Vert S-T\Vert +\Vert T\Vert`$. Taking the supremum gives $`\Vert S\Vert -\Vert T\Vert \le\Vert S-T\Vert`$. Interchanging $`S,T`$ gives the opposite signed inequality. Combining them proves $`\vert \Vert S\Vert -\Vert T\Vert \vert \le\Vert S-T\Vert`$.

#### Problem LA-15: stability of invertibility

Source: [Axler, Exercises 7F, 8(a,b)](https://linear.axler.net/LADR4e.pdf).

Problem. First show that $`\Vert I-T\Vert \lt 1`$ implies invertibility. Then prove that every sufficiently small operator-norm perturbation of an invertible $`S`$ remains invertible, giving an explicit radius.

Solution. If $`\Vert I-T\Vert \lt 1`$, the geometric operator series $`R=\sum_{k=0}^\infty(I-T)^k`$ converges in norm. Multiplying its finite partial sum by $`T`$ gives $`I-(I-T)^{N+1}`$, which tends to $`I`$. Hence $`RT=TR=I`$.

For invertible $`S`$, factor $`T=S[I+S^{-1}(T-S)]`$. The hypothesis $`\Vert T-S\Vert \lt 1/\Vert S^{-1}\Vert`$ makes the bracket's difference from $`I`$ have norm below one. Part (a) makes the bracket invertible, and the product is invertible. This gives a concrete neighborhood of $`S`$ containing no singular matrices.

#### Problem LA-16: positive definiteness under perturbation

Source: [Axler, Exercises 7F, 16](https://linear.axler.net/LADR4e.pdf).

Problem. Let $`S`$ be a positive-definite self-adjoint operator. Find an operator-norm neighborhood within the self-adjoint operators in which positivity is preserved.

Solution. Let $`S`$ be PD, let $`m=\lambda_{\min}(S)\gt 0`$, and let $`T`$ be symmetric with $`\Vert T-S\Vert _2\lt m`$. For nonzero $`x`$,

```math
x^\top Tx=x^\top Sx+x^\top(T-S)x
\ge(m-\Vert T-S\Vert _2)\Vert x\Vert ^2\gt 0.
```

Thus one can choose $`\delta=m`$, and in fact nearby $`T`$ is PD as well as PSD. Symmetry is needed to call $`T`$ a positive operator; a small nonsymmetric perturbation can preserve a positive quadratic form without producing a self-adjoint operator.

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

[Contents](README.md) · [Part 1](01-linear-algebra-part-1.md) · [Part 3](01-linear-algebra-part-3.md)
