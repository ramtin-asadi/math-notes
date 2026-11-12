# Linear Algebra — Part 3

These notes use real finite-dimensional spaces unless a different field is stated. Vectors are columns, $`A^\top`$ is the transpose, $`I`$ is the identity, and $`\Vert x\Vert _2=(x^\top x)^{1/2}`$. For complex vectors, replace transposes by conjugate transposes and use the Hermitian inner product.

Sources:

- Gilbert Strang, [MIT 18.06SC, Fall 2011](https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/), especially Units II and III.
- Sheldon Axler, [Linear Algebra Done Right, fourth edition](https://linear.axler.net/LADR4e.pdf), Chapters 2, 3, 6, and 7. Exercise numbers below refer to this edition.

The explanations and solutions are written for these notes. Source problems are briefly paraphrased and linked; examples and additional problems are constructed here.

[Contents](README.md) · [Part 1](01-linear-algebra-part-1.md) · [Part 2](01-linear-algebra-part-2.md)

Contents:

- [7. Further worked problems](#7-further-worked-problems)
- [7.1 Weighted projection onto a line](#71-weighted-projection-onto-a-line)
- [7.2 Least squares with a visible residual](#72-least-squares-with-a-visible-residual)
- [7.3 Common and opposite directions in a quadratic form](#73-common-and-opposite-directions-in-a-quadratic-form)
- [7.4 Rank-one least squares and minimum norm](#74-rank-one-least-squares-and-minimum-norm)
- [7.5 Curvature along a constraint](#75-curvature-along-a-constraint)
- [7.6 A small residual in a weak direction](#76-a-small-residual-in-a-weak-direction)

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

[Contents](README.md) · [Part 1](01-linear-algebra-part-1.md) · [Part 2](01-linear-algebra-part-2.md)
