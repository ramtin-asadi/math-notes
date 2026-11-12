# Linear Algebra — Part 1

These notes use real finite-dimensional spaces unless a different field is stated. Vectors are columns, $`A^\top`$ is the transpose, $`I`$ is the identity, and $`\Vert x\Vert _2=(x^\top x)^{1/2}`$. For complex vectors, replace transposes by conjugate transposes and use the Hermitian inner product.

Sources:

- Gilbert Strang, [MIT 18.06SC, Fall 2011](https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/), especially Units II and III.
- Sheldon Axler, [Linear Algebra Done Right, fourth edition](https://linear.axler.net/LADR4e.pdf), Chapters 2, 3, 6, and 7. Exercise numbers below refer to this edition.

The explanations and solutions are written for these notes. Source problems are briefly paraphrased and linked; examples and additional problems are constructed here.

[Contents](README.md) · [Part 2](01-linear-algebra-part-2.md) · [Part 3](01-linear-algebra-part-3.md)

Contents:

- [1. Vector spaces and linear maps](#1-vector-spaces-and-linear-maps)
- [2. Orthogonality and projections](#2-orthogonality-and-projections)
- [3. Eigenvalues and quadratic forms](#3-eigenvalues-and-quadratic-forms)

## 1. Vector spaces and linear maps

### 1.1 Theory

A vector space is a set closed under vector addition and scalar multiplication, with the usual associativity, distributivity, identity, and inverse rules. A subspace must contain zero and be closed under both operations. An affine plane such as $`x_1+x_2+x_3=1`$ is therefore not a subspace: it misses zero.

The span of $`v_1,\ldots,v_k`$ consists of all vectors $`\sum_i a_iv_i`$. The list is linearly independent when $`\sum_i a_iv_i=0`$ forces every $`a_i=0`$. A basis is an independent spanning list. Independence gives unique coordinates; spanning ensures that coordinates exist. All bases of a finite-dimensional space have the same length, called its dimension.

A map $`T:V\to W`$ is linear if $`T(ax+by)=aT(x)+bT(y)`$. Its kernel and range are

```math
\ker T=\lbrace v:T(v)=0\rbrace ,\qquad
\mathop{\mathrm{range}} T=\lbrace T(v):v\in V\rbrace .
```

Both are subspaces. The rank is the dimension of the range; the nullity is the dimension of the kernel. A matrix represents a map after bases are chosen. Changing the bases changes the matrix while leaving the underlying map unchanged.

For $`A\in\mathbb R^{m\times n}`$, the column space is its range. Solving $`Ax=b`$ asks whether $`b`$ lies in that space. If $`x_0`$ is one solution, every solution has the form $`x_0+z`$ with $`z\in\ker A`$: subtract two solutions to see why.

### 1.2 Derivations and proofs

#### Reading a linear combination

The equation $`b=c_1a_1+\cdots+c_na_n`$ asks whether the available directions $`a_i`$ can produce the target $`b`$. The coefficients can be positive, negative, or zero. A zero coefficient means that direction is unused; a negative one reverses it. The coefficients are the unknowns, while the vectors are the given information.

For example, if $`a_1=(1,1)`$ and $`a_2=(1,-1)`$, the equation $`(3,1)=c_1a_1+c_2a_2`$ becomes $`c_1+c_2=3`$, $`c_1-c_2=1`$. Adding and subtracting gives $`c_1=2`$, $`c_2=1`$. The target is a vector in the original space; $`(2,1)`$ is its coordinate vector relative to this particular basis. Using another basis changes its coordinates.

A spanning list may contain an unnecessary direction. If $`a_3=a_1+a_2`$, every use of $`a_3`$ can be replaced by uses of $`a_1`$ and $`a_2`$. A basis removes all such redundancy. To test whether a list is a basis, check both parts of the definition. A list can be independent but too short to span, or span but be dependent because it is too long.

#### Testing subspaces and linearity

For a proposed subspace, the closure test is: take arbitrary $`u,v`$ in the set and arbitrary real numbers $`a,b`$, and prove $`au+bv`$ is still in the set. This proves addition and scalar closure together. Check zero as well, since the empty set would otherwise pass a statement quantified over its elements.

For a set defined by homogeneous linear equations, the test is immediate. If $`Au=Av=0`$, then $`A(au+bv)=aAu+bAv=0`$. With a nonzero right-hand side $`c`$, the same computation gives $`(a+b)c`$, which need not equal $`c`$.

For a proposed linear map, testing two or three numerical inputs cannot establish linearity. Write the expression at general $`au+bv`$ and simplify. To disprove linearity, one counterexample suffices. In particular every linear map satisfies $`T(0)=0`$, because $`T(0)=T(0+0)=2T(0)`$.

The map $`T(x)=Ax+c`$ is linear exactly when $`c=0`$. When $`c\ne0`$ it is affine. The map $`T(x_1,x_2)=x_1x_2`$ sends zero to zero but is still nonlinear: $`T(e_1+e_2)=1`$, whereas $`T(e_1)+T(e_2)=0`$.

#### Extending an independent list to a basis

Start with an independent list in a finite-dimensional space. If it does not span, choose a vector outside its span and append it. The longer list remains independent: a dependence with a nonzero coefficient on the new vector would express that vector in the old span. If its new coefficient were zero, the original independence would force all remaining coefficients to vanish.

Continue until the list spans. The process stops because a finite-dimensional space cannot contain arbitrarily long independent lists. This construction explains the basis-extension step in rank-nullity. It also explains why the theorem applies to arbitrary linear maps, rather than only to matrices already in row-echelon form.

#### Four spaces belonging to one matrix

For $`A\in\mathbb R^{m\times n}`$, keep track of the ambient dimensions:

| Space | Ambient space | Dimension |
| --- | --- | --- |
| Column space $`\mathop{\mathrm{range}} A`$ | $`\mathbb R^m`$ | $`r`$ |
| Left null space $`\ker A^\top`$ | $`\mathbb R^m`$ | $`m-r`$ |
| Row space $`\mathop{\mathrm{range}} A^\top`$ | $`\mathbb R^n`$ | $`r`$ |
| Null space $`\ker A`$ | $`\mathbb R^n`$ | $`n-r`$ |

Column space and left null space are orthogonal complements. Row space and null space are orthogonal complements. A vector in $`\mathbb R^n`$ cannot generally be compared with a column-space vector in $`\mathbb R^m`$ when $`m\ne n`$; they live in different spaces.

The equality of row rank and column rank follows, for example, from elimination: multiplying by an invertible row-operation matrix preserves dependencies among columns, and a row-echelon matrix has one independent row and one pivot column for each pivot. Rank-nullity then supplies both null-space dimensions.

#### Rank-nullity

Let $`T:V\to W`$ be linear and $`\dim V=n`$. Choose a basis $`u_1,\ldots,u_k`$ of $`\ker T`$ and extend it to a basis

```math
u_1,\ldots,u_k,v_1,\ldots,v_{n-k}
```

of $`V`$. Every vector in $`V`$ has an expansion in this basis. Applying $`T`$ removes the $`u_i`$ terms, so $`Tv_1,\ldots,Tv_{n-k}`$ span the range.

To check independence, suppose $`\sum_j a_jTv_j=0`$. Then $`\sum_j a_jv_j\in\ker T`$, so it is also a linear combination of the $`u_i`$. Independence of the extended basis forces all $`a_j=0`$. The images are consequently a basis of the range. Thus

```math
\dim V=\dim\ker T+\dim\mathop{\mathrm{range}} T.
```

The proof also constructs a complement $`U=\mathop{\mathrm{span}} (v_1,\ldots,v_{n-k})`$ on which $`T`$ is one-to-one and maps onto its full range.

#### Row space and null space

The equation $`Ax=0`$ says that $`x`$ is perpendicular to every row of $`A`$. Hence $`\ker A=\mathop{\mathrm{range}} (A^\top)^\perp`$. In finite dimensions, taking orthogonal complements again gives

```math
(\ker A)^\perp=\mathop{\mathrm{range}} (A^\top).
```

This identifies the part of a parameter vector that the map can detect. Adding a null-space vector changes the input without changing the output.

#### Coordinates and a change of basis

A coordinate vector records coefficients in a chosen basis. It is not the underlying geometric vector independently of that basis. Let $`P`$ have the new domain basis vectors as columns in old coordinates, and let $`Q`$ do the same for the codomain. A vector with new domain coordinates $`z`$ has old coordinates $`Pz`$. Applying the old matrix $`A`$ gives $`APz`$, and converting the output to new coordinates gives $`Q^{-1}APz`$. The new matrix is therefore

```math
A_{\mathrm{new}}=Q^{-1}AP.
```

Both change-of-basis matrices are invertible because their columns are bases. Left and right multiplication by invertible maps preserve rank, so rank is independent of coordinates. The kernel's coordinate description changes, but its dimension does not.

For an operator from a space to itself with the same new basis on input and output, this becomes the similarity transformation $`P^{-1}AP`$. Similar matrices have the same eigenvalues: if $`Av=\lambda v`$, then $`P^{-1}v`$ is an eigenvector of the new matrix. A general two-sided coordinate transformation with unrelated $`P,Q`$ need not preserve eigenvalues.

There is another way to read the kernel. Two inputs $`x,y`$ give the same output exactly when $`x-y\in\ker T`$. Thus each attainable output corresponds to an entire translated kernel, $`x+\ker T`$. Rank-nullity counts the input directions that are distinguishable after applying the map and those that are lost. This interpretation also explains why adding null-space vectors cannot improve or worsen a least-squares fit.

#### Example LA-E19: the same operator in two bases

Let $`A=\mathop{\mathrm{diag}} (2,3)`$ in the standard basis of $`\mathbb R^2`$. Use the new basis $`p_1=(1,1)^\top`$, $`p_2=(1,-1)^\top`$.

Solution. The basis matrix is

```math
P=\begin{pmatrix}1&1\\[0pt] 1&-1\end{pmatrix}
```

, with inverse $`P/2`$. Multiplication gives

```math
P^{-1}AP=\frac12
\begin{pmatrix}5&-1\\[0pt] -1&5\end{pmatrix}.
```

The first column says that the old vector $`Ap_1=(2,3)^\top`$ has new coordinates $`(5/2,-1/2)^\top`$. Reconstructing it as $`(5/2)p_1-(1/2)p_2`$ checks the interpretation. The new matrix has eigenvalues two and three, as expected, even though it is no longer diagonal. A diagonal appearance depends on coordinates; the operator's eigenvalues do not.

### 1.3 Worked examples

#### Example LA-E01: a plane and its coordinates

Find a basis for $`S=\lbrace (x,y,z):x+y+z=0\rbrace`$.

Solution. Solve for $`z=-x-y`$:

```math
(x,y,z)=x(1,0,-1)+y(0,1,-1).
```

The two displayed vectors span $`S`$. If their linear combination is zero, its first two coordinates force both coefficients to vanish. They are a basis and $`\dim S=2`$. The equation defining the plane has removed one free coordinate from $`\mathbb R^3`$.

#### Example LA-E02: a kernel with a free variable

Let

```math
A=\begin{pmatrix}1&2&3\\[0pt] 0&1&1\end{pmatrix}.
```

Solution. The second row of $`Ax=0`$ gives $`x_2=-x_3`$; the first then gives $`x_1=-x_3`$. Thus $`\ker A=\mathop{\mathrm{span}} ((-1,-1,1))`$. The first two columns are independent, so $`\mathop{\mathrm{rank}} A=2`$. Rank-nullity reads $`3=1+2`$. For $`b=(4,1)^\top`$, setting $`x_3=t`$ gives all solutions $`(2-t,1-t,t)^\top`$.

#### Example LA-E03: differentiation as a linear map

Let $`D:\mathcal P_3\to\mathcal P_2`$ send $`p`$ to $`p'`$, where $`\mathcal P_k`$ means polynomials of degree at most $`k`$.

Solution. Differentiation respects sums and scalar multiples. Its kernel consists of constant polynomials, so its nullity is one. For any $`a+bx+cx^2`$, the polynomial $`ax+(b/2)x^2+(c/3)x^3`$ maps to it, so the range is all of $`\mathcal P_2`$ and the rank is three. Again $`4=1+3`$.

### 1.4 Selected exercises

#### Problem LA-01: a redundant spanning list

Source: [Axler, Exercises 2A, 1](https://linear.axler.net/LADR4e.pdf). Find four distinct vectors spanning the plane from Example LA-E01.

Solution. Use

```math
v_1=(1,0,-1),\quad v_2=(0,1,-1),\quad
v_3=(1,1,-2),\quad v_4=(1,-1,0).
```

Each coordinate sum is zero, so their span lies in the plane. Conversely the first two already span it, as shown in the example. The other two satisfy $`v_3=v_1+v_2`$ and $`v_4=v_1-v_2`$. The list has four members but rank two; a spanning list need not be a basis.

#### Problem LA-02: cumulative sums preserve span

Source: [Axler, Exercises 2A, 3](https://linear.axler.net/LADR4e.pdf). Compare the span of $`v_1,\ldots,v_m`$ with that of $`w_k=\sum_{i=1}^k v_i`$.

Solution. Every $`w_k`$ is a combination of the $`v_i`$, giving one containment. For the other, $`v_1=w_1`$ and $`v_k=w_k-w_{k-1}`$ for $`k\ge2`$. Each original vector is therefore in the new span. Both spans agree. The cumulative-sum transformation is triangular with diagonal entries one, and the differences explicitly give its inverse.

#### Problem LA-03: surjectivity from a specified kernel

Source: [Axler, Exercises 3B, 12](https://linear.axler.net/LADR4e.pdf). A map $`T:\mathbb F^4\to\mathbb F^2`$ has kernel determined by $`x_1=5x_2`$ and $`x_3=7x_4`$. Prove surjectivity.

Solution. Write every kernel vector as

```math
s(5,1,0,0)+t(0,0,7,1).
```

These two vectors are independent, so the nullity is two. Rank-nullity gives rank $`4-2=2`$. The range is a two-dimensional subspace of a two-dimensional codomain, and therefore equals that codomain. Merely knowing that two equations define the kernel would be insufficient if those equations were redundant; the explicit basis checks their independence.

#### Problem LA-04: decomposition for an idempotent

Source: [Axler, Exercises 3B, 27](https://linear.axler.net/LADR4e.pdf). For $`P^2=P`$, prove $`V=\ker P\oplus\mathop{\mathrm{range}} P`$.

Solution. For every $`v`$,

```math
v=(v-Pv)+Pv,\qquad P(v-Pv)=Pv-P^2v=0.
```

The first term is in the kernel and the second is in the range. If $`z`$ belongs to both, write $`z=Py`$ and use $`Pz=P^2y=Py=z`$. Kernel membership also says $`Pz=0`$, so $`z=0`$. The sum is direct. Idempotence alone does not make this an orthogonal decomposition; symmetry supplies that additional property.

## 2. Orthogonality and projections

### Further problems on spaces and maps

#### Problem LA-21: when changing one vector preserves independence

Source: [Axler, Exercises 2A, 9](https://linear.axler.net/LADR4e.pdf).

Given independent $`v_1,\ldots,v_m`$, replace $`v_1`$ by $`5v_1-4v_2`$, assuming $`m\ge2`$.

Solution. Suppose $`a_1(5v_1-4v_2)+a_2v_2+\cdots+a_mv_m=0`$. Collect coefficients of the original vectors:

```math
5a_1v_1+(a_2-4a_1)v_2+\sum_{j=3}^m a_jv_j=0.
```

Original independence gives $`5a_1=0`$, hence $`a_1=0`$, and then $`a_2=0`$ and all later coefficients zero. The replacement preserves independence over $`\mathbb R`$ or $`\mathbb C`$. The nonzero coefficient of the replaced vector makes the change reversible.

#### Problem LA-22: sums of independent lists

Source: [Axler, Exercises 2A, 11](https://linear.axler.net/LADR4e.pdf).

Problem. Determine whether adding corresponding vectors of two linearly independent lists must produce another independent list. Prove the claim or give a counterexample.

Solution. The assertion that coordinatewise sums of two independent lists must remain independent is false. Choose any independent $`v_1,\ldots,v_m`$ and let $`w_i=-v_i`$. The $`w_i`$ are independent, since a dependence among them is the negative of a dependence among the $`v_i`$. Yet each sum $`v_i+w_i=0`$, so the new list is dependent. Independence concerns relationships within one list, and gives no restriction on cancellation between two different lists.

#### Problem LA-23: preserving independent vectors

Source: [Axler, Exercises 3B, 9](https://linear.axler.net/LADR4e.pdf).

Problem. Prove that an injective linear map takes every independent list to an independent list, and examine the converse.

Solution. Let $`T`$ be injective and suppose $`\sum_i a_iTv_i=0`$. By linearity this is $`T(\sum_i a_iv_i)=0`$. Injectivity implies $`\sum_i a_iv_i=0`$, and independence of the $`v_i`$ implies $`a_i=0`$ for every $`i`$. Thus injective linear maps preserve independence. The converse also holds if a map preserves every independent one-vector list: a nonzero kernel vector would be sent to zero, violating that property.

#### Problem LA-24: rank of a composition

Source: [Axler, Exercises 3B, 23](https://linear.axler.net/LADR4e.pdf).

Problem. For compatible finite-dimensional linear maps $`T`$ and $`S`$, prove that the rank of $`ST`$ is at most the rank of either factor.

Solution. For $`T:U\to V`$ and $`S:V\to W`$, $`\mathop{\mathrm{range}} (ST)=S(\mathop{\mathrm{range}} T)`$. Its dimension cannot exceed $`\dim\mathop{\mathrm{range}} T`$, since images of a basis of that space span the image. It is also a subspace of $`\mathop{\mathrm{range}} S`$, so its dimension cannot exceed $`\mathop{\mathrm{rank}} S`$. Therefore

```math
\mathop{\mathrm{rank}} (ST)\le\min(\mathop{\mathrm{rank}} S,\mathop{\mathrm{rank}} T).
```

Equality need not hold: $`T`$ can map entirely into $`\ker S`$, in which case the composition is zero even though neither map is zero.

#### Problem LA-25: a parameter changes rank

Constructed problem. For

```math
A_t=\left(\begin{smallmatrix}1&1\\[0pt] 1&t\end{smallmatrix}\right)
```

, determine rank, kernel, and solvability for every $`t`$.

Solution. Subtract the first row from the second. The reduced second row is $`(0,t-1)`$. If $`t\ne1`$, both coordinates are pivots, the rank is two, the kernel is zero, and every right-hand side has a unique solution.

If $`t=1`$, the rank is one and the kernel is $`\mathop{\mathrm{span}} ((1,-1))`$. The equation $`A_1x=b`$ is solvable exactly when $`b_1=b_2`$. If that equality holds, all solutions are $`(b_1-s,s)`$ for $`s\in\mathbb R`$. Thus one parameter value changes both the number of admissible targets and the uniqueness of coefficients.

#### Problem LA-26: the matrix after a change of basis

Constructed problem. Let $`T(x,y)=(2x+y,x+2y)`$ and use the basis $`b_1=(1,1)`$, $`b_2=(1,-1)`$.

Solution. The basis matrix is

```math
B=\left(\begin{smallmatrix}1&1\\[0pt] 1&-1\end{smallmatrix}\right)
```

 and $`B^{-1}=\tfrac12B`$. For an input coordinate vector $`c`$, the original vector is $`Bc`$, its image is $`ABc`$, and the coordinates of the image are $`B^{-1}ABc`$. Hence the new matrix is $`B^{-1}AB`$.

Here $`Tb_1=3b_1`$ and $`Tb_2=b_2`$, so the new matrix is $`\mathop{\mathrm{diag}} (3,1)`$. The map has not changed; the new coordinates align with directions that it scales without mixing.

#### Problem LA-27: an inconsistent system

Constructed problem. Determine all least-squares targets that are exactly solvable for

```math
A=\left(\begin{smallmatrix}1&0\\[0pt] 0&1\\[0pt] 1&1\end{smallmatrix}\right)
```

.

Solution. The range consists of vectors $`(s,t,s+t)`$, so $`Ax=b`$ is exactly solvable if and only if $`b_3=b_1+b_2`$. A left-null vector is $`z=(-1,-1,1)`$, and the same condition is $`z^\top b=0`$. If the condition fails, the component of $`b`$ along $`z`$ is unreachable. Section 2 shows that least squares removes that component by orthogonal projection.

#### Problem LA-28: images and preimages

Constructed problem. Let $`T(x,y,z)=(x+y,y+z)`$ and find the preimage of the line $`\lbrace (s,s):s\in\mathbb R\rbrace`$.

Solution. Membership requires $`x+y=y+z`$, hence $`x=z`$. The preimage is $`\lbrace (x,y,x)\rbrace`$, with basis $`(1,0,1)`$ and $`(0,1,0)`$ and dimension two. The kernel of $`T`$ is generated by $`(-1,1,-1)`$, and its dimension is one. The target line contributes one additional independent image direction, giving $`2=1+1`$. A preimage under a linear map of a subspace is itself a subspace.

### 2.1 Theory

An inner product measures angles and lengths. In $`\mathbb R^n`$ it is $`x^\top y`$; orthogonality means $`x^\top y=0`$. The orthogonal complement $`S^\perp`$ contains the vectors perpendicular to every element of $`S`$.

An orthogonal projection $`P_Sb`$ is the element of $`S`$ closest to $`b`$. It is characterized by

```math
P_Sb\in S,\qquad b-P_Sb\in S^\perp.
```

If the columns of $`Q`$ are an orthonormal basis of $`S`$, then $`P_S=QQ^\top`$. If instead $`S`$ is represented by independent columns of $`A`$, then $`P_S=A(A^\top A)^{-1}A^\top`$. Neither formula requires $`S`$ to be all of the ambient space.

Gram-Schmidt constructs orthonormal directions from independent vectors:

```math
u_k=a_k-\sum_{j\lt k}(q_j^\top a_k)q_j,\qquad q_k=u_k/\Vert u_k\Vert .
```

Independence guarantees $`u_k\ne0`$. In numerical work, Householder QR or reorthogonalization is often more reliable than the unmodified procedure; see [Numerical Methods](07-numerical-methods-part-1.md#2-numerical-linear-algebra).

### 2.2 Derivations and proofs

#### Finding the closest point one coefficient at a time

For a line generated by nonzero $`a`$, every candidate is $`ca`$. Its squared error from $`b`$ is

```math
g(c)=\Vert b-ca\Vert ^2
=b^\top b-2c\,a^\top b+c^2a^\top a.
```

This is an upward-opening scalar quadratic because $`a^\top a\gt 0`$. Setting $`g'(c)=0`$ gives $`c=a^\top b/(a^\top a)`$. Substituting that coefficient into $`a^\top(b-ca)`$ gives zero. Thus the perpendicular residual arises from minimization; it is also a way to solve the minimization.

For several independent columns, the coefficients interact. If the columns are orthonormal, each coefficient is $`q_i^\top b`$ and can be found separately. If the columns are not orthogonal, the interaction is recorded by the Gram matrix $`A^\top A`$. Its $`(i,j)`$ entry is the inner product of columns $`i`$ and $`j`$.

The normal equations correct for that interaction. Writing $`\widehat x_i=a_i^\top b`$ for nonorthogonal columns counts overlapping directions more than once. Multiplication by $`(A^\top A)^{-1}`$ removes that overlap.

#### Cauchy-Schwarz and its equality case

For $`v\ne0`$, consider the nonnegative squared norm of $`u-tv`$ and choose $`t=\langle u,v\rangle/\Vert v\Vert ^2`$ in a real space. Expansion yields

```math
0\le\left\Vert u-\frac{\langle u,v\rangle}{\Vert v\Vert ^2}v\right\Vert ^2
=\Vert u\Vert ^2-\frac{\vert \langle u,v\rangle\vert ^2}{\Vert v\Vert ^2}.
```

Therefore $`\vert \langle u,v\rangle\vert \le\Vert u\Vert \Vert v\Vert`$. Equality holds exactly when the residual vector is zero, which means $`u`$ is a multiple of $`v`$. If $`v=0`$, both sides are zero. This proof is itself a projection calculation.

#### Existence in rank-deficient least squares

Dependent columns do not prevent a best fitted vector from existing. Their span is a finite-dimensional subspace and therefore closed. Choose an orthonormal basis $`Q`$ of that span and project $`b`$ to $`QQ^\top b`$. Because this fitted vector belongs to the column space, at least one coefficient vector produces it.

Two minimizing coefficient vectors differ by an element of $`\ker A`$. To choose one uniquely, require it to lie in $`(\ker A)^\perp`$. Every vector has a unique decomposition into these two orthogonal spaces; removing the kernel component leaves the fit unchanged and reduces the norm. This proves the minimum-norm characterization before any formula for a pseudoinverse is used.

#### Weighted least squares

Suppose $`W`$ is symmetric PD and squared error is measured by $`r^\top Wr`$ instead of $`r^\top r`$. An eigendecomposition gives a square root $`W^{1/2}`$, so

```math
(b-Ax)^\top W(b-Ax)=\Vert W^{1/2}b-W^{1/2}Ax\Vert ^2.
```

Applying ordinary least squares to the transformed system gives $`A^\top WA\widehat x=A^\top Wb`$. With independent columns,

```math
\widehat x=(A^\top WA)^{-1}A^\top Wb,\qquad
P_W=A(A^\top WA)^{-1}A^\top W.
```

The matrix $`P_W`$ is idempotent but generally not symmetric in Euclidean coordinates. Instead $`P_W^\top W=WP_W`$. Its residual is orthogonal to the fitted space under the weighted inner product $`\langle u,v\rangle_W=u^\top Wv`$. Always specify the geometry in which a projection is orthogonal.

#### Projection and least squares

Minimize $`\Vert b-Ax\Vert ^2`$ for $`A\in\mathbb R^{m\times n}`$. A minimizing residual must be perpendicular to each column: otherwise moving along that column decreases the residual length. Algebraically, differentiation gives

```math
A^\top(A\widehat x-b)=0.
```

These are the normal equations. If $`A`$ has independent columns, $`z^\top A^\top Az=\Vert Az\Vert ^2\gt 0`$ for $`z\ne0`$, so

```math
\widehat x=(A^\top A)^{-1}A^\top b,\qquad
P=A(A^\top A)^{-1}A^\top.
```

Transposing $`P`$ gives $`P^\top=P`$, because $`A^\top A`$ and its inverse are symmetric. Multiplying gives

```math
P^2=A(A^\top A)^{-1}(A^\top A)(A^\top A)^{-1}A^\top=P.
```

For any $`x`$, the orthogonal decomposition

```math
\Vert b-Ax\Vert ^2=\Vert b-Pb\Vert ^2+\Vert Pb-Ax\Vert ^2
```

proves global optimality. It also proves that the fitted vector is unique even when the coefficients are not. In the rank-deficient case use $`P=AA^+`$, where $`A^+`$ is the pseudoinverse introduced below.

#### Gram–Schmidt with every subtraction visible

Take independent vectors $`u_1=(1,1,0)^\top`$, $`u_2=(1,0,1)^\top`$, and $`u_3=(0,1,1)^\top`$. Normalize the first:

```math
q_1=\frac1{\sqrt2}(1,1,0)^\top.
```

Remove from $`u_2`$ the part already in the first direction:

```math
v_2=u_2-(q_1^\top u_2)q_1
=(1,0,1)^\top-\tfrac12(1,1,0)^\top
=\tfrac12(1,-1,2)^\top.
```

Its norm is $`\sqrt6/2`$, giving $`q_2=(1,-1,2)^\top/\sqrt6`$. For the third vector, the coefficients are $`q_1^\top u_3=1/\sqrt2`$ and $`q_2^\top u_3=1/\sqrt6`$. Therefore

```math
v_3=u_3-\tfrac12(1,1,0)^\top-\tfrac16(1,-1,2)^\top
=\tfrac23(-1,1,1)^\top.
```

Normalizing gives $`q_3=(-1,1,1)^\top/\sqrt3`$. Direct inner products check that the three $`q_i`$ are orthonormal. At each stage, the subtraction preserves the span built so far, while removing overlap with earlier directions.

Collecting $`Q=(q_1,q_2,q_3)`$ gives $`U=QR`$, with

```math
R=\begin{pmatrix}
\sqrt2&1/\sqrt2&1/\sqrt2\\[0pt]
0&\sqrt6/2&1/\sqrt6\\[0pt]
0&0&2/\sqrt3
\end{pmatrix}.
```

The entries above the diagonal are projection coefficients; the diagonal entries are residual norms. A zero residual would reveal dependence, so normalization at that step would be invalid.

#### Two successive projections

If orthogonal projection matrices $`P,Q`$ commute, their product is symmetric and idempotent:
$`(PQ)^\top=QP=PQ`$ and $`(PQ)^2=P^2Q^2=PQ`$. Its range is the intersection of the two ranges. Conversely, if their product is an orthogonal projection, symmetry forces $`PQ=QP`$.

For projection onto the first coordinate axis and projection onto the line $`(1,1)`$,

```math
P=\begin{pmatrix}1&0\\[0pt] 0&0\end{pmatrix},\qquad
Q=\tfrac12\begin{pmatrix}1&1\\[0pt] 1&1\end{pmatrix},\qquad
PQ=\tfrac12\begin{pmatrix}1&1\\[0pt] 0&0\end{pmatrix}.
```

The product is neither symmetric nor idempotent. Projecting onto one set and then another need not project onto their intersection in a single pass.

For a weighted inner product with $`W=C^\top C`$ invertible, the map $`z\mapsto Cz`$ converts weighted lengths into Euclidean lengths. Weighted least squares is ordinary least squares for $`CX`$ and $`Cy`$. Its projection in original coordinates is $`C^{-1}P_{CX}C`$, explaining why it is idempotent without generally being Euclidean-symmetric.

### 2.3 Worked examples

#### Example LA-E04: projection onto a line

Project $`b=(2,0,1)^\top`$ onto the line generated by $`a=(1,1,1)^\top`$.

Solution. The coefficient is $`(a^\top b)/(a^\top a)=3/3=1`$, so $`Pb=(1,1,1)^\top`$ and $`r=(1,-1,0)^\top`$. Directly $`a^\top r=0`$. The projection matrix is $`aa^\top/3`$, with every entry $`1/3`$. The squared distance is $`\Vert r\Vert ^2=2`$.

#### Example LA-E05: fitting a line to three observations

Fit $`y_i=a+bt_i`$ to $`(t,y)=(0,1),(1,2),(2,2)`$.

Solution. Here

```math
X=\begin{pmatrix}1&0\\[0pt] 1&1\\[0pt] 1&2\end{pmatrix},\quad
X^\top X=\begin{pmatrix}3&3\\[0pt] 3&5\end{pmatrix},\quad
X^\top y=\binom56.
```

Subtracting the first normal equation from the second gives $`2b=1`$, so $`b=1/2`$ and $`a=7/6`$. The fitted values are $`(7/6,5/3,13/6)`$ and the residuals are $`(-1/6,1/3,-1/6)`$. Their sum and their inner product with $`(0,1,2)`$ are both zero.

#### Example LA-E06: Gram-Schmidt in a plane

Orthonormalize $`a_1=(1,1,0)`$ and $`a_2=(1,0,1)`$.

Solution. First $`q_1=(1,1,0)/\sqrt2`$. Subtract the component $`(q_1^\top a_2)q_1=(1/2,1/2,0)`$:

```math
u_2=(1/2,-1/2,1),\qquad q_2=(1,-1,2)/\sqrt6.
```

Both vectors have unit length and their inner product is zero. The projection onto their span is $`q_1q_1^\top+q_2q_2^\top`$. A normal direction is $`(1,-1,-1)`$, so the same projection is $`I-nn^\top/3`$.

### 2.4 Selected exercises

#### Problem LA-05: projecting onto the orthogonal complement

Source: [Axler, Exercises 6C, 5](https://linear.axler.net/LADR4e.pdf).

Problem. For a subspace $`U`$ of a finite-dimensional inner-product space, express the orthogonal projection onto $`U^\perp`$ in terms of the projection onto $`U`$ and prove the identity.

Solution. Write $`v=u+w`$ with $`u\in U`$ and $`w\in U^\perp`$. By the projection characterization, $`P_Uv=u`$ and $`P_{U^\perp}v=w`$. Hence $`P_{U^\perp}v=v-P_Uv`$ for every $`v`$, and $`P_{U^\perp}=I-P_U`$. This also gives $`(I-P_U)^2=I-P_U`$ and $`(I-P_U)^\top=I-P_U`$.

#### Problem LA-06: recognizing an orthogonal projection

Source: [Axler, Exercises 6C, 9](https://linear.axler.net/LADR4e.pdf). Assume $`P^2=P`$ and $`\ker P\perp\mathop{\mathrm{range}} P`$.

Solution. Problem LA-04 gives $`v=(v-Pv)+Pv`$, with the terms in the kernel and range respectively. The assumed orthogonality makes this the orthogonal decomposition relative to $`U=\mathop{\mathrm{range}} P`$. Therefore $`Pv=P_Uv`$. For a direct minimization check, any $`u\in U`$ satisfies

```math
\Vert v-u\Vert ^2=\Vert v-Pv\Vert ^2+\Vert Pv-u\Vert ^2.
```

The second term is nonnegative and vanishes only at $`u=Pv`$. Thus $`P`$ selects the unique closest point.

#### Problem LA-07: a closest vector in four dimensions

Source: [Axler, Exercises 6C, 15](https://linear.axler.net/LADR4e.pdf). Minimize distance to $`(1,2,3,4)`$ over $`U=\mathop{\mathrm{span}} ((1,1,0,0),(1,1,1,2))`$.

Solution. Replace the generators by the orthogonal vectors $`a=(1,1,0,0)`$ and $`c=(0,0,1,2)`$; their span is unchanged. Then

```math
P_Ub=\frac{b^\top a}{a^\top a}a+
\frac{b^\top c}{c^\top c}c
=\frac32a+\frac{11}{5}c
=\left(\frac32,\frac32,\frac{11}{5},\frac{22}{5}\right).
```

The residual is $`(-1/2,1/2,4/5,-2/5)`$. Its inner products with $`a`$ and $`c`$ both vanish. The minimum squared distance is $`1/2+16/25+4/25=13/10`$. The Pythagorean identity above proves uniqueness.

## 3. Eigenvalues and quadratic forms

### Further problems on orthogonality and least squares

#### Problem LA-29: an inner product produced by a map

Source: [Axler, Exercises 6A, 2](https://linear.axler.net/LADR4e.pdf).

Problem. Let $`S`$ be a linear map between inner-product spaces. Determine exactly when $`\langle u,v\rangle_1=\langle Su,Sv\rangle`$ defines an inner product on its domain.

Solution. Define $`\langle u,v\rangle_1=\langle Su,Sv\rangle`$. Symmetry and linearity follow from the original inner product and linearity of $`S`$. Nonnegativity follows from $`\langle v,v\rangle_1=\Vert Sv\Vert ^2\ge0`$. This equals zero exactly when $`Sv=0`$. Thus positive definiteness holds exactly when $`\ker S=\lbrace 0\rbrace`$, that is, when $`S`$ is injective. For a matrix $`S`$, the induced Gram matrix is $`S^\top S`$, recovering the earlier full-column-rank criterion.

#### Problem LA-30: an orthogonal decomposition in two coordinates

Source: [Axler, Exercises 6A, 11](https://linear.axler.net/LADR4e.pdf).

Problem. Write $`(1,2)`$ as the sum of a vector parallel to $`(1,3)`$ and a vector perpendicular to $`(1,3)`$.

Solution. Decompose $`b=(1,2)`$ into a multiple of $`a=(1,3)`$ and a perpendicular residual. The projection coefficient is $`a^\top b/(a^\top a)=7/10`$. Therefore

```math
u=(7/10,21/10),\qquad v=b-u=(3/10,-1/10).
```

Their sum is $`b`$ and $`a^\top v=3/10-3/10=0`$. Also $`\Vert u\Vert ^2=49/10`$ and $`\Vert v\Vert ^2=1/10`$, whose sum equals $`\Vert b\Vert ^2=5`$. This checks both the arithmetic and the orthogonal decomposition.

#### Problem LA-31: Bessel's inequality and equality

Source: [Axler, Exercises 6B, 3](https://linear.axler.net/LADR4e.pdf).

Problem. For an orthonormal list $`e_1,\ldots,e_m`$, prove Bessel's inequality and characterize equality for a vector $`v`$.

Solution. For orthonormal $`e_1,\ldots,e_m`$, put $`p=\sum_i\langle v,e_i\rangle e_i`$ and $`r=v-p`$. Taking the inner product with any $`e_j`$ gives $`\langle r,e_j\rangle=0`$, so

```math
\Vert v\Vert ^2=\sum_i\vert \langle v,e_i\rangle\vert ^2+\Vert r\Vert ^2.
```

The sum of coefficient squares is at most $`\Vert v\Vert ^2`$. Equality holds exactly when $`r=0`$, or $`v\in\mathop{\mathrm{span}} (e_1,\ldots,e_m)`$. No assumption that the orthonormal list spans the whole ambient space is needed.

#### Problem LA-32: orthogonality between two subspaces

Source: [Axler, Exercises 6C, 7](https://linear.axler.net/LADR4e.pdf).

Problem. For subspaces $`X,Y`$, prove that every vector of $`X`$ is orthogonal to every vector of $`Y`$ exactly when $`P_XP_Y=0`$.

Solution. If every $`x\in X`$ is perpendicular to every $`y\in Y`$, then $`P_Yv\in Y\subseteq X^\perp`$, so $`P_XP_Yv=0`$ for every $`v`$.

Conversely suppose $`P_XP_Y=0`$. For $`y\in Y`$, $`P_Yy=y`$, hence $`P_Xy=0`$. The kernel of $`P_X`$ is $`X^\perp`$, so $`y\perp X`$. This holds for every $`y`$, establishing the claimed equivalence. It also implies $`P_YP_X=0`$ by transposition.

#### Problem LA-33: a weighted average

Constructed problem. Fit a common mean $`c`$ to observations $`1,3,5`$, with positive weights $`1,2,1`$.

Solution. Minimize $`(1-c)^2+2(3-c)^2+(5-c)^2`$. Differentiate:

```math
2(c-1)+4(c-3)+2(c-5)=8c-24.
```

Thus $`c=3`$. More generally, weights $`w_i\gt 0`$ give $`c=\sum_iw_iy_i/\sum_iw_i`$. The weighted residual condition is $`\sum_iw_i(y_i-c)=0`$. It need not make the unweighted residual sum zero when weights differ and the observations are not symmetric.

#### Problem LA-34: a projector that is not orthogonal

Constructed problem. Study

```math
P=\left(\begin{smallmatrix}1&2\\[0pt] 0&0\end{smallmatrix}\right)
```

.

Solution. Multiplication gives $`P^2=P`$. Its range is the horizontal axis and its kernel is $`\mathop{\mathrm{span}} ((-2,1))`$. These spaces are not perpendicular, and $`P^\top\ne P`$. For $`b=(0,1)`$, it returns $`(2,0)`$, whose distance from $`b`$ is $`\sqrt5`$. The closest horizontal-axis point is $`(0,0)`$ at distance one. Idempotence describes a projection along a complement; orthogonality is the condition that makes it a nearest-point map.

#### Problem LA-35: approximating a function by a straight line

Constructed problem. Find the best $`a+bx`$ approximation to $`x^2`$ in squared-integral error on $`[-1,1]`$.

Solution. Use the inner product $`\langle f,g\rangle=\int_{-1}^1f(x)g(x)\,dx`$. Orthogonality of the residual to $`1`$ and $`x`$ gives

```math
\int_{-1}^1(x^2-a-bx)\,dx=\frac23-2a=0,
```

```math
\int_{-1}^1x(x^2-a-bx)\,dx=-\frac23b=0.
```

Hence $`a=1/3`$, $`b=0`$. The minimum error is

```math
\int_{-1}^1(x^2-1/3)^2dx
=\frac25-\frac49+\frac29=\frac8{45}.
```

The odd integrals vanish by symmetry. Projection works in this polynomial subspace just as it does for a finite vector of observations.

#### Problem LA-36: fitting with duplicated columns

Constructed problem. Let $`X`$ have two identical columns, both equal to $`\mathbf1\in\mathbb R^n`$. Fit a response $`y`$.

Solution. The fitted vector is $`(b_1+b_2)\mathbf1`$. Minimizing over $`s=b_1+b_2`$ gives $`s=\bar y`$. Thus every $`(b_1,b_2)`$ with that sum minimizes squared error. Decompose any such pair as $`(\bar y/2,\bar y/2)+t(1,-1)`$. Its squared norm is $`\bar y^2/2+2t^2`$, uniquely minimized at $`t=0`$. The fitted vector is unique, the ordinary coefficient vector is not, and the minimum-norm coefficient vector is unique.

### 3.1 Theory

An eigenvector is a nonzero $`v`$ satisfying $`Av=\lambda v`$. A matrix is diagonalizable if it has a basis of eigenvectors: $`A=V\Lambda V^{-1}`$. A general matrix may fail to have enough eigenvectors. For example,

```math
\left(\begin{smallmatrix}1&1\\[0pt] 0&1\end{smallmatrix}\right)
```

 has only one independent eigenvector.

The real spectral theorem says that a symmetric real matrix has an orthonormal eigenbasis:

```math
A=Q\Lambda Q^\top,\qquad Q^\top Q=I.
```

The quadratic form is $`x^\top Ax`$. Symmetry can always be imposed for this purpose because $`x^\top Ax=x^\top[(A+A^\top)/2]x`$. A symmetric $`A`$ is positive semidefinite (PSD) when the form is nonnegative for every $`x`$, and positive definite (PD) when it is positive for every nonzero $`x`$.

For $`x\ne0`$, the Rayleigh quotient

```math
R_A(x)=\frac{x^\top Ax}{x^\top x}
```

lies between the smallest and largest eigenvalues. Eigenvectors for those endpoints attain the bounds.

### 3.2 Derivations and proofs

#### Working in eigenvector coordinates

The coordinate change $`y=Q^\top x`$ deserves a careful reading. It takes the same geometric vector $`x`$ and records its components along the orthonormal eigenvectors. The inverse change is $`x=Qy`$. Substituting into a symmetric quadratic gives

```math
x^\top Ax=(Qy)^\top(Q\Lambda Q^\top)(Qy)
=y^\top\Lambda y.
```

The middle factors cancel because $`Q^\top Q=I`$. No approximation has been made. In the new coordinates the cross terms disappear, so each eigenvalue determines the curvature along one independent direction.

A positive eigenvalue makes the form rise in that direction. A negative eigenvalue makes it fall. A zero eigenvalue makes it flat. A form with at least one positive and one negative eigenvalue is indefinite. For a symmetric PSD matrix, $`x^\top Ax=0`$ holds exactly when $`Ax=0`$, since the sum $`\sum_i\lambda_i y_i^2`$ can vanish only by setting every positive-eigenvalue coordinate to zero.

Repeated eigenvalues do not cause a failure of diagonalization for symmetric matrices. They make the choice of eigenvectors less unique. If the eigenspace for $`\lambda`$ has dimension two, any orthonormal basis of that plane is valid. The eigenspace is determined by the matrix even when the particular basis inside it is not.

#### Completing a block quadratic

For a symmetric block matrix

```math
M=\begin{pmatrix}B&C\\[0pt] C^\top&D\end{pmatrix},
```

assume $`B`$ is PD. Its quadratic form at $`(x,y)`$ can be rewritten as

```math
x^\top Bx+2x^\top Cy+y^\top Dy
=(x+B^{-1}Cy)^\top B(x+B^{-1}Cy)
+y^\top(D-C^\top B^{-1}C)y.
```

The expansion on the right produces the original cross term and cancels the added $`y^\top C^\top B^{-1}Cy`$. Therefore $`M`$ is PSD exactly when the Schur complement $`D-C^\top B^{-1}C`$ is PSD; it is PD exactly when that complement is PD. For necessity choose $`x=-B^{-1}Cy`$ to remove the first square. For sufficiency both displayed terms are nonnegative.

This derivation also describes minimization over $`x`$ with $`y`$ fixed. The best $`x`$ is $`-B^{-1}Cy`$, and the remaining cost is the Schur-complement quadratic in $`y`$.

#### Positive definite and invertible

For a symmetric PD matrix, $`Ax=0`$ would imply $`x^\top Ax=0`$, forcing $`x=0`$. Thus it is invertible. The converse is false: $`-I`$ and $`\mathop{\mathrm{diag}} (1,-1)`$ are invertible but not PD. For a symmetric PSD matrix the converse does hold, because invertibility removes every zero eigenvalue and PSD removes every negative one.

Its inverse is also PD: if $`A=Q\Lambda Q^\top`$ with all $`\lambda_i\gt 0`$, then $`A^{-1}=Q\Lambda^{-1}Q^\top`$ has positive eigenvalues $`1/\lambda_i`$. Its positive square root is $`A^{1/2}=Q\Lambda^{1/2}Q^\top`$. For PSD matrices the same square root exists, with zeros kept as zeros.

#### Spectral theorem in finite dimensions

For symmetric $`A`$, the continuous function $`x^\top Ax`$ has a maximum on the unit sphere, which is compact. Let $`q`$ be a maximizer. For every $`h\perp q`$, differentiate the form along $`(q+th)/\Vert q+th\Vert`$ at $`t=0`$. Its derivative is $`2h^\top Aq`$, so maximality gives $`h^\top Aq=0`$. Therefore $`Aq`$ is parallel to $`q`$ and $`Aq=\lambda q`$.

The orthogonal complement of $`q`$ is invariant: if $`x\perp q`$, then $`q^\top Ax=(Aq)^\top x=\lambda q^\top x=0`$. Restrict $`A`$ to that complement and repeat. Induction on dimension yields an orthonormal eigenbasis.

#### Positive semidefiniteness and eigenvalues

Put $`y=Q^\top x`$. Orthogonality preserves length and gives

```math
x^\top Ax=\sum_i\lambda_i y_i^2.
```

If every eigenvalue is nonnegative, the sum is nonnegative. Conversely, taking $`x=q_i`$ gives $`x^\top Ax=\lambda_i`$, so PSD forces every eigenvalue to be nonnegative. PD is equivalent to every eigenvalue being strictly positive. The Rayleigh bounds follow by dividing the sum by $`\sum_i y_i^2`$: the quotient is a weighted average of the eigenvalues.

#### Symmetry cannot be inferred from a quadratic expression

Every square real matrix decomposes as $`A=S+K`$, where $`S=(A+A^\top)/2`$ is symmetric and $`K=(A-A^\top)/2`$ is skew-symmetric. For any real vector,

```math
x^\top Kx=(x^\top Kx)^\top=x^\top K^\top x=-x^\top Kx,
```

so $`x^\top Kx=0`$. The quadratic form sees only $`S`$. In particular, a rotation matrix

```math
\begin{pmatrix}0&-1\\[0pt] 1&0\end{pmatrix}
```

 has quadratic expression zero for every $`x`$, but is not a symmetric positive-semidefinite operator. Whenever positivity is stated for an operator in these notes, symmetry is included.

Symmetry also ensures an orthonormal eigenbasis. A general matrix can lack any eigenbasis. For

```math
J=\begin{pmatrix}1&1\\[0pt] 0&1\end{pmatrix}
```

, the only eigenspace is the first coordinate axis, yet the matrix acts on a two-dimensional space. Its powers are

```math
J^n=\begin{pmatrix}1&n\\[0pt] 0&1\end{pmatrix}
```

. Knowing its eigenvalue is one does not describe the linear growth caused by the off-diagonal part.

#### A positive square root

If $`A=Q\mathop{\mathrm{diag}} (\lambda_i)Q^\top`$ is PSD, define

```math
A^{1/2}=Q\mathop{\mathrm{diag}} (\sqrt{\lambda_i})Q^\top.
```

This matrix is symmetric PSD and squares to $`A`$. If $`A`$ is PD, its square root is invertible and gives $`x^\top Ax=\Vert A^{1/2}x\Vert ^2`$. This identity turns a positive quadratic form into an ordinary squared length after a coordinate transformation.

The symmetric PSD square root is unique. Any such root $`B`$ commutes with $`A=B^2`$. Thus each eigenspace of $`A`$ is invariant under $`B`$. Diagonalize the symmetric restriction of $`B`$ on each one. Its eigenvalues must be nonnegative numbers whose squares equal the corresponding $`\lambda_i`$, forcing them to equal $`\sqrt{\lambda_i}`$.

A change of variables in a quadratic form gives $`P^\top AP`$, called a congruence, whereas an operator basis change gives $`P^{-1}AP`$, a similarity. Congruence preserves the signs available in the form when $`P`$ is invertible, but need not preserve eigenvalues. Taking $`A=I`$ and $`P=\mathop{\mathrm{diag}} (2,3)`$ gives form matrix $`\mathop{\mathrm{diag}} (4,9)`$.

#### Example LA-E20: whitening a covariance form

Let

```math
\Sigma=\begin{pmatrix}5&3\\[0pt] 3&5\end{pmatrix}
```

. Find a transformation making its covariance identity.

Solution. Along $`q_+=(1,1)^\top/\sqrt2`$ and $`q_-=(1,-1)^\top/\sqrt2`$, the variances are eight and two. Set $`C=\mathop{\mathrm{diag}} (1/\sqrt8,1/\sqrt2)Q^\top`$. Then $`C\Sigma C^\top=I`$. For a centered random vector $`X`$ with this covariance, $`CX`$ has uncorrelated unit-variance coordinates.

This does not imply independence without an additional distributional assumption. If $`X`$ is jointly Gaussian, the transformed vector is Gaussian and the identity covariance does imply independent coordinates. The linear algebra establishes covariance normalization; probability assumptions determine the stronger interpretation.

### 3.3 Worked examples

#### Example LA-E07: a symmetric matrix with coupled coordinates

Diagonalize

```math
A=\left(\begin{smallmatrix}2&1\\[0pt] 1&2\end{smallmatrix}\right)
```

.

Solution. The directions $`q_1=(1,1)/\sqrt2`$ and $`q_2=(1,-1)/\sqrt2`$ have eigenvalues $`3`$ and $`1`$. Thus $`Q=(q_1,q_2)`$ and $`A=Q\mathop{\mathrm{diag}} (3,1)Q^\top`$. For $`x=(2,0)`$ the eigenbasis coordinates are $`(\sqrt2,\sqrt2)`$, and $`x^\top Ax=3(2)+1(2)=8`$. Since both eigenvalues are positive, $`A`$ is PD.

#### Example LA-E08: an indefinite form

Classify

```math
B=\left(\begin{smallmatrix}1&2\\[0pt] 2&1\end{smallmatrix}\right)
```

.

Solution. Along $`(1,1)`$ the eigenvalue is $`3`$; along $`(1,-1)`$ it is $`-1`$. Hence $`x^\top Bx`$ can have either sign. Positive diagonal entries have not guaranteed PSD: the off-diagonal terms make some directions negative.

#### Example LA-E09: covariance rank

Let $`X`$ contain $`n`$ observations on $`p`$ variables, and let $`H=I-\mathbf1\mathbf1^\top/n`$ center each column. For $`n\gt 1`$, the sample covariance is $`S=X^\top HX/(n-1)`$.

Solution. Since $`H^\top=H^2=H`$,

```math
a^\top Sa=\Vert HXa\Vert ^2/(n-1)\ge0.
```

Also $`\mathop{\mathrm{rank}} S=\mathop{\mathrm{rank}} (HX)\le\min(p,n-1)`$. Therefore $`p\ge n`$ guarantees singularity. For $`p\lt n`$, singularity is still possible if a centered variable is an exact linear combination of others. Invertibility requires the centered columns to be independent.

### 3.4 Selected exercises

#### Problem LA-08: a tridiagonal positive operator

Source: [Axler, Exercises 7C, 2](https://linear.axler.net/LADR4e.pdf). The $`4\times4`$ matrix has diagonal $`2`$ and neighboring off-diagonals $`-1`$.

Solution. Expand the quadratic form:

```math
x^\top Ax=x_1^2+(x_1-x_2)^2+(x_2-x_3)^2+(x_3-x_4)^2+x_4^2.
```

Every term is nonnegative. If the sum is zero, $`x_1=x_4=0`$ and all neighboring differences vanish, forcing $`x=0`$. Thus the matrix is PD, so zero is not an eigenvalue and the matrix is invertible. For complex coordinates the same argument uses absolute squares and the conjugate transpose.

#### Problem LA-09: positive entries without positivity of the form

Source: [Axler, Exercises 7C, 4](https://linear.axler.net/LADR4e.pdf).

Problem. For every dimension $`n\gt 1`$, construct a symmetric matrix with all entries strictly positive whose quadratic form is not positive semidefinite.

Solution. For any $`n\gt 1`$, let $`J`$ be the all-ones matrix and put $`A=2J-I`$. Its diagonal entries are one and its off-diagonal entries are two, all positive. It is symmetric. Yet $`x=e_1-e_2`$ satisfies $`Jx=0`$, so $`x^\top Ax=-\Vert x\Vert ^2=-2`$. Thus $`A`$ is not PSD. Entrywise positivity and nonnegativity of a quadratic form are different properties.

#### Problem LA-10: shifting eigenvalues

Source: [Axler, Exercises 7C, 13(a,b)](https://linear.axler.net/LADR4e.pdf).

Problem. For a symmetric matrix $`A`$, characterize the scalars $`\alpha`$ for which $`A-\alpha I`$ or $`\alpha I-A`$ is positive semidefinite, using the extreme eigenvalues.

Solution. If $`Aq_i=\lambda_iq_i`$, then $`(A-\alpha I)q_i=(\lambda_i-\alpha)q_i`$. This matrix is PSD exactly when $`\alpha\le\lambda_{\min}(A)`$. Similarly, $`\alpha I-A`$ is PSD exactly when $`\alpha\ge\lambda_{\max}(A)`$. The endpoints allow zero eigenvalues. Strict inequalities make the respective matrix PD.

[Contents](README.md) · [Part 2](01-linear-algebra-part-2.md) · [Part 3](01-linear-algebra-part-3.md)
