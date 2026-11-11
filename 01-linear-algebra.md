# Linear Algebra

These notes use real finite-dimensional spaces unless a different field is stated. Vectors are columns, $A^\top$ is the transpose, $I$ is the identity, and $\|x\|_2=(x^\top x)^{1/2}$. For complex vectors, replace transposes by conjugate transposes and use the Hermitian inner product.

Sources:

- Gilbert Strang, [MIT 18.06SC, Fall 2011](https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/), especially Units II and III.
- Sheldon Axler, [Linear Algebra Done Right, fourth edition](https://linear.axler.net/LADR4e.pdf), Chapters 2, 3, 6, and 7. Exercise numbers below refer to this edition.

The explanations and solutions are written for these notes. Source problems are briefly paraphrased and linked; examples and additional problems are constructed here.

Contents:

1. [Vector spaces and linear maps](#1-vector-spaces-and-linear-maps)
2. [Orthogonality and projections](#2-orthogonality-and-projections)
3. [Eigenvalues and quadratic forms](#3-eigenvalues-and-quadratic-forms)
4. [SVD and low-rank approximation](#4-svd-and-low-rank-approximation)
5. [Conditioning and regularization](#5-conditioning-and-regularization)
6. [Matrix differentiation](#6-matrix-differentiation)
7. [Further worked problems](#7-further-worked-problems)

## 1. Vector spaces and linear maps

### 1.1 Theory

A vector space is a set closed under vector addition and scalar multiplication, with the usual associativity, distributivity, identity, and inverse rules. A subspace must contain zero and be closed under both operations. An affine plane such as $x_1+x_2+x_3=1$ is therefore not a subspace: it misses zero.

The span of $v_1,\ldots,v_k$ consists of all vectors $\sum_i a_iv_i$. The list is linearly independent when $\sum_i a_iv_i=0$ forces every $a_i=0$. A basis is an independent spanning list. Independence gives unique coordinates; spanning ensures that coordinates exist. All bases of a finite-dimensional space have the same length, called its dimension.

A map $T:V\to W$ is linear if $T(ax+by)=aT(x)+bT(y)$. Its kernel and range are

$$
\ker T=\{v:T(v)=0\},\qquad
\operatorname{range}T=\{T(v):v\in V\}.
$$

Both are subspaces. The rank is the dimension of the range; the nullity is the dimension of the kernel. A matrix represents a map after bases are chosen. Changing the bases changes the matrix while leaving the underlying map unchanged.

For $A\in\mathbb R^{m\times n}$, the column space is its range. Solving $Ax=b$ asks whether $b$ lies in that space. If $x_0$ is one solution, every solution has the form $x_0+z$ with $z\in\ker A$: subtract two solutions to see why.

### 1.2 Derivations and proofs

#### Reading a linear combination

The equation $b=c_1a_1+\cdots+c_na_n$ asks whether the available directions $a_i$ can produce the target $b$. The coefficients can be positive, negative, or zero. A zero coefficient means that direction is unused; a negative one reverses it. The coefficients are the unknowns, while the vectors are the given information.

For example, if $a_1=(1,1)$ and $a_2=(1,-1)$, the equation $(3,1)=c_1a_1+c_2a_2$ becomes $c_1+c_2=3$, $c_1-c_2=1$. Adding and subtracting gives $c_1=2$, $c_2=1$. The target is a vector in the original space; $(2,1)$ is its coordinate vector relative to this particular basis. Using another basis changes its coordinates.

A spanning list may contain an unnecessary direction. If $a_3=a_1+a_2$, every use of $a_3$ can be replaced by uses of $a_1$ and $a_2$. A basis removes all such redundancy. To test whether a list is a basis, check both parts of the definition. A list can be independent but too short to span, or span but be dependent because it is too long.

#### Testing subspaces and linearity

For a proposed subspace, the closure test is: take arbitrary $u,v$ in the set and arbitrary real numbers $a,b$, and prove $au+bv$ is still in the set. This proves addition and scalar closure together. Check zero as well, since the empty set would otherwise pass a statement quantified over its elements.

For a set defined by homogeneous linear equations, the test is immediate. If $Au=Av=0$, then $A(au+bv)=aAu+bAv=0$. With a nonzero right-hand side $c$, the same computation gives $(a+b)c$, which need not equal $c$.

For a proposed linear map, testing two or three numerical inputs cannot establish linearity. Write the expression at general $au+bv$ and simplify. To disprove linearity, one counterexample suffices. In particular every linear map satisfies $T(0)=0$, because $T(0)=T(0+0)=2T(0)$.

The map $T(x)=Ax+c$ is linear exactly when $c=0$. When $c\ne0$ it is affine. The map $T(x_1,x_2)=x_1x_2$ sends zero to zero but is still nonlinear: $T(e_1+e_2)=1$, whereas $T(e_1)+T(e_2)=0$.

#### Extending an independent list to a basis

Start with an independent list in a finite-dimensional space. If it does not span, choose a vector outside its span and append it. The longer list remains independent: a dependence with a nonzero coefficient on the new vector would express that vector in the old span. If its new coefficient were zero, the original independence would force all remaining coefficients to vanish.

Continue until the list spans. The process stops because a finite-dimensional space cannot contain arbitrarily long independent lists. This construction explains the basis-extension step in rank-nullity. It also explains why the theorem applies to arbitrary linear maps, rather than only to matrices already in row-echelon form.

#### Four spaces belonging to one matrix

For $A\in\mathbb R^{m\times n}$, keep track of the ambient dimensions:

| Space | Ambient space | Dimension |
| --- | --- | --- |
| Column space $\operatorname{range}A$ | $\mathbb R^m$ | $r$ |
| Left null space $\ker A^\top$ | $\mathbb R^m$ | $m-r$ |
| Row space $\operatorname{range}A^\top$ | $\mathbb R^n$ | $r$ |
| Null space $\ker A$ | $\mathbb R^n$ | $n-r$ |

Column space and left null space are orthogonal complements. Row space and null space are orthogonal complements. A vector in $\mathbb R^n$ cannot generally be compared with a column-space vector in $\mathbb R^m$ when $m\ne n$; they live in different spaces.

The equality of row rank and column rank follows, for example, from elimination: multiplying by an invertible row-operation matrix preserves dependencies among columns, and a row-echelon matrix has one independent row and one pivot column for each pivot. Rank-nullity then supplies both null-space dimensions.

#### Rank-nullity

Let $T:V\to W$ be linear and $\dim V=n$. Choose a basis $u_1,\ldots,u_k$ of $\ker T$ and extend it to a basis

$$
u_1,\ldots,u_k,v_1,\ldots,v_{n-k}
$$

of $V$. Every vector in $V$ has an expansion in this basis. Applying $T$ removes the $u_i$ terms, so $Tv_1,\ldots,Tv_{n-k}$ span the range.

To check independence, suppose $\sum_j a_jTv_j=0$. Then $\sum_j a_jv_j\in\ker T$, so it is also a linear combination of the $u_i$. Independence of the extended basis forces all $a_j=0$. The images are consequently a basis of the range. Thus

$$
\dim V=\dim\ker T+\dim\operatorname{range}T.
$$

The proof also constructs a complement $U=\operatorname{span}(v_1,\ldots,v_{n-k})$ on which $T$ is one-to-one and maps onto its full range.

#### Row space and null space

The equation $Ax=0$ says that $x$ is perpendicular to every row of $A$. Hence $\ker A=\operatorname{range}(A^\top)^\perp$. In finite dimensions, taking orthogonal complements again gives

$$
(\ker A)^\perp=\operatorname{range}(A^\top).
$$

This identifies the part of a parameter vector that the map can detect. Adding a null-space vector changes the input without changing the output.

#### Coordinates and a change of basis

A coordinate vector records coefficients in a chosen basis. It is not the underlying geometric vector independently of that basis. Let $P$ have the new domain basis vectors as columns in old coordinates, and let $Q$ do the same for the codomain. A vector with new domain coordinates $z$ has old coordinates $Pz$. Applying the old matrix $A$ gives $APz$, and converting the output to new coordinates gives $Q^{-1}APz$. The new matrix is therefore

$$
A_{\mathrm{new}}=Q^{-1}AP.
$$

Both change-of-basis matrices are invertible because their columns are bases. Left and right multiplication by invertible maps preserve rank, so rank is independent of coordinates. The kernel's coordinate description changes, but its dimension does not.

For an operator from a space to itself with the same new basis on input and output, this becomes the similarity transformation $P^{-1}AP$. Similar matrices have the same eigenvalues: if $Av=\lambda v$, then $P^{-1}v$ is an eigenvector of the new matrix. A general two-sided coordinate transformation with unrelated $P,Q$ need not preserve eigenvalues.

There is another way to read the kernel. Two inputs $x,y$ give the same output exactly when $x-y\in\ker T$. Thus each attainable output corresponds to an entire translated kernel, $x+\ker T$. Rank-nullity counts the input directions that are distinguishable after applying the map and those that are lost. This interpretation also explains why adding null-space vectors cannot improve or worsen a least-squares fit.

#### Example LA-E19: the same operator in two bases

Let $A=\operatorname{diag}(2,3)$ in the standard basis of $\mathbb R^2$. Use the new basis $p_1=(1,1)^\top$, $p_2=(1,-1)^\top$.

Solution. The basis matrix is $P=\begin{pmatrix}1&1\\1&-1\end{pmatrix}$, with inverse $P/2$. Multiplication gives

$$
P^{-1}AP=\frac12
\begin{pmatrix}5&-1\\-1&5\end{pmatrix}.
$$

The first column says that the old vector $Ap_1=(2,3)^\top$ has new coordinates $(5/2,-1/2)^\top$. Reconstructing it as $(5/2)p_1-(1/2)p_2$ checks the interpretation. The new matrix has eigenvalues two and three, as expected, even though it is no longer diagonal. A diagonal appearance depends on coordinates; the operator's eigenvalues do not.

### 1.3 Worked examples

#### Example LA-E01: a plane and its coordinates

Find a basis for $S=\{(x,y,z):x+y+z=0\}$.

Solution. Solve for $z=-x-y$:

$$
(x,y,z)=x(1,0,-1)+y(0,1,-1).
$$

The two displayed vectors span $S$. If their linear combination is zero, its first two coordinates force both coefficients to vanish. They are a basis and $\dim S=2$. The equation defining the plane has removed one free coordinate from $\mathbb R^3$.

#### Example LA-E02: a kernel with a free variable

Let

$$
A=\begin{pmatrix}1&2&3\\0&1&1\end{pmatrix}.
$$

Solution. The second row of $Ax=0$ gives $x_2=-x_3$; the first then gives $x_1=-x_3$. Thus $\ker A=\operatorname{span}((-1,-1,1))$. The first two columns are independent, so $\operatorname{rank}A=2$. Rank-nullity reads $3=1+2$. For $b=(4,1)^\top$, setting $x_3=t$ gives all solutions $(2-t,1-t,t)^\top$.

#### Example LA-E03: differentiation as a linear map

Let $D:\mathcal P_3\to\mathcal P_2$ send $p$ to $p'$, where $\mathcal P_k$ means polynomials of degree at most $k$.

Solution. Differentiation respects sums and scalar multiples. Its kernel consists of constant polynomials, so its nullity is one. For any $a+bx+cx^2$, the polynomial $ax+(b/2)x^2+(c/3)x^3$ maps to it, so the range is all of $\mathcal P_2$ and the rank is three. Again $4=1+3$.

### 1.4 Selected exercises

#### Problem LA-01: a redundant spanning list

Source: [Axler, Exercises 2A, 1](https://linear.axler.net/LADR4e.pdf). Find four distinct vectors spanning the plane from Example LA-E01.

Solution. Use

$$
v_1=(1,0,-1),\quad v_2=(0,1,-1),\quad
v_3=(1,1,-2),\quad v_4=(1,-1,0).
$$

Each coordinate sum is zero, so their span lies in the plane. Conversely the first two already span it, as shown in the example. The other two satisfy $v_3=v_1+v_2$ and $v_4=v_1-v_2$. The list has four members but rank two; a spanning list need not be a basis.

#### Problem LA-02: cumulative sums preserve span

Source: [Axler, Exercises 2A, 3](https://linear.axler.net/LADR4e.pdf). Compare the span of $v_1,\ldots,v_m$ with that of $w_k=\sum_{i=1}^k v_i$.

Solution. Every $w_k$ is a combination of the $v_i$, giving one containment. For the other, $v_1=w_1$ and $v_k=w_k-w_{k-1}$ for $k\ge2$. Each original vector is therefore in the new span. Both spans agree. The cumulative-sum transformation is triangular with diagonal entries one, and the differences explicitly give its inverse.

#### Problem LA-03: surjectivity from a specified kernel

Source: [Axler, Exercises 3B, 12](https://linear.axler.net/LADR4e.pdf). A map $T:\mathbb F^4\to\mathbb F^2$ has kernel determined by $x_1=5x_2$ and $x_3=7x_4$. Prove surjectivity.

Solution. Write every kernel vector as

$$
s(5,1,0,0)+t(0,0,7,1).
$$

These two vectors are independent, so the nullity is two. Rank-nullity gives rank $4-2=2$. The range is a two-dimensional subspace of a two-dimensional codomain, and therefore equals that codomain. Merely knowing that two equations define the kernel would be insufficient if those equations were redundant; the explicit basis checks their independence.

#### Problem LA-04: decomposition for an idempotent

Source: [Axler, Exercises 3B, 27](https://linear.axler.net/LADR4e.pdf). For $P^2=P$, prove $V=\ker P\oplus\operatorname{range}P$.

Solution. For every $v$,

$$
v=(v-Pv)+Pv,\qquad P(v-Pv)=Pv-P^2v=0.
$$

The first term is in the kernel and the second is in the range. If $z$ belongs to both, write $z=Py$ and use $Pz=P^2y=Py=z$. Kernel membership also says $Pz=0$, so $z=0$. The sum is direct. Idempotence alone does not make this an orthogonal decomposition; symmetry supplies that additional property.

## 2. Orthogonality and projections

### Further problems on spaces and maps

#### Problem LA-21: when changing one vector preserves independence

Source: [Axler, Exercises 2A, 9](https://linear.axler.net/LADR4e.pdf).

Given independent $v_1,\ldots,v_m$, replace $v_1$ by $5v_1-4v_2$, assuming $m\ge2$.

Solution. Suppose $a_1(5v_1-4v_2)+a_2v_2+\cdots+a_mv_m=0$. Collect coefficients of the original vectors:

$$
5a_1v_1+(a_2-4a_1)v_2+\sum_{j=3}^m a_jv_j=0.
$$

Original independence gives $5a_1=0$, hence $a_1=0$, and then $a_2=0$ and all later coefficients zero. The replacement preserves independence over $\mathbb R$ or $\mathbb C$. The nonzero coefficient of the replaced vector makes the change reversible.

#### Problem LA-22: sums of independent lists

Source: [Axler, Exercises 2A, 11](https://linear.axler.net/LADR4e.pdf).

Problem. Determine whether adding corresponding vectors of two linearly independent lists must produce another independent list. Prove the claim or give a counterexample.

Solution. The assertion that coordinatewise sums of two independent lists must remain independent is false. Choose any independent $v_1,\ldots,v_m$ and let $w_i=-v_i$. The $w_i$ are independent, since a dependence among them is the negative of a dependence among the $v_i$. Yet each sum $v_i+w_i=0$, so the new list is dependent. Independence concerns relationships within one list, and gives no restriction on cancellation between two different lists.

#### Problem LA-23: preserving independent vectors

Source: [Axler, Exercises 3B, 9](https://linear.axler.net/LADR4e.pdf).

Problem. Prove that an injective linear map takes every independent list to an independent list, and examine the converse.

Solution. Let $T$ be injective and suppose $\sum_i a_iTv_i=0$. By linearity this is $T(\sum_i a_iv_i)=0$. Injectivity implies $\sum_i a_iv_i=0$, and independence of the $v_i$ implies $a_i=0$ for every $i$. Thus injective linear maps preserve independence. The converse also holds if a map preserves every independent one-vector list: a nonzero kernel vector would be sent to zero, violating that property.

#### Problem LA-24: rank of a composition

Source: [Axler, Exercises 3B, 23](https://linear.axler.net/LADR4e.pdf).

Problem. For compatible finite-dimensional linear maps $T$ and $S$, prove that the rank of $ST$ is at most the rank of either factor.

Solution. For $T:U\to V$ and $S:V\to W$, $\operatorname{range}(ST)=S(\operatorname{range}T)$. Its dimension cannot exceed $\dim\operatorname{range}T$, since images of a basis of that space span the image. It is also a subspace of $\operatorname{range}S$, so its dimension cannot exceed $\operatorname{rank}S$. Therefore

$$
\operatorname{rank}(ST)\le\min(\operatorname{rank}S,\operatorname{rank}T).
$$

Equality need not hold: $T$ can map entirely into $\ker S$, in which case the composition is zero even though neither map is zero.

#### Problem LA-25: a parameter changes rank

Constructed problem. For $A_t=\left(\begin{smallmatrix}1&1\\1&t\end{smallmatrix}\right)$, determine rank, kernel, and solvability for every $t$.

Solution. Subtract the first row from the second. The reduced second row is $(0,t-1)$. If $t\ne1$, both coordinates are pivots, the rank is two, the kernel is zero, and every right-hand side has a unique solution.

If $t=1$, the rank is one and the kernel is $\operatorname{span}((1,-1))$. The equation $A_1x=b$ is solvable exactly when $b_1=b_2$. If that equality holds, all solutions are $(b_1-s,s)$ for $s\in\mathbb R$. Thus one parameter value changes both the number of admissible targets and the uniqueness of coefficients.

#### Problem LA-26: the matrix after a change of basis

Constructed problem. Let $T(x,y)=(2x+y,x+2y)$ and use the basis $b_1=(1,1)$, $b_2=(1,-1)$.

Solution. The basis matrix is $B=\left(\begin{smallmatrix}1&1\\1&-1\end{smallmatrix}\right)$ and $B^{-1}=\tfrac12B$. For an input coordinate vector $c$, the original vector is $Bc$, its image is $ABc$, and the coordinates of the image are $B^{-1}ABc$. Hence the new matrix is $B^{-1}AB$.

Here $Tb_1=3b_1$ and $Tb_2=b_2$, so the new matrix is $\operatorname{diag}(3,1)$. The map has not changed; the new coordinates align with directions that it scales without mixing.

#### Problem LA-27: an inconsistent system

Constructed problem. Determine all least-squares targets that are exactly solvable for $A=\left(\begin{smallmatrix}1&0\\0&1\\1&1\end{smallmatrix}\right)$.

Solution. The range consists of vectors $(s,t,s+t)$, so $Ax=b$ is exactly solvable if and only if $b_3=b_1+b_2$. A left-null vector is $z=(-1,-1,1)$, and the same condition is $z^\top b=0$. If the condition fails, the component of $b$ along $z$ is unreachable. Section 2 shows that least squares removes that component by orthogonal projection.

#### Problem LA-28: images and preimages

Constructed problem. Let $T(x,y,z)=(x+y,y+z)$ and find the preimage of the line $\{(s,s):s\in\mathbb R\}$.

Solution. Membership requires $x+y=y+z$, hence $x=z$. The preimage is $\{(x,y,x)\}$, with basis $(1,0,1)$ and $(0,1,0)$ and dimension two. The kernel of $T$ is generated by $(-1,1,-1)$, and its dimension is one. The target line contributes one additional independent image direction, giving $2=1+1$. A preimage under a linear map of a subspace is itself a subspace.

### 2.1 Theory

An inner product measures angles and lengths. In $\mathbb R^n$ it is $x^\top y$; orthogonality means $x^\top y=0$. The orthogonal complement $S^\perp$ contains the vectors perpendicular to every element of $S$.

An orthogonal projection $P_Sb$ is the element of $S$ closest to $b$. It is characterized by

$$
P_Sb\in S,\qquad b-P_Sb\in S^\perp.
$$

If the columns of $Q$ are an orthonormal basis of $S$, then $P_S=QQ^\top$. If instead $S$ is represented by independent columns of $A$, then $P_S=A(A^\top A)^{-1}A^\top$. Neither formula requires $S$ to be all of the ambient space.

Gram-Schmidt constructs orthonormal directions from independent vectors:

$$
u_k=a_k-\sum_{j<k}(q_j^\top a_k)q_j,\qquad q_k=u_k/\|u_k\|.
$$

Independence guarantees $u_k\ne0$. In numerical work, Householder QR or reorthogonalization is often more reliable than the unmodified procedure; see [Numerical Methods](07-numerical-methods.md#2-numerical-linear-algebra).

### 2.2 Derivations and proofs

#### Finding the closest point one coefficient at a time

For a line generated by nonzero $a$, every candidate is $ca$. Its squared error from $b$ is

$$
g(c)=\|b-ca\|^2
=b^\top b-2c\,a^\top b+c^2a^\top a.
$$

This is an upward-opening scalar quadratic because $a^\top a>0$. Setting $g'(c)=0$ gives $c=a^\top b/(a^\top a)$. Substituting that coefficient into $a^\top(b-ca)$ gives zero. Thus the perpendicular residual arises from minimization; it is also a way to solve the minimization.

For several independent columns, the coefficients interact. If the columns are orthonormal, each coefficient is $q_i^\top b$ and can be found separately. If the columns are not orthogonal, the interaction is recorded by the Gram matrix $A^\top A$. Its $(i,j)$ entry is the inner product of columns $i$ and $j$.

The normal equations correct for that interaction. Writing $\widehat x_i=a_i^\top b$ for nonorthogonal columns counts overlapping directions more than once. Multiplication by $(A^\top A)^{-1}$ removes that overlap.

#### Cauchy-Schwarz and its equality case

For $v\ne0$, consider the nonnegative squared norm of $u-tv$ and choose $t=\langle u,v\rangle/\|v\|^2$ in a real space. Expansion yields

$$
0\le\left\|u-\frac{\langle u,v\rangle}{\|v\|^2}v\right\|^2
=\|u\|^2-\frac{|\langle u,v\rangle|^2}{\|v\|^2}.
$$

Therefore $|\langle u,v\rangle|\le\|u\|\|v\|$. Equality holds exactly when the residual vector is zero, which means $u$ is a multiple of $v$. If $v=0$, both sides are zero. This proof is itself a projection calculation.

#### Existence in rank-deficient least squares

Dependent columns do not prevent a best fitted vector from existing. Their span is a finite-dimensional subspace and therefore closed. Choose an orthonormal basis $Q$ of that span and project $b$ to $QQ^\top b$. Because this fitted vector belongs to the column space, at least one coefficient vector produces it.

Two minimizing coefficient vectors differ by an element of $\ker A$. To choose one uniquely, require it to lie in $(\ker A)^\perp$. Every vector has a unique decomposition into these two orthogonal spaces; removing the kernel component leaves the fit unchanged and reduces the norm. This proves the minimum-norm characterization before any formula for a pseudoinverse is used.

#### Weighted least squares

Suppose $W$ is symmetric PD and squared error is measured by $r^\top Wr$ instead of $r^\top r$. An eigendecomposition gives a square root $W^{1/2}$, so

$$
(b-Ax)^\top W(b-Ax)=\|W^{1/2}b-W^{1/2}Ax\|^2.
$$

Applying ordinary least squares to the transformed system gives $A^\top WA\widehat x=A^\top Wb$. With independent columns,

$$
\widehat x=(A^\top WA)^{-1}A^\top Wb,\qquad
P_W=A(A^\top WA)^{-1}A^\top W.
$$

The matrix $P_W$ is idempotent but generally not symmetric in Euclidean coordinates. Instead $P_W^\top W=WP_W$. Its residual is orthogonal to the fitted space under the weighted inner product $\langle u,v\rangle_W=u^\top Wv$. Always specify the geometry in which a projection is orthogonal.

#### Projection and least squares

Minimize $\|b-Ax\|^2$ for $A\in\mathbb R^{m\times n}$. A minimizing residual must be perpendicular to each column: otherwise moving along that column decreases the residual length. Algebraically, differentiation gives

$$
A^\top(A\widehat x-b)=0.
$$

These are the normal equations. If $A$ has independent columns, $z^\top A^\top Az=\|Az\|^2>0$ for $z\ne0$, so

$$
\widehat x=(A^\top A)^{-1}A^\top b,\qquad
P=A(A^\top A)^{-1}A^\top.
$$

Transposing $P$ gives $P^\top=P$, because $A^\top A$ and its inverse are symmetric. Multiplying gives

$$
P^2=A(A^\top A)^{-1}(A^\top A)(A^\top A)^{-1}A^\top=P.
$$

For any $x$, the orthogonal decomposition

$$
\|b-Ax\|^2=\|b-Pb\|^2+\|Pb-Ax\|^2
$$

proves global optimality. It also proves that the fitted vector is unique even when the coefficients are not. In the rank-deficient case use $P=AA^+$, where $A^+$ is the pseudoinverse introduced below.

#### Gram–Schmidt with every subtraction visible

Take independent vectors $u_1=(1,1,0)^\top$, $u_2=(1,0,1)^\top$, and $u_3=(0,1,1)^\top$. Normalize the first:

$$
q_1=\frac1{\sqrt2}(1,1,0)^\top.
$$

Remove from $u_2$ the part already in the first direction:

$$
v_2=u_2-(q_1^\top u_2)q_1
=(1,0,1)^\top-\tfrac12(1,1,0)^\top
=\tfrac12(1,-1,2)^\top.
$$

Its norm is $\sqrt6/2$, giving $q_2=(1,-1,2)^\top/\sqrt6$. For the third vector, the coefficients are $q_1^\top u_3=1/\sqrt2$ and $q_2^\top u_3=1/\sqrt6$. Therefore

$$
v_3=u_3-\tfrac12(1,1,0)^\top-\tfrac16(1,-1,2)^\top
=\tfrac23(-1,1,1)^\top.
$$

Normalizing gives $q_3=(-1,1,1)^\top/\sqrt3$. Direct inner products check that the three $q_i$ are orthonormal. At each stage, the subtraction preserves the span built so far, while removing overlap with earlier directions.

Collecting $Q=(q_1,q_2,q_3)$ gives $U=QR$, with

$$
R=\begin{pmatrix}
\sqrt2&1/\sqrt2&1/\sqrt2\\
0&\sqrt6/2&1/\sqrt6\\
0&0&2/\sqrt3
\end{pmatrix}.
$$

The entries above the diagonal are projection coefficients; the diagonal entries are residual norms. A zero residual would reveal dependence, so normalization at that step would be invalid.

#### Two successive projections

If orthogonal projection matrices $P,Q$ commute, their product is symmetric and idempotent:
$(PQ)^\top=QP=PQ$ and $(PQ)^2=P^2Q^2=PQ$. Its range is the intersection of the two ranges. Conversely, if their product is an orthogonal projection, symmetry forces $PQ=QP$.

For projection onto the first coordinate axis and projection onto the line $(1,1)$,

$$
P=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
Q=\tfrac12\begin{pmatrix}1&1\\1&1\end{pmatrix},\qquad
PQ=\tfrac12\begin{pmatrix}1&1\\0&0\end{pmatrix}.
$$

The product is neither symmetric nor idempotent. Projecting onto one set and then another need not project onto their intersection in a single pass.

For a weighted inner product with $W=C^\top C$ invertible, the map $z\mapsto Cz$ converts weighted lengths into Euclidean lengths. Weighted least squares is ordinary least squares for $CX$ and $Cy$. Its projection in original coordinates is $C^{-1}P_{CX}C$, explaining why it is idempotent without generally being Euclidean-symmetric.

### 2.3 Worked examples

#### Example LA-E04: projection onto a line

Project $b=(2,0,1)^\top$ onto the line generated by $a=(1,1,1)^\top$.

Solution. The coefficient is $(a^\top b)/(a^\top a)=3/3=1$, so $Pb=(1,1,1)^\top$ and $r=(1,-1,0)^\top$. Directly $a^\top r=0$. The projection matrix is $aa^\top/3$, with every entry $1/3$. The squared distance is $\|r\|^2=2$.

#### Example LA-E05: fitting a line to three observations

Fit $y_i=a+bt_i$ to $(t,y)=(0,1),(1,2),(2,2)$.

Solution. Here

$$
X=\begin{pmatrix}1&0\\1&1\\1&2\end{pmatrix},\quad
X^\top X=\begin{pmatrix}3&3\\3&5\end{pmatrix},\quad
X^\top y=\binom56.
$$

Subtracting the first normal equation from the second gives $2b=1$, so $b=1/2$ and $a=7/6$. The fitted values are $(7/6,5/3,13/6)$ and the residuals are $(-1/6,1/3,-1/6)$. Their sum and their inner product with $(0,1,2)$ are both zero.

#### Example LA-E06: Gram-Schmidt in a plane

Orthonormalize $a_1=(1,1,0)$ and $a_2=(1,0,1)$.

Solution. First $q_1=(1,1,0)/\sqrt2$. Subtract the component $(q_1^\top a_2)q_1=(1/2,1/2,0)$:

$$
u_2=(1/2,-1/2,1),\qquad q_2=(1,-1,2)/\sqrt6.
$$

Both vectors have unit length and their inner product is zero. The projection onto their span is $q_1q_1^\top+q_2q_2^\top$. A normal direction is $(1,-1,-1)$, so the same projection is $I-nn^\top/3$.

### 2.4 Selected exercises

#### Problem LA-05: projecting onto the orthogonal complement

Source: [Axler, Exercises 6C, 5](https://linear.axler.net/LADR4e.pdf).

Problem. For a subspace $U$ of a finite-dimensional inner-product space, express the orthogonal projection onto $U^\perp$ in terms of the projection onto $U$ and prove the identity.

Solution. Write $v=u+w$ with $u\in U$ and $w\in U^\perp$. By the projection characterization, $P_Uv=u$ and $P_{U^\perp}v=w$. Hence $P_{U^\perp}v=v-P_Uv$ for every $v$, and $P_{U^\perp}=I-P_U$. This also gives $(I-P_U)^2=I-P_U$ and $(I-P_U)^\top=I-P_U$.

#### Problem LA-06: recognizing an orthogonal projection

Source: [Axler, Exercises 6C, 9](https://linear.axler.net/LADR4e.pdf). Assume $P^2=P$ and $\ker P\perp\operatorname{range}P$.

Solution. Problem LA-04 gives $v=(v-Pv)+Pv$, with the terms in the kernel and range respectively. The assumed orthogonality makes this the orthogonal decomposition relative to $U=\operatorname{range}P$. Therefore $Pv=P_Uv$. For a direct minimization check, any $u\in U$ satisfies

$$
\|v-u\|^2=\|v-Pv\|^2+\|Pv-u\|^2.
$$

The second term is nonnegative and vanishes only at $u=Pv$. Thus $P$ selects the unique closest point.

#### Problem LA-07: a closest vector in four dimensions

Source: [Axler, Exercises 6C, 15](https://linear.axler.net/LADR4e.pdf). Minimize distance to $(1,2,3,4)$ over $U=\operatorname{span}((1,1,0,0),(1,1,1,2))$.

Solution. Replace the generators by the orthogonal vectors $a=(1,1,0,0)$ and $c=(0,0,1,2)$; their span is unchanged. Then

$$
P_Ub=\frac{b^\top a}{a^\top a}a+
\frac{b^\top c}{c^\top c}c
=\frac32a+\frac{11}{5}c
=\left(\frac32,\frac32,\frac{11}{5},\frac{22}{5}\right).
$$

The residual is $(-1/2,1/2,4/5,-2/5)$. Its inner products with $a$ and $c$ both vanish. The minimum squared distance is $1/2+16/25+4/25=13/10$. The Pythagorean identity above proves uniqueness.

## 3. Eigenvalues and quadratic forms

### Further problems on orthogonality and least squares

#### Problem LA-29: an inner product produced by a map

Source: [Axler, Exercises 6A, 2](https://linear.axler.net/LADR4e.pdf).

Problem. Let $S$ be a linear map between inner-product spaces. Determine exactly when $\langle u,v\rangle_1=\langle Su,Sv\rangle$ defines an inner product on its domain.

Solution. Define $\langle u,v\rangle_1=\langle Su,Sv\rangle$. Symmetry and linearity follow from the original inner product and linearity of $S$. Nonnegativity follows from $\langle v,v\rangle_1=\|Sv\|^2\ge0$. This equals zero exactly when $Sv=0$. Thus positive definiteness holds exactly when $\ker S=\{0\}$, that is, when $S$ is injective. For a matrix $S$, the induced Gram matrix is $S^\top S$, recovering the earlier full-column-rank criterion.

#### Problem LA-30: an orthogonal decomposition in two coordinates

Source: [Axler, Exercises 6A, 11](https://linear.axler.net/LADR4e.pdf).

Problem. Write $(1,2)$ as the sum of a vector parallel to $(1,3)$ and a vector perpendicular to $(1,3)$.

Solution. Decompose $b=(1,2)$ into a multiple of $a=(1,3)$ and a perpendicular residual. The projection coefficient is $a^\top b/(a^\top a)=7/10$. Therefore

$$
u=(7/10,21/10),\qquad v=b-u=(3/10,-1/10).
$$

Their sum is $b$ and $a^\top v=3/10-3/10=0$. Also $\|u\|^2=49/10$ and $\|v\|^2=1/10$, whose sum equals $\|b\|^2=5$. This checks both the arithmetic and the orthogonal decomposition.

#### Problem LA-31: Bessel's inequality and equality

Source: [Axler, Exercises 6B, 3](https://linear.axler.net/LADR4e.pdf).

Problem. For an orthonormal list $e_1,\ldots,e_m$, prove Bessel's inequality and characterize equality for a vector $v$.

Solution. For orthonormal $e_1,\ldots,e_m$, put $p=\sum_i\langle v,e_i\rangle e_i$ and $r=v-p$. Taking the inner product with any $e_j$ gives $\langle r,e_j\rangle=0$, so

$$
\|v\|^2=\sum_i|\langle v,e_i\rangle|^2+\|r\|^2.
$$

The sum of coefficient squares is at most $\|v\|^2$. Equality holds exactly when $r=0$, or $v\in\operatorname{span}(e_1,\ldots,e_m)$. No assumption that the orthonormal list spans the whole ambient space is needed.

#### Problem LA-32: orthogonality between two subspaces

Source: [Axler, Exercises 6C, 7](https://linear.axler.net/LADR4e.pdf).

Problem. For subspaces $X,Y$, prove that every vector of $X$ is orthogonal to every vector of $Y$ exactly when $P_XP_Y=0$.

Solution. If every $x\in X$ is perpendicular to every $y\in Y$, then $P_Yv\in Y\subseteq X^\perp$, so $P_XP_Yv=0$ for every $v$.

Conversely suppose $P_XP_Y=0$. For $y\in Y$, $P_Yy=y$, hence $P_Xy=0$. The kernel of $P_X$ is $X^\perp$, so $y\perp X$. This holds for every $y$, establishing the claimed equivalence. It also implies $P_YP_X=0$ by transposition.

#### Problem LA-33: a weighted average

Constructed problem. Fit a common mean $c$ to observations $1,3,5$, with positive weights $1,2,1$.

Solution. Minimize $(1-c)^2+2(3-c)^2+(5-c)^2$. Differentiate:

$$
2(c-1)+4(c-3)+2(c-5)=8c-24.
$$

Thus $c=3$. More generally, weights $w_i>0$ give $c=\sum_iw_iy_i/\sum_iw_i$. The weighted residual condition is $\sum_iw_i(y_i-c)=0$. It need not make the unweighted residual sum zero when weights differ and the observations are not symmetric.

#### Problem LA-34: a projector that is not orthogonal

Constructed problem. Study $P=\left(\begin{smallmatrix}1&2\\0&0\end{smallmatrix}\right)$.

Solution. Multiplication gives $P^2=P$. Its range is the horizontal axis and its kernel is $\operatorname{span}((-2,1))$. These spaces are not perpendicular, and $P^\top\ne P$. For $b=(0,1)$, it returns $(2,0)$, whose distance from $b$ is $\sqrt5$. The closest horizontal-axis point is $(0,0)$ at distance one. Idempotence describes a projection along a complement; orthogonality is the condition that makes it a nearest-point map.

#### Problem LA-35: approximating a function by a straight line

Constructed problem. Find the best $a+bx$ approximation to $x^2$ in squared-integral error on $[-1,1]$.

Solution. Use the inner product $\langle f,g\rangle=\int_{-1}^1f(x)g(x)\,dx$. Orthogonality of the residual to $1$ and $x$ gives

$$
\int_{-1}^1(x^2-a-bx)\,dx=\frac23-2a=0,
$$

$$
\int_{-1}^1x(x^2-a-bx)\,dx=-\frac23b=0.
$$

Hence $a=1/3$, $b=0$. The minimum error is

$$
\int_{-1}^1(x^2-1/3)^2dx
=\frac25-\frac49+\frac29=\frac8{45}.
$$

The odd integrals vanish by symmetry. Projection works in this polynomial subspace just as it does for a finite vector of observations.

#### Problem LA-36: fitting with duplicated columns

Constructed problem. Let $X$ have two identical columns, both equal to $\mathbf1\in\mathbb R^n$. Fit a response $y$.

Solution. The fitted vector is $(b_1+b_2)\mathbf1$. Minimizing over $s=b_1+b_2$ gives $s=\bar y$. Thus every $(b_1,b_2)$ with that sum minimizes squared error. Decompose any such pair as $(\bar y/2,\bar y/2)+t(1,-1)$. Its squared norm is $\bar y^2/2+2t^2$, uniquely minimized at $t=0$. The fitted vector is unique, the ordinary coefficient vector is not, and the minimum-norm coefficient vector is unique.

### 3.1 Theory

An eigenvector is a nonzero $v$ satisfying $Av=\lambda v$. A matrix is diagonalizable if it has a basis of eigenvectors: $A=V\Lambda V^{-1}$. A general matrix may fail to have enough eigenvectors. For example, $\left(\begin{smallmatrix}1&1\\0&1\end{smallmatrix}\right)$ has only one independent eigenvector.

The real spectral theorem says that a symmetric real matrix has an orthonormal eigenbasis:

$$
A=Q\Lambda Q^\top,\qquad Q^\top Q=I.
$$

The quadratic form is $x^\top Ax$. Symmetry can always be imposed for this purpose because $x^\top Ax=x^\top[(A+A^\top)/2]x$. A symmetric $A$ is positive semidefinite (PSD) when the form is nonnegative for every $x$, and positive definite (PD) when it is positive for every nonzero $x$.

For $x\ne0$, the Rayleigh quotient

$$
R_A(x)=\frac{x^\top Ax}{x^\top x}
$$

lies between the smallest and largest eigenvalues. Eigenvectors for those endpoints attain the bounds.

### 3.2 Derivations and proofs

#### Working in eigenvector coordinates

The coordinate change $y=Q^\top x$ deserves a careful reading. It takes the same geometric vector $x$ and records its components along the orthonormal eigenvectors. The inverse change is $x=Qy$. Substituting into a symmetric quadratic gives

$$
x^\top Ax=(Qy)^\top(Q\Lambda Q^\top)(Qy)
=y^\top\Lambda y.
$$

The middle factors cancel because $Q^\top Q=I$. No approximation has been made. In the new coordinates the cross terms disappear, so each eigenvalue determines the curvature along one independent direction.

A positive eigenvalue makes the form rise in that direction. A negative eigenvalue makes it fall. A zero eigenvalue makes it flat. A form with at least one positive and one negative eigenvalue is indefinite. For a symmetric PSD matrix, $x^\top Ax=0$ holds exactly when $Ax=0$, since the sum $\sum_i\lambda_i y_i^2$ can vanish only by setting every positive-eigenvalue coordinate to zero.

Repeated eigenvalues do not cause a failure of diagonalization for symmetric matrices. They make the choice of eigenvectors less unique. If the eigenspace for $\lambda$ has dimension two, any orthonormal basis of that plane is valid. The eigenspace is determined by the matrix even when the particular basis inside it is not.

#### Completing a block quadratic

For a symmetric block matrix

$$
M=\begin{pmatrix}B&C\\C^\top&D\end{pmatrix},
$$

assume $B$ is PD. Its quadratic form at $(x,y)$ can be rewritten as

$$
x^\top Bx+2x^\top Cy+y^\top Dy
=(x+B^{-1}Cy)^\top B(x+B^{-1}Cy)
+y^\top(D-C^\top B^{-1}C)y.
$$

The expansion on the right produces the original cross term and cancels the added $y^\top C^\top B^{-1}Cy$. Therefore $M$ is PSD exactly when the Schur complement $D-C^\top B^{-1}C$ is PSD; it is PD exactly when that complement is PD. For necessity choose $x=-B^{-1}Cy$ to remove the first square. For sufficiency both displayed terms are nonnegative.

This derivation also describes minimization over $x$ with $y$ fixed. The best $x$ is $-B^{-1}Cy$, and the remaining cost is the Schur-complement quadratic in $y$.

#### Positive definite and invertible

For a symmetric PD matrix, $Ax=0$ would imply $x^\top Ax=0$, forcing $x=0$. Thus it is invertible. The converse is false: $-I$ and $\operatorname{diag}(1,-1)$ are invertible but not PD. For a symmetric PSD matrix the converse does hold, because invertibility removes every zero eigenvalue and PSD removes every negative one.

Its inverse is also PD: if $A=Q\Lambda Q^\top$ with all $\lambda_i>0$, then $A^{-1}=Q\Lambda^{-1}Q^\top$ has positive eigenvalues $1/\lambda_i$. Its positive square root is $A^{1/2}=Q\Lambda^{1/2}Q^\top$. For PSD matrices the same square root exists, with zeros kept as zeros.

#### Spectral theorem in finite dimensions

For symmetric $A$, the continuous function $x^\top Ax$ has a maximum on the unit sphere, which is compact. Let $q$ be a maximizer. For every $h\perp q$, differentiate the form along $(q+th)/\|q+th\|$ at $t=0$. Its derivative is $2h^\top Aq$, so maximality gives $h^\top Aq=0$. Therefore $Aq$ is parallel to $q$ and $Aq=\lambda q$.

The orthogonal complement of $q$ is invariant: if $x\perp q$, then $q^\top Ax=(Aq)^\top x=\lambda q^\top x=0$. Restrict $A$ to that complement and repeat. Induction on dimension yields an orthonormal eigenbasis.

#### Positive semidefiniteness and eigenvalues

Put $y=Q^\top x$. Orthogonality preserves length and gives

$$
x^\top Ax=\sum_i\lambda_i y_i^2.
$$

If every eigenvalue is nonnegative, the sum is nonnegative. Conversely, taking $x=q_i$ gives $x^\top Ax=\lambda_i$, so PSD forces every eigenvalue to be nonnegative. PD is equivalent to every eigenvalue being strictly positive. The Rayleigh bounds follow by dividing the sum by $\sum_i y_i^2$: the quotient is a weighted average of the eigenvalues.

#### Symmetry cannot be inferred from a quadratic expression

Every square real matrix decomposes as $A=S+K$, where $S=(A+A^\top)/2$ is symmetric and $K=(A-A^\top)/2$ is skew-symmetric. For any real vector,

$$
x^\top Kx=(x^\top Kx)^\top=x^\top K^\top x=-x^\top Kx,
$$

so $x^\top Kx=0$. The quadratic form sees only $S$. In particular, a rotation matrix $\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ has quadratic expression zero for every $x$, but is not a symmetric positive-semidefinite operator. Whenever positivity is stated for an operator in these notes, symmetry is included.

Symmetry also ensures an orthonormal eigenbasis. A general matrix can lack any eigenbasis. For $J=\begin{pmatrix}1&1\\0&1\end{pmatrix}$, the only eigenspace is the first coordinate axis, yet the matrix acts on a two-dimensional space. Its powers are $J^n=\begin{pmatrix}1&n\\0&1\end{pmatrix}$. Knowing its eigenvalue is one does not describe the linear growth caused by the off-diagonal part.

#### A positive square root

If $A=Q\operatorname{diag}(\lambda_i)Q^\top$ is PSD, define

$$
A^{1/2}=Q\operatorname{diag}(\sqrt{\lambda_i})Q^\top.
$$

This matrix is symmetric PSD and squares to $A$. If $A$ is PD, its square root is invertible and gives $x^\top Ax=\|A^{1/2}x\|^2$. This identity turns a positive quadratic form into an ordinary squared length after a coordinate transformation.

The symmetric PSD square root is unique. Any such root $B$ commutes with $A=B^2$. Thus each eigenspace of $A$ is invariant under $B$. Diagonalize the symmetric restriction of $B$ on each one. Its eigenvalues must be nonnegative numbers whose squares equal the corresponding $\lambda_i$, forcing them to equal $\sqrt{\lambda_i}$.

A change of variables in a quadratic form gives $P^\top AP$, called a congruence, whereas an operator basis change gives $P^{-1}AP$, a similarity. Congruence preserves the signs available in the form when $P$ is invertible, but need not preserve eigenvalues. Taking $A=I$ and $P=\operatorname{diag}(2,3)$ gives form matrix $\operatorname{diag}(4,9)$.

#### Example LA-E20: whitening a covariance form

Let $\Sigma=\begin{pmatrix}5&3\\3&5\end{pmatrix}$. Find a transformation making its covariance identity.

Solution. Along $q_+=(1,1)^\top/\sqrt2$ and $q_-=(1,-1)^\top/\sqrt2$, the variances are eight and two. Set $C=\operatorname{diag}(1/\sqrt8,1/\sqrt2)Q^\top$. Then $C\Sigma C^\top=I$. For a centered random vector $X$ with this covariance, $CX$ has uncorrelated unit-variance coordinates.

This does not imply independence without an additional distributional assumption. If $X$ is jointly Gaussian, the transformed vector is Gaussian and the identity covariance does imply independent coordinates. The linear algebra establishes covariance normalization; probability assumptions determine the stronger interpretation.

### 3.3 Worked examples

#### Example LA-E07: a symmetric matrix with coupled coordinates

Diagonalize $A=\left(\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right)$.

Solution. The directions $q_1=(1,1)/\sqrt2$ and $q_2=(1,-1)/\sqrt2$ have eigenvalues $3$ and $1$. Thus $Q=(q_1,q_2)$ and $A=Q\operatorname{diag}(3,1)Q^\top$. For $x=(2,0)$ the eigenbasis coordinates are $(\sqrt2,\sqrt2)$, and $x^\top Ax=3(2)+1(2)=8$. Since both eigenvalues are positive, $A$ is PD.

#### Example LA-E08: an indefinite form

Classify $B=\left(\begin{smallmatrix}1&2\\2&1\end{smallmatrix}\right)$.

Solution. Along $(1,1)$ the eigenvalue is $3$; along $(1,-1)$ it is $-1$. Hence $x^\top Bx$ can have either sign. Positive diagonal entries have not guaranteed PSD: the off-diagonal terms make some directions negative.

#### Example LA-E09: covariance rank

Let $X$ contain $n$ observations on $p$ variables, and let $H=I-\mathbf1\mathbf1^\top/n$ center each column. For $n>1$, the sample covariance is $S=X^\top HX/(n-1)$.

Solution. Since $H^\top=H^2=H$,

$$
a^\top Sa=\|HXa\|^2/(n-1)\ge0.
$$

Also $\operatorname{rank}S=\operatorname{rank}(HX)\le\min(p,n-1)$. Therefore $p\ge n$ guarantees singularity. For $p<n$, singularity is still possible if a centered variable is an exact linear combination of others. Invertibility requires the centered columns to be independent.

### 3.4 Selected exercises

#### Problem LA-08: a tridiagonal positive operator

Source: [Axler, Exercises 7C, 2](https://linear.axler.net/LADR4e.pdf). The $4\times4$ matrix has diagonal $2$ and neighboring off-diagonals $-1$.

Solution. Expand the quadratic form:

$$
x^\top Ax=x_1^2+(x_1-x_2)^2+(x_2-x_3)^2+(x_3-x_4)^2+x_4^2.
$$

Every term is nonnegative. If the sum is zero, $x_1=x_4=0$ and all neighboring differences vanish, forcing $x=0$. Thus the matrix is PD, so zero is not an eigenvalue and the matrix is invertible. For complex coordinates the same argument uses absolute squares and the conjugate transpose.

#### Problem LA-09: positive entries without positivity of the form

Source: [Axler, Exercises 7C, 4](https://linear.axler.net/LADR4e.pdf).

Problem. For every dimension $n>1$, construct a symmetric matrix with all entries strictly positive whose quadratic form is not positive semidefinite.

Solution. For any $n>1$, let $J$ be the all-ones matrix and put $A=2J-I$. Its diagonal entries are one and its off-diagonal entries are two, all positive. It is symmetric. Yet $x=e_1-e_2$ satisfies $Jx=0$, so $x^\top Ax=-\|x\|^2=-2$. Thus $A$ is not PSD. Entrywise positivity and nonnegativity of a quadratic form are different properties.

#### Problem LA-10: shifting eigenvalues

Source: [Axler, Exercises 7C, 13(a,b)](https://linear.axler.net/LADR4e.pdf).

Problem. For a symmetric matrix $A$, characterize the scalars $\alpha$ for which $A-\alpha I$ or $\alpha I-A$ is positive semidefinite, using the extreme eigenvalues.

Solution. If $Aq_i=\lambda_iq_i$, then $(A-\alpha I)q_i=(\lambda_i-\alpha)q_i$. This matrix is PSD exactly when $\alpha\le\lambda_{\min}(A)$. Similarly, $\alpha I-A$ is PSD exactly when $\alpha\ge\lambda_{\max}(A)$. The endpoints allow zero eigenvalues. Strict inequalities make the respective matrix PD.

## 4. SVD and low-rank approximation

### Further problems on eigenvalues and quadratic forms

#### Problem LA-37: an all-ones matrix

Source: [Axler, Exercises 7C, 3](https://linear.axler.net/LADR4e.pdf).

Problem. Prove positivity of the $n\times n$ all-ones matrix and determine its rank, null space, and eigenvalues.

Solution. For $J=\mathbf1\mathbf1^\top$,

$$
x^\top Jx=(\mathbf1^\top x)^2\ge0.
$$

It is symmetric and PSD. The vector $\mathbf1$ has eigenvalue $n$, since $J\mathbf1=n\mathbf1$. Every vector whose coordinates sum to zero has eigenvalue zero. Thus the rank is one, except in the vacuous zero-dimensional setting, and the nullity is $n-1$. Positive entries throughout the matrix have not prevented a large null space.

#### Problem LA-38: adding a positive semidefinite matrix

Source: [Axler, Exercises 7C, 7](https://linear.axler.net/LADR4e.pdf).

Problem. Show that adding a positive-semidefinite operator to a positive-definite operator preserves positive definiteness and invertibility.

Solution. Let $S$ be PD and $T$ PSD. For nonzero $x$, $x^\top(S+T)x=x^\top Sx+x^\top Tx>0+0$. Hence $S+T$ is PD and invertible. A quantitative form is $\lambda_{\min}(S+T)\ge\lambda_{\min}(S)$, obtained by minimizing the Rayleigh quotient over unit vectors. Commutativity is unnecessary for this lower bound.

#### Problem LA-39: changing variables in a positive form

Source: [Axler, Exercises 7C, 9](https://linear.axler.net/LADR4e.pdf).

Problem. For PSD $T$ and any compatible linear map $S$, prove positivity of $S^\top TS$. When $T$ is PD, determine the condition for the new form to be PD.

Solution. If $T$ is PSD, then $(S^\top TS)^\top=S^\top TS$ and

$$
x^\top S^\top TSx=(Sx)^\top T(Sx)\ge0.
$$

Thus $S^\top TS$ is PSD. If $T$ is PD, this new matrix is PD precisely when $S$ is injective: any nonzero kernel vector of $S$ makes the form zero, while injectivity ensures $Sx\ne0$ for nonzero $x$. This is the algebra behind both weighted least squares and covariance transformations.

#### Problem LA-40: products of positive operators

Source: [Axler, Exercises 7C, 18](https://linear.axler.net/LADR4e.pdf).

Problem. For self-adjoint PSD operators $S,T$ on the same space, prove that $ST$ is positive exactly when $S$ and $T$ commute.

Solution. If $ST$ is positive, it is self-adjoint, so $ST=(ST)^\top=TS$. Conversely suppose PSD matrices $S,T$ commute. Each eigenspace of $S$ is invariant under $T$, because $S(Tv)=T(Sv)=\lambda Tv$. Diagonalize the symmetric restriction of $T$ within each such eigenspace. This produces a common orthonormal eigenbasis in which both matrices are diagonal with nonnegative entries. Their product has products of those entries on the diagonal, all nonnegative. Thus $ST$ is PSD. Separate positivity does not suffice unless the product is also symmetric.

#### Problem LA-41: a three-dimensional diagonalization

Constructed problem. Diagonalize

$$
A=\begin{pmatrix}2&1&0\\1&2&0\\0&0&4\end{pmatrix}
$$

and minimize its quadratic form on the unit sphere.

Solution. Choose $q_1=(1,1,0)/\sqrt2$, $q_2=(1,-1,0)/\sqrt2$, $q_3=(0,0,1)$. Their eigenvalues are $3,1,4$. For $x=a q_1+b q_2+c q_3$, the unit constraint is $a^2+b^2+c^2=1$ and the form is $3a^2+b^2+4c^2$. Its minimum is one, achieved only when $a=c=0$ and $b=\pm1$. The maximum is four at $\pm q_3$.

#### Problem LA-42: a rank-one covariance update

Constructed problem. Show that $S+vv^\top$ is PSD if $S$ is PSD, and determine when a null direction disappears.

Solution. For every $x$,

$$
x^\top(S+vv^\top)x=x^\top Sx+(v^\top x)^2\ge0.
$$

Since both terms are nonnegative, this expression is zero exactly when $x\in\ker S$ and $v^\top x=0$. Hence $\ker(S+vv^\top)=\ker S\cap v^\perp$. A one-dimensional null space generated by $z$ disappears if $v^\top z\ne0$. A null space of dimension at least two cannot be removed completely by one rank-one update, because at least one independent null direction remains perpendicular to $v$.

#### Problem LA-43: a two-by-two positivity criterion

Constructed problem. Characterize PSD matrices $A=\left(\begin{smallmatrix}a&c\\c&b\end{smallmatrix}\right)$.

Solution. Necessity gives $a,b\ge0$ by testing the coordinate vectors. If $a>0$, complete the square:

$$
ax^2+2cxy+by^2=a(x+cy/a)^2+(b-c^2/a)y^2.
$$

Nonnegativity for every $(x,y)$ requires and is implied by $ab-c^2\ge0$. If $a=0$, nonnegativity forces $c=0$: otherwise fix $y=1$ and send $x$ in the direction making $2cx+b$ negative. Thus the full PSD criterion is $a\ge0$, $b\ge0$, and $ab-c^2\ge0$. PD requires $a>0$ and $ab-c^2>0$.

#### Problem LA-44: a diagonalizable matrix without an orthogonal eigenbasis

Constructed problem. Study $A=\left(\begin{smallmatrix}1&1\\0&2\end{smallmatrix}\right)$.

Solution. The eigenvector directions are $v_1=(1,0)$ for eigenvalue one and $v_2=(1,1)$ for eigenvalue two. They are independent, so $A=V\operatorname{diag}(1,2)V^{-1}$. Their inner product is one, and scaling either vector cannot make it zero. Thus this diagonalizable matrix has no orthogonal eigenbasis. Indeed $A^\top A=\left(\begin{smallmatrix}1&1\\1&5\end{smallmatrix}\right)$ differs from $AA^\top=\left(\begin{smallmatrix}2&2\\2&4\end{smallmatrix}\right)$, so it is not normal.

### 4.1 Theory

Every real $m\times n$ matrix admits a singular value decomposition

$$
A=U\Sigma V^\top.
$$

Here $U$ and $V$ are orthogonal, and the rectangular diagonal matrix $\Sigma$ has entries $\sigma_1\ge\cdots\ge\sigma_r>0$ followed by zeros. The rank is $r$. Singular values are the square roots of the eigenvalues of $A^\top A$.

The SVD separates a map into input directions, stretches, and output directions:

$$
Av_i=\sigma_i u_i,\qquad
A=\sum_{i=1}^r\sigma_i u_iv_i^\top.
$$

The Moore-Penrose pseudoinverse reverses the nonzero stretches:

$$
A^+=\sum_{i=1}^r\sigma_i^{-1}v_iu_i^\top.
$$

Then $AA^+$ projects onto $\operatorname{range}A$, and $A^+A$ projects onto $(\ker A)^\perp$. The vector $A^+b$ has the smallest Euclidean norm among all least-squares solutions.

### 4.2 Derivations and proofs

#### Singular vectors belong to two different spaces

For a rectangular matrix, right singular vectors $v_i$ belong to the input space $\mathbb R^n$, while left singular vectors $u_i$ belong to the output space $\mathbb R^m$. The map takes $v_i$ to the scaled output direction $\sigma_i u_i$. These vectors should not be identified merely because both lists are orthonormal.

The matrix $A^\top A$ acts on the input space and identifies which input directions are stretched strongly. The matrix $AA^\top$ acts on the output space and identifies the matching output directions. Both have the same positive eigenvalues $\sigma_i^2$, but they may have different numbers of zero eigenvalues when $m\ne n$.

For $b$, first compute its output coordinates $u_i^\top b$. Divide each attainable coordinate by its singular value, then place it along $v_i$. Coordinates in the orthogonal complement of the output range cannot be produced and are discarded. This sequence is exactly the formula for $A^+b$.

#### All four Moore-Penrose identities

The SVD immediately gives

$$
AA^+A=A,\qquad A^+AA^+=A^+,\qquad
(AA^+)^\top=AA^+,\qquad(A^+A)^\top=A^+A.
$$

For the first two identities, a positive singular value is multiplied by its reciprocal and then itself again. A zero singular value remains zero. For the symmetry identities, $AA^+=\sum_{i=1}^r u_iu_i^\top$ and $A^+A=\sum_{i=1}^r v_iv_i^\top$ are sums of symmetric rank-one projectors.

To see the smallest-norm property explicitly, let $x_0=A^+b$. Every other least-squares solution is $x_0+z$ for $z\in\ker A$. The SVD puts $x_0$ in the span of the positive right singular directions, which is perpendicular to $\ker A$. Therefore $\|x_0+z\|^2=\|x_0\|^2+\|z\|^2$. The unique minimum norm occurs at $z=0$.

#### Reconstruction error in PCA

For a centered observation vector $x$, reconstruction using a unit direction $w$ is $ww^\top x$. Orthogonal decomposition gives

$$
\|x-ww^\top x\|^2=\|x\|^2-(w^\top x)^2.
$$

Average this identity over the sample. The first term is fixed, so minimizing reconstruction error is equivalent to maximizing projected variance. For $k$ orthonormal columns of $W$, the same calculation gives residual sum of squares $\|X\|_F^2-\operatorname{tr}(W^\top X^\top XW)$. The spectral argument selects the leading $k$ eigenvectors.

Centering is part of this statement. Without centering, the leading singular direction can describe the mean offset instead of variation around the mean. If a leading eigenvalue is repeated, the best principal subspace may be unique at a larger dimension even though an individual first principal direction is not.

#### Constructing the SVD

Diagonalize $A^\top A$ using the spectral theorem. It is PSD because $x^\top A^\top Ax=\|Ax\|^2$. For each positive eigenvalue $\sigma_i^2$, define $u_i=Av_i/\sigma_i$. Then

$$
u_i^\top u_j=\frac{v_i^\top A^\top Av_j}{\sigma_i\sigma_j}=\delta_{ij}.
$$

For a zero eigenvalue, $\|Av_i\|^2=0$ implies $Av_i=0$. Extend the positive-direction $u_i$ to an orthonormal basis of the output space. Applying both sides to the input basis proves the decomposition.

#### Best low-rank approximation

The truncated matrix $A_k=\sum_{i=1}^k\sigma_i u_iv_i^\top$ has rank at most $k$ and error $\|A-A_k\|_2=\sigma_{k+1}$. For any matrix $B$ of rank at most $k$, its restriction to $\operatorname{span}(v_1,\ldots,v_{k+1})$ has a nonzero kernel vector. Normalize one to $z$. Then $Bz=0$ and

$$
\|A-B\|_2\ge\|(A-B)z\|=\|Az\|\ge\sigma_{k+1}.
$$

This proves spectral-norm optimality. In Frobenius norm, orthogonal projection of $A$ onto a fixed output subspace $S$ is $P_SA$, with residual squared $\|A\|_F^2-\operatorname{tr}(P_SAA^\top)$. Among spaces of dimension at most $k$, the trace is largest for the first $k$ left singular vectors: in an eigenbasis it is $\sum_i\sigma_i^2p_i$, where $0\le p_i\le1$ and $\sum_i p_i\le k$. Allocating those weights to the largest eigenvalues gives

$$
\min_{\operatorname{rank}B\le k}\|A-B\|_F^2=\sum_{i>k}\sigma_i^2.
$$

#### PCA as maximum variance

For centered data with covariance $S$, the variance along a unit direction $w$ is $w^\top Sw$. Maximizing it subject to $w^\top w=1$ gives $Sw=\lambda w$ by a Lagrange multiplier. The Rayleigh bound identifies the largest eigenvalue as the maximum. Subsequent directions maximize variance under orthogonality to earlier ones. For centered $X=U\Sigma V^\top$, the principal directions are the columns of $V$ and the variance explained by direction $i$ is $\sigma_i^2/(n-1)$.

#### Centering, covariance rank, and reconstruction

Let $X$ be an $n\times p$ data matrix with observations in rows. The centering matrix is $H=I-\mathbf1\mathbf1^\top/n$. It is symmetric and idempotent, kills the constant vector, and has rank $n-1$. The centered matrix $Z=HX$ has zero column means.

The sample covariance is $S=Z^\top Z/(n-1)$ for $n>1$. Since a Gram matrix has the same rank as its underlying matrix,

$$
\operatorname{rank}(S)=\operatorname{rank}(Z)
\le\min(p,n-1).
$$

Thus more features than independent centered observations necessarily produce a singular covariance. Even with many observations, exact column relations can lower the rank further.

For a unit direction $q$, the rank-one reconstruction of all centered rows is $Zqq^\top$. Expanding its squared error gives

$$
\|Z-Zqq^\top\|_F^2
=\operatorname{tr}(Z^\top Z)-q^\top Z^\top Zq.
$$

Minimizing reconstruction error is therefore equivalent to maximizing projected variance. With $k$ orthonormal directions collected in $Q_k$, the captured sum is $\operatorname{tr}(Q_k^\top Z^\top ZQ_k)$.

To justify choosing the top eigenvectors, let $v_i$ be an orthonormal eigenbasis with descending eigenvalues $\lambda_i$. Put $w_i=\|Q_k^\top v_i\|^2$. These numbers satisfy $0\le w_i\le1$ and $\sum_iw_i=k$. The captured variance is $\sum_i\lambda_iw_i$, which cannot exceed $\sum_{i=1}^k\lambda_i$: moving any weight from a lower eigenvalue to an available higher one cannot decrease the sum. Equality is attained by the leading eigenvectors.

If the cutoff eigenvalue is tied, different subspaces within the tied eigenspace can be equally good. Individual eigenvector signs never change the reconstruction. PCA identifies directions or subspaces, rather than a unique signed list in every dataset.

#### Example LA-E21: covariance PCA and rescaled PCA

Use the four observations $(2,0),(-2,0),(0,1),(0,-1)$.

Solution. The mean is zero and the scatter matrix is $Z^\top Z=\operatorname{diag}(8,2)$. The sample covariance eigenvalues are $8/3$ and $2/3$. The first coordinate is the leading principal direction and captures $8/(8+2)=80\%$ of the variance.

Keeping only that coordinate reconstructs the first two observations exactly and maps the other two to zero. Their squared errors are one each, totaling two, equal to the discarded scatter eigenvalue.

If each coordinate is first divided by its sample standard deviation, the covariance becomes identity. Both standardized directions now have equal variance; there is no unique leading direction. Standardizing changed the metric in which reconstruction errors are measured. A PCA result therefore needs its centering and scaling convention stated along with the computed eigenvectors.

### 4.3 Worked examples

#### Example LA-E10: an explicit rectangular SVD

Let $A=\left(\begin{smallmatrix}3&0\\0&2\\0&0\end{smallmatrix}\right)$.

Solution. The coordinate directions already give the SVD, with singular values $3,2$. The best rank-one approximation keeps only the entry $3$. Its spectral and Frobenius errors are both $2$. The pseudoinverse is $\left(\begin{smallmatrix}1/3&0&0\\0&1/2&0\end{smallmatrix}\right)$, so $A^+(6,6,7)^\top=(2,3)^\top$ and the residual $(0,0,7)^\top$ cannot be fitted.

#### Example LA-E11: rank-deficient least squares

Minimize $\|Ax-b\|^2$ for $A=\left(\begin{smallmatrix}1&1\\2&2\end{smallmatrix}\right)$ and $b=(1,0)^\top$.

Solution. Only $s=x_1+x_2$ affects the fit. Minimize $(s-1)^2+4s^2$, whose derivative is $10s-2$. Hence $s=1/5$. All minimizers are $(t,1/5-t)$; the smallest norm occurs at equal coordinates, $x^+=(1/10,1/10)$. The fitted vector is $(1/5,2/5)$, and the residual $(4/5,-2/5)$ is perpendicular to $(1,2)$.

#### Example LA-E12: two principal directions

Let $S=\left(\begin{smallmatrix}5&4\\4&5\end{smallmatrix}\right)$.

Solution. The unit eigenvectors are $(1,1)/\sqrt2$ and $(1,-1)/\sqrt2$, with eigenvalues $9$ and $1$. The first component explains $9/(9+1)=90\%$ of total variance. A one-dimensional reconstruction keeps the common movement of the coordinates and removes their difference. If the variables have different units, standardizing them changes the covariance matrix and can change these directions.

### 4.4 Selected exercises

#### Problem LA-11: pseudoinverse of a surjective map

Source: [Axler, Exercises 6C, 21(a,b,c)](https://linear.axler.net/LADR4e.pdf). Let $T(a,b,c)=(a+b+c,2b+3c)$.

Solution. Its matrix is $A=\left(\begin{smallmatrix}1&1&1\\0&2&3\end{smallmatrix}\right)$. Independent rows make it surjective. Therefore

$$
A^+=A^\top(AA^\top)^{-1},\quad
AA^\top=\begin{pmatrix}3&5\\5&13\end{pmatrix},\quad
(AA^\top)^{-1}=\frac1{14}\begin{pmatrix}13&-5\\-5&3\end{pmatrix}.
$$

Multiplying yields

$$
A^+\binom xy=\frac1{14}
\begin{pmatrix}13x-5y\\3x+y\\-2x+4y\end{pmatrix}.
$$

Direct multiplication gives $AA^+=I_2$, the projection onto the full range. A kernel basis is $n=(1,-3,2)^\top$, with $\|n\|^2=14$. Another multiplication gives

$$
A^+A=I_3-\frac{nn^\top}{14}.
$$

This is the projection onto the kernel's orthogonal complement, establishing both requested identities and the minimum-norm property.

#### Problem LA-12: singular values without eigenvalue calculations

Source: [Axler, Exercises 7E, 5](https://linear.axler.net/LADR4e.pdf). For $T(x,y)=(-4y,x)$, find the singular values.

Solution. Its matrix satisfies $A^\top A=\operatorname{diag}(1,16)$, so the singular values in decreasing order are $4,1$. One SVD uses $v_1=e_2,v_2=e_1$ and $u_1=-e_1,u_2=e_2$. Then $Av_1=4u_1$ and $Av_2=u_2$. The eigenvalues of $A$ itself are imaginary, which does not affect the real nonnegative singular values.

#### Problem LA-13: singular values of a self-adjoint operator

Source: [Axler, Exercises 7E, 7](https://linear.axler.net/LADR4e.pdf).

Problem. Relate the singular values of a self-adjoint operator to its eigenvalues. Extend the argument to a complex normal operator.

Solution. For a real symmetric matrix, $A=Q\Lambda Q^\top$ gives $A^\top A=Q\Lambda^2Q^\top$. The singular values are therefore $|\lambda_i|$, sorted in decreasing order. For a complex normal operator, the spectral theorem gives $A=Q\Lambda Q^*$ and $A^*A=Q|\Lambda|^2Q^*$, proving the same conclusion. Self-adjointness or normality is needed: a nonzero nilpotent matrix can have only zero eigenvalues and still have a positive singular value.

## 5. Conditioning and regularization

### Further problems on SVD and approximation

#### Problem LA-45: transposition preserves positive singular values

Source: [Axler, Exercises 7E, 9](https://linear.axler.net/LADR4e.pdf).

Problem. Prove that a matrix and its transpose have the same positive singular values, including the rectangular case.

Solution. Transposing $A=\sum_{i=1}^r\sigma_i u_iv_i^\top$ gives $A^\top=\sum_{i=1}^r\sigma_i v_iu_i^\top$. Both orthonormal lists simply exchange roles, leaving the positive singular values unchanged. If the matrix is rectangular, one of $A^\top A$ and $AA^\top$ has additional zeros, which is why the claim is stated for positive singular values.

#### Problem LA-46: singular values of an inverse

Source: [Axler, Exercises 7E, 10](https://linear.axler.net/LADR4e.pdf).

Problem. Express the singular values of an invertible matrix's inverse in terms of those of the original matrix, with the correct ordering.

Solution. If $A$ is square and invertible, all singular values are positive. Its inverse is $A^{-1}=V\Sigma^{-1}U^\top$. Thus the singular values are the reciprocals, reordered as $1/\sigma_n\ge\cdots\ge1/\sigma_1$. In particular, the inverse amplifies most strongly the output direction associated with the smallest original singular value.

#### Problem LA-47: a lower singular-value bound

Source: [Axler, Exercises 7E, 14](https://linear.axler.net/LADR4e.pdf).

Problem. Prove $\|Ax\|_2\ge\sigma_{\min}\|x\|_2$, including zero singular directions when the map is not injective.

Solution. Expand $x=\sum_i c_iv_i$ in the full right singular basis, including zero singular directions. Then

$$
\|Ax\|^2=\sum_i\sigma_i^2c_i^2
\ge\sigma_{\min}^2\sum_i c_i^2
=\sigma_{\min}^2\|x\|^2.
$$

Taking square roots proves the bound. When $\sigma_{\min}>0$, the map is injective and cannot collapse any nonzero direction. When $\sigma_{\min}=0$, the inequality still holds but provides no positive lower bound.

#### Problem LA-48: transposing a pseudoinverse

Source: [Axler, Exercises 7E, 16](https://linear.axler.net/LADR4e.pdf).

Problem. Prove $(A^\top)^+=(A^+)^\top$ for a matrix that may be rank deficient.

Solution. From the SVD, $A^+=\sum_i\sigma_i^{-1}v_iu_i^\top$. Its transpose is $\sum_i\sigma_i^{-1}u_iv_i^\top$. Applying the pseudoinverse construction to $A^\top=\sum_i\sigma_i v_iu_i^\top$ produces that same expression. Hence $(A^\top)^+=(A^+)^\top$. The argument includes rank-deficient matrices because the sums only use positive singular values.

#### Problem LA-49: a non-diagonal rank-one approximation

Constructed problem. Find the best rank-one approximation to $A=\left(\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right)$.

Solution. Its symmetric eigendecomposition has eigenvalues $3,1$ with unit directions $q_1=(1,1)/\sqrt2$, $q_2=(1,-1)/\sqrt2$. Both eigenvalues are positive, so this is also an SVD. Keep the first term:

$$
A_1=3q_1q_1^\top
=\frac32\begin{pmatrix}1&1\\1&1\end{pmatrix}.
$$

The residual is $q_2q_2^\top=\tfrac12\left(\begin{smallmatrix}1&-1\\-1&1\end{smallmatrix}\right)$, with spectral and Frobenius norm one. The best-approximation theorem proves that no rank-one matrix can reduce either error below one.

#### Problem LA-50: a rank-one pseudoinverse

Constructed problem. For nonzero $a\in\mathbb R^m$ and $b\in\mathbb R^n$, find $(ab^\top)^+$.

Solution. Normalize the vectors: $u=a/\|a\|$, $v=b/\|b\|$. Then $ab^\top=(\|a\|\|b\|)uv^\top$ is a one-term SVD. Reversing its only nonzero singular value gives

$$
(ab^\top)^+=\frac{ba^\top}{\|a\|^2\|b\|^2}.
$$

Multiplication gives $AA^+=aa^\top/\|a\|^2$ and $A^+A=bb^\top/\|b\|^2$, the appropriate line projectors. If either original vector is zero, the matrix and its pseudoinverse are zero, and the displayed division formula must not be used.

#### Problem LA-51: PCA from actual observations

Constructed problem. The four centered observations are $(2,1),(-2,-1),(1,2),(-1,-2)$. Find both principal directions and explained-variance shares.

Solution. The column sums vanish. Summing outer products gives

$$
X^\top X=\begin{pmatrix}10&8\\8&10\end{pmatrix},\qquad
S=\frac13\begin{pmatrix}10&8\\8&10\end{pmatrix}.
$$

The eigenvectors are $(1,1)/\sqrt2$ and $(1,-1)/\sqrt2$, with covariance eigenvalues $6$ and $2/3$. Total variance is $20/3$, so the shares are $9/10$ and $1/10$. The leading scores are $3/\sqrt2,-3/\sqrt2,3/\sqrt2,-3/\sqrt2$. Reconstructing each observation from its leading score produces $(3/2,3/2)$ or its negative. The residual sum of squares is the discarded eigenvalue of $X^\top X$, namely two.

#### Problem LA-52: nonuniqueness of a best rank-one approximation

Constructed problem. Find two distinct best rank-one approximations to $I_2$.

Solution. Both $\operatorname{diag}(1,0)$ and $\operatorname{diag}(0,1)$ have rank one and leave a residual with spectral norm one and Frobenius norm one. The two singular values of $I_2$ are tied at one, so either direction can be retained. More generally $ww^\top$ is best for every unit $w$. Its residual is the orthogonal projector onto the perpendicular line. A spectral gap between the last retained and first discarded singular values removes this particular ambiguity in the optimal subspace.

### 5.1 Theory

The induced matrix norm is $\|A\|=\sup_{\|x\|=1}\|Ax\|$. For the Euclidean norm it equals the largest singular value. For invertible $A$,

$$
\kappa_2(A)=\|A\|_2\|A^{-1}\|_2
=\frac{\sigma_{\max}}{\sigma_{\min}}\ge1.
$$

Conditioning describes sensitivity of the mathematical problem to changes in its inputs. Stability describes how an algorithm handles rounding. A stable algorithm can still produce a sensitive answer when the problem is ill-conditioned.

Ridge regularization minimizes $\|Ax-b\|^2+\lambda\|x\|^2$ for $\lambda>0$. It adds curvature in directions where the data give little or no information. It also changes the objective being solved; reducing sensitivity generally introduces shrinkage.

### 5.2 Derivations and proofs

#### Absolute and relative sensitivity

An absolute perturbation bound compares $\|\delta x\|$ with the size of the input change. A relative bound divides by the original sizes. These answer different questions when the true output is small. Relative error is undefined at a zero exact output, so an absolute error or a scale chosen independently of that output is needed there.

For fixed invertible $A$, $x=A^{-1}b$ and $\delta x=A^{-1}\delta b$ exactly. Thus the absolute condition number with respect to $b$ in the chosen norm is $\|A^{-1}\|$. At a particular nonzero $b$, the relative condition number is $\|A^{-1}\|\|b\|/\|A^{-1}b\|$, which is at most $\kappa(A)$. The usual matrix condition number is a bound over right-hand-side directions; it need not equal the sensitivity at every particular $b$.

#### Normal equations square the condition number

For a full-column-rank matrix, the positive eigenvalues of $A^\top A$ are $\sigma_i^2$. Hence

$$
\kappa_2(A^\top A)=\frac{\sigma_{\max}^2}{\sigma_{\min}^2}
=\kappa_2(A)^2.
$$

The normal equations are exact algebra, but explicitly forming $A^\top A$ in finite precision can lose information in small singular directions. QR works with orthogonal transformations whose Euclidean condition number is one. It avoids this squaring at the level of the matrix being solved, although the underlying least-squares problem can still be sensitive.

#### Choosing a regularization scale

The ridge transition occurs around $\sigma_i^2\approx\lambda$. If $\sigma_i^2\gg\lambda$, then $\sigma_i/(\sigma_i^2+\lambda)$ is close to $1/\sigma_i$. If $\sigma_i^2\ll\lambda$, then the factor is approximately $\sigma_i/\lambda$, which strongly suppresses that direction.

This also means that units affect the penalty. Rescaling one column of $A$ changes the numerical size of its coefficient without changing the fitted values available to the unpenalized model. A penalty $\lambda\|x\|^2$ is applied to coefficients in the chosen units, so rescaling columns can change the regularized fitted result. Standardization or a deliberately chosen matrix penalty $x^\top Rx$ makes that choice explicit.

#### Perturbing a linear system

Suppose $Ax=b$ and $(A+E)(x+\delta x)=b+e$. Subtracting gives

$$
(A+E)\delta x=e-Ex.
$$

If $\|A^{-1}E\|<1$, the series $\sum_{j\ge0}(-A^{-1}E)^j$ converges and inverts $I+A^{-1}E$. Thus

$$
\|\delta x\|\le
\frac{\|A^{-1}\|}{1-\|A^{-1}E\|}
(\|e\|+\|E\|\|x\|).
$$

For $b\ne0$, $x\ne0$, and $\kappa(A)\|E\|/\|A\|<1$, this implies

$$
\frac{\|\delta x\|}{\|x\|}
\le
\frac{\kappa(A)}
{1-\kappa(A)\|E\|/\|A\|}
\left(\frac{\|e\|}{\|b\|}+\frac{\|E\|}{\|A\|}\right).
$$

We used $\|b\|\le\|A\|\|x\|$ to express the right-hand-side perturbation relatively. This is an upper bound; particular perturbation directions may have much smaller effects.

#### Ridge in singular directions

Differentiating the ridge objective gives $(A^\top A+\lambda I)x=A^\top b$. Its Hessian is $2(A^\top A+\lambda I)$ and is PD, so the unique solution is

$$
x_\lambda=\sum_{i=1}^r
\frac{\sigma_i}{\sigma_i^2+\lambda}(u_i^\top b)v_i.
$$

Unregularized inversion uses $1/\sigma_i$; ridge replaces it by the bounded factor $\sigma_i/(\sigma_i^2+\lambda)$. The Hessian eigenvalues are $2(\sigma_i^2+\lambda)$, including $2\lambda$ in null directions. Consequently its condition number is $(\sigma_{\max}^2+\lambda)/(\sigma_{\min}^2+\lambda)$, with $\sigma_{\min}=0$ if the column rank is deficient.

#### General quadratic regularization

Consider $\min_x\|Ax-b\|^2+\lambda\|Lx\|^2$ with $\lambda>0$. The matrix $L$ can penalize selected coordinates or differences between neighboring coordinates. Its normal equation is

$$
(A^\top A+\lambda L^\top L)x=A^\top b.
$$

For any direction $h$, the associated quadratic form is $\|Ah\|^2+\lambda\|Lh\|^2$. It vanishes exactly when $h$ belongs to both kernels. Consequently the objective has a unique minimizer precisely when

$$
\ker A\cap\ker L=\{0\}.
$$

Penalizing some directions does not cure an unidentified direction left unpenalized by $L$. The standard ridge choice $L=I$ always gives a positive-definite system, but a penalty excluding an intercept requires checking the joint-kernel condition.

When $L=I$ and $A=U\Sigma V^\top$, ridge uses filter coefficients $\sigma_i/(\sigma_i^2+\lambda)$ on the left-singular data coordinates. The pseudoinverse instead uses $1/\sigma_i$ on positive singular values. Ridge suppresses amplification in small-singular-value directions by accepting some residual error.

#### Bias and variance are separate effects

Suppose $b=Ax_0+\epsilon$ with $E[\epsilon]=0$ and covariance $\sigma_\epsilon^2I$. The ridge estimator has expectation

$$
E[\hat x_\lambda]
=(A^\top A+\lambda I)^{-1}A^\top Ax_0.
$$

In a right-singular direction, the true coefficient is multiplied by $\sigma_i^2/(\sigma_i^2+\lambda)$. This produces shrinkage bias. The sampling variance in that direction is

$$
\sigma_\epsilon^2\frac{\sigma_i^2}{(\sigma_i^2+\lambda)^2}.
$$

The mean squared error is squared bias plus variance, summed over orthogonal directions. A null-space component of $x_0$ is not identified by the data; ridge sets it to zero and retains its squared contribution as bias. Increasing the penalty reduces variance but need not reduce total error for every possible true vector.

#### Example LA-E22: penalizing differences instead of levels

Let $A=(1,1)$ and $L=(1,-1)$. Minimize $(x_1+x_2-b)^2+\lambda(x_1-x_2)^2$.

Solution. The data identify the sum, while the penalty controls the difference. Their kernels are the opposite lines and intersect only at zero. The unique solution is $x_1=x_2=b/2$ for every positive $\lambda$, with objective zero. The normal matrix is

$$
\begin{pmatrix}1+\lambda&1-\lambda\\1-\lambda&1+\lambda\end{pmatrix},
$$

whose eigenvalues are two in the common direction and $2\lambda$ in the difference direction. If $L$ were instead equal to $A$, both terms would ignore the difference direction and uniqueness would fail. The direction penalized matters as much as the size of the penalty.

### 5.3 Worked examples

#### Example LA-E13: a sensitive direction

For $A=\operatorname{diag}(1,10^{-6})$ and $b=(1,10^{-6})$, the solution is $(1,1)$. Change only the second coordinate of $b$ by $10^{-6}$.

Solution. The new solution is $(1,2)$. An input perturbation of norm $10^{-6}$ has produced an output perturbation of norm one. The condition number is $10^6$. With $\lambda=10^{-4}$, the second ridge coordinate is $10^{-12}/(10^{-12}+10^{-4})$, about $10^{-8}$. It is insensitive but strongly shrunk relative to the unregularized value.

#### Example LA-E14: covariance shrinkage

Let $S=Q\operatorname{diag}(4,0)Q^\top$. Add $\lambda I$.

Solution. The eigenvectors stay the same and the eigenvalues become $4+\lambda,\lambda$. For $\lambda=1/10$, the condition number is $4.1/0.1=41$. Adding a scalar multiple of the identity makes a PSD covariance matrix invertible. Shrinkage toward a different target need not preserve the original eigenvectors.

#### Example LA-E15: a small residual with a large error

Use $A=\operatorname{diag}(1,\varepsilon)$, $b=(1,0)$, and an approximate solution $\widetilde x=(1,1)$.

Solution. The exact solution is $(1,0)$. The error has norm one, but the residual $b-A\widetilde x=(0,-\varepsilon)$ has norm $\varepsilon$. The relation $x-\widetilde x=A^{-1}r$ explains this. Reporting a tiny residual without the conditioning can overstate accuracy.

### 5.4 Selected exercises

#### Problem LA-14: reverse triangle inequality for operators

Source: [Axler, Exercises 7F, 1](https://linear.axler.net/LADR4e.pdf).

Problem. Prove the reverse triangle inequality $\lvert\|S\|-\|T\|\rvert\le\|S-T\|$ for induced operator norms.

Solution. For unit $x$, $\|Sx\|\le\|(S-T)x\|+\|Tx\|\le\|S-T\|+\|T\|$. Taking the supremum gives $\|S\|-\|T\|\le\|S-T\|$. Interchanging $S,T$ gives the opposite signed inequality. Combining them proves $|\|S\|-\|T\||\le\|S-T\|$.

#### Problem LA-15: stability of invertibility

Source: [Axler, Exercises 7F, 8(a,b)](https://linear.axler.net/LADR4e.pdf).

Problem. First show that $\|I-T\|<1$ implies invertibility. Then prove that every sufficiently small operator-norm perturbation of an invertible $S$ remains invertible, giving an explicit radius.

Solution. If $\|I-T\|<1$, the geometric operator series $R=\sum_{k=0}^\infty(I-T)^k$ converges in norm. Multiplying its finite partial sum by $T$ gives $I-(I-T)^{N+1}$, which tends to $I$. Hence $RT=TR=I$.

For invertible $S$, factor $T=S[I+S^{-1}(T-S)]$. The hypothesis $\|T-S\|<1/\|S^{-1}\|$ makes the bracket's difference from $I$ have norm below one. Part (a) makes the bracket invertible, and the product is invertible. This gives a concrete neighborhood of $S$ containing no singular matrices.

#### Problem LA-16: positive definiteness under perturbation

Source: [Axler, Exercises 7F, 16](https://linear.axler.net/LADR4e.pdf).

Problem. Let $S$ be a positive-definite self-adjoint operator. Find an operator-norm neighborhood within the self-adjoint operators in which positivity is preserved.

Solution. Let $S$ be PD, let $m=\lambda_{\min}(S)>0$, and let $T$ be symmetric with $\|T-S\|_2<m$. For nonzero $x$,

$$
x^\top Tx=x^\top Sx+x^\top(T-S)x
\ge(m-\|T-S\|_2)\|x\|^2>0.
$$

Thus one can choose $\delta=m$, and in fact nearby $T$ is PD as well as PSD. Symmetry is needed to call $T$ a positive operator; a small nonsymmetric perturbation can preserve a positive quadratic form without producing a self-adjoint operator.

## 6. Matrix differentiation

### Further problems on sensitivity and regularization

#### Problem LA-53: submultiplicativity

Source: [Axler, Exercises 7F, 5](https://linear.axler.net/LADR4e.pdf).

Problem. Prove submultiplicativity $\|ST\|\le\|S\|\|T\|$ for compatible induced operator norms.

Solution. For every unit vector $x$, $\|STx\|\le\|S\|\|Tx\|\le\|S\|\|T\|$. Taking the supremum over $x$ gives $\|ST\|\le\|S\|\|T\|$. If $Tx=0$, the first inequality still holds. The proof uses induced norms with compatible input, intermediate, and output vector norms; arbitrary matrix-entry norms need not obey the same formula without an additional constant.

#### Problem LA-54: the norm of a Gram matrix

Source: [Axler, Exercises 7F, 19](https://linear.axler.net/LADR4e.pdf).

Problem. Prove $\|A^\top A\|_2=\|A\|_2^2$ and explain the role of the Gram matrix's positivity.

Solution. The eigenvalues of $A^\top A$ are $\sigma_i^2$, all nonnegative, so its operator norm is $\sigma_{\max}^2=\|A\|_2^2$. Alternatively, the Rayleigh quotient gives $\|A^\top A\|_2=\max_{\|x\|=1}x^\top A^\top Ax=\max_{\|x\|=1}\|Ax\|^2$. The second equality relies on symmetry and PSD of the Gram matrix.

#### Problem LA-55: the smallest possible condition number

Source: [Axler, Exercises 7F, 24](https://linear.axler.net/LADR4e.pdf).

Problem. Characterize the invertible matrices satisfying $\|A^{-1}\|_2=\|A\|_2^{-1}$, and relate the result to condition number one.

Solution. For invertible $A$, $\|A^{-1}\|_2=\|A\|_2^{-1}$ exactly when $\sigma_{\max}/\sigma_{\min}=1$. Since every singular value lies between those two, they must all equal a common positive number $s=\|A\|_2$. The SVD then gives $A=sUV^\top$, so $A/\|A\|_2$ is orthogonal. Conversely a positive scalar times an orthogonal matrix has all singular values equal and condition number one.

#### Problem LA-56: product order and operator norm

Source: [Axler, Exercises 7F, 6](https://linear.axler.net/LADR4e.pdf).

Problem. Determine whether $\|ST\|_2=\|TS\|_2$ must hold for two square matrices. Supply a proof or a counterexample.

Solution. The assertion $\|ST\|=\|TS\|$ is false. Set

$$
S=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
T=\begin{pmatrix}0&1\\0&0\end{pmatrix}.
$$

Then $ST=T$, with norm one, but $TS=0$, with norm zero. Matrix products act in order: $TS$ first removes the second input coordinate, leaving nothing for $T$ to transfer into the first coordinate.

#### Problem LA-57: a nearly singular pair of equations

Constructed problem. Solve $x+y=2$, $x+(1+\varepsilon)y=2+\varepsilon$ for $\varepsilon>0$, and then perturb the second right-hand side by $\delta$.

Solution. Subtraction gives $\varepsilon y=\varepsilon$, so the exact solution is $(1,1)$. After perturbation, $\varepsilon\widetilde y=\varepsilon+\delta$, so $\widetilde y=1+\delta/\varepsilon$ and $\widetilde x=1-\delta/\varepsilon$. The error norm is $\sqrt2|\delta|/\varepsilon$. The equations become nearly duplicates as $\varepsilon\downarrow0$, and the subtraction needed to identify $y$ divides a small difference by another small number.

#### Problem LA-58: a residual-based error bound

Constructed problem. If $A$ is PD with $\lambda_{\min}(A)=m>0$, bound the error of an approximate solution using its residual.

Solution. With $r=b-A\widetilde x$ and exact solution $x$, $A(x-\widetilde x)=r$. The inverse norm is $1/m$, so $\|x-\widetilde x\|_2\le\|r\|_2/m$. Equality is possible when $r$ lies along an eigenvector for $m$. To guarantee solution error at most $\eta$, the residual should therefore be at most $m\eta$. A stopping rule based only on a fixed residual tolerance can fail to control error uniformly across differently conditioned systems.

#### Problem LA-59: ridge bias in each direction

Constructed problem. Suppose $b=Ax_0$ exactly and compute $x_\lambda-x_0$.

Solution. Expand $x_0=\sum_i c_iv_i$ in a full right singular basis. For each positive singular direction,

$$
(x_\lambda)_i=\frac{\sigma_i^2}{\sigma_i^2+\lambda}c_i,
\qquad
(x_\lambda-x_0)_i=-\frac{\lambda}{\sigma_i^2+\lambda}c_i.
$$

For a zero singular direction, $x_\lambda$ has coefficient zero and the error is $-c_i$. Therefore exact data do not make ridge unbiased relative to an arbitrary true coefficient vector. As $\lambda\downarrow0$, the solution approaches $A^+b$, which recovers only the component of $x_0$ perpendicular to the kernel.

#### Problem LA-60: shrinkage toward a general positive target

Constructed problem. Let $S$ be PSD, $T$ PD, and $S_\alpha=(1-\alpha)S+\alpha T$ for $0<\alpha\le1$. Prove invertibility and give a lower eigenvalue bound.

Solution. For a unit vector $x$,

$$
x^\top S_\alpha x=(1-\alpha)x^\top Sx+\alpha x^\top Tx
\ge\alpha\lambda_{\min}(T)>0.
$$

Thus $S_\alpha$ is PD and $\lambda_{\min}(S_\alpha)\ge\alpha\lambda_{\min}(T)$. Its largest eigenvalue is at most $(1-\alpha)\lambda_{\max}(S)+\alpha\lambda_{\max}(T)$, giving a condition-number upper bound by dividing by the earlier lower bound. Unlike an identity target, a general $T$ need not share eigenvectors with $S$.

### 6.1 Theory

For a differentiable scalar function $f:\mathbb R^n\to\mathbb R$, the gradient is defined by

$$
df=\nabla f(x)^\top dx.
$$

The Hessian is the matrix of second partial derivatives. For twice continuously differentiable $f$, it is symmetric. A positive-definite Hessian everywhere on a convex domain makes $f$ strictly convex; PSD gives convexity.

Useful identities are

$$
\nabla_x(a^\top x)=a,\qquad
\nabla_x(x^\top Ax)=(A+A^\top)x,\qquad
\nabla_x\|Ax-b\|^2=2A^\top(Ax-b).
$$

For a matrix variable, the Frobenius inner product is $\langle G,H\rangle_F=\operatorname{tr}(G^\top H)$. Define the matrix gradient by $df=\operatorname{tr}((\nabla_Xf)^\top dX)$. This convention avoids ambiguous numerator-layout and denominator-layout derivative tables.

### 6.2 Derivations and proofs

#### Dimensions before derivatives

For $A\in\mathbb R^{m\times n}$ and $x\in\mathbb R^n$, the residual $Ax-b$ has $m$ coordinates. Its gradient with respect to $x$ must have $n$ coordinates. The expression $A^\top(Ax-b)$ has the right size; $A(Ax-b)$ may not even be defined. Checking dimensions catches many derivative mistakes before any algebra.

For $X\in\mathbb R^{n\times k}$, a gradient with respect to $X$ has the same shape as $X$. The differential convention $\operatorname{tr}(G^\top dX)$ makes this visible. A trace is a scalar, and cyclic rearrangement is allowed only when the resulting products are defined. One cannot arbitrarily swap factors: $\operatorname{tr}(ABC)=\operatorname{tr}(BCA)$, but usually $\operatorname{tr}(ABC)\ne\operatorname{tr}(ACB)$.

#### A quadratic minimum along feasible directions

Positive definiteness of $H$ on all of $\mathbb R^n$ is sufficient for the constrained quadratic, but can be stronger than needed. Let the columns of $Z$ form a basis of $\ker C^\top$, and let $x_p$ be one feasible point. Every feasible $x$ is $x_p+Zz$. Substitution turns the problem into

$$
\frac12 z^\top(Z^\top HZ)z+
z^\top Z^\top(Hx_p-g)+\text{constant}.
$$

If $Z^\top HZ$ is PD, this reduced problem has a unique minimizer even if $H$ is singular or indefinite in infeasible directions. Its stationarity equation is

$$
Z^\top HZ\,z=Z^\top(g-Hx_p).
$$

If the reduced Hessian is merely PSD, use the range condition from Problem LA-17. A linear term outside its range makes the objective unbounded below along an allowed flat direction. This approach handles singular covariance examples without an unjustified matrix inverse.

#### Differentiating an inverse and a log determinant

For an invertible matrix $X$, differentiate $XX^{-1}=I$:

$$
(dX)X^{-1}+X\,d(X^{-1})=0,
\qquad d(X^{-1})=-X^{-1}(dX)X^{-1}.
$$

For PD $X$, the determinant expansion $\det(I+\varepsilon H)=1+\varepsilon\operatorname{tr}H+O(\varepsilon^2)$ gives

$$
d\log\det X=\operatorname{tr}(X^{-1}dX).
$$

To justify that expansion, use the permutation definition of a determinant. The linear terms come from choosing a single diagonal $\varepsilon H_{ii}$ and all other diagonal ones; any nonidentity permutation requires at least two off-diagonal factors and contributes at order at least two. The gradient of $\log\det X$ is $X^{-\top}$, which equals $X^{-1}$ on symmetric PD matrices.

#### Differentiating quadratic forms

Expand in a direction $h$:

$$
(x+h)^\top A(x+h)-x^\top Ax
=h^\top Ax+x^\top Ah+h^\top Ah.
$$

The linear term is $[(A+A^\top)x]^\top h$ and the remainder is $O(\|h\|^2)$. Thus the gradient is $(A+A^\top)x$ and the Hessian is $A+A^\top$. For $f(x)=\tfrac12x^\top Hx-g^\top x+c$ with symmetric PD $H$,

$$
\nabla f=Hx-g,\quad x^*=H^{-1}g,\quad
f(x)-f(x^*)=\tfrac12(x-x^*)^\top H(x-x^*).
$$

The last identity proves global optimality and uniqueness.

#### Equality-constrained quadratic

Minimize $\tfrac12x^\top Hx-g^\top x$ subject to $C^\top x=d$, where $H$ is PD and $C$ has independent columns. With multiplier $\nu$,

$$
L=\tfrac12x^\top Hx-g^\top x+\nu^\top(C^\top x-d).
$$

Stationarity gives $Hx-g+C\nu=0$, so $x=H^{-1}(g-C\nu)$. Substitution yields

$$
\nu=(C^\top H^{-1}C)^{-1}(C^\top H^{-1}g-d).
$$

The inverse exists because $z^\top C^\top H^{-1}Cz>0$ for nonzero $z$. For a feasible displacement $h$ with $C^\top h=0$, the linear objective change vanishes at this point and the remaining change is $\tfrac12h^\top Hh>0$. This proves uniqueness on the feasible set.

#### Matrix-valued unknowns and the Frobenius gradient

For a matrix variable $X$, a differential is a scalar linear function of a perturbation $H$. The Frobenius gradient $G$ is defined by $Df(X)[H]=\operatorname{tr}(G^\top H)$. This convention specifies both the dimensions and transpose placement.

Let $A\in\mathbb R^{m\times p}$, $X\in\mathbb R^{p\times q}$, $B\in\mathbb R^{q\times n}$, and $C\in\mathbb R^{m\times n}$. For

$$
f(X)=\frac12\|AXB-C\|_F^2,
$$

put $R=AXB-C$. Expanding $f(X+hH)$ to first order gives $Df(X)[H]=\operatorname{tr}(R^\top AHB)$. Cyclic movement inside the trace yields

$$
\nabla_Xf=A^\top(AXB-C)B^\top.
$$

The gradient has size $p\times q$, matching $X$. Its Hessian acts on a perturbation as $H\mapsto A^\top AHBB^\top$. The second directional derivative is $\|AHB\|_F^2\ge0$, proving convexity without vectorizing the matrix.

The objective is strictly convex when the map $H\mapsto AHB$ has no nonzero kernel. Full column rank of $A$ and full row rank of $B$ ensure this. Rank deficiency can leave whole perturbation directions invisible to the objective.

For the separate expression $f(X)=\operatorname{tr}(X^\top AXB)$ with square compatible $A,B$, differentiating both copies of $X$ gives

$$
\nabla_Xf=AXB+A^\top XB^\top.
$$

Only when both fixed matrices are symmetric does this simplify to $2AXB$. Treating a quadratic-looking expression as though its coefficients were symmetric is a common transpose error.

#### Example LA-E23: a rectangular-objective derivative

Take $A=\operatorname{diag}(1,2)$, $B=\operatorname{diag}(1,3)$, and let $C$ be the two-by-two all-ones matrix. Minimize $\tfrac12\|AXB-C\|_F^2$.

Solution. Both fixed matrices are invertible, so exact fit is possible and unique:

$$
X^*=A^{-1}CB^{-1}
=\begin{pmatrix}1&1/3\\1/2&1/6\end{pmatrix}.
$$

Multiplying $AX^*B$ gives $C$. The residual and the Frobenius gradient are both zero. For any perturbation $H\ne0$, invertibility gives $AHB\ne0$ and the loss increases by $\tfrac12\|AHB\|_F^2$. This checks uniqueness through the objective itself.

For symmetric positive-definite $X$, another useful differential is $D\log\det X[H]=\operatorname{tr}(X^{-1}H)$. In a symmetric direction its second derivative is $-\operatorname{tr}(X^{-1}HX^{-1}H)$. Writing $K=X^{-1/2}HX^{-1/2}$ makes this $-\operatorname{tr}(K^2)\le0$, proving concavity of log determinant on the positive-definite cone. The domain is essential: an arbitrary real matrix need not have a positive determinant or a real logarithm.

### 6.3 Worked examples

#### Example LA-E16: differentiating a nonsymmetric expression

Let $A=\left(\begin{smallmatrix}1&3\\0&2\end{smallmatrix}\right)$.

Solution. Since $x^\top Ax=x_1^2+3x_1x_2+2x_2^2$, the gradient is $(2x_1+3x_2,3x_1+4x_2)^\top$. It equals $(A+A^\top)x$, not $2Ax$. The Hessian has determinant $8-9=-1$, so the quadratic is indefinite.

#### Example LA-E17: a budget-constrained minimum

Minimize $w^\top\Sigma w$ subject to $\mathbf1^\top w=1$, with $\Sigma$ PD.

Solution. The Lagrangian $L=w^\top\Sigma w-\eta(\mathbf1^\top w-1)$ gives $2\Sigma w=\eta\mathbf1$. Enforcing the budget gives

$$
w^*=\frac{\Sigma^{-1}\mathbf1}{\mathbf1^\top\Sigma^{-1}\mathbf1}.
$$

The denominator is positive because $\Sigma^{-1}$ is PD. Every other feasible $w=w^*+h$ has $\mathbf1^\top h=0$, and the objective difference is $h^\top\Sigma h>0$ unless $h=0$. With $\Sigma=\operatorname{diag}(1,4)$, the weights are $(4/5,1/5)$ and the minimum variance is $4/5$.

#### Example LA-E18: a matrix gradient

For $F(X)=\tfrac12\|AX-B\|_F^2$, find $\nabla_XF$.

Solution. Write $R=AX-B$. Then

$$
dF=\operatorname{tr}(R^\top dR)
=\operatorname{tr}(R^\top A\,dX)
=\operatorname{tr}([A^\top R]^\top dX).
$$

Hence $\nabla_XF=A^\top(AX-B)$. Each column is an ordinary least-squares gradient, consistent with the sum-of-column-errors interpretation.

### 6.4 Additional solved problems

These problems are constructed for these notes.

#### Problem LA-17: a singular quadratic objective

For symmetric PSD $H$, characterize the minimizers of $f(x)=\tfrac12x^\top Hx-g^\top x$.

Solution. Decompose $g=g_R+g_N$ into $\operatorname{range}H$ and $\ker H$. If $g_N\ne0$, set $x=tg_N$. The quadratic term is zero and $f(x)=-t\|g_N\|^2\to-\infty$, so no minimum exists. If $g_N=0$, then $g$ lies in the range and the equation $Hx=g$ is solvable. Its solutions are $H^+g+z$, $z\in\ker H$. Expanding around any such solution gives a nonnegative difference $\tfrac12h^\top Hh$. All these solutions are minimizers; uniqueness occurs exactly when $H$ is PD.

#### Problem LA-18: a singular covariance with a budget constraint

Set $\Sigma=\operatorname{diag}(1,0)$ and minimize $w^\top\Sigma w$ under $w_1+w_2=1$.

Solution. Substitution gives the objective $w_1^2$, uniquely minimized at $(0,1)$ with value zero. The inverse-based formula cannot be used because $\Sigma$ is singular. Replacing the inverse by $\Sigma^+$ in that formula would give $(1,0)$, which is wrong. A null direction satisfying the budget can dominate every positive-variance direction. In general solve the constrained problem using its null space or its KKT system, and check curvature along feasible directions.

#### Problem LA-19: an active nonnegativity constraint

Minimize $\tfrac12[(x_1+1)^2+(x_2-2)^2]$ subject to $x_1,x_2\ge0$.

Solution. The unconstrained minimizer $(-1,2)$ is infeasible. Use constraints $-x_i\le0$ and multipliers $\mu_i\ge0$. Stationarity reads $x_1+1-\mu_1=0$ and $x_2-2-\mu_2=0$. At $(0,2)$ choose $(\mu_1,\mu_2)=(1,0)$. Feasibility, stationarity, and $\mu_ix_i=0$ all hold. Strict convexity gives the unique global minimum, whose value is $1/2$.

#### Problem LA-20: ridge with an unpenalized intercept

Minimize $\|y-a\mathbf1-Xb\|^2+\lambda\|b\|^2$, with $\lambda>0$.

Solution. Differentiating in $a$ gives $a=\bar y-\bar x^\top b$, where $\bar x^\top=\mathbf1^\top X/n$. Substituting centers both sides:

$$
\min_b\|Hy-HXb\|^2+\lambda\|b\|^2,\qquad
b^*=(X^\top HX+\lambda I)^{-1}X^\top Hy.
$$

The intercept follows from the earlier equation. Since $X^\top HX$ is PSD, adding $\lambda I$ makes the slope solution unique. Penalizing the intercept as well would define a different problem and would no longer force the residuals to sum to zero.

### Further problems on matrix differentiation

#### Problem LA-61: weighted ridge regression

Constructed problem. Minimize $(y-Xb)^\top W(y-Xb)+\lambda b^\top Rb$, with $W$ PD, $R$ PSD, and $\lambda>0$.

Solution. Differentiate the two quadratic terms:

$$
\nabla f(b)=-2X^\top W(y-Xb)+2\lambda Rb.
$$

Stationarity gives $(X^\top WX+\lambda R)b=X^\top Wy$. For any $z$, the coefficient matrix has quadratic form $\|W^{1/2}Xz\|^2+\lambda z^\top Rz$. Its kernel is $\ker X\cap\ker R$. Thus the solution is unique exactly when this intersection is zero. Merely calling the penalty ridge is insufficient to guarantee uniqueness if $R$ leaves some unidentified directions unpenalized.

#### Problem LA-62: sensitivity of a quadratic optimum

Constructed problem. Let $x^*(g)=H^{-1}g$ for fixed PD $H$. Differentiate the optimum and the minimized value with respect to $g$.

Solution. A perturbation $dg$ changes the optimum by $dx^*=H^{-1}dg$. The minimized value is $v(g)=c-\tfrac12g^\top H^{-1}g$, so $\nabla_gv=-H^{-1}g=-x^*$. Its Hessian is $-H^{-1}$, hence the value is concave in $g$. The decision is most sensitive along the smallest-eigenvalue directions of $H$, agreeing with the conditioning calculation.

#### Problem LA-63: optimizing a trace expression

Constructed problem. For symmetric PD $H$ and a given $G\in\mathbb R^{n\times k}$, minimize $F(X)=\tfrac12\operatorname{tr}(X^\top HX)-\operatorname{tr}(G^\top X)$.

Solution. Write $X$ and $G$ in columns. The trace expression is $\sum_{j=1}^k[\tfrac12x_j^\top Hx_j-g_j^\top x_j]$. Each column has the unique minimum $x_j=H^{-1}g_j$. Therefore $X^*=H^{-1}G$. Equivalently the matrix gradient is $HX-G$, and

$$
F(X)-F(X^*)=\tfrac12\operatorname{tr}((X-X^*)^\top H(X-X^*))\ge0,
$$

with equality only at $X=X^*$.

#### Problem LA-64: a quadratic with two equality constraints

Constructed problem. Minimize $x_1^2+x_2^2+x_3^2$ subject to $x_1+x_2+x_3=1$ and $x_1-x_3=0$.

Solution. The second equality sets $x_3=x_1=t$, and the first gives $x_2=1-2t$. Substitute:

$$
f(t)=2t^2+(1-2t)^2=6t^2-4t+1.
$$

Its derivative is $12t-4$, so $t=1/3$ and $x=(1/3,1/3,1/3)$. The second derivative is $12>0$, proving uniqueness on the feasible line. The minimum is $1/3$. Eliminating constraints first has reduced a three-variable constrained problem to an ordinary scalar quadratic.

#### Problem LA-65: curvature only on a feasible line

Constructed problem. Minimize $x_1^2-x_2^2$ subject to $x_2=0$, then subject to $x_1=0$.

Solution. Under the first constraint the objective becomes $x_1^2$, so zero is the unique minimizer even though the full Hessian is indefinite. Under the second constraint it becomes $-x_2^2$, which tends to $-\infty$ as $|x_2|\to\infty$. The two feasible tangent spaces select opposite curvature directions. A constrained second-order test must inspect feasible directions, not merely classify the unconstrained Hessian.

#### Problem LA-66: the derivative of an inverse quadratic

Constructed problem. For fixed $a\ne0$ and PD $X$, differentiate $f(X)=a^\top X^{-1}a$.

Solution. The inverse differential gives

$$
df=-a^\top X^{-1}(dX)X^{-1}a
=-\operatorname{tr}(X^{-1}aa^\top X^{-1}dX).
$$

Thus $\nabla_Xf=-X^{-1}aa^\top X^{-1}$ on symmetric matrices. Along a PSD direction $H$, $df[H]=-(X^{-1}a)^\top H(X^{-1}a)\le0$. Increasing the matrix in the PSD order can only decrease this inverse quadratic, as long as the path stays PD.

#### Problem LA-67: a log-determinant optimum

Constructed problem. For PD $S$, minimize $F(X)=\operatorname{tr}(SX)-\log\det X$ over PD $X$.

Solution. Stationarity gives $S-X^{-1}=0$, so $X^*=S^{-1}$. To verify global optimality without assuming it, set $Y=S^{1/2}XS^{1/2}$. Then $\operatorname{tr}(SX)=\operatorname{tr}Y$ and $\log\det X=\log\det Y-\log\det S$. Hence, up to a constant,

$$
F(X)=\sum_i(\lambda_i(Y)-\log\lambda_i(Y)).
$$

For a positive scalar $t$, $t-\log t$ has unique minimum one at $t=1$. Thus all eigenvalues of $Y$ must equal one, giving $Y=I$ and $X=S^{-1}$. This proves existence and uniqueness.

#### Problem LA-68: differentiating a normalized quadratic solution

Constructed problem. Let $w(\Sigma)=\Sigma^{-1}\mathbf1/(\mathbf1^\top\Sigma^{-1}\mathbf1)$ for PD $\Sigma$. Find its first-order change under a symmetric perturbation $E$.

Solution. Write $a=\Sigma^{-1}\mathbf1$ and $s=\mathbf1^\top a>0$. The inverse formula gives $da=-\Sigma^{-1}Ea$ and $ds=-a^\top Ea$. Applying the scalar quotient rule componentwise,

$$
dw=-\frac{\Sigma^{-1}Ea}{s}
+\frac{a(a^\top Ea)}{s^2}.
$$

Multiplying by $\mathbf1^\top$ gives zero: the first term sums to $-a^\top Ea/s$, and the second to $a^\top Ea/s$. Thus the first-order perturbation respects the budget constraint. If $\Sigma$ is nearly singular, the inverse factors can make this change large even for small $E$.

<!-- extended-practice-start -->

## 7. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 7.1 Weighted projection onto a line

For an ordinary Euclidean projection, all coordinate errors receive the same weight. With a positive diagonal matrix $W$, the loss is $(y-\hat y)^\top W(y-\hat y)$. A large diagonal entry makes the corresponding error more expensive. The fitted vector still lies in the chosen subspace, but the residual is orthogonal using the $W$ inner product. These problems compute both objects and check that orthogonality.

#### Problem LA-69: weighted line fit, case 1

Let $x=\begin{pmatrix}1\\2\end{pmatrix}$, $y=\begin{pmatrix}4\\1\end{pmatrix}$, and $W=\begin{pmatrix}3&0\\0&1\end{pmatrix}$. Find the scalar $a_*$ minimizing $g(a)=(y-ax)^\top W(y-ax)$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

$$
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
$$

Here $x^\top Wx=7>0$ and $x^\top Wy=14$. The derivative is $g'(a)=-2x^\top Wy+2a x^\top Wx$, so

$$
a_*=\frac{x^\top Wy}{x^\top Wx}=2,\qquad
\hat y=a_*x=\begin{pmatrix}2\\4\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}2\\-3\end{pmatrix}.
$$

The second derivative is $2x^\top Wx=14>0$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $x$.

For an arbitrary input vector $z$, the same derivation gives $P_Wz=x(x^\top Wz)/(x^\top Wx)$. Therefore

$$
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{3}{7}&\frac{2}{7}\\\frac{6}{7}&\frac{4}{7}\end{pmatrix}.
$$

Multiplying this matrix by $y$ recovers the fitted vector above. Also $P_Wx=x$, so applying the projection twice changes nothing: $P_W^2=P_W$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $P_W^\top W=WP_W$.

The residual check is $x^\top Wr=x^\top Wy-a_*x^\top Wx=0$. The minimum loss is

$$
r^\top Wr=21.
$$

Finally, writing $a=a_*+h$ and expanding the loss shows $g(a)-g(a_*)=7h^2$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-70: weighted line fit, case 2

Let $x=\begin{pmatrix}2\\-1\end{pmatrix}$, $y=\begin{pmatrix}1\\3\end{pmatrix}$, and $W=\begin{pmatrix}1&0\\0&4\end{pmatrix}$. Find the scalar $a_*$ minimizing $g(a)=(y-ax)^\top W(y-ax)$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

$$
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
$$

Here $x^\top Wx=8>0$ and $x^\top Wy=-10$. The derivative is $g'(a)=-2x^\top Wy+2a x^\top Wx$, so

$$
a_*=\frac{x^\top Wy}{x^\top Wx}=-\frac{5}{4},\qquad
\hat y=a_*x=\begin{pmatrix}-\frac{5}{2}\\\frac{5}{4}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{7}{2}\\\frac{7}{4}\end{pmatrix}.
$$

The second derivative is $2x^\top Wx=16>0$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $x$.

For an arbitrary input vector $z$, the same derivation gives $P_Wz=x(x^\top Wz)/(x^\top Wx)$. Therefore

$$
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{1}{2}&-1\\-\frac{1}{4}&\frac{1}{2}\end{pmatrix}.
$$

Multiplying this matrix by $y$ recovers the fitted vector above. Also $P_Wx=x$, so applying the projection twice changes nothing: $P_W^2=P_W$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $P_W^\top W=WP_W$.

The residual check is $x^\top Wr=x^\top Wy-a_*x^\top Wx=0$. The minimum loss is

$$
r^\top Wr=\frac{49}{2}.
$$

Finally, writing $a=a_*+h$ and expanding the loss shows $g(a)-g(a_*)=8h^2$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-71: weighted line fit, case 3

Let $x=\begin{pmatrix}1\\1\end{pmatrix}$, $y=\begin{pmatrix}-1\\4\end{pmatrix}$, and $W=\begin{pmatrix}2&0\\0&5\end{pmatrix}$. Find the scalar $a_*$ minimizing $g(a)=(y-ax)^\top W(y-ax)$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

$$
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
$$

Here $x^\top Wx=7>0$ and $x^\top Wy=18$. The derivative is $g'(a)=-2x^\top Wy+2a x^\top Wx$, so

$$
a_*=\frac{x^\top Wy}{x^\top Wx}=\frac{18}{7},\qquad
\hat y=a_*x=\begin{pmatrix}\frac{18}{7}\\\frac{18}{7}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}-\frac{25}{7}\\\frac{10}{7}\end{pmatrix}.
$$

The second derivative is $2x^\top Wx=14>0$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $x$.

For an arbitrary input vector $z$, the same derivation gives $P_Wz=x(x^\top Wz)/(x^\top Wx)$. Therefore

$$
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{2}{7}&\frac{5}{7}\\\frac{2}{7}&\frac{5}{7}\end{pmatrix}.
$$

Multiplying this matrix by $y$ recovers the fitted vector above. Also $P_Wx=x$, so applying the projection twice changes nothing: $P_W^2=P_W$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $P_W^\top W=WP_W$.

The residual check is $x^\top Wr=x^\top Wy-a_*x^\top Wx=0$. The minimum loss is

$$
r^\top Wr=\frac{250}{7}.
$$

Finally, writing $a=a_*+h$ and expanding the loss shows $g(a)-g(a_*)=7h^2$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-72: weighted line fit, case 4

Let $x=\begin{pmatrix}3\\1\end{pmatrix}$, $y=\begin{pmatrix}2\\-2\end{pmatrix}$, and $W=\begin{pmatrix}1&0\\0&2\end{pmatrix}$. Find the scalar $a_*$ minimizing $g(a)=(y-ax)^\top W(y-ax)$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

$$
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
$$

Here $x^\top Wx=11>0$ and $x^\top Wy=2$. The derivative is $g'(a)=-2x^\top Wy+2a x^\top Wx$, so

$$
a_*=\frac{x^\top Wy}{x^\top Wx}=\frac{2}{11},\qquad
\hat y=a_*x=\begin{pmatrix}\frac{6}{11}\\\frac{2}{11}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{16}{11}\\-\frac{24}{11}\end{pmatrix}.
$$

The second derivative is $2x^\top Wx=22>0$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $x$.

For an arbitrary input vector $z$, the same derivation gives $P_Wz=x(x^\top Wz)/(x^\top Wx)$. Therefore

$$
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{9}{11}&\frac{6}{11}\\\frac{3}{11}&\frac{2}{11}\end{pmatrix}.
$$

Multiplying this matrix by $y$ recovers the fitted vector above. Also $P_Wx=x$, so applying the projection twice changes nothing: $P_W^2=P_W$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $P_W^\top W=WP_W$.

The residual check is $x^\top Wr=x^\top Wy-a_*x^\top Wx=0$. The minimum loss is

$$
r^\top Wr=\frac{128}{11}.
$$

Finally, writing $a=a_*+h$ and expanding the loss shows $g(a)-g(a_*)=11h^2$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-73: weighted line fit, case 5

Let $x=\begin{pmatrix}1\\-2\end{pmatrix}$, $y=\begin{pmatrix}3\\2\end{pmatrix}$, and $W=\begin{pmatrix}4&0\\0&1\end{pmatrix}$. Find the scalar $a_*$ minimizing $g(a)=(y-ax)^\top W(y-ax)$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

$$
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
$$

Here $x^\top Wx=8>0$ and $x^\top Wy=8$. The derivative is $g'(a)=-2x^\top Wy+2a x^\top Wx$, so

$$
a_*=\frac{x^\top Wy}{x^\top Wx}=1,\qquad
\hat y=a_*x=\begin{pmatrix}1\\-2\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}2\\4\end{pmatrix}.
$$

The second derivative is $2x^\top Wx=16>0$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $x$.

For an arbitrary input vector $z$, the same derivation gives $P_Wz=x(x^\top Wz)/(x^\top Wx)$. Therefore

$$
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{1}{2}&-\frac{1}{4}\\-1&\frac{1}{2}\end{pmatrix}.
$$

Multiplying this matrix by $y$ recovers the fitted vector above. Also $P_Wx=x$, so applying the projection twice changes nothing: $P_W^2=P_W$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $P_W^\top W=WP_W$.

The residual check is $x^\top Wr=x^\top Wy-a_*x^\top Wx=0$. The minimum loss is

$$
r^\top Wr=32.
$$

Finally, writing $a=a_*+h$ and expanding the loss shows $g(a)-g(a_*)=8h^2$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

#### Problem LA-74: weighted line fit, case 6

Let $x=\begin{pmatrix}2\\3\end{pmatrix}$, $y=\begin{pmatrix}-1\\5\end{pmatrix}$, and $W=\begin{pmatrix}3&0\\0&2\end{pmatrix}$. Find the scalar $a_*$ minimizing $g(a)=(y-ax)^\top W(y-ax)$, its fitted vector, the projection matrix, and the minimum loss.

Solution. First expand the scalar quadratic rather than trying to invert a rectangular vector:

$$
g(a)=y^\top Wy-2a\,x^\top Wy+a^2x^\top Wx.
$$

Here $x^\top Wx=30>0$ and $x^\top Wy=24$. The derivative is $g'(a)=-2x^\top Wy+2a x^\top Wx$, so

$$
a_*=\frac{x^\top Wy}{x^\top Wx}=\frac{4}{5},\qquad
\hat y=a_*x=\begin{pmatrix}\frac{8}{5}\\\frac{12}{5}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}-\frac{13}{5}\\\frac{13}{5}\end{pmatrix}.
$$

The second derivative is $2x^\top Wx=60>0$. Thus the stationary point is the unique minimum, even when the fitted coefficient is negative. There is no restriction requiring a projection to point along the positive direction of $x$.

For an arbitrary input vector $z$, the same derivation gives $P_Wz=x(x^\top Wz)/(x^\top Wx)$. Therefore

$$
P_W=\frac{xx^\top W}{x^\top Wx}=\begin{pmatrix}\frac{2}{5}&\frac{2}{5}\\\frac{3}{5}&\frac{3}{5}\end{pmatrix}.
$$

Multiplying this matrix by $y$ recovers the fitted vector above. Also $P_Wx=x$, so applying the projection twice changes nothing: $P_W^2=P_W$. The matrix need not be symmetric in ordinary coordinates; its weighted symmetry is $P_W^\top W=WP_W$.

The residual check is $x^\top Wr=x^\top Wy-a_*x^\top Wx=0$. The minimum loss is

$$
r^\top Wr=\frac{169}{5}.
$$

Finally, writing $a=a_*+h$ and expanding the loss shows $g(a)-g(a_*)=30h^2$. The cross term vanishes by weighted orthogonality. This gives a direct check of global optimality and measures the extra loss from moving the coefficient away from its optimum.

### 7.2 Least squares with a visible residual

At the design points $-1,0,1$, the intercept column and slope column are orthogonal. This makes the normal equations diagonal without changing the meaning of least squares. The middle observation is especially useful: a straight line can fit all three observations exactly only if the second difference $y_0-2y_1+y_2$ vanishes.

#### Problem LA-75: three-observation regression, case 1

Fit $y_i=a+bt_i+r_i$ at $t=(-1,0,1)$ to $y=\begin{pmatrix}1\\2\\5\end{pmatrix}$ by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

$$
X=\begin{pmatrix}1&-1\\1&0\\1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\0&2\end{pmatrix}.
$$

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $2X^\top X$. The normal equations consequently have one solution:

$$
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
$$

Substitution gives $a=\frac{8}{3}$ and $b=2$. Evaluating the fitted line at each of the original design points gives

$$
\hat y=\begin{pmatrix}\frac{2}{3}\\\frac{8}{3}\\\frac{14}{3}\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}\frac{1}{3}\\-\frac{2}{3}\\\frac{1}{3}\end{pmatrix}.
$$

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $-r_0+r_2=0$. These are the two entries of $X^\top r=0$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $d=y_0-2y_1+y_2=2$. Direct substitution into the fitted values gives

$$
r=\frac d6\begin{pmatrix}1\\-2\\1\end{pmatrix},\qquad
\operatorname{SSE}=\frac{d^2}6=\frac{2}{3}.
$$

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $d$.

For any proposed change $h$ in the two coefficients, orthogonality gives $\|y-X(\hat\beta+h)\|_2^2=\|r\|_2^2+\|Xh\|_2^2$. The second term is strictly positive for nonzero $h$ because $X$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-76: three-observation regression, case 2

Fit $y_i=a+bt_i+r_i$ at $t=(-1,0,1)$ to $y=\begin{pmatrix}4\\1\\2\end{pmatrix}$ by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

$$
X=\begin{pmatrix}1&-1\\1&0\\1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\0&2\end{pmatrix}.
$$

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $2X^\top X$. The normal equations consequently have one solution:

$$
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
$$

Substitution gives $a=\frac{7}{3}$ and $b=-1$. Evaluating the fitted line at each of the original design points gives

$$
\hat y=\begin{pmatrix}\frac{10}{3}\\\frac{7}{3}\\\frac{4}{3}\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}\frac{2}{3}\\-\frac{4}{3}\\\frac{2}{3}\end{pmatrix}.
$$

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $-r_0+r_2=0$. These are the two entries of $X^\top r=0$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $d=y_0-2y_1+y_2=4$. Direct substitution into the fitted values gives

$$
r=\frac d6\begin{pmatrix}1\\-2\\1\end{pmatrix},\qquad
\operatorname{SSE}=\frac{d^2}6=\frac{8}{3}.
$$

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $d$.

For any proposed change $h$ in the two coefficients, orthogonality gives $\|y-X(\hat\beta+h)\|_2^2=\|r\|_2^2+\|Xh\|_2^2$. The second term is strictly positive for nonzero $h$ because $X$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-77: three-observation regression, case 3

Fit $y_i=a+bt_i+r_i$ at $t=(-1,0,1)$ to $y=\begin{pmatrix}-2\\0\\4\end{pmatrix}$ by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

$$
X=\begin{pmatrix}1&-1\\1&0\\1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\0&2\end{pmatrix}.
$$

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $2X^\top X$. The normal equations consequently have one solution:

$$
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
$$

Substitution gives $a=\frac{2}{3}$ and $b=3$. Evaluating the fitted line at each of the original design points gives

$$
\hat y=\begin{pmatrix}-\frac{7}{3}\\\frac{2}{3}\\\frac{11}{3}\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}\frac{1}{3}\\-\frac{2}{3}\\\frac{1}{3}\end{pmatrix}.
$$

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $-r_0+r_2=0$. These are the two entries of $X^\top r=0$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $d=y_0-2y_1+y_2=2$. Direct substitution into the fitted values gives

$$
r=\frac d6\begin{pmatrix}1\\-2\\1\end{pmatrix},\qquad
\operatorname{SSE}=\frac{d^2}6=\frac{2}{3}.
$$

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $d$.

For any proposed change $h$ in the two coefficients, orthogonality gives $\|y-X(\hat\beta+h)\|_2^2=\|r\|_2^2+\|Xh\|_2^2$. The second term is strictly positive for nonzero $h$ because $X$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-78: three-observation regression, case 4

Fit $y_i=a+bt_i+r_i$ at $t=(-1,0,1)$ to $y=\begin{pmatrix}3\\3\\3\end{pmatrix}$ by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

$$
X=\begin{pmatrix}1&-1\\1&0\\1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\0&2\end{pmatrix}.
$$

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $2X^\top X$. The normal equations consequently have one solution:

$$
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
$$

Substitution gives $a=3$ and $b=0$. Evaluating the fitted line at each of the original design points gives

$$
\hat y=\begin{pmatrix}3\\3\\3\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}0\\0\\0\end{pmatrix}.
$$

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $-r_0+r_2=0$. These are the two entries of $X^\top r=0$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $d=y_0-2y_1+y_2=0$. Direct substitution into the fitted values gives

$$
r=\frac d6\begin{pmatrix}1\\-2\\1\end{pmatrix},\qquad
\operatorname{SSE}=\frac{d^2}6=0.
$$

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $d$.

For any proposed change $h$ in the two coefficients, orthogonality gives $\|y-X(\hat\beta+h)\|_2^2=\|r\|_2^2+\|Xh\|_2^2$. The second term is strictly positive for nonzero $h$ because $X$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-79: three-observation regression, case 5

Fit $y_i=a+bt_i+r_i$ at $t=(-1,0,1)$ to $y=\begin{pmatrix}2\\-1\\2\end{pmatrix}$ by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

$$
X=\begin{pmatrix}1&-1\\1&0\\1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\0&2\end{pmatrix}.
$$

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $2X^\top X$. The normal equations consequently have one solution:

$$
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
$$

Substitution gives $a=1$ and $b=0$. Evaluating the fitted line at each of the original design points gives

$$
\hat y=\begin{pmatrix}1\\1\\1\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}1\\-2\\1\end{pmatrix}.
$$

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $-r_0+r_2=0$. These are the two entries of $X^\top r=0$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $d=y_0-2y_1+y_2=6$. Direct substitution into the fitted values gives

$$
r=\frac d6\begin{pmatrix}1\\-2\\1\end{pmatrix},\qquad
\operatorname{SSE}=\frac{d^2}6=6.
$$

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $d$.

For any proposed change $h$ in the two coefficients, orthogonality gives $\|y-X(\hat\beta+h)\|_2^2=\|r\|_2^2+\|Xh\|_2^2$. The second term is strictly positive for nonzero $h$ because $X$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

#### Problem LA-80: three-observation regression, case 6

Fit $y_i=a+bt_i+r_i$ at $t=(-1,0,1)$ to $y=\begin{pmatrix}5\\0\\-3\end{pmatrix}$ by least squares. Find the residual and prove that the answer is unique.

Solution. The design matrix and its Gram matrix are

$$
X=\begin{pmatrix}1&-1\\1&0\\1&1\end{pmatrix},\qquad
X^\top X=\begin{pmatrix}3&0\\0&2\end{pmatrix}.
$$

Both diagonal entries are positive. The two columns are independent, and the quadratic loss has a positive-definite Hessian $2X^\top X$. The normal equations consequently have one solution:

$$
3a=y_0+y_1+y_2,\qquad 2b=-y_0+y_2.
$$

Substitution gives $a=\frac{2}{3}$ and $b=-4$. Evaluating the fitted line at each of the original design points gives

$$
\hat y=\begin{pmatrix}\frac{14}{3}\\\frac{2}{3}\\-\frac{10}{3}\end{pmatrix},\qquad r=y-\hat y=\begin{pmatrix}\frac{1}{3}\\-\frac{2}{3}\\\frac{1}{3}\end{pmatrix}.
$$

The residuals sum to zero because the model contains an intercept. Their weighted sum using the design points is also zero: $-r_0+r_2=0$. These are the two entries of $X^\top r=0$. They are useful arithmetic checks after solving the normal equations.

There is a more informative expression for the entire residual. Define the second difference $d=y_0-2y_1+y_2=2$. Direct substitution into the fitted values gives

$$
r=\frac d6\begin{pmatrix}1\\-2\\1\end{pmatrix},\qquad
\operatorname{SSE}=\frac{d^2}6=\frac{2}{3}.
$$

The residual direction is perpendicular to both design columns and spans the one-dimensional orthogonal complement of the fitted plane. Thus all the failure of a line to fit three observations is summarized by one number, $d$.

For any proposed change $h$ in the two coefficients, orthogonality gives $\|y-X(\hat\beta+h)\|_2^2=\|r\|_2^2+\|Xh\|_2^2$. The second term is strictly positive for nonzero $h$ because $X$ has full column rank. This proves uniqueness without relying only on a derivative calculation.

### 7.3 Common and opposite directions in a quadratic form

For a symmetric matrix with equal diagonal entries, the directions $(1,1)$ and $(1,-1)$ separate common movement from opposite movement. A covariance matrix must be nonnegative in both directions. A matrix can have positive diagonal entries and still fail this requirement.

#### Problem LA-81: equal-diagonal quadratic form, case 1

Diagonalize $A=\begin{pmatrix}5&3\\3&5\end{pmatrix}$. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

$$
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\-1\end{pmatrix}.
$$

Their inner product is zero and each has norm one. Multiplication gives $Aq_+=8q_+$ and $Aq_-=2q_-$. Hence, with $Q=(q_+,q_-)$,

$$
A=Q\begin{pmatrix}8&0\\0&2\end{pmatrix}Q^\top.
$$

Every vector has a unique expansion $x=u q_++w q_-$, where $u=(x_1+x_2)/\sqrt2$ and $w=(x_1-x_2)/\sqrt2$. The quadratic form becomes

$$
x^\top Ax=8u^2+(2)w^2.
$$

This expression proves that $A$ is positive definite. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

Both eigenvalues are positive, so the inverse exists and is $\begin{pmatrix}\frac{5}{16}&-\frac{3}{16}\\-\frac{3}{16}&\frac{5}{16}\end{pmatrix}$.

As independent checks, the trace is $10=8+(2)$ and the determinant is $16=(8)(2)$. For nonzero $x$, the Rayleigh quotient is a weighted average of the two eigenvalues:

$$
\frac{x^\top Ax}{x^\top x}
=\frac{8u^2+(2)w^2}{u^2+w^2}.
$$

It lies between $2$ and $8$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-82: equal-diagonal quadratic form, case 2

Diagonalize $A=\begin{pmatrix}2&2\\2&2\end{pmatrix}$. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

$$
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\-1\end{pmatrix}.
$$

Their inner product is zero and each has norm one. Multiplication gives $Aq_+=4q_+$ and $Aq_-=0q_-$. Hence, with $Q=(q_+,q_-)$,

$$
A=Q\begin{pmatrix}4&0\\0&0\end{pmatrix}Q^\top.
$$

Every vector has a unique expansion $x=u q_++w q_-$, where $u=(x_1+x_2)/\sqrt2$ and $w=(x_1-x_2)/\sqrt2$. The quadratic form becomes

$$
x^\top Ax=4u^2+(0)w^2.
$$

This expression proves that $A$ is positive semidefinite and singular. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

The zero eigenvalue gives a nonzero direction with zero quadratic cost. The matrix has no inverse; a covariance interpretation is possible, with an exact linear relation between the two centered coordinates.

As independent checks, the trace is $4=4+(0)$ and the determinant is $0=(4)(0)$. For nonzero $x$, the Rayleigh quotient is a weighted average of the two eigenvalues:

$$
\frac{x^\top Ax}{x^\top x}
=\frac{4u^2+(0)w^2}{u^2+w^2}.
$$

It lies between $0$ and $4$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-83: equal-diagonal quadratic form, case 3

Diagonalize $A=\begin{pmatrix}2&3\\3&2\end{pmatrix}$. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

$$
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\-1\end{pmatrix}.
$$

Their inner product is zero and each has norm one. Multiplication gives $Aq_+=5q_+$ and $Aq_-=-1q_-$. Hence, with $Q=(q_+,q_-)$,

$$
A=Q\begin{pmatrix}5&0\\0&-1\end{pmatrix}Q^\top.
$$

Every vector has a unique expansion $x=u q_++w q_-$, where $u=(x_1+x_2)/\sqrt2$ and $w=(x_1-x_2)/\sqrt2$. The quadratic form becomes

$$
x^\top Ax=5u^2+(-1)w^2.
$$

This expression proves that $A$ is indefinite. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

The negative eigenvalue gives a direction with negative quadratic cost. Such a matrix cannot be a covariance matrix, since a variance cannot be negative.

As independent checks, the trace is $4=5+(-1)$ and the determinant is $-5=(5)(-1)$. For nonzero $x$, the Rayleigh quotient is a weighted average of the two eigenvalues:

$$
\frac{x^\top Ax}{x^\top x}
=\frac{5u^2+(-1)w^2}{u^2+w^2}.
$$

It lies between $-1$ and $5$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-84: equal-diagonal quadratic form, case 4

Diagonalize $A=\begin{pmatrix}4&-3\\-3&4\end{pmatrix}$. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

$$
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\-1\end{pmatrix}.
$$

Their inner product is zero and each has norm one. Multiplication gives $Aq_+=1q_+$ and $Aq_-=7q_-$. Hence, with $Q=(q_+,q_-)$,

$$
A=Q\begin{pmatrix}1&0\\0&7\end{pmatrix}Q^\top.
$$

Every vector has a unique expansion $x=u q_++w q_-$, where $u=(x_1+x_2)/\sqrt2$ and $w=(x_1-x_2)/\sqrt2$. The quadratic form becomes

$$
x^\top Ax=1u^2+(7)w^2.
$$

This expression proves that $A$ is positive definite. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

Both eigenvalues are positive, so the inverse exists and is $\begin{pmatrix}\frac{4}{7}&\frac{3}{7}\\\frac{3}{7}&\frac{4}{7}\end{pmatrix}$.

As independent checks, the trace is $8=1+(7)$ and the determinant is $7=(1)(7)$. For nonzero $x$, the Rayleigh quotient is a weighted average of the two eigenvalues:

$$
\frac{x^\top Ax}{x^\top x}
=\frac{1u^2+(7)w^2}{u^2+w^2}.
$$

It lies between $1$ and $7$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-85: equal-diagonal quadratic form, case 5

Diagonalize $A=\begin{pmatrix}3&0\\0&3\end{pmatrix}$. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

$$
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\-1\end{pmatrix}.
$$

Their inner product is zero and each has norm one. Multiplication gives $Aq_+=3q_+$ and $Aq_-=3q_-$. Hence, with $Q=(q_+,q_-)$,

$$
A=Q\begin{pmatrix}3&0\\0&3\end{pmatrix}Q^\top.
$$

Every vector has a unique expansion $x=u q_++w q_-$, where $u=(x_1+x_2)/\sqrt2$ and $w=(x_1-x_2)/\sqrt2$. The quadratic form becomes

$$
x^\top Ax=3u^2+(3)w^2.
$$

This expression proves that $A$ is positive definite. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

Both eigenvalues are positive, so the inverse exists and is $\begin{pmatrix}\frac{1}{3}&0\\0&\frac{1}{3}\end{pmatrix}$.

As independent checks, the trace is $6=3+(3)$ and the determinant is $9=(3)(3)$. For nonzero $x$, the Rayleigh quotient is a weighted average of the two eigenvalues:

$$
\frac{x^\top Ax}{x^\top x}
=\frac{3u^2+(3)w^2}{u^2+w^2}.
$$

It lies between $3$ and $3$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

#### Problem LA-86: equal-diagonal quadratic form, case 6

Diagonalize $A=\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$. Classify its quadratic form and determine whether it could be a covariance matrix.

Solution. Normalize the two natural directions:

$$
q_+=\frac1{\sqrt2}\begin{pmatrix}1\\1\end{pmatrix},
\qquad
q_-=\frac1{\sqrt2}\begin{pmatrix}1\\-1\end{pmatrix}.
$$

Their inner product is zero and each has norm one. Multiplication gives $Aq_+=0q_+$ and $Aq_-=2q_-$. Hence, with $Q=(q_+,q_-)$,

$$
A=Q\begin{pmatrix}0&0\\0&2\end{pmatrix}Q^\top.
$$

Every vector has a unique expansion $x=u q_++w q_-$, where $u=(x_1+x_2)/\sqrt2$ and $w=(x_1-x_2)/\sqrt2$. The quadratic form becomes

$$
x^\top Ax=0u^2+(2)w^2.
$$

This expression proves that $A$ is positive semidefinite and singular. It also identifies the direction responsible for each part of the classification; checking only the two diagonal entries would miss this information.

The zero eigenvalue gives a nonzero direction with zero quadratic cost. The matrix has no inverse; a covariance interpretation is possible, with an exact linear relation between the two centered coordinates.

As independent checks, the trace is $2=0+(2)$ and the determinant is $0=(0)(2)$. For nonzero $x$, the Rayleigh quotient is a weighted average of the two eigenvalues:

$$
\frac{x^\top Ax}{x^\top x}
=\frac{0u^2+(2)w^2}{u^2+w^2}.
$$

It lies between $0$ and $2$, with equality at the corresponding eigenvector. This turns the eigenvalue calculation into a complete statement about the size and sign of the quadratic form in every direction.

### 7.4 Rank-one least squares and minimum norm

A rank-one matrix can produce only multiples of one vector. Fitting data first projects the target onto that output line. Recovering coefficients then leaves an entire null-space direction undetermined. The pseudoinverse chooses the coefficient vector perpendicular to that null space.

#### Problem LA-87: rank-one pseudoinverse, case 1

Let $a=\begin{pmatrix}1\\2\end{pmatrix}$, $b=\begin{pmatrix}1\\1\end{pmatrix}$, $A=ab^\top$, and $y=\begin{pmatrix}3\\1\end{pmatrix}$. Find every least-squares solution of $Ax\approx y$ and identify the one with smallest Euclidean norm.

Solution. Since $Ax=a(b^\top x)$, the scalar $b^\top x$ determines the entire output. The column space of $A=\begin{pmatrix}1&1\\2&2\end{pmatrix}$ is the line spanned by $a$. Projecting $y$ onto this line gives

$$
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}1\\2\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}2\\-1\end{pmatrix}.
$$

Here $a^\top a=5$, $b^\top b=2$, and $a^\top y=5$. The check $a^\top r=0$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $b^\top x=1$.

The vector parallel to $b$ satisfying that equation is

$$
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{1}{2}\\\frac{1}{2}\end{pmatrix}.
$$

Consequently the pseudoinverse is

$$
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{1}{10}&\frac{1}{5}\\\frac{1}{10}&\frac{1}{5}\end{pmatrix}.
$$

The null space consists of vectors perpendicular to $b$. In two dimensions it is spanned by $z=\begin{pmatrix}-1\\1\end{pmatrix}$. All least-squares solutions are thus $x=x_++tz$ for real $t$. They give exactly the same fitted vector because $Az=0$.

The minimum-norm claim requires one further check. Since $x_+$ is parallel to $b$ while $z$ is perpendicular to $b$, their inner product is zero. Hence

$$
\|x_++tz\|_2^2=\|x_+\|_2^2+t^2\|z\|_2^2.
$$

The minimum occurs uniquely at $t=0$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-88: rank-one pseudoinverse, case 2

Let $a=\begin{pmatrix}2\\-1\end{pmatrix}$, $b=\begin{pmatrix}1\\-2\end{pmatrix}$, $A=ab^\top$, and $y=\begin{pmatrix}1\\4\end{pmatrix}$. Find every least-squares solution of $Ax\approx y$ and identify the one with smallest Euclidean norm.

Solution. Since $Ax=a(b^\top x)$, the scalar $b^\top x$ determines the entire output. The column space of $A=\begin{pmatrix}2&-4\\-1&2\end{pmatrix}$ is the line spanned by $a$. Projecting $y$ onto this line gives

$$
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}-\frac{4}{5}\\\frac{2}{5}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{9}{5}\\\frac{18}{5}\end{pmatrix}.
$$

Here $a^\top a=5$, $b^\top b=5$, and $a^\top y=-2$. The check $a^\top r=0$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $b^\top x=-\frac{2}{5}$.

The vector parallel to $b$ satisfying that equation is

$$
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}-\frac{2}{25}\\\frac{4}{25}\end{pmatrix}.
$$

Consequently the pseudoinverse is

$$
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{2}{25}&-\frac{1}{25}\\-\frac{4}{25}&\frac{2}{25}\end{pmatrix}.
$$

The null space consists of vectors perpendicular to $b$. In two dimensions it is spanned by $z=\begin{pmatrix}2\\1\end{pmatrix}$. All least-squares solutions are thus $x=x_++tz$ for real $t$. They give exactly the same fitted vector because $Az=0$.

The minimum-norm claim requires one further check. Since $x_+$ is parallel to $b$ while $z$ is perpendicular to $b$, their inner product is zero. Hence

$$
\|x_++tz\|_2^2=\|x_+\|_2^2+t^2\|z\|_2^2.
$$

The minimum occurs uniquely at $t=0$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-89: rank-one pseudoinverse, case 3

Let $a=\begin{pmatrix}1\\3\end{pmatrix}$, $b=\begin{pmatrix}2\\1\end{pmatrix}$, $A=ab^\top$, and $y=\begin{pmatrix}-1\\2\end{pmatrix}$. Find every least-squares solution of $Ax\approx y$ and identify the one with smallest Euclidean norm.

Solution. Since $Ax=a(b^\top x)$, the scalar $b^\top x$ determines the entire output. The column space of $A=\begin{pmatrix}2&1\\6&3\end{pmatrix}$ is the line spanned by $a$. Projecting $y$ onto this line gives

$$
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}\frac{1}{2}\\\frac{3}{2}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}-\frac{3}{2}\\\frac{1}{2}\end{pmatrix}.
$$

Here $a^\top a=10$, $b^\top b=5$, and $a^\top y=5$. The check $a^\top r=0$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $b^\top x=\frac{1}{2}$.

The vector parallel to $b$ satisfying that equation is

$$
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{1}{5}\\\frac{1}{10}\end{pmatrix}.
$$

Consequently the pseudoinverse is

$$
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{1}{25}&\frac{3}{25}\\\frac{1}{50}&\frac{3}{50}\end{pmatrix}.
$$

The null space consists of vectors perpendicular to $b$. In two dimensions it is spanned by $z=\begin{pmatrix}-1\\2\end{pmatrix}$. All least-squares solutions are thus $x=x_++tz$ for real $t$. They give exactly the same fitted vector because $Az=0$.

The minimum-norm claim requires one further check. Since $x_+$ is parallel to $b$ while $z$ is perpendicular to $b$, their inner product is zero. Hence

$$
\|x_++tz\|_2^2=\|x_+\|_2^2+t^2\|z\|_2^2.
$$

The minimum occurs uniquely at $t=0$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-90: rank-one pseudoinverse, case 4

Let $a=\begin{pmatrix}3\\1\end{pmatrix}$, $b=\begin{pmatrix}1\\-1\end{pmatrix}$, $A=ab^\top$, and $y=\begin{pmatrix}2\\5\end{pmatrix}$. Find every least-squares solution of $Ax\approx y$ and identify the one with smallest Euclidean norm.

Solution. Since $Ax=a(b^\top x)$, the scalar $b^\top x$ determines the entire output. The column space of $A=\begin{pmatrix}3&-3\\1&-1\end{pmatrix}$ is the line spanned by $a$. Projecting $y$ onto this line gives

$$
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}\frac{33}{10}\\\frac{11}{10}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}-\frac{13}{10}\\\frac{39}{10}\end{pmatrix}.
$$

Here $a^\top a=10$, $b^\top b=2$, and $a^\top y=11$. The check $a^\top r=0$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $b^\top x=\frac{11}{10}$.

The vector parallel to $b$ satisfying that equation is

$$
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{11}{20}\\-\frac{11}{20}\end{pmatrix}.
$$

Consequently the pseudoinverse is

$$
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{3}{20}&\frac{1}{20}\\-\frac{3}{20}&-\frac{1}{20}\end{pmatrix}.
$$

The null space consists of vectors perpendicular to $b$. In two dimensions it is spanned by $z=\begin{pmatrix}1\\1\end{pmatrix}$. All least-squares solutions are thus $x=x_++tz$ for real $t$. They give exactly the same fitted vector because $Az=0$.

The minimum-norm claim requires one further check. Since $x_+$ is parallel to $b$ while $z$ is perpendicular to $b$, their inner product is zero. Hence

$$
\|x_++tz\|_2^2=\|x_+\|_2^2+t^2\|z\|_2^2.
$$

The minimum occurs uniquely at $t=0$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-91: rank-one pseudoinverse, case 5

Let $a=\begin{pmatrix}1\\-1\end{pmatrix}$, $b=\begin{pmatrix}3\\2\end{pmatrix}$, $A=ab^\top$, and $y=\begin{pmatrix}4\\1\end{pmatrix}$. Find every least-squares solution of $Ax\approx y$ and identify the one with smallest Euclidean norm.

Solution. Since $Ax=a(b^\top x)$, the scalar $b^\top x$ determines the entire output. The column space of $A=\begin{pmatrix}3&2\\-3&-2\end{pmatrix}$ is the line spanned by $a$. Projecting $y$ onto this line gives

$$
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}\frac{3}{2}\\-\frac{3}{2}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{5}{2}\\\frac{5}{2}\end{pmatrix}.
$$

Here $a^\top a=2$, $b^\top b=13$, and $a^\top y=3$. The check $a^\top r=0$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $b^\top x=\frac{3}{2}$.

The vector parallel to $b$ satisfying that equation is

$$
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{9}{26}\\\frac{3}{13}\end{pmatrix}.
$$

Consequently the pseudoinverse is

$$
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{3}{26}&-\frac{3}{26}\\\frac{1}{13}&-\frac{1}{13}\end{pmatrix}.
$$

The null space consists of vectors perpendicular to $b$. In two dimensions it is spanned by $z=\begin{pmatrix}-2\\3\end{pmatrix}$. All least-squares solutions are thus $x=x_++tz$ for real $t$. They give exactly the same fitted vector because $Az=0$.

The minimum-norm claim requires one further check. Since $x_+$ is parallel to $b$ while $z$ is perpendicular to $b$, their inner product is zero. Hence

$$
\|x_++tz\|_2^2=\|x_+\|_2^2+t^2\|z\|_2^2.
$$

The minimum occurs uniquely at $t=0$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

#### Problem LA-92: rank-one pseudoinverse, case 6

Let $a=\begin{pmatrix}2\\3\end{pmatrix}$, $b=\begin{pmatrix}2\\-3\end{pmatrix}$, $A=ab^\top$, and $y=\begin{pmatrix}1\\-2\end{pmatrix}$. Find every least-squares solution of $Ax\approx y$ and identify the one with smallest Euclidean norm.

Solution. Since $Ax=a(b^\top x)$, the scalar $b^\top x$ determines the entire output. The column space of $A=\begin{pmatrix}4&-6\\6&-9\end{pmatrix}$ is the line spanned by $a$. Projecting $y$ onto this line gives

$$
\hat y=a\frac{a^\top y}{a^\top a}=\begin{pmatrix}-\frac{8}{13}\\-\frac{12}{13}\end{pmatrix},\qquad
r=y-\hat y=\begin{pmatrix}\frac{21}{13}\\-\frac{14}{13}\end{pmatrix}.
$$

Here $a^\top a=13$, $b^\top b=13$, and $a^\top y=-4$. The check $a^\top r=0$ proves that no other output on the line can reduce the squared residual. Every least-squares coefficient vector must therefore satisfy the single equation $b^\top x=-\frac{4}{13}$.

The vector parallel to $b$ satisfying that equation is

$$
x_+=\frac{b\,a^\top y}{(a^\top a)(b^\top b)}=\begin{pmatrix}-\frac{8}{169}\\\frac{12}{169}\end{pmatrix}.
$$

Consequently the pseudoinverse is

$$
A^+=\frac{ba^\top}{(a^\top a)(b^\top b)}=\begin{pmatrix}\frac{4}{169}&\frac{6}{169}\\-\frac{6}{169}&-\frac{9}{169}\end{pmatrix}.
$$

The null space consists of vectors perpendicular to $b$. In two dimensions it is spanned by $z=\begin{pmatrix}3\\2\end{pmatrix}$. All least-squares solutions are thus $x=x_++tz$ for real $t$. They give exactly the same fitted vector because $Az=0$.

The minimum-norm claim requires one further check. Since $x_+$ is parallel to $b$ while $z$ is perpendicular to $b$, their inner product is zero. Hence

$$
\|x_++tz\|_2^2=\|x_+\|_2^2+t^2\|z\|_2^2.
$$

The minimum occurs uniquely at $t=0$. This explains both parts of the pseudoinverse solution: first minimize the output error, then minimize the coefficient norm among all vectors with that optimal output.

### 7.5 Curvature along a constraint

A quadratic form need not be positive definite on the whole plane to have a unique minimum on a line. Under $x_1+x_2=s$, the only permitted changes are multiples of $(1,-1)$. The relevant curvature is the quadratic form evaluated along that direction.

#### Problem LA-93: constrained diagonal quadratic, case 1

Minimize $f(x)=1x_1^2+(4)x_2^2$ subject to $x_1+x_2=1$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $x_2=1-x_1$. The objective on the feasible line is

$$
g(x_1)=1x_1^2+(4)(1-x_1)^2.
$$

Its derivative is $g'(x_1)=2(5)x_1-2(4)(1)$ and its second derivative is $2(5)>0$. Therefore

$$
x_1^*=\frac{h_2s}{h_1+h_2}=\frac{4}{5},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=\frac{1}{5}.
$$

These coordinates sum to $1$, as required. The objective value is $\frac{4}{5}$.

For a separate derivation, use the convention $L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)$. Stationarity gives $2h_1x_1=\lambda$ and $2h_2x_2=\lambda$. Together with the equality constraint, these equations yield the same coordinates and $\lambda=\frac{8}{5}$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $x^*+(t,-t)$. Expanding the objective and using stationarity cancels the linear terms:

$$
f(x^*+(t,-t))-f(x^*)=(1+4)t^2=5t^2.
$$

This difference is positive whenever $t\ne0$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-94: constrained diagonal quadratic, case 2

Minimize $f(x)=2x_1^2+(3)x_2^2$ subject to $x_1+x_2=2$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $x_2=2-x_1$. The objective on the feasible line is

$$
g(x_1)=2x_1^2+(3)(2-x_1)^2.
$$

Its derivative is $g'(x_1)=2(5)x_1-2(3)(2)$ and its second derivative is $2(5)>0$. Therefore

$$
x_1^*=\frac{h_2s}{h_1+h_2}=\frac{6}{5},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=\frac{4}{5}.
$$

These coordinates sum to $2$, as required. The objective value is $\frac{24}{5}$.

For a separate derivation, use the convention $L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)$. Stationarity gives $2h_1x_1=\lambda$ and $2h_2x_2=\lambda$. Together with the equality constraint, these equations yield the same coordinates and $\lambda=\frac{24}{5}$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $x^*+(t,-t)$. Expanding the objective and using stationarity cancels the linear terms:

$$
f(x^*+(t,-t))-f(x^*)=(2+3)t^2=5t^2.
$$

This difference is positive whenever $t\ne0$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-95: constrained diagonal quadratic, case 3

Minimize $f(x)=0x_1^2+(5)x_2^2$ subject to $x_1+x_2=1$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $x_2=1-x_1$. The objective on the feasible line is

$$
g(x_1)=0x_1^2+(5)(1-x_1)^2.
$$

Its derivative is $g'(x_1)=2(5)x_1-2(5)(1)$ and its second derivative is $2(5)>0$. Therefore

$$
x_1^*=\frac{h_2s}{h_1+h_2}=1,\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=0.
$$

These coordinates sum to $1$, as required. The objective value is $0$.

For a separate derivation, use the convention $L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)$. Stationarity gives $2h_1x_1=\lambda$ and $2h_2x_2=\lambda$. Together with the equality constraint, these equations yield the same coordinates and $\lambda=0$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $x^*+(t,-t)$. Expanding the objective and using stationarity cancels the linear terms:

$$
f(x^*+(t,-t))-f(x^*)=(0+5)t^2=5t^2.
$$

This difference is positive whenever $t\ne0$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-96: constrained diagonal quadratic, case 4

Minimize $f(x)=-1x_1^2+(3)x_2^2$ subject to $x_1+x_2=1$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $x_2=1-x_1$. The objective on the feasible line is

$$
g(x_1)=-1x_1^2+(3)(1-x_1)^2.
$$

Its derivative is $g'(x_1)=2(2)x_1-2(3)(1)$ and its second derivative is $2(2)>0$. Therefore

$$
x_1^*=\frac{h_2s}{h_1+h_2}=\frac{3}{2},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=-\frac{1}{2}.
$$

These coordinates sum to $1$, as required. The objective value is $-\frac{3}{2}$.

For a separate derivation, use the convention $L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)$. Stationarity gives $2h_1x_1=\lambda$ and $2h_2x_2=\lambda$. Together with the equality constraint, these equations yield the same coordinates and $\lambda=-3$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $x^*+(t,-t)$. Expanding the objective and using stationarity cancels the linear terms:

$$
f(x^*+(t,-t))-f(x^*)=(-1+3)t^2=2t^2.
$$

This difference is positive whenever $t\ne0$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-97: constrained diagonal quadratic, case 5

Minimize $f(x)=4x_1^2+(1)x_2^2$ subject to $x_1+x_2=-2$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $x_2=-2-x_1$. The objective on the feasible line is

$$
g(x_1)=4x_1^2+(1)(-2-x_1)^2.
$$

Its derivative is $g'(x_1)=2(5)x_1-2(1)(-2)$ and its second derivative is $2(5)>0$. Therefore

$$
x_1^*=\frac{h_2s}{h_1+h_2}=-\frac{2}{5},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=-\frac{8}{5}.
$$

These coordinates sum to $-2$, as required. The objective value is $\frac{16}{5}$.

For a separate derivation, use the convention $L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)$. Stationarity gives $2h_1x_1=\lambda$ and $2h_2x_2=\lambda$. Together with the equality constraint, these equations yield the same coordinates and $\lambda=-\frac{16}{5}$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $x^*+(t,-t)$. Expanding the objective and using stationarity cancels the linear terms:

$$
f(x^*+(t,-t))-f(x^*)=(4+1)t^2=5t^2.
$$

This difference is positive whenever $t\ne0$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

#### Problem LA-98: constrained diagonal quadratic, case 6

Minimize $f(x)=3x_1^2+(3)x_2^2$ subject to $x_1+x_2=1$. Derive the solution by elimination and by a multiplier.

Solution. Eliminate the equality by writing $x_2=1-x_1$. The objective on the feasible line is

$$
g(x_1)=3x_1^2+(3)(1-x_1)^2.
$$

Its derivative is $g'(x_1)=2(6)x_1-2(3)(1)$ and its second derivative is $2(6)>0$. Therefore

$$
x_1^*=\frac{h_2s}{h_1+h_2}=\frac{1}{2},\qquad
x_2^*=\frac{h_1s}{h_1+h_2}=\frac{1}{2}.
$$

These coordinates sum to $1$, as required. The objective value is $\frac{3}{2}$.

For a separate derivation, use the convention $L(x,\lambda)=h_1x_1^2+h_2x_2^2-\lambda(x_1+x_2-s)$. Stationarity gives $2h_1x_1=\lambda$ and $2h_2x_2=\lambda$. Together with the equality constraint, these equations yield the same coordinates and $\lambda=3$. Keeping the sign convention visible prevents a sign error in the multiplier.

To prove global optimality, consider any other feasible point $x^*+(t,-t)$. Expanding the objective and using stationarity cancels the linear terms:

$$
f(x^*+(t,-t))-f(x^*)=(3+3)t^2=6t^2.
$$

This difference is positive whenever $t\ne0$. It proves uniqueness on the feasible line, including cases where one diagonal coefficient is zero or negative. The unrestricted Hessian alone would not establish that conclusion. In a general equality-constrained quadratic problem, the analogous test restricts the Hessian to the null space of the constraint matrix.

### 7.6 A small residual in a weak direction

A residual measures how well a computed vector satisfies the supplied equations. A forward error measures its distance from the true solution. If the matrix almost removes a direction, a large forward error in that direction can produce a small residual. These cases keep the distinction visible using exact arithmetic.

#### Problem LA-99: weak-direction sensitivity, case 1

Let $A=\operatorname{diag}(1,\frac{1}{10})$ and $b=(1,\frac{1}{10})^\top$. Replace the second component of $b$ by $\frac{11}{100}$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $x_1=1$ and $\frac{1}{10}x_2=\frac{1}{10}$, so $x=(1,1)^\top$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

$$
\tilde x=\begin{pmatrix}1\\\frac{11}{10}\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\\frac{1}{10}\end{pmatrix}.
$$

The input perturbation has norm $|\delta|=\frac{1}{100}$, while the solution perturbation has norm $|\delta|/\varepsilon=\frac{1}{10}$. Division by $\varepsilon$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $r=b-A\tilde x$,

$$
r=\begin{pmatrix}0\\-\frac{1}{100}\end{pmatrix},\qquad
\|r\|_2=\frac{1}{100}.
$$

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $A$ are $1$ and $\varepsilon=\frac{1}{10}$. Therefore

$$
\|A\|_2=1,\qquad \|A^{-1}\|_2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=10.
$$

The general absolute bound $\|\Delta x\|_2\le\|A^{-1}\|_2\|\Delta b\|_2$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-100: weak-direction sensitivity, case 2

Let $A=\operatorname{diag}(1,\frac{1}{100})$ and $b=(1,\frac{1}{100})^\top$. Replace the second component of $b$ by $\frac{11}{1000}$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $x_1=1$ and $\frac{1}{100}x_2=\frac{1}{100}$, so $x=(1,1)^\top$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

$$
\tilde x=\begin{pmatrix}1\\\frac{11}{10}\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\\frac{1}{10}\end{pmatrix}.
$$

The input perturbation has norm $|\delta|=\frac{1}{1000}$, while the solution perturbation has norm $|\delta|/\varepsilon=\frac{1}{10}$. Division by $\varepsilon$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $r=b-A\tilde x$,

$$
r=\begin{pmatrix}0\\-\frac{1}{1000}\end{pmatrix},\qquad
\|r\|_2=\frac{1}{1000}.
$$

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $A$ are $1$ and $\varepsilon=\frac{1}{100}$. Therefore

$$
\|A\|_2=1,\qquad \|A^{-1}\|_2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=100.
$$

The general absolute bound $\|\Delta x\|_2\le\|A^{-1}\|_2\|\Delta b\|_2$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-101: weak-direction sensitivity, case 3

Let $A=\operatorname{diag}(1,\frac{1}{1000})$ and $b=(1,\frac{1}{1000})^\top$. Replace the second component of $b$ by $\frac{1}{500}$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $x_1=1$ and $\frac{1}{1000}x_2=\frac{1}{1000}$, so $x=(1,1)^\top$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

$$
\tilde x=\begin{pmatrix}1\\2\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\1\end{pmatrix}.
$$

The input perturbation has norm $|\delta|=\frac{1}{1000}$, while the solution perturbation has norm $|\delta|/\varepsilon=1$. Division by $\varepsilon$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $r=b-A\tilde x$,

$$
r=\begin{pmatrix}0\\-\frac{1}{1000}\end{pmatrix},\qquad
\|r\|_2=\frac{1}{1000}.
$$

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $A$ are $1$ and $\varepsilon=\frac{1}{1000}$. Therefore

$$
\|A\|_2=1,\qquad \|A^{-1}\|_2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=1000.
$$

The general absolute bound $\|\Delta x\|_2\le\|A^{-1}\|_2\|\Delta b\|_2$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-102: weak-direction sensitivity, case 4

Let $A=\operatorname{diag}(1,\frac{1}{10000})$ and $b=(1,\frac{1}{10000})^\top$. Replace the second component of $b$ by $\frac{101}{1000000}$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $x_1=1$ and $\frac{1}{10000}x_2=\frac{1}{10000}$, so $x=(1,1)^\top$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

$$
\tilde x=\begin{pmatrix}1\\\frac{101}{100}\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\\frac{1}{100}\end{pmatrix}.
$$

The input perturbation has norm $|\delta|=\frac{1}{1000000}$, while the solution perturbation has norm $|\delta|/\varepsilon=\frac{1}{100}$. Division by $\varepsilon$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $r=b-A\tilde x$,

$$
r=\begin{pmatrix}0\\-\frac{1}{1000000}\end{pmatrix},\qquad
\|r\|_2=\frac{1}{1000000}.
$$

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $A$ are $1$ and $\varepsilon=\frac{1}{10000}$. Therefore

$$
\|A\|_2=1,\qquad \|A^{-1}\|_2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=10000.
$$

The general absolute bound $\|\Delta x\|_2\le\|A^{-1}\|_2\|\Delta b\|_2$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-103: weak-direction sensitivity, case 5

Let $A=\operatorname{diag}(1,\frac{1}{1000000})$ and $b=(1,\frac{1}{1000000})^\top$. Replace the second component of $b$ by $\frac{1}{500000}$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $x_1=1$ and $\frac{1}{1000000}x_2=\frac{1}{1000000}$, so $x=(1,1)^\top$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

$$
\tilde x=\begin{pmatrix}1\\2\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\1\end{pmatrix}.
$$

The input perturbation has norm $|\delta|=\frac{1}{1000000}$, while the solution perturbation has norm $|\delta|/\varepsilon=1$. Division by $\varepsilon$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $r=b-A\tilde x$,

$$
r=\begin{pmatrix}0\\-\frac{1}{1000000}\end{pmatrix},\qquad
\|r\|_2=\frac{1}{1000000}.
$$

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $A$ are $1$ and $\varepsilon=\frac{1}{1000000}$. Therefore

$$
\|A\|_2=1,\qquad \|A^{-1}\|_2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=1000000.
$$

The general absolute bound $\|\Delta x\|_2\le\|A^{-1}\|_2\|\Delta b\|_2$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.

#### Problem LA-104: weak-direction sensitivity, case 6

Let $A=\operatorname{diag}(1,\frac{1}{1000000})$ and $b=(1,\frac{1}{1000000})^\top$. Replace the second component of $b$ by $-\frac{9}{1000000}$. Compute the old and new solutions, the original-system residual, and the condition number.

Solution. The original equations read $x_1=1$ and $\frac{1}{1000000}x_2=\frac{1}{1000000}$, so $x=(1,1)^\top$. Only the second right-hand-side component changes. Dividing that changed equation by the small diagonal entry gives

$$
\tilde x=\begin{pmatrix}1\\-9\end{pmatrix},\qquad
\Delta x=\begin{pmatrix}0\\-10\end{pmatrix}.
$$

The input perturbation has norm $|\delta|=\frac{1}{100000}$, while the solution perturbation has norm $|\delta|/\varepsilon=10$. Division by $\varepsilon$ is the source of the amplification.

Now evaluate the new solution in the original equations. With residual convention $r=b-A\tilde x$,

$$
r=\begin{pmatrix}0\\\frac{1}{100000}\end{pmatrix},\qquad
\|r\|_2=\frac{1}{100000}.
$$

The residual can be tiny even when the second solution coordinate changes substantially. This is an exact example, so the effect does not depend on a particular floating-point algorithm.

The singular values of $A$ are $1$ and $\varepsilon=\frac{1}{1000000}$. Therefore

$$
\|A\|_2=1,\qquad \|A^{-1}\|_2=\varepsilon^{-1},\qquad
\kappa_2(A)=\varepsilon^{-1}=1000000.
$$

The general absolute bound $\|\Delta x\|_2\le\|A^{-1}\|_2\|\Delta b\|_2$ is an equality here because the perturbation points exactly along the weakest singular direction. A perturbation in the first coordinate would have amplification one. Conditioning gives a worst-direction bound; it does not say that every error is amplified equally.
