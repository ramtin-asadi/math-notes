# Analysis and Optimization — Part 2

The starting point is the real number system. Later sections use Euclidean spaces with $`\Vert x\Vert =(\sum_i x_i^2)^{1/2}`$. A statement about all points of a set includes its boundary unless an interior is explicitly specified.

Sources:

- Jiří Lebl, [Basic Analysis I](https://www.jirka.org/ra/realanal.pdf) and [Basic Analysis II](https://www.jirka.org/ra/realanal2.pdf), version 6.3.
- [MIT 18.100A, Real Analysis, Fall 2020](https://ocw.mit.edu/courses/18-100a-real-analysis-fall-2020/), taught by Casey Rodriguez.
- Stephen Boyd and Lieven Vandenberghe, [Convex Optimization](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf), Chapters 2–5 and 9.

Source exercise statements are summarized; the linked texts contain the originals. Examples and the additional problems in Section 8 are constructed here.

[Contents](README.md) · [Part 1](02-analysis-and-optimization-part-1.md)

Contents:

- [9. Further worked problems](#9-further-worked-problems)

## 9. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 9.1 Finding an explicit convergence index

An epsilon proof must produce a single index that works for every later term. Computing a limit by dividing numerator and denominator by $`n`$ suggests the answer, but does not yet provide that index. The following calculations also show how the sign of an error determines monotonicity.

#### Problem AN-25: rational sequence and epsilon index, case 1

For $`x_n=(n+2)/(n+5)`$, $`n\ge1`$, prove convergence, determine monotonicity, and give an index for error less than $`\frac{1}{100}`$. Give a direct Cauchy bound as well.

Solution. Subtract the proposed limit before taking absolute values:

```math
x_n-1=\frac{-3}{n+5},\qquad
\vert x_n-1\vert =\frac{3}{n+5}.
```

Since the denominator grows without bound and is positive for every allowed index, the limit is one. To obtain a quantitative proof, the inequality $`\vert x_n-1\vert \lt \epsilon`$ is equivalent to $`n\gt 3/\epsilon-5`$. At the requested tolerance, one valid index is $`N=296`$. For every integer $`n\ge N`$ the strict inequality holds. Using a floor followed by one avoids an equality at the boundary.

For arbitrary $`\epsilon\gt 0`$ the same construction, with a maximum against one, proves the full definition of convergence. The chosen index may depend on $`\epsilon`$; it cannot depend on which later term is being checked.

Consecutive differences are

```math
x_{n+1}-x_n
=\frac{3}{(n+5)(n+1+5)}.
```

Thus the sequence is strictly increasing. Its terms remain on the lower side of one. This is consistent with convergence and gives a second route using monotone convergence.

For a direct Cauchy estimate, let $`m,n\ge M`$. The triangle inequality gives

```math
\vert x_n-x_m\vert
\le \vert x_n-1\vert +\vert x_m-1\vert
\le\frac{6}{M+5}.
```

Taking $`M=596`$ makes this less than the requested tolerance. This bound is sufficient even though it is not the smallest possible one. An epsilon proof needs a valid uniform bound; it does not need the optimal index.

#### Problem AN-26: rational sequence and epsilon index, case 2

For $`x_n=(n+7)/(n+1)`$, $`n\ge1`$, prove convergence, determine monotonicity, and give an index for error less than $`\frac{1}{50}`$. Give a direct Cauchy bound as well.

Solution. Subtract the proposed limit before taking absolute values:

```math
x_n-1=\frac{6}{n+1},\qquad
\vert x_n-1\vert =\frac{6}{n+1}.
```

Since the denominator grows without bound and is positive for every allowed index, the limit is one. To obtain a quantitative proof, the inequality $`\vert x_n-1\vert \lt \epsilon`$ is equivalent to $`n\gt 6/\epsilon-1`$. At the requested tolerance, one valid index is $`N=300`$. For every integer $`n\ge N`$ the strict inequality holds. Using a floor followed by one avoids an equality at the boundary.

For arbitrary $`\epsilon\gt 0`$ the same construction, with a maximum against one, proves the full definition of convergence. The chosen index may depend on $`\epsilon`$; it cannot depend on which later term is being checked.

Consecutive differences are

```math
x_{n+1}-x_n
=\frac{-6}{(n+1)(n+1+1)}.
```

Thus the sequence is strictly decreasing. Its terms remain on the upper side of one. This is consistent with convergence and gives a second route using monotone convergence.

For a direct Cauchy estimate, let $`m,n\ge M`$. The triangle inequality gives

```math
\vert x_n-x_m\vert
\le \vert x_n-1\vert +\vert x_m-1\vert
\le\frac{12}{M+1}.
```

Taking $`M=600`$ makes this less than the requested tolerance. This bound is sufficient even though it is not the smallest possible one. An epsilon proof needs a valid uniform bound; it does not need the optimal index.

#### Problem AN-27: rational sequence and epsilon index, case 3

For $`x_n=(n+-3)/(n+4)`$, $`n\ge1`$, prove convergence, determine monotonicity, and give an index for error less than $`\frac{1}{1000}`$. Give a direct Cauchy bound as well.

Solution. Subtract the proposed limit before taking absolute values:

```math
x_n-1=\frac{-7}{n+4},\qquad
\vert x_n-1\vert =\frac{7}{n+4}.
```

Since the denominator grows without bound and is positive for every allowed index, the limit is one. To obtain a quantitative proof, the inequality $`\vert x_n-1\vert \lt \epsilon`$ is equivalent to $`n\gt 7/\epsilon-4`$. At the requested tolerance, one valid index is $`N=6997`$. For every integer $`n\ge N`$ the strict inequality holds. Using a floor followed by one avoids an equality at the boundary.

For arbitrary $`\epsilon\gt 0`$ the same construction, with a maximum against one, proves the full definition of convergence. The chosen index may depend on $`\epsilon`$; it cannot depend on which later term is being checked.

Consecutive differences are

```math
x_{n+1}-x_n
=\frac{7}{(n+4)(n+1+4)}.
```

Thus the sequence is strictly increasing. Its terms remain on the lower side of one. This is consistent with convergence and gives a second route using monotone convergence.

For a direct Cauchy estimate, let $`m,n\ge M`$. The triangle inequality gives

```math
\vert x_n-x_m\vert
\le \vert x_n-1\vert +\vert x_m-1\vert
\le\frac{14}{M+4}.
```

Taking $`M=13997`$ makes this less than the requested tolerance. This bound is sufficient even though it is not the smallest possible one. An epsilon proof needs a valid uniform bound; it does not need the optimal index.

#### Problem AN-28: rational sequence and epsilon index, case 4

For $`x_n=(n+1)/(n+3)`$, $`n\ge1`$, prove convergence, determine monotonicity, and give an index for error less than $`\frac{1}{20}`$. Give a direct Cauchy bound as well.

Solution. Subtract the proposed limit before taking absolute values:

```math
x_n-1=\frac{-2}{n+3},\qquad
\vert x_n-1\vert =\frac{2}{n+3}.
```

Since the denominator grows without bound and is positive for every allowed index, the limit is one. To obtain a quantitative proof, the inequality $`\vert x_n-1\vert \lt \epsilon`$ is equivalent to $`n\gt 2/\epsilon-3`$. At the requested tolerance, one valid index is $`N=38`$. For every integer $`n\ge N`$ the strict inequality holds. Using a floor followed by one avoids an equality at the boundary.

For arbitrary $`\epsilon\gt 0`$ the same construction, with a maximum against one, proves the full definition of convergence. The chosen index may depend on $`\epsilon`$; it cannot depend on which later term is being checked.

Consecutive differences are

```math
x_{n+1}-x_n
=\frac{2}{(n+3)(n+1+3)}.
```

Thus the sequence is strictly increasing. Its terms remain on the lower side of one. This is consistent with convergence and gives a second route using monotone convergence.

For a direct Cauchy estimate, let $`m,n\ge M`$. The triangle inequality gives

```math
\vert x_n-x_m\vert
\le \vert x_n-1\vert +\vert x_m-1\vert
\le\frac{4}{M+3}.
```

Taking $`M=78`$ makes this less than the requested tolerance. This bound is sufficient even though it is not the smallest possible one. An epsilon proof needs a valid uniform bound; it does not need the optimal index.

### 9.2 A usable delta for continuity

In a continuity proof, an expression involving both $`x`$ and the base point must be bounded using information fixed in advance. The common first restriction $`\vert x-a\vert \lt 1`$ creates such a local bound. A second restriction then controls the desired error. Both restrictions are enforced by taking a minimum.

#### Problem AN-29: epsilon-delta proof for a square, case 1

Prove continuity of $`f(x)=x^2`$ at $`a=2`$ and exhibit a delta for $`\epsilon=\frac{1}{10}`$.

Solution. Factor the difference:

```math
\vert x^2-a^2\vert =\vert x-a\vert \vert x+a\vert .
```

The second factor still depends on $`x`$, so choosing $`\delta=\epsilon/\vert x+a\vert`$ would be circular. Instead first require $`\vert x-a\vert \lt 1`$. Since $`x+a=(x-a)+2a`$, this yields $`\vert x+a\vert \le\vert x-a\vert +2\vert a\vert \lt 1+2\vert a\vert =5`$.

Now choose

```math
\delta=\min\left(1,\frac{\epsilon}{5}\right).
```

Whenever $`\vert x-a\vert \lt \delta`$, both required bounds hold, and therefore

```math
\vert f(x)-f(a)\vert \lt  5\delta\le\epsilon.
```

For the supplied tolerance this choice is $`\delta=\frac{1}{50}`$. The base-point value is $`f(a)=4`$. Every point within the chosen delta has its squared value within epsilon of that number.

The order of the choices is part of the proof. The base point is fixed first, epsilon is supplied next, and delta is then chosen using only those two quantities. The point $`x`$ is allowed to vary afterwards.

This argument also gives a local Lipschitz bound $`\vert x^2-a^2\vert \le 5\vert x-a\vert`$ within distance one of this base point. It does not produce a global Lipschitz constant on the whole real line. Indeed, comparing $`x`$ and $`x+1`$ gives a difference $`\vert 2x+1\vert`$, which has no finite global upper bound. Continuity at every point and one global Lipschitz estimate are different assertions.

#### Problem AN-30: epsilon-delta proof for a square, case 2

Prove continuity of $`f(x)=x^2`$ at $`a=-3`$ and exhibit a delta for $`\epsilon=\frac{1}{100}`$.

Solution. Factor the difference:

```math
\vert x^2-a^2\vert =\vert x-a\vert \vert x+a\vert .
```

The second factor still depends on $`x`$, so choosing $`\delta=\epsilon/\vert x+a\vert`$ would be circular. Instead first require $`\vert x-a\vert \lt 1`$. Since $`x+a=(x-a)+2a`$, this yields $`\vert x+a\vert \le\vert x-a\vert +2\vert a\vert \lt 1+2\vert a\vert =7`$.

Now choose

```math
\delta=\min\left(1,\frac{\epsilon}{7}\right).
```

Whenever $`\vert x-a\vert \lt \delta`$, both required bounds hold, and therefore

```math
\vert f(x)-f(a)\vert \lt  7\delta\le\epsilon.
```

For the supplied tolerance this choice is $`\delta=\frac{1}{700}`$. The base-point value is $`f(a)=9`$. Every point within the chosen delta has its squared value within epsilon of that number.

The order of the choices is part of the proof. The base point is fixed first, epsilon is supplied next, and delta is then chosen using only those two quantities. The point $`x`$ is allowed to vary afterwards.

This argument also gives a local Lipschitz bound $`\vert x^2-a^2\vert \le 7\vert x-a\vert`$ within distance one of this base point. It does not produce a global Lipschitz constant on the whole real line. Indeed, comparing $`x`$ and $`x+1`$ gives a difference $`\vert 2x+1\vert`$, which has no finite global upper bound. Continuity at every point and one global Lipschitz estimate are different assertions.

#### Problem AN-31: epsilon-delta proof for a square, case 3

Prove continuity of $`f(x)=x^2`$ at $`a=0`$ and exhibit a delta for $`\epsilon=\frac{1}{1000}`$.

Solution. Factor the difference:

```math
\vert x^2-a^2\vert =\vert x-a\vert \vert x+a\vert .
```

The second factor still depends on $`x`$, so choosing $`\delta=\epsilon/\vert x+a\vert`$ would be circular. Instead first require $`\vert x-a\vert \lt 1`$. Since $`x+a=(x-a)+2a`$, this yields $`\vert x+a\vert \le\vert x-a\vert +2\vert a\vert \lt 1+2\vert a\vert =1`$.

Now choose

```math
\delta=\min\left(1,\frac{\epsilon}{1}\right).
```

Whenever $`\vert x-a\vert \lt \delta`$, both required bounds hold, and therefore

```math
\vert f(x)-f(a)\vert \lt  1\delta\le\epsilon.
```

For the supplied tolerance this choice is $`\delta=\frac{1}{1000}`$. The base-point value is $`f(a)=0`$. Every point within the chosen delta has its squared value within epsilon of that number.

The order of the choices is part of the proof. The base point is fixed first, epsilon is supplied next, and delta is then chosen using only those two quantities. The point $`x`$ is allowed to vary afterwards.

This argument also gives a local Lipschitz bound $`\vert x^2-a^2\vert \le 1\vert x-a\vert`$ within distance one of this base point. It does not produce a global Lipschitz constant on the whole real line. Indeed, comparing $`x`$ and $`x+1`$ gives a difference $`\vert 2x+1\vert`$, which has no finite global upper bound. Continuity at every point and one global Lipschitz estimate are different assertions.

#### Problem AN-32: epsilon-delta proof for a square, case 4

Prove continuity of $`f(x)=x^2`$ at $`a=5`$ and exhibit a delta for $`\epsilon=\frac{1}{20}`$.

Solution. Factor the difference:

```math
\vert x^2-a^2\vert =\vert x-a\vert \vert x+a\vert .
```

The second factor still depends on $`x`$, so choosing $`\delta=\epsilon/\vert x+a\vert`$ would be circular. Instead first require $`\vert x-a\vert \lt 1`$. Since $`x+a=(x-a)+2a`$, this yields $`\vert x+a\vert \le\vert x-a\vert +2\vert a\vert \lt 1+2\vert a\vert =11`$.

Now choose

```math
\delta=\min\left(1,\frac{\epsilon}{11}\right).
```

Whenever $`\vert x-a\vert \lt \delta`$, both required bounds hold, and therefore

```math
\vert f(x)-f(a)\vert \lt  11\delta\le\epsilon.
```

For the supplied tolerance this choice is $`\delta=\frac{1}{220}`$. The base-point value is $`f(a)=25`$. Every point within the chosen delta has its squared value within epsilon of that number.

The order of the choices is part of the proof. The base point is fixed first, epsilon is supplied next, and delta is then chosen using only those two quantities. The point $`x`$ is allowed to vary afterwards.

This argument also gives a local Lipschitz bound $`\vert x^2-a^2\vert \le 11\vert x-a\vert`$ within distance one of this base point. It does not produce a global Lipschitz constant on the whole real line. Indeed, comparing $`x`$ and $`x+1`$ gives a difference $`\vert 2x+1\vert`$, which has no finite global upper bound. Continuity at every point and one global Lipschitz estimate are different assertions.

### 9.3 Taylor approximation with a remainder

For $`\log(1+x)`$, the first two Taylor terms are easy to write. The useful question is the size of what was omitted. A derivative bound must hold throughout the interval between zero and the evaluation point, including when that point is negative.

#### Problem AN-33: logarithm with a certified error, case 1

Approximate $`\log(1+x)`$ at $`x=\frac{1}{10}`$ by its degree-two Taylor polynomial at zero. Derive a rigorous error bound and compare it with the actual numerical error.

Solution. The derivatives are $`f'(x)=1/(1+x)`$, $`f''(x)=-1/(1+x)^2`$, and $`f'''(x)=2/(1+x)^3`$. At zero, the constant term is zero, the first derivative is one, and the second derivative is minus one. Thus

```math
P_2(x)=x-\frac{x^2}2,\qquad P_2(\frac{1}{10})=\frac{19}{200}.
```

Every point $`\xi`$ between zero and the supplied $`x`$ lies in $`[-q,q]`$ for $`q=\frac{1}{10}\lt 1`$. On this whole interval, $`1+\xi\ge1-q\gt 0`$, so

```math
\vert f'''(\xi)\vert \le \frac2{(1-q)^3}.
```

Taylor's theorem with Lagrange remainder then gives

```math
\vert \log(1+x)-P_2(x)\vert
\le\frac{\vert x\vert ^3}{3(1-q)^3}
=\frac{1}{2187}.
```

Numerically the approximation is $`0.09500000`$, the logarithm is $`0.09531018`$, and the actual absolute error is about $`0.00031018`$. It is below the bound $`0.00045725`$.

The bound is deliberately uniform over the interval. It can be conservative because the largest derivative on the interval need not occur near the unknown remainder point. Nevertheless it certifies the approximation without knowing that point.

For $`x\gt 0`$, the third derivative is positive and the remainder has the sign of $`x^3`$, so the quadratic approximation is below the logarithm. For $`-1\lt x\lt 0`$ the sign reverses. Near $`x=-1`$ the derivative bound grows sharply; a Taylor approximation about zero becomes less useful there. The domain restriction $`x\gt -1`$ must be checked before any expansion is used.

#### Problem AN-34: logarithm with a certified error, case 2

Approximate $`\log(1+x)`$ at $`x=-\frac{1}{5}`$ by its degree-two Taylor polynomial at zero. Derive a rigorous error bound and compare it with the actual numerical error.

Solution. The derivatives are $`f'(x)=1/(1+x)`$, $`f''(x)=-1/(1+x)^2`$, and $`f'''(x)=2/(1+x)^3`$. At zero, the constant term is zero, the first derivative is one, and the second derivative is minus one. Thus

```math
P_2(x)=x-\frac{x^2}2,\qquad P_2(-\frac{1}{5})=-\frac{11}{50}.
```

Every point $`\xi`$ between zero and the supplied $`x`$ lies in $`[-q,q]`$ for $`q=\frac{1}{5}\lt 1`$. On this whole interval, $`1+\xi\ge1-q\gt 0`$, so

```math
\vert f'''(\xi)\vert \le \frac2{(1-q)^3}.
```

Taylor's theorem with Lagrange remainder then gives

```math
\vert \log(1+x)-P_2(x)\vert
\le\frac{\vert x\vert ^3}{3(1-q)^3}
=\frac{1}{192}.
```

Numerically the approximation is $`-0.22000000`$, the logarithm is $`-0.22314355`$, and the actual absolute error is about $`0.00314355`$. It is below the bound $`0.00520833`$.

The bound is deliberately uniform over the interval. It can be conservative because the largest derivative on the interval need not occur near the unknown remainder point. Nevertheless it certifies the approximation without knowing that point.

For $`x\gt 0`$, the third derivative is positive and the remainder has the sign of $`x^3`$, so the quadratic approximation is below the logarithm. For $`-1\lt x\lt 0`$ the sign reverses. Near $`x=-1`$ the derivative bound grows sharply; a Taylor approximation about zero becomes less useful there. The domain restriction $`x\gt -1`$ must be checked before any expansion is used.

#### Problem AN-35: logarithm with a certified error, case 3

Approximate $`\log(1+x)`$ at $`x=\frac{1}{4}`$ by its degree-two Taylor polynomial at zero. Derive a rigorous error bound and compare it with the actual numerical error.

Solution. The derivatives are $`f'(x)=1/(1+x)`$, $`f''(x)=-1/(1+x)^2`$, and $`f'''(x)=2/(1+x)^3`$. At zero, the constant term is zero, the first derivative is one, and the second derivative is minus one. Thus

```math
P_2(x)=x-\frac{x^2}2,\qquad P_2(\frac{1}{4})=\frac{7}{32}.
```

Every point $`\xi`$ between zero and the supplied $`x`$ lies in $`[-q,q]`$ for $`q=\frac{1}{4}\lt 1`$. On this whole interval, $`1+\xi\ge1-q\gt 0`$, so

```math
\vert f'''(\xi)\vert \le \frac2{(1-q)^3}.
```

Taylor's theorem with Lagrange remainder then gives

```math
\vert \log(1+x)-P_2(x)\vert
\le\frac{\vert x\vert ^3}{3(1-q)^3}
=\frac{1}{81}.
```

Numerically the approximation is $`0.21875000`$, the logarithm is $`0.22314355`$, and the actual absolute error is about $`0.00439355`$. It is below the bound $`0.01234568`$.

The bound is deliberately uniform over the interval. It can be conservative because the largest derivative on the interval need not occur near the unknown remainder point. Nevertheless it certifies the approximation without knowing that point.

For $`x\gt 0`$, the third derivative is positive and the remainder has the sign of $`x^3`$, so the quadratic approximation is below the logarithm. For $`-1\lt x\lt 0`$ the sign reverses. Near $`x=-1`$ the derivative bound grows sharply; a Taylor approximation about zero becomes less useful there. The domain restriction $`x\gt -1`$ must be checked before any expansion is used.

#### Problem AN-36: logarithm with a certified error, case 4

Approximate $`\log(1+x)`$ at $`x=-\frac{1}{3}`$ by its degree-two Taylor polynomial at zero. Derive a rigorous error bound and compare it with the actual numerical error.

Solution. The derivatives are $`f'(x)=1/(1+x)`$, $`f''(x)=-1/(1+x)^2`$, and $`f'''(x)=2/(1+x)^3`$. At zero, the constant term is zero, the first derivative is one, and the second derivative is minus one. Thus

```math
P_2(x)=x-\frac{x^2}2,\qquad P_2(-\frac{1}{3})=-\frac{7}{18}.
```

Every point $`\xi`$ between zero and the supplied $`x`$ lies in $`[-q,q]`$ for $`q=\frac{1}{3}\lt 1`$. On this whole interval, $`1+\xi\ge1-q\gt 0`$, so

```math
\vert f'''(\xi)\vert \le \frac2{(1-q)^3}.
```

Taylor's theorem with Lagrange remainder then gives

```math
\vert \log(1+x)-P_2(x)\vert
\le\frac{\vert x\vert ^3}{3(1-q)^3}
=\frac{1}{24}.
```

Numerically the approximation is $`-0.38888889`$, the logarithm is $`-0.40546511`$, and the actual absolute error is about $`0.01657622`$. It is below the bound $`0.04166667`$.

The bound is deliberately uniform over the interval. It can be conservative because the largest derivative on the interval need not occur near the unknown remainder point. Nevertheless it certifies the approximation without knowing that point.

For $`x\gt 0`$, the third derivative is positive and the remainder has the sign of $`x^3`$, so the quadratic approximation is below the logarithm. For $`-1\lt x\lt 0`$ the sign reverses. Near $`x=-1`$ the derivative bound grows sharply; a Taylor approximation about zero becomes less useful there. The domain restriction $`x\gt -1`$ must be checked before any expansion is used.

### 9.4 Exact minimization of a convex quadratic

For a quadratic, the Taylor expansion around a stationary point has no remainder. This makes it possible to prove global optimality with one exact identity. Checking positive definiteness remains necessary: solving the gradient equation alone only finds a stationary point.

#### Problem AN-37: two-variable convex quadratic, case 1

Minimize $`f(x)=\tfrac12x^\top Hx-b^\top x`$ for

```math
H=\begin{pmatrix}2&1\\[0pt] 1&3\end{pmatrix}
```

 and

```math
b=\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

. Prove existence and uniqueness and compute the minimum value.

Solution. The gradient is $`\nabla f(x)=Hx-b`$ because $`H`$ is symmetric. The Hessian is the constant matrix $`H`$. Its first leading principal minor is $`2\gt 0`$ and its determinant is $`5\gt 0`$. Completing a square in the quadratic form gives

```math
h^\top Hh
=2\left(h_1+\frac{1}{2}h_2\right)^2
+\frac{5}{2}h_2^2\gt 0
\quad(h\ne0).
```

Thus the objective is strictly convex. The stationary equation $`Hx=b`$ has the solution

```math
x^*=H^{-1}b=\begin{pmatrix}\frac{1}{5}\\[0pt] \frac{3}{5}\end{pmatrix},\qquad
H^{-1}=\begin{pmatrix}\frac{3}{5}&-\frac{1}{5}\\[0pt] -\frac{1}{5}&\frac{2}{5}\end{pmatrix}.
```

Multiplying $`H`$ by the proposed vector checks both component equations. At this vector, $`b=Hx^*`$, so the objective simplifies to

```math
f(x^*)=-\frac12b^\top x^*=-\frac{7}{10}.
```

Now take an arbitrary vector $`h`$ and expand $`f(x^*+h)`$. The terms linear in $`h`$ equal $`h^\top(Hx^*-b)=0`$. Therefore

```math
f(x^*+h)-f(x^*)=\frac12h^\top Hh.
```

The right side is positive for every nonzero $`h`$. This proves that the stationary solution exists, is a global minimizer, and is the only minimizer. It also gives the exact cost of any coefficient error.

Positive definiteness implies coercivity here: as $`\Vert x\Vert _2`$ grows, the positive quadratic term eventually dominates the linear term. A continuous coercive function attains a minimum on some sufficiently large closed ball. The explicit stationary solution agrees with that existence argument.

#### Problem AN-38: two-variable convex quadratic, case 2

Minimize $`f(x)=\tfrac12x^\top Hx-b^\top x`$ for

```math
H=\begin{pmatrix}4&-1\\[0pt] -1&2\end{pmatrix}
```

 and

```math
b=\begin{pmatrix}3\\[0pt] -2\end{pmatrix}
```

. Prove existence and uniqueness and compute the minimum value.

Solution. The gradient is $`\nabla f(x)=Hx-b`$ because $`H`$ is symmetric. The Hessian is the constant matrix $`H`$. Its first leading principal minor is $`4\gt 0`$ and its determinant is $`7\gt 0`$. Completing a square in the quadratic form gives

```math
h^\top Hh
=4\left(h_1+\frac{-1}{4}h_2\right)^2
+\frac{7}{4}h_2^2\gt 0
\quad(h\ne0).
```

Thus the objective is strictly convex. The stationary equation $`Hx=b`$ has the solution

```math
x^*=H^{-1}b=\begin{pmatrix}\frac{4}{7}\\[0pt] -\frac{5}{7}\end{pmatrix},\qquad
H^{-1}=\begin{pmatrix}\frac{2}{7}&\frac{1}{7}\\[0pt] \frac{1}{7}&\frac{4}{7}\end{pmatrix}.
```

Multiplying $`H`$ by the proposed vector checks both component equations. At this vector, $`b=Hx^*`$, so the objective simplifies to

```math
f(x^*)=-\frac12b^\top x^*=-\frac{11}{7}.
```

Now take an arbitrary vector $`h`$ and expand $`f(x^*+h)`$. The terms linear in $`h`$ equal $`h^\top(Hx^*-b)=0`$. Therefore

```math
f(x^*+h)-f(x^*)=\frac12h^\top Hh.
```

The right side is positive for every nonzero $`h`$. This proves that the stationary solution exists, is a global minimizer, and is the only minimizer. It also gives the exact cost of any coefficient error.

Positive definiteness implies coercivity here: as $`\Vert x\Vert _2`$ grows, the positive quadratic term eventually dominates the linear term. A continuous coercive function attains a minimum on some sufficiently large closed ball. The explicit stationary solution agrees with that existence argument.

#### Problem AN-39: two-variable convex quadratic, case 3

Minimize $`f(x)=\tfrac12x^\top Hx-b^\top x`$ for

```math
H=\begin{pmatrix}3&0\\[0pt] 0&5\end{pmatrix}
```

 and

```math
b=\begin{pmatrix}-1\\[0pt] 4\end{pmatrix}
```

. Prove existence and uniqueness and compute the minimum value.

Solution. The gradient is $`\nabla f(x)=Hx-b`$ because $`H`$ is symmetric. The Hessian is the constant matrix $`H`$. Its first leading principal minor is $`3\gt 0`$ and its determinant is $`15\gt 0`$. Completing a square in the quadratic form gives

```math
h^\top Hh
=3\left(h_1+\frac{0}{3}h_2\right)^2
+\frac{15}{3}h_2^2\gt 0
\quad(h\ne0).
```

Thus the objective is strictly convex. The stationary equation $`Hx=b`$ has the solution

```math
x^*=H^{-1}b=\begin{pmatrix}-\frac{1}{3}\\[0pt] \frac{4}{5}\end{pmatrix},\qquad
H^{-1}=\begin{pmatrix}\frac{1}{3}&0\\[0pt] 0&\frac{1}{5}\end{pmatrix}.
```

Multiplying $`H`$ by the proposed vector checks both component equations. At this vector, $`b=Hx^*`$, so the objective simplifies to

```math
f(x^*)=-\frac12b^\top x^*=-\frac{53}{30}.
```

Now take an arbitrary vector $`h`$ and expand $`f(x^*+h)`$. The terms linear in $`h`$ equal $`h^\top(Hx^*-b)=0`$. Therefore

```math
f(x^*+h)-f(x^*)=\frac12h^\top Hh.
```

The right side is positive for every nonzero $`h`$. This proves that the stationary solution exists, is a global minimizer, and is the only minimizer. It also gives the exact cost of any coefficient error.

Positive definiteness implies coercivity here: as $`\Vert x\Vert _2`$ grows, the positive quadratic term eventually dominates the linear term. A continuous coercive function attains a minimum on some sufficiently large closed ball. The explicit stationary solution agrees with that existence argument.

#### Problem AN-40: two-variable convex quadratic, case 4

Minimize $`f(x)=\tfrac12x^\top Hx-b^\top x`$ for

```math
H=\begin{pmatrix}5&2\\[0pt] 2&2\end{pmatrix}
```

 and

```math
b=\begin{pmatrix}2\\[0pt] 1\end{pmatrix}
```

. Prove existence and uniqueness and compute the minimum value.

Solution. The gradient is $`\nabla f(x)=Hx-b`$ because $`H`$ is symmetric. The Hessian is the constant matrix $`H`$. Its first leading principal minor is $`5\gt 0`$ and its determinant is $`6\gt 0`$. Completing a square in the quadratic form gives

```math
h^\top Hh
=5\left(h_1+\frac{2}{5}h_2\right)^2
+\frac{6}{5}h_2^2\gt 0
\quad(h\ne0).
```

Thus the objective is strictly convex. The stationary equation $`Hx=b`$ has the solution

```math
x^*=H^{-1}b=\begin{pmatrix}\frac{1}{3}\\[0pt] \frac{1}{6}\end{pmatrix},\qquad
H^{-1}=\begin{pmatrix}\frac{1}{3}&-\frac{1}{3}\\[0pt] -\frac{1}{3}&\frac{5}{6}\end{pmatrix}.
```

Multiplying $`H`$ by the proposed vector checks both component equations. At this vector, $`b=Hx^*`$, so the objective simplifies to

```math
f(x^*)=-\frac12b^\top x^*=-\frac{5}{12}.
```

Now take an arbitrary vector $`h`$ and expand $`f(x^*+h)`$. The terms linear in $`h`$ equal $`h^\top(Hx^*-b)=0`$. Therefore

```math
f(x^*+h)-f(x^*)=\frac12h^\top Hh.
```

The right side is positive for every nonzero $`h`$. This proves that the stationary solution exists, is a global minimizer, and is the only minimizer. It also gives the exact cost of any coefficient error.

Positive definiteness implies coercivity here: as $`\Vert x\Vert _2`$ grows, the positive quadratic term eventually dominates the linear term. A continuous coercive function attains a minimum on some sufficiently large closed ball. The explicit stationary solution agrees with that existence argument.

### 9.5 Projection onto a half-space

The closest point in a closed half-space is either the original point, if feasible, or a point on the boundary. The multiplier determines how far to move in the normal direction. Complementary slackness distinguishes these two cases without guessing that every inequality must bind.

#### Problem AN-41: half-space projection, case 1

Minimize $`\tfrac12\Vert x-z\Vert _2^2`$ subject to $`a^\top x\le b`$, with

```math
a=\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

,

```math
z=\begin{pmatrix}2\\[0pt] 2\end{pmatrix}
```

, and $`b=3`$.

Solution. The Lagrangian is $`L(x,\lambda)=\tfrac12\Vert x-z\Vert _2^2+\lambda(a^\top x-b)`$ with $`\lambda\ge0`$. Differentiation in $`x`$ gives $`x-z+\lambda a=0`$, hence $`x=z-\lambda a`$.

At the unconstrained minimum $`z`$, the constraint violation is $`a^\top z-b=3`$. If this is nonpositive, $`z`$ is feasible and no movement is needed. Otherwise the closest feasible point must lie on the boundary, and substituting stationarity into $`a^\top x=b`$ gives

```math
\lambda=\frac{a^\top z-b}{a^\top a}.
```

Both cases are summarized by

```math
\lambda^*=\max\left(0,\frac{a^\top z-b}{a^\top a}\right)
=\frac{3}{5},\qquad
x^*=z-\lambda^*a=\begin{pmatrix}\frac{7}{5}\\[0pt] \frac{4}{5}\end{pmatrix}.
```

For these data, $`a^\top a=5`$, $`a^\top x^*-b=0`$, and the minimum objective is $`\frac{9}{10}`$. The multiplier is nonnegative, stationarity holds, and $`\lambda^*(a^\top x^*-b)=0`$.

A direct geometric proof gives another check. For any feasible $`y`$,

```math
\tfrac12\Vert y-z\Vert _2^2-\tfrac12\Vert x^*-z\Vert _2^2
=\tfrac12\Vert y-x^*\Vert _2^2
-\lambda^*a^\top(y-x^*).
```

If the multiplier is positive, the boundary equality and feasibility make the last term nonnegative. If the multiplier is zero, that term vanishes. The difference is therefore nonnegative in either case, and is strictly positive for $`y\ne x^*`$. This proves global optimality and uniqueness.

#### Problem AN-42: half-space projection, case 2

Minimize $`\tfrac12\Vert x-z\Vert _2^2`$ subject to $`a^\top x\le b`$, with

```math
a=\begin{pmatrix}2\\[0pt] -1\end{pmatrix}
```

,

```math
z=\begin{pmatrix}0\\[0pt] 0\end{pmatrix}
```

, and $`b=1`$.

Solution. The Lagrangian is $`L(x,\lambda)=\tfrac12\Vert x-z\Vert _2^2+\lambda(a^\top x-b)`$ with $`\lambda\ge0`$. Differentiation in $`x`$ gives $`x-z+\lambda a=0`$, hence $`x=z-\lambda a`$.

At the unconstrained minimum $`z`$, the constraint violation is $`a^\top z-b=-1`$. If this is nonpositive, $`z`$ is feasible and no movement is needed. Otherwise the closest feasible point must lie on the boundary, and substituting stationarity into $`a^\top x=b`$ gives

```math
\lambda=\frac{a^\top z-b}{a^\top a}.
```

Both cases are summarized by

```math
\lambda^*=\max\left(0,\frac{a^\top z-b}{a^\top a}\right)
=0,\qquad
x^*=z-\lambda^*a=\begin{pmatrix}0\\[0pt] 0\end{pmatrix}.
```

For these data, $`a^\top a=5`$, $`a^\top x^*-b=-1`$, and the minimum objective is $`0`$. The multiplier is nonnegative, stationarity holds, and $`\lambda^*(a^\top x^*-b)=0`$.

A direct geometric proof gives another check. For any feasible $`y`$,

```math
\tfrac12\Vert y-z\Vert _2^2-\tfrac12\Vert x^*-z\Vert _2^2
=\tfrac12\Vert y-x^*\Vert _2^2
-\lambda^*a^\top(y-x^*).
```

If the multiplier is positive, the boundary equality and feasibility make the last term nonnegative. If the multiplier is zero, that term vanishes. The difference is therefore nonnegative in either case, and is strictly positive for $`y\ne x^*`$. This proves global optimality and uniqueness.

#### Problem AN-43: half-space projection, case 3

Minimize $`\tfrac12\Vert x-z\Vert _2^2`$ subject to $`a^\top x\le b`$, with

```math
a=\begin{pmatrix}1\\[0pt] 1\end{pmatrix}
```

,

```math
z=\begin{pmatrix}-2\\[0pt] 1\end{pmatrix}
```

, and $`b=-2`$.

Solution. The Lagrangian is $`L(x,\lambda)=\tfrac12\Vert x-z\Vert _2^2+\lambda(a^\top x-b)`$ with $`\lambda\ge0`$. Differentiation in $`x`$ gives $`x-z+\lambda a=0`$, hence $`x=z-\lambda a`$.

At the unconstrained minimum $`z`$, the constraint violation is $`a^\top z-b=1`$. If this is nonpositive, $`z`$ is feasible and no movement is needed. Otherwise the closest feasible point must lie on the boundary, and substituting stationarity into $`a^\top x=b`$ gives

```math
\lambda=\frac{a^\top z-b}{a^\top a}.
```

Both cases are summarized by

```math
\lambda^*=\max\left(0,\frac{a^\top z-b}{a^\top a}\right)
=\frac{1}{2},\qquad
x^*=z-\lambda^*a=\begin{pmatrix}-\frac{5}{2}\\[0pt] \frac{1}{2}\end{pmatrix}.
```

For these data, $`a^\top a=2`$, $`a^\top x^*-b=0`$, and the minimum objective is $`\frac{1}{4}`$. The multiplier is nonnegative, stationarity holds, and $`\lambda^*(a^\top x^*-b)=0`$.

A direct geometric proof gives another check. For any feasible $`y`$,

```math
\tfrac12\Vert y-z\Vert _2^2-\tfrac12\Vert x^*-z\Vert _2^2
=\tfrac12\Vert y-x^*\Vert _2^2
-\lambda^*a^\top(y-x^*).
```

If the multiplier is positive, the boundary equality and feasibility make the last term nonnegative. If the multiplier is zero, that term vanishes. The difference is therefore nonnegative in either case, and is strictly positive for $`y\ne x^*`$. This proves global optimality and uniqueness.

#### Problem AN-44: half-space projection, case 4

Minimize $`\tfrac12\Vert x-z\Vert _2^2`$ subject to $`a^\top x\le b`$, with

```math
a=\begin{pmatrix}3\\[0pt] 1\end{pmatrix}
```

,

```math
z=\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

, and $`b=0`$.

Solution. The Lagrangian is $`L(x,\lambda)=\tfrac12\Vert x-z\Vert _2^2+\lambda(a^\top x-b)`$ with $`\lambda\ge0`$. Differentiation in $`x`$ gives $`x-z+\lambda a=0`$, hence $`x=z-\lambda a`$.

At the unconstrained minimum $`z`$, the constraint violation is $`a^\top z-b=5`$. If this is nonpositive, $`z`$ is feasible and no movement is needed. Otherwise the closest feasible point must lie on the boundary, and substituting stationarity into $`a^\top x=b`$ gives

```math
\lambda=\frac{a^\top z-b}{a^\top a}.
```

Both cases are summarized by

```math
\lambda^*=\max\left(0,\frac{a^\top z-b}{a^\top a}\right)
=\frac{1}{2},\qquad
x^*=z-\lambda^*a=\begin{pmatrix}-\frac{1}{2}\\[0pt] \frac{3}{2}\end{pmatrix}.
```

For these data, $`a^\top a=10`$, $`a^\top x^*-b=0`$, and the minimum objective is $`\frac{5}{4}`$. The multiplier is nonnegative, stationarity holds, and $`\lambda^*(a^\top x^*-b)=0`$.

A direct geometric proof gives another check. For any feasible $`y`$,

```math
\tfrac12\Vert y-z\Vert _2^2-\tfrac12\Vert x^*-z\Vert _2^2
=\tfrac12\Vert y-x^*\Vert _2^2
-\lambda^*a^\top(y-x^*).
```

If the multiplier is positive, the boundary equality and feasibility make the last term nonnegative. If the multiplier is zero, that term vanishes. The difference is therefore nonnegative in either case, and is strictly positive for $`y\ne x^*`$. This proves global optimality and uniqueness.

### 9.6 A dual function with active and inactive cases

The dual function is found by minimizing the Lagrangian over the primal variable without retaining the primal constraint. The multiplier domain remains $`\lambda\ge0`$. Solving this small example in both primal and dual variables makes the role of that domain explicit.

#### Problem AN-45: quadratic dual and sensitivity, case 1

Minimize $`f(x)=\tfrac12(2)x^2+(-4)x`$ subject to $`x\ge 3`$. Derive and maximize the dual function and interpret its multiplier.

Solution. The unconstrained minimizer solves $`2x+(-4)=0`$, giving $`x_0=2`$. The positive second derivative $`2`$ means that the constrained minimizer is the closest feasible point to $`x_0`$:

```math
x^*=\max(b,-a/q)=3,\qquad f(x^*)=-3.
```

Write the constraint as $`b-x\le0`$. The Lagrangian and its unconstrained minimizer are

```math
L(x,\lambda)=\tfrac q2x^2+ax+\lambda(b-x),\qquad
x(\lambda)=\frac{\lambda-a}q.
```

Completing the square gives

```math
g(\lambda)=\inf_xL(x,\lambda)
=\lambda b-\frac{(a-\lambda)^2}{2q},
\qquad \lambda\ge0.
```

The derivative is $`g'(\lambda)=b+(a-\lambda)/q`$. The unrestricted maximum would occur at $`qb+a`$. Restricting to the nonnegative half-line gives $`\lambda^*=\max(0,qb+a)=2`$. Substituting it into $`x(\lambda)`$ recovers the primal solution. Direct substitution also gives $`g(\lambda^*)=f(x^*)=-3`$.

Thus the duality gap is zero. Complementary slackness $`\lambda^*(b-x^*)=0`$ says that a strictly slack lower bound has zero multiplier. A binding lower bound can also have zero multiplier exactly where it first meets the unconstrained optimum.

If the lower bound is treated as a parameter, the value function is constant below $`-a/q`$ and equals $`qb^2/2+ab`$ above it. Its derivative is $`\max(0,qb+a)`$. Increasing this bound tightens the feasible set, so a positive multiplier means a positive marginal increase in minimum cost. The sign follows from the chosen parameterization.

#### Problem AN-46: quadratic dual and sensitivity, case 2

Minimize $`f(x)=\tfrac12(3)x^2+(6)x`$ subject to $`x\ge -1`$. Derive and maximize the dual function and interpret its multiplier.

Solution. The unconstrained minimizer solves $`3x+(6)=0`$, giving $`x_0=-2`$. The positive second derivative $`3`$ means that the constrained minimizer is the closest feasible point to $`x_0`$:

```math
x^*=\max(b,-a/q)=-1,\qquad f(x^*)=-\frac{9}{2}.
```

Write the constraint as $`b-x\le0`$. The Lagrangian and its unconstrained minimizer are

```math
L(x,\lambda)=\tfrac q2x^2+ax+\lambda(b-x),\qquad
x(\lambda)=\frac{\lambda-a}q.
```

Completing the square gives

```math
g(\lambda)=\inf_xL(x,\lambda)
=\lambda b-\frac{(a-\lambda)^2}{2q},
\qquad \lambda\ge0.
```

The derivative is $`g'(\lambda)=b+(a-\lambda)/q`$. The unrestricted maximum would occur at $`qb+a`$. Restricting to the nonnegative half-line gives $`\lambda^*=\max(0,qb+a)=3`$. Substituting it into $`x(\lambda)`$ recovers the primal solution. Direct substitution also gives $`g(\lambda^*)=f(x^*)=-\frac{9}{2}`$.

Thus the duality gap is zero. Complementary slackness $`\lambda^*(b-x^*)=0`$ says that a strictly slack lower bound has zero multiplier. A binding lower bound can also have zero multiplier exactly where it first meets the unconstrained optimum.

If the lower bound is treated as a parameter, the value function is constant below $`-a/q`$ and equals $`qb^2/2+ab`$ above it. Its derivative is $`\max(0,qb+a)`$. Increasing this bound tightens the feasible set, so a positive multiplier means a positive marginal increase in minimum cost. The sign follows from the chosen parameterization.

#### Problem AN-47: quadratic dual and sensitivity, case 3

Minimize $`f(x)=\tfrac12(4)x^2+(-8)x`$ subject to $`x\ge 1`$. Derive and maximize the dual function and interpret its multiplier.

Solution. The unconstrained minimizer solves $`4x+(-8)=0`$, giving $`x_0=2`$. The positive second derivative $`4`$ means that the constrained minimizer is the closest feasible point to $`x_0`$:

```math
x^*=\max(b,-a/q)=2,\qquad f(x^*)=-8.
```

Write the constraint as $`b-x\le0`$. The Lagrangian and its unconstrained minimizer are

```math
L(x,\lambda)=\tfrac q2x^2+ax+\lambda(b-x),\qquad
x(\lambda)=\frac{\lambda-a}q.
```

Completing the square gives

```math
g(\lambda)=\inf_xL(x,\lambda)
=\lambda b-\frac{(a-\lambda)^2}{2q},
\qquad \lambda\ge0.
```

The derivative is $`g'(\lambda)=b+(a-\lambda)/q`$. The unrestricted maximum would occur at $`qb+a`$. Restricting to the nonnegative half-line gives $`\lambda^*=\max(0,qb+a)=0`$. Substituting it into $`x(\lambda)`$ recovers the primal solution. Direct substitution also gives $`g(\lambda^*)=f(x^*)=-8`$.

Thus the duality gap is zero. Complementary slackness $`\lambda^*(b-x^*)=0`$ says that a strictly slack lower bound has zero multiplier. A binding lower bound can also have zero multiplier exactly where it first meets the unconstrained optimum.

If the lower bound is treated as a parameter, the value function is constant below $`-a/q`$ and equals $`qb^2/2+ab`$ above it. Its derivative is $`\max(0,qb+a)`$. Increasing this bound tightens the feasible set, so a positive multiplier means a positive marginal increase in minimum cost. The sign follows from the chosen parameterization.

#### Problem AN-48: quadratic dual and sensitivity, case 4

Minimize $`f(x)=\tfrac12(1)x^2+(2)x`$ subject to $`x\ge -3`$. Derive and maximize the dual function and interpret its multiplier.

Solution. The unconstrained minimizer solves $`1x+(2)=0`$, giving $`x_0=-2`$. The positive second derivative $`1`$ means that the constrained minimizer is the closest feasible point to $`x_0`$:

```math
x^*=\max(b,-a/q)=-2,\qquad f(x^*)=-2.
```

Write the constraint as $`b-x\le0`$. The Lagrangian and its unconstrained minimizer are

```math
L(x,\lambda)=\tfrac q2x^2+ax+\lambda(b-x),\qquad
x(\lambda)=\frac{\lambda-a}q.
```

Completing the square gives

```math
g(\lambda)=\inf_xL(x,\lambda)
=\lambda b-\frac{(a-\lambda)^2}{2q},
\qquad \lambda\ge0.
```

The derivative is $`g'(\lambda)=b+(a-\lambda)/q`$. The unrestricted maximum would occur at $`qb+a`$. Restricting to the nonnegative half-line gives $`\lambda^*=\max(0,qb+a)=0`$. Substituting it into $`x(\lambda)`$ recovers the primal solution. Direct substitution also gives $`g(\lambda^*)=f(x^*)=-2`$.

Thus the duality gap is zero. Complementary slackness $`\lambda^*(b-x^*)=0`$ says that a strictly slack lower bound has zero multiplier. A binding lower bound can also have zero multiplier exactly where it first meets the unconstrained optimum.

If the lower bound is treated as a parameter, the value function is constant below $`-a/q`$ and equals $`qb^2/2+ab`$ above it. Its derivative is $`\max(0,qb+a)`$. Increasing this bound tightens the feasible set, so a positive multiplier means a positive marginal increase in minimum cost. The sign follows from the chosen parameterization.

### 9.7 Box uncertainty and sign-dependent worst cases

When coefficients vary independently in intervals, the worst coefficient for a positive decision coordinate is its upper endpoint. For a negative coordinate it is its lower endpoint. The support function formula records both choices through absolute values and avoids silently assuming nonnegative decisions.

#### Problem AN-49: robust scaling in a fixed direction, case 1

Choose $`x=td`$ with $`t\ge0`$,

```math
d=\begin{pmatrix}1\\[0pt] 1\end{pmatrix}
```

, as large as possible subject to $`a^\top x\le 12`$ for every $`a_i\in[a_{0i}-r_i,a_{0i}+r_i]`$. Here

```math
a_0=\begin{pmatrix}2\\[0pt] 1\end{pmatrix}
```

 and

```math
r=\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

.

Solution. Write $`a=a_0+u`$ with $`\vert u_i\vert \le r_i`$. For a fixed $`x`$,

```math
\sup_a a^\top x=a_0^\top x+\sum_i r_i\vert x_i\vert .
```

To prove the formula, bound each product by $`u_ix_i\le r_i\vert x_i\vert`$. Equality is attained by choosing $`u_i=r_i`$ when $`x_i\gt 0`$ and $`u_i=-r_i`$ when $`x_i\lt 0`$; a zero coordinate allows either endpoint.

Since $`x=td`$ with $`t\ge0`$, the worst-case left side is

```math
t\left(a_0^\top d+\sum_i r_i\vert d_i\vert \right)
=t(6).
```

The coefficient is strictly positive for the supplied data. The robust feasible interval is therefore $`0\le t\le 2`$, and the largest scaling is $`t^*=2`$. The resulting decision is

```math
x^*=\begin{pmatrix}2\\[0pt] 2\end{pmatrix}
```

.

One worst-case coefficient vector is

```math
a^*=\begin{pmatrix}3\\[0pt] 3\end{pmatrix}
```

. Checking its inner product gives $`(a^*)^\top x^*=12`$, so the robust constraint binds. Every other allowed coefficient vector has an inner product no larger, by the support-function inequality above.

The absolute values are essential when a coordinate of $`d`$ is negative. Replacing the uncertain vector by all upper endpoints would generally fail in that case. If the computed support coefficient had been zero or negative while the right side remained positive, this particular constraint would impose no finite upper bound on $`t`$; dividing by it as though positive would be invalid.

#### Problem AN-50: robust scaling in a fixed direction, case 2

Choose $`x=td`$ with $`t\ge0`$,

```math
d=\begin{pmatrix}1\\[0pt] -2\end{pmatrix}
```

, as large as possible subject to $`a^\top x\le 14`$ for every $`a_i\in[a_{0i}-r_i,a_{0i}+r_i]`$. Here

```math
a_0=\begin{pmatrix}3\\[0pt] -1\end{pmatrix}
```

 and

```math
r=\begin{pmatrix}1\\[0pt] 1\end{pmatrix}
```

.

Solution. Write $`a=a_0+u`$ with $`\vert u_i\vert \le r_i`$. For a fixed $`x`$,

```math
\sup_a a^\top x=a_0^\top x+\sum_i r_i\vert x_i\vert .
```

To prove the formula, bound each product by $`u_ix_i\le r_i\vert x_i\vert`$. Equality is attained by choosing $`u_i=r_i`$ when $`x_i\gt 0`$ and $`u_i=-r_i`$ when $`x_i\lt 0`$; a zero coordinate allows either endpoint.

Since $`x=td`$ with $`t\ge0`$, the worst-case left side is

```math
t\left(a_0^\top d+\sum_i r_i\vert d_i\vert \right)
=t(8).
```

The coefficient is strictly positive for the supplied data. The robust feasible interval is therefore $`0\le t\le \frac{7}{4}`$, and the largest scaling is $`t^*=\frac{7}{4}`$. The resulting decision is

```math
x^*=\begin{pmatrix}\frac{7}{4}\\[0pt] -\frac{7}{2}\end{pmatrix}
```

.

One worst-case coefficient vector is

```math
a^*=\begin{pmatrix}4\\[0pt] -2\end{pmatrix}
```

. Checking its inner product gives $`(a^*)^\top x^*=14`$, so the robust constraint binds. Every other allowed coefficient vector has an inner product no larger, by the support-function inequality above.

The absolute values are essential when a coordinate of $`d`$ is negative. Replacing the uncertain vector by all upper endpoints would generally fail in that case. If the computed support coefficient had been zero or negative while the right side remained positive, this particular constraint would impose no finite upper bound on $`t`$; dividing by it as though positive would be invalid.

#### Problem AN-51: robust scaling in a fixed direction, case 3

Choose $`x=td`$ with $`t\ge0`$,

```math
d=\begin{pmatrix}-1\\[0pt] 2\end{pmatrix}
```

, as large as possible subject to $`a^\top x\le 15`$ for every $`a_i\in[a_{0i}-r_i,a_{0i}+r_i]`$. Here

```math
a_0=\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

 and

```math
r=\begin{pmatrix}2\\[0pt] 1\end{pmatrix}
```

.

Solution. Write $`a=a_0+u`$ with $`\vert u_i\vert \le r_i`$. For a fixed $`x`$,

```math
\sup_a a^\top x=a_0^\top x+\sum_i r_i\vert x_i\vert .
```

To prove the formula, bound each product by $`u_ix_i\le r_i\vert x_i\vert`$. Equality is attained by choosing $`u_i=r_i`$ when $`x_i\gt 0`$ and $`u_i=-r_i`$ when $`x_i\lt 0`$; a zero coordinate allows either endpoint.

Since $`x=td`$ with $`t\ge0`$, the worst-case left side is

```math
t\left(a_0^\top d+\sum_i r_i\vert d_i\vert \right)
=t(7).
```

The coefficient is strictly positive for the supplied data. The robust feasible interval is therefore $`0\le t\le \frac{15}{7}`$, and the largest scaling is $`t^*=\frac{15}{7}`$. The resulting decision is

```math
x^*=\begin{pmatrix}-\frac{15}{7}\\[0pt] \frac{30}{7}\end{pmatrix}
```

.

One worst-case coefficient vector is

```math
a^*=\begin{pmatrix}-1\\[0pt] 3\end{pmatrix}
```

. Checking its inner product gives $`(a^*)^\top x^*=15`$, so the robust constraint binds. Every other allowed coefficient vector has an inner product no larger, by the support-function inequality above.

The absolute values are essential when a coordinate of $`d`$ is negative. Replacing the uncertain vector by all upper endpoints would generally fail in that case. If the computed support coefficient had been zero or negative while the right side remained positive, this particular constraint would impose no finite upper bound on $`t`$; dividing by it as though positive would be invalid.

#### Problem AN-52: robust scaling in a fixed direction, case 4

Choose $`x=td`$ with $`t\ge0`$,

```math
d=\begin{pmatrix}-2\\[0pt] 3\end{pmatrix}
```

, as large as possible subject to $`a^\top x\le 20`$ for every $`a_i\in[a_{0i}-r_i,a_{0i}+r_i]`$. Here

```math
a_0=\begin{pmatrix}2\\[0pt] 3\end{pmatrix}
```

 and

```math
r=\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

.

Solution. Write $`a=a_0+u`$ with $`\vert u_i\vert \le r_i`$. For a fixed $`x`$,

```math
\sup_a a^\top x=a_0^\top x+\sum_i r_i\vert x_i\vert .
```

To prove the formula, bound each product by $`u_ix_i\le r_i\vert x_i\vert`$. Equality is attained by choosing $`u_i=r_i`$ when $`x_i\gt 0`$ and $`u_i=-r_i`$ when $`x_i\lt 0`$; a zero coordinate allows either endpoint.

Since $`x=td`$ with $`t\ge0`$, the worst-case left side is

```math
t\left(a_0^\top d+\sum_i r_i\vert d_i\vert \right)
=t(13).
```

The coefficient is strictly positive for the supplied data. The robust feasible interval is therefore $`0\le t\le \frac{20}{13}`$, and the largest scaling is $`t^*=\frac{20}{13}`$. The resulting decision is

```math
x^*=\begin{pmatrix}-\frac{40}{13}\\[0pt] \frac{60}{13}\end{pmatrix}
```

.

One worst-case coefficient vector is

```math
a^*=\begin{pmatrix}1\\[0pt] 5\end{pmatrix}
```

. Checking its inner product gives $`(a^*)^\top x^*=20`$, so the robust constraint binds. Every other allowed coefficient vector has an inner product no larger, by the support-function inequality above.

The absolute values are essential when a coordinate of $`d`$ is negative. Replacing the uncertain vector by all upper endpoints would generally fail in that case. If the computed support coefficient had been zero or negative while the right side remained positive, this particular constraint would impose no finite upper bound on $`t`$; dividing by it as though positive would be invalid.

### 9.8 Maximizing a minimum of two affine functions

The lower of an increasing line and a decreasing line rises until the crossing and falls afterwards. An interval restriction can put the crossing outside the feasible set. The optimizer is then an endpoint. The auxiliary variable formulation gives the same result as this direct argument.

#### Problem AN-53: two-reward max-min problem, case 1

Maximize $`\min\lbrace 1x+(0),\,2-(2)x\rbrace`$ over $`0\le x\le 1`$.

Solution. Set the two affine rewards equal to locate their crossing:

```math
1x+(0)=2-(2)x,
\qquad x_c=\frac{v-u}{p+q}=\frac{2}{3}.
```

The slopes have opposite signs because $`p=1\gt 0`$ and $`q=2\gt 0`$. To the left of the crossing, the increasing reward is the smaller one. To the right, the decreasing reward is the smaller one. The objective is therefore increasing before $`x_c`$ and decreasing after it.

Restricting to the supplied interval gives the clipped solution

```math
x^*=\min(U,\max(L,x_c))=\frac{2}{3}.
```

The two rewards there are $`\frac{2}{3}`$ and $`\frac{2}{3}`$, so the maximum minimum reward is $`\frac{2}{3}`$. If the crossing is outside the interval, the two rewards need not be equal at the optimum. Equalizing them without first checking feasibility would produce an invalid decision.

The equivalent linear optimization problem introduces a reward floor $`z`$:

```math
\max_{x,z}z
\quad\text{subject to}\quad
z\le px+u,\quad z\le v-qx,\quad L\le x\le U.
```

For each fixed $`x`$, the largest feasible $`z`$ is exactly the smaller reward, proving equivalence. This formulation remains useful with more than two affine rewards even when a single crossing argument no longer describes the entire solution. Here the direct monotonicity proof already establishes global optimality and uniqueness.

#### Problem AN-54: two-reward max-min problem, case 2

Maximize $`\min\lbrace 2x+(1),\,7-(1)x\rbrace`$ over $`0\le x\le 1`$.

Solution. Set the two affine rewards equal to locate their crossing:

```math
2x+(1)=7-(1)x,
\qquad x_c=\frac{v-u}{p+q}=2.
```

The slopes have opposite signs because $`p=2\gt 0`$ and $`q=1\gt 0`$. To the left of the crossing, the increasing reward is the smaller one. To the right, the decreasing reward is the smaller one. The objective is therefore increasing before $`x_c`$ and decreasing after it.

Restricting to the supplied interval gives the clipped solution

```math
x^*=\min(U,\max(L,x_c))=1.
```

The two rewards there are $`3`$ and $`6`$, so the maximum minimum reward is $`3`$. If the crossing is outside the interval, the two rewards need not be equal at the optimum. Equalizing them without first checking feasibility would produce an invalid decision.

The equivalent linear optimization problem introduces a reward floor $`z`$:

```math
\max_{x,z}z
\quad\text{subject to}\quad
z\le px+u,\quad z\le v-qx,\quad L\le x\le U.
```

For each fixed $`x`$, the largest feasible $`z`$ is exactly the smaller reward, proving equivalence. This formulation remains useful with more than two affine rewards even when a single crossing argument no longer describes the entire solution. Here the direct monotonicity proof already establishes global optimality and uniqueness.

#### Problem AN-55: two-reward max-min problem, case 3

Maximize $`\min\lbrace 3x+(-2),\,3-(2)x\rbrace`$ over $`-1\le x\le 2`$.

Solution. Set the two affine rewards equal to locate their crossing:

```math
3x+(-2)=3-(2)x,
\qquad x_c=\frac{v-u}{p+q}=1.
```

The slopes have opposite signs because $`p=3\gt 0`$ and $`q=2\gt 0`$. To the left of the crossing, the increasing reward is the smaller one. To the right, the decreasing reward is the smaller one. The objective is therefore increasing before $`x_c`$ and decreasing after it.

Restricting to the supplied interval gives the clipped solution

```math
x^*=\min(U,\max(L,x_c))=1.
```

The two rewards there are $`1`$ and $`1`$, so the maximum minimum reward is $`1`$. If the crossing is outside the interval, the two rewards need not be equal at the optimum. Equalizing them without first checking feasibility would produce an invalid decision.

The equivalent linear optimization problem introduces a reward floor $`z`$:

```math
\max_{x,z}z
\quad\text{subject to}\quad
z\le px+u,\quad z\le v-qx,\quad L\le x\le U.
```

For each fixed $`x`$, the largest feasible $`z`$ is exactly the smaller reward, proving equivalence. This formulation remains useful with more than two affine rewards even when a single crossing argument no longer describes the entire solution. Here the direct monotonicity proof already establishes global optimality and uniqueness.

#### Problem AN-56: two-reward max-min problem, case 4

Maximize $`\min\lbrace 1x+(4),\,0-(3)x\rbrace`$ over $`0\le x\le 2`$.

Solution. Set the two affine rewards equal to locate their crossing:

```math
1x+(4)=0-(3)x,
\qquad x_c=\frac{v-u}{p+q}=-1.
```

The slopes have opposite signs because $`p=1\gt 0`$ and $`q=3\gt 0`$. To the left of the crossing, the increasing reward is the smaller one. To the right, the decreasing reward is the smaller one. The objective is therefore increasing before $`x_c`$ and decreasing after it.

Restricting to the supplied interval gives the clipped solution

```math
x^*=\min(U,\max(L,x_c))=0.
```

The two rewards there are $`4`$ and $`0`$, so the maximum minimum reward is $`0`$. If the crossing is outside the interval, the two rewards need not be equal at the optimum. Equalizing them without first checking feasibility would produce an invalid decision.

The equivalent linear optimization problem introduces a reward floor $`z`$:

```math
\max_{x,z}z
\quad\text{subject to}\quad
z\le px+u,\quad z\le v-qx,\quad L\le x\le U.
```

For each fixed $`x`$, the largest feasible $`z`$ is exactly the smaller reward, proving equivalence. This formulation remains useful with more than two affine rewards even when a single crossing argument no longer describes the entire solution. Here the direct monotonicity proof already establishes global optimality and uniqueness.

<!-- further-distinct-problems-AN -->

### 9.9 Further proofs and constrained problems

These additional problems are constructed for these notes.

#### Problem AN-57: uniform convergence and interchange of limits

For $`f_n(x)=x^n`$ on $`[0,1]`$, find the pointwise limit and determine uniform convergence on this interval and on $`[0,a]`$ for fixed $`a\lt 1`$. Prove that a uniform limit of continuous functions is continuous. Give a separate example showing that uniform convergence of functions does not ensure convergence of their derivatives to the derivative of the limit.

Solution. For $`0\le x\lt 1`$, powers tend to zero; at one they remain one. The pointwise limit is therefore zero except for value one at the right endpoint. Its discontinuity already suggests that convergence is not uniform.

Directly, the supremum of $`\vert f_n-f\vert`$ over points below one is one: values can be made arbitrarily close to one by choosing $`x\lt 1`$ sufficiently near one. This supremum does not decrease with $`n`$. On $`[0,a]`$, however, the supremum is $`a^n\to0`$, proving uniform convergence there.

To prove the continuous-limit theorem, suppose $`f_n\to f`$ uniformly and every $`f_n`$ is continuous. Fix a point $`x_0`$ and tolerance $`\epsilon\gt 0`$. Choose one index $`N`$ for which $`\sup_x\vert f_N(x)-f(x)\vert \lt \epsilon/3`$. Continuity of this fixed $`f_N`$ supplies a neighborhood in which $`\vert f_N(x)-f_N(x_0)\vert \lt \epsilon/3`$. The triangle inequality gives

```math
\vert f(x)-f(x_0)\vert
\le \vert f(x)-f_N(x)\vert +\vert f_N(x)-f_N(x_0)\vert +\vert f_N(x_0)-f(x_0)\vert
\lt \epsilon.
```

The uniform choice of $`N`$ controls both evaluation points. Pointwise convergence would not give the same control throughout the neighborhood.

For derivatives, take $`g_n(x)=\sin(nx)/n`$ on the real line. The bound $`\vert g_n(x)\vert \le1/n`$ proves uniform convergence to zero. Yet $`g_n'(0)=1`$ for every $`n`$, while the derivative of the limit is zero. Stronger assumptions, such as suitable uniform convergence of derivatives together with a value at one point, are needed for a differentiation-under-limit theorem.

For $`x^n`$, the integrals equal $`1/(n+1)`$ and converge to zero, agreeing with the integral of the pointwise limit. Nonuniform convergence prevents a particular general theorem from being applied; it does not force every possible interchange to fail.

#### Problem AN-58: a contraction equation with an explicit error bound

Prove that $`x=\cos x`$ has exactly one solution in $`[0,1]`$ and that iteration $`x_{n+1}=\cos x_n`$ converges from every point of that interval. Derive a computable error bound from successive iterates.

Solution. The cosine maps $`[0,1]`$ into $`[\cos1,1]`$, which lies within the same interval. Its derivative magnitude is $`\sin x\le q=\sin1\lt 1`$ there. The mean value theorem gives $`\vert \cos x-\cos y\vert \le q\vert x-y\vert`$.

Consecutive differences therefore satisfy $`\vert x_{n+1}-x_n\vert \le q^n\vert x_1-x_0\vert`$. For $`m\gt n`$, summing these differences gives

```math
\vert x_m-x_n\vert
\le\sum_{k=n}^{m-1}\vert x_{k+1}-x_k\vert
\le\frac{q^n}{1-q}\vert x_1-x_0\vert .
```

The right side tends to zero, so the sequence is Cauchy. Completeness gives a limit, and closedness of the interval keeps it in $`[0,1]`$. Continuity of cosine then makes the limit a fixed point.

If $`r,s`$ were two fixed points, $`\vert r-s\vert =\vert \cos r-\cos s\vert \le q\vert r-s\vert`$. Since $`q\lt 1`$, this forces equality of the points. The root is approximately $`0.7390851332`$.

Letting $`m\to\infty`$ in the tail estimate gives the a priori error bound $`\vert x_n-r\vert \le q^n\vert x_1-x_0\vert /(1-q)`$. Repeating the argument from the latest known difference gives

```math
\vert x_n-r\vert \le\frac q{1-q}\vert x_n-x_{n-1}\vert .
```

This is a valid stopping certificate on this invariant contraction interval. An arbitrary iteration with small successive steps need not have such a certificate. Here the derivative bound and invariance are what convert a step difference into a bound on the unseen remaining tail.

#### Problem AN-59: interior critical points and boundary extrema

Find every global minimum and maximum of $`f(x)=x^3-3x`$ on $`[-2,2]`$. Classify its interior stationary points and explain why the endpoints must be checked separately.

Solution. The polynomial is continuous on a compact nonempty interval, so both extrema exist. Interior extrema must satisfy $`f'(x)=3x^2-3=0`$, giving $`x=-1,1`$. The second derivative is $`6x`$, so minus one is a strict local maximum and one is a strict local minimum.

The candidate list includes both endpoints because an endpoint optimum need not have derivative zero. Evaluate all four values:

```math
f(-2)=-2,\qquad f(-1)=2,\qquad f(1)=-2,\qquad f(2)=2.
```

The global minimum is minus two, attained at $`-2`$ and $`1`$. The global maximum is two, attained at $`-1`$ and $`2`$.

For completeness, the derivative is positive where $`\vert x\vert \gt 1`$ and negative where $`\vert x\vert \lt 1`$. Thus the function increases on $`[-2,-1]`$, decreases on $`[-1,1]`$, and increases on $`[1,2]`$. These monotonic pieces also prove that the candidate list contains every extremum.

The second-derivative test classified the interior points locally. The value comparison established which ones are global and found the tied boundary extrema. Neither step alone would provide the complete answer.

#### Problem AN-60: an absolute-value penalty

For real $`a`$ and $`\lambda\ge0`$, minimize $`F(x)=\tfrac12(x-a)^2+\lambda\vert x\vert`$. Derive the optimizer for every $`a`$ and justify the case at zero.

Solution. On $`x\gt 0`$, differentiation gives $`F'(x)=x-a+\lambda`$. Its stationary point is $`a-\lambda`$ and belongs to that region only if $`a\gt \lambda`$. On $`x\lt 0`$, the derivative is $`x-a-\lambda`$, giving $`a+\lambda`$ only when $`a\lt -\lambda`$.

At zero, the left derivative is $`-a-\lambda`$ and the right derivative is $`-a+\lambda`$. A minimum there requires the left derivative to be nonpositive and the right derivative nonnegative, equivalent to $`\vert a\vert \le\lambda`$.

Thus the solution is

```math
x^*=\begin{cases}
a-\lambda,&a\gt \lambda,\\[0pt]
0,&\vert a\vert \le\lambda,\\[0pt]
a+\lambda,&a\lt -\lambda.
\end{cases}
```

The quadratic part is strictly convex and the absolute value is convex, so their sum is strictly convex. Every candidate meeting the appropriate one-sided condition is therefore the unique global minimum.

This can also be expressed using a subgradient. At zero, the subgradients of $`\vert x\vert`$ form $`[-1,1]`$ because $`\vert y\vert \ge sy`$ for every $`y`$ exactly when $`\vert s\vert \le1`$. The optimality condition $`0\in x-a+\lambda\partial\vert x\vert`$ gives all three cases in one statement.

With $`\lambda=1`$, inputs $`a=3,1/2,-3`$ produce solutions $`2,0,-2`$. An interval of inputs is mapped exactly to zero. For a vector problem with a sum of separate squared errors and absolute penalties, the objective separates into coordinates and this formula applies to each. Correlated design columns generally destroy that coordinatewise closed form.

#### Problem AN-61: a disk constraint and an exact dual certificate

Minimize $`(x-2)^2+(y-1)^2`$ subject to $`x^2+y^2\le1`$. Derive the primal solution, all needed KKT quantities, and a matching dual value.

Solution. Write $`z=(2,1)^\top`$, whose norm is $`\sqrt5\gt 1`$. The unconstrained minimizer lies outside the disk. With multiplier $`\lambda\ge0`$, the Lagrangian is

```math
L(w,\lambda)=(1+\lambda)\Vert w\Vert ^2-2z^\top w+5-\lambda.
```

Stationarity gives $`w=z/(1+\lambda)`$. If the constraint were slack, complementary slackness would give $`\lambda=0`$ and $`w=z`$, which is infeasible. Thus the constraint binds, yielding $`1+\lambda=\sqrt5`$.

The solution is $`w^*=(2,1)^\top/\sqrt5`$, with multiplier $`\lambda^*=\sqrt5-1\gt 0`$. Its cost is the squared radial distance $`(\sqrt5-1)^2=6-2\sqrt5`$. Feasibility, multiplier nonnegativity, stationarity, and complementary slackness all hold.

Minimizing the displayed Lagrangian over unrestricted $`w`$ gives

```math
g(\lambda)=5-\lambda-\frac5{1+\lambda}.
```

Its derivative is $`-1+5/(1+\lambda)^2`$ and its second derivative is $`-10/(1+\lambda)^3\lt 0`$. The unique dual maximizer is the same multiplier, with $`g(\lambda^*)=6-2\sqrt5`$.

The matching values certify optimality. The primal objective is strictly convex and the disk is convex, so the minimizer is unique. The center of the disk is strictly feasible, also verifying Slater's condition. The multiplier is attached to the squared-norm constraint as written; using $`\Vert w\Vert -1\le0`$ would produce a differently scaled multiplier while leaving the optimizer unchanged.

#### Problem AN-62: an equality with an active nonnegativity bound

Minimize $`\tfrac12(x-3)^2+\tfrac12(y+1)^2`$ subject to $`x+y=1`$ and $`x,y\ge0`$. Solve by elimination and give a KKT certificate with explicit multiplier signs.

Solution. Eliminate $`y=1-x`$. Nonnegativity restricts $`x`$ to $`[0,1]`$. The scalar objective is $`\tfrac12(x-3)^2+\tfrac12(2-x)^2`$ with derivative $`2x-5`$. It decreases throughout that interval, so $`x^*=1`$ and $`y^*=0`$. The minimum value is $`5/2`$.

Use $`L=f+\nu(x+y-1)-\alpha x-\gamma y`$, with $`\alpha,\gamma\ge0`$ and unrestricted equality multiplier $`\nu`$. Stationarity at the proposed solution gives

```math
-2+\nu-\alpha=0,\qquad
1+\nu-\gamma=0.
```

Since $`x^*\gt 0`$, complementary slackness requires $`\alpha=0`$. Therefore $`\nu=2`$ and $`\gamma=3`$. The second nonnegativity constraint binds and has positive multiplier. All signs match the convention of writing its inequality as $`-y\le0`$.

The objective Hessian is identity and the feasible set is convex, so KKT is sufficient and the solution is unique. Existence also follows directly from continuity on the compact line segment.

For a sensitivity check, replace the equality's right side by $`b`$ near one. The active set remains $`x=b,y=0`$, giving value $`p(b)=\tfrac12(b-3)^2+\tfrac12`$. Hence $`p'(1)=-2=-\nu`$. An equality multiplier is not constrained to the sign used for inequality multipliers; its sensitivity sign comes from the term $`-\nu b`$.

#### Problem AN-63: weighted worst-case squared errors

Minimize $`\max\lbrace x^2,4(x-3)^2\rbrace`$ over all real $`x`$. Prove a global lower bound and derive the multipliers in an auxiliary-variable formulation.

Solution. Let the worst-case value at $`x`$ be $`t`$. Then $`\vert x\vert \le\sqrt t`$ and $`\vert x-3\vert \le\sqrt t/2`$. The triangle inequality gives

```math
3\le\vert x\vert +\vert x-3\vert \le\frac32\sqrt t.
```

Thus $`t\ge4`$ for every real decision. At $`x=2`$, the two squared errors both equal four, so the bound is attained and the optimum is four.

To identify uniqueness, equality requires distances two from zero and one from three, with $`x`$ between the endpoints. This forces $`x=2`$.

The epigraph problem minimizes $`t`$ subject to $`x^2-t\le0`$ and $`4(x-3)^2-t\le0`$. Its Lagrangian is

```math
L=t+\lambda_1(x^2-t)+\lambda_2(4(x-3)^2-t).
```

Stationarity in $`t`$ gives $`\lambda_1+\lambda_2=1`$. At $`x=2`$, stationarity in $`x`$ gives $`4\lambda_1-8\lambda_2=0`$. Therefore $`\lambda_1=2/3`$ and $`\lambda_2=1/3`$. They are nonnegative, both inequalities bind, and all KKT conditions hold.

The auxiliary problem is convex and has a strictly feasible point, such as $`x=2,t=5`$. The multipliers can be read as weights on the two active error functions. They are not equal because the two functions have different slopes at the balancing point. Equal worst-case losses do not imply equal dual weights.

#### Problem AN-64: a continuous value with a discontinuous optimizer

Let $`p(a)=\min_{-1\le x\le1}ax`$. Find the value and all minimizers. Then add $`\epsilon x^2/2`$ with $`\epsilon\gt 0`$ and determine how the solution and value change.

Solution. For $`a\gt 0`$ the objective increases in $`x`$, so the minimizer is minus one. For $`a\lt 0`$ it decreases, so the minimizer is one. At $`a=0`$ every feasible point is optimal. Thus $`p(a)=-\vert a\vert`$. The value is continuous and Lipschitz with constant one, but is not differentiable at zero. Every optimizing selection must jump across zero.

For the regularized objective $`ax+\epsilon x^2/2`$, the unrestricted stationary point is $`-a/\epsilon`$. Restricting to the interval gives

```math
x_\epsilon(a)=\min(1,\max(-1,-a/\epsilon)).
```

Strict convexity makes this optimizer unique. Substitution gives

```math
p_\epsilon(a)=
\begin{cases}
-a^2/(2\epsilon),&\vert a\vert \le\epsilon,\\[0pt]
\epsilon/2-\vert a\vert ,&\vert a\vert \ge\epsilon.
\end{cases}
```

The expressions agree at the boundary. Their derivatives also agree, and $`p_\epsilon'(a)=x_\epsilon(a)`$. Thus the regularized value is differentiable everywhere.

The penalty lies between zero and $`\epsilon/2`$ throughout the feasible interval, so $`0\le p_\epsilon(a)-p(a)\le\epsilon/2`$ uniformly in $`a`$. Uniform convergence of smooth value functions has produced a nonsmooth limiting value. This agrees with the earlier counterexample about derivatives: convergence of functions alone does not justify interchanging differentiation and a limit.

#### Problem AN-65: weighted arithmetic and geometric means

For positive numbers $`x_i`$ and nonnegative weights $`w_i`$ summing to one, prove $`\prod_i x_i^{w_i}\le\sum_iw_ix_i`$ and characterize equality. Apply it to maximizing a product subject to a fixed positive sum.

Solution. The logarithm has second derivative $`-1/x^2\lt 0`$ on its positive domain, so it is strictly concave. Finite Jensen applied to it gives

```math
\sum_iw_i\log x_i\le
\log\left(\sum_iw_ix_i\right).
```

Exponentiation preserves order and turns the left side into the weighted product. Equality holds exactly when all $`x_i`$ carrying positive weights are equal; zero-weight entries do not affect either side.

For $`n`$ positive numbers with sum $`B\gt 0`$, take equal weights $`1/n`$. Then $`(\prod_i x_i)^{1/n}\le B/n`$, so the largest possible product is $`(B/n)^n`$, attained only at $`x_i=B/n`$.

Existence can also be checked by extending the product continuously to the compact nonnegative simplex. Every boundary point has product zero, while the equal positive allocation has positive product. Therefore the maximum is interior, making the equality condition applicable.

For a weighted log objective $`\sum_iw_i\log x_i`$ with strictly positive weights, the same budget instead gives stationarity $`w_i/x_i=\lambda`$ and hence $`x_i=w_iB`$. This is a different allocation problem from the unweighted product. The sum constraint binds by monotonicity, strict concavity gives uniqueness, and $`\lambda=1/B`$ follows from summing the allocations.

#### Problem AN-66: comparative statics from a unique stationary point

For $`a\gt 0`$ and real $`b`$, minimize $`f(x;a,b)=x^4/4+ax^2/2-bx`$. Prove existence and uniqueness, derive the optimizer's local derivatives with respect to both parameters, and evaluate them at $`a=1,b=2`$.

Solution. The quartic leading term makes the objective coercive, so continuity gives existence. Its second derivative is $`3x^2+a\gt 0`$, proving strict convexity and uniqueness.

The first-order equation is $`x^3+ax-b=0`$. Its left side is strictly increasing in $`x`$ and tends to opposite infinities at the two ends of the real line, so it has exactly one root. The derivative in $`x`$ is nonzero. The implicit-function theorem therefore gives a smooth local optimizer map, and differentiating the equation yields

```math
\frac{\partial x^*}{\partial b}=\frac1{3(x^*)^2+a}\gt 0,\qquad
\frac{\partial x^*}{\partial a}=-\frac{x^*}{3(x^*)^2+a}.
```

At $`a=1,b=2`$ the root is one, so the derivatives are $`1/4`$ and $`-1/4`$. The minimum value there is $`1/4+1/2-2=-5/4`$.

Let $`v(a,b)=f(x^*(a,b);a,b)`$. The chain-rule terms involving derivatives of $`x^*`$ disappear because $`f_x=0`$. Thus $`v_a=(x^*)^2/2`$ and $`v_b=-x^*`$. These are envelope derivatives.

The value is concave in $`(a,b)`$ on $`a\gt 0`$ because it is the infimum over $`x`$ of affine functions of those parameters. This concerns parameter dependence of the optimized value. The original objective is strictly convex in the decision variable. Convexity in one set of variables and concavity in another are compatible.

#### Problem AN-67: curvature along a circle

Find the global extrema of $`xy`$ on $`x^2+y^2=1`$. Compute equality multipliers and test the Lagrangian Hessian along tangent directions.

Solution. The identities $`(x+y)^2=1+2xy`$ and $`(x-y)^2=1-2xy`$ imply $`-1/2\le xy\le1/2`$. The upper value occurs at $`x=y=\pm1/\sqrt2`$, and the lower at $`x=-y=\pm1/\sqrt2`$.

Use $`L=xy+\nu(x^2+y^2-1)`$. Stationarity gives $`y+2\nu x=0`$ and $`x+2\nu y=0`$. At a minimum with $`y=-x`$, $`\nu=1/2`$. At a maximum with $`y=x`$, $`\nu=-1/2`$.

The Hessian in the decision variables is

```math
\nabla^2_{x,y}L=
\begin{pmatrix}2\nu&1\\[0pt] 1&2\nu\end{pmatrix}.
```

A unit tangent direction to the circle at $`(x,y)`$ is $`h=(-y,x)`$. At either minimum, its quadratic value is two, strictly positive. At either maximum it is minus two, strictly negative.

The Hessian need not be positive definite on the entire plane at a constrained minimum. Here its other eigenvalue is zero. Only tangent directions represent first-order feasible motion, and the multiplier incorporates the circle's curvature. The direct square identities establish global extrema; the tangent Hessian independently checks their local classification.

The constraint gradient is $`2(x,y)`$ and never vanishes on the circle, so the equality constraint qualification holds. Solving a multiplier system without this check could miss or mishandle a singular feasible point in another problem.

#### Problem AN-68: a polynomial whose infimum escapes to infinity

Analyze global minimization of $`f(x,y)=(xy-1)^2+x^2`$ over $`\mathbb R^2`$. Find its infimum, determine attainment, and classify its stationary point.

Solution. Both terms are nonnegative, so the infimum is at least zero. Along $`x_n=1/n,y_n=n`$, the product equals one and the objective equals $`1/n^2`$. Thus the infimum is zero.

Attaining zero would require $`x=0`$ from the second term and $`xy=1`$ from the first, an impossibility. There is no global minimizer. The domain is closed and the objective is continuous, but neither compactness nor coercivity is available.

For stationary points,

```math
f_x=2y(xy-1)+2x,\qquad f_y=2x(xy-1).
```

The second equation implies $`x=0`$ or $`xy=1`$. The latter makes the first equation $`2x=0`$ and is impossible with $`xy=1`$. The former makes the first equation $`-2y=0`$. Thus the only stationary point is $`(0,0)`$, with value one.

Its Hessian is

```math
\begin{pmatrix}2&-2\\[0pt] -2&0\end{pmatrix}
```

, whose determinant is negative. It is a saddle. Direct paths confirm this: $`f(t,0)=1+t^2\gt 1`$, while $`f(t,t)=1-t^2+t^4\lt 1`$ for sufficiently small nonzero $`t`$.

A search for stationary points would therefore find a saddle and never find a minimizer, because none exists. The near-optimal sequence sends one coordinate to zero and the other to infinity. Reporting the numerical location of a small objective would need to acknowledge that this sequence does not approach a finite optimizer.

[Contents](README.md) · [Part 1](02-analysis-and-optimization-part-1.md)
