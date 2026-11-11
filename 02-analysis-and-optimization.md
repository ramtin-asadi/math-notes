# Analysis and Optimization

The starting point is the real number system. Later sections use Euclidean spaces with $\|x\|=(\sum_i x_i^2)^{1/2}$. A statement about all points of a set includes its boundary unless an interior is explicitly specified.

Sources:

- Jiří Lebl, [Basic Analysis I](https://www.jirka.org/ra/realanal.pdf) and [Basic Analysis II](https://www.jirka.org/ra/realanal2.pdf), version 6.3.
- [MIT 18.100A, Real Analysis, Fall 2020](https://ocw.mit.edu/courses/18-100a-real-analysis-fall-2020/), taught by Casey Rodriguez.
- Stephen Boyd and Lieven Vandenberghe, [Convex Optimization](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf), Chapters 2–5 and 9.

Source exercise statements are summarized; the linked texts contain the originals. Examples and the additional problems in Section 8 are constructed here.

Contents:

1. [Reasoning and completeness](#1-reasoning-and-completeness)
2. [Sequences and convergence](#2-sequences-and-convergence)
3. [Continuity and compactness](#3-continuity-and-compactness)
4. [Differentiability and Taylor approximation](#4-differentiability-and-taylor-approximation)
5. [Convex sets and functions](#5-convex-sets-and-functions)
6. [Constrained optimization](#6-constrained-optimization)
7. [Lagrangian duality](#7-lagrangian-duality)
8. [Optimization with uncertainty](#8-optimization-with-uncertainty)
9. [Further worked problems](#9-further-worked-problems)

## 1. Reasoning and completeness

### 1.1 Theory

The order of quantifiers changes a statement. The assertion $\forall x\,\exists y:y>x$ is true over $\mathbb R$: choose $y=x+1$. The assertion $\exists y\,\forall x:y>x$ is false: choosing $x=y$ contradicts it. To negate $\forall x\,\exists y\,P(x,y)$, write $\exists x\,\forall y\,\neg P(x,y)$.

A direct proof starts from assumptions and derives the conclusion. A contradiction proof assumes the conclusion fails and derives an impossibility. Induction proves a base case and then proves that the statement for $n$ implies the statement for $n+1$. In each method, the assumptions must stay fixed during the argument.

An upper bound of $A\subset\mathbb R$ is a number $u$ with $a\le u$ for all $a\in A$. The supremum is the least upper bound. Completeness of $\mathbb R$ states that every nonempty set bounded above has a real supremum. A maximum must belong to the set; a supremum need not. For $A=(0,1)$, $\sup A=1$ but there is no maximum.

For $s=\sup A$, two conditions are useful:

$$
a\le s\quad(a\in A),\qquad
\text{for each }\varepsilon>0\text{ there is }a\in A
\text{ with }s-\varepsilon<a.
$$

The second condition prevents a smaller number from being an upper bound.

### 1.2 Derivations and proofs

#### A supremum argument

Let $A$ be nonempty and bounded above, with supremum $s$. For a fixed $c\in\mathbb R$, prove $\sup(A+c)=s+c$.

Every $a+c\le s+c$, so $s+c$ is an upper bound. Given $\varepsilon>0$, choose $a>s-\varepsilon$. Then $a+c>s+c-\varepsilon$, so no number below $s+c$ is an upper bound. This proves the claim without assuming the supremum belongs to $A$.

#### Bernoulli's inequality

For $x>-1$ and an integer $n\ge1$,

$$
(1+x)^n\ge1+nx.
$$

The case $n=1$ is equality. If the result holds for $n$, multiplication by the positive number $1+x$ preserves the inequality:

$$
(1+x)^{n+1}\ge(1+nx)(1+x)
=1+(n+1)x+nx^2\ge1+(n+1)x.
$$

This proves the induction step. Positivity of the multiplier is the reason for the assumption on $x$.

#### Completeness, large integers, and rational approximations

The completeness axiom concerns arbitrary bounded sets of real numbers. Some familiar properties of integers and rational numbers can be derived from it, so they need not be silently assumed in a proof about limits.

First, the natural numbers have no real upper bound. Suppose they did, and let $s=\sup\mathbb N$. The number $s-1$ is smaller than the least upper bound, so it cannot itself be an upper bound. There is an integer $n>s-1$. Then the integer $n+1>s$ contradicts the definition of $s$. This proves the Archimedean property: given any real number, some integer is larger.

It follows that for every $\epsilon>0$ there is an integer $n$ with $1/n<\epsilon$. The reciprocal is taken only after choosing a positive integer, so reversing the inequality is justified.

Now take $a<b$. Choose a positive integer $n$ large enough that $1/n<b-a$, and let $m=\lfloor na\rfloor+1$. The floor property gives $na<m\le na+1$. Dividing by $n$ yields

$$
a<\frac mn\le a+\frac1n<b.
$$

Thus every nonempty open real interval contains a rational number. This is density of the rationals. Density says that rational points approximate real points; completeness says that certain limiting points actually exist within the ambient system. A dense subset need not be complete.

#### Constructing a square root using a supremum

Let $A=\{x\in\mathbb R:x\ge0,\ x^2\le2\}$. The set contains one and is bounded above by two. Completeness supplies $s=\sup A$, with $1\le s\le2$. We show $s^2=2$ without assuming the square root already exists.

If $s^2<2$, choose $0<h<1$ so small that $h<(2-s^2)/(2s+1)$. Then

$$
(s+h)^2=s^2+2sh+h^2
\le s^2+(2s+1)h<2.
$$

So $s+h\in A$, contradicting that $s$ is an upper bound.

If $s^2>2$, choose $0<h<\min(s/2,(s^2-2)/(2s))$. Then $s-h>0$ and

$$
(s-h)^2=s^2-2sh+h^2>2.
$$

Every $x\in A$ must satisfy $x<s-h$, since otherwise its square would exceed two. Thus $s-h$ is an upper bound smaller than $s$, another contradiction. The only remaining possibility is $s^2=2$.

Uniqueness among nonnegative numbers follows from $u^2-v^2=(u-v)(u+v)$: if two positive numbers have the same square, their difference is zero. The supremum argument established existence; the factorization established uniqueness.

#### Example AN-E17: a supremum with a quantitative witness

Let $A=\{n/(n+1):n=1,2,\ldots\}$. Determine its infimum and supremum, and identify which bounds are attained.

Solution. Every term is $1-1/(n+1)$, so it is below one. For a proposed gap $\epsilon>0$, choose $n>1/\epsilon-1$. Then $n/(n+1)>1-\epsilon$. Therefore no number below one bounds the set, and $\sup A=1$. No term equals one.

The consecutive difference is $1/[(n+1)(n+2)]>0$, so the smallest term is the first, $1/2$. Hence $\inf A=\min A=1/2$. An infimum can be attained even when the supremum of the same set is not.

### 1.3 Worked examples

#### Example AN-E01: choosing witnesses

Prove that for every $\varepsilon>0$ there is $N$ with $1/n<\varepsilon$ whenever $n\ge N$.

Solution. The Archimedean property gives an integer $N>1/\varepsilon$. For $n\ge N$, $1/n\le1/N<\varepsilon$. The choice of $N$ depends on $\varepsilon$ but must work for every later $n$. Choosing a new $N$ separately for each $n$ would not establish the quantified statement.

#### Example AN-E02: a supremum outside the rational numbers

Let $A=\{q\in\mathbb Q:q\ge0,\ q^2<2\}$. Find its supremum in $\mathbb R$.

Solution. It is $\sqrt2$. Every element is below $\sqrt2$. For any $\varepsilon>0$, density of $\mathbb Q$ gives a rational $q$ between $\max(0,\sqrt2-\varepsilon)$ and $\sqrt2$, so $q\in A$ lies within $\varepsilon$ of that upper bound. In $\mathbb Q$, this set has no rational supremum. Completeness is a property of the ambient number system.

### 1.4 Selected exercises

#### Problem AN-01: an upper bound in the set

Source: [Lebl I, Exercise 1.1.5](https://www.jirka.org/ra/realanal.pdf).

Problem. If a set has an upper bound $b$ that also belongs to the set, prove that $b$ is its supremum and maximum.

Solution. Suppose $b$ is an upper bound of $A$ and $b\in A$. Any other upper bound $u$ must satisfy $b\le u$, because it bounds every member, including $b$. Thus $b$ is the least upper bound and also the maximum. The proof works in any ordered set where the stated comparisons make sense.

#### Problem AN-02: the arithmetic-geometric mean inequality

Source: [Lebl I, Exercise 1.2.7](https://www.jirka.org/ra/realanal.pdf).

Problem. For positive $x,y$, prove the arithmetic-geometric mean inequality $\sqrt{xy}\le(x+y)/2$ and determine its equality condition.

Solution. For $x,y>0$, $(\sqrt x-\sqrt y)^2\ge0$. Expanding gives $x+y\ge2\sqrt{xy}$, or

$$
\sqrt{xy}\le\frac{x+y}{2}.
$$

Equality holds exactly when $\sqrt x=\sqrt y$, hence when $x=y$. For a fixed sum, this also proves that the product is largest when the two terms are equal.

#### Problem AN-03: bounds for sums of sets

Source: [Lebl I, Exercise 1.2.9](https://www.jirka.org/ra/realanal.pdf).

Problem. For nonempty bounded sets $A,B\subset\mathbb R$, determine the infimum and supremum of $\{a+b:a\in A,b\in B\}$ and prove the formulas.

Solution. For nonempty bounded $A,B$, put $C=\{a+b:a\in A,b\in B\}$. If $\alpha=\inf A$, $\beta=\inf B$, $s=\sup A$, and $t=\sup B$, then $\alpha+\beta\le c\le s+t$ for every $c\in C$, proving boundedness.

Given $\varepsilon>0$, choose $a>s-\varepsilon/2$ and $b>t-\varepsilon/2$. Then $a+b>s+t-\varepsilon$, proving $\sup C=s+t$. Choose instead $a<\alpha+\varepsilon/2$ and $b<\beta+\varepsilon/2$ to obtain $\inf C=\alpha+\beta$. Approximation to the bounds replaces any unjustified assumption that they are attained.

## 2. Sequences and convergence

### 2.1 Theory

A sequence $x_n$ converges to $L$ if

$$
\forall\varepsilon>0\ \exists N\ \forall n\ge N:
|x_n-L|<\varepsilon.
$$

It is Cauchy if its late terms are close to each other:

$$
\forall\varepsilon>0\ \exists N\ \forall m,n\ge N:
|x_n-x_m|<\varepsilon.
$$

Every convergent sequence is bounded and Cauchy. In $\mathbb R$, every Cauchy sequence converges. In an incomplete space this implication can fail: rational approximations to $\sqrt2$ are Cauchy in $\mathbb Q$ but have no rational limit.

A subsequence $x_{n_k}$ uses strictly increasing indices $n_k$. Every subsequence of a convergent sequence has the same limit. Bolzano-Weierstrass states that every bounded real sequence has a convergent subsequence. Boundedness alone does not make the whole sequence converge.

### 2.2 Derivations and proofs

#### Convergence implies the Cauchy property

Suppose $x_n\to L$. Given $\varepsilon>0$, choose $N$ so that $|x_n-L|<\varepsilon/2$ for $n\ge N$. For $m,n\ge N$,

$$
|x_n-x_m|\le|x_n-L|+|x_m-L|<\varepsilon.
$$

The same $N$ controls both terms.

#### Monotone convergence from completeness

If $x_n$ is increasing and bounded above, let $s=\sup\{x_n:n\ge1\}$. For any $\varepsilon>0$, there is $N$ with $x_N>s-\varepsilon$. Monotonicity then gives $s-\varepsilon<x_N\le x_n\le s$ for all $n\ge N$. Therefore $x_n\to s$.

#### Bolzano-Weierstrass

Enclose a bounded sequence in a closed interval. Split that interval into two halves; at least one half contains infinitely many sequence terms. Keep such a half and repeat. The resulting nested closed intervals have lengths tending to zero. Completeness gives a unique common point $L$. Choose successively a term in each interval whose index exceeds all previously chosen indices. Its distance from $L$ is at most the interval length, so the selected subsequence converges.

To obtain convergence of a real Cauchy sequence, first prove it bounded by comparing its tail with one fixed late term. Bolzano-Weierstrass gives a convergent subsequence $x_{n_k}\to L$. For late $n$, compare $x_n$ with a late subsequence term and then with $L$; the triangle inequality proves $x_n\to L$.

#### Identifying a limit does not prove convergence

A common calculation begins with a recursion, assumes a limit exists, and solves an equation for that limit. This is a useful final step. It is not an existence proof. For $x_{n+1}=1-x_n$ with $x_1=0$, assuming convergence would give $L=1-L$ and hence $L=1/2$. The actual sequence alternates between zero and one.

A complete recursive-limit argument usually establishes three separate facts. The recursion stays in a chosen interval, its values move monotonically or satisfy a Cauchy estimate, and the defining function is continuous at a possible limit. The bounded monotone convergence theorem then provides the limit before continuity is used to identify it.

There is another useful criterion. Suppose a bounded real sequence has only one possible subsequential limit, $L$. Then the full sequence converges to $L$. If it did not, some $\epsilon_0>0$ would allow infinitely many terms outside $(L-\epsilon_0,L+\epsilon_0)$. Select those terms as a subsequence. Boundedness gives a convergent further subsequence, whose limit remains at least $\epsilon_0$ away from $L$. This contradicts uniqueness of the subsequential limit. Boundedness is essential to extract the further convergent subsequence.

#### Series as sequences of partial sums

An infinite sum $\sum_{k=1}^\infty a_k$ is defined as the limit of partial sums $S_n=\sum_{k=1}^n a_k$. The symbols do not assert that the limit exists. Since $S_m-S_n=\sum_{k=n+1}^m a_k$, completeness gives the Cauchy criterion

$$
\forall\epsilon>0\ \exists N\ \forall m>n\ge N:
\left|\sum_{k=n+1}^m a_k\right|<\epsilon.
$$

If $\sum |a_k|$ converges, its tails satisfy this criterion. The triangle inequality bounds every signed tail by the corresponding absolute tail, proving convergence of $\sum a_k$. This is the implication from absolute convergence to convergence.

The condition $a_n\to0$ is necessary: if $S_n\to S$ then $a_n=S_n-S_{n-1}\to0$. It is not sufficient, as the harmonic series shows. For every integer $j\ge0$,

$$
\sum_{k=2^j+1}^{2^{j+1}}\frac1k
\ge \frac{2^j}{2^{j+1}}=\frac12.
$$

Arbitrarily late tails therefore contain blocks with sum at least one half. The partial sums fail the Cauchy criterion even though the individual summands vanish.

#### Example AN-E18: Newton's recursion for a square root

Let $x_0=2$ and $x_{n+1}=(x_n+2/x_n)/2$. Prove convergence and explain the rapid improvement.

Solution. Put $r=\sqrt2$, whose existence was established above. For every positive $x_n$,

$$
x_{n+1}-r=\frac{(x_n-r)^2}{2x_n}\ge0.
$$

Thus every next iterate is at least $r$. If $x_n\ge r$, then $2/x_n\le x_n$, so $x_{n+1}\le x_n$. The sequence is decreasing and bounded below, hence converges. Taking limits in the continuous update gives $L=(L+2/L)/2$, so $L^2=2$ and positivity gives $L=r$.

The same identity describes the error, rather than only its limit. Near $r$, the next error is approximately the square of the current error divided by $2r$. The first iterates are $2$, $3/2$, $17/12$, and $577/408$. Once the error is small, squaring it roughly doubles the number of correct digits. The proof used positivity to keep the denominator away from zero; a similar-looking recursion without that control might fail.

### 2.3 Worked examples

#### Example AN-E03: bounded oscillation

For $x_n=(-1)^n$, the even subsequence equals $1$ and the odd subsequence equals $-1$. A convergent full sequence would force both subsequences to share its limit, which is impossible. Also late adjacent terms remain distance two apart, so the sequence is not Cauchy.

#### Example AN-E04: recursive convergence

Let $x_1=1$ and $x_{n+1}=\sqrt{2+x_n}$.

Solution. If $1\le x_n\le2$, then $1\le x_{n+1}\le2$. Also $2+x_n-x_n^2=(2-x_n)(x_n+1)\ge0$, so $x_{n+1}\ge x_n$. Induction proves bounded monotonicity. Let the limit be $L\in[1,2]$. Continuity of the square root gives $L=\sqrt{2+L}$, so $(L-2)(L+1)=0$. The interval excludes $-1$, leaving $L=2$.

### 2.4 Selected exercises

#### Problem AN-04: decreasing monotone convergence

Source: [Lebl I, Exercise 2.1.11](https://www.jirka.org/ra/realanal.pdf).

Problem. Prove that a decreasing real sequence bounded below converges, using completeness.

Solution. Let $x_n$ decrease and be bounded below. Set $s=\inf\{x_n\}$. For $\varepsilon>0$, choose $N$ with $x_N<s+\varepsilon$. Then $s\le x_n\le x_N<s+\varepsilon$ for $n\ge N$. Hence $|x_n-s|<\varepsilon$ and $x_n\to s$. Equivalently, apply increasing monotone convergence to $-x_n$.

#### Problem AN-05: convergent subsequences of an alternating sequence

Source: [Lebl I, Exercise 2.1.14](https://www.jirka.org/ra/realanal.pdf).

Problem. For $x_n=(-1)^n$, exhibit convergent subsequences and determine their possible limits.

Solution. Set $n_k=2k$. Then $(-1)^{n_k}=1$ for every $k$, so this subsequence converges to $1$. Taking $n_k=2k-1$ produces a subsequence converging to $-1$. A subsequence with infinitely many indices of both parities cannot converge, because its tail contains terms distance two apart.

#### Problem AN-06: the geometric series

Source: [Lebl I, Exercise 2.5.2](https://www.jirka.org/ra/realanal.pdf).

Problem. Prove convergence of the geometric series for $|r|<1$, find its sum, and bound its remaining tail.

Solution. For $|r|<1$, the finite identity

$$
(1-r)\sum_{n=0}^N r^n=1-r^{N+1}
$$

follows by cancellation of neighboring powers. Since $r^{N+1}\to0$,

$$
\sum_{n=0}^\infty r^n=\frac1{1-r}.
$$

The tail after $N$ is $r^{N+1}/(1-r)$, with absolute value at most $|r|^{N+1}/(1-|r|)$. This quantitative remainder is useful in contraction and discounting arguments.

## 3. Continuity and compactness

### 3.1 Theory

Continuity at $a$ means that for every $\varepsilon>0$ there exists $\delta>0$ such that $|x-a|<\delta$ within the domain implies $|f(x)-f(a)|<\varepsilon$. In metric spaces, replace absolute values by the domain and codomain distances. Sequential continuity, $x_n\to a\Rightarrow f(x_n)\to f(a)$, is equivalent.

A metric satisfies nonnegativity, identity of indiscernibles, symmetry, and the triangle inequality. An open ball is $B(x,r)=\{y:d(x,y)<r\}$. A set is closed when it contains limits of its convergent sequences. It is compact when every open cover has a finite subcover; in metric spaces this is equivalent to every sequence having a subsequence converging to a point of the set.

In $\mathbb R^n$, compactness is equivalent to closedness and boundedness. This characterization does not hold in every metric space. A continuous real function on a nonempty compact set attains its minimum and maximum and is uniformly continuous.

### 3.2 Derivations and proofs

#### Attainment of extrema

Let $f:K\to\mathbb R$ be continuous and $K$ nonempty compact. First $f$ is bounded: if not, there are $x_n\in K$ with $|f(x_n)|>n$. Compactness gives a subsequence tending to $x\in K$, while continuity gives $f(x_{n_k})\to f(x)$, contradicting unbounded growth.

Now $M=\sup_{x\in K}f(x)$ is finite. Choose $x_n$ with $f(x_n)>M-1/n$. A convergent subsequence has limit $x^*\in K$, and continuity gives $f(x^*)=M$. Apply the same argument to $-f$ for the minimum.

#### Uniform continuity on a compact set

If uniform continuity fails, some $\varepsilon_0>0$ permits pairs $x_n,y_n\in K$ with $\|x_n-y_n\|<1/n$ but $|f(x_n)-f(y_n)|\ge\varepsilon_0$. Extract $x_{n_k}\to x$. The shrinking pairwise distance implies $y_{n_k}\to x$ too. Continuity at $x$ forces both function values to converge to $f(x)$, a contradiction.

#### Sequential compactness in Euclidean space

The closed-and-bounded criterion has two ingredients. Boundedness prevents a sequence from escaping to infinity. Closedness prevents the limit of a selected subsequence from escaping through a missing boundary point.

For a bounded sequence $x_n\in\mathbb R^d$, its first coordinate is a bounded real sequence. Select a subsequence whose first coordinate converges. From that subsequence select a further one whose second coordinate converges. Continue through the finite list of coordinates. On the final subsequence every coordinate converges, and therefore the Euclidean norm of the difference from the coordinatewise limit tends to zero. If the original set is closed, that limit belongs to it.

This proves the sequential form of compactness for a closed bounded subset of $\mathbb R^d$. The finite dimension is used when the extraction ends after finitely many coordinates. The equivalence with the open-cover definition is a metric-space theorem; the elementary coordinate argument is enough for the optimization arguments here.

A compact subset of Euclidean space must be bounded, since the increasing balls $B(0,n)$ cover the space and a finite subcover gives one finite radius. It must also be closed: a convergent sequence in the set has a subsequence converging to a point in the set by compactness, and uniqueness of limits makes that point the original limit.

#### Continuous images and minimization

If $K$ is compact and $f$ continuous, every sequence $f(x_n)$ has a subsequence $f(x_{n_k})$ converging to $f(x)$ for some $x\in K$. This proves sequential compactness of the image. In particular, a continuous real image is closed and bounded. Its finite supremum and infimum belong to the image, which is another proof of attainment of extrema.

The domain of an optimization problem need not itself be bounded. Suppose $f:\mathbb R^d\to\mathbb R$ is continuous and coercive:

$$
\|x\|\to\infty\quad\Longrightarrow\quad f(x)\to+\infty.
$$

Fix a reference point $x_0$. Coercivity gives a radius $R$ such that points outside the ball have value greater than $f(x_0)$. Enlarge the radius to include $x_0$. A minimum on the closed ball exists by compactness and is a global minimum, since no exterior point improves on $x_0$.

These are sufficient conditions, not necessary ones. A constant function attains its minimum on an unbounded domain without being coercive. Strict convexity addresses uniqueness if a minimum exists, while coercivity or compactness can supply existence.

#### Example AN-E19: existence before solving the derivative

Minimize $f(x)=x^2+e^{-x}$ over all real $x$.

Solution. The function is continuous and $f(x)\ge x^2$, so it is coercive and has a global minimizer. Its derivatives are

$$
f'(x)=2x-e^{-x},\qquad f''(x)=2+e^{-x}>0.
$$

Thus it is strictly convex and the minimizer is unique. Since $f'(0)=-1$ and $f'(1)=2-e^{-1}>0$, continuity of the derivative gives a root in $(0,1)$. The positive second derivative makes the root unique and proves that it is the minimizer. Numerically it is about $0.351734$.

No closed-form root is needed for the existence or uniqueness conclusions. A numerical root finder only locates a point whose mathematical role has already been established. Conversely, a numerical routine returning a plausible number would not by itself rule out another minimum elsewhere.

### 3.3 Worked examples

#### Example AN-E05: a maximum that escapes the domain

On $(0,1)$, $f(x)=x$ has supremum one but no maximum. Closing the interval adds the missing point and gives a maximum at $1$. On $\mathbb R$, $f(x)=e^{-x}$ has infimum zero but no minimum. Continuity without a compact domain does not ensure attainment.

#### Example AN-E06: strict convexity without existence

The function $f(x)=e^x$ is strictly convex on $\mathbb R$, since $f''(x)>0$. Its infimum is zero, approached as $x\to-\infty$, and it has no minimizer. In contrast, a constant function on $[0,1]$ attains a minimum everywhere. Existence and uniqueness require separate arguments.

### 3.4 Selected exercises

#### Problem AN-07: continuity of a square

Source: [Lebl I, Exercise 3.2.1](https://www.jirka.org/ra/realanal.pdf).

Problem. Use the epsilon-delta definition to prove that $f(x)=x^2$ is continuous at every real point.

Solution. Fix $a\in\mathbb R$ and $\varepsilon>0$. Choose $\delta=\min(1,\varepsilon/(2|a|+1))$. If $|x-a|<\delta$, then $|x|\le|a|+1$, so

$$
|x^2-a^2|=|x-a||x+a|
<\delta(2|a|+1)\le\varepsilon.
$$

This proves continuity at arbitrary $a$, hence everywhere. Bounding the factor $|x+a|$ locally is the step that allows a fixed $\delta$ to work.

#### Problem AN-08: a bounded function with neither extremum

Source: [Lebl I, Exercise 3.3.2](https://www.jirka.org/ra/realanal.pdf).

Problem. Construct a bounded function on $[0,1]$ that attains neither a maximum nor a minimum. Explain which extreme-value hypothesis fails.

Solution. Define $f(x)=x$ for $0<x<1$, and $f(0)=f(1)=1/2$. Its image is exactly $(0,1)$, which is bounded and has neither a largest nor a smallest member. The domain is compact, but the function is discontinuous at both endpoints, so the extreme-value theorem does not apply.

#### Problem AN-09: zero endpoint limits

Source: [Lebl I, Exercise 3.3.3](https://www.jirka.org/ra/realanal.pdf).

Problem. Let $f$ be continuous on $(0,1)$ and have limit zero at both endpoints. Prove that it attains at least one of its two extrema in the open interval, and determine whether both must be attained there.

Solution. Extend $f:(0,1)\to\mathbb R$ continuously by $f(0)=f(1)=0$. The extension attains both extrema on $[0,1]$. If $f$ is identically zero, every interior point attains both. Otherwise choose an interior point with nonzero value. If that value is positive, the global maximum is positive and cannot occur only at the zero endpoints, so it occurs in $(0,1)$. A negative value similarly forces an interior minimum. Both extrema need not be attained in the open interval: $x(1-x)$ has an interior maximum and infimum zero reached only at the excluded endpoints.

## 4. Differentiability and Taylor approximation

### 4.1 Theory

Differentiability at $x$ means the function has a linear approximation with a smaller-order remainder:

$$
f(x+h)=f(x)+Df(x)h+o(\|h\|).
$$

For scalar-valued functions, $Df(x)h=\nabla f(x)^\top h$. The existence of all partial derivatives at a point does not alone guarantee differentiability there. Continuous first partial derivatives in a neighborhood do.

The mean value theorem says that if $f$ is continuous on $[a,b]$ and differentiable inside, then $f(b)-f(a)=f'(c)(b-a)$ for some $c\in(a,b)$. It follows from Rolle's theorem applied to the difference between $f$ and its secant line.

At an unconstrained interior local minimum of a differentiable function, the gradient is zero. If the function is twice differentiable, its Hessian there is PSD. A PD Hessian at a stationary point is sufficient for a strict local minimum; a singular PSD Hessian requires further analysis.

### 4.2 Derivations and proofs

#### Multivariable Taylor expansion along a line

Let $f$ be twice continuously differentiable near the segment from $x$ to $x+h$, and set $\phi(t)=f(x+th)$. The chain rule gives $\phi'(t)=\nabla f(x+th)^\top h$ and $\phi''(t)=h^\top\nabla^2f(x+th)h$. Integrating twice,

$$
f(x+h)=f(x)+\nabla f(x)^\top h+
\int_0^1(1-t)h^\top\nabla^2f(x+th)h\,dt.
$$

Continuity of the Hessian gives

$$
f(x+h)=f(x)+\nabla f(x)^\top h+
\tfrac12h^\top\nabla^2f(x)h+o(\|h\|^2).
$$

If $\|\nabla^2f\|\le M$ along the segment, the first-order remainder is bounded by $M\|h\|^2/2$. If the Hessian is Lipschitz with constant $L$, the second-order remainder is bounded by $L\|h\|^3/6$.

For the necessary conditions, fix any direction $v$ and consider $f(x+tv)$. A local minimum gives derivative zero and second derivative nonnegative at $t=0$. Since $v$ was arbitrary, $\nabla f=0$ and $v^\top\nabla^2f\,v\ge0$.

#### Differentiability is a linear approximation

At $a\in\mathbb R^d$, differentiability means that there is a linear map $L$ with

$$
f(a+h)=f(a)+Lh+r(h),\qquad
\frac{|r(h)|}{\|h\|}\longrightarrow0.
$$

For a scalar-valued differentiable function, $Lh=\nabla f(a)^\top h$. The remainder condition must hold as $h$ approaches zero in every direction, with no direction chosen in advance.

Existence of coordinate partial derivatives alone does not prove this. Define $f(x,y)=xy/\sqrt{x^2+y^2}$ away from zero and $f(0,0)=0$. Both coordinate partial derivatives at zero are zero because the function vanishes on each coordinate axis. It is continuous there, since $|xy|\le(x^2+y^2)/2$ gives $|f(x,y)|\le\sqrt{x^2+y^2}/2$.

If it were differentiable at zero, its derivative would be zero. Along $h=(t,t)$, however,

$$
\frac{|f(t,t)|}{\sqrt{2t^2}}=\frac12.
$$

The required ratio does not tend to zero. The function is continuous and has both partial derivatives at the point, yet lacks a linear approximation there. Continuous partial derivatives in a neighborhood are a standard sufficient condition for differentiability.

#### Mean-value estimates and integral remainders

For a real function continuous on $[a,b]$ and differentiable inside, subtract the line through its endpoint values. The resulting function has equal endpoint values, so Rolle's theorem supplies a point with derivative zero. Undoing the subtraction proves the mean value theorem:

$$
f(b)-f(a)=f'(\xi)(b-a)
$$

for some interior $\xi$. In several variables, apply the one-dimensional theorem to $g(t)=f(a+th)$. A gradient bound $\|\nabla f\|\le M$ throughout the segment implies $|f(a+h)-f(a)|\le M\|h\|$.

When $f$ is twice continuously differentiable along that segment, two integrations of the fundamental theorem of calculus give

$$
f(a+h)=f(a)+\nabla f(a)^\top h+
\int_0^1(1-t)h^\top\nabla^2f(a+th)h\,dt.
$$

If the Hessian is Lipschitz with constant $L$ on the segment, subtracting $\tfrac12h^\top\nabla^2f(a)h$ bounds the remaining error by

$$
\int_0^1(1-t)Lt\|h\|^3dt=\frac L6\|h\|^3.
$$

This shows where both the power of the error and its constant come from.

#### Example AN-E20: a mixed two-variable approximation

Expand $f(x,y)=e^x\cos y$ about zero to degree two.

Solution. At zero the value is one, the gradient is $(1,0)$, and the Hessian is $\operatorname{diag}(1,-1)$. The mixed derivative is $-e^x\sin y$ and vanishes there. Thus

$$
P_2(x,y)=1+x+\frac{x^2-y^2}{2}.
$$

At $(0.1,0.2)$ this gives $1.085$. To bound the remainder, write $g(t)=e^{tx}\cos(ty)=\operatorname{Re}e^{t(x+iy)}$. Its third derivative has magnitude at most $e^{|x|}(x^2+y^2)^{3/2}$ for $0\le t\le1$. The one-dimensional Taylor theorem therefore bounds the error by $e^{|x|}(x^2+y^2)^{3/2}/6$, about $0.002060$ at this point. This bound covers the approximation error without treating the two coordinates as independent one-variable expansions.

### 4.3 Worked examples

#### Example AN-E07: a two-variable Taylor approximation

For $f(x,y)=e^x\log(1+y)$ near $(0,0)$, compute the quadratic approximation.

Solution. At the origin, $f=0$, $\nabla f=(0,1)$, and

$$
\nabla^2f(0,0)=\begin{pmatrix}0&1\\1&-1\end{pmatrix}.
$$

Therefore $f(x,y)=y+xy-y^2/2+o(x^2+y^2)$. The domain restriction $y>-1$ is satisfied in a sufficiently small neighborhood. At $(0.01,0.02)$, the approximation is $0.02+0.0002-0.0002=0.02$.

#### Example AN-E08: a zero Hessian

At zero, both $f(x)=x^4$ and $g(x)=-x^4$ have first and second derivatives zero. The first has a strict minimum and the second a strict maximum. A semidefinite second-order test is inconclusive when it has zero directions.

### 4.4 Selected exercises

#### Problem AN-10: the product rule

Source: [Lebl I, Exercise 4.1.1](https://www.jirka.org/ra/realanal.pdf).

Problem. Starting from difference quotients, prove the product rule for functions differentiable at a common point.

Solution. For $h\ne0$, add and subtract $f(a+h)g(a)$:

$$
\frac{f(a+h)g(a+h)-f(a)g(a)}h
=f(a+h)\frac{g(a+h)-g(a)}h+
g(a)\frac{f(a+h)-f(a)}h.
$$

Differentiability implies continuity, so $f(a+h)\to f(a)$. Taking limits yields $(fg)'(a)=f(a)g'(a)+g(a)f'(a)$. Each limit exists under the assumption that both functions are differentiable at $a$.

#### Problem AN-11: the quotient rule

Source: [Lebl I, Exercise 4.1.2](https://www.jirka.org/ra/realanal.pdf).

Problem. Prove the reciprocal and quotient differentiation rules, including the nonzero-denominator requirement.

Solution. Assume $g(a)\ne0$. Continuity makes $g$ nonzero near $a$. Direct subtraction gives

$$
\frac{1/g(a+h)-1/g(a)}h
=-\frac{g(a+h)-g(a)}{h\,g(a+h)g(a)}
\longrightarrow-\frac{g'(a)}{g(a)^2}.
$$

Apply the product rule to $f(1/g)$:

$$
\left(\frac fg\right)'(a)
=\frac{f'(a)g(a)-f(a)g'(a)}{g(a)^2}.
$$

The nonzero-denominator assumption is part of both the domain and the proof.

## 5. Convex sets and functions

### 5.1 Theory

A set $C$ is convex if it contains the segment between any two of its points. Its convex hull consists of finite convex combinations. Halfspaces, affine spaces, norm balls, and intersections of convex sets are convex.

A function on a convex domain is convex when

$$
f(\theta x+(1-\theta)y)
\le\theta f(x)+(1-\theta)f(y),\qquad0\le\theta\le1.
$$

Strict convexity replaces the inequality by a strict one for $x\ne y$ and $0<\theta<1$. The epigraph $\{(x,t):t\ge f(x)\}$ is convex exactly when $f$ is convex.

A differentiable convex function lies above every tangent plane:

$$
f(y)\ge f(x)+\nabla f(x)^\top(y-x).
$$

This follows by applying the convexity inequality at $x+t(y-x)$, dividing by $t>0$, and taking $t\downarrow0$. It makes every stationary point a global minimizer.

A quasiconvex function has convex sublevel sets $\{x:f(x)\le a\}$. Convexity implies quasiconvexity; the converse fails. For example, $x^3$ has interval sublevel sets on $\mathbb R$, but its second derivative changes sign.

### 5.2 Derivations and proofs

#### Finite Jensen inequality

For nonnegative $\theta_i$ summing to one,

$$
f\left(\sum_{i=1}^k\theta_ix_i\right)
\le\sum_{i=1}^k\theta_if(x_i).
$$

Prove this by induction. The two-point case is the definition. For the next case, let $s=1-\theta_k$. If $s=0$, the claim is equality. Otherwise write the argument as $s\sum_{i<k}(\theta_i/s)x_i+\theta_kx_k$. Apply two-point convexity, then the induction hypothesis to the normalized first $k-1$ weights. Zero weights cause no difficulty.

#### First-order and second-order convexity tests

Let $f$ be differentiable on an open convex set. If $f$ is convex, apply the convexity inequality at $x+t(y-x)$, subtract $f(x)$, divide by $t>0$, and let $t$ decrease to zero. This gives

$$
f(y)\ge f(x)+\nabla f(x)^\top(y-x).
$$

Thus each tangent affine function lies below the graph. Conversely, assume this inequality holds at every pair of points. Set $z=(1-t)x+ty$ and apply the inequality at $z$ first to $x$ and then to $y$. Multiply by $1-t$ and $t$ and add. The gradient terms cancel because $(1-t)(x-z)+t(y-z)=0$. The result is $f(z)\le(1-t)f(x)+tf(y)$, proving convexity.

If $f$ is twice continuously differentiable, restriction to every line has second derivative $h^\top\nabla^2f(x)h$. Nonnegativity for every direction is equivalent to a positive-semidefinite Hessian and gives convexity on a convex domain. Domain convexity is needed to keep the line segment within the function's domain.

A uniform bound $\nabla^2f(x)\succeq\mu I$ with $\mu>0$ gives the stronger inequality

$$
f(y)\ge f(x)+\nabla f(x)^\top(y-x)+\frac\mu2\|y-x\|^2.
$$

This is strong convexity. It supplies a quadratic lower bound on error around an optimum. Strict convexity alone does not supply a fixed positive $\mu$. For example, $x^4$ is strictly convex on the real line, but its second derivative vanishes at zero.

#### Quasiconvexity through sublevel sets

A function is quasiconvex when $f((1-t)x+ty)\le\max(f(x),f(y))$. This is equivalent to convexity of every sublevel set $\{x:f(x)\le a\}$. If the function inequality holds, a segment joining two points of one sublevel set remains in that set. Conversely, choose the sublevel height $a=\max(f(x),f(y))$ to recover the inequality.

Every convex function is quasiconvex, but $f(x)=x^3$ on the real line is a simple converse failure. It is increasing, so its sublevel sets are half-lines. Its second derivative is negative for negative $x$, so it is not convex on the whole line. Quasiconvexity controls the shape of sublevel sets and does not provide the full tangent-plane inequality.

#### Example AN-E21: a smooth maximum

For $x\in\mathbb R^d$ define $F(x)=\log\sum_i e^{x_i}$. Prove convexity and compare it with the maximum coordinate.

Solution. Put $p_i=e^{x_i}/\sum_j e^{x_j}$. Differentiation gives $\nabla F=p$ and $\nabla^2F=\operatorname{diag}(p)-pp^\top$. For any vector $h$,

$$
h^\top\nabla^2Fh
=\sum_i p_ih_i^2-\left(\sum_i p_ih_i\right)^2\ge0,
$$

because this is a variance under probabilities $p_i$. Hence $F$ is convex. The all-ones direction has zero curvature, so it is not strictly convex on all of $\mathbb R^d$.

With $m=\max_i x_i$, factor out $e^m$. The remaining sum lies between one and $d$, giving $m\le F(x)\le m+\log d$. More generally, $\gamma\log\sum_i e^{x_i/\gamma}$ approximates the maximum with error at most $\gamma\log d$ for $\gamma>0$. The same factorization avoids overflow in numerical evaluation.

### 5.3 Worked examples

#### Example AN-E09: convexity of an ellipsoid

Let $H$ be PSD and $C=\{x:x^\top Hx\le1\}$.

Solution. The function $f(x)=x^\top Hx$ has Hessian $2H$, so it is convex. For $x,y\in C$, convexity gives $f(\theta x+(1-\theta)y)\le\theta+(1-\theta)=1$, proving the set convex. If $H$ is singular, the set can be unbounded in its null directions.

#### Example AN-E10: a Hessian test with a parameter

For $f(x,y)=x^2+2axy+y^2$, the Hessian is $2\left(\begin{smallmatrix}1&a\\a&1\end{smallmatrix}\right)$, with eigenvalues $2(1+a)$ and $2(1-a)$. Thus $f$ is convex exactly when $|a|\le1$, and strictly convex when $|a|<1$. At $a=1$, it equals $(x+y)^2$ and is flat along $x=-y$.

### 5.4 Selected exercises

#### Problem AN-12: convex combinations of any finite length

Source: [Boyd and Vandenberghe, Exercise 2.1](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf).

Problem. Prove that a convex set contains every finite convex combination of its points.

Solution. For two points this is the definition of a convex set. Assume the claim for $k-1$ points. Given weights $\theta_1,\ldots,\theta_k$, if $\theta_k=1$ the combination is $x_k$. Otherwise the normalized combination $z=\sum_{i<k}\theta_ix_i/(1-\theta_k)$ belongs to $C$ by induction. The full combination $(1-\theta_k)z+\theta_kx_k$ belongs to $C$ by two-point convexity. This establishes the claim for every finite $k$.

#### Problem AN-13: testing convexity on lines

Source: [Boyd and Vandenberghe, Exercise 2.2](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf).

Problem. Prove that a set is convex exactly when its intersection with every line is convex. Establish the corresponding affine-set statement.

Solution. If $C$ is convex, intersecting it with a line preserves convexity because both sets contain the segment between any common points. Conversely, any two distinct $x,y\in C$ lie on a line; convexity of its intersection with $C$ supplies the entire segment. The same reasoning with all real coefficients $\theta$, rather than $\theta\in[0,1]$, proves the affine version. Empty intersections are convex and affine under the usual vacuous convention.

#### Problem AN-14: closed midpoint-convex sets

Source: [Boyd and Vandenberghe, Exercise 2.3](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf).

Problem. Prove that a closed set containing the midpoint of every pair of its points is convex. Explain why closedness cannot simply be dropped.

Solution. Repeated midpoint formation shows that $(1-k/2^m)x+(k/2^m)y\in C$ for every integer $0\le k\le2^m$. For any $\theta\in[0,1]$, choose dyadic rationals $k_m/2^m\to\theta$. The corresponding points of $C$ converge to $(1-\theta)x+\theta y$. Closedness places the limit in $C$, proving convexity. Without closedness, $\mathbb Q\subset\mathbb R$ is midpoint convex but misses many points on the segment from $0$ to $1$.

#### Problem AN-15: running averages

Source: [Boyd and Vandenberghe, Exercise 3.5](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf).

Problem. If $f$ is a finite convex function on a domain containing $[0,\infty)$, show that its running average $F(x)=x^{-1}\int_0^xf(t)dt$ is convex for $x>0$.

Solution. For convex $f$ on a domain containing $[0,\infty)$,

$$
F(x)=\frac1x\int_0^x f(t)\,dt=\int_0^1f(sx)\,ds,\qquad x>0.
$$

For each $s\in[0,1]$, convexity gives $f(s[\theta x+(1-\theta)y])\le\theta f(sx)+(1-\theta)f(sy)$. Integrating preserves the inequality and proves convexity of $F$. Convex finite functions are continuous in the interior; their possible endpoint behavior here is bounded by nearby finite secants, so these finite-interval integrals exist.

#### Problem AN-16: a convex function bounded above on all space

Source: [Boyd and Vandenberghe, Exercise 3.7](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf).

Problem. Prove that a finite convex function on all of $\mathbb R^n$ which is bounded above must be constant.

Solution. Suppose $f(y)>f(x)$ for two points. For $t>1$, let $z_t=x+t(y-x)$. Since $y=(1-1/t)x+(1/t)z_t$, convexity implies

$$
f(z_t)\ge f(x)+t[f(y)-f(x)]\longrightarrow\infty.
$$

This contradicts the assumed upper bound. Reversing $x,y$ rules out the other strict inequality. Thus every pair has equal function values and $f$ is constant. The full domain matters because the ray $z_t$ must remain available.

#### Problem AN-17: the Hessian criterion

Source: [Boyd and Vandenberghe, Exercise 3.8](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf). Work on an open convex domain with a twice continuously differentiable function.

Solution. If $f$ is convex, its restriction $\phi(t)=f(x+tv)$ is convex for every line through an interior point. The midpoint inequality gives $\phi(h)-2\phi(0)+\phi(-h)\ge0$. Divide by $h^2$ and let $h\to0$ to get $v^\top\nabla^2f(x)v\ge0$.

Conversely, if the Hessian is PSD throughout the domain, every line restriction has nonnegative second derivative. Its first derivative is increasing, so the function lies below its secants and is convex on each segment. This gives the multivariable convexity inequality. Openness ensures that the directional derivatives used in the first implication are defined on both sides.

## 6. Constrained optimization

### 6.1 Theory

Continuity on a nonempty compact feasible set ensures existence. On a closed unbounded set, coercivity, $f(x)\to\infty$ as $\|x\|\to\infty$, can replace boundedness: any minimizing sequence lies eventually in a compact sublevel set. Strict convexity on a convex feasible set ensures at most one minimizer.

For

$$
\min_x f_0(x)\quad\text{subject to }g_i(x)\le0,\quad h_j(x)=0,
$$

the Lagrangian is $L=f_0+\sum_i\lambda_i g_i+\sum_j\nu_jh_j$. The KKT conditions are primal feasibility, $\lambda_i\ge0$, stationarity $\nabla_xL=0$, and complementary slackness $\lambda_i g_i(x)=0$.

For a differentiable convex problem with affine equalities, any KKT point is globally optimal. Necessity requires a suitable constraint qualification, such as Slater's condition in the convex setting. A constraint with $g_i(x)=0$ is active; an inactive constraint has zero multiplier. An active constraint may also have zero multiplier.

### 6.2 Derivations and proofs

#### Where multipliers come from

For smooth equality constraints $h(x)=0$, assume their gradients are independent at a local optimum. The implicit-function theorem identifies the feasible tangent directions with $\ker Dh(x)$. The objective derivative must vanish along these directions, so $\nabla f_0$ lies in the orthogonal complement of that kernel, which is the range of $Dh(x)^\top$. Thus $\nabla f_0+Dh^\top\nu=0$ for some $\nu$.

With active inequalities and the linear-independence constraint qualification, the linearized feasible directions satisfy $\nabla g_i^\top d\le0$ and $Dh\,d=0$. No such direction can have $\nabla f_0^\top d<0$. The finite-dimensional theorem of alternatives expresses this as stationarity with nonnegative coefficients on active inequality gradients. Set inactive multipliers to zero; complementary slackness follows. This invokes the theorem of alternatives rather than proving general multiplier theory.

For convex sufficiency, convexity and stationarity make $x^*$ a global minimizer of $L(\cdot,\lambda^*,\nu^*)$. For any feasible $x$,

$$
f_0(x)\ge L(x,\lambda^*,\nu^*)
\ge L(x^*,\lambda^*,\nu^*)=f_0(x^*).
$$

The first inequality uses nonnegative multipliers; the last uses feasibility and complementary slackness.

#### Necessary conditions and sufficient certificates

For an unconstrained differentiable objective at an interior local minimum, both small positive and negative steps in every coordinate are feasible. Their directional derivatives force $\nabla f(x^*)=0$. At a boundary, some directions are forbidden, so setting every derivative to zero is generally incorrect.

For differentiable convex $f$ on a convex feasible set $C$, a useful optimality condition is

$$
\nabla f(x^*)^\top(y-x^*)\ge0\qquad(y\in C).
$$

It is necessary by differentiating the objective along each feasible segment from $x^*$. It is sufficient because the first-order convexity inequality then gives $f(y)\ge f(x^*)$. The condition says that every feasible first-order move has nonnegative cost.

For smooth inequality and equality descriptions, KKT multipliers express this geometric condition through the constraint gradients. Necessity requires a constraint qualification. For example, linear independence of active inequality gradients together with equality gradients is a common smooth qualification. In a convex problem with affine equalities, Slater's strict-feasibility condition supports existence of optimal multipliers under the usual finite-optimum assumptions.

Sufficiency has a direct proof that does not need a constraint qualification. Suppose $f,g_i$ are convex, equalities are affine, and a feasible $x^*$ satisfies stationarity, $\lambda_i\ge0$, and complementary slackness. The Lagrangian is convex in $x$ and has zero gradient at $x^*$, so $L(x,\lambda,\nu)\ge L(x^*,\lambda,\nu)$ for every $x$. For feasible $x$, $L(x,\lambda,\nu)\le f(x)$. At $x^*$, complementary slackness makes $L(x^*,\lambda,\nu)=f(x^*)$. Combining gives $f(x)\ge f(x^*)$.

For a nonconvex problem, stationarity and complementary slackness alone do not supply this global certificate. Second-order tests must use feasible tangent directions, and even a strict local minimum can fail to be global.

#### Example AN-E22: allocating a positive resource

Minimize $-\log x-\log y$ subject to $x+y\le B$, with $B>0$ and $x,y>0$.

Solution. A feasible point with a strict resource inequality exists, such as $x=y=B/4$. The objective is strictly convex because its Hessian is $\operatorname{diag}(1/x^2,1/y^2)$.

With multiplier $\lambda\ge0$ on $x+y-B$, stationarity gives $-1/x+\lambda=0$ and $-1/y+\lambda=0$. Thus $x=y=1/\lambda$ and the multiplier must be positive. Complementary slackness then makes the resource constraint bind. Therefore

$$
x^*=y^*=\frac B2,\qquad \lambda^*=\frac2B,
\qquad f(x^*,y^*)=-2\log(B/2).
$$

All KKT conditions hold, so convexity proves global optimality. Strict convexity proves uniqueness. The positive domain is part of the problem, not a pair of constraints on which a zero allocation could be allowed: the logarithm is undefined at zero and the objective diverges there. Increasing $B$ lowers minimum cost at rate $-2/B=-\lambda^*$.

### 6.3 Worked examples

#### Example AN-E11: a strictly convex unconstrained quadratic

Minimize $f(x,y)=x^2+xy+y^2-4x-2y$.

Solution. Stationarity gives $2x+y=4$ and $x+2y=2$, so $(x,y)=(2,0)$. The Hessian $\left(\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right)$ has eigenvalues $3,1$, proving strict convexity. Hence this stationary point is the unique global minimum, with value $-4$.

#### Example AN-E12: a boundary optimum

Minimize $(x-2)^2$ subject to $x\le1$.

Solution. The Lagrangian is $(x-2)^2+\lambda(x-1)$. The point $x=1$ satisfies stationarity with $\lambda=2$, and complementary slackness holds because $\lambda(x-1)=0$. The objective is strictly convex, so this is the unique minimum. The positive multiplier measures the first-order benefit of relaxing the bound.

### 6.4 Selected exercises

#### Problem AN-18: five objectives on one feasible set

Source: [Boyd and Vandenberghe, Exercise 4.1(a–e)](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf), analytical solutions.

The common constraints are $2x+y\ge1$, $x+3y\ge1$, and $x,y\ge0$.

Solution. The lower boundary joins $(0,1)$, $(2/5,1/5)$, and $(1,0)$; the feasible set extends upward and rightward.

(a) To minimize $x+y$, combine the first two inequalities with weights $2/5$ and $1/5$:

$$
x+y=\tfrac25(2x+y)+\tfrac15(x+3y)\ge\tfrac35.
$$

Equality forces both constraints to bind, so $(2/5,1/5)$ is the unique minimizer.

(b) The objective $-x-y$ tends to $-\infty$ along $(t,t)$ as $t\to\infty$, so no minimum exists.

(c) The minimum of $x$ is zero, attained at all $(0,y)$ with $y\ge1$.

(d) If $t=\max(x,y)$, then $2x+y\le3t$, so feasibility gives $t\ge1/3$. The point $(1/3,1/3)$ is feasible and attains the bound. Equality in $2x+y\le3t$ forces $x=y=t$, so it is unique.

(e) For $x^2+9y^2$, use $\nabla f=(2x,18y)$. If only $2x+y=1$ binds, stationarity gives $x=\lambda$, $y=\lambda/18$, and hence $\lambda=18/37$. Thus $(x,y)=(18/37,1/37)$ would violate $x+3y\ge1$. At the intersection $(2/5,1/5)$, solve

$$
\binom{4/5}{18/5}
=\lambda_1\binom21+\lambda_2\binom13.
$$

This gives $\lambda_1=-6/25$, so that point cannot satisfy KKT with nonnegative multipliers. Instead bind $x+3y=1$: stationarity gives $2x=\lambda$, $18y=3\lambda$, hence $x=1/2$, $y=1/6$. The first constraint is slack because $2x+y=7/6>1$. The multiplier is $\lambda=1$, all other multipliers vanish, and KKT holds. Strict convexity gives the unique minimum $1/4+9/36=1/2$.

## 7. Lagrangian duality

### 7.1 Theory

The dual function is

$$
q(\lambda,\nu)=\inf_xL(x,\lambda,\nu),\qquad \lambda\ge0.
$$

It is concave in the multipliers because it is an infimum of affine functions. Maximizing it gives the dual problem. Weak duality, $d^*\le p^*$, follows for every primal-feasible $x$ and dual-feasible $(\lambda,\nu)$ from $q(\lambda,\nu)\le L(x,\lambda,\nu)\le f_0(x)$.

Strong duality means equality of optimal values. A standard sufficient condition is a convex problem with affine equalities and a point in the relative interior of the objective domain satisfying every nonlinear inequality strictly. This is Slater's condition. In the usual finite-dimensional setting, a finite optimum then also has dual multipliers. Strong duality alone does not imply primal attainment.

For a perturbation $g_i(x)\le u_i$, the dual expression subtracts $\lambda^\top u$. When the value function is differentiable and the usual sensitivity assumptions hold, $\partial p^*/\partial u_i=-\lambda_i^*$. More generally the multiplier gives a supporting affine bound.

### 7.2 Derivations and proofs

#### A quadratic-program dual

For $H$ PD, consider $\min \tfrac12x^\top Hx+c^\top x$ subject to $Ax=b$. Completing the square in $x$ gives

$$
x(\nu)=-H^{-1}(c+A^\top\nu),
$$

$$
q(\nu)=-\tfrac12(c+A^\top\nu)^\top H^{-1}(c+A^\top\nu)-b^\top\nu.
$$

This is a concave quadratic. If $A$ has independent rows, its stationarity equation determines a unique multiplier. The resulting $x$ satisfies the equality constraint, and the primal and dual values agree by substitution. This establishes strong duality directly for this problem.

#### Reading a dual bound in the correct direction

For a minimization problem with $g_i(x)\le0$, a nonnegative multiplier makes $\lambda_i g_i(x)\le0$ at every feasible point. Therefore $L(x,\lambda,\nu)\le f(x)$ at feasible $x$. Minimizing the Lagrangian over all $x$ can only decrease it further. Every dual-feasible multiplier consequently gives a lower bound on every feasible primal cost:

$$
g(\lambda,\nu)=\inf_x L(x,\lambda,\nu)\le p^*.
$$

The dual maximizes these lower bounds. A positive duality gap means that the best available bound is still below the minimum cost. It does not mean that a dual point violates a primal constraint, since the two points belong to different spaces.

The dual function is concave even when the original problem is nonconvex: it is the pointwise infimum of affine functions of the multipliers. For a convex combination of multiplier vectors, each fixed-$x$ Lagrangian value is the corresponding affine combination. Taking the infimum gives the concavity inequality.

#### An equality-constrained quadratic with matching values

Consider $\min \tfrac12x^2+2y^2$ subject to $x+y=b$. Use the Lagrangian

$$
L(x,y,\nu)=\tfrac12x^2+2y^2+\nu(x+y-b).
$$

For any fixed real $\nu$, its unique unrestricted minimum has $x=-\nu$ and $y=-\nu/4$. Substituting gives

$$
g(\nu)=-\frac58\nu^2-b\nu.
$$

This is a concave quadratic. Its maximum occurs at $\nu^*=-4b/5$, producing $x^*=4b/5$, $y^*=b/5$. Both primal and dual values equal $2b^2/5$. For $b=1$, the solution is $(4/5,1/5)$ and the common value is $2/5$.

The equality multiplier is unrestricted in sign. Imposing $\nu\ge0$ would be an extra condition that destroys the correct dual solution for positive $b$.

#### Sensitivity from an optimal multiplier

Let $p(u)$ be the minimum cost when the inequalities are relaxed to $g_i(x)\le u_i$, while the other constraints remain fixed. Suppose strong duality holds at $u=0$ and an optimal multiplier $\lambda^*$ exists. Using it at the perturbed problem gives

$$
p(u)\ge p(0)-(\lambda^*)^\top u.
$$

This is a global supporting bound. If the value function is differentiable at zero, its gradient there is $-\lambda^*$. Differentiability is an additional property, not a consequence of having written down a multiplier.

For the equality example, $p(b)=2b^2/5$, so $p'(b)=4b/5=-\nu^*(b)$. The minus sign comes from the term $-\nu b$ in the Lagrangian. Reversing the equality convention would reverse the multiplier sign while leaving the derivative of the actual value unchanged.

Strong duality concerns equality of optimal values. It should be distinguished from attainment of the dual maximum. For $\min x$ subject to $x^2\le0$, the only feasible point is zero and the primal value is zero. For $\lambda>0$, the dual value is $-1/(4\lambda)$; its supremum is zero but no finite multiplier attains it. Strict feasibility fails. The limiting lower bounds become exact without yielding a finite KKT multiplier.

### 7.3 Worked examples

#### Example AN-E13: a dual certificate

Minimize $x^2$ subject to $x\ge1$.

Solution. With $g(x)=1-x$, $L=x^2+\lambda(1-x)$. Its minimizer is $x=\lambda/2$, so $q(\lambda)=\lambda-\lambda^2/4$. The maximum over $\lambda\ge0$ occurs at $\lambda=2$ and equals one. The feasible point $x=1$ also has value one, certifying both optima.

#### Example AN-E14: failure of a constraint qualification

Minimize $x$ subject to $x^2\le0$.

Solution. The only feasible point is zero and the primal value is zero. At this point, stationarity would require $1+2\lambda x=1=0$, impossible for any finite multiplier. For $\lambda>0$, $q(\lambda)=-1/(4\lambda)$; its supremum is zero but is never attained. Slater fails. The example separates equality of values from the existence of KKT multipliers.

### 7.4 Selected exercises

#### Problem AN-19: a complete scalar dual and sensitivity calculation

Source: [Boyd and Vandenberghe, Exercise 5.1](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf), parts (a), (c), (d), and the algebraic part of (b).

Solution. Minimize $x^2+1$ subject to $(x-2)(x-4)\le0$. The feasible interval is $[2,4]$, so $x^*=2$ and $p^*=5$. The Lagrangian is

$$
L(x,\lambda)=(1+\lambda)x^2-6\lambda x+8\lambda+1.
$$

For $\lambda\ge0$, its minimizer is $x=3\lambda/(1+\lambda)$, giving

$$
q(\lambda)=1+8\lambda-\frac{9\lambda^2}{1+\lambda}.
$$

Its derivative is $q'(\lambda)=-1+9/(1+\lambda)^2$ and its second derivative is $-18/(1+\lambda)^3<0$. Hence $\lambda^*=2$ uniquely maximizes it and $q(2)=5$, establishing strong duality. For every feasible $x$, $q(\lambda)\le L(x,\lambda)\le x^2+1$, verifying the lower-bound property.

For $(x-2)(x-4)\le u$, rewrite $(x-3)^2\le1+u$. It is infeasible for $u<-1$. For $-1\le u<8$, the feasible interval stays to the right of zero and its closest point to zero is $3-\sqrt{1+u}$. For $u\ge8$, zero is feasible. Thus

$$
p^*(u)=
\begin{cases}
+\infty,&u<-1,\\
1+(3-\sqrt{1+u})^2,&-1\le u<8,\\
1,&u\ge8.
\end{cases}
$$

At zero, the derivative is $-(3-\sqrt1)/\sqrt1=-2=-\lambda^*$.

#### Problem AN-20: weak duality with infinite values

Source: [Boyd and Vandenberghe, Exercise 5.2](https://www.seas.ucla.edu/~vandenbe/cvxbook/bv_cvxbook.pdf).

Problem. Use weak duality to determine what follows when a primal optimum is $-\infty$ or a dual optimum is $+\infty$, under the usual finite-valued objective convention.

Solution. If $p^*=-\infty$, there are feasible $x_k$ whose objective tends to $-\infty$. Every dual-feasible multiplier gives $q(\lambda,\nu)\le f_0(x_k)$ for all $k$, forcing $q=-\infty$. Therefore $d^*=-\infty$, including the convention for an empty dual-feasible set.

If $d^*=+\infty$ and there were a primal-feasible point with finite objective, that value would bound every dual value above, a contradiction. In the usual finite-valued objective formulation there is consequently no primal-feasible point, and $p^*=+\infty$.

## 8. Optimization with uncertainty

### 8.1 Theory

A robust constraint requires feasibility for every uncertainty realization:

$$
g(x,u)\le0\quad\text{for every }u\in\mathcal U.
$$

Its feasible set is an intersection over $u$. If each $g(\cdot,u)$ is convex, that intersection is convex. A worst-case objective $\sup_{u\in\mathcal U}f(x,u)$ is convex in $x$ whenever every $f(\cdot,u)$ is convex.

To maximize the worst of finitely many rewards $r_i(x)$, introduce $t$ and solve $\max t$ subject to $t\le r_i(x)$ for every $i$. With concave rewards and a convex decision set this is a convex optimization problem in maximization form.

Existence still needs checking. A supremum over an open uncertainty set need not be attained by any realization. A minimax interchange also needs assumptions; changing $\min_x\max_u$ to $\max_u\min_x$ without justification changes the order in which decisions are made.

### 8.2 Derivations and proofs

#### A norm-bounded uncertain linear constraint

Suppose $(a_0+u)^\top x\le b$ must hold for every $\|u\|_2\le\rho$. Cauchy-Schwarz gives $u^\top x\le\rho\|x\|_2$. For $x\ne0$, equality is attained at $u=\rho x/\|x\|_2$; for $x=0$, it is zero for all $u$. Therefore the exact robust counterpart is

$$
a_0^\top x+\rho\|x\|_2\le b.
$$

The extra term is the support function of the uncertainty ball. It is a convex function of $x$.

#### Robust feasible sets as intersections

Suppose the decision $x$ is chosen before an uncertain coefficient $a$ is known. Requiring $a^\top x\le b$ for every $a\in\mathcal U$ gives the feasible set

$$
C_{\mathrm{rob}}=\bigcap_{a\in\mathcal U}\{x:a^\top x\le b\}.
$$

Each set in the intersection is a half-space. Their intersection is convex, regardless of whether the uncertainty set itself is convex. Computing a compact reformulation is useful, but convexity already follows from this intersection description.

If uncertainty sets are nested, $\mathcal U_1\subseteq\mathcal U_2$, the corresponding robust feasible sets are nested in the reverse direction. Thus enlarging uncertainty cannot lower the minimum of a fixed objective, and cannot raise the maximum. The optimizer need not move smoothly, because a new worst-case coefficient or a new active constraint can take over.

For ellipsoidal uncertainty $a=a_0+Lu$ with $\|u\|_2\le\rho$,

$$
\sup_a a^\top x
=a_0^\top x+\sup_{\|u\|\le\rho}u^\top L^\top x
=a_0^\top x+\rho\|L^\top x\|_2.
$$

Cauchy–Schwarz proves the upper bound. If $L^\top x\ne0$, equality is attained at $u=\rho L^\top x/\|L^\top x\|$; if it is zero, every allowed $u$ has the same effect. The formula remains valid for singular $L$.

#### The order of optimization

For any real payoff $F(x,u)$, one always has the weak minimax inequality

$$
\sup_u\inf_x F(x,u)\le\inf_x\sup_u F(x,u).
$$

Indeed, for each fixed pair $x,u$, the infimum over decisions is no larger than $F(x,u)$, which is no larger than the supremum over uncertainty. Taking the remaining bounds gives the inequality. Equality requires additional structure.

The pure two-action payoff $xu$ with $x,u\in\{-1,1\}$ gives a strict gap. If both players may randomize independently, with probabilities $p,q$ of choosing one, the expected payoff becomes $(2p-1)(2q-1)$. Choosing $p=1/2$ guarantees zero for the minimizing player, and choosing $q=1/2$ guarantees zero for the maximizing player. The mixed game therefore has value zero. Randomization changed the feasible objects and removed this particular gap.

#### Example AN-E23: the cost of coefficient uncertainty

Maximize $x_1+x_2$ with $x_1,x_2\ge0$ subject to $a^\top x\le1$ for every $a$ in a Euclidean ball of radius $\rho\ge0$ centered at $(1,1)$.

Solution. The robust constraint is $x_1+x_2+\rho\sqrt{x_1^2+x_2^2}\le1$. Put $t=x_1+x_2$. Cauchy–Schwarz gives $\sqrt{x_1^2+x_2^2}\ge t/\sqrt2$, with equality at equal coordinates. Every feasible decision therefore satisfies

$$
t\left(1+\frac\rho{\sqrt2}\right)\le1.
$$

The upper bound is attained by setting

$$
t^*=\frac1{1+\rho/\sqrt2},\qquad
x_1^*=x_2^*=\frac1{2(1+\rho/\sqrt2)}.
$$

For positive uncertainty radius, equal allocation is required to attain the bound. At $\rho=0$, every nonnegative allocation with sum one is optimal. The optimal value decreases continuously with the radius, while uniqueness changes at zero. Both facts follow from the norm inequality and its equality condition.

### 8.3 Worked examples

#### Example AN-E15: a robust upper bound

Maximize $x\ge0$ subject to $ax\le2$ for all $a\in[1,3]$.

Solution. For nonnegative $x$, the worst coefficient is $a=3$, so the constraint is $3x\le2$. The maximum is $x^*=2/3$. If negative $x$ were allowed, selecting the worst coefficient would depend on the sign; maximizing over the uncertainty set must be done before discarding that dependence.

#### Example AN-E16: balancing two rewards

Maximize $\min(x,2-2x)$ over $0\le x\le1$.

Solution. The epigraph-style maximization form is $\max t$ with $t\le x$ and $t\le2-2x$. The first bound increases and the second decreases; their intersection is $x=2/3$, $t=2/3$. For $x$ below this point, the first reward is smaller; for $x$ above it, the second is smaller. This proves global optimality.

### 8.4 Additional solved problems

#### Problem AN-21: a robust minimum-distance decision

Minimize $(x-1)^2$ subject to $ax\le1$ for all $a\in[1,2]$.

Solution. If $x\ge0$, the strongest inequality is $2x\le1$. If $x<0$, all constraints hold automatically. Thus the full feasible set is $(-\infty,1/2]$. The unconstrained minimum at one is outside this set, so the closest feasible point is $x^*=1/2$ and the minimum value is $1/4$. Reducing the uncertainty interval to $[1,1]$ would allow the unconstrained optimum.

#### Problem AN-22: minimum worst-case squared error

Choose $x\in\mathbb R$ to minimize $\max\{(x-1)^2,(x+3)^2\}$.

Solution. The two reference points are distance four apart. Any $x$ must be at least distance two from one of them; otherwise the triangle inequality would make their distance less than four. At their midpoint $x=-1$, both distances equal two. The minimum worst-case squared error is therefore four. Equivalently solve $\min t$ subject to $(x-1)^2\le t$ and $(x+3)^2\le t$.

#### Problem AN-23: a failed minimax interchange

Let $x,u\in\{-1,1\}$ and payoff $xu$. Compare the two orders.

Solution. For fixed $x$, the maximizing player chooses $u=x$, giving payoff one. Hence $\min_x\max_u xu=1$. For fixed $u$, the minimizing player chooses $x=-u$, giving payoff minus one. Hence $\max_u\min_x xu=-1$. The nonconvex decision sets do not meet the usual convex-concave minimax assumptions. Allowing mixed decisions gives a different problem.

#### Problem AN-24: a value function with a corner

Minimize $x^2$ subject to $x\ge b$, with $b$ a real parameter.

Solution. If $b\le0$, zero is feasible and optimal. If $b>0$, the solution is $b$. Thus $p(b)=\max(b,0)^2$, which is differentiable at zero but has a change in second derivative. The multiplier for $b-x\le0$ is $\lambda(b)=2\max(b,0)$ and $p'(b)=\lambda(b)$. The sign differs from the relaxation convention $g(x)\le u$ because increasing $b$ tightens this constraint.

<!-- extended-practice-start -->

## 9. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 9.1 Finding an explicit convergence index

An epsilon proof must produce a single index that works for every later term. Computing a limit by dividing numerator and denominator by $n$ suggests the answer, but does not yet provide that index. The following calculations also show how the sign of an error determines monotonicity.

#### Problem AN-25: rational sequence and epsilon index, case 1

For $x_n=(n+2)/(n+5)$, $n\ge1$, prove convergence, determine monotonicity, and give an index for error less than $\frac{1}{100}$. Give a direct Cauchy bound as well.

Solution. Subtract the proposed limit before taking absolute values:

$$
x_n-1=\frac{-3}{n+5},\qquad
|x_n-1|=\frac{3}{n+5}.
$$

Since the denominator grows without bound and is positive for every allowed index, the limit is one. To obtain a quantitative proof, the inequality $|x_n-1|<\epsilon$ is equivalent to $n>3/\epsilon-5$. At the requested tolerance, one valid index is $N=296$. For every integer $n\ge N$ the strict inequality holds. Using a floor followed by one avoids an equality at the boundary.

For arbitrary $\epsilon>0$ the same construction, with a maximum against one, proves the full definition of convergence. The chosen index may depend on $\epsilon$; it cannot depend on which later term is being checked.

Consecutive differences are

$$
x_{n+1}-x_n
=\frac{3}{(n+5)(n+1+5)}.
$$

Thus the sequence is strictly increasing. Its terms remain on the lower side of one. This is consistent with convergence and gives a second route using monotone convergence.

For a direct Cauchy estimate, let $m,n\ge M$. The triangle inequality gives

$$
|x_n-x_m|
\le |x_n-1|+|x_m-1|
\le\frac{6}{M+5}.
$$

Taking $M=596$ makes this less than the requested tolerance. This bound is sufficient even though it is not the smallest possible one. An epsilon proof needs a valid uniform bound; it does not need the optimal index.

#### Problem AN-26: rational sequence and epsilon index, case 2

For $x_n=(n+7)/(n+1)$, $n\ge1$, prove convergence, determine monotonicity, and give an index for error less than $\frac{1}{50}$. Give a direct Cauchy bound as well.

Solution. Subtract the proposed limit before taking absolute values:

$$
x_n-1=\frac{6}{n+1},\qquad
|x_n-1|=\frac{6}{n+1}.
$$

Since the denominator grows without bound and is positive for every allowed index, the limit is one. To obtain a quantitative proof, the inequality $|x_n-1|<\epsilon$ is equivalent to $n>6/\epsilon-1$. At the requested tolerance, one valid index is $N=300$. For every integer $n\ge N$ the strict inequality holds. Using a floor followed by one avoids an equality at the boundary.

For arbitrary $\epsilon>0$ the same construction, with a maximum against one, proves the full definition of convergence. The chosen index may depend on $\epsilon$; it cannot depend on which later term is being checked.

Consecutive differences are

$$
x_{n+1}-x_n
=\frac{-6}{(n+1)(n+1+1)}.
$$

Thus the sequence is strictly decreasing. Its terms remain on the upper side of one. This is consistent with convergence and gives a second route using monotone convergence.

For a direct Cauchy estimate, let $m,n\ge M$. The triangle inequality gives

$$
|x_n-x_m|
\le |x_n-1|+|x_m-1|
\le\frac{12}{M+1}.
$$

Taking $M=600$ makes this less than the requested tolerance. This bound is sufficient even though it is not the smallest possible one. An epsilon proof needs a valid uniform bound; it does not need the optimal index.

#### Problem AN-27: rational sequence and epsilon index, case 3

For $x_n=(n+-3)/(n+4)$, $n\ge1$, prove convergence, determine monotonicity, and give an index for error less than $\frac{1}{1000}$. Give a direct Cauchy bound as well.

Solution. Subtract the proposed limit before taking absolute values:

$$
x_n-1=\frac{-7}{n+4},\qquad
|x_n-1|=\frac{7}{n+4}.
$$

Since the denominator grows without bound and is positive for every allowed index, the limit is one. To obtain a quantitative proof, the inequality $|x_n-1|<\epsilon$ is equivalent to $n>7/\epsilon-4$. At the requested tolerance, one valid index is $N=6997$. For every integer $n\ge N$ the strict inequality holds. Using a floor followed by one avoids an equality at the boundary.

For arbitrary $\epsilon>0$ the same construction, with a maximum against one, proves the full definition of convergence. The chosen index may depend on $\epsilon$; it cannot depend on which later term is being checked.

Consecutive differences are

$$
x_{n+1}-x_n
=\frac{7}{(n+4)(n+1+4)}.
$$

Thus the sequence is strictly increasing. Its terms remain on the lower side of one. This is consistent with convergence and gives a second route using monotone convergence.

For a direct Cauchy estimate, let $m,n\ge M$. The triangle inequality gives

$$
|x_n-x_m|
\le |x_n-1|+|x_m-1|
\le\frac{14}{M+4}.
$$

Taking $M=13997$ makes this less than the requested tolerance. This bound is sufficient even though it is not the smallest possible one. An epsilon proof needs a valid uniform bound; it does not need the optimal index.

#### Problem AN-28: rational sequence and epsilon index, case 4

For $x_n=(n+1)/(n+3)$, $n\ge1$, prove convergence, determine monotonicity, and give an index for error less than $\frac{1}{20}$. Give a direct Cauchy bound as well.

Solution. Subtract the proposed limit before taking absolute values:

$$
x_n-1=\frac{-2}{n+3},\qquad
|x_n-1|=\frac{2}{n+3}.
$$

Since the denominator grows without bound and is positive for every allowed index, the limit is one. To obtain a quantitative proof, the inequality $|x_n-1|<\epsilon$ is equivalent to $n>2/\epsilon-3$. At the requested tolerance, one valid index is $N=38$. For every integer $n\ge N$ the strict inequality holds. Using a floor followed by one avoids an equality at the boundary.

For arbitrary $\epsilon>0$ the same construction, with a maximum against one, proves the full definition of convergence. The chosen index may depend on $\epsilon$; it cannot depend on which later term is being checked.

Consecutive differences are

$$
x_{n+1}-x_n
=\frac{2}{(n+3)(n+1+3)}.
$$

Thus the sequence is strictly increasing. Its terms remain on the lower side of one. This is consistent with convergence and gives a second route using monotone convergence.

For a direct Cauchy estimate, let $m,n\ge M$. The triangle inequality gives

$$
|x_n-x_m|
\le |x_n-1|+|x_m-1|
\le\frac{4}{M+3}.
$$

Taking $M=78$ makes this less than the requested tolerance. This bound is sufficient even though it is not the smallest possible one. An epsilon proof needs a valid uniform bound; it does not need the optimal index.

### 9.2 A usable delta for continuity

In a continuity proof, an expression involving both $x$ and the base point must be bounded using information fixed in advance. The common first restriction $|x-a|<1$ creates such a local bound. A second restriction then controls the desired error. Both restrictions are enforced by taking a minimum.

#### Problem AN-29: epsilon-delta proof for a square, case 1

Prove continuity of $f(x)=x^2$ at $a=2$ and exhibit a delta for $\epsilon=\frac{1}{10}$.

Solution. Factor the difference:

$$
|x^2-a^2|=|x-a||x+a|.
$$

The second factor still depends on $x$, so choosing $\delta=\epsilon/|x+a|$ would be circular. Instead first require $|x-a|<1$. Since $x+a=(x-a)+2a$, this yields $|x+a|\le|x-a|+2|a|<1+2|a|=5$.

Now choose

$$
\delta=\min\left(1,\frac{\epsilon}{5}\right).
$$

Whenever $|x-a|<\delta$, both required bounds hold, and therefore

$$
|f(x)-f(a)|< 5\delta\le\epsilon.
$$

For the supplied tolerance this choice is $\delta=\frac{1}{50}$. The base-point value is $f(a)=4$. Every point within the chosen delta has its squared value within epsilon of that number.

The order of the choices is part of the proof. The base point is fixed first, epsilon is supplied next, and delta is then chosen using only those two quantities. The point $x$ is allowed to vary afterwards.

This argument also gives a local Lipschitz bound $|x^2-a^2|\le 5|x-a|$ within distance one of this base point. It does not produce a global Lipschitz constant on the whole real line. Indeed, comparing $x$ and $x+1$ gives a difference $|2x+1|$, which has no finite global upper bound. Continuity at every point and one global Lipschitz estimate are different assertions.

#### Problem AN-30: epsilon-delta proof for a square, case 2

Prove continuity of $f(x)=x^2$ at $a=-3$ and exhibit a delta for $\epsilon=\frac{1}{100}$.

Solution. Factor the difference:

$$
|x^2-a^2|=|x-a||x+a|.
$$

The second factor still depends on $x$, so choosing $\delta=\epsilon/|x+a|$ would be circular. Instead first require $|x-a|<1$. Since $x+a=(x-a)+2a$, this yields $|x+a|\le|x-a|+2|a|<1+2|a|=7$.

Now choose

$$
\delta=\min\left(1,\frac{\epsilon}{7}\right).
$$

Whenever $|x-a|<\delta$, both required bounds hold, and therefore

$$
|f(x)-f(a)|< 7\delta\le\epsilon.
$$

For the supplied tolerance this choice is $\delta=\frac{1}{700}$. The base-point value is $f(a)=9$. Every point within the chosen delta has its squared value within epsilon of that number.

The order of the choices is part of the proof. The base point is fixed first, epsilon is supplied next, and delta is then chosen using only those two quantities. The point $x$ is allowed to vary afterwards.

This argument also gives a local Lipschitz bound $|x^2-a^2|\le 7|x-a|$ within distance one of this base point. It does not produce a global Lipschitz constant on the whole real line. Indeed, comparing $x$ and $x+1$ gives a difference $|2x+1|$, which has no finite global upper bound. Continuity at every point and one global Lipschitz estimate are different assertions.

#### Problem AN-31: epsilon-delta proof for a square, case 3

Prove continuity of $f(x)=x^2$ at $a=0$ and exhibit a delta for $\epsilon=\frac{1}{1000}$.

Solution. Factor the difference:

$$
|x^2-a^2|=|x-a||x+a|.
$$

The second factor still depends on $x$, so choosing $\delta=\epsilon/|x+a|$ would be circular. Instead first require $|x-a|<1$. Since $x+a=(x-a)+2a$, this yields $|x+a|\le|x-a|+2|a|<1+2|a|=1$.

Now choose

$$
\delta=\min\left(1,\frac{\epsilon}{1}\right).
$$

Whenever $|x-a|<\delta$, both required bounds hold, and therefore

$$
|f(x)-f(a)|< 1\delta\le\epsilon.
$$

For the supplied tolerance this choice is $\delta=\frac{1}{1000}$. The base-point value is $f(a)=0$. Every point within the chosen delta has its squared value within epsilon of that number.

The order of the choices is part of the proof. The base point is fixed first, epsilon is supplied next, and delta is then chosen using only those two quantities. The point $x$ is allowed to vary afterwards.

This argument also gives a local Lipschitz bound $|x^2-a^2|\le 1|x-a|$ within distance one of this base point. It does not produce a global Lipschitz constant on the whole real line. Indeed, comparing $x$ and $x+1$ gives a difference $|2x+1|$, which has no finite global upper bound. Continuity at every point and one global Lipschitz estimate are different assertions.

#### Problem AN-32: epsilon-delta proof for a square, case 4

Prove continuity of $f(x)=x^2$ at $a=5$ and exhibit a delta for $\epsilon=\frac{1}{20}$.

Solution. Factor the difference:

$$
|x^2-a^2|=|x-a||x+a|.
$$

The second factor still depends on $x$, so choosing $\delta=\epsilon/|x+a|$ would be circular. Instead first require $|x-a|<1$. Since $x+a=(x-a)+2a$, this yields $|x+a|\le|x-a|+2|a|<1+2|a|=11$.

Now choose

$$
\delta=\min\left(1,\frac{\epsilon}{11}\right).
$$

Whenever $|x-a|<\delta$, both required bounds hold, and therefore

$$
|f(x)-f(a)|< 11\delta\le\epsilon.
$$

For the supplied tolerance this choice is $\delta=\frac{1}{220}$. The base-point value is $f(a)=25$. Every point within the chosen delta has its squared value within epsilon of that number.

The order of the choices is part of the proof. The base point is fixed first, epsilon is supplied next, and delta is then chosen using only those two quantities. The point $x$ is allowed to vary afterwards.

This argument also gives a local Lipschitz bound $|x^2-a^2|\le 11|x-a|$ within distance one of this base point. It does not produce a global Lipschitz constant on the whole real line. Indeed, comparing $x$ and $x+1$ gives a difference $|2x+1|$, which has no finite global upper bound. Continuity at every point and one global Lipschitz estimate are different assertions.

### 9.3 Taylor approximation with a remainder

For $\log(1+x)$, the first two Taylor terms are easy to write. The useful question is the size of what was omitted. A derivative bound must hold throughout the interval between zero and the evaluation point, including when that point is negative.

#### Problem AN-33: logarithm with a certified error, case 1

Approximate $\log(1+x)$ at $x=\frac{1}{10}$ by its degree-two Taylor polynomial at zero. Derive a rigorous error bound and compare it with the actual numerical error.

Solution. The derivatives are $f'(x)=1/(1+x)$, $f''(x)=-1/(1+x)^2$, and $f'''(x)=2/(1+x)^3$. At zero, the constant term is zero, the first derivative is one, and the second derivative is minus one. Thus

$$
P_2(x)=x-\frac{x^2}2,\qquad P_2(\frac{1}{10})=\frac{19}{200}.
$$

Every point $\xi$ between zero and the supplied $x$ lies in $[-q,q]$ for $q=\frac{1}{10}<1$. On this whole interval, $1+\xi\ge1-q>0$, so

$$
|f'''(\xi)|\le \frac2{(1-q)^3}.
$$

Taylor's theorem with Lagrange remainder then gives

$$
|\log(1+x)-P_2(x)|
\le\frac{|x|^3}{3(1-q)^3}
=\frac{1}{2187}.
$$

Numerically the approximation is $0.09500000$, the logarithm is $0.09531018$, and the actual absolute error is about $0.00031018$. It is below the bound $0.00045725$.

The bound is deliberately uniform over the interval. It can be conservative because the largest derivative on the interval need not occur near the unknown remainder point. Nevertheless it certifies the approximation without knowing that point.

For $x>0$, the third derivative is positive and the remainder has the sign of $x^3$, so the quadratic approximation is below the logarithm. For $-1<x<0$ the sign reverses. Near $x=-1$ the derivative bound grows sharply; a Taylor approximation about zero becomes less useful there. The domain restriction $x>-1$ must be checked before any expansion is used.

#### Problem AN-34: logarithm with a certified error, case 2

Approximate $\log(1+x)$ at $x=-\frac{1}{5}$ by its degree-two Taylor polynomial at zero. Derive a rigorous error bound and compare it with the actual numerical error.

Solution. The derivatives are $f'(x)=1/(1+x)$, $f''(x)=-1/(1+x)^2$, and $f'''(x)=2/(1+x)^3$. At zero, the constant term is zero, the first derivative is one, and the second derivative is minus one. Thus

$$
P_2(x)=x-\frac{x^2}2,\qquad P_2(-\frac{1}{5})=-\frac{11}{50}.
$$

Every point $\xi$ between zero and the supplied $x$ lies in $[-q,q]$ for $q=\frac{1}{5}<1$. On this whole interval, $1+\xi\ge1-q>0$, so

$$
|f'''(\xi)|\le \frac2{(1-q)^3}.
$$

Taylor's theorem with Lagrange remainder then gives

$$
|\log(1+x)-P_2(x)|
\le\frac{|x|^3}{3(1-q)^3}
=\frac{1}{192}.
$$

Numerically the approximation is $-0.22000000$, the logarithm is $-0.22314355$, and the actual absolute error is about $0.00314355$. It is below the bound $0.00520833$.

The bound is deliberately uniform over the interval. It can be conservative because the largest derivative on the interval need not occur near the unknown remainder point. Nevertheless it certifies the approximation without knowing that point.

For $x>0$, the third derivative is positive and the remainder has the sign of $x^3$, so the quadratic approximation is below the logarithm. For $-1<x<0$ the sign reverses. Near $x=-1$ the derivative bound grows sharply; a Taylor approximation about zero becomes less useful there. The domain restriction $x>-1$ must be checked before any expansion is used.

#### Problem AN-35: logarithm with a certified error, case 3

Approximate $\log(1+x)$ at $x=\frac{1}{4}$ by its degree-two Taylor polynomial at zero. Derive a rigorous error bound and compare it with the actual numerical error.

Solution. The derivatives are $f'(x)=1/(1+x)$, $f''(x)=-1/(1+x)^2$, and $f'''(x)=2/(1+x)^3$. At zero, the constant term is zero, the first derivative is one, and the second derivative is minus one. Thus

$$
P_2(x)=x-\frac{x^2}2,\qquad P_2(\frac{1}{4})=\frac{7}{32}.
$$

Every point $\xi$ between zero and the supplied $x$ lies in $[-q,q]$ for $q=\frac{1}{4}<1$. On this whole interval, $1+\xi\ge1-q>0$, so

$$
|f'''(\xi)|\le \frac2{(1-q)^3}.
$$

Taylor's theorem with Lagrange remainder then gives

$$
|\log(1+x)-P_2(x)|
\le\frac{|x|^3}{3(1-q)^3}
=\frac{1}{81}.
$$

Numerically the approximation is $0.21875000$, the logarithm is $0.22314355$, and the actual absolute error is about $0.00439355$. It is below the bound $0.01234568$.

The bound is deliberately uniform over the interval. It can be conservative because the largest derivative on the interval need not occur near the unknown remainder point. Nevertheless it certifies the approximation without knowing that point.

For $x>0$, the third derivative is positive and the remainder has the sign of $x^3$, so the quadratic approximation is below the logarithm. For $-1<x<0$ the sign reverses. Near $x=-1$ the derivative bound grows sharply; a Taylor approximation about zero becomes less useful there. The domain restriction $x>-1$ must be checked before any expansion is used.

#### Problem AN-36: logarithm with a certified error, case 4

Approximate $\log(1+x)$ at $x=-\frac{1}{3}$ by its degree-two Taylor polynomial at zero. Derive a rigorous error bound and compare it with the actual numerical error.

Solution. The derivatives are $f'(x)=1/(1+x)$, $f''(x)=-1/(1+x)^2$, and $f'''(x)=2/(1+x)^3$. At zero, the constant term is zero, the first derivative is one, and the second derivative is minus one. Thus

$$
P_2(x)=x-\frac{x^2}2,\qquad P_2(-\frac{1}{3})=-\frac{7}{18}.
$$

Every point $\xi$ between zero and the supplied $x$ lies in $[-q,q]$ for $q=\frac{1}{3}<1$. On this whole interval, $1+\xi\ge1-q>0$, so

$$
|f'''(\xi)|\le \frac2{(1-q)^3}.
$$

Taylor's theorem with Lagrange remainder then gives

$$
|\log(1+x)-P_2(x)|
\le\frac{|x|^3}{3(1-q)^3}
=\frac{1}{24}.
$$

Numerically the approximation is $-0.38888889$, the logarithm is $-0.40546511$, and the actual absolute error is about $0.01657622$. It is below the bound $0.04166667$.

The bound is deliberately uniform over the interval. It can be conservative because the largest derivative on the interval need not occur near the unknown remainder point. Nevertheless it certifies the approximation without knowing that point.

For $x>0$, the third derivative is positive and the remainder has the sign of $x^3$, so the quadratic approximation is below the logarithm. For $-1<x<0$ the sign reverses. Near $x=-1$ the derivative bound grows sharply; a Taylor approximation about zero becomes less useful there. The domain restriction $x>-1$ must be checked before any expansion is used.

### 9.4 Exact minimization of a convex quadratic

For a quadratic, the Taylor expansion around a stationary point has no remainder. This makes it possible to prove global optimality with one exact identity. Checking positive definiteness remains necessary: solving the gradient equation alone only finds a stationary point.

#### Problem AN-37: two-variable convex quadratic, case 1

Minimize $f(x)=\tfrac12x^\top Hx-b^\top x$ for $H=\begin{pmatrix}2&1\\1&3\end{pmatrix}$ and $b=\begin{pmatrix}1\\2\end{pmatrix}$. Prove existence and uniqueness and compute the minimum value.

Solution. The gradient is $\nabla f(x)=Hx-b$ because $H$ is symmetric. The Hessian is the constant matrix $H$. Its first leading principal minor is $2>0$ and its determinant is $5>0$. Completing a square in the quadratic form gives

$$
h^\top Hh
=2\left(h_1+\frac{1}{2}h_2\right)^2
+\frac{5}{2}h_2^2>0
\quad(h\ne0).
$$

Thus the objective is strictly convex. The stationary equation $Hx=b$ has the solution

$$
x^*=H^{-1}b=\begin{pmatrix}\frac{1}{5}\\\frac{3}{5}\end{pmatrix},\qquad
H^{-1}=\begin{pmatrix}\frac{3}{5}&-\frac{1}{5}\\-\frac{1}{5}&\frac{2}{5}\end{pmatrix}.
$$

Multiplying $H$ by the proposed vector checks both component equations. At this vector, $b=Hx^*$, so the objective simplifies to

$$
f(x^*)=-\frac12b^\top x^*=-\frac{7}{10}.
$$

Now take an arbitrary vector $h$ and expand $f(x^*+h)$. The terms linear in $h$ equal $h^\top(Hx^*-b)=0$. Therefore

$$
f(x^*+h)-f(x^*)=\frac12h^\top Hh.
$$

The right side is positive for every nonzero $h$. This proves that the stationary solution exists, is a global minimizer, and is the only minimizer. It also gives the exact cost of any coefficient error.

Positive definiteness implies coercivity here: as $\|x\|_2$ grows, the positive quadratic term eventually dominates the linear term. A continuous coercive function attains a minimum on some sufficiently large closed ball. The explicit stationary solution agrees with that existence argument.

#### Problem AN-38: two-variable convex quadratic, case 2

Minimize $f(x)=\tfrac12x^\top Hx-b^\top x$ for $H=\begin{pmatrix}4&-1\\-1&2\end{pmatrix}$ and $b=\begin{pmatrix}3\\-2\end{pmatrix}$. Prove existence and uniqueness and compute the minimum value.

Solution. The gradient is $\nabla f(x)=Hx-b$ because $H$ is symmetric. The Hessian is the constant matrix $H$. Its first leading principal minor is $4>0$ and its determinant is $7>0$. Completing a square in the quadratic form gives

$$
h^\top Hh
=4\left(h_1+\frac{-1}{4}h_2\right)^2
+\frac{7}{4}h_2^2>0
\quad(h\ne0).
$$

Thus the objective is strictly convex. The stationary equation $Hx=b$ has the solution

$$
x^*=H^{-1}b=\begin{pmatrix}\frac{4}{7}\\-\frac{5}{7}\end{pmatrix},\qquad
H^{-1}=\begin{pmatrix}\frac{2}{7}&\frac{1}{7}\\\frac{1}{7}&\frac{4}{7}\end{pmatrix}.
$$

Multiplying $H$ by the proposed vector checks both component equations. At this vector, $b=Hx^*$, so the objective simplifies to

$$
f(x^*)=-\frac12b^\top x^*=-\frac{11}{7}.
$$

Now take an arbitrary vector $h$ and expand $f(x^*+h)$. The terms linear in $h$ equal $h^\top(Hx^*-b)=0$. Therefore

$$
f(x^*+h)-f(x^*)=\frac12h^\top Hh.
$$

The right side is positive for every nonzero $h$. This proves that the stationary solution exists, is a global minimizer, and is the only minimizer. It also gives the exact cost of any coefficient error.

Positive definiteness implies coercivity here: as $\|x\|_2$ grows, the positive quadratic term eventually dominates the linear term. A continuous coercive function attains a minimum on some sufficiently large closed ball. The explicit stationary solution agrees with that existence argument.

#### Problem AN-39: two-variable convex quadratic, case 3

Minimize $f(x)=\tfrac12x^\top Hx-b^\top x$ for $H=\begin{pmatrix}3&0\\0&5\end{pmatrix}$ and $b=\begin{pmatrix}-1\\4\end{pmatrix}$. Prove existence and uniqueness and compute the minimum value.

Solution. The gradient is $\nabla f(x)=Hx-b$ because $H$ is symmetric. The Hessian is the constant matrix $H$. Its first leading principal minor is $3>0$ and its determinant is $15>0$. Completing a square in the quadratic form gives

$$
h^\top Hh
=3\left(h_1+\frac{0}{3}h_2\right)^2
+\frac{15}{3}h_2^2>0
\quad(h\ne0).
$$

Thus the objective is strictly convex. The stationary equation $Hx=b$ has the solution

$$
x^*=H^{-1}b=\begin{pmatrix}-\frac{1}{3}\\\frac{4}{5}\end{pmatrix},\qquad
H^{-1}=\begin{pmatrix}\frac{1}{3}&0\\0&\frac{1}{5}\end{pmatrix}.
$$

Multiplying $H$ by the proposed vector checks both component equations. At this vector, $b=Hx^*$, so the objective simplifies to

$$
f(x^*)=-\frac12b^\top x^*=-\frac{53}{30}.
$$

Now take an arbitrary vector $h$ and expand $f(x^*+h)$. The terms linear in $h$ equal $h^\top(Hx^*-b)=0$. Therefore

$$
f(x^*+h)-f(x^*)=\frac12h^\top Hh.
$$

The right side is positive for every nonzero $h$. This proves that the stationary solution exists, is a global minimizer, and is the only minimizer. It also gives the exact cost of any coefficient error.

Positive definiteness implies coercivity here: as $\|x\|_2$ grows, the positive quadratic term eventually dominates the linear term. A continuous coercive function attains a minimum on some sufficiently large closed ball. The explicit stationary solution agrees with that existence argument.

#### Problem AN-40: two-variable convex quadratic, case 4

Minimize $f(x)=\tfrac12x^\top Hx-b^\top x$ for $H=\begin{pmatrix}5&2\\2&2\end{pmatrix}$ and $b=\begin{pmatrix}2\\1\end{pmatrix}$. Prove existence and uniqueness and compute the minimum value.

Solution. The gradient is $\nabla f(x)=Hx-b$ because $H$ is symmetric. The Hessian is the constant matrix $H$. Its first leading principal minor is $5>0$ and its determinant is $6>0$. Completing a square in the quadratic form gives

$$
h^\top Hh
=5\left(h_1+\frac{2}{5}h_2\right)^2
+\frac{6}{5}h_2^2>0
\quad(h\ne0).
$$

Thus the objective is strictly convex. The stationary equation $Hx=b$ has the solution

$$
x^*=H^{-1}b=\begin{pmatrix}\frac{1}{3}\\\frac{1}{6}\end{pmatrix},\qquad
H^{-1}=\begin{pmatrix}\frac{1}{3}&-\frac{1}{3}\\-\frac{1}{3}&\frac{5}{6}\end{pmatrix}.
$$

Multiplying $H$ by the proposed vector checks both component equations. At this vector, $b=Hx^*$, so the objective simplifies to

$$
f(x^*)=-\frac12b^\top x^*=-\frac{5}{12}.
$$

Now take an arbitrary vector $h$ and expand $f(x^*+h)$. The terms linear in $h$ equal $h^\top(Hx^*-b)=0$. Therefore

$$
f(x^*+h)-f(x^*)=\frac12h^\top Hh.
$$

The right side is positive for every nonzero $h$. This proves that the stationary solution exists, is a global minimizer, and is the only minimizer. It also gives the exact cost of any coefficient error.

Positive definiteness implies coercivity here: as $\|x\|_2$ grows, the positive quadratic term eventually dominates the linear term. A continuous coercive function attains a minimum on some sufficiently large closed ball. The explicit stationary solution agrees with that existence argument.

### 9.5 Projection onto a half-space

The closest point in a closed half-space is either the original point, if feasible, or a point on the boundary. The multiplier determines how far to move in the normal direction. Complementary slackness distinguishes these two cases without guessing that every inequality must bind.

#### Problem AN-41: half-space projection, case 1

Minimize $\tfrac12\|x-z\|_2^2$ subject to $a^\top x\le b$, with $a=\begin{pmatrix}1\\2\end{pmatrix}$, $z=\begin{pmatrix}2\\2\end{pmatrix}$, and $b=3$.

Solution. The Lagrangian is $L(x,\lambda)=\tfrac12\|x-z\|_2^2+\lambda(a^\top x-b)$ with $\lambda\ge0$. Differentiation in $x$ gives $x-z+\lambda a=0$, hence $x=z-\lambda a$.

At the unconstrained minimum $z$, the constraint violation is $a^\top z-b=3$. If this is nonpositive, $z$ is feasible and no movement is needed. Otherwise the closest feasible point must lie on the boundary, and substituting stationarity into $a^\top x=b$ gives

$$
\lambda=\frac{a^\top z-b}{a^\top a}.
$$

Both cases are summarized by

$$
\lambda^*=\max\left(0,\frac{a^\top z-b}{a^\top a}\right)
=\frac{3}{5},\qquad
x^*=z-\lambda^*a=\begin{pmatrix}\frac{7}{5}\\\frac{4}{5}\end{pmatrix}.
$$

For these data, $a^\top a=5$, $a^\top x^*-b=0$, and the minimum objective is $\frac{9}{10}$. The multiplier is nonnegative, stationarity holds, and $\lambda^*(a^\top x^*-b)=0$.

A direct geometric proof gives another check. For any feasible $y$,

$$
\tfrac12\|y-z\|_2^2-\tfrac12\|x^*-z\|_2^2
=\tfrac12\|y-x^*\|_2^2
-\lambda^*a^\top(y-x^*).
$$

If the multiplier is positive, the boundary equality and feasibility make the last term nonnegative. If the multiplier is zero, that term vanishes. The difference is therefore nonnegative in either case, and is strictly positive for $y\ne x^*$. This proves global optimality and uniqueness.

#### Problem AN-42: half-space projection, case 2

Minimize $\tfrac12\|x-z\|_2^2$ subject to $a^\top x\le b$, with $a=\begin{pmatrix}2\\-1\end{pmatrix}$, $z=\begin{pmatrix}0\\0\end{pmatrix}$, and $b=1$.

Solution. The Lagrangian is $L(x,\lambda)=\tfrac12\|x-z\|_2^2+\lambda(a^\top x-b)$ with $\lambda\ge0$. Differentiation in $x$ gives $x-z+\lambda a=0$, hence $x=z-\lambda a$.

At the unconstrained minimum $z$, the constraint violation is $a^\top z-b=-1$. If this is nonpositive, $z$ is feasible and no movement is needed. Otherwise the closest feasible point must lie on the boundary, and substituting stationarity into $a^\top x=b$ gives

$$
\lambda=\frac{a^\top z-b}{a^\top a}.
$$

Both cases are summarized by

$$
\lambda^*=\max\left(0,\frac{a^\top z-b}{a^\top a}\right)
=0,\qquad
x^*=z-\lambda^*a=\begin{pmatrix}0\\0\end{pmatrix}.
$$

For these data, $a^\top a=5$, $a^\top x^*-b=-1$, and the minimum objective is $0$. The multiplier is nonnegative, stationarity holds, and $\lambda^*(a^\top x^*-b)=0$.

A direct geometric proof gives another check. For any feasible $y$,

$$
\tfrac12\|y-z\|_2^2-\tfrac12\|x^*-z\|_2^2
=\tfrac12\|y-x^*\|_2^2
-\lambda^*a^\top(y-x^*).
$$

If the multiplier is positive, the boundary equality and feasibility make the last term nonnegative. If the multiplier is zero, that term vanishes. The difference is therefore nonnegative in either case, and is strictly positive for $y\ne x^*$. This proves global optimality and uniqueness.

#### Problem AN-43: half-space projection, case 3

Minimize $\tfrac12\|x-z\|_2^2$ subject to $a^\top x\le b$, with $a=\begin{pmatrix}1\\1\end{pmatrix}$, $z=\begin{pmatrix}-2\\1\end{pmatrix}$, and $b=-2$.

Solution. The Lagrangian is $L(x,\lambda)=\tfrac12\|x-z\|_2^2+\lambda(a^\top x-b)$ with $\lambda\ge0$. Differentiation in $x$ gives $x-z+\lambda a=0$, hence $x=z-\lambda a$.

At the unconstrained minimum $z$, the constraint violation is $a^\top z-b=1$. If this is nonpositive, $z$ is feasible and no movement is needed. Otherwise the closest feasible point must lie on the boundary, and substituting stationarity into $a^\top x=b$ gives

$$
\lambda=\frac{a^\top z-b}{a^\top a}.
$$

Both cases are summarized by

$$
\lambda^*=\max\left(0,\frac{a^\top z-b}{a^\top a}\right)
=\frac{1}{2},\qquad
x^*=z-\lambda^*a=\begin{pmatrix}-\frac{5}{2}\\\frac{1}{2}\end{pmatrix}.
$$

For these data, $a^\top a=2$, $a^\top x^*-b=0$, and the minimum objective is $\frac{1}{4}$. The multiplier is nonnegative, stationarity holds, and $\lambda^*(a^\top x^*-b)=0$.

A direct geometric proof gives another check. For any feasible $y$,

$$
\tfrac12\|y-z\|_2^2-\tfrac12\|x^*-z\|_2^2
=\tfrac12\|y-x^*\|_2^2
-\lambda^*a^\top(y-x^*).
$$

If the multiplier is positive, the boundary equality and feasibility make the last term nonnegative. If the multiplier is zero, that term vanishes. The difference is therefore nonnegative in either case, and is strictly positive for $y\ne x^*$. This proves global optimality and uniqueness.

#### Problem AN-44: half-space projection, case 4

Minimize $\tfrac12\|x-z\|_2^2$ subject to $a^\top x\le b$, with $a=\begin{pmatrix}3\\1\end{pmatrix}$, $z=\begin{pmatrix}1\\2\end{pmatrix}$, and $b=0$.

Solution. The Lagrangian is $L(x,\lambda)=\tfrac12\|x-z\|_2^2+\lambda(a^\top x-b)$ with $\lambda\ge0$. Differentiation in $x$ gives $x-z+\lambda a=0$, hence $x=z-\lambda a$.

At the unconstrained minimum $z$, the constraint violation is $a^\top z-b=5$. If this is nonpositive, $z$ is feasible and no movement is needed. Otherwise the closest feasible point must lie on the boundary, and substituting stationarity into $a^\top x=b$ gives

$$
\lambda=\frac{a^\top z-b}{a^\top a}.
$$

Both cases are summarized by

$$
\lambda^*=\max\left(0,\frac{a^\top z-b}{a^\top a}\right)
=\frac{1}{2},\qquad
x^*=z-\lambda^*a=\begin{pmatrix}-\frac{1}{2}\\\frac{3}{2}\end{pmatrix}.
$$

For these data, $a^\top a=10$, $a^\top x^*-b=0$, and the minimum objective is $\frac{5}{4}$. The multiplier is nonnegative, stationarity holds, and $\lambda^*(a^\top x^*-b)=0$.

A direct geometric proof gives another check. For any feasible $y$,

$$
\tfrac12\|y-z\|_2^2-\tfrac12\|x^*-z\|_2^2
=\tfrac12\|y-x^*\|_2^2
-\lambda^*a^\top(y-x^*).
$$

If the multiplier is positive, the boundary equality and feasibility make the last term nonnegative. If the multiplier is zero, that term vanishes. The difference is therefore nonnegative in either case, and is strictly positive for $y\ne x^*$. This proves global optimality and uniqueness.

### 9.6 A dual function with active and inactive cases

The dual function is found by minimizing the Lagrangian over the primal variable without retaining the primal constraint. The multiplier domain remains $\lambda\ge0$. Solving this small example in both primal and dual variables makes the role of that domain explicit.

#### Problem AN-45: quadratic dual and sensitivity, case 1

Minimize $f(x)=\tfrac12(2)x^2+(-4)x$ subject to $x\ge 3$. Derive and maximize the dual function and interpret its multiplier.

Solution. The unconstrained minimizer solves $2x+(-4)=0$, giving $x_0=2$. The positive second derivative $2$ means that the constrained minimizer is the closest feasible point to $x_0$:

$$
x^*=\max(b,-a/q)=3,\qquad f(x^*)=-3.
$$

Write the constraint as $b-x\le0$. The Lagrangian and its unconstrained minimizer are

$$
L(x,\lambda)=\tfrac q2x^2+ax+\lambda(b-x),\qquad
x(\lambda)=\frac{\lambda-a}q.
$$

Completing the square gives

$$
g(\lambda)=\inf_xL(x,\lambda)
=\lambda b-\frac{(a-\lambda)^2}{2q},
\qquad \lambda\ge0.
$$

The derivative is $g'(\lambda)=b+(a-\lambda)/q$. The unrestricted maximum would occur at $qb+a$. Restricting to the nonnegative half-line gives $\lambda^*=\max(0,qb+a)=2$. Substituting it into $x(\lambda)$ recovers the primal solution. Direct substitution also gives $g(\lambda^*)=f(x^*)=-3$.

Thus the duality gap is zero. Complementary slackness $\lambda^*(b-x^*)=0$ says that a strictly slack lower bound has zero multiplier. A binding lower bound can also have zero multiplier exactly where it first meets the unconstrained optimum.

If the lower bound is treated as a parameter, the value function is constant below $-a/q$ and equals $qb^2/2+ab$ above it. Its derivative is $\max(0,qb+a)$. Increasing this bound tightens the feasible set, so a positive multiplier means a positive marginal increase in minimum cost. The sign follows from the chosen parameterization.

#### Problem AN-46: quadratic dual and sensitivity, case 2

Minimize $f(x)=\tfrac12(3)x^2+(6)x$ subject to $x\ge -1$. Derive and maximize the dual function and interpret its multiplier.

Solution. The unconstrained minimizer solves $3x+(6)=0$, giving $x_0=-2$. The positive second derivative $3$ means that the constrained minimizer is the closest feasible point to $x_0$:

$$
x^*=\max(b,-a/q)=-1,\qquad f(x^*)=-\frac{9}{2}.
$$

Write the constraint as $b-x\le0$. The Lagrangian and its unconstrained minimizer are

$$
L(x,\lambda)=\tfrac q2x^2+ax+\lambda(b-x),\qquad
x(\lambda)=\frac{\lambda-a}q.
$$

Completing the square gives

$$
g(\lambda)=\inf_xL(x,\lambda)
=\lambda b-\frac{(a-\lambda)^2}{2q},
\qquad \lambda\ge0.
$$

The derivative is $g'(\lambda)=b+(a-\lambda)/q$. The unrestricted maximum would occur at $qb+a$. Restricting to the nonnegative half-line gives $\lambda^*=\max(0,qb+a)=3$. Substituting it into $x(\lambda)$ recovers the primal solution. Direct substitution also gives $g(\lambda^*)=f(x^*)=-\frac{9}{2}$.

Thus the duality gap is zero. Complementary slackness $\lambda^*(b-x^*)=0$ says that a strictly slack lower bound has zero multiplier. A binding lower bound can also have zero multiplier exactly where it first meets the unconstrained optimum.

If the lower bound is treated as a parameter, the value function is constant below $-a/q$ and equals $qb^2/2+ab$ above it. Its derivative is $\max(0,qb+a)$. Increasing this bound tightens the feasible set, so a positive multiplier means a positive marginal increase in minimum cost. The sign follows from the chosen parameterization.

#### Problem AN-47: quadratic dual and sensitivity, case 3

Minimize $f(x)=\tfrac12(4)x^2+(-8)x$ subject to $x\ge 1$. Derive and maximize the dual function and interpret its multiplier.

Solution. The unconstrained minimizer solves $4x+(-8)=0$, giving $x_0=2$. The positive second derivative $4$ means that the constrained minimizer is the closest feasible point to $x_0$:

$$
x^*=\max(b,-a/q)=2,\qquad f(x^*)=-8.
$$

Write the constraint as $b-x\le0$. The Lagrangian and its unconstrained minimizer are

$$
L(x,\lambda)=\tfrac q2x^2+ax+\lambda(b-x),\qquad
x(\lambda)=\frac{\lambda-a}q.
$$

Completing the square gives

$$
g(\lambda)=\inf_xL(x,\lambda)
=\lambda b-\frac{(a-\lambda)^2}{2q},
\qquad \lambda\ge0.
$$

The derivative is $g'(\lambda)=b+(a-\lambda)/q$. The unrestricted maximum would occur at $qb+a$. Restricting to the nonnegative half-line gives $\lambda^*=\max(0,qb+a)=0$. Substituting it into $x(\lambda)$ recovers the primal solution. Direct substitution also gives $g(\lambda^*)=f(x^*)=-8$.

Thus the duality gap is zero. Complementary slackness $\lambda^*(b-x^*)=0$ says that a strictly slack lower bound has zero multiplier. A binding lower bound can also have zero multiplier exactly where it first meets the unconstrained optimum.

If the lower bound is treated as a parameter, the value function is constant below $-a/q$ and equals $qb^2/2+ab$ above it. Its derivative is $\max(0,qb+a)$. Increasing this bound tightens the feasible set, so a positive multiplier means a positive marginal increase in minimum cost. The sign follows from the chosen parameterization.

#### Problem AN-48: quadratic dual and sensitivity, case 4

Minimize $f(x)=\tfrac12(1)x^2+(2)x$ subject to $x\ge -3$. Derive and maximize the dual function and interpret its multiplier.

Solution. The unconstrained minimizer solves $1x+(2)=0$, giving $x_0=-2$. The positive second derivative $1$ means that the constrained minimizer is the closest feasible point to $x_0$:

$$
x^*=\max(b,-a/q)=-2,\qquad f(x^*)=-2.
$$

Write the constraint as $b-x\le0$. The Lagrangian and its unconstrained minimizer are

$$
L(x,\lambda)=\tfrac q2x^2+ax+\lambda(b-x),\qquad
x(\lambda)=\frac{\lambda-a}q.
$$

Completing the square gives

$$
g(\lambda)=\inf_xL(x,\lambda)
=\lambda b-\frac{(a-\lambda)^2}{2q},
\qquad \lambda\ge0.
$$

The derivative is $g'(\lambda)=b+(a-\lambda)/q$. The unrestricted maximum would occur at $qb+a$. Restricting to the nonnegative half-line gives $\lambda^*=\max(0,qb+a)=0$. Substituting it into $x(\lambda)$ recovers the primal solution. Direct substitution also gives $g(\lambda^*)=f(x^*)=-2$.

Thus the duality gap is zero. Complementary slackness $\lambda^*(b-x^*)=0$ says that a strictly slack lower bound has zero multiplier. A binding lower bound can also have zero multiplier exactly where it first meets the unconstrained optimum.

If the lower bound is treated as a parameter, the value function is constant below $-a/q$ and equals $qb^2/2+ab$ above it. Its derivative is $\max(0,qb+a)$. Increasing this bound tightens the feasible set, so a positive multiplier means a positive marginal increase in minimum cost. The sign follows from the chosen parameterization.

### 9.7 Box uncertainty and sign-dependent worst cases

When coefficients vary independently in intervals, the worst coefficient for a positive decision coordinate is its upper endpoint. For a negative coordinate it is its lower endpoint. The support function formula records both choices through absolute values and avoids silently assuming nonnegative decisions.

#### Problem AN-49: robust scaling in a fixed direction, case 1

Choose $x=td$ with $t\ge0$, $d=\begin{pmatrix}1\\1\end{pmatrix}$, as large as possible subject to $a^\top x\le 12$ for every $a_i\in[a_{0i}-r_i,a_{0i}+r_i]$. Here $a_0=\begin{pmatrix}2\\1\end{pmatrix}$ and $r=\begin{pmatrix}1\\2\end{pmatrix}$.

Solution. Write $a=a_0+u$ with $|u_i|\le r_i$. For a fixed $x$,

$$
\sup_a a^\top x=a_0^\top x+\sum_i r_i|x_i|.
$$

To prove the formula, bound each product by $u_ix_i\le r_i|x_i|$. Equality is attained by choosing $u_i=r_i$ when $x_i>0$ and $u_i=-r_i$ when $x_i<0$; a zero coordinate allows either endpoint.

Since $x=td$ with $t\ge0$, the worst-case left side is

$$
t\left(a_0^\top d+\sum_i r_i|d_i|\right)
=t(6).
$$

The coefficient is strictly positive for the supplied data. The robust feasible interval is therefore $0\le t\le 2$, and the largest scaling is $t^*=2$. The resulting decision is $x^*=\begin{pmatrix}2\\2\end{pmatrix}$.

One worst-case coefficient vector is $a^*=\begin{pmatrix}3\\3\end{pmatrix}$. Checking its inner product gives $(a^*)^\top x^*=12$, so the robust constraint binds. Every other allowed coefficient vector has an inner product no larger, by the support-function inequality above.

The absolute values are essential when a coordinate of $d$ is negative. Replacing the uncertain vector by all upper endpoints would generally fail in that case. If the computed support coefficient had been zero or negative while the right side remained positive, this particular constraint would impose no finite upper bound on $t$; dividing by it as though positive would be invalid.

#### Problem AN-50: robust scaling in a fixed direction, case 2

Choose $x=td$ with $t\ge0$, $d=\begin{pmatrix}1\\-2\end{pmatrix}$, as large as possible subject to $a^\top x\le 14$ for every $a_i\in[a_{0i}-r_i,a_{0i}+r_i]$. Here $a_0=\begin{pmatrix}3\\-1\end{pmatrix}$ and $r=\begin{pmatrix}1\\1\end{pmatrix}$.

Solution. Write $a=a_0+u$ with $|u_i|\le r_i$. For a fixed $x$,

$$
\sup_a a^\top x=a_0^\top x+\sum_i r_i|x_i|.
$$

To prove the formula, bound each product by $u_ix_i\le r_i|x_i|$. Equality is attained by choosing $u_i=r_i$ when $x_i>0$ and $u_i=-r_i$ when $x_i<0$; a zero coordinate allows either endpoint.

Since $x=td$ with $t\ge0$, the worst-case left side is

$$
t\left(a_0^\top d+\sum_i r_i|d_i|\right)
=t(8).
$$

The coefficient is strictly positive for the supplied data. The robust feasible interval is therefore $0\le t\le \frac{7}{4}$, and the largest scaling is $t^*=\frac{7}{4}$. The resulting decision is $x^*=\begin{pmatrix}\frac{7}{4}\\-\frac{7}{2}\end{pmatrix}$.

One worst-case coefficient vector is $a^*=\begin{pmatrix}4\\-2\end{pmatrix}$. Checking its inner product gives $(a^*)^\top x^*=14$, so the robust constraint binds. Every other allowed coefficient vector has an inner product no larger, by the support-function inequality above.

The absolute values are essential when a coordinate of $d$ is negative. Replacing the uncertain vector by all upper endpoints would generally fail in that case. If the computed support coefficient had been zero or negative while the right side remained positive, this particular constraint would impose no finite upper bound on $t$; dividing by it as though positive would be invalid.

#### Problem AN-51: robust scaling in a fixed direction, case 3

Choose $x=td$ with $t\ge0$, $d=\begin{pmatrix}-1\\2\end{pmatrix}$, as large as possible subject to $a^\top x\le 15$ for every $a_i\in[a_{0i}-r_i,a_{0i}+r_i]$. Here $a_0=\begin{pmatrix}1\\2\end{pmatrix}$ and $r=\begin{pmatrix}2\\1\end{pmatrix}$.

Solution. Write $a=a_0+u$ with $|u_i|\le r_i$. For a fixed $x$,

$$
\sup_a a^\top x=a_0^\top x+\sum_i r_i|x_i|.
$$

To prove the formula, bound each product by $u_ix_i\le r_i|x_i|$. Equality is attained by choosing $u_i=r_i$ when $x_i>0$ and $u_i=-r_i$ when $x_i<0$; a zero coordinate allows either endpoint.

Since $x=td$ with $t\ge0$, the worst-case left side is

$$
t\left(a_0^\top d+\sum_i r_i|d_i|\right)
=t(7).
$$

The coefficient is strictly positive for the supplied data. The robust feasible interval is therefore $0\le t\le \frac{15}{7}$, and the largest scaling is $t^*=\frac{15}{7}$. The resulting decision is $x^*=\begin{pmatrix}-\frac{15}{7}\\\frac{30}{7}\end{pmatrix}$.

One worst-case coefficient vector is $a^*=\begin{pmatrix}-1\\3\end{pmatrix}$. Checking its inner product gives $(a^*)^\top x^*=15$, so the robust constraint binds. Every other allowed coefficient vector has an inner product no larger, by the support-function inequality above.

The absolute values are essential when a coordinate of $d$ is negative. Replacing the uncertain vector by all upper endpoints would generally fail in that case. If the computed support coefficient had been zero or negative while the right side remained positive, this particular constraint would impose no finite upper bound on $t$; dividing by it as though positive would be invalid.

#### Problem AN-52: robust scaling in a fixed direction, case 4

Choose $x=td$ with $t\ge0$, $d=\begin{pmatrix}-2\\3\end{pmatrix}$, as large as possible subject to $a^\top x\le 20$ for every $a_i\in[a_{0i}-r_i,a_{0i}+r_i]$. Here $a_0=\begin{pmatrix}2\\3\end{pmatrix}$ and $r=\begin{pmatrix}1\\2\end{pmatrix}$.

Solution. Write $a=a_0+u$ with $|u_i|\le r_i$. For a fixed $x$,

$$
\sup_a a^\top x=a_0^\top x+\sum_i r_i|x_i|.
$$

To prove the formula, bound each product by $u_ix_i\le r_i|x_i|$. Equality is attained by choosing $u_i=r_i$ when $x_i>0$ and $u_i=-r_i$ when $x_i<0$; a zero coordinate allows either endpoint.

Since $x=td$ with $t\ge0$, the worst-case left side is

$$
t\left(a_0^\top d+\sum_i r_i|d_i|\right)
=t(13).
$$

The coefficient is strictly positive for the supplied data. The robust feasible interval is therefore $0\le t\le \frac{20}{13}$, and the largest scaling is $t^*=\frac{20}{13}$. The resulting decision is $x^*=\begin{pmatrix}-\frac{40}{13}\\\frac{60}{13}\end{pmatrix}$.

One worst-case coefficient vector is $a^*=\begin{pmatrix}1\\5\end{pmatrix}$. Checking its inner product gives $(a^*)^\top x^*=20$, so the robust constraint binds. Every other allowed coefficient vector has an inner product no larger, by the support-function inequality above.

The absolute values are essential when a coordinate of $d$ is negative. Replacing the uncertain vector by all upper endpoints would generally fail in that case. If the computed support coefficient had been zero or negative while the right side remained positive, this particular constraint would impose no finite upper bound on $t$; dividing by it as though positive would be invalid.

### 9.8 Maximizing a minimum of two affine functions

The lower of an increasing line and a decreasing line rises until the crossing and falls afterwards. An interval restriction can put the crossing outside the feasible set. The optimizer is then an endpoint. The auxiliary variable formulation gives the same result as this direct argument.

#### Problem AN-53: two-reward max-min problem, case 1

Maximize $\min\{1x+(0),\,2-(2)x\}$ over $0\le x\le 1$.

Solution. Set the two affine rewards equal to locate their crossing:

$$
1x+(0)=2-(2)x,
\qquad x_c=\frac{v-u}{p+q}=\frac{2}{3}.
$$

The slopes have opposite signs because $p=1>0$ and $q=2>0$. To the left of the crossing, the increasing reward is the smaller one. To the right, the decreasing reward is the smaller one. The objective is therefore increasing before $x_c$ and decreasing after it.

Restricting to the supplied interval gives the clipped solution

$$
x^*=\min(U,\max(L,x_c))=\frac{2}{3}.
$$

The two rewards there are $\frac{2}{3}$ and $\frac{2}{3}$, so the maximum minimum reward is $\frac{2}{3}$. If the crossing is outside the interval, the two rewards need not be equal at the optimum. Equalizing them without first checking feasibility would produce an invalid decision.

The equivalent linear optimization problem introduces a reward floor $z$:

$$
\max_{x,z}z
\quad\text{subject to}\quad
z\le px+u,\quad z\le v-qx,\quad L\le x\le U.
$$

For each fixed $x$, the largest feasible $z$ is exactly the smaller reward, proving equivalence. This formulation remains useful with more than two affine rewards even when a single crossing argument no longer describes the entire solution. Here the direct monotonicity proof already establishes global optimality and uniqueness.

#### Problem AN-54: two-reward max-min problem, case 2

Maximize $\min\{2x+(1),\,7-(1)x\}$ over $0\le x\le 1$.

Solution. Set the two affine rewards equal to locate their crossing:

$$
2x+(1)=7-(1)x,
\qquad x_c=\frac{v-u}{p+q}=2.
$$

The slopes have opposite signs because $p=2>0$ and $q=1>0$. To the left of the crossing, the increasing reward is the smaller one. To the right, the decreasing reward is the smaller one. The objective is therefore increasing before $x_c$ and decreasing after it.

Restricting to the supplied interval gives the clipped solution

$$
x^*=\min(U,\max(L,x_c))=1.
$$

The two rewards there are $3$ and $6$, so the maximum minimum reward is $3$. If the crossing is outside the interval, the two rewards need not be equal at the optimum. Equalizing them without first checking feasibility would produce an invalid decision.

The equivalent linear optimization problem introduces a reward floor $z$:

$$
\max_{x,z}z
\quad\text{subject to}\quad
z\le px+u,\quad z\le v-qx,\quad L\le x\le U.
$$

For each fixed $x$, the largest feasible $z$ is exactly the smaller reward, proving equivalence. This formulation remains useful with more than two affine rewards even when a single crossing argument no longer describes the entire solution. Here the direct monotonicity proof already establishes global optimality and uniqueness.

#### Problem AN-55: two-reward max-min problem, case 3

Maximize $\min\{3x+(-2),\,3-(2)x\}$ over $-1\le x\le 2$.

Solution. Set the two affine rewards equal to locate their crossing:

$$
3x+(-2)=3-(2)x,
\qquad x_c=\frac{v-u}{p+q}=1.
$$

The slopes have opposite signs because $p=3>0$ and $q=2>0$. To the left of the crossing, the increasing reward is the smaller one. To the right, the decreasing reward is the smaller one. The objective is therefore increasing before $x_c$ and decreasing after it.

Restricting to the supplied interval gives the clipped solution

$$
x^*=\min(U,\max(L,x_c))=1.
$$

The two rewards there are $1$ and $1$, so the maximum minimum reward is $1$. If the crossing is outside the interval, the two rewards need not be equal at the optimum. Equalizing them without first checking feasibility would produce an invalid decision.

The equivalent linear optimization problem introduces a reward floor $z$:

$$
\max_{x,z}z
\quad\text{subject to}\quad
z\le px+u,\quad z\le v-qx,\quad L\le x\le U.
$$

For each fixed $x$, the largest feasible $z$ is exactly the smaller reward, proving equivalence. This formulation remains useful with more than two affine rewards even when a single crossing argument no longer describes the entire solution. Here the direct monotonicity proof already establishes global optimality and uniqueness.

#### Problem AN-56: two-reward max-min problem, case 4

Maximize $\min\{1x+(4),\,0-(3)x\}$ over $0\le x\le 2$.

Solution. Set the two affine rewards equal to locate their crossing:

$$
1x+(4)=0-(3)x,
\qquad x_c=\frac{v-u}{p+q}=-1.
$$

The slopes have opposite signs because $p=1>0$ and $q=3>0$. To the left of the crossing, the increasing reward is the smaller one. To the right, the decreasing reward is the smaller one. The objective is therefore increasing before $x_c$ and decreasing after it.

Restricting to the supplied interval gives the clipped solution

$$
x^*=\min(U,\max(L,x_c))=0.
$$

The two rewards there are $4$ and $0$, so the maximum minimum reward is $0$. If the crossing is outside the interval, the two rewards need not be equal at the optimum. Equalizing them without first checking feasibility would produce an invalid decision.

The equivalent linear optimization problem introduces a reward floor $z$:

$$
\max_{x,z}z
\quad\text{subject to}\quad
z\le px+u,\quad z\le v-qx,\quad L\le x\le U.
$$

For each fixed $x$, the largest feasible $z$ is exactly the smaller reward, proving equivalence. This formulation remains useful with more than two affine rewards even when a single crossing argument no longer describes the entire solution. Here the direct monotonicity proof already establishes global optimality and uniqueness.

<!-- further-distinct-problems-AN -->

### 9.9 Further proofs and constrained problems

These additional problems are constructed for these notes.

#### Problem AN-57: uniform convergence and interchange of limits

For $f_n(x)=x^n$ on $[0,1]$, find the pointwise limit and determine uniform convergence on this interval and on $[0,a]$ for fixed $a<1$. Prove that a uniform limit of continuous functions is continuous. Give a separate example showing that uniform convergence of functions does not ensure convergence of their derivatives to the derivative of the limit.

Solution. For $0\le x<1$, powers tend to zero; at one they remain one. The pointwise limit is therefore zero except for value one at the right endpoint. Its discontinuity already suggests that convergence is not uniform.

Directly, the supremum of $|f_n-f|$ over points below one is one: values can be made arbitrarily close to one by choosing $x<1$ sufficiently near one. This supremum does not decrease with $n$. On $[0,a]$, however, the supremum is $a^n\to0$, proving uniform convergence there.

To prove the continuous-limit theorem, suppose $f_n\to f$ uniformly and every $f_n$ is continuous. Fix a point $x_0$ and tolerance $\epsilon>0$. Choose one index $N$ for which $\sup_x|f_N(x)-f(x)|<\epsilon/3$. Continuity of this fixed $f_N$ supplies a neighborhood in which $|f_N(x)-f_N(x_0)|<\epsilon/3$. The triangle inequality gives

$$
|f(x)-f(x_0)|
\le |f(x)-f_N(x)|+|f_N(x)-f_N(x_0)|+|f_N(x_0)-f(x_0)|
<\epsilon.
$$

The uniform choice of $N$ controls both evaluation points. Pointwise convergence would not give the same control throughout the neighborhood.

For derivatives, take $g_n(x)=\sin(nx)/n$ on the real line. The bound $|g_n(x)|\le1/n$ proves uniform convergence to zero. Yet $g_n'(0)=1$ for every $n$, while the derivative of the limit is zero. Stronger assumptions, such as suitable uniform convergence of derivatives together with a value at one point, are needed for a differentiation-under-limit theorem.

For $x^n$, the integrals equal $1/(n+1)$ and converge to zero, agreeing with the integral of the pointwise limit. Nonuniform convergence prevents a particular general theorem from being applied; it does not force every possible interchange to fail.

#### Problem AN-58: a contraction equation with an explicit error bound

Prove that $x=\cos x$ has exactly one solution in $[0,1]$ and that iteration $x_{n+1}=\cos x_n$ converges from every point of that interval. Derive a computable error bound from successive iterates.

Solution. The cosine maps $[0,1]$ into $[\cos1,1]$, which lies within the same interval. Its derivative magnitude is $\sin x\le q=\sin1<1$ there. The mean value theorem gives $|\cos x-\cos y|\le q|x-y|$.

Consecutive differences therefore satisfy $|x_{n+1}-x_n|\le q^n|x_1-x_0|$. For $m>n$, summing these differences gives

$$
|x_m-x_n|
\le\sum_{k=n}^{m-1}|x_{k+1}-x_k|
\le\frac{q^n}{1-q}|x_1-x_0|.
$$

The right side tends to zero, so the sequence is Cauchy. Completeness gives a limit, and closedness of the interval keeps it in $[0,1]$. Continuity of cosine then makes the limit a fixed point.

If $r,s$ were two fixed points, $|r-s|=|\cos r-\cos s|\le q|r-s|$. Since $q<1$, this forces equality of the points. The root is approximately $0.7390851332$.

Letting $m\to\infty$ in the tail estimate gives the a priori error bound $|x_n-r|\le q^n|x_1-x_0|/(1-q)$. Repeating the argument from the latest known difference gives

$$
|x_n-r|\le\frac q{1-q}|x_n-x_{n-1}|.
$$

This is a valid stopping certificate on this invariant contraction interval. An arbitrary iteration with small successive steps need not have such a certificate. Here the derivative bound and invariance are what convert a step difference into a bound on the unseen remaining tail.

#### Problem AN-59: interior critical points and boundary extrema

Find every global minimum and maximum of $f(x)=x^3-3x$ on $[-2,2]$. Classify its interior stationary points and explain why the endpoints must be checked separately.

Solution. The polynomial is continuous on a compact nonempty interval, so both extrema exist. Interior extrema must satisfy $f'(x)=3x^2-3=0$, giving $x=-1,1$. The second derivative is $6x$, so minus one is a strict local maximum and one is a strict local minimum.

The candidate list includes both endpoints because an endpoint optimum need not have derivative zero. Evaluate all four values:

$$
f(-2)=-2,\qquad f(-1)=2,\qquad f(1)=-2,\qquad f(2)=2.
$$

The global minimum is minus two, attained at $-2$ and $1$. The global maximum is two, attained at $-1$ and $2$.

For completeness, the derivative is positive where $|x|>1$ and negative where $|x|<1$. Thus the function increases on $[-2,-1]$, decreases on $[-1,1]$, and increases on $[1,2]$. These monotonic pieces also prove that the candidate list contains every extremum.

The second-derivative test classified the interior points locally. The value comparison established which ones are global and found the tied boundary extrema. Neither step alone would provide the complete answer.

#### Problem AN-60: an absolute-value penalty

For real $a$ and $\lambda\ge0$, minimize $F(x)=\tfrac12(x-a)^2+\lambda|x|$. Derive the optimizer for every $a$ and justify the case at zero.

Solution. On $x>0$, differentiation gives $F'(x)=x-a+\lambda$. Its stationary point is $a-\lambda$ and belongs to that region only if $a>\lambda$. On $x<0$, the derivative is $x-a-\lambda$, giving $a+\lambda$ only when $a<-\lambda$.

At zero, the left derivative is $-a-\lambda$ and the right derivative is $-a+\lambda$. A minimum there requires the left derivative to be nonpositive and the right derivative nonnegative, equivalent to $|a|\le\lambda$.

Thus the solution is

$$
x^*=\begin{cases}
a-\lambda,&a>\lambda,\\
0,&|a|\le\lambda,\\
a+\lambda,&a<-\lambda.
\end{cases}
$$

The quadratic part is strictly convex and the absolute value is convex, so their sum is strictly convex. Every candidate meeting the appropriate one-sided condition is therefore the unique global minimum.

This can also be expressed using a subgradient. At zero, the subgradients of $|x|$ form $[-1,1]$ because $|y|\ge sy$ for every $y$ exactly when $|s|\le1$. The optimality condition $0\in x-a+\lambda\partial|x|$ gives all three cases in one statement.

With $\lambda=1$, inputs $a=3,1/2,-3$ produce solutions $2,0,-2$. An interval of inputs is mapped exactly to zero. For a vector problem with a sum of separate squared errors and absolute penalties, the objective separates into coordinates and this formula applies to each. Correlated design columns generally destroy that coordinatewise closed form.

#### Problem AN-61: a disk constraint and an exact dual certificate

Minimize $(x-2)^2+(y-1)^2$ subject to $x^2+y^2\le1$. Derive the primal solution, all needed KKT quantities, and a matching dual value.

Solution. Write $z=(2,1)^\top$, whose norm is $\sqrt5>1$. The unconstrained minimizer lies outside the disk. With multiplier $\lambda\ge0$, the Lagrangian is

$$
L(w,\lambda)=(1+\lambda)\|w\|^2-2z^\top w+5-\lambda.
$$

Stationarity gives $w=z/(1+\lambda)$. If the constraint were slack, complementary slackness would give $\lambda=0$ and $w=z$, which is infeasible. Thus the constraint binds, yielding $1+\lambda=\sqrt5$.

The solution is $w^*=(2,1)^\top/\sqrt5$, with multiplier $\lambda^*=\sqrt5-1>0$. Its cost is the squared radial distance $(\sqrt5-1)^2=6-2\sqrt5$. Feasibility, multiplier nonnegativity, stationarity, and complementary slackness all hold.

Minimizing the displayed Lagrangian over unrestricted $w$ gives

$$
g(\lambda)=5-\lambda-\frac5{1+\lambda}.
$$

Its derivative is $-1+5/(1+\lambda)^2$ and its second derivative is $-10/(1+\lambda)^3<0$. The unique dual maximizer is the same multiplier, with $g(\lambda^*)=6-2\sqrt5$.

The matching values certify optimality. The primal objective is strictly convex and the disk is convex, so the minimizer is unique. The center of the disk is strictly feasible, also verifying Slater's condition. The multiplier is attached to the squared-norm constraint as written; using $\|w\|-1\le0$ would produce a differently scaled multiplier while leaving the optimizer unchanged.

#### Problem AN-62: an equality with an active nonnegativity bound

Minimize $\tfrac12(x-3)^2+\tfrac12(y+1)^2$ subject to $x+y=1$ and $x,y\ge0$. Solve by elimination and give a KKT certificate with explicit multiplier signs.

Solution. Eliminate $y=1-x$. Nonnegativity restricts $x$ to $[0,1]$. The scalar objective is $\tfrac12(x-3)^2+\tfrac12(2-x)^2$ with derivative $2x-5$. It decreases throughout that interval, so $x^*=1$ and $y^*=0$. The minimum value is $5/2$.

Use $L=f+\nu(x+y-1)-\alpha x-\gamma y$, with $\alpha,\gamma\ge0$ and unrestricted equality multiplier $\nu$. Stationarity at the proposed solution gives

$$
-2+\nu-\alpha=0,\qquad
1+\nu-\gamma=0.
$$

Since $x^*>0$, complementary slackness requires $\alpha=0$. Therefore $\nu=2$ and $\gamma=3$. The second nonnegativity constraint binds and has positive multiplier. All signs match the convention of writing its inequality as $-y\le0$.

The objective Hessian is identity and the feasible set is convex, so KKT is sufficient and the solution is unique. Existence also follows directly from continuity on the compact line segment.

For a sensitivity check, replace the equality's right side by $b$ near one. The active set remains $x=b,y=0$, giving value $p(b)=\tfrac12(b-3)^2+\tfrac12$. Hence $p'(1)=-2=-\nu$. An equality multiplier is not constrained to the sign used for inequality multipliers; its sensitivity sign comes from the term $-\nu b$.

#### Problem AN-63: weighted worst-case squared errors

Minimize $\max\{x^2,4(x-3)^2\}$ over all real $x$. Prove a global lower bound and derive the multipliers in an auxiliary-variable formulation.

Solution. Let the worst-case value at $x$ be $t$. Then $|x|\le\sqrt t$ and $|x-3|\le\sqrt t/2$. The triangle inequality gives

$$
3\le|x|+|x-3|\le\frac32\sqrt t.
$$

Thus $t\ge4$ for every real decision. At $x=2$, the two squared errors both equal four, so the bound is attained and the optimum is four.

To identify uniqueness, equality requires distances two from zero and one from three, with $x$ between the endpoints. This forces $x=2$.

The epigraph problem minimizes $t$ subject to $x^2-t\le0$ and $4(x-3)^2-t\le0$. Its Lagrangian is

$$
L=t+\lambda_1(x^2-t)+\lambda_2(4(x-3)^2-t).
$$

Stationarity in $t$ gives $\lambda_1+\lambda_2=1$. At $x=2$, stationarity in $x$ gives $4\lambda_1-8\lambda_2=0$. Therefore $\lambda_1=2/3$ and $\lambda_2=1/3$. They are nonnegative, both inequalities bind, and all KKT conditions hold.

The auxiliary problem is convex and has a strictly feasible point, such as $x=2,t=5$. The multipliers can be read as weights on the two active error functions. They are not equal because the two functions have different slopes at the balancing point. Equal worst-case losses do not imply equal dual weights.

#### Problem AN-64: a continuous value with a discontinuous optimizer

Let $p(a)=\min_{-1\le x\le1}ax$. Find the value and all minimizers. Then add $\epsilon x^2/2$ with $\epsilon>0$ and determine how the solution and value change.

Solution. For $a>0$ the objective increases in $x$, so the minimizer is minus one. For $a<0$ it decreases, so the minimizer is one. At $a=0$ every feasible point is optimal. Thus $p(a)=-|a|$. The value is continuous and Lipschitz with constant one, but is not differentiable at zero. Every optimizing selection must jump across zero.

For the regularized objective $ax+\epsilon x^2/2$, the unrestricted stationary point is $-a/\epsilon$. Restricting to the interval gives

$$
x_\epsilon(a)=\min(1,\max(-1,-a/\epsilon)).
$$

Strict convexity makes this optimizer unique. Substitution gives

$$
p_\epsilon(a)=
\begin{cases}
-a^2/(2\epsilon),&|a|\le\epsilon,\\
\epsilon/2-|a|,&|a|\ge\epsilon.
\end{cases}
$$

The expressions agree at the boundary. Their derivatives also agree, and $p_\epsilon'(a)=x_\epsilon(a)$. Thus the regularized value is differentiable everywhere.

The penalty lies between zero and $\epsilon/2$ throughout the feasible interval, so $0\le p_\epsilon(a)-p(a)\le\epsilon/2$ uniformly in $a$. Uniform convergence of smooth value functions has produced a nonsmooth limiting value. This agrees with the earlier counterexample about derivatives: convergence of functions alone does not justify interchanging differentiation and a limit.

#### Problem AN-65: weighted arithmetic and geometric means

For positive numbers $x_i$ and nonnegative weights $w_i$ summing to one, prove $\prod_i x_i^{w_i}\le\sum_iw_ix_i$ and characterize equality. Apply it to maximizing a product subject to a fixed positive sum.

Solution. The logarithm has second derivative $-1/x^2<0$ on its positive domain, so it is strictly concave. Finite Jensen applied to it gives

$$
\sum_iw_i\log x_i\le
\log\left(\sum_iw_ix_i\right).
$$

Exponentiation preserves order and turns the left side into the weighted product. Equality holds exactly when all $x_i$ carrying positive weights are equal; zero-weight entries do not affect either side.

For $n$ positive numbers with sum $B>0$, take equal weights $1/n$. Then $(\prod_i x_i)^{1/n}\le B/n$, so the largest possible product is $(B/n)^n$, attained only at $x_i=B/n$.

Existence can also be checked by extending the product continuously to the compact nonnegative simplex. Every boundary point has product zero, while the equal positive allocation has positive product. Therefore the maximum is interior, making the equality condition applicable.

For a weighted log objective $\sum_iw_i\log x_i$ with strictly positive weights, the same budget instead gives stationarity $w_i/x_i=\lambda$ and hence $x_i=w_iB$. This is a different allocation problem from the unweighted product. The sum constraint binds by monotonicity, strict concavity gives uniqueness, and $\lambda=1/B$ follows from summing the allocations.

#### Problem AN-66: comparative statics from a unique stationary point

For $a>0$ and real $b$, minimize $f(x;a,b)=x^4/4+ax^2/2-bx$. Prove existence and uniqueness, derive the optimizer's local derivatives with respect to both parameters, and evaluate them at $a=1,b=2$.

Solution. The quartic leading term makes the objective coercive, so continuity gives existence. Its second derivative is $3x^2+a>0$, proving strict convexity and uniqueness.

The first-order equation is $x^3+ax-b=0$. Its left side is strictly increasing in $x$ and tends to opposite infinities at the two ends of the real line, so it has exactly one root. The derivative in $x$ is nonzero. The implicit-function theorem therefore gives a smooth local optimizer map, and differentiating the equation yields

$$
\frac{\partial x^*}{\partial b}=\frac1{3(x^*)^2+a}>0,\qquad
\frac{\partial x^*}{\partial a}=-\frac{x^*}{3(x^*)^2+a}.
$$

At $a=1,b=2$ the root is one, so the derivatives are $1/4$ and $-1/4$. The minimum value there is $1/4+1/2-2=-5/4$.

Let $v(a,b)=f(x^*(a,b);a,b)$. The chain-rule terms involving derivatives of $x^*$ disappear because $f_x=0$. Thus $v_a=(x^*)^2/2$ and $v_b=-x^*$. These are envelope derivatives.

The value is concave in $(a,b)$ on $a>0$ because it is the infimum over $x$ of affine functions of those parameters. This concerns parameter dependence of the optimized value. The original objective is strictly convex in the decision variable. Convexity in one set of variables and concavity in another are compatible.

#### Problem AN-67: curvature along a circle

Find the global extrema of $xy$ on $x^2+y^2=1$. Compute equality multipliers and test the Lagrangian Hessian along tangent directions.

Solution. The identities $(x+y)^2=1+2xy$ and $(x-y)^2=1-2xy$ imply $-1/2\le xy\le1/2$. The upper value occurs at $x=y=\pm1/\sqrt2$, and the lower at $x=-y=\pm1/\sqrt2$.

Use $L=xy+\nu(x^2+y^2-1)$. Stationarity gives $y+2\nu x=0$ and $x+2\nu y=0$. At a minimum with $y=-x$, $\nu=1/2$. At a maximum with $y=x$, $\nu=-1/2$.

The Hessian in the decision variables is

$$
\nabla^2_{x,y}L=
\begin{pmatrix}2\nu&1\\1&2\nu\end{pmatrix}.
$$

A unit tangent direction to the circle at $(x,y)$ is $h=(-y,x)$. At either minimum, its quadratic value is two, strictly positive. At either maximum it is minus two, strictly negative.

The Hessian need not be positive definite on the entire plane at a constrained minimum. Here its other eigenvalue is zero. Only tangent directions represent first-order feasible motion, and the multiplier incorporates the circle's curvature. The direct square identities establish global extrema; the tangent Hessian independently checks their local classification.

The constraint gradient is $2(x,y)$ and never vanishes on the circle, so the equality constraint qualification holds. Solving a multiplier system without this check could miss or mishandle a singular feasible point in another problem.

#### Problem AN-68: a polynomial whose infimum escapes to infinity

Analyze global minimization of $f(x,y)=(xy-1)^2+x^2$ over $\mathbb R^2$. Find its infimum, determine attainment, and classify its stationary point.

Solution. Both terms are nonnegative, so the infimum is at least zero. Along $x_n=1/n,y_n=n$, the product equals one and the objective equals $1/n^2$. Thus the infimum is zero.

Attaining zero would require $x=0$ from the second term and $xy=1$ from the first, an impossibility. There is no global minimizer. The domain is closed and the objective is continuous, but neither compactness nor coercivity is available.

For stationary points,

$$
f_x=2y(xy-1)+2x,\qquad f_y=2x(xy-1).
$$

The second equation implies $x=0$ or $xy=1$. The latter makes the first equation $2x=0$ and is impossible with $xy=1$. The former makes the first equation $-2y=0$. Thus the only stationary point is $(0,0)$, with value one.

Its Hessian is $\begin{pmatrix}2&-2\\-2&0\end{pmatrix}$, whose determinant is negative. It is a saddle. Direct paths confirm this: $f(t,0)=1+t^2>1$, while $f(t,t)=1-t^2+t^4<1$ for sufficiently small nonzero $t$.

A search for stationary points would therefore find a saddle and never find a minimizer, because none exists. The near-optimal sequence sends one coordinate to zero and the other to infinity. Reporting the numerical location of a small objective would need to acknowledge that this sequence does not approach a finite optimizer.
