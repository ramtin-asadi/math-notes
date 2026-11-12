# Stochastic Calculus

Brownian motion has continuous paths but does not behave like a differentiable function. Its increments are small enough to have zero mean and large enough that their squares accumulate. Stochastic calculus keeps track of that accumulated second-order contribution.

Sources:

- Steven Lalley, [Brownian Motion and Stochastic Calculus, University of Chicago](https://galton.uchicago.edu/~lalley/Courses/385/2014/index.html), including the linked lecture notes and the assignment PDFs currently dated 2016. Citations identify the actual document and problem rather than assuming its date from the older index.
- [MIT 18.642, Topics in Mathematics with Applications in Finance, Fall 2024](https://ocw.mit.edu/courses/18-642-topics-in-mathematics-with-applications-in-finance-fall-2024/), especially the stochastic-calculus and pricing lectures.
- The prerequisite treatment of [conditional expectation and Gaussian variables](03-probability-and-statistics.md).

Unless otherwise stated, $`W`$ is a standard Brownian motion relative to a filtration $`(\mathcal F_t)`$ on a probability space $`(\Omega,\mathcal F,P)`$. The filtration contains the information already observed and satisfies the usual completeness and right-continuity conditions. Equalities between stochastic processes are understood up to the appropriate almost-sure versions.

Contents:

1. [Information and martingales](#1-information-and-martingales)
2. [Brownian motion](#2-brownian-motion)
3. [Stochastic integration](#3-stochastic-integration)
4. [Itô's formula](#4-itôs-formula)
5. [Stochastic differential equations](#5-stochastic-differential-equations)
6. [Risk-neutral valuation](#6-risk-neutral-valuation)
7. [Pricing PDEs and Feynman-Kac](#7-pricing-pdes-and-feynman-kac)
8. [Further worked problems](#8-further-worked-problems)

## 1. Information and martingales

### 1.1 Theory

A filtration is an increasing family of sigma-algebras: $`\mathcal F_s\subseteq\mathcal F_t`$ for $`s\le t`$. It records which events can be determined by each date. A process $`X_t`$ is adapted if $`X_t`$ is $`\mathcal F_t`$-measurable. A decision made at date $`t`$ must use information available then, rather than a future observation.

An integrable adapted process $`M_t`$ is a martingale if $`E[M_t\mid\mathcal F_s]=M_s`$ for every $`s\le t`$. Integrability, $`E\vert M_t\vert \lt \infty`$, is part of the definition. Constant unconditional means alone are insufficient: the conditional requirement must hold after any admissible past information.

A stopping time $`\tau`$ satisfies $`\lbrace \tau\le t\rbrace \in\mathcal F_t`$ for every $`t`$. A first hitting time of a closed level by a continuous adapted process is a standard example. The time at which a path attains its maximum over a future fixed interval is generally not a stopping time, since later observations can change whether the current maximum is final.

Optional stopping is an assumption-dependent extension of the martingale mean property. For a right-continuous martingale and a bounded stopping time, $`E[M_\tau]=E[M_0]`$. For unbounded stopping times, one needs conditions such as uniform integrability of the stopped family, or other hypotheses that justify passing to the limit.

### 1.2 Derivations and proofs

#### A conditional-expectation martingale

For an integrable terminal variable $`Z`$, define $`M_t=E[Z\mid\mathcal F_t]`$. It is adapted and integrable. The tower property gives $`E[M_t\mid\mathcal F_s]=E[E[Z\mid\mathcal F_t]\mid\mathcal F_s]=E[Z\mid\mathcal F_s]=M_s`$. This construction produces martingales from future payoffs without claiming that the future payoff itself is known early.

#### Optional stopping for a bounded discrete time

Let $`\tau\le N`$ be an integer-valued stopping time. Write

```math
M_\tau=M_0+\sum_{k=1}^N
\mathbf1_{\lbrace \tau\ge k\rbrace }(M_k-M_{k-1}).
```

The event $`\lbrace \tau\ge k\rbrace =\lbrace \tau\gt k-1\rbrace`$ belongs to $`\mathcal F_{k-1}`$. Condition each summand on that sigma-algebra: its expectation is zero because the indicator is already known and the martingale increment has conditional mean zero. Summing finitely many integrable terms proves $`E[M_\tau]=E[M_0]`$.

#### A finite-partition conditional expectation

If information distinguishes only atoms $`A_1,\ldots,A_k`$, then the conditional mean on each positive-probability atom is the average of $`X`$ restricted to it. This is a concrete model of information: within one atom no further outcome distinctions are observable.

#### What information a martingale statement uses

A martingale is always defined relative to a probability measure and a filtration. The statement $`E[M_t\mid\mathcal F_s]=M_s`$ requires that $`M_t`$ be integrable and that $`M_s`$ be measurable with respect to the information available at date $`s`$. It does not say that a sample path is constant, increasing, or equally likely to move up and down.

For a discrete process, the one-step conditions imply the multi-step conditions by repeated use of the tower property. If $`E[M_{n+1}\mid\mathcal F_n]=M_n`$, then

```math
E[M_{n+2}\mid\mathcal F_n]
=E[E[M_{n+2}\mid\mathcal F_{n+1}]\mid\mathcal F_n]
=E[M_{n+1}\mid\mathcal F_n]=M_n.
```

Iterating proves the assertion for all later dates. An adapted integrable process with constant unconditional expectation need not be a martingale: its conditional expectation can still depend on current information in the wrong way.

A stopping time allows the decision to stop at date $`n`$ to be made from $`\mathcal F_n`$. The event $`\lbrace \tau\le n\rbrace`$ must therefore be known at that date. A first hitting time is a standard example. A last visit before a future deadline generally is not, because verifying that a visit was the last one requires observing the remaining path.

#### A complete stopped-walk calculation

Let $`S_n`$ be a simple symmetric random walk started at integer $`i`$ with $`0\lt i\lt m`$. Let $`\tau`$ be its first visit to zero or $`m`$. We derive both the probability of reaching the upper boundary and the expected duration.

The processes $`S_n`$ and $`S_n^2-n`$ are martingales. For the bounded stopping time $`\tau\wedge N`$, optional stopping gives

```math
E[S_{\tau\wedge N}]=i,\qquad
E[\tau\wedge N]=E[S_{\tau\wedge N}^2]-i^2\le m^2.
```

The stopped position always belongs to $`[0,m]`$. Monotone convergence applied to $`\tau\wedge N`$ gives $`E[\tau]\le m^2\lt \infty`$, which in particular proves $`\tau\lt \infty`$ almost surely. We can now pass to the terminal position. Dominated convergence applies to $`S_{\tau\wedge N}`$ and its square because they are bounded by $`m`$ and $`m^2`$.

Writing $`p=P(S_\tau=m)`$, the first identity becomes $`mp=i`$. The second becomes

```math
E[\tau]=m^2p-i^2=i(m-i).
```

For $`i=2`$ and $`m=5`$, the upper-boundary probability is $`2/5`$ and the expected duration is six steps. The proof did not apply optional stopping directly at an unbounded time without justification. It first stopped at a finite truncation and then checked the limit.

#### Example SC-E21: a last visit uses future information

For a symmetric random walk with $`S_0=0`$, define $`L=\max\lbrace n\in\lbrace 0,1,2\rbrace :S_n=0\rbrace`$. Is $`L`$ a stopping time?

Solution. The walk cannot be zero at date one. It returns at date two exactly when the two increments have opposite signs. Thus $`L=2`$ on that event and $`L=0`$ otherwise. The event $`\lbrace L\le0\rbrace =\lbrace S_2\ne0\rbrace`$ has probability one half and is not measurable with respect to the initial, trivial information. Hence $`L`$ is not a stopping time. Knowing the current position at zero does not tell us whether another zero will occur later.

### 1.3 Worked examples

#### Example SC-E01: a random-walk martingale

Let independent increments $`\xi_i`$ equal $`\pm1`$ with equal probabilities and $`S_n=\sum_{i=1}^n\xi_i`$. With the history filtration, $`E[S_{n+1}\mid\mathcal F_n]=S_n+E[\xi_{n+1}\mid\mathcal F_n]=S_n`$. The process is integrable because $`\vert S_n\vert \le n`$. Also $`S_n^2-n`$ is a martingale, since the conditional expansion of $`(S_n+\xi_{n+1})^2`$ adds exactly one.

#### Example SC-E02: bounded stopped gambling

For the same walk, stop at the first hit of $`+2`$ or $`-1`$, or at date ten if neither level is reached. The stopping date is bounded, so $`ES_\tau=0`$. If the date-ten cap were removed, the conclusion would need a further argument. It happens to remain valid for exit from a finite interval because the stopped walk is bounded, but that justification must be supplied.

### 1.4 Selected exercises

#### Problem SC-01: conditioning on a finite information partition

Source: [Lalley, Conditional Expectation, Exercise 2](https://galton.uchicago.edu/~lalley/Courses/385/ConditionalExpectation.pdf).

Problem. For a finite measurable partition and an integrable variable $`X`$, derive the conditional expectation on the partition's sigma-algebra and check the $`L^2`$ bound when $`X`$ is square integrable.

Solution. For atoms $`A_i`$ with positive probabilities define $`m=\sum_i[E(X\mathbf1_{A_i})/P(A_i)]\mathbf1_{A_i}`$. Set arbitrary values on probability-zero atoms. The function is measurable with respect to the partition. For an event $`B`$ in its sigma-algebra, $`B`$ is a union of atoms, and

```math
E[m\mathbf1_B]
=\sum_{A_i\subseteq B}E[X\mathbf1_{A_i}]
=E[X\mathbf1_B].
```

Thus it satisfies the defining conditional-expectation identity. For $`X\in L^2`$, Jensen on each atom gives $`E[m^2]\le E[X^2]`$, so it belongs to the required square-integrable space as well.

#### Problem SC-02: minima and maxima of stopping times

Source: [Lalley, Brownian Motion, Exercise 3(a)](https://galton.uchicago.edu/~lalley/Courses/385/BrownianMotion.pdf).

Problem. Prove that the minimum and maximum of two stopping times are stopping times.

Solution. For stopping times $`\tau,\sigma`$,

```math
\lbrace \tau\wedge\sigma\le t\rbrace =\lbrace \tau\le t\rbrace \cup\lbrace \sigma\le t\rbrace ,
\quad
\lbrace \tau\vee\sigma\le t\rbrace =\lbrace \tau\le t\rbrace \cap\lbrace \sigma\le t\rbrace .
```

Both events belong to $`\mathcal F_t`$ by closure of a sigma-algebra under unions and intersections. Hence both random times are stopping times. This proof directly checks the information requirement rather than relying on the names of the operations.

#### Problem SC-03: rounding a stopping time upward

Source: [Lalley, Brownian Motion, Exercise 3(b)](https://galton.uchicago.edu/~lalley/Courses/385/BrownianMotion.pdf).

Problem. For a finite stopping time $`\tau`$, show that its upward dyadic rounding $`2^{-n}\lceil2^n\tau\rceil`$ is a stopping time and decreases to $`\tau`$.

Solution. Let $`\tau_n=2^{-n}\lceil2^n\tau\rceil`$. For fixed $`t`$, put $`q_n(t)=2^{-n}\lfloor2^nt\rfloor`$. Then $`\lbrace \tau_n\le t\rbrace =\lbrace \tau\le q_n(t)\rbrace \in\mathcal F_{q_n(t)}\subseteq\mathcal F_t`$. Thus $`\tau_n`$ is a stopping time. Also $`\tau\le\tau_n\lt \tau+2^{-n}`$, so $`\tau_n\downarrow\tau`$ for finite $`\tau`$. Upward rounding permits a discrete-grid approximation using only information available at the rounded date.

#### Problem SC-04: the stopping sigma-algebra

Source: [Lalley, Brownian Motion, Exercise 4(a)](https://galton.uchicago.edu/~lalley/Courses/385/BrownianMotion.pdf).

Problem. Prove that the events observable at a stopping time form a sigma-algebra under the displayed definition of $`\mathcal F_\tau`$.

Solution. Define $`\mathcal F_\tau=\lbrace A\in\mathcal F:A\cap\lbrace \tau\le t\rbrace \in\mathcal F_t\text{ for every }t\rbrace`$. The whole space belongs because $`\tau`$ is a stopping time. If $`A`$ belongs, then $`A^c\cap\lbrace \tau\le t\rbrace =\lbrace \tau\le t\rbrace \setminus[A\cap\lbrace \tau\le t\rbrace ]`$ belongs to $`\mathcal F_t`$. For countably many $`A_i`$, intersecting their union with $`\lbrace \tau\le t\rbrace`$ gives the union of the corresponding measurable intersections. Thus $`\mathcal F_\tau`$ is a sigma-algebra.

#### Problem SC-05: earlier stopping reveals less information

Source: [Lalley, Brownian Motion, Exercise 4(b)](https://galton.uchicago.edu/~lalley/Courses/385/BrownianMotion.pdf).

Problem. If stopping times satisfy $`\sigma\le\tau`$, prove $`\mathcal F_\sigma\subseteq\mathcal F_\tau`$.

Solution. If $`\sigma\le\tau`$ and $`A\in\mathcal F_\sigma`$, then

```math
A\cap\lbrace \tau\le t\rbrace
=[A\cap\lbrace \sigma\le t\rbrace ]\cap\lbrace \tau\le t\rbrace .
```

Both factors belong to $`\mathcal F_t`$, so $`A\in\mathcal F_\tau`$. Therefore $`\mathcal F_\sigma\subseteq\mathcal F_\tau`$. The inequality of random times, holding pointwise or after a suitable almost-sure version, is what permits the event identity.

### 1.5 Additional solved problems

#### Problem SC-06: a stopped walk with no integrable terminal passage time

Constructed problem. Explain the failure of naively stopping a fair random walk at its first hit of $`+1`$.

Solution. The hitting time is almost surely finite, and the terminal walk value is one. If one applied optional stopping without further assumptions, it would predict mean zero and contradict that value. The bounded approximations $`\tau\wedge n`$ do have mean-zero stopped values, but they are not uniformly integrable in the way needed to pass their expectations to the limiting terminal value. Almost-sure finiteness alone is insufficient.

#### Problem SC-07: deterministic drift

Constructed problem. Is $`W_t+t`$ a martingale?

Solution. For $`s\lt t`$, $`E[W_t+t\mid\mathcal F_s]=W_s+t`$, whereas the process at $`s`$ is $`W_s+s`$. The difference is $`t-s`$. It is an integrable adapted submartingale, not a martingale. Subtracting its known drift returns $`W_t`$.

## 2. Brownian motion

### 2.1 Theory

A standard Brownian motion starts at zero, has continuous paths, and has independent Gaussian increments: $`W_t-W_s\sim N(0,t-s)`$ for $`s\lt t`$, independently of $`\mathcal F_s`$. Its increment standard deviation is $`\sqrt{t-s}`$, which is larger than the time length at small scales.

The covariance is $`E[W_sW_t]=\min(s,t)`$. To derive this for $`s\le t`$, write $`W_t=W_s+(W_t-W_s)`$ and use independence and zero mean of the increment. Finite-dimensional distributions are multivariate Gaussian because their coordinates are sums of independent Gaussian increments.

For $`c\gt 0`$, $`c^{-1/2}W_{ct}`$ is another standard Brownian motion in law. Its increments have variance $`c^{-1}c(t-s)=t-s`$, and independence and continuity are preserved. This is a statement about process distributions, not equality of individual paths.

Quadratic variation measures accumulated squared increments. Along deterministic partitions whose mesh tends to zero, Brownian squared increments sum to elapsed time. For smooth finite-variation paths the same sums tend to zero. This difference creates the extra term in Itô's formula.

### 2.2 Derivations and proofs

#### Quadratic variation in $`L^2`$

For partition $`0=t_0\lt \cdots\lt t_n=T`$, let $`Q_\pi=\sum_i(\Delta W_i)^2`$. A centered normal increment of variance $`h_i`$ has squared mean $`h_i`$ and fourth moment $`3h_i^2`$, so its square has variance $`2h_i^2`$. Independence gives

```math
EQ_\pi=T,\qquad
E[(Q_\pi-T)^2]=2\sum_i h_i^2
\le2T\vert \pi\vert \to0.
```

Thus $`Q_\pi\to T`$ in $`L^2`$ and probability. For dyadic equal partitions, the variance bounds are summable over refinement levels. Chebyshev and Borel-Cantelli give almost-sure convergence along that sequence too.

#### Infinite variation

For a continuous finite-variation path $`f`$,

```math
\sum_i(\Delta f_i)^2
\le\max_i\vert \Delta f_i\vert \sum_i\vert \Delta f_i\vert
\le\max_i\vert \Delta f_i\vert \,\mathop{\mathrm{TV}} _{[0,T]}(f)\to0,
```

by uniform continuity. Brownian motion has dyadic squared-increment sums tending almost surely to $`T\gt 0`$, so its total variation cannot be finite on that interval. This rules out ordinary pathwise Riemann-Stieltjes integration against Brownian motion under the usual finite-variation integrator framework.

#### Overlapping increments and Brownian scaling

For two Brownian increments, the covariance equals the length of overlap of their time intervals. For example,

```math
\mathop{\mathrm{Cov}} (W_3-W_1,W_4-W_2)=1,
```

because the shared increment over $`[2,3]`$ has variance one and the remaining disjoint pieces are independent. More generally, write each increment as an integral of an interval indicator against Brownian motion; their covariance is the integral of the product of the indicators.

The familiar covariance $`E[W_sW_t]=\min(s,t)`$ follows by taking intervals starting at zero. Every finite vector of Brownian values is jointly Gaussian. Its mean and covariance determine its distribution, so a linear combination $`\sum_i a_iW_{t_i}`$ is centered normal with variance $`\sum_{i,j}a_ia_j\min(t_i,t_j)`$.

For $`c\gt 0`$, the process $`B_t=W_{ct}/\sqrt c`$ has continuous paths, starts at zero, and has independent centered Gaussian increments with variance equal to their time length. It is therefore Brownian motion. This is equality in process law, not an assertion that the two processes have identical paths on the original sample space.

Brownian increments have a typical size proportional to the square root of the elapsed time. Dividing an increment by its time length gives a random slope whose variance grows as the interval shrinks. A full proof that Brownian paths are almost surely nowhere differentiable requires more than this variance observation; it is a path theorem in the cited notes. The quadratic-variation argument above already proves infinite total variation, which is enough to rule out ordinary finite-variation integration against Brownian paths.

#### Reflection and the maximum

For $`a\gt 0`$, reflect a Brownian path after its first hit of $`a`$. Symmetry and the strong Markov property imply that the reflected process has the same Brownian law. The strong Markov property, stated here as a theorem, extends independent increments to a finite stopping time: the subsequent increments form a fresh Brownian motion independent of the stopped information.

Among paths that hit $`a`$ by time $`T`$, reflection pairs terminal values below $`a`$ with terminal values above it. Consequently

```math
P\left(\max_{0\le t\le T}W_t\ge a\right)
=2P(W_T\ge a)
=2\left[1-\Phi\left(\frac a{\sqrt T}\right)\right].
```

The first-passage time $`\tau_a`$ therefore has CDF equal to this expression as a function of $`T\gt 0`$. Differentiating gives density

```math
f_{\tau_a}(t)=\frac a{\sqrt{2\pi t^3}}e^{-a^2/(2t)},\qquad t\gt 0.
```

It hits eventually with probability one, but the tail is of order $`t^{-1/2}`$ and its mean is infinite. Almost-sure finiteness and finite expected duration are different properties.

#### Example SC-E22: a conditional barrier crossing

Let $`a\gt 0`$ and condition on $`W_T=b\lt a`$. Find the probability that the path crossed $`a`$ before $`T`$.

Solution. Reflection gives joint terminal density $`\phi_T(2a-b)`$ for paths that crossed the barrier and ended at $`b`$, where $`\phi_T(x)=(2\pi T)^{-1/2}e^{-x^2/(2T)}`$. Dividing by the unconditional terminal density gives

```math
P\left(\max_{t\le T}W_t\ge a\mid W_T=b\right)
=\frac{\phi_T(2a-b)}{\phi_T(b)}
=\exp\left(-\frac{2a(a-b)}T\right).
```

For $`a=2`$, $`b=1`$, $`T=1`$, this is $`e^{-4}`$. If $`b\ge a`$, continuity makes the crossing probability one. Checking only simulated endpoints would miss some barrier crossings when both endpoints lie below the barrier; the conditional formula quantifies that missing path information.

### 2.3 Worked examples

#### Example SC-E03: a conditional increment

Given $`\mathcal F_2`$, $`W_5-W_2\sim N(0,3)`$ and is independent of the past. Hence $`E[W_5\mid\mathcal F_2]=W_2`$ and $`\mathop{\mathrm{Var}} (W_5\mid\mathcal F_2)=3`$. Conditioning only on $`W_2=1`$ gives $`W_5\sim N(1,3)`$ under its conditional law.

#### Example SC-E04: a correlated pair of dates

The covariance of $`(W_1,W_3)`$ is

```math
\left(\begin{smallmatrix}1&1\\[0pt] 1&3\end{smallmatrix}\right)
```

. The correlation is $`1/\sqrt3`$. The variables are not independent, although $`W_1`$ and $`W_3-W_1`$ are independent. Independent increments should not be confused with independent levels.

#### Example SC-E05: cross-variation of independent coordinates

For independent Brownian motions $`W,B`$, the sum $`\sum_i\Delta W_i\Delta B_i`$ has mean zero and variance $`\sum_i h_i^2\to0`$. Thus their cross-variation is zero. For $`B=\rho W+\sqrt{1-\rho^2}Z`$ with independent Brownian $`Z`$, the cross-variation is instead $`\rho T`$.

### 2.4 Selected exercise and additional problems

#### Problem SC-08: a complex exponential martingale

Source: [Lalley, Homework 2, Problem 1(a)](https://galton.uchicago.edu/~lalley/Courses/385/HW2.pdf).

Problem. For independent Brownian coordinates $`X,Y`$ and real $`\theta`$, prove that $`\exp(i\theta X_t+\theta Y_t)`$ is a complex-valued martingale.

Solution. Let $`X,Y`$ be independent standard Brownian coordinates and $`M_t=\exp(i\theta X_t+\theta Y_t)`$ for real $`\theta`$. For $`s\lt t`$, factor out $`M_s`$ and use independent increments:

```math
E[M_t\mid\mathcal F_s]
=M_s e^{-\theta^2(t-s)/2}e^{\theta^2(t-s)/2}=M_s.
```

Integrability holds because $`\vert M_t\vert =e^{\theta Y_t}`$ has finite expectation $`e^{\theta^2t/2}`$. Complex martingale statements mean that both real and imaginary parts obey the conditional mean property.

#### Problem SC-09: a Brownian bridge at one date

Constructed problem. Find the conditional law of $`W_s`$ given $`W_T=w`$, for $`0\lt s\lt T`$.

Solution. The covariance of $`(W_s,W_T)`$ is

```math
\left(\begin{smallmatrix}s&s\\[0pt] s&T\end{smallmatrix}\right)
```

. Conditional Gaussian formulas give mean $`(s/T)w`$ and variance $`s-s^2/T=s(T-s)/T`$. Knowing the future endpoint reduces variance relative to the unconditioned value $`s`$. This conditional model is a bridge; it does not make the original forward Brownian increments independent after the endpoint is revealed.

#### Problem SC-10: a squared Brownian martingale

Constructed problem. Prove $`W_t^2-t`$ is a martingale.

Solution. Let $`\Delta=W_t-W_s`$, independent of $`\mathcal F_s`$, with mean zero and variance $`t-s`$. Expanding,

```math
E[W_t^2-t\mid\mathcal F_s]
=W_s^2+2W_sE\Delta+E\Delta^2-t=W_s^2-s.
```

Second moments are finite, so integrability is satisfied. The subtracted time compensates exactly for variance accumulated after $`s`$.

#### Problem SC-11: distributional scaling of a first passage time

Constructed problem. Let $`\tau_a=\inf\lbrace t:W_t=a\rbrace`$ for $`a\gt 0`$. Determine its scale relative to $`\tau_1`$.

Solution. The process $`a^{-1}W_{a^2t}`$ is standard Brownian. It reaches one exactly when the original path reaches $`a`$ at time $`a^2t`$. Hence $`\tau_a\overset d=a^2\tau_1`$. Scaling a distance by $`a`$ scales the hitting time by $`a^2`$. No finite mean is implied: the one-sided Brownian first-passage time has infinite expectation.

## 3. Stochastic integration

### 3.1 Theory

Start with a simple predictable process $`H_t=\sum_iH_i\mathbf1_{(t_i,t_{i+1}]}(t)`$, where $`H_i`$ is $`\mathcal F_{t_i}`$-measurable. Define

```math
\int_0^T H_t\,dW_t=\sum_iH_i(W_{t_{i+1}}-W_{t_i}).
```

The left-endpoint information condition prevents choosing a stake after seeing the increment it multiplies. If $`E\int_0^TH_t^2dt\lt \infty`$, the integral extends by $`L^2`$ limits to predictable square-integrable integrands. Continuous adapted integrands are predictable in the usual setup.

The resulting integral process is a square-integrable martingale, starts at zero, and has expectation zero. If only $`\int_0^TH_t^2dt\lt \infty`$ almost surely without its expectation being finite, one can define a local integral by localization, but its local-martingale status does not automatically supply zero expectation.

### 3.2 Derivations and proofs

#### Itô isometry

For a simple integrand, expand the square of its integral. Each diagonal term has expectation $`E[H_i^2](t_{i+1}-t_i)`$ because the coefficient is known before the independent increment. For $`i\lt j`$, condition the cross term on $`\mathcal F_{t_j}`$. All earlier factors and $`H_j`$ are then known, while the last increment has mean zero. Thus cross terms vanish:

```math
E\left[\left(\int_0^TH_t\,dW_t\right)^2\right]
=E\int_0^TH_t^2dt.
```

The same calculation for the difference of two simple integrands shows that approximation in the time-probability $`L^2`$ norm produces a Cauchy sequence of integrals. Completeness defines the limit and makes it independent of the approximating sequence. That is the construction's role for the isometry.

#### Integrating Brownian motion itself

For a partition,

```math
W_T^2=\sum_i[W_{t_{i+1}}^2-W_{t_i}^2]
=2\sum_iW_{t_i}\Delta W_i+\sum_i(\Delta W_i)^2.
```

The first sum converges to twice the Itô integral and the second to $`T`$. Therefore

```math
\int_0^T W_t\,dW_t=\frac12(W_T^2-T).
```

Ordinary differentiation of $`W_t^2/2`$ would miss the accumulated squared increments. The correction has a nonzero deterministic size even as each increment becomes small.

#### Why the left endpoint appears in the definition

For an elementary integrand $`H_t=\sum_i H_i1_{(t_i,t_{i+1}]}(t)`$, the coefficient $`H_i`$ is known at $`t_i`$. The stochastic integral is $`\sum_i H_i(W_{t_{i+1}}-W_{t_i})`$. This matches a strategy that chooses its holding before the price increment occurs.

Conditioning on $`\mathcal F_{t_i}`$ gives expectation zero for each term, because the coefficient is known while the future Brownian increment is independent and centered. For two different intervals, the cross product also has expectation zero: condition at the start of the later interval, when all the earlier factors and the later coefficient are already known.

The diagonal terms satisfy

```math
E[H_i^2(W_{t_{i+1}}-W_{t_i})^2]
=(t_{i+1}-t_i)E[H_i^2].
```

Summing proves the isometry. Completion in the norm $`(E\int_0^T H_t^2dt)^{1/2}`$ then defines the integral for predictable square-integrable integrands. Continuous adapted processes with this integrability are an important included class. Measurability in time and information is part of the construction; saying only that each separate $`H_t`$ is adapted does not specify every technical detail for an arbitrary irregular process.

The integral process is a martingale. For elementary integrands, future increments have conditional mean zero by the same conditioning argument. For general square-integrable integrands, approximate by elementary ones and use the contraction property of conditional expectation in $`L^2`$ to pass the martingale identity to the limit.

#### Left sums, right sums, and symmetric sums

For the integrand $`W_t`$, the right-endpoint sum differs from the left-endpoint sum by

```math
\sum_i W_{t_{i+1}}\Delta W_i-\sum_i W_{t_i}\Delta W_i
=\sum_i(\Delta W_i)^2.
```

The difference converges to $`T`$, not zero. Thus the choice of evaluation point survives refinement. The left limit is the Itô integral $`(W_T^2-T)/2`$. The symmetric endpoint-average sum is the average of left and right sums and converges to $`W_T^2/2`$.

The symmetric limit is the Stratonovich integral in this example. For a smooth state function $`f(W_t)`$, its relation to the Itô integral includes $`\tfrac12\int f'(W_t)dt`$. This is another convention with its own calculus. It does not make a right-endpoint coefficient an admissible left-endpoint Itô trading strategy.

#### Example SC-E23: an integral that is a product of normals

Let $`I=\int_1^3 W_1\,dW_t`$. Find its conditional law given $`\mathcal F_1`$, its variance, and whether it is Gaussian unconditionally.

Solution. The coefficient is already known when the interval begins, so $`I=W_1(W_3-W_1)`$. Conditional on $`\mathcal F_1`$, it is normal with mean zero and variance $`2W_1^2`$. Unconditionally it is a mixture with a random variance.

Independence of $`W_1`$ and $`W_3-W_1`$ gives $`E[I^2]=1\cdot2=2`$. Their fourth moments are three and twelve, so $`E[I^4]=36`$. A centered normal with variance two would instead have fourth moment twelve. Hence $`I`$ is not Gaussian. The isometry gives its variance correctly without implying a Gaussian distribution for a random integrand.

### 3.3 Worked examples

#### Example SC-E06: a deterministic integrand

For $`H_t=t`$, the integral is centered Gaussian with variance $`\int_0^Tt^2dt=T^3/3`$. Gaussianity follows by approximating with deterministic linear combinations of Brownian increments and taking their distributional limit. Its mean is zero; the ordinary integral $`\int_0^Tt\,dt=T^2/2`$ is a different operation.

#### Example SC-E07: a coefficient learned at date one

For $`H_t=W_1\mathbf1_{(1,2]}(t)`$, the integral over $`[0,2]`$ is $`W_1(W_2-W_1)`$. Its mean is zero, and variance is $`E[W_1^2]E[(W_2-W_1)^2]=1`$. The coefficient is random but known when its integration interval begins.

#### Example SC-E08: two deterministic integrals

For $`I=\int_0^T f(t)dW_t`$ and $`J=\int_0^Tg(t)dW_t`$, polarization of the isometry gives $`\mathop{\mathrm{Cov}} (I,J)=\int_0^Tf(t)g(t)dt`$. With $`f=1`$ and $`g=t`$, the covariance is $`T^2/2`$. The integrals share the same Brownian driver and are generally correlated.

### 3.4 Additional solved problems

#### Problem SC-12: a random integrand and a variance

Constructed problem. Find the variance of $`\int_0^TW_t\,dW_t`$.

Solution. The isometry gives $`\int_0^TE[W_t^2]dt=\int_0^Tt\,dt=T^2/2`$. The explicit formula gives the same answer: since $`W_T\sim N(0,T)`$, $`\mathop{\mathrm{Var}} (W_T^2)=2T^2`$, and dividing by four gives $`T^2/2`$. These two routes check both the stochastic identity and the moment calculation.

#### Problem SC-13: a deterministic exponential kernel

Constructed problem. Find the law of $`I_T=\int_0^Te^{-\kappa(T-s)}dW_s`$ for $`\kappa\gt 0`$.

Solution. It is a deterministic-integrand Gaussian integral with mean zero. Its variance is $`\int_0^Te^{-2\kappa(T-s)}ds=(1-e^{-2\kappa T})/(2\kappa)`$. As $`T\to\infty`$, the variance approaches $`1/(2\kappa)`$. This kernel is the random component of the Ornstein-Uhlenbeck solution.

#### Problem SC-14: the information condition cannot be skipped

Constructed problem. On a one-step interval $`[0,T]`$, compare the adapted constant coefficient $`H=1`$ with the future-dependent coefficient $`H=W_T`$.

Solution. The first gives integral $`W_T`$ with mean zero. Naively multiplying the second coefficient by the increment gives $`W_T^2`$, with expectation $`T`$. That second expression does not belong to the simple predictable definition because $`W_T`$ is unavailable at time zero. Applying the mean-zero Itô-integral theorem to it would violate the theorem's assumptions.

## 4. Itô's formula

### 4.1 Theory

Let $`dX_t=b_tdt+\sigma_tdW_t`$ and let $`f(t,x)`$ have one continuous time derivative and two continuous space derivatives. Under local integrability and the usual localization conditions,

```math
df(t,X_t)=
\left(f_t+b_tf_x+\frac12\sigma_t^2f_{xx}\right)dt
+\sigma_tf_xdW_t.
```

Every derivative is evaluated at $`(t,X_t)`$. The notation is an integral identity: integrating the right-hand side gives $`f(T,X_T)-f(0,X_0)`$. It does not define an ordinary derivative $`dW_t/dt`$.

The mnemonic $`(dW)^2=dt`$, $`dt\,dW=0`$, and $`(dt)^2=0`$ summarizes orders in this limiting calculation. These are quadratic-variation rules, rather than algebraic equalities between ordinary infinitesimals.

For two Itô processes driven by correlated Brownian motions, the product rule includes their covariation:

```math
d(XY)=X\,dY+Y\,dX+d[X,Y].
```

If both use the same Brownian motion with diffusion coefficients $`\sigma_X,\sigma_Y`$, then $`d[X,Y]=\sigma_X\sigma_Ydt`$.

### 4.2 Derivations and proofs

#### The second-order contribution

For Brownian input and a smooth $`f`$, Taylor expansion of each small spatial increment contains $`f_x\Delta W+\tfrac12f_{xx}(\Delta W)^2`$. The first sum converges to the stochastic integral. For bounded adapted coefficients, the difference between the weighted square sum and its time-step counterpart has variance bounded by a constant times $`\sum_i(\Delta t_i)^2`$, so it vanishes as the mesh shrinks.

Third-order spatial remainders vanish under bounded third derivatives, since $`E\sum_i\vert \Delta W_i\vert ^3`$ is of order $`T\sqrt{\vert \pi\vert }`$. Mixed time-space terms and squared time increments also vanish. Time variation contributes $`\int f_tdt`$. Localizing to bounded state regions removes global derivative bounds, and standard approximation extends the argument to the $`C^{1,2}`$ regularity used in the theorem.

Thus the Taylor second-order term survives because Brownian quadratic variation is nonzero. For a differentiable input path, its square-increment sum would vanish and the ordinary chain rule would result.

#### Deriving the product rule from Itô

For $`f(x,y)=xy`$, the pure second derivatives vanish while both mixed derivatives equal one. Multidimensional Itô contributes $`\tfrac12(f_{xy}+f_{yx})d[X,Y]=d[X,Y]`$, in addition to $`Y\,dX+X\,dY`$. The factor one half and two equal mixed derivatives explain why the cross term appears once.

#### The differential multiplication rules

For an Itô process $`dX=b\,dt+\sigma\,dW`$, the rules $`(dW)^2=dt`$, $`dt\,dW=0`$, and $`(dt)^2=0`$ summarize limits of sums, not ordinary multiplication of infinitesimal real numbers. A squared Brownian increment contributes order $`\Delta t`$ on average, and these contributions accumulate. Products involving one time increment and one Brownian increment vanish in the limiting sums under the usual regularity conditions.

With two Brownian motions satisfying $`d[W^1,W^2]_t=\rho\,dt`$, the cross rule is $`dW^1dW^2=\rho\,dt`$. For processes with diffusions $`\sigma_X`$ and $`\sigma_Y`$, their cross variation contributes $`\rho\sigma_X\sigma_Ydt`$.

The two-dimensional formula contains

```math
\frac12f_{xx}\,d[X]+\frac12f_{yy}\,d[Y]+f_{xy}\,d[X,Y].
```

The mixed term has no factor one half because the two symmetric mixed derivatives combine. Taking $`f(x,y)=xy`$ gives the product rule. Taking a time-dependent function also adds $`f_tdt`$.

These rules should always be used after identifying the drift, diffusion, covariance structure, and domain of the function. For $`\log X`$ or $`1/X`$, positivity or a stopping argument away from zero is necessary. Applying a formula with singular derivatives at a point the process may reach requires further analysis.

#### Example SC-E24: the ratio of correlated geometric processes

Let positive processes satisfy

```math
\frac{dX_t}{X_t}=\mu_Xdt+\sigma_XdW_t^1,\qquad
\frac{dY_t}{Y_t}=\mu_Ydt+\sigma_YdW_t^2,
\qquad d[W^1,W^2]_t=\rho\,dt.
```

Find the dynamics of $`R_t=X_t/Y_t`$.

Solution. For $`f(x,y)=x/y`$, the derivatives are $`f_x=1/y`$, $`f_y=-x/y^2`$, $`f_{xx}=0`$, $`f_{xy}=-1/y^2`$, and $`f_{yy}=2x/y^3`$. Substitution gives

```math
\frac{dR_t}{R_t}
=\left(\mu_X-\mu_Y+\sigma_Y^2-\rho\sigma_X\sigma_Y\right)dt
+\sigma_XdW_t^1-\sigma_YdW_t^2.
```

The instantaneous variance rate of the noise is $`\sigma_X^2+\sigma_Y^2-2\rho\sigma_X\sigma_Y`$. It is nonnegative for $`\vert \rho\vert \le1`$.

There is an independent check. Subtract the log equations for $`X`$ and $`Y`$:

```math
d\log R_t
=\left(\mu_X-\mu_Y-\tfrac12\sigma_X^2+\tfrac12\sigma_Y^2\right)dt
+\sigma_XdW_t^1-\sigma_YdW_t^2.
```

Applying Itô's exponential formula to this log ratio adds half its variance rate to the level drift, reproducing the term $`\sigma_Y^2-\rho\sigma_X\sigma_Y`$. Omitting either the reciprocal correction or the cross variation would fail this check.

#### Time-dependent martingale checks

For $`f(t,W_t)`$, the drift is $`f_t+\tfrac12f_{xx}`$. A smooth solution of the backward heat equation therefore gives a local martingale after composition with Brownian motion. To claim a true martingale, check integrability of the stochastic integral or another sufficient condition. Polynomials on a finite horizon have the needed finite Gaussian moments. A merely vanishing formal drift is not a general proof of the true martingale property.

### 4.3 Worked examples

#### Example SC-E09: a cubic integral

Apply Itô to $`f(x)=x^3`$:

```math
d(W_t^3)=3W_t^2dW_t+3W_tdt.
```

Therefore $`\int_0^TW_t^2dW_t=W_T^3/3-\int_0^TW_tdt`$. Both terms are random. Subtracting a deterministic multiple of $`T`$ would not be the correct correction for the cubic.

#### Example SC-E10: an exponential

For $`M_t=e^{\theta W_t-\theta^2t/2}`$, the time derivative is $`-\theta^2M_t/2`$ and the second spatial derivative is $`\theta^2M_t`$. Their drift contributions cancel, leaving $`dM_t=\theta M_tdW_t`$. Independent Gaussian increments verify that it is a true martingale on each finite horizon, beyond the local differential calculation.

#### Example SC-E11: a product with a common shock

If $`dX=a\,dt+\sigma\,dW`$ and $`dY=b\,dt+\eta\,dW`$, then

```math
d(XY)=(aY+bX+\sigma\eta)dt+(\sigma Y+\eta X)dW.
```

The extra $`\sigma\eta`$ drift measures the shared instantaneous noise. With independent drivers the cross term vanishes; with correlation $`\rho`$ it becomes $`\rho\sigma\eta`$.

### 4.4 Selected exercises

#### Problem SC-15: normal even moments from Itô

Source: [Lalley, Homework 4, Problem 1(A)](https://galton.uchicago.edu/~lalley/Courses/385/HW4.pdf).

Problem. Use Itô's formula and Brownian scaling to derive the even-moment recursion for a standard normal variable and evaluate all even moments.

Solution. Apply Itô to $`W_t^{2n}`$:

```math
dW_t^{2n}=2nW_t^{2n-1}dW_t+n(2n-1)W_t^{2n-2}dt.
```

Gaussian moments ensure the finite-horizon stochastic integral is integrable with zero mean; localization is another route. Brownian scaling gives $`E[W_t^{2n}]=t^nm_{2n}`$, where $`m_{2n}=E[W_1^{2n}]`$. Taking expectations and integrating gives

```math
t^nm_{2n}=n(2n-1)m_{2n-2}\int_0^ts^{n-1}ds
=(2n-1)m_{2n-2}t^n.
```

Starting from $`m_0=1`$, obtain $`m_{2n}=(2n-1)(2n-3)\cdots1`$. Odd moments vanish by symmetry.

#### Problem SC-16: pairings and the same product

Source: [Lalley, Homework 4, Problem 1(B)](https://galton.uchicago.edu/~lalley/Courses/385/HW4.pdf).

Problem. Count partitions of $`2n`$ labeled objects into unordered pairs and relate the count to the Gaussian even-moment product.

Solution. To pair $`2n`$ labeled objects, the partner of object one has $`2n-1`$ choices. After fixing that pair, the smallest unpaired object has $`2n-3`$ partner choices, and so on. Multiplying gives $`(2n-1)(2n-3)\cdots1`$. No overcount occurs because selecting the smallest remaining label fixes a unique order in which pairs are formed. This proves the counting statement independently of the Gaussian moment recursion.

#### Problem SC-17: an exponential concentration bound

Source: [Lalley, Homework 4, Problem 2](https://galton.uchicago.edu/~lalley/Courses/385/HW4.pdf).

Problem. Let a continuous Itô-integral local martingale start at zero and satisfy $`[M]_\infty\le1`$. Prove $`P(\sup_tM_t\ge x)\le e^{-x^2/2}`$ for $`x\ge0`$.

Solution. Let $`M_0=0`$ be a continuous Itô-integral local martingale with $`[M]_\infty\le1`$. For $`\theta\gt 0`$, the exponential $`Z_t=\exp(\theta M_t-\theta^2[M]_t/2)`$ is a nonnegative local martingale and hence a supermartingale. Stop at the first time $`M`$ reaches $`x\ge0`$, capped at a finite horizon and localized if necessary. On the hitting event, $`Z\ge e^{\theta x-\theta^2/2}`$.

The supermartingale bound gives $`P(\sup_{t\le T}M_t\ge x)\le e^{-\theta x+\theta^2/2}`$. Let $`T\to\infty`$ using increasing events and choose $`\theta=x`$ to obtain $`e^{-x^2/2}`$. For $`x=0`$ the bound is one, obtained by continuity or directly. Localization and nonnegativity justify the expectation inequality without assuming an unrestricted exponential is automatically a true martingale.

### 4.5 Additional solved problems

#### Problem SC-18: the logarithm of a positive Itô process

Constructed problem. If $`dX_t=\mu X_tdt+\sigma X_tdW_t`$ and $`X_t\gt 0`$, compute $`d\log X_t`$.

Solution. With $`f'=1/x`$ and $`f''=-1/x^2`$, Itô gives $`(\mu-\sigma^2/2)dt+\sigma dW_t`$. The correction is negative because logarithm is concave. This calculation is valid before any potential hit of zero; the explicit solution in Section 5 then proves positivity for positive initial value and justifies its global use on finite horizons.

#### Problem SC-19: a quartic martingale

Constructed problem. Find a polynomial martingale of form $`W_t^4+a\,tW_t^2+b\,t^2`$.

Solution. Itô gives $`dW^4=4W^3dW+6W^2dt`$ and $`d(tW^2)=W^2dt+t(2WdW+dt)`$. The total drift is $`(6+a)W^2+(a+2b)t`$. Set $`a=-6`$, $`b=3`$. The stochastic differential becomes $`(4W^3-12tW)dW`$, whose finite-horizon second moment is integrable by Gaussian moments. Thus $`W_t^4-6tW_t^2+3t^2`$ is a true martingale.

#### Problem SC-20: correlated squared coordinates

Constructed problem. Let $`W,B`$ have correlation $`\rho`$. Find $`E[W_t^2B_t^2]`$.

Solution. Write $`B_t=\rho W_t+\sqrt{1-\rho^2}Z_t`$ with independent Brownian $`Z`$. Expanding and using zero odd moments gives $`\rho^2E[W_t^4]+(1-\rho^2)E[W_t^2]E[Z_t^2]`$. This equals $`3\rho^2t^2+(1-\rho^2)t^2=(1+2\rho^2)t^2`$. At $`\rho=0`$ it factors; at $`\rho=\pm1`$ it becomes the fourth normal moment.

## 5. Stochastic differential equations

### 5.1 Theory

An SDE $`dX_t=b(t,X_t)dt+\sigma(t,X_t)dW_t`$ means

```math
X_t=X_0+\int_0^tb(s,X_s)ds+\int_0^t\sigma(s,X_s)dW_s.
```

A strong solution is adapted to a specified Brownian driver and initial information on a given probability space. A weak solution includes a choice of probability space and Brownian motion supporting the required law. Pathwise uniqueness compares solutions driven by the same Brownian motion; uniqueness in law compares their distributions. These are separate definitions, although standard theorems connect them under appropriate hypotheses.

Globally Lipschitz coefficients with linear growth are standard sufficient conditions for a unique nonexplosive strong solution with suitable moment bounds. Locally Lipschitz coefficients give a local solution up to potential explosion or boundary hitting; those possibilities require further analysis.

Geometric Brownian motion has proportional drift and diffusion. Ornstein-Uhlenbeck dynamics have linear mean reversion and additive noise. They illustrate why the noise coefficient matters: proportional noise preserves positivity in the explicit GBM solution, while an OU process can take negative values.

### 5.2 Derivations and proofs

#### Geometric Brownian motion

For $`dS_t=\mu S_tdt+\sigma S_tdW_t`$, apply the log calculation and integrate:

```math
S_t=S_0\exp[(\mu-\tfrac12\sigma^2)t+\sigma W_t].
```

For $`S_0\gt 0`$, this formula is strictly positive at every finite time. Applying Itô to the exponential verifies that it solves the original SDE. The coefficients are globally Lipschitz and of linear growth, so the standard strong existence and uniqueness theorem applies too.

Gaussian exponential moments give

```math
ES_t=S_0e^{\mu t},\qquad
\mathop{\mathrm{Var}} (S_t)=S_0^2e^{2\mu t}(e^{\sigma^2t}-1).
```

The drift of the log is $`\mu-\sigma^2/2`$, while the growth rate of the arithmetic mean is $`\mu`$. Confusing them loses the Itô correction.

#### Ornstein-Uhlenbeck process

For $`dX_t=\kappa(\theta-X_t)dt+\sigma dW_t`$, $`\kappa\gt 0`$, multiply by $`e^{\kappa t}`$. Since that factor has finite variation, the product rule gives

```math
d[e^{\kappa t}(X_t-\theta)]=\sigma e^{\kappa t}dW_t.
```

Therefore

```math
X_t=\theta+(X_0-\theta)e^{-\kappa t}
+\sigma\int_0^te^{-\kappa(t-s)}dW_s.
```

Conditional on $`\mathcal F_s`$, the future integral is independent Gaussian with variance $`\sigma^2(1-e^{-2\kappa(t-s)})/(2\kappa)`$. Its conditional mean is zero. This gives the full conditional normal law. A stationary initialization is $`X_0\sim N(\theta,\sigma^2/(2\kappa))`$, independent of future noise. Starting from a fixed value approaches that law but is not stationary at the initial dates.

#### Existence, uniqueness, and the role of growth conditions

A strong solution is constructed on a space with a specified Brownian motion and is adapted to its given information, together with the initial condition. A weak solution includes a choice of probability space and Brownian motion as part of the solution. These notions concern the stochastic model, not the numerical accuracy of an approximation.

A standard sufficient theorem gives a unique global strong solution when drift and diffusion are globally Lipschitz in the state, uniformly on each finite time interval, with a suitable linear-growth bound and square-integrable initial state. A proof constructs successive integral approximations and controls their differences using the Itô isometry and deterministic integral inequalities.

The uniqueness estimate can be seen directly. Let $`X,Y`$ solve the same equation with the same Brownian motion and initial state, and suppose both coefficients have Lipschitz constant $`L`$. On $`[0,T]`$, using $`\vert a+b\vert ^2\le2\vert a\vert ^2+2\vert b\vert ^2`$, Cauchy–Schwarz for the drift integral, and the isometry for the stochastic integral gives

```math
E\vert X_t-Y_t\vert ^2
\le 2L^2(T+1)\int_0^t E\vert X_s-Y_s\vert ^2ds.
```

Gronwall's inequality forces the difference to be zero. The linear-growth condition separately prevents moments from blowing up in finite time.

Local regularity alone is insufficient for a global solution. The deterministic equation $`dX=X^2dt`$ from $`X_0\gt 0`$ has the locally smooth solution $`X_t=X_0/(1-X_0t)`$, which explodes at $`1/X_0`$. Conversely, lack of Lipschitz continuity can allow nonuniqueness even without noise: $`dX=\sqrt{\vert X\vert }dt`$ from zero has solutions that remain zero until any chosen date $`c`$ and then equal $`(t-c)^2/4`$.

GBM has globally Lipschitz state coefficients and an explicit positive solution. OU coefficients are also globally Lipschitz. The nonlinear positive example solved earlier needs its explicit transformation to establish positivity and nonexplosion; it is not covered by a global Lipschitz assertion at zero.

#### Example SC-E25: a time-varying drift with additive noise

Solve $`dX_t=(t-X_t)dt+2dW_t`$ with $`X_0=x`$.

Solution. Multiply by $`e^t`$. The product rule gives $`d(e^tX_t)=te^tdt+2e^tdW_t`$. Integrating and multiplying back yields

```math
X_t=xe^{-t}+\int_0^t e^{-(t-s)}s\,ds
+2\int_0^t e^{-(t-s)}dW_s.
```

The ordinary integral equals $`t-1+e^{-t}`$, by integration of $`se^s`$. The stochastic integral has a deterministic kernel. Therefore

```math
X_t\sim N\left(xe^{-t}+t-1+e^{-t},\ 2(1-e^{-2t})\right).
```

The mean satisfies $`m'(t)=t-m(t)`$ and $`m(0)=x`$, giving a deterministic check. The variance satisfies $`v'(t)=-2v(t)+4`$ and $`v(0)=0`$, giving another. The process has a bounded limiting variance but no stationary mean, since its mean grows approximately as $`t-1`$. Mean-reverting dependence on the current state alone does not imply stationarity when the forcing term changes with time.

### 5.3 Worked examples

#### Example SC-E12: a GBM conditional distribution

Given $`S_t=s\gt 0`$, the log ratio over horizon $`h`$ is normal with mean $`(\mu-\sigma^2/2)h`$ and variance $`\sigma^2h`$. Thus $`P(S_{t+h}\gt K\mid\mathcal F_t)=1-\Phi([\log(K/s)-(\mu-\sigma^2/2)h]/[\sigma\sqrt h])`$ for $`\sigma\gt 0`$. The conditional statement uses Brownian increments after $`t`$.

#### Example SC-E13: mean reversion with a current state

Let $`\kappa=2`$, $`\theta=1`$, $`\sigma=3`$, and $`X_t=4`$. At horizon $`h`$, the conditional mean is $`1+3e^{-2h}`$ and variance $`9(1-e^{-4h})/4`$. The mean approaches one, while the variance approaches $`9/4`$. Mean reversion does not make long-horizon uncertainty disappear.

#### Example SC-E14: an additive-noise SDE

For $`dX=a\,dt+\sigma dW`$, the solution is $`X_t=X_0+at+\sigma W_t`$. At deterministic initial value, its law is normal with mean $`X_0+at`$ and variance $`\sigma^2t`$. There is no Itô drift correction in $`X`$ itself because the solution is linear in $`W`$; nonlinear functions of this solution still acquire corrections.

### 5.4 Selected exercises

#### Problem SC-21: a coupled linear SDE

Source: [Lalley, Homework 5, Problem 1](https://galton.uchicago.edu/~lalley/Courses/385/HW5.pdf).

Problem. Solve $`dX_1=X_2dt+\alpha dW_1`$ and $`dX_2=X_1dt+\beta dW_2`$ for independent Brownian coordinates and given initial values, and examine dependence after diagonalization.

Solution. The equations are $`dX_1=X_2dt+\alpha dW_1`$ and $`dX_2=X_1dt+\beta dW_2`$, with independent Brownian coordinates. Put $`U=X_1+X_2`$ and $`V=X_1-X_2`$. Then

```math
dU=Udt+\alpha dW_1+\beta dW_2,\qquad
dV=-Vdt+\alpha dW_1-\beta dW_2.
```

Integrating factors give

```math
U_t=e^tU_0+\int_0^te^{t-s}(\alpha dW_{1s}+\beta dW_{2s}),
```

```math
V_t=e^{-t}V_0+\int_0^te^{-(t-s)}(\alpha dW_{1s}-\beta dW_{2s}).
```

Recover $`X_1=(U+V)/2`$ and $`X_2=(U-V)/2`$. These formulas completely specify the jointly Gaussian solution. The transformed noises have cross-variation $`(\alpha^2-\beta^2)dt`$, so $`U,V`$ are not generally independent. The positive drift eigenvalue in $`U`$ also shows that the system is not mean-reverting as a whole.

#### Problem SC-22: a stochastic integrating factor

Source: [Lalley, Homework 5, Problem 2(A)](https://galton.uchicago.edu/~lalley/Courses/385/HW5.pdf).

Problem. For $`dX=f(t,X)dt+a(t)X\,dW`$, derive a stochastic integrating factor and the corresponding random ordinary equation, under the usual finite-horizon integrability assumptions.

Solution. For $`dX=f(t,X)dt+a(t)X\,dW`$, define

```math
G_t=\exp\left[-\int_0^ta(s)dW_s+\frac12\int_0^ta(s)^2ds\right].
```

Itô gives $`dG=-aGdW+a^2Gdt`$. The product rule then gives

```math
d(GX)=Gf\,dt+aGX\,dW-aGX\,dW+a^2GXdt-a^2GXdt
=Gf(t,X)dt.
```

The cross-variation cancels the extra integrating-factor drift. Defining $`Y=GX`$ therefore gives the random ordinary equation $`dY=Gf(t,G^{-1}Y)dt`$. This fixes the transformed-variable notation explicitly; it is $`GX`$, not $`G^{-1}X`$, for the displayed choice of $`G`$.

#### Problem SC-23: a positive nonlinear SDE

Source: [Lalley, Homework 5, Problem 2(C)](https://galton.uchicago.edu/~lalley/Courses/385/HW5.pdf).

Problem. Solve $`dX=X^{-1}dt+\alpha XdW`$ from $`X_0\gt 0`$ and justify positivity and nonexplosion on finite horizons.

Solution. For $`dX=X^{-1}dt+\alpha XdW`$ and $`X_0=x_0\gt 0`$, use $`G_t=e^{-\alpha W_t+\alpha^2t/2}`$. From SC-22, $`Y=GX`$ satisfies $`dY=G^2Y^{-1}dt`$. Thus $`d(Y^2)=2G^2dt`$ and

```math
X_t=G_t^{-1}
\left[x_0^2+2\int_0^tG_s^2ds\right]^{1/2}.
```

The square-root argument stays strictly positive, and $`G`$ is positive and continuous, so the formula gives a finite positive solution on every finite horizon. Local uniqueness follows from local Lipschitz coefficients on $`(0,\infty)`$; the explicit expression rules out both hitting zero and finite-time explosion. Applying the product calculation backward verifies the SDE.

### 5.5 Additional solved problems

#### Problem SC-24: a stationary OU autocovariance

Constructed problem. Initialize OU in its stationary normal distribution and compute $`\mathop{\mathrm{Cov}} (X_t,X_{t+h})`$.

Solution. The conditional future expression is $`X_{t+h}-\theta=e^{-\kappa h}(X_t-\theta)`$ plus independent centered future noise. Multiplying by $`X_t-\theta`$ and taking expectations gives $`\sigma^2e^{-\kappa h}/(2\kappa)`$. The correlation is $`e^{-\kappa h}`$. Discretely sampled stationary OU is therefore an AR(1) with coefficient $`e^{-\kappa h}`$ and Gaussian innovation variance $`\sigma^2(1-e^{-2\kappa h})/(2\kappa)`$.

#### Problem SC-25: moments of a stochastic exponential

Constructed problem. Find $`E[M_t^p]`$ for $`M_t=e^{\theta W_t-\theta^2t/2}`$ and real $`p`$.

Solution. Gaussian exponential moments give

```math
E[M_t^p]
=e^{-p\theta^2t/2}e^{p^2\theta^2t/2}
=e^{p(p-1)\theta^2t/2}.
```

For $`p=1`$ the mean is one. For $`p=2`$ the second moment is $`e^{\theta^2t}`$, finite on finite horizons but growing over time. Finite-horizon martingale properties do not automatically imply a uniformly integrable family over an infinite horizon.

#### Problem SC-26: a reciprocal GBM

Constructed problem. Find the SDE for $`Y=1/S`$ when $`S`$ follows positive GBM.

Solution. With $`f'=-1/s^2`$, $`f''=2/s^3`$, Itô gives

```math
dY=(-\mu+\sigma^2)Ydt-\sigma YdW.
```

The reciprocal does not merely negate the original drift. Its explicit formula is $`S_0^{-1}e^{(-\mu+\sigma^2/2)t-\sigma W_t}`$, consistent with the GBM solution for drift $`-\mu+\sigma^2`$ and volatility $`-\sigma`$.

#### Problem SC-27: a deterministic limit

Constructed problem. Let the OU volatility tend to zero with deterministic initial value.

Solution. The stochastic integral has variance $`\sigma^2(1-e^{-2\kappa t})/(2\kappa)`$, tending to zero. The solution converges in $`L^2`$ at each fixed time to $`\theta+(X_0-\theta)e^{-\kappa t}`$, the deterministic mean-reversion ODE solution. A stochastic formula should reproduce this zero-noise limit.

## 6. Risk-neutral valuation

### 6.1 Theory

Let the bank account satisfy $`dB_t=rB_tdt`$, $`B_0=1`$, and a nondividend-paying stock satisfy $`dS_t=\mu S_tdt+\sigma S_tdW_t^P`$, with constant $`\sigma\gt 0`$. Under the physical measure $`P`$, $`\mu`$ describes the modeled expected stock return. A pricing measure changes probabilities so discounted traded prices become martingales under the required integrability assumptions.

Equivalent measures have the same probability-zero events on the specified finite-horizon sigma-algebra. They can assign different probabilities to other events. Equivalence preserves almost-sure path properties such as quadratic variation; it does not preserve drift or expectations.

The market price of Brownian risk is $`\theta=(\mu-r)/\sigma`$ in this one-factor model. A Girsanov measure change can shift the Brownian drift while leaving the diffusion coefficient unchanged. For a progressively measurable $`\theta_t`$, a standard sufficient condition for the density exponential to be a true martingale is Novikov's condition,

```math
E_P\exp\left(\frac12\int_0^T\theta_t^2dt\right)\lt \infty.
```

Girsanov's theorem is stated and applied here, not proved in full generality. Under its hypotheses, density

```math
Z_T=\exp\left[-\int_0^T\theta_tdW_t^P-\frac12\int_0^T\theta_t^2dt\right]
```

defines $`Q`$ by $`dQ/dP=Z_T`$, and $`W_t^Q=W_t^P+\int_0^t\theta_sds`$ is Brownian under $`Q`$ on that horizon.

### 6.2 Derivations and proofs

#### Constant-parameter change of measure

For constant $`\theta`$, Novikov's expression is the finite constant $`e^{\theta^2T/2}`$. Gaussian integration gives $`E_PZ_T=1`$, and $`Z_T\gt 0`$ almost surely, so $`Q`$ is equivalent to $`P`$ on $`\mathcal F_T`$.

Substitute $`dW^P=dW^Q-\theta dt`$:

```math
dS_t=(\mu-\sigma\theta)S_tdt+\sigma S_tdW_t^Q
=rS_tdt+\sigma S_tdW_t^Q.
```

Then

```math
e^{-rt}S_t=S_0e^{-\sigma^2t/2+\sigma W_t^Q},
```

a true finite-horizon martingale by its explicit Gaussian conditional expectation. Showing zero drift alone would establish only a local-martingale statement without the integrability check.

#### Pricing a replicated payoff

For a self-financing strategy with integrable discounted gain process and terminal value $`H`$, the discounted portfolio is a $`Q`$-martingale. Therefore its value at date $`t`$ is

```math
V_t=E_Q[e^{-r(T-t)}H\mid\mathcal F_t].
```

In the standard Brownian stock-bank model with nonzero volatility and the Brownian filtration, suitable square-integrable claims can be replicated using martingale representation. That theorem is invoked under these assumptions. In an incomplete market, equivalent pricing measures may be nonunique, and the expectation formula needs a chosen measure or an additional pricing principle.

#### A one-period replication before continuous time

Let a stock cost $`S_0`$ and end at $`uS_0`$ or $`dS_0`$, with $`d\lt R\lt u`$, where $`R`$ is the money-market gross return. A claim pays $`H_u,H_d`$. Holding $`\Delta`$ shares and an initial cash position $`B_0`$ gives the two equations

```math
\Delta uS_0+RB_0=H_u,\qquad
\Delta dS_0+RB_0=H_d.
```

Subtracting gives $`\Delta=(H_u-H_d)/(S_0(u-d))`$. Substitution determines $`B_0`$ and the initial replication cost:

```math
v_0=\Delta S_0+B_0
=\frac1R[qH_u+(1-q)H_d],
\qquad q=\frac{R-d}{u-d}.
```

The inequalities make $`0\lt q\lt 1`$. Under this probability, the discounted stock has unchanged expectation. Replication derived the price before introducing this probability; the probability is a convenient representation of the resulting linear pricing rule.

For $`S_0=100`$, $`u=1.2`$, $`d=0.8`$, $`R=1.05`$, and call strike one hundred, the payoffs are twenty and zero. The hedge is half a share, initial cash is $`-40/1.05`$, and price is $`12.5/1.05\approx11.904762`$. The pricing probability is $`q=0.625`$. The physical probability of an up move was not needed for this complete replication.

If the physical up probability is $`p\in(0,1)`$, the pricing density takes values $`q/p`$ and $`(1-q)/(1-p)`$. It is positive and has physical expectation one. This finite model illustrates what an equivalent measure does: it reweights the same possible outcomes.

#### Discounting, dividends, and completeness

With money-market account $`B_t=\exp(\int_0^t r_sds)`$, its reciprocal has finite variation. Applying the product rule gives

```math
d(S_t/B_t)=B_t^{-1}dS_t-r_tB_t^{-1}S_tdt.
```

Under a pricing measure, a no-dividend stock with drift $`r_tS_t`$ therefore has zero discounted drift. If the stock pays dividends, the discounted cumulative gains process includes those dividends. The ex-dividend price alone generally has a different drift.

In the standard one-Brownian Black–Scholes market with nonzero volatility and the Brownian-generated filtration, a martingale representation theorem underlies replication of suitable integrable claims. This is an advanced theorem stated and used under its assumptions. If the information includes an independent untraded Brownian shock, a stock driven only by the first shock cannot replicate every payoff depending on the second. Equivalent martingale measures can then be nonunique.

Thus constructing one pricing measure does not by itself prove market completeness or a unique price for every conceivable claim. Replication, admissible strategies, the information set, and the number of independent traded risk exposures all matter.

#### Example SC-E26: zero drift after the wrong discounting

Under $`Q`$ let $`dS=rSdt+\sigma SdW`$. Is $`e^{-rt}S_t^2`$ a martingale?

Solution. Itô gives $`d(S^2)=(2r+\sigma^2)S^2dt+2\sigma S^2dW`$. Multiplying by $`e^{-rt}`$ leaves drift $`(r+\sigma^2)e^{-rt}S^2dt`$. It is generally nonzero. Discounting every function of a traded price by the money-market account does not automatically make that function a martingale.

The process $`e^{-(2r+\sigma^2)t}S_t^2`$ does have zero drift and is a true martingale on finite horizons by its explicit exponential-normal form. This altered factor is a moment normalization, not the money-market discount factor used to price arbitrary payoffs.

### 6.3 Worked examples

#### Example SC-E15: computing the drift shift

For $`\mu=0.10`$, $`r=0.04`$, and $`\sigma=0.20`$, the risk price is $`\theta=0.3`$. The physical log drift is $`0.08`$, while the risk-neutral log drift is $`0.02`$. The instantaneous variance remains $`0.04`$ per unit time. These measures model different expectations for the same finite-horizon path events.

#### Example SC-E16: a forward payoff

For $`H=S_T-K`$, risk-neutral expectation gives $`V_t=S_t-Ke^{-r(T-t)}`$. A zero-cost forward therefore has delivery price $`K=S_te^{r(T-t)}`$. The physical drift $`\mu`$ does not enter because the payoff is replicated by holding the stock and financing the delivery payment.

#### Example SC-E17: a dividend yield

If the stock pays continuous proportional yield $`q`$, its ex-dividend risk-neutral drift is $`r-q`$. The discounted stock price alone is generally not the relevant martingale; discounted gains include the dividends. The forward delivery price becomes $`S_te^{(r-q)(T-t)}`$ under the constant-yield model.

### 6.4 Selected exercise and additional problems

#### Problem SC-28: first passage to a declining line

Source: [Lalley, Brownian Motion, Exercise 8](https://galton.uchicago.edu/~lalley/Courses/385/BrownianMotion.pdf), the Laplace-transform option.

Problem. For Brownian motion, let $`\tau`$ be the first hit of $`a-bt`$ with $`a,b\gt 0`$. Derive its Laplace transform and use it to obtain its passage probability, mean, and variance.

Solution. The hitting condition $`W_t=a-bt`$ is equivalent to $`X_t=W_t+bt`$ reaching $`a`$, with $`a,b\gt 0`$. For $`\lambda\gt 0`$, choose $`q=\sqrt{b^2+2\lambda}-b\gt 0`$, so $`bq+q^2/2=\lambda`$. Then $`M_t=e^{qX_t-\lambda t}=e^{qW_t-q^2t/2}`$ is a martingale.

Before first passage, $`X_t\lt a`$, so $`M_{\tau\wedge T}\le e^{qa}`$. Optional stopping at the bounded time gives expectation one. On $`\lbrace \tau\gt T\rbrace`$, the terminal value is at most $`e^{qa-\lambda T}`$, tending to zero. Dominated convergence gives

```math
E[e^{-\lambda\tau};\tau\lt \infty]
=e^{-qa}=\exp[a(b-\sqrt{b^2+2\lambda})].
```

Letting $`\lambda\downarrow0`$ gives passage probability one. Differentiating the transform at zero yields $`E\tau=a/b`$ and $`\mathop{\mathrm{Var}} (\tau)=a/b^3`$. This uses an exponential martingale with a checked stopping bound rather than an unrestricted optional-stopping assertion.

#### Problem SC-29: a failed density candidate

Constructed problem. Explain why a nonnegative local martingale density is insufficient unless its expectation is one.

Solution. A probability measure defined by $`Q(A)=E_P[Z_T\mathbf1_A]`$ must satisfy $`Q(\Omega)=E_PZ_T=1`$. A strict local martingale can have expectation below its initial value. It would then define a subprobability measure through that formula. A condition such as Novikov's, or another valid true-martingale argument, is needed before using it as a change-of-probability density.

#### Problem SC-30: a discounted martingale check by conditioning

Constructed problem. Verify the martingale property of discounted risk-neutral GBM directly for $`s\lt t`$.

Solution. Factor

```math
e^{-rt}S_t=e^{-rs}S_s
\exp[-\sigma^2(t-s)/2+\sigma(W_t^Q-W_s^Q)].
```

The second factor is independent of $`\mathcal F_s`$ and has mean one by the normal moment-generating function. Conditioning gives $`E_Q[e^{-rt}S_t\mid\mathcal F_s]=e^{-rs}S_s`$. This verifies the conditional property and finite expectations at once.

#### Problem SC-31: an extra untraded shock

Constructed problem. A claim depends on an independent Brownian driver not spanned by the stock. Explain the pricing implication.

Solution. One stock diffusion and a bank account span only the stock's Brownian risk direction. A payoff with an independent noise component generally cannot be replicated by those assets alone. Conditions making the discounted stock a martingale constrain its risk premium but can leave the independent shock's pricing adjustment unspecified. Equivalent martingale measures can then assign different expectations to that claim.

## 7. Pricing PDEs and Feynman-Kac

### 7.1 Theory

A pricing PDE describes a value as a function of current time and state. A risk-neutral expectation describes the same value as an average over future paths, when the relevant existence, smoothness, and integrability conditions hold. Terminal and boundary conditions are part of the problem.

For $`dX=b(t,X)dt+\sigma(t,X)dW`$ under a chosen measure, a smooth solution of

```math
u_t+b u_x+\tfrac12\sigma^2u_{xx}-c\,u=0,\qquad u(T,x)=g(x),
```

can be represented as an expected discounted terminal payoff. Sufficient classical hypotheses include well-posed nonexplosive state dynamics, regular coefficients, a $`C^{1,2}`$ solution for $`t\lt T`$, suitable growth bounds, and integrability that makes the derived stochastic integral a true martingale. Less regular payoffs require an appropriate limiting or weak/viscosity formulation.

### 7.2 Derivations and proofs

#### Black-Scholes replication

Let the claim value be $`v(t,S_t)`$, smooth before maturity. Itô under the physical measure gives

```math
dv=(v_t+\mu Sv_S+\tfrac12\sigma^2S^2v_{SS})dt
+\sigma Sv_SdW^P.
```

A self-financing portfolio with $`\Delta`$ shares and bank position $`v-\Delta S`$ has

```math
dv=\Delta\,dS+r(v-\Delta S)dt.
```

Match diffusion terms to obtain $`\Delta=v_S`$. Matching drifts then cancels the physical-drift term $`\mu Sv_S`$ and yields

```math
v_t+\frac12\sigma^2S^2v_{SS}+rSv_S-rv=0.
```

The bank position finances changes in stock holdings; self-financing is not the same as holding both positions constant. The terminal condition for a call is $`(S-K)^+`$.

#### Feynman-Kac verification

Fix starting time $`t`$ and define $`D_s=\exp[-\int_t^s c(v,X_v)dv]`$. Itô and the finite-variation product rule give

```math
d[D_su(s,X_s)]
=D_s(u_t+bu_x+\tfrac12\sigma^2u_{xx}-cu)ds
+D_s\sigma u_xdW_s.
```

The PDE makes the drift zero. If the stochastic integral is a true martingale and the terminal limit is integrable, condition at time $`t`$ to obtain

```math
u(t,x)=E_{t,x}\left[e^{-\int_t^Tc(s,X_s)ds}g(X_T)\right].
```

This proves the representation for an already regular solution. Proving that the expectation is a classical solution requires additional regularity results; the stochastic calculation alone does not establish all PDE existence theory.

#### A European call from the Gaussian integral

Under risk-neutral GBM, let $`\tau=T-t\gt 0`$ and $`S_T=s\exp[(r-\sigma^2/2)\tau+\sigma\sqrt\tau Z]`$, $`Z\sim N(0,1)`$. The call payoff is positive when $`Z\gt -d_2`$, where

```math
d_2=\frac{\log(s/K)+(r-\sigma^2/2)\tau}{\sigma\sqrt\tau},
\qquad d_1=d_2+\sigma\sqrt\tau.
```

Split the discounted payoff expectation into $`e^{-r\tau}E[S_T\mathbf1_{\lbrace Z\gt -d_2\rbrace }]-Ke^{-r\tau}P(Z\gt -d_2)`$. The second probability is $`\Phi(d_2)`$. In the first integral, completing the normal square gives

```math
e^{a z-a^2/2}\varphi(z)=\varphi(z-a),
\qquad a=\sigma\sqrt\tau.
```

Thus the first term is $`s\Phi(d_1)`$ and

```math
C(t,s)=s\Phi(d_1)-Ke^{-r\tau}\Phi(d_2).
```

The payoff kink at maturity does not prevent smoothness for $`s\gt 0`$ and $`\tau\gt 0`$ in this nondegenerate model. Terminal convergence and growth conditions connect the formula to the original payoff.

#### Feynman–Kac with intermediate cash flows

Let $`dX_s=b(s,X_s)ds+\sigma(s,X_s)dW_s`$ and write $`\mathcal Lv=bv_x+\tfrac12\sigma^2v_{xx}`$. Suppose a sufficiently regular candidate satisfies

```math
v_t+\mathcal Lv-rv+f=0,\qquad v(T,x)=g(x).
```

Here $`f`$ is a running payment rate and $`g`$ the terminal payment. Set $`D_s=\exp(-\int_t^s r(u,X_u)du)`$. The product rule and Itô's formula show that

```math
M_s=D_sv(s,X_s)+\int_t^sD_uf(u,X_u)du
```

has zero drift. If its stochastic-integral term is a true martingale and the remaining terms are integrable, taking conditional expectations gives

```math
v(t,x)=E_{t,x}\left[
D_Tg(X_T)+\int_t^TD_sf(s,X_s)ds
\right].
```

Sufficient conditions can be supplied by bounded coefficients and derivatives, or by polynomial-growth bounds together with the corresponding finite state moments. In a particular application these must actually be checked. If the stochastic integral is only local, taking expectations without localization and a justified limit can fail.

The same argument with the difference of two candidate solutions gives uniqueness within a class where the martingale and terminal limits are valid. A PDE and terminal condition written without a growth or integrability class need not identify a unique solution.

For a process killed on leaving a domain, stop at the earlier of exit and maturity. The representation includes any specified payment on the exit boundary, discounted to the exit date. Boundary conditions therefore encode an economic or probabilistic rule, rather than merely helping a numerical solver.

#### The call hedge from its closed-form price

For a no-dividend European call, $`v=s\Phi(d_1)-Ke^{-r\tau}\Phi(d_2)`$ with $`d_2=d_1-\sigma\sqrt\tau`$. The normal-density identity

```math
s\phi(d_1)=Ke^{-r\tau}\phi(d_2)
```

follows by substituting the definitions and comparing the exponents. Differentiating the price in $`s`$ creates two density terms. The identity cancels them, leaving

```math
\Delta=\Phi(d_1),\qquad
\Gamma=\frac{\phi(d_1)}{s\sigma\sqrt\tau}.
```

Delta lies between zero and one, and gamma is positive for positive $`s,\sigma,\tau`$. These are useful checks on a derived hedge. At maturity the payoff has a corner at the strike, so classical second derivatives need not persist there even though the price is smooth for earlier dates.

For the call, the zero-stock boundary value is zero and the large-stock behavior is $`s-Ke^{-r\tau}`$. A digital payoff has a discontinuous terminal condition and different limiting behavior. Smoothness for $`t\lt T`$ does not make every terminal payoff smooth.

#### Example SC-E27: continuous coupons and a terminal payment

With constant rate $`r`$, a claim pays continuously at rate $`c`$ until $`T`$ and one unit at $`T`$. Find its value.

Solution. There is no remaining uncertainty about the payments. Discounting gives

```math
v(t)=e^{-r\tau}+c\int_0^\tau e^{-rs}ds
=e^{-r\tau}+\frac c r(1-e^{-r\tau})
```

for $`r\ne0`$, with limit $`1+c\tau`$ for $`r=0`$. Its equation is $`v_t-rv+c=0`$ and terminal value one. Direct differentiation verifies both. The plus sign before the payment rate in the PDE matches the positive discounted cash-flow integral. Using the homogeneous pricing PDE would omit the coupons.

### 7.3 Worked examples

#### Example SC-E18: a quadratic payoff

For terminal payoff $`S_T^2`$, the risk-neutral second moment gives $`v(t,s)=s^2e^{(r+\sigma^2)(T-t)}`$. Differentiate: $`v_t=-(r+\sigma^2)v`$, $`sv_s=2v`$, and $`s^2v_{ss}=2v`$. Substitution into the Black-Scholes PDE gives $`[-r-\sigma^2+\sigma^2+2r-r]v=0`$, verifying it directly.

#### Example SC-E19: a heat equation

For $`u_t+\tfrac12u_{xx}=0`$ with terminal $`g(x)=x^2`$, Brownian expectation gives $`u(t,x)=x^2+T-t`$. Its derivatives give $`u_t=-1`$ and $`u_{xx}=2`$, satisfying the equation. The added time term is the future Brownian variance.

#### Example SC-E20: a deterministic short rate

If $`r(s)`$ is deterministic, a zero-coupon bond paying one at $`T`$ has price $`\exp[-\int_t^Tr(s)ds]`$. With stochastic rates, the bond is $`E_Q[e^{-\int_t^Tr_sds}\mid\mathcal F_t]`$, which generally differs from the exponential of the expected integral.

### 7.4 Selected exercises

For the next three problems, let Brownian motion start at $`x\in[a,b]`$ and let $`T`$ be its first exit from $`(a,b)`$.

#### Problem SC-32: a stopped Feynman-Kac martingale

Source: [Lalley, Homework 4, Problem 5(A)](https://galton.uchicago.edu/~lalley/Courses/385/HW4.pdf).

Problem. Let $`T`$ be Brownian exit from $`(a,b)`$ and suppose $`\phi''/2=k\phi`$ with continuous $`k\ge0`$ and $`\phi\in C^2`$ on a neighborhood of the closed interval. Show that the discounted, stopped $`\phi(W)`$ process is a bounded martingale.

Solution. Assume $`\phi\in C^2`$, $`k`$ continuous and nonnegative, and $`\phi''/2=k\phi`$. Set $`Z_t=\phi(W_t)e^{-\int_0^tk(W_s)ds}`$. Itô and the product rule leave

```math
dZ_t=e^{-\int_0^tk(W_s)ds}\phi'(W_t)dW_t.
```

Stop at $`T`$. On the compact interval $`[a,b]`$, $`\phi`$ and $`\phi'`$ are bounded, and the exponential is at most one. On each finite horizon the stopped integrand is square-integrable, so the stopped process is a true martingale. Also $`\vert Z_{t\wedge T}\vert \le\max_{[a,b]}\vert \phi\vert`$, a uniform bound giving $`L^2`$ boundedness and uniform integrability.

#### Problem SC-33: terminal boundary representation

Source: [Lalley, Homework 4, Problem 5(B)](https://galton.uchicago.edu/~lalley/Courses/385/HW4.pdf).

Problem. Use the stopped martingale of SC-32 to derive its boundary expectation representation at the almost-surely finite interval exit.

Solution. Brownian exit from a bounded interval is almost surely finite. For example, on every unit interval an increment larger in magnitude than $`b-a`$ forces exit; independence supplies a fixed positive conditional chance and a geometric survival bound. By continuity $`W_T\in\lbrace a,b\rbrace`$. The bounded stopped martingale from SC-32 converges almost surely, and dominated convergence preserves its mean:

```math
\phi(x)=E_x\left[\phi(W_T)
\exp\left(-\int_0^Tk(W_s)ds\right)\right].
```

The finite-interval boundedness and $`k\ge0`$ justify this limit. An unbounded state domain or a negative potential would require a different integrability argument.

#### Problem SC-34: the Laplace transform of interval exit

Source: [Lalley, Homework 4, Problem 5(C)](https://galton.uchicago.edu/~lalley/Courses/385/HW4.pdf).

Problem. For exit time $`T`$ from $`(a,b)`$, solve the constant-potential boundary equation to find $`E_xe^{-\lambda T}`$. Specialize to symmetric endpoints and compute two exit-time moments.

Solution. Choose constant $`k=\lambda\gt 0`$ and boundary values $`\phi(a)=\phi(b)=1`$. The ODE is $`\phi''=2\lambda\phi`$. Put $`m=(a+b)/2`$, $`L=(b-a)/2`$, and $`q=\sqrt{2\lambda}`$. Its boundary solution is

```math
\phi(x)=\frac{\cosh(q(x-m))}{\cosh(qL)}.
```

Direct differentiation checks the ODE and both endpoints. SC-33 gives $`E_x[e^{-\lambda T}]=\phi(x)`$. For symmetric endpoints $`[-a,a]`$ and start zero, this is $`\mathop{\mathrm{sech}} (a\sqrt{2\lambda})`$. Expanding around zero yields $`ET=a^2`$ and $`ET^2=5a^4/3`$, so $`\mathop{\mathrm{Var}} (T)=2a^4/3`$.

### 7.5 Additional solved problems

#### Problem SC-35: an OU short-rate bond

Constructed problem. Under $`Q`$, let $`dr=\kappa(\theta-r)dt+\sigma dW^Q`$, with current rate $`r_t=r`$. Price a unit payment at $`t+\tau`$.

Solution. Integrate the explicit OU solution over the horizon. The integral $`I=\int_t^{t+\tau}r_sds`$ is Gaussian with mean

```math
m=\theta\tau+(r-\theta)B,\qquad
B=\frac{1-e^{-\kappa\tau}}{\kappa},
```

and variance

```math
v=\frac{\sigma^2}{\kappa^2}
\left[\tau-\frac{2(1-e^{-\kappa\tau})}{\kappa}
+\frac{1-e^{-2\kappa\tau}}{2\kappa}\right].
```

The variance follows by stochastic Fubini, whose square-integrability conditions hold for the bounded deterministic kernel on a finite horizon. Gaussian exponential moments give the bond price $`e^{-m+v/2}`$. Replacing it by $`e^{-m}`$ would omit the randomness correction.

#### Problem SC-36: put-call parity

Constructed problem. Derive the European put-call relation under the nondividend stock model.

Solution. Pointwise $`(S_T-K)^+-(K-S_T)^+=S_T-K`$. Linearity of risk-neutral valuation gives $`C-P=S_t-Ke^{-r\tau}`$. The same relation follows by replicating that terminal payoff with one stock and a bank liability. It requires matching strikes, maturities, and European exercise rights.

#### Problem SC-37: a digital payoff

Constructed problem. Price $`\mathbf1_{\lbrace S_T\gt K\rbrace }`$ under constant-parameter risk-neutral GBM.

Solution. The exercise event is the same normal threshold used in the call derivation, with probability $`\Phi(d_2)`$. Discounting gives $`e^{-r\tau}\Phi(d_2)`$. Its payoff is discontinuous at the strike, but for $`\tau\gt 0`$ and $`s\gt 0`$ the normal transition density produces a smooth value. At maturity, convergence at the discontinuity requires a convention; the expectation representation remains unambiguous because a positive-time lognormal law assigns zero probability to a single strike.

#### Problem SC-38: a price boundary

Constructed problem. State sensible call-price boundary behavior under Black-Scholes with $`K\gt 0`$.

Solution. At $`s=0`$, the multiplicative stock model remains at zero and the call is worthless, so $`C(t,0)=0`$. As $`s\to\infty`$ at fixed positive remaining time, both normal probabilities approach one, giving $`C(t,s)\sim s-Ke^{-r\tau}`$. At maturity the terminal payoff is $`(s-K)^+`$. These conditions help specify the PDE problem and its numerical truncation.

#### Problem SC-39: the difference between a local and true martingale in pricing

Constructed problem. Explain the remaining step after obtaining a zero drift in a discounted candidate price.

Solution. The Itô differential shows a local stochastic-integral representation. To take conditional expectations and conclude the current value equals an expected future payoff, prove integrability and the true-martingale property, perhaps by a square-integrable integrand, uniform integrability, bounded stopping, or an explicit conditional expectation. Without this step, the expectation can decrease for a nonnegative strict local martingale, so the desired equality need not hold.

#### Problem SC-40: a nontradable state in Feynman-Kac

Constructed problem. Does a Feynman-Kac expectation for a state diffusion automatically give a unique financial price?

Solution. Feynman-Kac links a PDE and an expectation under a specified diffusion law and discount. Financial pricing also needs a justified pricing measure and assumptions about tradability or market completion. If an unhedgeable state risk admits several pricing adjustments, several different PDE coefficients can be consistent with traded-asset restrictions. The analytical representation does not by itself choose among them.

<!-- extended-practice-start -->

## 8. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 8.1 Conditioning Brownian motion on two observations

An intermediate Brownian value conditioned on both endpoints has less uncertainty than a value conditioned only on the past. The conditional process is a Brownian bridge. This is a conditional distribution calculation; knowledge of a future endpoint changes the information available.

#### Problem SC-41: Brownian bridge at an intermediate date, case 1

For standard Brownian motion, take $`s=0\lt t=\frac{1}{2}\lt T=1`$ and condition on $`W_s=0`$, $`W_T=2`$. Find the conditional law of $`W_t`$ and the probability that $`W_t\gt 1`$.

Solution. Set $`A=W_t-W_s`$ and $`B=W_T-W_t`$. Independent increments give independent normals with variances $`t-s`$ and $`T-t`$. Observing both endpoints fixes $`A+B=b-a`$.

The pair $`(A,A+B)`$ is jointly Gaussian, with covariance $`t-s`$ and second-coordinate variance $`T-s`$. Gaussian conditioning therefore gives

```math
E[A\mid A+B=b-a]=\frac{t-s}{T-s}(b-a),
```

and conditional variance

```math
(t-s)-\frac{(t-s)^2}{T-s}
=\frac{(t-s)(T-t)}{T-s}.
```

Adding the known value $`W_s=a`$ yields

```math
W_t\mid(W_s=0,W_T=2)
\sim N\left(1,\frac{1}{4}\right).
```

The requested conditional probability is

```math
1-\Phi\left(\frac{1-(1)}{\sqrt{\frac{1}{4}}}\right)
\approx 0.500000.
```

Conditioning only on $`W_s`$ would leave variance $`t-s=\frac{1}{2}`$, which is larger. The future endpoint removes uncertainty about the total increment but leaves uncertainty about how that increment is split.

The conditional mean linearly interpolates the observed endpoints. The conditional variance tends to zero as the intermediate date approaches either endpoint. For endpoints that have probability zero as exact observations, this expression is understood as a regular Gaussian conditional distribution. It does not justify treating a trader who knows $`W_T`$ at time $`t\lt T`$ as adapted to the original Brownian filtration.

#### Problem SC-42: Brownian bridge at an intermediate date, case 2

For standard Brownian motion, take $`s=1\lt t=2\lt T=4`$ and condition on $`W_s=1`$, $`W_T=-1`$. Find the conditional law of $`W_t`$ and the probability that $`W_t\gt 0`$.

Solution. Set $`A=W_t-W_s`$ and $`B=W_T-W_t`$. Independent increments give independent normals with variances $`t-s`$ and $`T-t`$. Observing both endpoints fixes $`A+B=b-a`$.

The pair $`(A,A+B)`$ is jointly Gaussian, with covariance $`t-s`$ and second-coordinate variance $`T-s`$. Gaussian conditioning therefore gives

```math
E[A\mid A+B=b-a]=\frac{t-s}{T-s}(b-a),
```

and conditional variance

```math
(t-s)-\frac{(t-s)^2}{T-s}
=\frac{(t-s)(T-t)}{T-s}.
```

Adding the known value $`W_s=a`$ yields

```math
W_t\mid(W_s=1,W_T=-1)
\sim N\left(\frac{1}{3},\frac{2}{3}\right).
```

The requested conditional probability is

```math
1-\Phi\left(\frac{0-(\frac{1}{3})}{\sqrt{\frac{2}{3}}}\right)
\approx 0.658454.
```

Conditioning only on $`W_s`$ would leave variance $`t-s=1`$, which is larger. The future endpoint removes uncertainty about the total increment but leaves uncertainty about how that increment is split.

The conditional mean linearly interpolates the observed endpoints. The conditional variance tends to zero as the intermediate date approaches either endpoint. For endpoints that have probability zero as exact observations, this expression is understood as a regular Gaussian conditional distribution. It does not justify treating a trader who knows $`W_T`$ at time $`t\lt T`$ as adapted to the original Brownian filtration.

#### Problem SC-43: Brownian bridge at an intermediate date, case 3

For standard Brownian motion, take $`s=0\lt t=1\lt T=3`$ and condition on $`W_s=0`$, $`W_T=0`$. Find the conditional law of $`W_t`$ and the probability that $`W_t\gt 1`$.

Solution. Set $`A=W_t-W_s`$ and $`B=W_T-W_t`$. Independent increments give independent normals with variances $`t-s`$ and $`T-t`$. Observing both endpoints fixes $`A+B=b-a`$.

The pair $`(A,A+B)`$ is jointly Gaussian, with covariance $`t-s`$ and second-coordinate variance $`T-s`$. Gaussian conditioning therefore gives

```math
E[A\mid A+B=b-a]=\frac{t-s}{T-s}(b-a),
```

and conditional variance

```math
(t-s)-\frac{(t-s)^2}{T-s}
=\frac{(t-s)(T-t)}{T-s}.
```

Adding the known value $`W_s=a`$ yields

```math
W_t\mid(W_s=0,W_T=0)
\sim N\left(0,\frac{2}{3}\right).
```

The requested conditional probability is

```math
1-\Phi\left(\frac{1-(0)}{\sqrt{\frac{2}{3}}}\right)
\approx 0.110336.
```

Conditioning only on $`W_s`$ would leave variance $`t-s=1`$, which is larger. The future endpoint removes uncertainty about the total increment but leaves uncertainty about how that increment is split.

The conditional mean linearly interpolates the observed endpoints. The conditional variance tends to zero as the intermediate date approaches either endpoint. For endpoints that have probability zero as exact observations, this expression is understood as a regular Gaussian conditional distribution. It does not justify treating a trader who knows $`W_T`$ at time $`t\lt T`$ as adapted to the original Brownian filtration.

#### Problem SC-44: Brownian bridge at an intermediate date, case 4

For standard Brownian motion, take $`s=2\lt t=3\lt T=5`$ and condition on $`W_s=-2`$, $`W_T=1`$. Find the conditional law of $`W_t`$ and the probability that $`W_t\gt 0`$.

Solution. Set $`A=W_t-W_s`$ and $`B=W_T-W_t`$. Independent increments give independent normals with variances $`t-s`$ and $`T-t`$. Observing both endpoints fixes $`A+B=b-a`$.

The pair $`(A,A+B)`$ is jointly Gaussian, with covariance $`t-s`$ and second-coordinate variance $`T-s`$. Gaussian conditioning therefore gives

```math
E[A\mid A+B=b-a]=\frac{t-s}{T-s}(b-a),
```

and conditional variance

```math
(t-s)-\frac{(t-s)^2}{T-s}
=\frac{(t-s)(T-t)}{T-s}.
```

Adding the known value $`W_s=a`$ yields

```math
W_t\mid(W_s=-2,W_T=1)
\sim N\left(-1,\frac{2}{3}\right).
```

The requested conditional probability is

```math
1-\Phi\left(\frac{0-(-1)}{\sqrt{\frac{2}{3}}}\right)
\approx 0.110336.
```

Conditioning only on $`W_s`$ would leave variance $`t-s=1`$, which is larger. The future endpoint removes uncertainty about the total increment but leaves uncertainty about how that increment is split.

The conditional mean linearly interpolates the observed endpoints. The conditional variance tends to zero as the intermediate date approaches either endpoint. For endpoints that have probability zero as exact observations, this expression is understood as a regular Gaussian conditional distribution. It does not justify treating a trader who knows $`W_T`$ at time $`t\lt T`$ as adapted to the original Brownian filtration.

### 8.2 Polynomial stochastic integrals

Itô's correction changes the ordinary chain rule. For a polynomial it produces an identity expressing a stochastic integral through a terminal power and a time integral. The Itô isometry gives the integral's second moment without needing a full distribution.

#### Problem SC-45: integrating the Brownian power 1

For a standard Brownian motion started at zero, express $`I_T=\int_0^T W_t^{1}\,dW_t`$ in terms of ordinary quantities. Find its expectation and variance.

Solution. Apply Itô's formula to $`f(x)=x^m`$ with $`m=2`$. The first derivative is $`m x^{m-1}`$ and the second is $`m(m-1)x^{m-2}`$. Since $`d[W]_t=dt`$,

```math
d(W_t^m)=mW_t^{m-1}\,dW_t
+\frac{m(m-1)}2W_t^{m-2}\,dt.
```

Integrate from zero to $`T`$ and use $`W_0=0`$:

```math
I_T=\frac{W_T^m}m-\frac{m-1}2
\int_0^T W_t^{m-2}\,dt,
\qquad m=2.
```

The time-integral correction would be absent in the ordinary differentiable-path chain rule. Brownian quadratic variation is precisely what produces it.

For finite $`T`$, Gaussian moments ensure $`E\int_0^T W_t^{2m-2}dt\lt \infty`$. The stochastic integral is therefore a square-integrable martingale evaluated at $`T`$, so $`E[I_T]=0`$. The isometry gives

```math
E[I_T^2]=\int_0^T E[W_t^{2m-2}]\,dt.
```

For a centered normal with variance $`t`$, the even moment of order $`2m-2`$ is $`(2m-3)!!\,t^{m-1}`$. Here $`(2m-3)!!=1`$. Integration yields

```math
\mathop{\mathrm{Var}} (I_T)=
\frac{(2m-3)!!}mT^m=\frac{1}{2}T^{2}.
```

This also checks the integrability condition used above. A stochastic integral with an adapted integrand is not automatically Gaussian. In particular, the displayed polynomial identity can produce a visibly non-Gaussian terminal variable. Gaussianity follows directly for deterministic integrands, as in the next group, but generally fails for nonlinear Brownian integrands.

#### Problem SC-46: integrating the Brownian power 2

For a standard Brownian motion started at zero, express $`I_T=\int_0^T W_t^{2}\,dW_t`$ in terms of ordinary quantities. Find its expectation and variance.

Solution. Apply Itô's formula to $`f(x)=x^m`$ with $`m=3`$. The first derivative is $`m x^{m-1}`$ and the second is $`m(m-1)x^{m-2}`$. Since $`d[W]_t=dt`$,

```math
d(W_t^m)=mW_t^{m-1}\,dW_t
+\frac{m(m-1)}2W_t^{m-2}\,dt.
```

Integrate from zero to $`T`$ and use $`W_0=0`$:

```math
I_T=\frac{W_T^m}m-\frac{m-1}2
\int_0^T W_t^{m-2}\,dt,
\qquad m=3.
```

The time-integral correction would be absent in the ordinary differentiable-path chain rule. Brownian quadratic variation is precisely what produces it.

For finite $`T`$, Gaussian moments ensure $`E\int_0^T W_t^{2m-2}dt\lt \infty`$. The stochastic integral is therefore a square-integrable martingale evaluated at $`T`$, so $`E[I_T]=0`$. The isometry gives

```math
E[I_T^2]=\int_0^T E[W_t^{2m-2}]\,dt.
```

For a centered normal with variance $`t`$, the even moment of order $`2m-2`$ is $`(2m-3)!!\,t^{m-1}`$. Here $`(2m-3)!!=3`$. Integration yields

```math
\mathop{\mathrm{Var}} (I_T)=
\frac{(2m-3)!!}mT^m=1T^{3}.
```

This also checks the integrability condition used above. A stochastic integral with an adapted integrand is not automatically Gaussian. In particular, the displayed polynomial identity can produce a visibly non-Gaussian terminal variable. Gaussianity follows directly for deterministic integrands, as in the next group, but generally fails for nonlinear Brownian integrands.

#### Problem SC-47: integrating the Brownian power 3

For a standard Brownian motion started at zero, express $`I_T=\int_0^T W_t^{3}\,dW_t`$ in terms of ordinary quantities. Find its expectation and variance.

Solution. Apply Itô's formula to $`f(x)=x^m`$ with $`m=4`$. The first derivative is $`m x^{m-1}`$ and the second is $`m(m-1)x^{m-2}`$. Since $`d[W]_t=dt`$,

```math
d(W_t^m)=mW_t^{m-1}\,dW_t
+\frac{m(m-1)}2W_t^{m-2}\,dt.
```

Integrate from zero to $`T`$ and use $`W_0=0`$:

```math
I_T=\frac{W_T^m}m-\frac{m-1}2
\int_0^T W_t^{m-2}\,dt,
\qquad m=4.
```

The time-integral correction would be absent in the ordinary differentiable-path chain rule. Brownian quadratic variation is precisely what produces it.

For finite $`T`$, Gaussian moments ensure $`E\int_0^T W_t^{2m-2}dt\lt \infty`$. The stochastic integral is therefore a square-integrable martingale evaluated at $`T`$, so $`E[I_T]=0`$. The isometry gives

```math
E[I_T^2]=\int_0^T E[W_t^{2m-2}]\,dt.
```

For a centered normal with variance $`t`$, the even moment of order $`2m-2`$ is $`(2m-3)!!\,t^{m-1}`$. Here $`(2m-3)!!=15`$. Integration yields

```math
\mathop{\mathrm{Var}} (I_T)=
\frac{(2m-3)!!}mT^m=\frac{15}{4}T^{4}.
```

This also checks the integrability condition used above. A stochastic integral with an adapted integrand is not automatically Gaussian. In particular, the displayed polynomial identity can produce a visibly non-Gaussian terminal variable. Gaussianity follows directly for deterministic integrands, as in the next group, but generally fails for nonlinear Brownian integrands.

#### Problem SC-48: integrating the Brownian power 4

For a standard Brownian motion started at zero, express $`I_T=\int_0^T W_t^{4}\,dW_t`$ in terms of ordinary quantities. Find its expectation and variance.

Solution. Apply Itô's formula to $`f(x)=x^m`$ with $`m=5`$. The first derivative is $`m x^{m-1}`$ and the second is $`m(m-1)x^{m-2}`$. Since $`d[W]_t=dt`$,

```math
d(W_t^m)=mW_t^{m-1}\,dW_t
+\frac{m(m-1)}2W_t^{m-2}\,dt.
```

Integrate from zero to $`T`$ and use $`W_0=0`$:

```math
I_T=\frac{W_T^m}m-\frac{m-1}2
\int_0^T W_t^{m-2}\,dt,
\qquad m=5.
```

The time-integral correction would be absent in the ordinary differentiable-path chain rule. Brownian quadratic variation is precisely what produces it.

For finite $`T`$, Gaussian moments ensure $`E\int_0^T W_t^{2m-2}dt\lt \infty`$. The stochastic integral is therefore a square-integrable martingale evaluated at $`T`$, so $`E[I_T]=0`$. The isometry gives

```math
E[I_T^2]=\int_0^T E[W_t^{2m-2}]\,dt.
```

For a centered normal with variance $`t`$, the even moment of order $`2m-2`$ is $`(2m-3)!!\,t^{m-1}`$. Here $`(2m-3)!!=105`$. Integration yields

```math
\mathop{\mathrm{Var}} (I_T)=
\frac{(2m-3)!!}mT^m=21T^{5}.
```

This also checks the integrability condition used above. A stochastic integral with an adapted integrand is not automatically Gaussian. In particular, the displayed polynomial identity can produce a visibly non-Gaussian terminal variable. Gaussianity follows directly for deterministic integrands, as in the next group, but generally fails for nonlinear Brownian integrands.

### 8.3 Deterministic integrands and endpoint information

For deterministic square-integrable integrands, stochastic integration preserves Gaussianity. Joint laws with the Brownian endpoint follow from inner products of the deterministic integrands. An integrand can have a nonzero stochastic integral that is independent of the endpoint.

#### Problem SC-49: linear deterministic integrand, case 1

Let $`I=\int_0^T(a+bt)\,dW_t`$ with $`a=1`$, $`b=1`$, $`T=1`$. Find its law, its covariance with $`W_T`$, and its conditional law given $`W_T=w`$.

Solution. Approximate the deterministic integrand by step functions. Each step integral is a linear combination of jointly Gaussian increments, hence Gaussian. The $`L^2`$ limit remains centered Gaussian, with variance given by the isometry:

```math
v=\int_0^T(a+bt)^2dt
=a^2T+abT^2+\frac{b^2T^3}3=\frac{7}{3}.
```

Thus $`I\sim N(0,\frac{7}{3})`$. The Brownian endpoint is itself $`W_T=\int_0^T1\,dW_t`$. The bilinear isometry gives

```math
\mathop{\mathrm{Cov}} (I,W_T)=\int_0^T(a+bt)\,dt
=aT+\frac{bT^2}2=\frac{3}{2}.
```

Since $`\mathop{\mathrm{Var}} (W_T)=T`$, the Gaussian regression coefficient is $`c=\mathop{\mathrm{Cov}} (I,W_T)/T=\frac{3}{2}`$. The residual $`R=I-cW_T`$ is jointly Gaussian with $`W_T`$ and has zero covariance with it, so the two are independent.

Its variance is

```math
\mathop{\mathrm{Var}} (R)=v-\frac{\mathop{\mathrm{Cov}} (I,W_T)^2}T
=\frac{b^2T^3}{12}=\frac{1}{12}.
```

Consequently

```math
I\mid W_T=w\sim N\left((\frac{3}{2})w,\frac{1}{12}\right).
```

The stochastic integral depends on the full path, while $`W_T`$ records only its total increment. The residual represents path information not captured by that total. If the covariance is zero, joint Gaussianity proves independence even though both quantities are formed from the same Brownian path. For a constant integrand, the residual variance would instead be zero and the integral would be completely determined by the endpoint.

#### Problem SC-50: linear deterministic integrand, case 2

Let $`I=\int_0^T(a+bt)\,dW_t`$ with $`a=2`$, $`b=-1`$, $`T=2`$. Find its law, its covariance with $`W_T`$, and its conditional law given $`W_T=w`$.

Solution. Approximate the deterministic integrand by step functions. Each step integral is a linear combination of jointly Gaussian increments, hence Gaussian. The $`L^2`$ limit remains centered Gaussian, with variance given by the isometry:

```math
v=\int_0^T(a+bt)^2dt
=a^2T+abT^2+\frac{b^2T^3}3=\frac{8}{3}.
```

Thus $`I\sim N(0,\frac{8}{3})`$. The Brownian endpoint is itself $`W_T=\int_0^T1\,dW_t`$. The bilinear isometry gives

```math
\mathop{\mathrm{Cov}} (I,W_T)=\int_0^T(a+bt)\,dt
=aT+\frac{bT^2}2=2.
```

Since $`\mathop{\mathrm{Var}} (W_T)=T`$, the Gaussian regression coefficient is $`c=\mathop{\mathrm{Cov}} (I,W_T)/T=1`$. The residual $`R=I-cW_T`$ is jointly Gaussian with $`W_T`$ and has zero covariance with it, so the two are independent.

Its variance is

```math
\mathop{\mathrm{Var}} (R)=v-\frac{\mathop{\mathrm{Cov}} (I,W_T)^2}T
=\frac{b^2T^3}{12}=\frac{2}{3}.
```

Consequently

```math
I\mid W_T=w\sim N\left((1)w,\frac{2}{3}\right).
```

The stochastic integral depends on the full path, while $`W_T`$ records only its total increment. The residual represents path information not captured by that total. If the covariance is zero, joint Gaussianity proves independence even though both quantities are formed from the same Brownian path. For a constant integrand, the residual variance would instead be zero and the integral would be completely determined by the endpoint.

#### Problem SC-51: linear deterministic integrand, case 3

Let $`I=\int_0^T(a+bt)\,dW_t`$ with $`a=0`$, $`b=2`$, $`T=1`$. Find its law, its covariance with $`W_T`$, and its conditional law given $`W_T=w`$.

Solution. Approximate the deterministic integrand by step functions. Each step integral is a linear combination of jointly Gaussian increments, hence Gaussian. The $`L^2`$ limit remains centered Gaussian, with variance given by the isometry:

```math
v=\int_0^T(a+bt)^2dt
=a^2T+abT^2+\frac{b^2T^3}3=\frac{4}{3}.
```

Thus $`I\sim N(0,\frac{4}{3})`$. The Brownian endpoint is itself $`W_T=\int_0^T1\,dW_t`$. The bilinear isometry gives

```math
\mathop{\mathrm{Cov}} (I,W_T)=\int_0^T(a+bt)\,dt
=aT+\frac{bT^2}2=1.
```

Since $`\mathop{\mathrm{Var}} (W_T)=T`$, the Gaussian regression coefficient is $`c=\mathop{\mathrm{Cov}} (I,W_T)/T=1`$. The residual $`R=I-cW_T`$ is jointly Gaussian with $`W_T`$ and has zero covariance with it, so the two are independent.

Its variance is

```math
\mathop{\mathrm{Var}} (R)=v-\frac{\mathop{\mathrm{Cov}} (I,W_T)^2}T
=\frac{b^2T^3}{12}=\frac{1}{3}.
```

Consequently

```math
I\mid W_T=w\sim N\left((1)w,\frac{1}{3}\right).
```

The stochastic integral depends on the full path, while $`W_T`$ records only its total increment. The residual represents path information not captured by that total. If the covariance is zero, joint Gaussianity proves independence even though both quantities are formed from the same Brownian path. For a constant integrand, the residual variance would instead be zero and the integral would be completely determined by the endpoint.

#### Problem SC-52: linear deterministic integrand, case 4

Let $`I=\int_0^T(a+bt)\,dW_t`$ with $`a=1`$, $`b=-2`$, $`T=1`$. Find its law, its covariance with $`W_T`$, and its conditional law given $`W_T=w`$.

Solution. Approximate the deterministic integrand by step functions. Each step integral is a linear combination of jointly Gaussian increments, hence Gaussian. The $`L^2`$ limit remains centered Gaussian, with variance given by the isometry:

```math
v=\int_0^T(a+bt)^2dt
=a^2T+abT^2+\frac{b^2T^3}3=\frac{1}{3}.
```

Thus $`I\sim N(0,\frac{1}{3})`$. The Brownian endpoint is itself $`W_T=\int_0^T1\,dW_t`$. The bilinear isometry gives

```math
\mathop{\mathrm{Cov}} (I,W_T)=\int_0^T(a+bt)\,dt
=aT+\frac{bT^2}2=0.
```

Since $`\mathop{\mathrm{Var}} (W_T)=T`$, the Gaussian regression coefficient is $`c=\mathop{\mathrm{Cov}} (I,W_T)/T=0`$. The residual $`R=I-cW_T`$ is jointly Gaussian with $`W_T`$ and has zero covariance with it, so the two are independent.

Its variance is

```math
\mathop{\mathrm{Var}} (R)=v-\frac{\mathop{\mathrm{Cov}} (I,W_T)^2}T
=\frac{b^2T^3}{12}=\frac{1}{3}.
```

Consequently

```math
I\mid W_T=w\sim N\left((0)w,\frac{1}{3}\right).
```

The stochastic integral depends on the full path, while $`W_T`$ records only its total increment. The residual represents path information not captured by that total. If the covariance is zero, joint Gaussianity proves independence even though both quantities are formed from the same Brownian path. For a constant integrand, the residual variance would instead be zero and the integral would be completely determined by the endpoint.

### 8.4 Geometric Brownian moments and threshold events

Geometric Brownian motion is lognormal, so its mean, median, variance, and threshold probabilities can all be calculated explicitly. The drift of the logarithm differs from the drift of the level. Confusing them changes every probability involving a future price.

#### Problem SC-53: lognormal forecast, case 1

Let $`dS_t=\mu S_tdt+\sigma S_tdW_t`$ with $`S_0=100`$, $`\mu=\frac{1}{10}`$, $`\sigma=\frac{1}{5}`$. At $`T=1`$, find the mean, median, variance, and probability $`P(S_T\gt 110)`$.

Solution. Itô's formula for the logarithm gives

```math
d\log S_t=\left(\mu-\frac{\sigma^2}2\right)dt+\sigma dW_t.
```

Integrating yields $`S_T=S_0\exp((\mu-\sigma^2/2)T+\sigma W_T)`$. Hence $`\log(S_T/S_0)`$ is normal with mean $`(\mu-\sigma^2/2)T=\frac{2}{25}`$ and variance $`\sigma^2T=\frac{1}{25}`$.

Using $`E[e^{aW_T}]=e^{a^2T/2}`$ gives

```math
E[S_T]=S_0e^{\mu T}\approx 110.517092,\qquad
\mathop{\mathrm{median}} (S_T)=S_0e^{(\mu-\sigma^2/2)T}\approx 108.328707.
```

The mean exceeds the median because the positive tail contributes to the average. The second moment is $`S_0^2e^{(2\mu+\sigma^2)T}`$, so

```math
\mathop{\mathrm{Var}} (S_T)=S_0^2e^{2\mu T}(e^{\sigma^2T}-1)
\approx 498.463922.
```

For the threshold event, take logarithms, subtract the log mean, and divide by the positive log standard deviation:

```math
P(S_T\gt 110)=
1-\Phi\left(
\frac{\log(110/100)-(\mu-\sigma^2/2)T}{\sigma\sqrt T}
\right)
\approx 0.469490.
```

These are probabilities under the measure in which the supplied drift is $`\mu`$. For derivative pricing under a risk-neutral measure in this simple no-dividend model, the drift is replaced by the risk-free rate. A change of measure changes the drift and probabilities; it does not change the diffusion coefficient or turn the physical mean into a pricing formula.

#### Problem SC-54: lognormal forecast, case 2

Let $`dS_t=\mu S_tdt+\sigma S_tdW_t`$ with $`S_0=50`$, $`\mu=\frac{1}{20}`$, $`\sigma=\frac{3}{10}`$. At $`T=2`$, find the mean, median, variance, and probability $`P(S_T\gt 50)`$.

Solution. Itô's formula for the logarithm gives

```math
d\log S_t=\left(\mu-\frac{\sigma^2}2\right)dt+\sigma dW_t.
```

Integrating yields $`S_T=S_0\exp((\mu-\sigma^2/2)T+\sigma W_T)`$. Hence $`\log(S_T/S_0)`$ is normal with mean $`(\mu-\sigma^2/2)T=\frac{1}{100}`$ and variance $`\sigma^2T=\frac{9}{50}`$.

Using $`E[e^{aW_T}]=e^{a^2T/2}`$ gives

```math
E[S_T]=S_0e^{\mu T}\approx 55.258546,\qquad
\mathop{\mathrm{median}} (S_T)=S_0e^{(\mu-\sigma^2/2)T}\approx 50.502508.
```

The mean exceeds the median because the positive tail contributes to the average. The second moment is $`S_0^2e^{(2\mu+\sigma^2)T}`$, so

```math
\mathop{\mathrm{Var}} (S_T)=S_0^2e^{2\mu T}(e^{\sigma^2T}-1)
\approx 602.204578.
```

For the threshold event, take logarithms, subtract the log mean, and divide by the positive log standard deviation:

```math
P(S_T\gt 50)=
1-\Phi\left(
\frac{\log(50/50)-(\mu-\sigma^2/2)T}{\sigma\sqrt T}
\right)
\approx 0.509402.
```

These are probabilities under the measure in which the supplied drift is $`\mu`$. For derivative pricing under a risk-neutral measure in this simple no-dividend model, the drift is replaced by the risk-free rate. A change of measure changes the drift and probabilities; it does not change the diffusion coefficient or turn the physical mean into a pricing formula.

#### Problem SC-55: lognormal forecast, case 3

Let $`dS_t=\mu S_tdt+\sigma S_tdW_t`$ with $`S_0=80`$, $`\mu=-\frac{1}{20}`$, $`\sigma=\frac{1}{4}`$. At $`T=1`$, find the mean, median, variance, and probability $`P(S_T\gt 70)`$.

Solution. Itô's formula for the logarithm gives

```math
d\log S_t=\left(\mu-\frac{\sigma^2}2\right)dt+\sigma dW_t.
```

Integrating yields $`S_T=S_0\exp((\mu-\sigma^2/2)T+\sigma W_T)`$. Hence $`\log(S_T/S_0)`$ is normal with mean $`(\mu-\sigma^2/2)T=-\frac{13}{160}`$ and variance $`\sigma^2T=\frac{1}{16}`$.

Using $`E[e^{aW_T}]=e^{a^2T/2}`$ gives

```math
E[S_T]=S_0e^{\mu T}\approx 76.098354,\qquad
\mathop{\mathrm{median}} (S_T)=S_0e^{(\mu-\sigma^2/2)T}\approx 73.757054.
```

The mean exceeds the median because the positive tail contributes to the average. The second moment is $`S_0^2e^{(2\mu+\sigma^2)T}`$, so

```math
\mathop{\mathrm{Var}} (S_T)=S_0^2e^{2\mu T}(e^{\sigma^2T}-1)
\approx 373.484798.
```

For the threshold event, take logarithms, subtract the log mean, and divide by the positive log standard deviation:

```math
P(S_T\gt 70)=
1-\Phi\left(
\frac{\log(70/80)-(\mu-\sigma^2/2)T}{\sigma\sqrt T}
\right)
\approx 0.582825.
```

These are probabilities under the measure in which the supplied drift is $`\mu`$. For derivative pricing under a risk-neutral measure in this simple no-dividend model, the drift is replaced by the risk-free rate. A change of measure changes the drift and probabilities; it does not change the diffusion coefficient or turn the physical mean into a pricing formula.

#### Problem SC-56: lognormal forecast, case 4

Let $`dS_t=\mu S_tdt+\sigma S_tdW_t`$ with $`S_0=120`$, $`\mu=0`$, $`\sigma=\frac{2}{5}`$. At $`T=\frac{1}{2}`$, find the mean, median, variance, and probability $`P(S_T\gt 130)`$.

Solution. Itô's formula for the logarithm gives

```math
d\log S_t=\left(\mu-\frac{\sigma^2}2\right)dt+\sigma dW_t.
```

Integrating yields $`S_T=S_0\exp((\mu-\sigma^2/2)T+\sigma W_T)`$. Hence $`\log(S_T/S_0)`$ is normal with mean $`(\mu-\sigma^2/2)T=-\frac{1}{25}`$ and variance $`\sigma^2T=\frac{2}{25}`$.

Using $`E[e^{aW_T}]=e^{a^2T/2}`$ gives

```math
E[S_T]=S_0e^{\mu T}\approx 120.000000,\qquad
\mathop{\mathrm{median}} (S_T)=S_0e^{(\mu-\sigma^2/2)T}\approx 115.294733.
```

The mean exceeds the median because the positive tail contributes to the average. The second moment is $`S_0^2e^{(2\mu+\sigma^2)T}`$, so

```math
\mathop{\mathrm{Var}} (S_T)=S_0^2e^{2\mu T}(e^{\sigma^2T}-1)
\approx 1199.333775.
```

For the threshold event, take logarithms, subtract the log mean, and divide by the positive log standard deviation:

```math
P(S_T\gt 130)=
1-\Phi\left(
\frac{\log(130/120)-(\mu-\sigma^2/2)T}{\sigma\sqrt T}
\right)
\approx 0.335632.
```

These are probabilities under the measure in which the supplied drift is $`\mu`$. For derivative pricing under a risk-neutral measure in this simple no-dividend model, the drift is replaced by the risk-free rate. A change of measure changes the drift and probabilities; it does not change the diffusion coefficient or turn the physical mean into a pricing formula.

### 8.5 An Ornstein–Uhlenbeck transition

The integrating factor solves the mean-reverting linear SDE exactly. It separates a decaying initial-condition effect from a Gaussian innovation. Its transition law permits an exact simulation step without Euler discretization error.

#### Problem SC-57: exact mean-reverting transition, case 1

Solve $`dX_t=\kappa(\theta-X_t)dt+\sigma dW_t`$ from $`X_0=5`$, with $`\kappa=1`$, $`\theta=2`$, $`\sigma=1`$. Find its law at $`T=1`$ and its stationary variance.

Solution. Multiply the centered process by $`e^{\kappa t}`$. The ordinary factor has finite variation, so there is no cross quadratic-variation term:

```math
d[e^{\kappa t}(X_t-\theta)]
=\sigma e^{\kappa t}dW_t.
```

After integration and division by the factor,

```math
X_T=\theta+(x-\theta)e^{-\kappa T}
+\sigma\int_0^T e^{-\kappa(T-s)}dW_s.
```

The last term is centered Gaussian because its integrand is deterministic. Its variance follows from the isometry:

```math
\sigma^2\int_0^T e^{-2\kappa(T-s)}ds
=\frac{\sigma^2}{2\kappa}(1-e^{-2\kappa T})
\approx 0.432332.
```

Thus $`X_T\sim N(3.103638,0.432332)`$. The conditional mean is pulled from the initial value toward $`\theta`$, while the variance grows from zero toward $`\sigma^2/(2\kappa)=\frac{1}{2}`$.

If the initial state instead has the stationary law $`N(\theta,\sigma^2/(2\kappa))`$ independently of future increments, this variance is preserved at every date. For lag $`h\ge0`$, the stationary covariance is $`[\sigma^2/(2\kappa)]e^{-\kappa h}`$.

An exact simulation step over length $`h`$ is

```math
X_{t+h}=\theta+(X_t-\theta)e^{-\kappa h}
+\sqrt{\frac{\sigma^2}{2\kappa}(1-e^{-2\kappa h})}\,Z,
\qquad Z\sim N(0,1),
```

using an independent $`Z`$ at each step. Euler simulation replaces the exponential decay and exact variance by first-order approximations. It can be useful for more general equations, but is unnecessary for this transition and can be unstable if the step is too large.

#### Problem SC-58: exact mean-reverting transition, case 2

Solve $`dX_t=\kappa(\theta-X_t)dt+\sigma dW_t`$ from $`X_0=-1`$, with $`\kappa=\frac{1}{2}`$, $`\theta=0`$, $`\sigma=2`$. Find its law at $`T=2`$ and its stationary variance.

Solution. Multiply the centered process by $`e^{\kappa t}`$. The ordinary factor has finite variation, so there is no cross quadratic-variation term:

```math
d[e^{\kappa t}(X_t-\theta)]
=\sigma e^{\kappa t}dW_t.
```

After integration and division by the factor,

```math
X_T=\theta+(x-\theta)e^{-\kappa T}
+\sigma\int_0^T e^{-\kappa(T-s)}dW_s.
```

The last term is centered Gaussian because its integrand is deterministic. Its variance follows from the isometry:

```math
\sigma^2\int_0^T e^{-2\kappa(T-s)}ds
=\frac{\sigma^2}{2\kappa}(1-e^{-2\kappa T})
\approx 3.458659.
```

Thus $`X_T\sim N(-0.367879,3.458659)`$. The conditional mean is pulled from the initial value toward $`\theta`$, while the variance grows from zero toward $`\sigma^2/(2\kappa)=4`$.

If the initial state instead has the stationary law $`N(\theta,\sigma^2/(2\kappa))`$ independently of future increments, this variance is preserved at every date. For lag $`h\ge0`$, the stationary covariance is $`[\sigma^2/(2\kappa)]e^{-\kappa h}`$.

An exact simulation step over length $`h`$ is

```math
X_{t+h}=\theta+(X_t-\theta)e^{-\kappa h}
+\sqrt{\frac{\sigma^2}{2\kappa}(1-e^{-2\kappa h})}\,Z,
\qquad Z\sim N(0,1),
```

using an independent $`Z`$ at each step. Euler simulation replaces the exponential decay and exact variance by first-order approximations. It can be useful for more general equations, but is unnecessary for this transition and can be unstable if the step is too large.

#### Problem SC-59: exact mean-reverting transition, case 3

Solve $`dX_t=\kappa(\theta-X_t)dt+\sigma dW_t`$ from $`X_0=0`$, with $`\kappa=2`$, $`\theta=1`$, $`\sigma=\frac{1}{2}`$. Find its law at $`T=\frac{1}{2}`$ and its stationary variance.

Solution. Multiply the centered process by $`e^{\kappa t}`$. The ordinary factor has finite variation, so there is no cross quadratic-variation term:

```math
d[e^{\kappa t}(X_t-\theta)]
=\sigma e^{\kappa t}dW_t.
```

After integration and division by the factor,

```math
X_T=\theta+(x-\theta)e^{-\kappa T}
+\sigma\int_0^T e^{-\kappa(T-s)}dW_s.
```

The last term is centered Gaussian because its integrand is deterministic. Its variance follows from the isometry:

```math
\sigma^2\int_0^T e^{-2\kappa(T-s)}ds
=\frac{\sigma^2}{2\kappa}(1-e^{-2\kappa T})
\approx 0.054042.
```

Thus $`X_T\sim N(0.632121,0.054042)`$. The conditional mean is pulled from the initial value toward $`\theta`$, while the variance grows from zero toward $`\sigma^2/(2\kappa)=\frac{1}{16}`$.

If the initial state instead has the stationary law $`N(\theta,\sigma^2/(2\kappa))`$ independently of future increments, this variance is preserved at every date. For lag $`h\ge0`$, the stationary covariance is $`[\sigma^2/(2\kappa)]e^{-\kappa h}`$.

An exact simulation step over length $`h`$ is

```math
X_{t+h}=\theta+(X_t-\theta)e^{-\kappa h}
+\sqrt{\frac{\sigma^2}{2\kappa}(1-e^{-2\kappa h})}\,Z,
\qquad Z\sim N(0,1),
```

using an independent $`Z`$ at each step. Euler simulation replaces the exponential decay and exact variance by first-order approximations. It can be useful for more general equations, but is unnecessary for this transition and can be unstable if the step is too large.

#### Problem SC-60: exact mean-reverting transition, case 4

Solve $`dX_t=\kappa(\theta-X_t)dt+\sigma dW_t`$ from $`X_0=3`$, with $`\kappa=\frac{1}{4}`$, $`\theta=3`$, $`\sigma=1`$. Find its law at $`T=4`$ and its stationary variance.

Solution. Multiply the centered process by $`e^{\kappa t}`$. The ordinary factor has finite variation, so there is no cross quadratic-variation term:

```math
d[e^{\kappa t}(X_t-\theta)]
=\sigma e^{\kappa t}dW_t.
```

After integration and division by the factor,

```math
X_T=\theta+(x-\theta)e^{-\kappa T}
+\sigma\int_0^T e^{-\kappa(T-s)}dW_s.
```

The last term is centered Gaussian because its integrand is deterministic. Its variance follows from the isometry:

```math
\sigma^2\int_0^T e^{-2\kappa(T-s)}ds
=\frac{\sigma^2}{2\kappa}(1-e^{-2\kappa T})
\approx 1.729329.
```

Thus $`X_T\sim N(3.000000,1.729329)`$. The conditional mean is pulled from the initial value toward $`\theta`$, while the variance grows from zero toward $`\sigma^2/(2\kappa)=2`$.

If the initial state instead has the stationary law $`N(\theta,\sigma^2/(2\kappa))`$ independently of future increments, this variance is preserved at every date. For lag $`h\ge0`$, the stationary covariance is $`[\sigma^2/(2\kappa)]e^{-\kappa h}`$.

An exact simulation step over length $`h`$ is

```math
X_{t+h}=\theta+(X_t-\theta)e^{-\kappa h}
+\sqrt{\frac{\sigma^2}{2\kappa}(1-e^{-2\kappa h})}\,Z,
\qquad Z\sim N(0,1),
```

using an independent $`Z`$ at each step. Euler simulation replaces the exponential decay and exact variance by first-order approximations. It can be useful for more general equations, but is unnecessary for this transition and can be unstable if the step is too large.

### 8.6 A constant change of Brownian drift

For a constant market price of risk, the density is an exponential of a normal variable. Its expectation, the shifted Brownian mean, and the transformed price drift can all be checked directly. These checks clarify the sign convention in Girsanov's theorem.

#### Problem SC-61: risk-neutral drift conversion, case 1

Under $`P`$, let $`dS_t=\mu S_tdt+\sigma S_tdW_t^P`$ with $`\mu=\frac{1}{10}`$, $`r=\frac{1}{20}`$, $`\sigma=\frac{1}{5}`$. Construct the equivalent measure on $`[0,T]`$ for $`T=1`$ that gives price drift $`r`$.

Solution. Define the constant $`\theta=(\mu-r)/\sigma=\frac{1}{4}`$ and the density process

```math
Z_t=\exp\left(-\theta W_t^P-\frac12\theta^2t\right).
```

The sign in front of the Brownian term determines the sign of the drift shift. The normal moment-generating function gives $`E_P[Z_t]=e^{\theta^2t/2}e^{-\theta^2t/2}=1`$. The constant-integrand Novikov expression is $`\exp(\theta^2T/2)\lt \infty`$, so $`Z`$ is a true martingale on this finite interval.

Set $`dQ/dP=Z_T`$. Girsanov's theorem gives a $`Q`$ Brownian motion

```math
W_t^Q=W_t^P+\theta t.
```

As a direct terminal check, exponential tilting changes $`W_T^P`$ from mean zero to mean $`-\theta T=-\frac{1}{4}`$ and leaves its variance $`T`$ unchanged. Therefore $`W_T^Q`$ has mean zero and variance $`T`$ under $`Q`$.

Substitute $`dW_t^P=dW_t^Q-\theta dt`$ into the price equation:

```math
dS_t=(\mu-\sigma\theta)S_tdt+\sigma S_tdW_t^Q
=rS_tdt+\sigma S_tdW_t^Q.
```

The discounted price is a true martingale in this constant-coefficient model because its exponential-normal factor has conditional expectation one. Thus $`E_Q[e^{-rT}S_T]=S_0`$.

The measures are equivalent on the stated finite horizon because the density is strictly positive. A finite-horizon construction should not be silently extended to equivalence on an infinite-horizon path sigma-field. If $`\theta`$ were a random process, the martingale property of the density would also require an appropriate condition; positivity of a stochastic exponential alone would not suffice.

#### Problem SC-62: risk-neutral drift conversion, case 2

Under $`P`$, let $`dS_t=\mu S_tdt+\sigma S_tdW_t^P`$ with $`\mu=\frac{1}{5}`$, $`r=\frac{1}{25}`$, $`\sigma=\frac{2}{5}`$. Construct the equivalent measure on $`[0,T]`$ for $`T=2`$ that gives price drift $`r`$.

Solution. Define the constant $`\theta=(\mu-r)/\sigma=\frac{2}{5}`$ and the density process

```math
Z_t=\exp\left(-\theta W_t^P-\frac12\theta^2t\right).
```

The sign in front of the Brownian term determines the sign of the drift shift. The normal moment-generating function gives $`E_P[Z_t]=e^{\theta^2t/2}e^{-\theta^2t/2}=1`$. The constant-integrand Novikov expression is $`\exp(\theta^2T/2)\lt \infty`$, so $`Z`$ is a true martingale on this finite interval.

Set $`dQ/dP=Z_T`$. Girsanov's theorem gives a $`Q`$ Brownian motion

```math
W_t^Q=W_t^P+\theta t.
```

As a direct terminal check, exponential tilting changes $`W_T^P`$ from mean zero to mean $`-\theta T=-\frac{4}{5}`$ and leaves its variance $`T`$ unchanged. Therefore $`W_T^Q`$ has mean zero and variance $`T`$ under $`Q`$.

Substitute $`dW_t^P=dW_t^Q-\theta dt`$ into the price equation:

```math
dS_t=(\mu-\sigma\theta)S_tdt+\sigma S_tdW_t^Q
=rS_tdt+\sigma S_tdW_t^Q.
```

The discounted price is a true martingale in this constant-coefficient model because its exponential-normal factor has conditional expectation one. Thus $`E_Q[e^{-rT}S_T]=S_0`$.

The measures are equivalent on the stated finite horizon because the density is strictly positive. A finite-horizon construction should not be silently extended to equivalence on an infinite-horizon path sigma-field. If $`\theta`$ were a random process, the martingale property of the density would also require an appropriate condition; positivity of a stochastic exponential alone would not suffice.

#### Problem SC-63: risk-neutral drift conversion, case 3

Under $`P`$, let $`dS_t=\mu S_tdt+\sigma S_tdW_t^P`$ with $`\mu=-\frac{1}{20}`$, $`r=\frac{1}{20}`$, $`\sigma=\frac{1}{4}`$. Construct the equivalent measure on $`[0,T]`$ for $`T=1`$ that gives price drift $`r`$.

Solution. Define the constant $`\theta=(\mu-r)/\sigma=-\frac{2}{5}`$ and the density process

```math
Z_t=\exp\left(-\theta W_t^P-\frac12\theta^2t\right).
```

The sign in front of the Brownian term determines the sign of the drift shift. The normal moment-generating function gives $`E_P[Z_t]=e^{\theta^2t/2}e^{-\theta^2t/2}=1`$. The constant-integrand Novikov expression is $`\exp(\theta^2T/2)\lt \infty`$, so $`Z`$ is a true martingale on this finite interval.

Set $`dQ/dP=Z_T`$. Girsanov's theorem gives a $`Q`$ Brownian motion

```math
W_t^Q=W_t^P+\theta t.
```

As a direct terminal check, exponential tilting changes $`W_T^P`$ from mean zero to mean $`-\theta T=\frac{2}{5}`$ and leaves its variance $`T`$ unchanged. Therefore $`W_T^Q`$ has mean zero and variance $`T`$ under $`Q`$.

Substitute $`dW_t^P=dW_t^Q-\theta dt`$ into the price equation:

```math
dS_t=(\mu-\sigma\theta)S_tdt+\sigma S_tdW_t^Q
=rS_tdt+\sigma S_tdW_t^Q.
```

The discounted price is a true martingale in this constant-coefficient model because its exponential-normal factor has conditional expectation one. Thus $`E_Q[e^{-rT}S_T]=S_0`$.

The measures are equivalent on the stated finite horizon because the density is strictly positive. A finite-horizon construction should not be silently extended to equivalence on an infinite-horizon path sigma-field. If $`\theta`$ were a random process, the martingale property of the density would also require an appropriate condition; positivity of a stochastic exponential alone would not suffice.

#### Problem SC-64: risk-neutral drift conversion, case 4

Under $`P`$, let $`dS_t=\mu S_tdt+\sigma S_tdW_t^P`$ with $`\mu=\frac{3}{25}`$, $`r=\frac{1}{50}`$, $`\sigma=\frac{1}{2}`$. Construct the equivalent measure on $`[0,T]`$ for $`T=\frac{1}{2}`$ that gives price drift $`r`$.

Solution. Define the constant $`\theta=(\mu-r)/\sigma=\frac{1}{5}`$ and the density process

```math
Z_t=\exp\left(-\theta W_t^P-\frac12\theta^2t\right).
```

The sign in front of the Brownian term determines the sign of the drift shift. The normal moment-generating function gives $`E_P[Z_t]=e^{\theta^2t/2}e^{-\theta^2t/2}=1`$. The constant-integrand Novikov expression is $`\exp(\theta^2T/2)\lt \infty`$, so $`Z`$ is a true martingale on this finite interval.

Set $`dQ/dP=Z_T`$. Girsanov's theorem gives a $`Q`$ Brownian motion

```math
W_t^Q=W_t^P+\theta t.
```

As a direct terminal check, exponential tilting changes $`W_T^P`$ from mean zero to mean $`-\theta T=-\frac{1}{10}`$ and leaves its variance $`T`$ unchanged. Therefore $`W_T^Q`$ has mean zero and variance $`T`$ under $`Q`$.

Substitute $`dW_t^P=dW_t^Q-\theta dt`$ into the price equation:

```math
dS_t=(\mu-\sigma\theta)S_tdt+\sigma S_tdW_t^Q
=rS_tdt+\sigma S_tdW_t^Q.
```

The discounted price is a true martingale in this constant-coefficient model because its exponential-normal factor has conditional expectation one. Thus $`E_Q[e^{-rT}S_T]=S_0`$.

The measures are equivalent on the stated finite horizon because the density is strictly positive. A finite-horizon construction should not be silently extended to equivalence on an infinite-horizon path sigma-field. If $`\theta`$ were a random process, the martingale property of the density would also require an appropriate condition; positivity of a stochastic exponential alone would not suffice.

### 8.7 Pricing and hedging a quadratic terminal payoff

A polynomial payoff gives a useful pricing check because its expectation is explicit and its derivatives can be inserted directly into the pricing PDE. The hedge generally differs from the number of shares suggested by expanding the payoff at the current stock price.

#### Problem SC-65: quadratic payoff valuation, case 1

In a no-dividend Black–Scholes model under $`Q`$, price $`(S_T-K)^2`$. Take $`K=100`$, $`r=\frac{1}{20}`$, $`\sigma=\frac{1}{5}`$, time to maturity $`\tau=1`$, and current $`S_t=100`$. Find price, delta, and gamma and check the PDE.

Solution. Write $`v(t,s)=e^{-r\tau}E_Q[(S_T-K)^2\mid S_t=s]`$ with $`\tau=T-t`$. Under $`Q`$,

```math
E_Q[S_T\mid S_t=s]=se^{r\tau},\qquad
E_Q[S_T^2\mid S_t=s]=s^2e^{(2r+\sigma^2)\tau}.
```

Expanding the payoff and discounting each expectation gives

```math
v(t,s)=s^2e^{(r+\sigma^2)\tau}-2Ks+K^2e^{-r\tau}.
```

At the supplied state, this is approximately $`454.037082`$. Differentiating with respect to the current stock price gives

```math
\Delta=v_s=2se^{(r+\sigma^2)\tau}-2K
\approx 18.834857,\qquad
\Gamma=v_{ss}=2e^{(r+\sigma^2)\tau}\approx 2.188349.
```

To check the pricing PDE, put $`A=e^{(r+\sigma^2)\tau}`$ and $`B=e^{-r\tau}`$. Since $`\partial_t=-\partial_\tau`$,

```math
v_t=-(r+\sigma^2)s^2A+rK^2B.
```

Substitute into $`v_t+rsv_s+\tfrac12\sigma^2s^2v_{ss}-rv`$. The coefficients of $`s^2A`$, $`Ks`$, and $`K^2B`$ each cancel separately, giving zero. At maturity, $`A=B=1`$ and the expression reduces to $`(s-K)^2`$.

Holding $`\Delta`$ shares and placing $`v-\Delta s`$ in the money-market account gives the replicating portfolio under the usual self-financing strategy. The payoff has quadratic growth, but geometric Brownian motion has finite moments of all fixed orders on a finite horizon, so its discounted expectation is finite. This integrability check supports the expectation and PDE calculations.

#### Problem SC-66: quadratic payoff valuation, case 2

In a no-dividend Black–Scholes model under $`Q`$, price $`(S_T-K)^2`$. Take $`K=50`$, $`r=0`$, $`\sigma=\frac{3}{10}`$, time to maturity $`\tau=2`$, and current $`S_t=60`$. Find price, delta, and gamma and check the PDE.

Solution. Write $`v(t,s)=e^{-r\tau}E_Q[(S_T-K)^2\mid S_t=s]`$ with $`\tau=T-t`$. Under $`Q`$,

```math
E_Q[S_T\mid S_t=s]=se^{r\tau},\qquad
E_Q[S_T^2\mid S_t=s]=s^2e^{(2r+\sigma^2)\tau}.
```

Expanding the payoff and discounting each expectation gives

```math
v(t,s)=s^2e^{(r+\sigma^2)\tau}-2Ks+K^2e^{-r\tau}.
```

At the supplied state, this is approximately $`809.982507`$. Differentiating with respect to the current stock price gives

```math
\Delta=v_s=2se^{(r+\sigma^2)\tau}-2K
\approx 43.666084,\qquad
\Gamma=v_{ss}=2e^{(r+\sigma^2)\tau}\approx 2.394435.
```

To check the pricing PDE, put $`A=e^{(r+\sigma^2)\tau}`$ and $`B=e^{-r\tau}`$. Since $`\partial_t=-\partial_\tau`$,

```math
v_t=-(r+\sigma^2)s^2A+rK^2B.
```

Substitute into $`v_t+rsv_s+\tfrac12\sigma^2s^2v_{ss}-rv`$. The coefficients of $`s^2A`$, $`Ks`$, and $`K^2B`$ each cancel separately, giving zero. At maturity, $`A=B=1`$ and the expression reduces to $`(s-K)^2`$.

Holding $`\Delta`$ shares and placing $`v-\Delta s`$ in the money-market account gives the replicating portfolio under the usual self-financing strategy. The payoff has quadratic growth, but geometric Brownian motion has finite moments of all fixed orders on a finite horizon, so its discounted expectation is finite. This integrability check supports the expectation and PDE calculations.

#### Problem SC-67: quadratic payoff valuation, case 3

In a no-dividend Black–Scholes model under $`Q`$, price $`(S_T-K)^2`$. Take $`K=80`$, $`r=\frac{1}{25}`$, $`\sigma=\frac{1}{4}`$, time to maturity $`\tau=\frac{1}{2}`$, and current $`S_t=70`$. Find price, delta, and gamma and check the PDE.

Solution. Write $`v(t,s)=e^{-r\tau}E_Q[(S_T-K)^2\mid S_t=s]`$ with $`\tau=T-t`$. Under $`Q`$,

```math
E_Q[S_T\mid S_t=s]=se^{r\tau},\qquad
E_Q[S_T^2\mid S_t=s]=s^2e^{(2r+\sigma^2)\tau}.
```

Expanding the payoff and discounting each expectation gives

```math
v(t,s)=s^2e^{(r+\sigma^2)\tau}-2Ks+K^2e^{-r\tau}.
```

At the supplied state, this is approximately $`230.942943`$. Differentiating with respect to the current stock price gives

```math
\Delta=v_s=2se^{(r+\sigma^2)\tau}-2K
\approx -12.637959,\qquad
\Gamma=v_{ss}=2e^{(r+\sigma^2)\tau}\approx 2.105172.
```

To check the pricing PDE, put $`A=e^{(r+\sigma^2)\tau}`$ and $`B=e^{-r\tau}`$. Since $`\partial_t=-\partial_\tau`$,

```math
v_t=-(r+\sigma^2)s^2A+rK^2B.
```

Substitute into $`v_t+rsv_s+\tfrac12\sigma^2s^2v_{ss}-rv`$. The coefficients of $`s^2A`$, $`Ks`$, and $`K^2B`$ each cancel separately, giving zero. At maturity, $`A=B=1`$ and the expression reduces to $`(s-K)^2`$.

Holding $`\Delta`$ shares and placing $`v-\Delta s`$ in the money-market account gives the replicating portfolio under the usual self-financing strategy. The payoff has quadratic growth, but geometric Brownian motion has finite moments of all fixed orders on a finite horizon, so its discounted expectation is finite. This integrability check supports the expectation and PDE calculations.

#### Problem SC-68: quadratic payoff valuation, case 4

In a no-dividend Black–Scholes model under $`Q`$, price $`(S_T-K)^2`$. Take $`K=120`$, $`r=\frac{1}{10}`$, $`\sigma=\frac{2}{5}`$, time to maturity $`\tau=1`$, and current $`S_t=130`$. Find price, delta, and gamma and check the PDE.

Solution. Write $`v(t,s)=e^{-r\tau}E_Q[(S_T-K)^2\mid S_t=s]`$ with $`\tau=T-t`$. Under $`Q`$,

```math
E_Q[S_T\mid S_t=s]=se^{r\tau},\qquad
E_Q[S_T^2\mid S_t=s]=s^2e^{(2r+\sigma^2)\tau}.
```

Expanding the payoff and discounting each expectation gives

```math
v(t,s)=s^2e^{(r+\sigma^2)\tau}-2Ks+K^2e^{-r\tau}.
```

At the supplied state, this is approximately $`3747.777284`$. Differentiating with respect to the current stock price gives

```math
\Delta=v_s=2se^{(r+\sigma^2)\tau}-2K
\approx 97.201823,\qquad
\Gamma=v_{ss}=2e^{(r+\sigma^2)\tau}\approx 2.593860.
```

To check the pricing PDE, put $`A=e^{(r+\sigma^2)\tau}`$ and $`B=e^{-r\tau}`$. Since $`\partial_t=-\partial_\tau`$,

```math
v_t=-(r+\sigma^2)s^2A+rK^2B.
```

Substitute into $`v_t+rsv_s+\tfrac12\sigma^2s^2v_{ss}-rv`$. The coefficients of $`s^2A`$, $`Ks`$, and $`K^2B`$ each cancel separately, giving zero. At maturity, $`A=B=1`$ and the expression reduces to $`(s-K)^2`$.

Holding $`\Delta`$ shares and placing $`v-\Delta s`$ in the money-market account gives the replicating portfolio under the usual self-financing strategy. The payoff has quadratic growth, but geometric Brownian motion has finite moments of all fixed orders on a finite horizon, so its discounted expectation is finite. This integrability check supports the expectation and PDE calculations.

<!-- further-distinct-problems-SC -->

### 8.8 Further martingale and SDE problems

These additional problems are constructed for these notes.

#### Problem SC-69: Brownian exit probabilities and duration

Let $`\tau`$ be the first exit of standard Brownian motion started at zero from $`(-a,b)`$, with $`a,b\gt 0`$. Find the probability of upper exit and the expected exit time, justifying the stopping limits.

Solution. For the bounded time $`\tau\wedge n`$, the martingales $`W_t`$ and $`W_t^2-t`$ give $`E[W_{\tau\wedge n}]=0`$ and $`E[\tau\wedge n]=E[W_{\tau\wedge n}^2]\le\max(a,b)^2`$. Monotone convergence proves finite expected exit time, hence almost-sure exit.

Continuity gives $`W_\tau=-a`$ or $`b`$. The stopped positions are uniformly bounded, so dominated convergence passes both position expectations to the exit. If $`p=P(W_\tau=b)`$, then

```math
0=bp-a(1-p),\qquad p=\frac a{a+b}.
```

The square identity gives

```math
E[\tau]=b^2p+a^2(1-p)=ab.
```

For $`a=1,b=2`$, the upper-exit probability is $`1/3`$ and expected duration is two. The unequal boundary distances explain why the more distant upper endpoint is less likely.

The bounded stopped-position argument is what allows an unbounded exit time to be used in the final identity. Without it, writing $`E[W_\tau]=0`$ from a martingale property alone would not be a proof. On an unbounded interval, the same domination can fail, and a one-sided Brownian hitting time has infinite expectation.

#### Problem SC-70: the integral of a Brownian path

Let $`A_T=\int_0^T W_s\,ds`$. Find its joint Gaussian law with $`W_T`$ and its conditional law given $`W_T=w`$.

Solution. Use the product rule on $`tW_t`$:

```math
d(tW_t)=W_tdt+t\,dW_t.
```

Integration gives $`A_T=TW_T-\int_0^T s\,dW_s=\int_0^T(T-s)dW_s`$. The integrand is deterministic, so $`A_T`$ and $`W_T=\int_0^T1\,dW_s`$ are jointly Gaussian.

The isometry and its bilinear version yield

```math
\mathop{\mathrm{Var}} (A_T)=\int_0^T(T-s)^2ds=\frac{T^3}3,\qquad
\mathop{\mathrm{Cov}} (A_T,W_T)=\int_0^T(T-s)ds=\frac{T^2}2.
```

Both means are zero, and $`\mathop{\mathrm{Var}} (W_T)=T`$. Therefore Gaussian conditioning gives

```math
A_T\mid W_T=w\sim N\left(\frac T2w,\frac{T^3}{12}\right).
```

The residual $`A_T-(T/2)W_T`$ is independent of $`W_T`$, with variance $`T^3/12`$. Knowing the final Brownian value therefore determines the average linear interpolation component but leaves random fluctuations around that interpolation.

For $`T=2`$ and endpoint one, the conditional mean is one and variance $`2/3`$. The time average $`A_T/T`$ has conditional variance $`T/12`$. It cannot be recovered from the terminal value alone, even though both are Gaussian functionals of the same path.

#### Problem SC-71: integrated squared Brownian motion

For $`B_T=\int_0^T W_s^2ds`$, compute its mean and variance and its conditional expectation given $`W_T=w`$.

Solution. Nonnegativity allows expectation to pass through the integral. Since $`E[W_s^2]=s`$, $`E[B_T]=T^2/2`$.

For $`s\le t`$, write $`W_t=W_s+(W_t-W_s)`$. The increment is centered independent normal with variance $`t-s`$. Expanding gives

```math
E[W_s^2W_t^2]=E[W_s^4]+E[W_s^2](t-s)
=3s^2+s(t-s)=st+2s^2.
```

Thus $`\mathop{\mathrm{Cov}} (W_s^2,W_t^2)=2\min(s,t)^2`$. Finite fourth moments justify integrating this covariance:

```math
\mathop{\mathrm{Var}} (B_T)
=2\int_0^T\int_0^T\min(s,t)^2dsdt
=4\int_0^T\int_0^t s^2dsdt
=\frac{T^4}3.
```

Conditional on the endpoint, the bridge value at $`s`$ has mean $`sw/T`$ and variance $`s(T-s)/T`$. Its conditional second moment is the variance plus the squared mean. Integrating yields

```math
E[B_T\mid W_T=w]
=\int_0^T\left[\frac{s(T-s)}T+\frac{s^2w^2}{T^2}\right]ds
=\frac{T^2}6+\frac{Tw^2}3.
```

Averaging this expression over $`W_T`$, whose second moment is $`T`$, recovers $`T^2/2`$. This checks the conditioning calculation. The conditional expectation is not obtained by simply integrating the squared conditional mean; the bridge's conditional variance contributes the extra $`T^2/6`$.

#### Problem SC-72: a correlated product martingale

Let standard Brownian motions $`X,Y`$ have constant correlation $`\rho\in[-1,1]`$. Prove that $`M_t=X_tY_t-\rho t`$ is a true martingale on every finite horizon and compute its variance.

Solution. The product rule gives

```math
d(X_tY_t)=X_tdY_t+Y_tdX_t+\rho\,dt.
```

Subtracting $`\rho t`$ leaves a stochastic integral. Its integrand has finite expected squared norm on finite horizons, so it is square integrable and is a true martingale, not merely a formal zero-drift process.

For a moment calculation, represent $`Y_t=\rho X_t+\sqrt{1-\rho^2}Z_t`$ with independent Brownian $`Z`$. Expanding and using the centered fourth normal moment gives

```math
E[X_t^2Y_t^2]
=3\rho^2t^2+(1-\rho^2)t^2
=(1+2\rho^2)t^2.
```

The product's mean is $`\rho t`$, so $`\mathop{\mathrm{Var}} (M_t)=(1+\rho^2)t^2`$.

There is an isometry check. The quadratic variation of the martingale is

```math
[M]_t=\int_0^t(X_s^2+Y_s^2+2\rho X_sY_s)ds.
```

Taking expectations gives $`\int_0^t(2s+2\rho^2s)ds=(1+\rho^2)t^2`$, matching the direct result. At $`\rho=0`$ the uncorrected product is a martingale. At $`\rho=1`$ it reduces to $`X_t^2-t`$. The endpoint cases remain valid despite singular joint covariance.

#### Problem SC-73: a bounded random market price of risk

On a finite horizon let $`\theta_t=k1_{\lbrace W_t\gt 0\rbrace }`$ for a real constant $`k`$. Determine whether $`Z_t=\exp(-\int_0^t\theta_s\,dW_s-\tfrac12\int_0^t\theta_s^2ds)`$ defines a valid density. State the Brownian motion under the resulting measure and give a second-moment bound for the density.

Solution. The bounded progressively measurable integrand is admissible for stochastic integration. Its quadratic integral satisfies $`\int_0^T\theta_s^2ds\le k^2T`$. Hence the Novikov expectation is bounded by $`e^{k^2T/2}`$, which is finite.

The stochastic exponential is consequently a positive true martingale with expectation one. Defining $`dQ/dP=Z_T`$ gives an equivalent measure on this horizon. Girsanov gives

```math
W_t^Q=W_t+\int_0^t\theta_sds
=W_t+k\int_0^t1_{\lbrace W_s\gt 0\rbrace }ds.
```

The original $`W`$ therefore has drift $`-k1_{\lbrace W_t\gt 0\rbrace }`$ when expressed using $`W^Q`$. The shifted process, rather than the original one, is the Brownian motion under $`Q`$.

For the second moment, factor

```math
Z_T^2=
\exp\left(-2\int_0^T\theta_s dW_s-2\int_0^T\theta_s^2ds\right)
\exp\left(\int_0^T\theta_s^2ds\right).
```

The first factor is the stochastic exponential for the bounded integrand $`-2\theta`$ and has expectation one. The second is bounded by $`e^{k^2T}`$, giving $`E[Z_T^2]\le e^{k^2T}`$.

The random-integrand log density is not generally a normal variable. The constant-integrand moment-generating-function shortcut cannot simply be reused. Boundedness supplied a direct martingale criterion and moment bound without requiring an explicit distribution.

#### Problem SC-74: a digital claim under two probability measures

Let $`S_0=K=100`$, maturity one, volatility $`0.2`$, physical drift $`0.10`$, and risk-free rate $`0.02`$. Compare the physical probability of $`S_T\gt K`$ with the arbitrage price of a one-unit cash payment on that event. Derive its delta for positive remaining maturity.

Solution. Under the physical measure, the log return has mean $`0.10-0.2^2/2=0.08`$ and standard deviation $`0.2`$. Therefore $`P(S_T\gt K)=\Phi(0.4)\approx0.655422`$.

Under the risk-neutral measure, the log mean is $`0.02-0.2^2/2=0`$. The event probability is one half. Discounting gives price $`e^{-0.02}/2\approx0.490099`$. The physical event probability is a forecasting quantity; the replication price uses the risk-neutral law in this complete model.

For current stock $`s`$ and remaining time $`\tau\gt 0`$, the price is

```math
v(t,s)=e^{-r\tau}\Phi(d_2),\qquad
d_2=\frac{\log(s/K)+(r-\sigma^2/2)\tau}{\sigma\sqrt\tau}.
```

Differentiating gives

```math
\Delta=e^{-r\tau}\frac{\phi(d_2)}{s\sigma\sqrt\tau}.
```

This hedge becomes concentrated around the strike as maturity approaches. The discontinuous terminal payoff does not require a discontinuous price at earlier times: the lognormal transition smooths the conditional expectation.

Multiplying the physical probability by the money-market discount factor would give a different number. That calculation is an expected discounted payment under physical probabilities and does not satisfy the replication argument used to derive this model's arbitrage price.

#### Problem SC-75: integrating a geometric process against its driver

Let $`dS_t=\mu S_tdt+\sigma S_tdW_t`$ with $`S_0=s\gt 0`$ and $`\sigma\ne0`$. Express $`I_T=\int_0^TS_t\,dW_t`$ using $`S_T`$ and a time integral. Compute its variance and covariance with $`W_T`$.

Solution. Integrating the defining SDE gives $`S_T-s=\mu\int_0^TS_tdt+\sigma I_T`$. Therefore

```math
I_T=\frac{S_T-s-\mu\int_0^TS_tdt}{\sigma}.
```

Finite lognormal moments imply square integrability of the integrand on every finite horizon, so $`E[I_T]=0`$. The isometry gives

```math
\mathop{\mathrm{Var}} (I_T)
=\int_0^TE[S_t^2]dt
=s^2\int_0^Te^{(2\mu+\sigma^2)t}dt.
```

This is $`s^2(e^{(2\mu+\sigma^2)T}-1)/(2\mu+\sigma^2)`$ when the denominator is nonzero, and $`s^2T`$ when it is zero. The latter is the continuous limit of the former.

The bilinear isometry with $`W_T=\int_0^T1\,dW_t`$ gives

```math
\mathop{\mathrm{Cov}} (I_T,W_T)
=E\int_0^TS_tdt
=s\int_0^Te^{\mu t}dt.
```

Thus the covariance is $`s(e^{\mu T}-1)/\mu`$ for nonzero $`\mu`$, and $`sT`$ for zero $`\mu`$.

In the special case $`\mu=0`$ and $`\sigma\gt 0`$, $`I_T=(S_T-s)/\sigma`$ is bounded below but nonconstant. It cannot be Gaussian. This provides an explicit random-integrand example whose mean and variance are obtained exactly by stochastic integration rules even though its terminal law is a shifted lognormal.

#### Problem SC-76: an affine short-rate bond and its PDE

Under the pricing measure, let $`dr_t=\kappa(\theta-r_t)dt+\sigma dW_t`$ with $`\kappa\gt 0`$. Derive an exponential-affine zero-coupon bond price and verify its coefficient equations.

Solution. The price is $`P(\tau,r)=E_r[\exp(-\int_0^\tau r_sds)]`$. Seek $`P=\exp(A(\tau)-B(\tau)r)`$ with $`A(0)=B(0)=0`$. The time-to-maturity PDE is

```math
P_\tau=\kappa(\theta-r)P_r+\tfrac12\sigma^2P_{rr}-rP.
```

The derivatives satisfy $`P_r=-BP`$ and $`P_{rr}=B^2P`$. Comparing the coefficient of $`r`$ and the constant term gives

```math
B'=1-\kappa B,\qquad
A'=-\kappa\theta B+\tfrac12\sigma^2B^2.
```

The first equation gives $`B=(1-e^{-\kappa\tau})/\kappa`$. Since $`\kappa\int_0^\tau B(s)ds=\tau-B(\tau)`$, integration of the second yields

```math
A=-\theta(\tau-B)+\frac{\sigma^2}{2\kappa^2}
\left[\tau-2B+\frac{1-e^{-2\kappa\tau}}{2\kappa}\right].
```

The integrated short rate is Gaussian with mean $`rB+\theta(\tau-B)`$ and variance equal to twice the positive variance term in $`A`$. Its exponential-normal expectation reproduces the same formula, independently checking the PDE calculation.

The Gaussian short rate can be negative. Consequently a bond price above one is possible in some states and does not contradict this model's discounting rule. Excluding negative rates would require a different state model, not clipping values inside the derived Gaussian expectation.

#### Problem SC-77: an inhomogeneous Feynman–Kac calculation

For Brownian motion $`X_s=x+W_s`$ started at $`x`$, compute $`v(t,x)=E_x[\int_0^{T-t}X_s^2ds]`$. Verify its terminal condition and PDE, and explain the martingale used in the representation.

Solution. For each elapsed time $`s`$, $`E_x[X_s^2]=x^2+s`$. Nonnegativity allows integration of expectations, giving

```math
v(t,x)=x^2\tau+\frac{\tau^2}{2},\qquad \tau=T-t.
```

The derivatives are $`v_t=-x^2-\tau`$ and $`v_{xx}=2\tau`$. Therefore

```math
v_t+\tfrac12v_{xx}+x^2=0,\qquad v(T,x)=0.
```

Along a Brownian state process, Itô's formula gives $`dv(s,X_s)=-X_s^2ds+v_x(s,X_s)dW_s`$. Adding the accumulated running payment $`\int_t^sX_u^2du`$ cancels the drift.

The resulting stochastic integrand is $`2(T-s)X_s`$. Its expected squared integral is finite on the bounded horizon because $`E[X_s^2]`$ is finite and grows only linearly in elapsed time. Thus the corrected process is a true martingale, and taking expectations at maturity is justified.

This example keeps the source term's sign visible. The running payment enters the backward PDE with a plus sign and the representation with a positive integral. A homogeneous terminal-value equation would give zero for the zero terminal payoff and would miss the accumulated payments entirely.

#### Problem SC-78: one-sided passage with a negative drift

Let $`X_t=W_t+bt`$ and let $`\tau_a`$ be its first hit of $`a\gt 0`$, with any real drift $`b`$. Find the eventual hitting probability and, when $`b\lt 0`$, the conditional mean given a hit.

Solution. For $`\lambda\gt 0`$, the same bounded exponential-stopping argument as in SC-28 gives

```math
E[e^{-\lambda\tau_a};\tau_a\lt \infty]
=\exp\left[a\left(b-\sqrt{b^2+2\lambda}\right)\right].
```

Alternatively, this solves the stopped generator equation $`u''/2+bu'=\lambda u`$ below the barrier with boundary value one at the barrier and decay at the lower infinite endpoint for positive $`\lambda`$.

Let $`\lambda`$ decrease to zero. The exponential on the finite-hitting event increases to one, so monotone convergence gives

```math
P(\tau_a\lt \infty)=
\begin{cases}
1,&b\ge0,\\[0pt]
e^{2ab},&b\lt 0.
\end{cases}
```

For negative drift, divide the Laplace transform by this probability. The conditional transform becomes $`\exp[a(\vert b\vert -\sqrt{b^2+2\lambda})]`$. Its derivative at zero gives

```math
E[\tau_a\mid\tau_a\lt \infty]=\frac a{\vert b\vert }.
```

For $`a=2,b=-1`$, the hitting probability is $`e^{-4}`$ and the conditional expected duration is two. The unconditional expectation of the extended hitting time is infinite because the time equals infinity with positive probability. It must not be confused with either the conditional mean or $`E[\tau_a1_{\lbrace \tau_a\lt \infty\rbrace }]`$, the latter being $`e^{2ab}a/\vert b\vert`$.

At zero drift, hitting is almost sure but the mean is still infinite. The transform's nonfinite derivative at zero agrees with the first-passage density's heavy tail.
