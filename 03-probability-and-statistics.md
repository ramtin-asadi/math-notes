# Probability and Statistics

Probability describes a model before observations arrive. Statistics uses observations to learn about an unknown feature of a model. These notes develop both, starting with events and ending with inference for dependent observations.

Sources:

- Joe Blitzstein and Jessica Hwang, [Harvard Stat 110 and its probability text](https://stat110.hsites.harvard.edu/). The [strategic practice collection](https://stat110.hsites.harvard.edu/strategic-practice-problems) has its own section numbering within each practice sheet.
- Philippe Rigollet, [MIT 18.650, Statistics for Applications, Fall 2016](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/), especially Problem Sets 1–6 and 10.

The source links identify the selected problems; their statements are paraphrased. Additional problems and worked examples are constructed for these notes. Throughout, iid means independent and identically distributed. The notation $N(\mu,\sigma^2)$ uses variance as its second argument, and $\operatorname{Exp}(\lambda)$ uses a rate, with mean $1/\lambda$.

Contents:

1. [Probability foundations](#1-probability-foundations)
2. [Random variables and distributions](#2-random-variables-and-distributions)
3. [Conditional expectation](#3-conditional-expectation)
4. [Multivariate probability](#4-multivariate-probability)
5. [Convergence and asymptotic theory](#5-convergence-and-asymptotic-theory)
6. [Statistical estimation](#6-statistical-estimation)
7. [Statistical inference](#7-statistical-inference)
8. [Resampling and dependence](#8-resampling-and-dependence)
9. [Further worked problems](#9-further-worked-problems)

## 1. Probability foundations

### 1.1 Theory

A probability model is a triple $(\Omega,\mathcal F,P)$. The sample space $\Omega$ contains possible outcomes. The sigma-algebra $\mathcal F$ contains the events to which probabilities are assigned. It contains $\Omega$, is closed under complements, and is closed under countable unions. These rules also imply closure under countable intersections by De Morgan's law.

The probability measure satisfies $P(A)\ge0$, $P(\Omega)=1$, and countable additivity: for pairwise disjoint events $A_i$, $P(\bigcup_iA_i)=\sum_iP(A_i)$. Disjointness means the events cannot occur together. Independence means learning one occurred does not change the other's probability. They are different properties: two disjoint events with positive probabilities cannot be independent.

For $P(B)>0$, conditional probability is

$$
P(A\mid B)=\frac{P(A\cap B)}{P(B)}.
$$

The denominator restricts attention to outcomes in $B$ and rescales their total probability to one. For fixed $B$, the function $A\mapsto P(A\mid B)$ is itself a probability measure. If $P(B)=0$, this ratio is undefined. Conditional distributions given continuous random variables require a density or conditional-expectation construction instead.

Events $A,B$ are independent when $P(A\cap B)=P(A)P(B)$. For several events, mutual independence requires the corresponding product equality for every subcollection. Pairwise equalities do not control triple intersections.

### 1.2 Derivations and proofs

#### Total probability and Bayes' theorem

Let $B_1,\ldots,B_k$ partition $\Omega$, with positive probabilities. The sets $A\cap B_i$ are disjoint and their union is $A$, so

$$
P(A)=\sum_iP(A\mid B_i)P(B_i).
$$

When $P(A)>0$, substitute $P(A\cap B_j)=P(A\mid B_j)P(B_j)$ into the definition:

$$
P(B_j\mid A)=
\frac{P(A\mid B_j)P(B_j)}
{\sum_iP(A\mid B_i)P(B_i)}.
$$

The likelihood $P(A\mid B_j)$ and prior weight $P(B_j)$ both affect the posterior probability. A high probability of observing evidence under a hypothesis does not by itself imply a high posterior probability for that hypothesis.

#### Inclusion-exclusion

Split $A\cup B$ into $A$ and $B\setminus A$. Since $P(B)=P(B\setminus A)+P(A\cap B)$,

$$
P(A\cup B)=P(A)+P(B)-P(A\cap B).
$$

For three events, add their individual probabilities, subtract the three pairwise intersections, and add the triple intersection back. An outcome in all three events was first counted three times, then removed three times, and must finally be counted once.

#### Events, random variables, and information

A probability space consists of outcomes $\Omega$, a sigma-algebra $\mathcal F$ of events, and a probability measure $P$. A sigma-algebra contains the whole space, is closed under complements, and is closed under countable unions. Countability matters because limits of sequences of events can then be described as events.

A random variable is a measurable map from outcomes to real values. Measurability means that events such as $\{X\le x\}$ belong to $\mathcal F$. The CDF assigns probabilities to these events. In a finite outcome space with every subset measurable, this condition is automatic; in continuous models it is part of the definition.

Countable additivity concerns disjoint events. For arbitrary events, $P(A\cup B)=P(A)+P(B)-P(A\cap B)$. Omitting the intersection double-counts outcomes in both events. Independence instead says $P(A\cap B)=P(A)P(B)$ and is an extra property of the model.

Conditional probability $P(A\mid B)=P(A\cap B)/P(B)$ needs $P(B)>0$. Conditioning on the exact value of a continuous variable cannot usually use that ratio because $P(X=x)=0$. A conditional density or regular conditional distribution supplies the appropriate meaning. Its values on conditioning events of probability zero can be chosen in different versions without changing the probabilistic identities.

#### Independence beyond a pair of events

Mutual independence requires factorization for every finite subcollection, not only every pair. If $X,Y$ are independent fair bits and $Z$ is their parity, each pair is independent, but knowing any two determines the third. Pairwise checks therefore do not justify factoring a three-way probability.

For independent random variables, measurable functions of separate variables remain independent. For example, $X^2$ and $e^Y$ are independent if $X,Y$ are independent. Functions of the same variable need not be independent merely because their covariance vanishes.

#### Example PR-E21: repeated evidence needs a model

Suppose a latent event has prior probability $p$, and two tests each have sensitivity $s$ and false-positive probability $f$. Find the posterior after two positives when the tests are conditionally independent given the latent event and its complement.

Solution. Under the stated conditional independence, $P(++\mid D)=s^2$ and $P(++\mid D^c)=f^2$. Bayes gives

$$
P(D\mid++)=\frac{ps^2}{ps^2+(1-p)f^2}.
$$

If the second test simply repeats the first result, these squared probabilities are wrong: the pair provides only one observation. Marginally identical sensitivities do not determine the joint test behavior. The conditional-independence assumption supplies the missing information.

### 1.3 Worked examples

#### Example PR-E01: verifying a finite model

Take $\Omega=\{1,2,3\}$, $\mathcal F=2^\Omega$, and assign masses $1/2,1/3,1/6$.

Solution. Every subset is measurable. The masses are nonnegative and sum to one. Defining $P(A)$ as the sum of masses in $A$ gives additivity because disjoint events contain no outcome twice. Thus this is a probability model. The collection $\{\varnothing,\{1\},\Omega\}$ would not be a sigma-algebra, since it omits the complement $\{2,3\}$.

#### Example PR-E02: a diagnostic calculation

Suppose $P(D)=0.01$, $P(+\mid D)=0.9$, and $P(+\mid D^c)=0.05$.

Solution. Total probability gives $P(+)=0.9(0.01)+0.05(0.99)=0.0585$. Hence $P(D\mid+)=0.009/0.0585=2/13$, about $0.154$. Among ten thousand modeled cases, the expected counts are ninety positive cases from $D$ and 495 from $D^c$. The posterior fraction is $90/(90+495)$, the same result expressed as counts.

### 1.4 Selected exercises

#### Problem PR-01: independence from itself

Source: [Harvard Strategic Practice 2, Section 2, Problem 1](https://stat110.hsites.harvard.edu/sites/g/files/omnuum10111/files/stat110/files/strategic_practice_and_homework_2.pdf).

Problem. Determine exactly when an event can be independent of itself.

Solution. Independence would require $P(A)=P(A\cap A)=P(A)^2$. Writing $p=P(A)$ gives $p(1-p)=0$, so $p=0$ or $p=1$. Both values also satisfy the equality, so the condition is necessary and sufficient. Probability-zero events can be nonempty; the result concerns probability rather than set emptiness.

#### Problem PR-02: complements of independent events

Source: [Harvard Strategic Practice 2, Section 2, Problem 2](https://stat110.hsites.harvard.edu/sites/g/files/omnuum10111/files/stat110/files/strategic_practice_and_homework_2.pdf).

Problem. Prove that complementing either or both of two independent events preserves independence.

Solution. If $P(A\cap B)=P(A)P(B)$, inclusion-exclusion gives

$$
P(A^c\cap B^c)=1-P(A)-P(B)+P(A)P(B)
=(1-P(A))(1-P(B)).
$$

This equals $P(A^c)P(B^c)$, proving independence. Also $P(A\cap B^c)=P(A)-P(A\cap B)=P(A)(1-P(B))$, so complementing only one event preserves independence too.

#### Problem PR-03: pairwise independence without mutual independence

Source: [Harvard Strategic Practice 2, Section 2, Problem 3](https://stat110.hsites.harvard.edu/sites/g/files/omnuum10111/files/stat110/files/strategic_practice_and_homework_2.pdf).

Problem. Construct three events that are pairwise independent but not mutually independent, verifying every required probability.

Solution. Toss two fair independent coins. Let $A$ mean the first is heads, $B$ mean the second is heads, and $C$ mean they agree. Each has probability $1/2$. Each pairwise intersection consists of the single outcome HH, with probability $1/4$, so every pair is independent. Yet $P(A\cap B\cap C)=1/4$, whereas $P(A)P(B)P(C)=1/8$. Knowing the first two events determines whether the third occurs.

#### Problem PR-04: a phrase-based spam probability

Source: [Harvard Strategic Practice 2, Section 3, Problem 2](https://stat110.hsites.harvard.edu/sites/g/files/omnuum10111/files/stat110/files/strategic_practice_and_homework_2.pdf).

Problem. Suppose 80% of messages are spam, a phrase occurs in 10% of spam messages and 1% of other messages. Find the probability of spam given the phrase and determine whether the phrase and spam events are independent.

Solution. Let $S$ mean spam and $F$ the specified phrase. The data are $P(S)=0.8$, $P(F\mid S)=0.1$, and $P(F\mid S^c)=0.01$. Thus

$$
P(S\mid F)=\frac{0.1(0.8)}{0.1(0.8)+0.01(0.2)}
=\frac{40}{41}.
$$

The events are not independent, since $P(F\mid S)=0.1$ differs from $P(F)=0.082$. Both the phrase frequency and the high initial prevalence contribute to the posterior.

### 1.5 Additional solved problems

#### Problem PR-05: a coarse sigma-algebra

Constructed problem. Let $\Omega=\{1,2,3,4\}$ and $\mathcal F=\{\varnothing,\{1,2\},\{3,4\},\Omega\}$. Verify its sigma-algebra properties and describe its information.

Solution. The complement of either two-point set is the other, and the complements of $\Omega$ and $\varnothing$ exchange them. Every union of sets from this collection is one of the four listed sets, including countable unions because only four distinct sets are available. It is a sigma-algebra. Observing an event in this collection can distinguish the pair $\{1,2\}$ from $\{3,4\}$ but cannot distinguish one from two. A random variable measurable with respect to it must be constant within each pair.

#### Problem PR-06: two-stage updating

Constructed problem. Show that conditioning first on $B$ and then on $C$ agrees with conditioning on $B\cap C$, when $P(B\cap C)>0$.

Solution. Define $P_B(A)=P(A\mid B)$. Then

$$
P_B(A\mid C)=\frac{P_B(A\cap C)}{P_B(C)}
=\frac{P(A\cap B\cap C)/P(B)}{P(B\cap C)/P(B)}
=P(A\mid B\cap C).
$$

The second update must use the already-conditioned probability model in both numerator and denominator. Multiplying unrelated marginal likelihoods would not generally produce this identity.

## 2. Random variables and distributions

### 2.1 Theory

A real random variable $X$ is a measurable function from outcomes to real numbers. Its cumulative distribution function is $F_X(x)=P(X\le x)$. Every CDF is nondecreasing, right-continuous, and has limits zero and one at the two infinities.

A discrete variable has a probability mass function $p(x)=P(X=x)$, whose nonnegative values sum to one. An absolutely continuous variable has a density $f$, with $F(x)=\int_{-\infty}^xf(t)\,dt$. A density value is not a point probability: for a continuous variable, $P(X=x)=0$ even if $f(x)>0$. Densities may exceed one, provided their area is one.

For integrable $X$, expectation is a weighted average. In the discrete and continuous cases,

$$
E[g(X)]=\sum_xg(x)p(x),\qquad
E[g(X)]=\int g(x)f(x)\,dx.
$$

An expectation is finite when $E|g(X)|<\infty$. The variance is $\operatorname{Var}(X)=E[(X-EX)^2]=E[X^2]-(EX)^2$, provided $E[X^2]<\infty$.

A joint density describes several variables at once. Marginalizing means integrating out unwanted coordinates. The integration region must respect the support; rectangular limits cannot be assumed for a triangular or curved support.

### 2.2 Derivations and proofs

#### Transforming a continuous variable

For a strictly increasing differentiable map $Y=g(X)$ with inverse $h$,

$$
F_Y(y)=P(X\le h(y))=F_X(h(y)).
$$

Differentiation gives $f_Y(y)=f_X(h(y))|h'(y)|$. The absolute value handles decreasing maps as well. For a non-one-to-one map, split the support into monotone branches and add their contributions. A transformation can produce atoms if it maps a positive-probability region to a single value.

#### The tail integral formula

For a nonnegative $X$, pointwise $X=\int_0^\infty\mathbf1_{\{X>t\}}dt$. Interchanging two nonnegative integrals gives

$$
EX=\int_0^\infty P(X>t)\,dt.
$$

The interchange is valid even if the answer is infinite, by Tonelli's theorem. For an exponential variable, integrating $e^{-\lambda t}$ immediately gives $EX=1/\lambda$.

#### Expectations and the order of integration

For a discrete variable, expectation is a probability-weighted sum. For a variable with density, it is the corresponding integral. A finite real expectation requires integrability of the absolute value, $E|X|<\infty$. A cancellation between infinite positive and negative contributions does not define a finite expectation.

For a nonnegative random variable, the tail formula is

$$
E[X]=\int_0^\infty P(X>t)\,dt.
$$

To see this, write $X(\omega)=\int_0^\infty1_{\{X(\omega)>t\}}dt$ and interchange the two nonnegative integrals. Tonelli's theorem permits this even when the answer is infinite. For signed integrands, a Fubini interchange needs an integrability condition.

Variance is $E[(X-EX)^2]$ and requires a finite second moment. Expanding gives $E[X^2]-(EX)^2$ when these terms are finite. The expression should not be used as an undefined subtraction of infinities.

#### Nonmonotone transformations need every branch

For $Y=g(X)$, a monotone differentiable branch gives density $f_X(g^{-1}(y))|(g^{-1})'(y)|$. If several distinct input points map to the same output, sum the contributions of all branches. A missing branch generally produces a density integrating to less than one.

The CDF method avoids this bookkeeping when the event $\{g(X)\le y\}$ has a simple form. It also handles transformations with atoms, for which a density alone would be incomplete.

#### Example PR-E22: squaring a standard normal

Let $Z\sim N(0,1)$ and $Y=Z^2$. Find its density and expectation.

Solution. For $y>0$, the two inverse points are $\sqrt y$ and $-\sqrt y$. Each has inverse-derivative magnitude $1/(2\sqrt y)$. Symmetry gives

$$
f_Y(y)=\frac{\phi(\sqrt y)+\phi(-\sqrt y)}{2\sqrt y}
=\frac{e^{-y/2}}{\sqrt{2\pi y}},\qquad y>0.
$$

It is zero for negative $y$ and has no atom at zero. This is the chi-squared density with one degree of freedom. Its mean is $E[Z^2]=1$. Its variance is $E[Z^4]-1=3-1=2$.

Integrating only the positive inverse branch would yield total mass one half. The singular density near zero is integrable; the change of variables concentrates many small normal values into a short interval of squared values.

### 2.3 Worked examples

#### Example PR-E03: a monotone transformation

Let $U\sim\operatorname{Unif}(0,1)$ and $Y=-\log U$.

Solution. For $y\ge0$, $P(Y\le y)=P(U\ge e^{-y})=1-e^{-y}$. Thus $f_Y(y)=e^{-y}$ on $y>0$: $Y$ is exponential with rate one. The decreasing transformation reverses the inequality before the probability is evaluated.

#### Example PR-E04: squaring a symmetric variable

Let $X\sim\operatorname{Unif}(-1,1)$ and $Y=X^2$.

Solution. For $0\le y\le1$, $P(Y\le y)=P(-\sqrt y\le X\le\sqrt y)=\sqrt y$. Therefore $f_Y(y)=1/(2\sqrt y)$ on $(0,1)$. The density diverges near zero but integrates to one. Using only the positive inverse branch would lose half the probability.

### 2.4 Selected exercises

#### Problem PR-07: marginals inside a ball

Source: [Harvard Strategic Practice 7, Section 1, Problem 1(a–c)](https://stat110.hsites.harvard.edu/sites/g/files/omnuum10111/files/stat110/files/strategic_practice_and_homework_7.pdf).

Problem. Choose a point uniformly from the three-dimensional unit ball. Find its joint density, the density of its first two coordinates, and the marginal density of one coordinate.

Solution. The unit ball in three dimensions has volume $4\pi/3$, so the joint density is $3/(4\pi)$ inside it and zero elsewhere. For fixed $(x,y)$ inside the unit disk, the allowed $z$ interval has length $2\sqrt{1-x^2-y^2}$. Hence

$$
f_{X,Y}(x,y)=\frac{3}{2\pi}\sqrt{1-x^2-y^2},
\qquad x^2+y^2<1.
$$

For fixed $x$, integrate over $-\sqrt{1-x^2}<y<\sqrt{1-x^2}$. The integral of $\sqrt{a^2-y^2}$ over $[-a,a]$ is the area of a semicircle, $\pi a^2/2$. Therefore

$$
f_X(x)=\frac34(1-x^2),\qquad -1<x<1.
$$

Its integral is $\tfrac34(2-2/3)=1$. A coordinate of a point uniform in a ball is not uniform on its coordinate interval; central slices contain more volume.

#### Problem PR-08: distance between two uniforms

Source: [Harvard Strategic Practice 7, Section 1, Problem 2](https://stat110.hsites.harvard.edu/sites/g/files/omnuum10111/files/stat110/files/strategic_practice_and_homework_7.pdf).

Problem. For independent uniform variables on $[0,1]$, find the mean and standard deviation of their absolute difference.

Solution. For independent uniforms $X,Y$, symmetry across the diagonal of the unit square gives

$$
E|X-Y|=2\int_0^1\int_0^y(y-x)\,dx\,dy
=2\int_0^1y^2/2\,dy=\frac13.
$$

Also $E|X-Y|^2=E(X-Y)^2=\operatorname{Var}(X)+\operatorname{Var}(Y)=1/6$, since their means agree and they are independent. The variance of the distance is $1/6-1/9=1/18$ and its standard deviation is $1/(3\sqrt2)$.

#### Problem PR-09: the smallest and largest of three uniforms

Source: [Harvard Strategic Practice 7, Section 1, Problem 3(a,b)](https://stat110.hsites.harvard.edu/sites/g/files/omnuum10111/files/stat110/files/strategic_practice_and_homework_7.pdf).

Problem. For three independent uniforms on $[0,1]$, derive the maximum's distribution, the joint law of minimum and maximum, and the conditional density of the maximum given the minimum.

Solution. For $M=\max(U_1,U_2,U_3)$, independence gives $F_M(m)=m^3$ on $[0,1]$, hence $f_M(m)=3m^2$. For $0\le l\le m\le1$,

$$
P(L>l,M\le m)=(m-l)^3,
$$

so $F_{L,M}(l,m)=m^3-(m-l)^3$. Differentiating first in $l$, then in $m$, gives $f_{L,M}(l,m)=6(m-l)$ on $0<l<m<1$. Outside this region the joint density is zero.

The full CDF is zero for $l<0$ or $m<0$; equals $m^3$ when $0\le m\le1$ and $l\ge m$; and equals $1-(1-l)^3$ for $0\le l\le1$ and $m\ge1$. It equals one when both arguments are at least one.

Since $P(L>l)=(1-l)^3$, $f_L(l)=3(1-l)^2$. Dividing gives

$$
f_{M\mid L}(m\mid l)=\frac{2(m-l)}{(1-l)^2},
\qquad l<m<1.
$$

Its integral over $m$ equals one. Once the minimum is fixed at $l$, the two remaining observations occupy $(l,1)$, and their maximum has this conditional density.

### 2.5 Additional solved problems

#### Problem PR-10: a triangular joint density

Constructed problem. Let $f(x,y)=c$ on $0<x<y<1$. Normalize it and compute $E[X\mid Y=y]$.

Solution. The triangle has area $1/2$, so $c=2$. Integrating $x$ from zero to $y$ gives $f_Y(y)=2y$. Consequently $f_{X\mid Y}(x\mid y)=1/y$ on $(0,y)$, a uniform distribution. Its mean is $y/2$. Integrating the other coordinate gives $f_X(x)=2(1-x)$ on $(0,1)$, so $EX=\int_0^12x(1-x)dx=1/3$. The tower property later gives the same answer from $E(Y/2)=(1/2)(2/3)$.

#### Problem PR-11: a mixed distribution

Constructed problem. Let $U$ be uniform on $(-1,1)$ and set $Y=\max(U,0)$.

Solution. Half the original interval is mapped to zero, so $P(Y=0)=1/2$. For $0<y<1$, $F_Y(y)=P(U\le y)=(y+1)/2$, giving a continuous density $1/2$ on that interval. The CDF jumps from zero to $1/2$ at zero. A description using only a density would omit the atom. The mean is $\int_0^1y/2\,dy=1/4$, with the atom contributing zero.

## 3. Conditional expectation

### 3.1 Theory

For discrete $Y$, $E[X\mid Y=y]$ averages $X$ under its conditional distribution given that value. The result as a function of the observed $Y$ is the random variable $E[X\mid Y]$. Before observing $Y$ it is random; after observing $Y=y$ it becomes a number.

More generally, for a sigma-algebra $\mathcal G$, conditional expectation $m=E[X\mid\mathcal G]$ is $\mathcal G$-measurable and satisfies $E[m\mathbf1_A]=E[X\mathbf1_A]$ for every $A\in\mathcal G$. This definition expresses exactly the information that $\mathcal G$ permits us to distinguish. It determines $m$ up to events of probability zero.

The tower property says $E[E[X\mid\mathcal G]]=EX$. If $\mathcal H\subseteq\mathcal G$, then $E[E[X\mid\mathcal G]\mid\mathcal H]=E[X\mid\mathcal H]$. A known factor can be pulled out: $E[ZX\mid\mathcal G]=ZE[X\mid\mathcal G]$ for suitable $\mathcal G$-measurable $Z$ and integrable products.

Independence of $X$ from $\mathcal G$ gives $E[X\mid\mathcal G]=EX$. Zero covariance gives much less: it controls only one linear average of joint behavior.

### 3.2 Derivations and proofs

#### Conditional expectation minimizes mean squared error

Assume $X\in L^2$ and put $m=E[X\mid\mathcal G]$. For any square-integrable $\mathcal G$-measurable predictor $g$, write $X-g=(X-m)+(m-g)$. Expansion gives

$$
E[(X-g)^2]=E[(X-m)^2]+E[(m-g)^2]
+2E[(X-m)(m-g)].
$$

Condition the cross term on $\mathcal G$. The factor $m-g$ is known there, and $E[X-m\mid\mathcal G]=0$, so the cross term vanishes. Thus

$$
E[(X-g)^2]=E[(X-m)^2]+E[(m-g)^2].
$$

The second term is nonnegative and is zero only when $g=m$ almost surely. This proves both optimality and uniqueness in the probabilistic sense.

#### Total variance

Apply the same orthogonal decomposition with $g=EX$. The residual contribution is $E[\operatorname{Var}(X\mid\mathcal G)]$, while the predictor contribution is $\operatorname{Var}(E[X\mid\mathcal G])$. Therefore

$$
\operatorname{Var}(X)=E[\operatorname{Var}(X\mid\mathcal G)]
+\operatorname{Var}(E[X\mid\mathcal G]).
$$

The first term is uncertainty left after observing the information. The second is variation in the conditional mean across the information states.

#### Conditional expectation as an information-based variable

For integrable $Y$ and information sigma-algebra $\mathcal G$, a conditional expectation $m=E[Y\mid\mathcal G]$ is a $\mathcal G$-measurable integrable variable satisfying

$$
E[m1_A]=E[Y1_A]\qquad(A\in\mathcal G).
$$

The measurable variable summarizes what can be predicted using the available information. Its defining equalities identify it up to events of probability zero. For a finite partition, these equalities reduce to within-cell averages.

If $\mathcal H\subseteq\mathcal G$, the tower property follows from the definition: for every $A\in\mathcal H$, integrating $E[Y\mid\mathcal G]$ over $A$ gives the same value as integrating $Y$. Therefore

$$
E[E[Y\mid\mathcal G]\mid\mathcal H]=E[Y\mid\mathcal H].
$$

The nesting is necessary. Conditioning in an arbitrary order on two unrelated information sets does not generally return either of the original conditional expectations.

For a bounded $\mathcal G$-measurable variable $Z$, one can take it outside the conditional expectation: $E[ZY\mid\mathcal G]=ZE[Y\mid\mathcal G]$. The identity extends under suitable integrability conditions. It does not apply to an unobserved future factor merely because its symbol is being treated as fixed in a calculation.

#### More information and average prediction error

Let $Y\in L^2$ and $\mathcal H\subseteq\mathcal G$. The projection identity gives

$$
E[(Y-E[Y\mid\mathcal H])^2]
=E[(Y-E[Y\mid\mathcal G])^2]
+E[(E[Y\mid\mathcal G]-E[Y\mid\mathcal H])^2].
$$

To derive it, split the first prediction error into the latter two terms. Their cross product has expectation zero after conditioning on $\mathcal G$. Thus more information cannot increase the minimum achievable mean squared prediction error. This is an average statement and does not compare every individual conditional variance.

#### Example PR-E23: one observation explains nonlinear variation

Let $X$ be uniform on $[-1,1]$ and $Y=X^2$. Find the best constant predictor and the best predictor based on $X$.

Solution. The best constant is $E[Y]=1/3$. Its squared prediction error is $\operatorname{Var}(Y)=E[X^4]-(E[X^2])^2=1/5-1/9=4/45$. Knowing $X$ determines $Y$ exactly, so $E[Y\mid X]=X^2$ and its error is zero.

Nevertheless $\operatorname{Cov}(X,Y)=E[X^3]-E[X]E[X^2]=0$. A linear regression on $X$ with an intercept gives the constant $1/3$ and misses the nonlinear relation. Conditional expectation ranges over all square-integrable functions of the information, while linear projection restricts the class of predictors.

### 3.3 Worked examples

#### Example PR-E05: a two-group mixture

Let $P(G=1)=1/4$, $E[X\mid G=1]=4$, and $E[X\mid G=0]=0$, with conditional variance one in both groups.

Solution. The mean is $(1/4)4+(3/4)0=1$. The average conditional variance is one. The conditional mean takes values four and zero, so its variance is $(1/4)(4-1)^2+(3/4)(0-1)^2=3$. Hence $\operatorname{Var}(X)=4$. Ignoring group differences would miss three quarters of the total variance.

#### Example PR-E06: dependence with zero covariance

Let $X$ be uniform on $(-1,1)$ and $Y=X^2$.

Solution. Symmetry gives $EX=EX^3=0$, so $\operatorname{Cov}(X,Y)=EX^3-EX\,EX^2=0$. Yet $Y$ is completely determined by $X$. In particular $E[Y\mid X]=X^2$, whereas $EY=1/3$. Zero covariance has not made the conditional mean constant.

### 3.4 Selected exercise and additional problems

#### Problem PR-12: shared opponent strength

Source: [Harvard Strategic Practice 2, Section 3, Problem 5(a–c)](https://stat110.hsites.harvard.edu/sites/g/files/omnuum10111/files/stat110/files/strategic_practice_and_homework_2.pdf).

Problem. An opponent is equally likely to be one of three types with win probabilities $0.9,0.5,0.3$. Games are independent conditional on type. Find the first-game win probability, the probability of two wins, and the second-game win probability after a first win.

Solution. Opponent types have equal probabilities and conditional win probabilities $0.9,0.5,0.3$. The first-game win probability is their average, $17/30$. Conditional on type, the games are independent, so the probability of two wins is $(0.9^2+0.5^2+0.3^2)/3=23/60$. Consequently

$$
P(W_2\mid W_1)=\frac{23/60}{17/30}=\frac{23}{34}.
$$

This exceeds $17/30$. A first win shifts posterior weight toward easier opponents. Conditional independence allows this common hidden type to create unconditional dependence. Indeed $P(W_1\cap W_2)=23/60\ne(17/30)^2$.

#### Problem PR-13: conditional Jensen inequality

Constructed problem. Prove $f(E[X\mid\mathcal G])\le E[f(X)\mid\mathcal G]$ when $X$ has finite support and $f$ is convex.

Solution. Let the possible values be $x_1,\ldots,x_k$ and put $p_i=P(X=x_i\mid\mathcal G)$. Almost surely the $p_i$ are nonnegative and sum to one. Then $E[X\mid\mathcal G]=\sum_ip_ix_i$ and $E[f(X)\mid\mathcal G]=\sum_ip_if(x_i)$. Apply the finite Jensen inequality pointwise to those weights. In particular $E[X\mid\mathcal G]^2\le E[X^2\mid\mathcal G]$, showing directly that conditional expectation does not increase the $L^2$ norm after taking expectations.

#### Problem PR-14: prediction from a coarse observation

Constructed problem. Let $X$ be uniform on $(0,1)$ and observe only $B=\mathbf1_{\{X>1/2\}}$. Find the best squared-error predictor and its MSE.

Solution. Given $B=0$, $X$ is uniform on $(0,1/2)$ with mean $1/4$. Given $B=1$, it is uniform on $(1/2,1)$ with mean $3/4$. The predictor is $m=1/4+B/2$. Each conditional variance is $(1/2)^2/12=1/48$, so the MSE is $1/48$. The unconditional variance is $1/12$. Observing the half-interval removes $1/12-1/48=1/16$ of variance, exactly the variance of the conditional mean.

#### Problem PR-15: conditioning on more information

Constructed problem. For $\mathcal H\subseteq\mathcal G$ and $X\in L^2$, compare the two optimal prediction errors.

Solution. Every $\mathcal H$-measurable predictor is also $\mathcal G$-measurable. The minimum over the larger predictor class cannot be greater. More precisely, let $m_H=E[X\mid\mathcal H]$ and $m_G=E[X\mid\mathcal G]$. Orthogonality gives

$$
E[(X-m_H)^2]=E[(X-m_G)^2]+E[(m_G-m_H)^2].
$$

The reduction in prediction error is exactly the squared distance between the two conditional means. Extra information gives no improvement precisely when the conditional means already agree almost surely.

## 4. Multivariate probability

### 4.1 Theory

For a random vector $X\in\mathbb R^p$ with finite second moments, define $\mu=EX$ and

$$
\Sigma=E[(X-\mu)(X-\mu)^\top].
$$

Its $(i,j)$ entry is $\operatorname{Cov}(X_i,X_j)$. The diagonal entries are variances. For a fixed matrix $A$ and vector $b$, $E[AX+b]=A\mu+b$ and $\operatorname{Cov}(AX+b)=A\Sigma A^\top$.

A jointly Gaussian vector is one for which every scalar linear combination is normally distributed, permitting degenerate zero-variance normals. It is characterized by its mean and covariance. A collection of variables can each be marginally normal without being jointly normal. Zero covariance implies independence for jointly Gaussian coordinates; it does not do so for general variables.

If $\Sigma$ is PD, the joint density is

$$
f(x)=\frac{\exp[-\tfrac12(x-\mu)^\top\Sigma^{-1}(x-\mu)]}
{(2\pi)^{p/2}\det(\Sigma)^{1/2}}.
$$

A singular Gaussian vector still exists, but is supported on a lower-dimensional affine subspace and has no full-dimensional density of this form.

### 4.2 Derivations and proofs

#### Conditional Gaussian distribution

Partition $X=(X_1^\top,X_2^\top)^\top$ and its covariance into blocks $\Sigma_{11},\Sigma_{12},\Sigma_{21},\Sigma_{22}$, with $\Sigma_{22}$ invertible. Set $B=\Sigma_{12}\Sigma_{22}^{-1}$ and

$$
R=X_1-\mu_1-B(X_2-\mu_2).
$$

Then $\operatorname{Cov}(R,X_2)=\Sigma_{12}-B\Sigma_{22}=0$. The pair $(R,X_2)$ is jointly Gaussian because it is a linear transformation of $X$, so zero cross-covariance makes the two blocks independent. Its residual covariance is $\Sigma_{11}-\Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21}$. Therefore

$$
X_1\mid X_2=x_2\sim
N\left(\mu_1+B(x_2-\mu_2),
\Sigma_{11}-\Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21}\right).
$$

The covariance is PSD because it is the covariance of the explicitly constructed residual. The independence step is the place where Gaussianity is essential.

#### Expected quadratic forms

For symmetric $A$, write $X=\mu+Z$ with $EZ=0$. Expand $X^\top AX$. The cross term has expectation zero, and

$$
E[Z^\top AZ]=E[\operatorname{tr}(AZZ^\top)]
=\operatorname{tr}(A\Sigma).
$$

Hence $E[X^\top AX]=\mu^\top A\mu+\operatorname{tr}(A\Sigma)$. No Gaussian assumption is needed.

#### Singular Gaussian distributions

For a random vector $X$, a covariance matrix is always PSD because $a^\top\Sigma a=\operatorname{Var}(a^\top X)\ge0$. It is positive definite exactly when no nonzero linear combination of the centered coordinates is almost surely constant.

A Gaussian vector can have singular covariance. For example, $(Z,2Z)$ is jointly Gaussian but lives on one line in the plane. It has no two-dimensional Lebesgue density. The familiar normal density involving $\det\Sigma$ and $\Sigma^{-1}$ applies only to positive-definite covariance.

If $X\sim N(\mu,\Sigma)$ and $A$ is fixed, then $AX+b$ is Gaussian with mean $A\mu+b$ and covariance $A\Sigma A^\top$. This remains true in singular cases and is often safer than trying to manipulate a density that does not exist.

#### Quadratic forms in independent normals

Let $Z\sim N(0,I_d)$ and $P$ be an orthogonal projection of rank $r$. Its eigenvalues are $r$ ones and the remaining zeros. Orthogonal invariance of the standard normal and diagonalization give

$$
Z^\top PZ=\sum_{i=1}^r U_i^2\sim\chi_r^2.
$$

If $P,Q$ are orthogonal projections with $PQ=0$, their projected Gaussian vectors have zero cross covariance and are independent. Their squared lengths are therefore independent too. Gaussianity is the assumption that converts geometric orthogonality into statistical independence.

#### Example PR-E24: why the sample variance uses a chi-square law

For independent $X_i\sim N(\mu,\sigma^2)$, write $X=\mu\mathbf1+\sigma Z$. Let $P=\mathbf1\mathbf1^\top/n$ and $M=I-P$. The sample mean depends on $PZ$, while centered residuals equal $\sigma MZ$.

The ranks are one and $n-1$, and $PM=0$. Thus

$$
\frac{\sum_i(X_i-\bar X)^2}{\sigma^2}
=Z^\top MZ\sim\chi_{n-1}^2,
$$

independently of $\bar X$. Dividing by $n-1$ makes $s^2$ unbiased because a chi-square with that many degrees of freedom has that expectation.

The t pivot combines the independent normal mean and chi-square variance estimate. These exact laws arise from the Gaussian vector's orthogonal decomposition. With non-Gaussian data, residual orthogonality still holds algebraically, but it does not in general give independence or an exact chi-square distribution.

### 4.3 Worked examples

#### Example PR-E07: a Gaussian spread

Suppose $(X,Y)$ is jointly Gaussian with means $1,2$, variances $4,9$, and covariance $3$. Find the law of $2X-Y$.

Solution. The mean is $2(1)-2=0$. Its variance is $4(4)+9-4(3)=13$, so $2X-Y\sim N(0,13)$. The covariance term has coefficient twice the product of the two linear coefficients, namely $2(2)(-1)=-4$.

#### Example PR-E08: a conditional Gaussian calculation

Using the same variables, $E[X\mid Y=y]=1+(3/9)(y-2)$ and $\operatorname{Var}(X\mid Y)=4-9/9=3$. At $Y=5$, the conditional mean is two. The variance falls from four to three; it does not fall to zero because knowing $Y$ leaves an independent residual component in $X$.

### 4.4 Selected exercise and additional problems

#### Problem PR-16: covariance positivity and linear transformations

Source: [MIT 18.650, Problem Set 10, Problem 3, Questions 3–5](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/ff133adcf6955bcaefcec1ea49014008_MIT18_650F16_PSet10.pdf).

Problem. Prove that population and sample covariance matrices are PSD, derive covariance under a linear transformation, and explain the rank obstruction when the output dimension exceeds the input dimension.

Solution. For any fixed $a$, $a^\top\Sigma a=E[(a^\top(X-\mu))^2]=\operatorname{Var}(a^\top X)\ge0$, proving PSD. For a sample covariance, replace expectation by the average of centered outer products; the same sum-of-squares argument works with either positive denominator convention.

For $Y=AX$, center it as $Y-EY=A(X-\mu)$ and take outer products, yielding $\operatorname{Cov}(Y)=A\Sigma A^\top$. If $A$ is $q\times p$ with $q>p$, its covariance has rank at most $p<q$ and cannot be invertible. Finally $\operatorname{Var}(u^\top X)=u^\top\Sigma u$ by the same transformation formula with a one-row matrix.

#### Problem PR-17: normal marginals without a Gaussian pair

Constructed problem. Let $Z\sim N(0,1)$ and let an independent sign $S$ equal $\pm1$ with equal probabilities. Put $X=Z$, $Y=SZ$.

Solution. Both marginals are standard normal because changing the sign of a centered normal preserves its law. Also $E[XY]=E[S]E[Z^2]=0$. Nevertheless $Y^2=X^2$ with probability one, so the pair is dependent. Its support lies on two crossing lines, whereas a jointly Gaussian pair with identity covariance has a two-dimensional density. Thus the zero-covariance independence rule cannot be applied merely from normal marginal distributions.

#### Problem PR-18: correlated normals from independent ones

Constructed problem. For independent standard normals $Z_1,Z_2$, define $X=Z_1$ and $Y=\rho Z_1+\sqrt{1-\rho^2}Z_2$ for $|\rho|<1$.

Solution. The pair is jointly Gaussian as a linear transformation of independent normals. Both variances equal one and the covariance is $\rho$. Given $X=x$, the remaining $Z_2$ is independent and unchanged, so $Y\mid X=x\sim N(\rho x,1-\rho^2)$. This construction directly explains both the conditional mean and variance. At $|\rho|=1$ the residual disappears and the covariance matrix becomes singular.

#### Problem PR-19: a chi-square quadratic

Constructed problem. If $X\sim N_p(\mu,\Sigma)$ with $\Sigma$ PD, find the law of $(X-\mu)^\top\Sigma^{-1}(X-\mu)$.

Solution. Let $Z=\Sigma^{-1/2}(X-\mu)$. This is jointly Gaussian with mean zero and covariance identity, so its coordinates are independent standard normals. The quadratic becomes $Z^\top Z=\sum_{j=1}^pZ_j^2$, which by definition has the $\chi_p^2$ distribution. Its mean is $p$ and variance $2p$. The invertibility assumption is necessary for this exact expression; a rank-$r$ singular covariance with its pseudoinverse yields $\chi_r^2$ on the supported subspace.

## 5. Convergence and asymptotic theory

### 5.1 Theory

There are several ways for random variables to approach a limit. Convergence in probability, $X_n\to_p X$, means $P(|X_n-X|>\varepsilon)\to0$ for each fixed $\varepsilon>0$. Convergence in $L^2$ means $E[(X_n-X)^2]\to0$. Almost-sure convergence means that the numerical sequence $X_n(\omega)$ converges for all outcomes outside a probability-zero set.

Convergence in distribution, $X_n\Rightarrow X$, means $F_{X_n}(x)\to F_X(x)$ at every continuity point of $F_X$. It compares laws rather than requiring the variables to become close on a common probability space. Convergence in probability implies convergence in distribution; the converse generally fails. Convergence in $L^2$ implies convergence in probability by Markov's inequality applied to squared error.

The weak law of large numbers says an iid sample mean approaches its expectation under suitable integrability assumptions. With finite variance the proof is short. The central limit theorem goes further: for iid observations with mean $\mu$ and finite positive variance $\sigma^2$,

$$
\sqrt n(\bar X_n-\mu)\Rightarrow N(0,\sigma^2).
$$

The CLT describes the size and shape of the remaining error, at scale $n^{-1/2}$. The limiting normal approximation is to the standardized mean, not necessarily to an individual observation.

The continuous mapping theorem preserves convergence under functions continuous at the relevant limits. Slutsky's theorem allows estimated quantities converging in probability to constants to be inserted into a distributional limit. It does not say that two separately converging marginal distributions may be combined without joint information.

### 5.2 Derivations and proofs

#### Chebyshev's inequality and the weak law

For a finite-variance variable, apply Markov's inequality to $(X-EX)^2$:

$$
P(|X-EX|\ge a)\le\frac{\operatorname{Var}(X)}{a^2}.
$$

For iid $X_i$ with variance $\sigma^2$, independence removes covariance terms and $\operatorname{Var}(\bar X_n)=\sigma^2/n$. Therefore

$$
P(|\bar X_n-\mu|\ge\varepsilon)\le\frac{\sigma^2}{n\varepsilon^2}\to0.
$$

The same calculation gives $E[(\bar X_n-\mu)^2]=\sigma^2/n$, proving $L^2$ convergence. With serial dependence, the variance contains additional covariance terms and the iid proof must be changed.

#### Delta method

Suppose $\sqrt n(\widehat\theta-\theta)\Rightarrow N(0,V)$ and $g$ is differentiable at $\theta$. Taylor expansion gives

$$
g(\widehat\theta)-g(\theta)
=g'(\theta)(\widehat\theta-\theta)+o_p(n^{-1/2}).
$$

Hence $\sqrt n[g(\widehat\theta)-g(\theta)]\Rightarrow N(0,[g'(\theta)]^2V)$. In several dimensions, replace the derivative by the row gradient: the variance is $\nabla g(\theta)^\top V\nabla g(\theta)$ when gradients are columns. The remainder is negligible because the estimator is consistent and its scaled error stays bounded in probability.

For a ratio $g(a,b)=a/b$, with $b\ne0$, the gradient is $(1/b,-a/b^2)^\top$. If the joint limit covariance of $(\widehat a,\widehat b)$ is $V$, both diagonal variance terms and their cross-covariance contribute to the ratio's variance.

#### A short route to the iid CLT

A full CLT proof is beyond the elementary inequalities above, but its finite-variance mechanism can be stated precisely. For centered variance-one $Z$, its characteristic function satisfies $E[e^{itZ}]=1-t^2/2+o(t^2)$ near zero. Independence makes the characteristic function of $n^{-1/2}\sum_iZ_i$ equal to $[1-t^2/(2n)+o(1/n)]^n\to e^{-t^2/2}$. The characteristic-function continuity theorem identifies the limit as standard normal. The expansion and that continuity theorem are the non-elementary ingredients invoked here.

#### Comparing modes of convergence

Almost-sure convergence controls paths except on a null event. Convergence in probability controls $P(|X_n-X|>\epsilon)$ for every fixed positive tolerance. Convergence in $L^2$ controls $E|X_n-X|^2$. Convergence in distribution controls CDF limits at continuity points and may refer only to laws rather than variables on one common space.

The implication $L^2\Rightarrow$ probability follows from Markov's inequality applied to the squared error. Almost-sure convergence implies convergence in probability by dominated convergence applied to the error-event indicators. Convergence in probability implies convergence in distribution, but the converse generally fails.

For a rare-spike example, let $X_n=\sqrt n$ with probability $1/n$ and zero otherwise. Then $P(|X_n|>\epsilon)\to0$ for fixed $\epsilon>0$ and $E|X_n|=1/\sqrt n\to0$. Yet $E[X_n^2]=1$ for every $n$, so there is no $L^2$ convergence to zero. Rare large errors can disappear in probability while retaining a nonvanishing second moment.

For another example, let $Z\sim N(0,1)$ and $X_n=(-1)^nZ$. Every $X_n$ has the standard normal law, so the distributions converge to that law. The sequence does not converge in probability to $Z$: on odd indices, its difference from $Z$ is $-2Z$. A limiting marginal law does not determine a pathwise coupling.

#### The multivariable delta method for a ratio

Suppose $\sqrt n((\hat a,\hat b)-(a,b))\Rightarrow N(0,\Sigma)$ with $b\ne0$. For $g(a,b)=a/b$, the gradient is $(1/b,-a/b^2)^\top$. A first-order Taylor expansion and consistency give

$$
\sqrt n\left(\frac{\hat a}{\hat b}-\frac ab\right)
\Rightarrow N(0,\nabla g^\top\Sigma\nabla g).
$$

Writing $R=a/b$, the asymptotic variance is

$$
\frac{\Sigma_{11}-2R\Sigma_{12}+R^2\Sigma_{22}}{b^2}.
$$

The covariance term is essential when the numerator and denominator use the same observations. Independence cannot be inferred from two separate marginal CLTs. If the population denominator is zero, this differentiable-map argument does not apply.

#### Example PR-E25: a ratio with correlated estimation error

Take $a=2$, $b=1$, and asymptotic covariance $\Sigma=\begin{pmatrix}4&1\\1&1\end{pmatrix}$. The ratio is two and the delta variance is $4-2(2)(1)+4(1)=4$. Thus its large-sample standard error is approximately $2/\sqrt n$.

Ignoring the covariance would give variance eight and a different interval. A denominator near zero would make the variance large and a symmetric normal interval unreliable. Consistency of the denominator away from zero is part of the approximation, rather than an arithmetic detail.

### 5.3 Worked examples

#### Example PR-E09: an average of exponentials

For iid $\operatorname{Exp}(2)$ observations, estimate $P(\bar X_{100}>0.6)$.

Solution. The mean is $1/2$ and variance $1/4$, so the sample-mean standard deviation is $0.05$. The standardized threshold is $(0.6-0.5)/0.05=2$. The CLT approximation is $1-\Phi(2)\approx0.0228$. This is an approximation, not an exact identity; an exact answer could instead use the gamma distribution of the sum.

#### Example PR-E10: estimating a logarithm

If $X_i$ are iid Bernoulli with $0<p<1$, then $\sqrt n(\bar X-p)\Rightarrow N(0,p(1-p))$. For $g(p)=\log p$, the derivative is $1/p$, so the asymptotic variance of $\sqrt n(\log\bar X-\log p)$ is $(1-p)/p$. The event $\bar X=0$ has probability $(1-p)^n\to0$; define the estimator on that event separately before making a finite-sample statement.

#### Example PR-E11: a zero derivative changes the rate

Suppose $\sqrt n\widehat\theta\Rightarrow N(0,\sigma^2)$ at true $\theta=0$, and estimate $\theta^2$ by $\widehat\theta^2$.

Solution. The ordinary first-order delta variance is zero because $g'(0)=0$. At scale $n$, however, $n\widehat\theta^2=(\sqrt n\widehat\theta)^2\Rightarrow\sigma^2\chi_1^2$. The faster rate and nonnormal limit require a different confidence-interval calculation.

### 5.4 Selected exercises

#### Problem PR-20: rare large values

Source: [MIT 18.650, Problem Set 1, Problem 1, Question 1](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/1aea12bd54eb19c821501f730cc247d9_MIT18_650F16_PSet1.pdf).

Problem. Let $X_n=1/n$ with probability $1-1/n^2$ and $X_n=n$ with probability $1/n^2$. Compare convergence in probability and in $L^2$.

Solution. Let $P(X_n=1/n)=1-1/n^2$ and $P(X_n=n)=1/n^2$. For any fixed $\varepsilon>0$, once $1/n<\varepsilon$, the probability of exceeding $\varepsilon$ is at most $1/n^2$, so $X_n\to_p0$. But

$$
E[X_n^2]=(1-1/n^2)/n^2+n^2/n^2
=1+1/n^2-1/n^4\to1.
$$

Thus it does not converge to zero in $L^2$. It cannot have any different $L^2$ limit, because an $L^2$ limit would also be a probability limit, and probability limits are unique almost surely. Rare large values vanish in probability while retaining a nonvanishing second-moment contribution.

#### Problem PR-21: Bernoulli sample means

Source: [MIT 18.650, Problem Set 1, Problem 1, Question 2(a,b)](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/1aea12bd54eb19c821501f730cc247d9_MIT18_650F16_PSet1.pdf).

Problem. For independent Bernoulli observations with success probability $p$, identify the law of the sample sum and prove $L^2$ convergence of the sample mean to $p$.

Solution. The sum $n\bar X_n$ counts successes in $n$ independent Bernoulli trials, so it has the $\operatorname{Bin}(n,p)$ distribution. Its mean is $np$ and variance $np(1-p)$. Dividing by $n$ gives $E\bar X_n=p$ and $\operatorname{Var}(\bar X_n)=p(1-p)/n$. Therefore $E[(\bar X_n-p)^2]\to0$, proving $L^2$ convergence, including the degenerate endpoints $p=0,1$.

#### Problem PR-22: marginal distributional limits do not determine a sum

Source: [MIT 18.650, Problem Set 1, Problem 2, Question 3](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/1aea12bd54eb19c821501f730cc247d9_MIT18_650F16_PSet1.pdf).

Problem. Determine whether separate convergence in distribution of $X_n$ and $Y_n$ implies convergence of their sum to the sum of specified limiting variables. Give a counterexample if not.

Solution. Let $Z\sim N(0,1)$, $X_n=Z$, and $Y_n=-Z$ for every $n$. Choose limit variables $X=Z$ and $Y=Z$. Both marginal convergences $X_n\Rightarrow X$ and $Y_n\Rightarrow Y$ hold, since $Z$ and $-Z$ have identical laws. Yet $X_n+Y_n=0$ always, while $X+Y=2Z$ is nondegenerate. Thus the asserted sum convergence fails. Joint convergence to the specified pair would make continuous mapping applicable, but separate marginal convergence does not supply it.

### 5.5 Additional solved problems

#### Problem PR-23: a ratio of sample averages

Constructed problem. Suppose iid pairs $(X_i,Y_i)$ have mean $(a,b)$, $b\ne0$, variances $\sigma_X^2,\sigma_Y^2$, and covariance $\sigma_{XY}$. Find the asymptotic variance of $\bar X/\bar Y$.

Solution. The joint CLT and the ratio gradient give

$$
\sqrt n\left(\frac{\bar X}{\bar Y}-\frac ab\right)
\Rightarrow N\left(0,\frac{\sigma_X^2}{b^2}
+\frac{a^2\sigma_Y^2}{b^4}
-\frac{2a\sigma_{XY}}{b^3}\right).
$$

The denominator stays away from zero with probability tending to one because $\bar Y\to_p b\ne0$. The covariance term cannot be discarded unless its coefficient vanishes or the pairs have zero cross-covariance. Separate independence across observations does not imply independence of the two components within each pair.

#### Problem PR-24: a sample-size guarantee

Constructed problem. With known bound $\operatorname{Var}(X)\le4$, find a sufficient iid sample size for $P(|\bar X-\mu|\ge0.1)\le0.05$.

Solution. Chebyshev gives a bound $4/(0.01n)=400/n$. Requiring $400/n\le0.05$ gives $n\ge8000$. This is a finite-sample guarantee under only the variance bound and independence. A normal approximation can suggest a smaller size, but that approximation supplies different information and does not follow from Chebyshev alone.

#### Problem PR-25: why finite variance matters for the mean CLT

Constructed problem. For iid standard Cauchy variables, explain why the usual sample-mean LLN and finite-variance CLT do not apply.

Solution. A standard Cauchy has neither an integrable absolute first moment nor a finite variance. Its characteristic function is $e^{-|t|}$. The sample mean has characteristic function $\prod_{i=1}^ne^{-|t|/n}=e^{-|t|}$, so its distribution is the same Cauchy law for every $n$. It does not concentrate around a constant. Increasing the sample size does not repair a theorem whose moment assumptions fail.

## 6. Statistical estimation

### 6.1 Theory

A statistical model is a family $\{P_\theta:\theta\in\Theta\}$ of possible sampling distributions. Identification means different parameter values produce different observable distributions. An estimator $\widehat\theta$ is a function of the observed sample. Its sampling distribution describes how it varies under repeated samples with the same true parameter.

Bias is $E_\theta\widehat\theta-\theta$. Mean squared error decomposes as

$$
E_\theta[(\widehat\theta-\theta)^2]
=\operatorname{Var}_\theta(\widehat\theta)
+[\operatorname{Bias}_\theta(\widehat\theta)]^2.
$$

Unbiasedness is a finite-sample mean statement. Consistency is convergence to the true parameter as sample size increases. Neither property alone establishes small error at a particular sample size.

Method of moments equates sample moments with their model expectations. Maximum likelihood chooses a parameter making the observed sample density or mass largest. For iid data, the likelihood is $\prod_i f_\theta(x_i)$ and its log is $\ell_n(\theta)=\sum_i\log f_\theta(x_i)$. The observations are fixed when maximizing this expression; the parameter is varied.

The score is $s_\theta(X)=\partial_\theta\log f_\theta(X)$. In a regular model, Fisher information per observation is $I(\theta)=E_\theta[s_\theta(X)^2]=-E_\theta[\partial_\theta^2\log f_\theta(X)]$. Interchanging derivatives and integrals, a common support, and the required integrability justify these equalities. Parameter-dependent support can invalidate them.

### 6.2 Derivations and proofs

#### The information lower bound

For a regular scalar model and an unbiased estimator $T$ of $\theta$, differentiating $E_\theta T=\theta$ gives $E_\theta[T S_n]=1$, where $S_n=\sum_i s_\theta(X_i)$. Since $ES_n=0$, this is $\operatorname{Cov}(T,S_n)=1$. Cauchy-Schwarz gives

$$
1\le\operatorname{Var}(T)\operatorname{Var}(S_n)
=\operatorname{Var}(T)\,nI(\theta).
$$

Thus $\operatorname{Var}(T)\ge1/[nI(\theta)]$. The differentiation step is an assumption-dependent part of the argument. This bound is for unbiased estimators in the stated regular setting, not for every possible biased procedure.

#### A regular MLE's asymptotic variance

If an interior consistent MLE solves $\ell_n'(\widehat\theta)=0$, expand around the truth:

$$
0=\ell_n'(\theta)+\ell_n''(\widetilde\theta)(\widehat\theta-\theta).
$$

A score CLT gives $n^{-1/2}\ell_n'(\theta)\Rightarrow N(0,I(\theta))$, and a uniform law of large numbers under regularity gives $-n^{-1}\ell_n''(\widetilde\theta)\to_p I(\theta)$. Rearrangement and Slutsky yield asymptotic variance $I(\theta)^{-1}$. Consistency, interiority, differentiability, and the laws of large numbers are requirements, not consequences of merely writing the likelihood.

#### Information bounds require regularity

For a regular scalar model with score $S_\theta=\partial_\theta\log f_\theta(X)$, differentiation under the integral gives $E_\theta[S_\theta]=0$. If $T$ is unbiased for $\theta$ and the same interchange is justified, differentiating $E_\theta[T]=\theta$ gives

$$
E_\theta[(T-\theta)S_\theta]=1.
$$

Cauchy–Schwarz then yields $\operatorname{Var}(T)I(\theta)\ge1$. For independent observations, information adds, giving the Cramér–Rao lower bound $1/(nI_1(\theta))$ for regular unbiased estimators.

The proof identifies the assumptions behind the bound. Parameter-dependent support can contribute boundary terms and invalidate the score identities. A smaller variance in such a model is not a contradiction.

#### Example PR-E26: a nonregular uniform endpoint

Let $X_i$ be independent uniform variables on $(0,\theta)$ with $\theta>0$. Find the maximum-likelihood estimate and an unbiased correction.

Solution. For positive observations, the likelihood is $\theta^{-n}1_{\{\theta\ge M\}}$ with $M=\max_iX_i$. It decreases for $\theta\ge M$, so the MLE is $M$. The support depends on the parameter, and the maximum occurs at its boundary.

For $0<m<\theta$, $P(M\le m)=(m/\theta)^n$. Differentiation and integration give

$$
E[M]=\frac n{n+1}\theta,\qquad
\operatorname{Var}(M)=
\frac{n\theta^2}{(n+1)^2(n+2)}.
$$

Thus $T=(n+1)M/n$ is unbiased, with variance $\theta^2/[n(n+2)]$. Its variance shrinks at order $n^{-2}$. The regular score calculation would wrongly omit the moving support boundary. It cannot be used to claim a usual $n^{-1}$ lower bound here.

#### Bias, consistency, and efficiency

Bias is a finite-sample expectation difference. Consistency is a limiting probability statement. An estimator can be biased for every finite sample and still be consistent, as the normal variance MLE illustrates. It can also be unbiased but inconsistent when its variance fails to shrink.

Mean squared error decomposes as $E[(T-\theta)^2]=\operatorname{Var}(T)+\operatorname{Bias}(T)^2$ whenever the second moment exists. Comparing only unbiasedness can therefore select a less accurate estimator under squared error. “Efficient” must specify the comparison class and whether it concerns finite samples or an asymptotic covariance.

A likelihood is a function of the parameter with the observations fixed. Its numerical value is not the probability that a continuous parameter equals a candidate. Maximum likelihood chooses the parameter that best supports the observed data under the specified model; it does not verify the model itself.

### 6.3 Worked examples

#### Example PR-E12: method of moments for an exponential rate

For iid exponentials, $EX=1/\lambda$. Equating $\bar X$ to this mean gives $\widehat\lambda=1/\bar X$. The LLN and continuity of reciprocal at $1/\lambda>0$ prove consistency. Applying the delta method with derivative $-\lambda^2$ gives asymptotic variance $\lambda^2$ for the $\sqrt n$ error.

#### Example PR-E13: an unbiased estimator can be noisy

For iid observations with mean $\mu$, both $X_1$ and $\bar X_n$ are unbiased estimators of $\mu$. If the variance is $\sigma^2>0$, their variances are $\sigma^2$ and $\sigma^2/n$. The first ignores all new observations and is generally not consistent; the second uses the expanding sample. Unbiasedness alone does not choose between them.

#### Example PR-E14: information lost by thresholding

If $Y\sim\operatorname{Exp}(\lambda)$ but only $X=\mathbf1_{\{Y>z\}}$ is observed, then $X$ is Bernoulli with $p=e^{-\lambda z}$. Its information for $\lambda$ is

$$
I_X(\lambda)=\frac{(dp/d\lambda)^2}{p(1-p)}
=\frac{z^2e^{-\lambda z}}{1-e^{-\lambda z}}.
$$

The full observation has information $1/\lambda^2$. Thresholding keeps only whether the waiting time crossed one cutoff and discards its exact length.

### 6.4 Selected exercises

#### Problem PR-26: correcting a variance estimator's bias

Source: [MIT 18.650, Problem Set 2, Problem 1, Questions 1–4](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/680f2b6e1bd22ea4ce020accf752cbe9_MIT18_650F16_PSet2.pdf).

Problem. For a Bernoulli sample, examine $\bar X(1-\bar X)$ as an estimator of the population variance: establish consistency, compute bias, and construct an unbiased correction.

Solution. Bernoulli $X$ satisfies $X^2=X$, so $\operatorname{Var}(X)=p-p^2$. Since $\bar X\to_p p$, the continuous mapping theorem makes $\bar X(1-\bar X)$ consistent for $p(1-p)$. Its expectation is

$$
E[\bar X-\bar X^2]
=p-\left(p^2+\frac{p(1-p)}n\right)
=\left(1-\frac1n\right)p(1-p).
$$

The bias is $-p(1-p)/n$. For $n>1$, multiplying by $n/(n-1)$ gives an unbiased estimator. For Bernoulli observations this equals the ordinary sample variance with denominator $n-1$.

#### Problem PR-27: Bernoulli likelihood and information

Source: [MIT 18.650, Problem Set 4, Problem 1, Question 1](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/b6bf22e01213547707d09e588819e94f_MIT18_650F16_PSet4.pdf).

Problem. Derive the Bernoulli MLE, score, and Fisher information. Treat all-success and all-failure samples and distinguish open from closed parameter spaces.

Solution. With $k=\sum_i x_i$, $\ell(p)=k\log p+(n-k)\log(1-p)$. For $0<k<n$, stationarity gives $k/p=(n-k)/(1-p)$, hence $\widehat p=k/n$. The second derivative is negative, proving the unique maximum.

If $k=0$ or $k=n$ and the parameter space is strictly $(0,1)$, the likelihood has a supremum at the missing endpoint but no MLE inside the parameter space. On the closed model $[0,1]$, the endpoint estimate is allowed. This finite-sample boundary event must be distinguished from interior asymptotic regularity.

For one observation, the score is $(X-p)/[p(1-p)]$. Its mean is zero and its second moment is $I(p)=1/[p(1-p)]$. The sample mean's variance $p(1-p)/n$ attains the unbiased information bound.

#### Problem PR-28: exponential likelihood and information

Source: [MIT 18.650, Problem Set 4, Problem 1, Question 3](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/b6bf22e01213547707d09e588819e94f_MIT18_650F16_PSet4.pdf).

Problem. For an independent exponential sample parameterized by rate $\lambda>0$, derive the MLE, score, Fisher information, and finite-sample expectation of the MLE.

Solution. The log likelihood is $n\log\lambda-\lambda\sum_i x_i$, with derivative $n/\lambda-\sum_i x_i$. For a positive observed sum, the maximum is $\widehat\lambda=n/\sum_i x_i=1/\bar X$, and the negative second derivative proves uniqueness. The per-observation score is $1/\lambda-X$, with variance $1/\lambda^2$, so $I(\lambda)=1/\lambda^2$. The rate estimator is consistent but not exactly unbiased: for $n>1$, the gamma law of the sum gives $E\widehat\lambda=n\lambda/(n-1)$.

#### Problem PR-29: the normal mean and variance

Source: [MIT 18.650, Problem Set 4, Problem 1, Question 4](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/b6bf22e01213547707d09e588819e94f_MIT18_650F16_PSet4.pdf).

Problem. For a normal sample with both mean and strictly positive variance unknown, derive the MLEs and information matrix, and examine the all-equal observation case.

Solution. Write $v=\sigma^2>0$. The log likelihood is $-\tfrac n2\log(2\pi v)-\sum_i(x_i-\mu)^2/(2v)$. At fixed $v$, completing the square in $\mu$ gives $\widehat\mu=\bar X$. Then differentiating in $v$ gives $\widehat v=n^{-1}\sum_i(X_i-\bar X)^2$, assuming a nonzero residual sum.

The score components for one observation are $(X-\mu)/v$ and $-1/(2v)+(X-\mu)^2/(2v^2)$. Centered Gaussian odd moments vanish, so their cross-information is zero. Using $E(X-\mu)^4=3v^2$ gives

$$
I(\mu,v)=\begin{pmatrix}1/v&0\\0&1/(2v^2)\end{pmatrix}.
$$

For an all-equal sample, the likelihood increases without bound as $v\downarrow0$ and no MLE exists with strictly positive variance. With $n\ge2$ under a genuine continuous normal model, this event has probability zero.

### 6.5 Additional solved problems

#### Problem PR-30: identifying parameters after keeping only signs

Constructed problem. If $X\sim N(\mu,\sigma^2)$ and only $\mathbf1_{\{X>0\}}$ is observed, can both parameters be identified?

Solution. The observable success probability is $\Phi(\mu/\sigma)$. Thus the sign distribution identifies the ratio $\mu/\sigma$ but not the pair. For any $c>0$, $(c\mu,c\sigma)$ produces the same observations. No amount of sign-only sampling distinguishes these pairs. The issue is identification, so improving an estimation algorithm cannot resolve it.

#### Problem PR-31: a uniform maximum

Constructed problem. For iid $U(0,\theta)$ observations, find the mean, bias, and a correction for $M_n=\max_iX_i$.

Solution. Its CDF is $(m/\theta)^n$ for $0\le m\le\theta$, so its density is $nm^{n-1}/\theta^n$. Integration gives $EM_n=n\theta/(n+1)$. Its bias is $-\theta/(n+1)$ and $(n+1)M_n/n$ is unbiased. Also $P(M_n\le\theta-\varepsilon)=(1-\varepsilon/\theta)^n\to0$ for $0<\varepsilon<\theta$, proving consistency. The parameter-dependent support places this model outside the usual regular score argument.

#### Problem PR-32: bias-variance tradeoff for a normal mean

Constructed problem. Let $X_i\sim N(\mu,\sigma^2)$ iid and estimate $\mu$ by $a\bar X$ for fixed $a$.

Solution. The bias is $(a-1)\mu$ and variance is $a^2\sigma^2/n$. Thus $\operatorname{MSE}=a^2\sigma^2/n+(a-1)^2\mu^2$. Differentiating in $a$ gives the oracle choice $a=\mu^2/(\mu^2+\sigma^2/n)$. It depends on the unknown $\mu$, so it is not directly an available rule. The calculation nevertheless explains why some shrinkage can reduce MSE near zero while adding bias.

#### Problem PR-33: information under reparameterization

Constructed problem. For an exponential distribution, replace rate $\lambda$ by mean $m=1/\lambda$. Find information for $m$.

Solution. The chain rule gives $s_m=s_\lambda(d\lambda/dm)$, so $I_m(m)=I_\lambda(1/m)(-1/m^2)^2$. Since $I_\lambda(\lambda)=1/\lambda^2$, this becomes $m^2/m^4=1/m^2$. The MLE for $m$ is $\bar X$, with variance $m^2/n$. Information has units tied to the chosen parameter; its numerical value changes when the parameter is rescaled.

## 7. Statistical inference

### 7.1 Theory

A confidence interval is a random set $C(X_1,\ldots,X_n)$ whose repeated-sampling coverage is $P_\theta(\theta\in C)\ge1-\alpha$ or approaches that level asymptotically. In a frequentist model, $\theta$ is fixed and the interval is random before sampling. A realized interval either contains it or does not; a posterior probability statement requires a different model including a prior.

A test rejects a null hypothesis on a rejection region. Its Type I error probability is rejection under the null. Its Type II error is failure to reject at an alternative. Power is one minus Type II error. For a composite null, level $\alpha$ means the rejection probability is at most $\alpha$ at every null parameter.

A p-value is a null-calibrated measure of how extreme the observed statistic is in the specified direction. It is not the probability the null is true. Failing to reject does not establish equality or practical equivalence.

A Wald test standardizes an estimator's distance from the null. A likelihood-ratio test compares the best fit under the null with the unrestricted best fit. Under regular interior null hypotheses, the log-likelihood ratio multiplied by two has an asymptotic chi-square law, with degrees of freedom equal to the number of restrictions. Boundary and unidentified-nuisance cases require separate limits.

### 7.2 Derivations and proofs

#### From an asymptotic pivot to an interval

If $\sqrt n(\widehat\theta-\theta)/\widehat s\Rightarrow N(0,1)$ and $\widehat s\to_p s>0$, then

$$
P_\theta\left(
\widehat\theta-z_{1-\alpha/2}\widehat s/\sqrt n
\le\theta\le
\widehat\theta+z_{1-\alpha/2}\widehat s/\sqrt n
\right)\to1-\alpha.
$$

The event is obtained by solving $-z\le\sqrt n(\widehat\theta-\theta)/\widehat s\le z$ for $\theta$. The consistency of $\widehat s$ and positivity of $s$ justify studentization.

For iid normal data with unknown variance, an exact pivot is $\sqrt n(\bar X-\mu)/S\sim t_{n-1}$, where $S^2=\sum_i(X_i-\bar X)^2/(n-1)$. Orthogonal Gaussian projections make the sample mean independent of the residual sum of squares, and the latter divided by $\sigma^2$ has $\chi_{n-1}^2$ law. The ratio defines the Student distribution.

#### Error rates, power, and a rejection rule

For a test function $\varphi$ equal to one on rejection, a level-$\alpha$ test satisfies $\sup_{\theta\in H_0}E_\theta[\varphi]\le\alpha$. Power at an alternative is $E_\theta[\varphi]$. Type II error at that alternative is one minus power. An observed p-value is not the probability that the null hypothesis is true.

For independent Gaussian observations with known $\sigma$, test $H_0:\mu\le\mu_0$ against $\mu>\mu_0$ by rejecting when

$$
\bar X>\mu_0+z_{1-\alpha}\frac{\sigma}{\sqrt n}.
$$

The largest rejection probability within the null occurs at the boundary $\mu_0$. At alternative $\mu_1>\mu_0$, power is

$$
1-\Phi\left(z_{1-\alpha}
-\frac{\sqrt n(\mu_1-\mu_0)}{\sigma}\right).
$$

Solving for target power $1-\beta$ gives the sufficient sample size $n\ge[(z_{1-\alpha}+z_{1-\beta})\sigma/(\mu_1-\mu_0)]^2$.

#### Example PR-E27: planning a one-sided mean test

Take $\alpha=0.05$, target power $0.8$, known $\sigma=2$, and an effect $\mu_1-\mu_0=0.5$. The normal quantiles are approximately $1.644854$ and $0.841621$. The required continuous sample size is about $98.92$, so use at least ninety-nine observations.

This calculation assumes the supplied variance, effect, direction, and Gaussian sampling model. A two-sided test uses a different critical value and generally needs more observations for the same specified alternative. Estimating variance or accounting for dependence changes the design calculation.

#### Three large-sample tests for a Bernoulli probability

For an interior null $p_0$ and sample proportion $\hat p$, the likelihood-ratio statistic is

$$
2n\left[
\hat p\log\frac{\hat p}{p_0}
+(1-\hat p)\log\frac{1-\hat p}{1-p_0}
\right].
$$

Use the continuous convention $0\log0=0$ at an endpoint estimate. A Wald statistic uses $n(\hat p-p_0)^2/[\hat p(1-\hat p)]$, while the score statistic uses the null variance $n(\hat p-p_0)^2/[p_0(1-p_0)]$. Under fixed interior $p_0$ and large samples, all have a chi-square limit with one degree of freedom.

At a boundary estimate, the Wald denominator can vanish even though the exact binomial problem is well defined. Exact binomial tail probabilities avoid that asymptotic denominator. Discrete tests often have rejection probability below the nominal level because no integer cutoff gives exactly the desired tail mass.

### 7.3 Worked examples

#### Example PR-E15: a known-variance normal interval

For $n=100$, observed mean $3$, and known $\sigma=2$, the exact normal 95% interval is $3\pm1.96(2/10)$, or $[2.608,3.392]$. Exactness uses normal sampling and known variance. With a general finite-variance population, the same form is asymptotic; with estimated variance, an additional consistency argument is needed.

#### Example PR-E16: a Bernoulli Wald interval

For $n=100$ and sixty successes, $\widehat p=0.6$ and estimated standard error is $\sqrt{0.6(0.4)/100}\approx0.0490$. The Wald interval is approximately $[0.504,0.696]$. Near zero or one, this approximation can have poor coverage and can extend outside $[0,1]$. Wilson or exact binomial inversion uses a different construction.

#### Example PR-E17: power of a one-sided normal test

With known $\sigma$, reject $H_0:\mu\le\mu_0$ if $\sqrt n(\bar X-\mu_0)/\sigma>z_{1-\alpha}$. At true $\mu$, power is

$$
1-\Phi\left(z_{1-\alpha}-\frac{\sqrt n(\mu-\mu_0)}{\sigma}\right).
$$

At the boundary it equals $\alpha$; farther inside the null it is smaller. At a fixed alternative above the boundary it approaches one as $n$ grows.

### 7.4 Selected exercises

#### Problem PR-34: tests for an exponential rate

Source: [MIT 18.650, Problem Set 6, Problem 1, Questions 1–3](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/99e0952a99c4d37c82f15f3728fdd7a8_MIT18_650F16_PSet6.pdf).

Problem. Construct asymptotic two-sided and upper one-sided tests for an exponential rate. Evaluate the upper one-sided test of rate at most one for fifty observations with sample mean $0.98$.

Solution. Since $\widehat\lambda=1/\bar X$ has $\sqrt n$ asymptotic variance $\lambda^2$, under $\lambda=\lambda_0$,

$$
Z_n=\frac{\sqrt n(\widehat\lambda-\lambda_0)}{\widehat\lambda}
\Rightarrow N(0,1).
$$

Reject equality against a two-sided alternative when $|Z_n|>z_{1-\alpha/2}$. For $H_0:\lambda\le\lambda_0$ against $\lambda>\lambda_0$, reject when $Z_n>z_{1-\alpha}$. At the boundary the limiting rejection probability is $\alpha$; for any fixed $\lambda<\lambda_0$, the statistic tends to $-\infty$, so the rejection probability tends to zero.

For fifty interarrival observations with $\bar X=0.98$ and $\lambda_0=1$, $\widehat\lambda=1/0.98$ and $Z_n=\sqrt{50}(1-0.98)\approx0.1414$. The one-sided asymptotic p-value is $1-\Phi(0.1414)\approx0.444$. This does not reject at 5%. The calculation does not establish that the arrival rate is at most one; the observed difference is small relative to sampling uncertainty.

#### Problem PR-35: a Poisson interval and its inverted test

Source: [MIT 18.650, Problem Set 5, Problem 1, Questions 1–2](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/fc95ad9117580c46a8fddce6df5fb1ec_MIT18_650F16_PSet5.pdf).

Problem. For independent Poisson observations, construct a plug-in asymptotic confidence interval for the positive rate and the corresponding inverted two-sided test. Discuss the zero-count boundary.

Solution. For fixed $\lambda>0$, the CLT gives $\sqrt n(\bar X-\lambda)\Rightarrow N(0,\lambda)$ and the LLN gives $\bar X\to_p\lambda$. Hence

$$
C_n=\left[\bar X-z_{1-\alpha/2}\sqrt{\bar X/n},
\bar X+z_{1-\alpha/2}\sqrt{\bar X/n}\right]
$$

has asymptotic coverage $1-\alpha$. Intersecting it with $[0,\infty)$ preserves coverage for positive $\lambda$. Reject $H_0:\lambda=\lambda_0$ when $\lambda_0\notin C_n$, equivalently when $\sqrt n|\bar X-\lambda_0|/\sqrt{\bar X}>z_{1-\alpha/2}$ for positive $\bar X$. The zero-count event requires an explicit convention; its probability $e^{-n\lambda}$ vanishes at fixed positive $\lambda$, but a small-sample exact Poisson interval handles it more effectively.

### 7.5 Additional solved problems

#### Problem PR-36: likelihood-ratio and Wald tests for a Bernoulli parameter

Constructed problem. Test $p=p_0\in(0,1)$ with $k$ successes in $n$ trials.

Solution. Put $\widehat p=k/n$. The likelihood-ratio statistic is

$$
LR=2\left[
k\log\frac{\widehat p}{p_0}
+(n-k)\log\frac{1-\widehat p}{1-p_0}\right],
$$

with $0\log0$ interpreted by continuity. The Wald statistic is $W=n(\widehat p-p_0)^2/[\widehat p(1-\widehat p)]$ when the estimate is interior. Both approach $\chi_1^2$ under a fixed interior null. A second-order Taylor expansion of the Bernoulli log likelihood around $p_0$ explains their asymptotic agreement. Their finite-sample values differ, and the Wald denominator is zero at endpoint estimates while the LR remains well-defined.

#### Problem PR-37: an exact binomial tail test

Constructed problem. For $H_0:p\le1/2$ against $p>1/2$, explain an exact nonrandomized test using $K\sim\operatorname{Bin}(n,p)$.

Solution. Choose the smallest integer $c$ such that $P_{1/2}(K\ge c)\le\alpha$, and reject when $K\ge c$. Couple Bernoulli trials as $\mathbf1_{\{U_i\le p\}}$ using common independent uniforms. Their sum increases with $p$, so the upper-tail probability is largest within the null at $p=1/2$. Thus the test has level at most $\alpha$. The discrete tail may make the actual level strictly smaller.

#### Problem PR-38: a confidence interval after many comparisons

Constructed problem. Construct simultaneous intervals for $m$ parameters from individually valid intervals.

Solution. If each interval misses its parameter with probability at most $\alpha/m$, the union bound gives

$$
P(\text{at least one miss})\le\sum_{j=1}^m\alpha/m=\alpha.
$$

Thus using individual level $1-\alpha/m$ yields joint coverage at least $1-\alpha$. Independence among intervals is unnecessary. This Bonferroni argument can be conservative, especially when the estimates are strongly correlated.

#### Problem PR-39: size and power are different probabilities

Constructed problem. Suppose a level-5% test has power 80% at a specified alternative. In a population of repeated studies, 10% use that alternative and 90% use the null. Find the fraction of rejections arising under the alternative.

Solution. The model gives rejection probability $0.8(0.1)+0.05(0.9)=0.125$. The alternative-and-rejection probability is $0.08$. Their ratio is $0.08/0.125=0.64$. This calculation uses the extra mixture model for study types. The original test's 5% level alone is not a posterior error probability after rejection.

#### Problem PR-40: equivalence requires a specified margin

Constructed problem. Explain why failing to reject $\mu=0$ does not establish $|\mu|<\delta$.

Solution. A low-precision sample can produce a wide interval containing zero and values well beyond $\pm\delta$. To establish equivalence at a chosen margin, test both $H_{01}:\mu\le-\delta$ and $H_{02}:\mu\ge\delta$ and reject both using the appropriate one-sided procedures. Equivalently, in a regular normal setting, an associated interval must lie wholly inside $(-\delta,\delta)$. The margin must be specified in the units of the outcome before interpreting practical closeness.

## 8. Resampling and dependence

### 8.1 Theory

The nonparametric bootstrap treats the empirical distribution $F_n=n^{-1}\sum_i\delta_{X_i}$ as an estimate of the population distribution. Conditional on the data, draw a new sample of size $n$ independently from $F_n$ and recompute the estimator. Repeating this approximates a sampling distribution when the bootstrap is valid for the particular statistic and data-generating assumptions.

For a regular iid mean, the conditional bootstrap distribution of $\sqrt n(\bar X^*-\bar X)$ approximates that of $\sqrt n(\bar X-\mu)$. The centering is essential: bootstrap randomness is generated around the observed empirical mean.

A percentile interval uses bootstrap estimator quantiles directly. A basic interval reflects centered-error quantiles around the observed estimate. A studentized interval resamples a standardized statistic, requiring a standard-error estimate in each resample. These constructions need not give the same endpoints, and none is universally valid.

Under serial dependence, independent resampling destroys the sequence of lag relationships. Moving-block bootstrap samples consecutive blocks and joins them until the desired sample length is reached. Under appropriate stationarity, mixing, moment, and block-length conditions, this preserves enough local dependence for many smooth statistics. Typically the block length $\ell$ must grow while $\ell/n\to0$.

### 8.2 Derivations and proofs

#### Exact conditional variance of a bootstrap mean

Conditional on the observed data, each bootstrap draw has mean $\bar X$ and variance $v_n=n^{-1}\sum_i(X_i-\bar X)^2$. The draws are conditionally independent, so

$$
E_*[\bar X^*]=\bar X,\qquad
\operatorname{Var}_*(\bar X^*)=\frac{v_n}{n}.
$$

For iid finite-variance observations, $v_n\to_p\sigma^2$, giving the correct first-order variance scale. For dependent observations, $v_n$ estimates only lag-zero variance, while the mean depends on covariances at all lags.

#### Variance of a stationary sample mean

For a weakly stationary series with autocovariance $\gamma_k$,

$$
\operatorname{Var}(\bar X_n)
=\frac1{n^2}\sum_{t=1}^n\sum_{s=1}^n\gamma_{t-s}
=\frac1n\left[\gamma_0+
2\sum_{k=1}^{n-1}\left(1-\frac kn\right)\gamma_k\right].
$$

There are $n-k$ pairs at each positive lag. If autocovariances are absolutely summable, multiplying by $n$ and taking a limit gives the long-run variance $\Omega=\gamma_0+2\sum_{k\ge1}\gamma_k$. An iid bootstrap instead targets $\gamma_0$. Positive serial correlation can therefore cause substantial understatement of uncertainty.

#### The bootstrap probability is conditional on the data

Given observations $x_1,\ldots,x_n$, an ordinary bootstrap draw selects each observed value with probability $1/n$, independently across bootstrap positions. The starred expectation $E_*$ refers to this conditional resampling distribution, with the original data held fixed.

For the resampled mean,

$$
E_*[\bar X^*]=\bar x,\qquad
\operatorname{Var}_*(\bar X^*)=
\frac1{n^2}\sum_{i=1}^n(x_i-\bar x)^2.
$$

The second formula is the empirical-distribution variance divided by $n$. Its denominator differs slightly from using the unbiased sample variance divided by $n$. The bootstrap imitates sampling from the empirical distribution; it does not automatically apply every finite-sample bias correction.

A vector of bootstrap counts has a multinomial distribution. Counts $k_1,\ldots,k_n$ summing to $n$ have probability $n!/(k_1!\cdots k_n!)(1/n)^n$. This permits exact enumeration for very small samples, without simulation noise.

#### Example PR-E28: a small sample with degenerate resamples

Take data $(1,2,6)$. The mean is three and the sum of squared centered observations is fourteen. The conditional bootstrap variance of the mean is $14/9$. There are twenty-seven equally likely ordered resamples.

Three resamples contain only one repeated value: $(1,1,1)$, $(2,2,2)$, and $(6,6,6)$. Their total probability is $1/9$. Each has zero resampled standard deviation, so a studentized bootstrap statistic dividing by that standard deviation is undefined on a nonnegligible part of this tiny resampling space.

The bootstrap principle does not guarantee that every chosen statistic behaves regularly in every small sample. A percentile interval, a studentized interval, and an exact confidence procedure are different constructions with different assumptions.

#### Blocks retain only a finite stretch of dependence

For a moving-block bootstrap of a time series, form overlapping length-$\ell$ segments, draw whole segments, concatenate, and truncate to the desired sample length. Dependence within a selected segment is retained. Dependence across newly joined segments is broken.

Choosing $\ell=1$ reduces to independent resampling. Choosing one block as long as the sample supplies very little resampling variation. Standard consistency results for weakly dependent stationary sequences typically require $\ell\to\infty$ and $\ell/n\to0$, together with moment and dependence conditions. These rate conditions alone do not cover every process.

Circular blocks wrap the end of the series to the beginning to reduce unequal endpoint representation, while stationary-bootstrap variants use random block lengths. These are design choices, not universal repairs. Structural breaks, long memory, and nonstationary levels can require other methods.

Increasing the number of bootstrap replications reduces Monte Carlo noise in a quantile estimate. It does not repair an invalid resampling model. The sample size, dependence assumptions, statistic regularity, and block design determine whether the resampling distribution approximates the intended sampling distribution.

### 8.3 Worked examples

#### Example PR-E18: enumerating every resample

The observed sample is $(0,2)$, and the statistic is its mean.

Solution. The four equally likely ordered bootstrap samples are $(0,0),(0,2),(2,0),(2,2)$, with means $0,1,1,2$. Thus the bootstrap mean has probabilities $1/4,1/2,1/4$ at those values, conditional mean one, and variance $1/2$. The formula above gives $v_n=[(0-1)^2+(2-1)^2]/2=1$ and $v_n/n=1/2$.

#### Example PR-E19: serial correlation changes the uncertainty

For a stationary AR(1) with coefficient $\rho=0.8$, $\gamma_k=\gamma_0\rho^{|k|}$. The long-run variance is

$$
\Omega=\gamma_0\left(1+\frac{2\rho}{1-\rho}\right)
=9\gamma_0.
$$

For a long sample, the true sample-mean standard error is about three times the iid formula. Resampling individual observations would reproduce the marginal distribution and still miss this factor.

#### Example PR-E20: overlapping blocks

For data $(x_1,\ldots,x_6)$ and block length three, the moving blocks are $(x_1,x_2,x_3)$, $(x_2,x_3,x_4)$, $(x_3,x_4,x_5)$, $(x_4,x_5,x_6)$. Draw two blocks independently and concatenate them to obtain a six-observation resample. Within each block the original order is preserved; joins introduce artificial boundaries. Longer blocks preserve longer lags but leave fewer effectively distinct blocks.

### 8.4 Additional solved problems

#### Problem PR-41: percentile and basic intervals

Constructed problem. An estimate is $\widehat\theta=4$, and bootstrap 2.5% and 97.5% quantiles are $3.2$ and $5.1$. Compute two interval constructions.

Solution. The percentile interval is $[3.2,5.1]$. The basic interval reflects the estimator quantiles about the observed estimate:

$$
[2\widehat\theta-5.1,\ 2\widehat\theta-3.2]=[2.9,4.8].
$$

The difference arises because the bootstrap distribution is asymmetric around four. These endpoints alone do not prove 95% coverage; that requires validity of the resampling approximation for the estimator under the sampling model.

#### Problem PR-42: studentization

Constructed problem. If bootstrap quantiles of $T^*=(\widehat\theta^*-\widehat\theta)/\widehat{\operatorname{se}}^*$ are $q_L,q_U$, derive the interval.

Solution. Approximate the event $q_L\le(\widehat\theta-\theta)/\widehat{\operatorname{se}}\le q_U$. Solving for $\theta$ gives

$$
[\widehat\theta-q_U\widehat{\operatorname{se}},
\widehat\theta-q_L\widehat{\operatorname{se}}].
$$

The order reverses when solving the two inequalities because the target is subtracted in the numerator. Every resample needs its own usable standard error. With the two-point sample from Example PR-E18, some resamples have zero sample variance, so ordinary studentization can be undefined.

#### Problem PR-43: an effective sample size for a persistent mean

Constructed problem. For a long stationary AR(1) sample, find $n_{\mathrm{eff}}$ such that $\operatorname{Var}(\bar X_n)\approx\gamma_0/n_{\mathrm{eff}}$.

Solution. The long-run variance calculation gives $\operatorname{Var}(\bar X_n)\approx\gamma_0(1+\rho)/[n(1-\rho)]$. Therefore $n_{\mathrm{eff}}\approx n(1-\rho)/(1+\rho)$. At $\rho=0.8$ this is $n/9$. Negative correlation can give an effective size above $n$ for estimating a mean; the quantity summarizes a variance comparison and is not an actual count of independent observations.

#### Problem PR-44: a sample maximum and bootstrap failure

Constructed problem. Explain a difficulty with the ordinary bootstrap for the endpoint estimator from a uniform sample.

Solution. A bootstrap sample drawn from observed values can never exceed the observed maximum $M_n$. It also repeats that exact endpoint with conditional probability $1-(1-1/n)^n\to1-e^{-1}$. The true continuous-sample maximum has no atom at the old observed maximum. Thus the usual scaled maximum bootstrap distribution has a persistent atom absent from its target limit. Smooth-mean bootstrap intuition does not automatically transfer to boundary estimators.

#### Problem PR-45: pairs versus residual resampling

Constructed problem. In regression, distinguish resampling observation pairs from resampling residuals.

Solution. Pairs resampling draws $(X_i,Y_i)$ together and preserves their empirical joint relationship, appropriate to an iid random-design sampling model under suitable regularity. Residual resampling keeps regressors fixed, draws centered fitted residuals, and constructs $Y_i^*=X_i^\top\widehat\beta+u_i^*$. Its simplest form assumes identically distributed errors independent of the design. Under heteroskedasticity, pooling residuals across all regressor values can impose a false constant error distribution; a suitably designed wild bootstrap is a common alternative.

#### Problem PR-46: block length extremes

Constructed problem. Explain the consequences of choosing block length one or the full sample length.

Solution. Length one reduces moving-block resampling to independent observation resampling, losing all nonzero-lag dependence. Length $n$ leaves only one ordinary overlapping block, the original full sample, so repeated selection produces no useful sampling variability. Intermediate lengths trade dependence preservation against the number of distinct blocks. Conditions such as $\ell\to\infty$ and $\ell/n\to0$ exclude both extremes asymptotically.

#### Problem PR-47: nonstationarity is not repaired by blocking

Constructed problem. Suppose $X_t=t/n+\varepsilon_t$ with iid centered errors. Explain the problem with treating ordinary blocks as exchangeable.

Solution. Blocks near the beginning have systematically smaller means than blocks near the end. Randomly rearranging them changes the deterministic time pattern, so the resampling distribution no longer represents the original fixed trend. A procedure might first model and remove the trend, then resample a justified residual process, or use a resampling method designed for the relevant local nonstationarity. The series should not be declared stationary merely because each block is short.

#### Problem PR-48: simulation error in a bootstrap quantile

Constructed problem. Distinguish increasing the original sample size $n$ from increasing the number of bootstrap draws $B$.

Solution. Increasing $B$ approximates the conditional bootstrap distribution more accurately. It does not add new original information or remove a mismatch between that distribution and the true sampling law. For a fixed tail probability $p$, the Monte Carlo proportion of bootstrap statistics below a threshold has conditional standard deviation $\sqrt{p(1-p)/B}$. Quantile error also depends on the density near the target quantile. Increasing $n$ changes the statistical information; increasing $B$ reduces only simulation noise conditional on that information.

<!-- extended-practice-start -->

## 9. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 9.1 Bayes calculations with rare events

A positive result can be common among true cases and still have a modest posterior probability when true cases are rare. The denominator must include positive results from both populations. Computing their separate probability masses makes the calculation easier to inspect.

#### Problem PR-49: posterior after an imperfect signal, case 1

Let $D$ be an event with prior probability $\frac{1}{100}$. A signal satisfies $P(+\mid D)=\frac{9}{10}$ and $P(+\mid D^c)=\frac{1}{20}$. Find the posterior after a positive and after a negative signal.

Solution. Separate the positive-signal event according to whether $D$ occurred:

$$
P(D\cap +)=P(D)P(+\mid D)=\frac{9}{1000},\qquad
P(D^c\cap +)=P(D^c)P(+\mid D^c)=\frac{99}{2000}.
$$

The two events are disjoint and exhaust all positive signals. Therefore $P(+)=\frac{117}{2000}$. The conditional probability is

$$
P(D\mid +)=\frac{P(D\cap +)}{P(+)}=\frac{2}{13}
\approx 0.153846.
$$

The sensitivity $P(+\mid D)$ and this posterior condition in opposite directions. They answer different questions and need not be close.

For a negative signal, replace sensitivity and false-positive probability by their complements. The total negative probability is $1-P(+)=\frac{1883}{2000}$, while $P(D\cap -)=\frac{1}{1000}$. Hence

$$
P(D\mid -)=\frac{2}{1883}\approx 0.001062.
$$

An odds calculation gives an independent check of the positive posterior:

$$
\frac{P(D\mid +)}{1-P(D\mid +)}
=\frac{P(D)}{1-P(D)}\,
\frac{P(+\mid D)}{P(+\mid D^c)}.
$$

The first factor is prior odds and the second is the likelihood ratio of the observed signal. Multiplying and converting odds back to probability recovers $\frac{2}{13}$. Neither calculation assumes that successive signals would be conditionally independent. If several observations are combined, that additional assumption must be stated or their joint likelihood must be supplied.

#### Problem PR-50: posterior after an imperfect signal, case 2

Let $D$ be an event with prior probability $\frac{1}{1000}$. A signal satisfies $P(+\mid D)=\frac{99}{100}$ and $P(+\mid D^c)=\frac{1}{100}$. Find the posterior after a positive and after a negative signal.

Solution. Separate the positive-signal event according to whether $D$ occurred:

$$
P(D\cap +)=P(D)P(+\mid D)=\frac{99}{100000},\qquad
P(D^c\cap +)=P(D^c)P(+\mid D^c)=\frac{999}{100000}.
$$

The two events are disjoint and exhaust all positive signals. Therefore $P(+)=\frac{549}{50000}$. The conditional probability is

$$
P(D\mid +)=\frac{P(D\cap +)}{P(+)}=\frac{11}{122}
\approx 0.090164.
$$

The sensitivity $P(+\mid D)$ and this posterior condition in opposite directions. They answer different questions and need not be close.

For a negative signal, replace sensitivity and false-positive probability by their complements. The total negative probability is $1-P(+)=\frac{49451}{50000}$, while $P(D\cap -)=\frac{1}{100000}$. Hence

$$
P(D\mid -)=\frac{1}{98902}\approx 0.000010.
$$

An odds calculation gives an independent check of the positive posterior:

$$
\frac{P(D\mid +)}{1-P(D\mid +)}
=\frac{P(D)}{1-P(D)}\,
\frac{P(+\mid D)}{P(+\mid D^c)}.
$$

The first factor is prior odds and the second is the likelihood ratio of the observed signal. Multiplying and converting odds back to probability recovers $\frac{11}{122}$. Neither calculation assumes that successive signals would be conditionally independent. If several observations are combined, that additional assumption must be stated or their joint likelihood must be supplied.

#### Problem PR-51: posterior after an imperfect signal, case 3

Let $D$ be an event with prior probability $\frac{1}{5}$. A signal satisfies $P(+\mid D)=\frac{4}{5}$ and $P(+\mid D^c)=\frac{1}{10}$. Find the posterior after a positive and after a negative signal.

Solution. Separate the positive-signal event according to whether $D$ occurred:

$$
P(D\cap +)=P(D)P(+\mid D)=\frac{4}{25},\qquad
P(D^c\cap +)=P(D^c)P(+\mid D^c)=\frac{2}{25}.
$$

The two events are disjoint and exhaust all positive signals. Therefore $P(+)=\frac{6}{25}$. The conditional probability is

$$
P(D\mid +)=\frac{P(D\cap +)}{P(+)}=\frac{2}{3}
\approx 0.666667.
$$

The sensitivity $P(+\mid D)$ and this posterior condition in opposite directions. They answer different questions and need not be close.

For a negative signal, replace sensitivity and false-positive probability by their complements. The total negative probability is $1-P(+)=\frac{19}{25}$, while $P(D\cap -)=\frac{1}{25}$. Hence

$$
P(D\mid -)=\frac{1}{19}\approx 0.052632.
$$

An odds calculation gives an independent check of the positive posterior:

$$
\frac{P(D\mid +)}{1-P(D\mid +)}
=\frac{P(D)}{1-P(D)}\,
\frac{P(+\mid D)}{P(+\mid D^c)}.
$$

The first factor is prior odds and the second is the likelihood ratio of the observed signal. Multiplying and converting odds back to probability recovers $\frac{2}{3}$. Neither calculation assumes that successive signals would be conditionally independent. If several observations are combined, that additional assumption must be stated or their joint likelihood must be supplied.

#### Problem PR-52: posterior after an imperfect signal, case 4

Let $D$ be an event with prior probability $\frac{1}{50}$. A signal satisfies $P(+\mid D)=\frac{19}{20}$ and $P(+\mid D^c)=\frac{1}{100}$. Find the posterior after a positive and after a negative signal.

Solution. Separate the positive-signal event according to whether $D$ occurred:

$$
P(D\cap +)=P(D)P(+\mid D)=\frac{19}{1000},\qquad
P(D^c\cap +)=P(D^c)P(+\mid D^c)=\frac{49}{5000}.
$$

The two events are disjoint and exhaust all positive signals. Therefore $P(+)=\frac{18}{625}$. The conditional probability is

$$
P(D\mid +)=\frac{P(D\cap +)}{P(+)}=\frac{95}{144}
\approx 0.659722.
$$

The sensitivity $P(+\mid D)$ and this posterior condition in opposite directions. They answer different questions and need not be close.

For a negative signal, replace sensitivity and false-positive probability by their complements. The total negative probability is $1-P(+)=\frac{607}{625}$, while $P(D\cap -)=\frac{1}{1000}$. Hence

$$
P(D\mid -)=\frac{5}{4856}\approx 0.001030.
$$

An odds calculation gives an independent check of the positive posterior:

$$
\frac{P(D\mid +)}{1-P(D\mid +)}
=\frac{P(D)}{1-P(D)}\,
\frac{P(+\mid D)}{P(+\mid D^c)}.
$$

The first factor is prior odds and the second is the likelihood ratio of the observed signal. Multiplying and converting odds back to probability recovers $\frac{95}{144}$. Neither calculation assumes that successive signals would be conditionally independent. If several observations are combined, that additional assumption must be stated or their joint likelihood must be supplied.

### 9.2 Deriving a transformed distribution

A monotone transformation can be handled by the CDF without guessing a density. The support is determined first. Differentiating the CDF then gives the density, and moments can be checked independently by integrating with respect to the original uniform variable.

#### Problem PR-53: power of a uniform variable, case 1

Let $U\sim\operatorname{Uniform}(0,1)$ and $Y=U^p$ with $p=2$. Find the support, CDF, density, mean, variance, and median.

Solution. Because $p>0$, the transformation is strictly increasing and maps the unit interval into itself. Thus $Y$ has support $[0,1]$. For $0<y<1$,

$$
P(Y\le y)=P(U^p\le y)=P(U\le y^{1/p})=y^{1/p}.
$$

The complete CDF is zero for $y\le0$, equals $y^{1/p}$ inside the interval, and is one for $y\ge1$. There are no point masses at the endpoints because $U$ has none.

Differentiating on the interior gives

$$
f_Y(y)=\frac1p y^{1/p-1},\qquad 0<y<1.
$$

Its integral is $[y^{1/p}]_0^1=1$. A density can grow without bound near zero and still integrate to one; boundedness of a density is not a requirement.

For any real $k$ with $pk>-1$, direct integration over $U$ gives

$$
E[Y^k]=E[U^{pk}]=\int_0^1u^{pk}\,du=\frac1{pk+1}.
$$

In particular, $E[Y]=\frac{1}{3}$ and $E[Y^2]=\frac{1}{5}$. Subtracting the square of the mean gives $\operatorname{Var}(Y)=\frac{4}{45}$. The same values follow by integrating $y f_Y(y)$ and $y^2f_Y(y)$, which checks the Jacobian factor in the density.

The median $m$ satisfies $m^{1/p}=1/2$, so $m=2^{-p}=2^{-2}$. It generally differs from the mean. The transformation preserves the ordering of observations, but changes the distances between them and therefore changes averages and variances.

#### Problem PR-54: power of a uniform variable, case 2

Let $U\sim\operatorname{Uniform}(0,1)$ and $Y=U^p$ with $p=3$. Find the support, CDF, density, mean, variance, and median.

Solution. Because $p>0$, the transformation is strictly increasing and maps the unit interval into itself. Thus $Y$ has support $[0,1]$. For $0<y<1$,

$$
P(Y\le y)=P(U^p\le y)=P(U\le y^{1/p})=y^{1/p}.
$$

The complete CDF is zero for $y\le0$, equals $y^{1/p}$ inside the interval, and is one for $y\ge1$. There are no point masses at the endpoints because $U$ has none.

Differentiating on the interior gives

$$
f_Y(y)=\frac1p y^{1/p-1},\qquad 0<y<1.
$$

Its integral is $[y^{1/p}]_0^1=1$. A density can grow without bound near zero and still integrate to one; boundedness of a density is not a requirement.

For any real $k$ with $pk>-1$, direct integration over $U$ gives

$$
E[Y^k]=E[U^{pk}]=\int_0^1u^{pk}\,du=\frac1{pk+1}.
$$

In particular, $E[Y]=\frac{1}{4}$ and $E[Y^2]=\frac{1}{7}$. Subtracting the square of the mean gives $\operatorname{Var}(Y)=\frac{9}{112}$. The same values follow by integrating $y f_Y(y)$ and $y^2f_Y(y)$, which checks the Jacobian factor in the density.

The median $m$ satisfies $m^{1/p}=1/2$, so $m=2^{-p}=2^{-3}$. It generally differs from the mean. The transformation preserves the ordering of observations, but changes the distances between them and therefore changes averages and variances.

#### Problem PR-55: power of a uniform variable, case 3

Let $U\sim\operatorname{Uniform}(0,1)$ and $Y=U^p$ with $p=\frac{1}{2}$. Find the support, CDF, density, mean, variance, and median.

Solution. Because $p>0$, the transformation is strictly increasing and maps the unit interval into itself. Thus $Y$ has support $[0,1]$. For $0<y<1$,

$$
P(Y\le y)=P(U^p\le y)=P(U\le y^{1/p})=y^{1/p}.
$$

The complete CDF is zero for $y\le0$, equals $y^{1/p}$ inside the interval, and is one for $y\ge1$. There are no point masses at the endpoints because $U$ has none.

Differentiating on the interior gives

$$
f_Y(y)=\frac1p y^{1/p-1},\qquad 0<y<1.
$$

Its integral is $[y^{1/p}]_0^1=1$. A density can grow without bound near zero and still integrate to one; boundedness of a density is not a requirement.

For any real $k$ with $pk>-1$, direct integration over $U$ gives

$$
E[Y^k]=E[U^{pk}]=\int_0^1u^{pk}\,du=\frac1{pk+1}.
$$

In particular, $E[Y]=\frac{2}{3}$ and $E[Y^2]=\frac{1}{2}$. Subtracting the square of the mean gives $\operatorname{Var}(Y)=\frac{1}{18}$. The same values follow by integrating $y f_Y(y)$ and $y^2f_Y(y)$, which checks the Jacobian factor in the density.

The median $m$ satisfies $m^{1/p}=1/2$, so $m=2^{-p}=2^{-\frac{1}{2}}$. It generally differs from the mean. The transformation preserves the ordering of observations, but changes the distances between them and therefore changes averages and variances.

#### Problem PR-56: power of a uniform variable, case 4

Let $U\sim\operatorname{Uniform}(0,1)$ and $Y=U^p$ with $p=\frac{1}{3}$. Find the support, CDF, density, mean, variance, and median.

Solution. Because $p>0$, the transformation is strictly increasing and maps the unit interval into itself. Thus $Y$ has support $[0,1]$. For $0<y<1$,

$$
P(Y\le y)=P(U^p\le y)=P(U\le y^{1/p})=y^{1/p}.
$$

The complete CDF is zero for $y\le0$, equals $y^{1/p}$ inside the interval, and is one for $y\ge1$. There are no point masses at the endpoints because $U$ has none.

Differentiating on the interior gives

$$
f_Y(y)=\frac1p y^{1/p-1},\qquad 0<y<1.
$$

Its integral is $[y^{1/p}]_0^1=1$. A density can grow without bound near zero and still integrate to one; boundedness of a density is not a requirement.

For any real $k$ with $pk>-1$, direct integration over $U$ gives

$$
E[Y^k]=E[U^{pk}]=\int_0^1u^{pk}\,du=\frac1{pk+1}.
$$

In particular, $E[Y]=\frac{3}{4}$ and $E[Y^2]=\frac{3}{5}$. Subtracting the square of the mean gives $\operatorname{Var}(Y)=\frac{3}{80}$. The same values follow by integrating $y f_Y(y)$ and $y^2f_Y(y)$, which checks the Jacobian factor in the density.

The median $m$ satisfies $m^{1/p}=1/2$, so $m=2^{-p}=2^{-\frac{1}{3}}$. It generally differs from the mean. The transformation preserves the ordering of observations, but changes the distances between them and therefore changes averages and variances.

### 9.3 Conditional expectation from a joint table

For binary variables, a four-cell table determines every marginal and conditional quantity. The conditional mean $E[Y\mid X]$ is itself a random variable with one value for each possible $X$. The law of total variance separates average within-group uncertainty from variation between those conditional means.

#### Problem PR-57: binary table and total variance, case 1

The rows $X=0,1$ and columns $Y=0,1$ have joint probabilities $\begin{pmatrix}\frac{2}{5}&\frac{1}{10}\\\frac{1}{5}&\frac{3}{10}\end{pmatrix}$. Find $E[Y\mid X]$, covariance, and both parts of the law of total variance.

Solution. Sum each row to obtain $P(X=0)=\frac{1}{2}$ and $P(X=1)=\frac{1}{2}$. Sum the second column to obtain $P(Y=1)=\frac{2}{5}$. All row probabilities are positive, so both conditional distributions are defined.

Since $Y$ is binary, its conditional expectation is its conditional probability of one:

$$
E[Y\mid X=0]=\frac{P(X=0,Y=1)}{P(X=0)}=\frac{1}{5},\qquad
E[Y\mid X=1]=\frac{3}{5}.
$$

Thus the random variable $E[Y\mid X]$ equals the first number when $X=0$ and the second when $X=1$. Averaging with row probabilities gives $E[E[Y\mid X]]=\frac{2}{5}=E[Y]$, which checks the conditioning calculation.

Also $E[X]=\frac{1}{2}$ and $E[XY]=P(X=1,Y=1)=\frac{3}{10}$. Therefore

$$
\operatorname{Cov}(X,Y)=E[XY]-E[X]E[Y]=\frac{1}{10}.
$$

Within each row, $Y$ is Bernoulli with success probability $m_x$. Its conditional variance is $m_x(1-m_x)$. Averaging these two variances gives

$$
E[\operatorname{Var}(Y\mid X)]=\frac{1}{5}.
$$

The variance of the conditional mean, computed around its mean $E[Y]$, is

$$
\operatorname{Var}(E[Y\mid X])=\frac{1}{25}.
$$

Their sum is $\frac{6}{25}$, exactly $P(Y=1)(1-P(Y=1))$. The equality explains the decomposition: learning $X$ removes the between-group part of uncertainty on average, while the within-group part remains. Individual conditional variances need not all be smaller than the unconditional variance.

#### Problem PR-58: binary table and total variance, case 2

The rows $X=0,1$ and columns $Y=0,1$ have joint probabilities $\begin{pmatrix}\frac{1}{10}&\frac{3}{10}\\\frac{1}{5}&\frac{2}{5}\end{pmatrix}$. Find $E[Y\mid X]$, covariance, and both parts of the law of total variance.

Solution. Sum each row to obtain $P(X=0)=\frac{2}{5}$ and $P(X=1)=\frac{3}{5}$. Sum the second column to obtain $P(Y=1)=\frac{7}{10}$. All row probabilities are positive, so both conditional distributions are defined.

Since $Y$ is binary, its conditional expectation is its conditional probability of one:

$$
E[Y\mid X=0]=\frac{P(X=0,Y=1)}{P(X=0)}=\frac{3}{4},\qquad
E[Y\mid X=1]=\frac{2}{3}.
$$

Thus the random variable $E[Y\mid X]$ equals the first number when $X=0$ and the second when $X=1$. Averaging with row probabilities gives $E[E[Y\mid X]]=\frac{7}{10}=E[Y]$, which checks the conditioning calculation.

Also $E[X]=\frac{3}{5}$ and $E[XY]=P(X=1,Y=1)=\frac{2}{5}$. Therefore

$$
\operatorname{Cov}(X,Y)=E[XY]-E[X]E[Y]=-\frac{1}{50}.
$$

Within each row, $Y$ is Bernoulli with success probability $m_x$. Its conditional variance is $m_x(1-m_x)$. Averaging these two variances gives

$$
E[\operatorname{Var}(Y\mid X)]=\frac{5}{24}.
$$

The variance of the conditional mean, computed around its mean $E[Y]$, is

$$
\operatorname{Var}(E[Y\mid X])=\frac{1}{600}.
$$

Their sum is $\frac{21}{100}$, exactly $P(Y=1)(1-P(Y=1))$. The equality explains the decomposition: learning $X$ removes the between-group part of uncertainty on average, while the within-group part remains. Individual conditional variances need not all be smaller than the unconditional variance.

#### Problem PR-59: binary table and total variance, case 3

The rows $X=0,1$ and columns $Y=0,1$ have joint probabilities $\begin{pmatrix}\frac{1}{2}&\frac{1}{5}\\\frac{1}{10}&\frac{1}{5}\end{pmatrix}$. Find $E[Y\mid X]$, covariance, and both parts of the law of total variance.

Solution. Sum each row to obtain $P(X=0)=\frac{7}{10}$ and $P(X=1)=\frac{3}{10}$. Sum the second column to obtain $P(Y=1)=\frac{2}{5}$. All row probabilities are positive, so both conditional distributions are defined.

Since $Y$ is binary, its conditional expectation is its conditional probability of one:

$$
E[Y\mid X=0]=\frac{P(X=0,Y=1)}{P(X=0)}=\frac{2}{7},\qquad
E[Y\mid X=1]=\frac{2}{3}.
$$

Thus the random variable $E[Y\mid X]$ equals the first number when $X=0$ and the second when $X=1$. Averaging with row probabilities gives $E[E[Y\mid X]]=\frac{2}{5}=E[Y]$, which checks the conditioning calculation.

Also $E[X]=\frac{3}{10}$ and $E[XY]=P(X=1,Y=1)=\frac{1}{5}$. Therefore

$$
\operatorname{Cov}(X,Y)=E[XY]-E[X]E[Y]=\frac{2}{25}.
$$

Within each row, $Y$ is Bernoulli with success probability $m_x$. Its conditional variance is $m_x(1-m_x)$. Averaging these two variances gives

$$
E[\operatorname{Var}(Y\mid X)]=\frac{22}{105}.
$$

The variance of the conditional mean, computed around its mean $E[Y]$, is

$$
\operatorname{Var}(E[Y\mid X])=\frac{16}{525}.
$$

Their sum is $\frac{6}{25}$, exactly $P(Y=1)(1-P(Y=1))$. The equality explains the decomposition: learning $X$ removes the between-group part of uncertainty on average, while the within-group part remains. Individual conditional variances need not all be smaller than the unconditional variance.

#### Problem PR-60: binary table and total variance, case 4

The rows $X=0,1$ and columns $Y=0,1$ have joint probabilities $\begin{pmatrix}\frac{1}{5}&\frac{1}{5}\\\frac{1}{2}&\frac{1}{10}\end{pmatrix}$. Find $E[Y\mid X]$, covariance, and both parts of the law of total variance.

Solution. Sum each row to obtain $P(X=0)=\frac{2}{5}$ and $P(X=1)=\frac{3}{5}$. Sum the second column to obtain $P(Y=1)=\frac{3}{10}$. All row probabilities are positive, so both conditional distributions are defined.

Since $Y$ is binary, its conditional expectation is its conditional probability of one:

$$
E[Y\mid X=0]=\frac{P(X=0,Y=1)}{P(X=0)}=\frac{1}{2},\qquad
E[Y\mid X=1]=\frac{1}{6}.
$$

Thus the random variable $E[Y\mid X]$ equals the first number when $X=0$ and the second when $X=1$. Averaging with row probabilities gives $E[E[Y\mid X]]=\frac{3}{10}=E[Y]$, which checks the conditioning calculation.

Also $E[X]=\frac{3}{5}$ and $E[XY]=P(X=1,Y=1)=\frac{1}{10}$. Therefore

$$
\operatorname{Cov}(X,Y)=E[XY]-E[X]E[Y]=-\frac{2}{25}.
$$

Within each row, $Y$ is Bernoulli with success probability $m_x$. Its conditional variance is $m_x(1-m_x)$. Averaging these two variances gives

$$
E[\operatorname{Var}(Y\mid X)]=\frac{11}{60}.
$$

The variance of the conditional mean, computed around its mean $E[Y]$, is

$$
\operatorname{Var}(E[Y\mid X])=\frac{2}{75}.
$$

Their sum is $\frac{21}{100}$, exactly $P(Y=1)(1-P(Y=1))$. The equality explains the decomposition: learning $X$ removes the between-group part of uncertainty on average, while the within-group part remains. Individual conditional variances need not all be smaller than the unconditional variance.

### 9.4 Gaussian conditioning with correlation

For a jointly Gaussian pair, subtracting the linear prediction produces a Gaussian residual uncorrelated with the conditioning variable. Joint Gaussianity then upgrades zero covariance to independence. Without that distributional assumption, the same calculation gives a linear projection rather than necessarily the conditional mean.

#### Problem PR-61: conditional Gaussian prediction, case 1

Let $(X,Y)$ be jointly Gaussian with mean $(1,2)$ and covariance $\begin{pmatrix}4&3\\3&9\end{pmatrix}$. Find the law of $Y\mid X=3$ and the conditional probability that $Y>2$.

Solution. The covariance matrix is valid and nonsingular: its first diagonal entry is $4>0$ and its determinant is $27>0$. The correlation is $0.500000$, which lies strictly between minus one and one.

Choose $b=\operatorname{Cov}(X,Y)/\operatorname{Var}(X)=\frac{3}{4}$ and define

$$
R=Y-\mu_Y-b(X-\mu_X).
$$

Then $E[R]=0$ and $\operatorname{Cov}(R,X)=\operatorname{Cov}(Y,X)-b\operatorname{Var}(X)=0$. Because $(R,X)$ is jointly Gaussian, this zero covariance proves independence. Expanding its variance gives

$$
\operatorname{Var}(R)=\sigma_Y^2-\frac{\sigma_{XY}^2}{\sigma_X^2}
=\frac{27}{4}>0.
$$

Conditional on $X=3$, the linear term becomes a constant while $R$ keeps its unconditional Gaussian law. Therefore

$$
Y\mid X=3\sim
N\left(\frac{7}{2},\frac{27}{4}\right).
$$

Standardizing the requested threshold gives

$$
P(Y>2\mid X=3)
=1-\Phi\left(\frac{2-(\frac{7}{2})}{\sqrt{\frac{27}{4}}}\right)
\approx 0.718149.
$$

The conditional variance does not depend on the observed value of $X$ in this Gaussian model. It equals the unconditional variance minus an explained component. A singular covariance matrix would require separate treatment: the residual variance could be zero and the conditional law could be a point mass, making division by its standard deviation invalid.

#### Problem PR-62: conditional Gaussian prediction, case 2

Let $(X,Y)$ be jointly Gaussian with mean $(0,-1)$ and covariance $\begin{pmatrix}1&-1\\-1&4\end{pmatrix}$. Find the law of $Y\mid X=2$ and the conditional probability that $Y>-1$.

Solution. The covariance matrix is valid and nonsingular: its first diagonal entry is $1>0$ and its determinant is $3>0$. The correlation is $-0.500000$, which lies strictly between minus one and one.

Choose $b=\operatorname{Cov}(X,Y)/\operatorname{Var}(X)=-1$ and define

$$
R=Y-\mu_Y-b(X-\mu_X).
$$

Then $E[R]=0$ and $\operatorname{Cov}(R,X)=\operatorname{Cov}(Y,X)-b\operatorname{Var}(X)=0$. Because $(R,X)$ is jointly Gaussian, this zero covariance proves independence. Expanding its variance gives

$$
\operatorname{Var}(R)=\sigma_Y^2-\frac{\sigma_{XY}^2}{\sigma_X^2}
=3>0.
$$

Conditional on $X=2$, the linear term becomes a constant while $R$ keeps its unconditional Gaussian law. Therefore

$$
Y\mid X=2\sim
N\left(-3,3\right).
$$

Standardizing the requested threshold gives

$$
P(Y>-1\mid X=2)
=1-\Phi\left(\frac{-1-(-3)}{\sqrt{3}}\right)
\approx 0.124107.
$$

The conditional variance does not depend on the observed value of $X$ in this Gaussian model. It equals the unconditional variance minus an explained component. A singular covariance matrix would require separate treatment: the residual variance could be zero and the conditional law could be a point mass, making division by its standard deviation invalid.

#### Problem PR-63: conditional Gaussian prediction, case 3

Let $(X,Y)$ be jointly Gaussian with mean $(-2,1)$ and covariance $\begin{pmatrix}9&3\\3&4\end{pmatrix}$. Find the law of $Y\mid X=1$ and the conditional probability that $Y>1$.

Solution. The covariance matrix is valid and nonsingular: its first diagonal entry is $9>0$ and its determinant is $27>0$. The correlation is $0.500000$, which lies strictly between minus one and one.

Choose $b=\operatorname{Cov}(X,Y)/\operatorname{Var}(X)=\frac{1}{3}$ and define

$$
R=Y-\mu_Y-b(X-\mu_X).
$$

Then $E[R]=0$ and $\operatorname{Cov}(R,X)=\operatorname{Cov}(Y,X)-b\operatorname{Var}(X)=0$. Because $(R,X)$ is jointly Gaussian, this zero covariance proves independence. Expanding its variance gives

$$
\operatorname{Var}(R)=\sigma_Y^2-\frac{\sigma_{XY}^2}{\sigma_X^2}
=3>0.
$$

Conditional on $X=1$, the linear term becomes a constant while $R$ keeps its unconditional Gaussian law. Therefore

$$
Y\mid X=1\sim
N\left(2,3\right).
$$

Standardizing the requested threshold gives

$$
P(Y>1\mid X=1)
=1-\Phi\left(\frac{1-(2)}{\sqrt{3}}\right)
\approx 0.718149.
$$

The conditional variance does not depend on the observed value of $X$ in this Gaussian model. It equals the unconditional variance minus an explained component. A singular covariance matrix would require separate treatment: the residual variance could be zero and the conditional law could be a point mass, making division by its standard deviation invalid.

#### Problem PR-64: conditional Gaussian prediction, case 4

Let $(X,Y)$ be jointly Gaussian with mean $(3,0)$ and covariance $\begin{pmatrix}4&1\\1&1\end{pmatrix}$. Find the law of $Y\mid X=1$ and the conditional probability that $Y>0$.

Solution. The covariance matrix is valid and nonsingular: its first diagonal entry is $4>0$ and its determinant is $3>0$. The correlation is $0.500000$, which lies strictly between minus one and one.

Choose $b=\operatorname{Cov}(X,Y)/\operatorname{Var}(X)=\frac{1}{4}$ and define

$$
R=Y-\mu_Y-b(X-\mu_X).
$$

Then $E[R]=0$ and $\operatorname{Cov}(R,X)=\operatorname{Cov}(Y,X)-b\operatorname{Var}(X)=0$. Because $(R,X)$ is jointly Gaussian, this zero covariance proves independence. Expanding its variance gives

$$
\operatorname{Var}(R)=\sigma_Y^2-\frac{\sigma_{XY}^2}{\sigma_X^2}
=\frac{3}{4}>0.
$$

Conditional on $X=1$, the linear term becomes a constant while $R$ keeps its unconditional Gaussian law. Therefore

$$
Y\mid X=1\sim
N\left(-\frac{1}{2},\frac{3}{4}\right).
$$

Standardizing the requested threshold gives

$$
P(Y>0\mid X=1)
=1-\Phi\left(\frac{0-(-\frac{1}{2})}{\sqrt{\frac{3}{4}}}\right)
\approx 0.281851.
$$

The conditional variance does not depend on the observed value of $X$ in this Gaussian model. It equals the unconditional variance minus an explained component. A singular covariance matrix would require separate treatment: the residual variance could be zero and the conditional law could be a point mass, making division by its standard deviation invalid.

### 9.5 Finite-sample bounds and normal approximations

Chebyshev's inequality uses only a variance and gives a guaranteed upper bound. The central limit theorem gives an approximation whose finite-sample accuracy depends on the distribution and threshold. For Bernoulli samples, the binomial distribution allows an exact comparison with both.

#### Problem PR-65: Bernoulli mean error probability, case 1

Let $X_1,\ldots,X_n$ be independent Bernoulli variables with $p=\frac{1}{2}$ and $n=100$. Compare a variance bound, a normal approximation, and the exact probability of $|\bar X-p|\ge \frac{1}{10}$.

Solution. Independence removes cross-covariances, so

$$
E[\bar X]=p,\qquad
\operatorname{Var}(\bar X)=\frac{p(1-p)}n=\frac{1}{400}.
$$

Chebyshev's inequality therefore gives

$$
P(|\bar X-p|\ge\epsilon)
\le \min\left(1,\frac{p(1-p)}{n\epsilon^2}\right)
=\frac{1}{4}.
$$

The minimum with one records the trivial probability bound when the variance inequality is uninformative. This is a guaranteed finite-sample statement under the given independence and variance assumptions.

The normal approximation standardizes by $\sqrt{p(1-p)/n}$. The threshold is $z=2.000000$ standard deviations, giving the uncorrected approximation $2[1-\Phi(z)]\approx 0.045500$.

For an exact answer, let $S=n\bar X$. Then $S\sim\operatorname{Binomial}(n,p)$, and the event is

$$
S\le 40\quad\text{or}\quad S\ge 60.
$$

The floor and ceiling follow from the discrete support and the non-strict inequality in the event. Summing the binomial masses in these two tails gives approximately $0.056888$.

The exact value is below the Chebyshev bound, as required. It need not equal the normal approximation. A continuity correction can improve the latter by moving the integer cutoffs to half-integer boundaries. The central limit theorem itself is a limit statement for fixed standardized thresholds; it does not guarantee accurate probabilities for every far-tail event in a small sample.

#### Problem PR-66: Bernoulli mean error probability, case 2

Let $X_1,\ldots,X_n$ be independent Bernoulli variables with $p=\frac{1}{10}$ and $n=200$. Compare a variance bound, a normal approximation, and the exact probability of $|\bar X-p|\ge \frac{1}{20}$.

Solution. Independence removes cross-covariances, so

$$
E[\bar X]=p,\qquad
\operatorname{Var}(\bar X)=\frac{p(1-p)}n=\frac{9}{20000}.
$$

Chebyshev's inequality therefore gives

$$
P(|\bar X-p|\ge\epsilon)
\le \min\left(1,\frac{p(1-p)}{n\epsilon^2}\right)
=\frac{9}{50}.
$$

The minimum with one records the trivial probability bound when the variance inequality is uninformative. This is a guaranteed finite-sample statement under the given independence and variance assumptions.

The normal approximation standardizes by $\sqrt{p(1-p)/n}$. The threshold is $z=2.357023$ standard deviations, giving the uncorrected approximation $2[1-\Phi(z)]\approx 0.018422$.

For an exact answer, let $S=n\bar X$. Then $S\sim\operatorname{Binomial}(n,p)$, and the event is

$$
S\le 10\quad\text{or}\quad S\ge 30.
$$

The floor and ceiling follow from the discrete support and the non-strict inequality in the event. Summing the binomial masses in these two tails gives approximately $0.024398$.

The exact value is below the Chebyshev bound, as required. It need not equal the normal approximation. A continuity correction can improve the latter by moving the integer cutoffs to half-integer boundaries. The central limit theorem itself is a limit statement for fixed standardized thresholds; it does not guarantee accurate probabilities for every far-tail event in a small sample.

#### Problem PR-67: Bernoulli mean error probability, case 3

Let $X_1,\ldots,X_n$ be independent Bernoulli variables with $p=\frac{1}{4}$ and $n=80$. Compare a variance bound, a normal approximation, and the exact probability of $|\bar X-p|\ge \frac{1}{8}$.

Solution. Independence removes cross-covariances, so

$$
E[\bar X]=p,\qquad
\operatorname{Var}(\bar X)=\frac{p(1-p)}n=\frac{3}{1280}.
$$

Chebyshev's inequality therefore gives

$$
P(|\bar X-p|\ge\epsilon)
\le \min\left(1,\frac{p(1-p)}{n\epsilon^2}\right)
=\frac{3}{20}.
$$

The minimum with one records the trivial probability bound when the variance inequality is uninformative. This is a guaranteed finite-sample statement under the given independence and variance assumptions.

The normal approximation standardizes by $\sqrt{p(1-p)/n}$. The threshold is $z=2.581989$ standard deviations, giving the uncorrected approximation $2[1-\Phi(z)]\approx 0.009823$.

For an exact answer, let $S=n\bar X$. Then $S\sim\operatorname{Binomial}(n,p)$, and the event is

$$
S\le 10\quad\text{or}\quad S\ge 30.
$$

The floor and ceiling follow from the discrete support and the non-strict inequality in the event. Summing the binomial masses in these two tails gives approximately $0.013599$.

The exact value is below the Chebyshev bound, as required. It need not equal the normal approximation. A continuity correction can improve the latter by moving the integer cutoffs to half-integer boundaries. The central limit theorem itself is a limit statement for fixed standardized thresholds; it does not guarantee accurate probabilities for every far-tail event in a small sample.

#### Problem PR-68: Bernoulli mean error probability, case 4

Let $X_1,\ldots,X_n$ be independent Bernoulli variables with $p=\frac{4}{5}$ and $n=150$. Compare a variance bound, a normal approximation, and the exact probability of $|\bar X-p|\ge \frac{1}{10}$.

Solution. Independence removes cross-covariances, so

$$
E[\bar X]=p,\qquad
\operatorname{Var}(\bar X)=\frac{p(1-p)}n=\frac{2}{1875}.
$$

Chebyshev's inequality therefore gives

$$
P(|\bar X-p|\ge\epsilon)
\le \min\left(1,\frac{p(1-p)}{n\epsilon^2}\right)
=\frac{8}{75}.
$$

The minimum with one records the trivial probability bound when the variance inequality is uninformative. This is a guaranteed finite-sample statement under the given independence and variance assumptions.

The normal approximation standardizes by $\sqrt{p(1-p)/n}$. The threshold is $z=3.061862$ standard deviations, giving the uncorrected approximation $2[1-\Phi(z)]\approx 0.002200$.

For an exact answer, let $S=n\bar X$. Then $S\sim\operatorname{Binomial}(n,p)$, and the event is

$$
S\le 105\quad\text{or}\quad S\ge 135.
$$

The floor and ceiling follow from the discrete support and the non-strict inequality in the event. Summing the binomial masses in these two tails gives approximately $0.003070$.

The exact value is below the Chebyshev bound, as required. It need not equal the normal approximation. A continuity correction can improve the latter by moving the integer cutoffs to half-integer boundaries. The central limit theorem itself is a limit statement for fixed standardized thresholds; it does not guarantee accurate probabilities for every far-tail event in a small sample.

### 9.6 Exponential likelihood and an exact interval

An exponential rate is the reciprocal of a population mean. The maximum-likelihood estimate consequently involves the reciprocal of the sample mean. The sum has a gamma law, which also gives an exact confidence interval and avoids relying entirely on a large-sample approximation.

#### Problem PR-69: rate estimation from positive observations, case 1

The observations $1, 2, 3, 2, 2$ are modeled as independent exponential variables with density $f(x;\lambda)=\lambda e^{-\lambda x}$ for $x>0$ and $\lambda>0$. Find the MLE, estimated asymptotic standard error, and an exact 95% interval.

Solution. Write $S=\sum_i x_i=10$ and $n=5$. Ignoring factors independent of the rate, the likelihood and log-likelihood are

$$
L(\lambda)=\lambda^n e^{-\lambda S},\qquad
\ell(\lambda)=n\log\lambda-\lambda S.
$$

The score equation is $n/\lambda-S=0$, giving $\hat\lambda=n/S=\frac{1}{2}$. Since $\ell''(\lambda)=-n/\lambda^2<0$, the log-likelihood is strictly concave. It tends to minus infinity at both ends of the positive parameter line for these data, proving that this is the unique global maximizer.

The Fisher information is $I_n(\lambda)=n/\lambda^2$. Substituting the estimate gives asymptotic standard error

$$
\widehat{\operatorname{se}}(\hat\lambda)
=\frac{\hat\lambda}{\sqrt n}\approx 0.223607.
$$

This is a large-sample approximation, not an exact small-sample sampling standard deviation.

For an exact interval, the sum of $n$ independent exponential variables has a gamma distribution. Equivalently, $2\lambda S\sim\chi^2_{2n}$. Thus

$$
P\left(
\frac{\chi^2_{2n,0.025}}{2S}
\le\lambda\le
\frac{\chi^2_{2n,0.975}}{2S}
\right)=0.95.
$$

The observed interval is approximately $[0.162349,1.024159]$. Its endpoints are positive and it is generally asymmetric around the MLE. The coverage statement refers to repeated samples from this exponential model. Model misspecification, dependence, or selected observations would require a different sampling argument.

#### Problem PR-70: rate estimation from positive observations, case 2

The observations $\frac{1}{2}, 1, \frac{3}{2}, 2, 1, \frac{1}{2}$ are modeled as independent exponential variables with density $f(x;\lambda)=\lambda e^{-\lambda x}$ for $x>0$ and $\lambda>0$. Find the MLE, estimated asymptotic standard error, and an exact 95% interval.

Solution. Write $S=\sum_i x_i=\frac{13}{2}$ and $n=6$. Ignoring factors independent of the rate, the likelihood and log-likelihood are

$$
L(\lambda)=\lambda^n e^{-\lambda S},\qquad
\ell(\lambda)=n\log\lambda-\lambda S.
$$

The score equation is $n/\lambda-S=0$, giving $\hat\lambda=n/S=\frac{12}{13}$. Since $\ell''(\lambda)=-n/\lambda^2<0$, the log-likelihood is strictly concave. It tends to minus infinity at both ends of the positive parameter line for these data, proving that this is the unique global maximizer.

The Fisher information is $I_n(\lambda)=n/\lambda^2$. Substituting the estimate gives asymptotic standard error

$$
\widehat{\operatorname{se}}(\hat\lambda)
=\frac{\hat\lambda}{\sqrt n}\approx 0.376845.
$$

This is a large-sample approximation, not an exact small-sample sampling standard deviation.

For an exact interval, the sum of $n$ independent exponential variables has a gamma distribution. Equivalently, $2\lambda S\sim\chi^2_{2n}$. Thus

$$
P\left(
\frac{\chi^2_{2n,0.025}}{2S}
\le\lambda\le
\frac{\chi^2_{2n,0.975}}{2S}
\right)=0.95.
$$

The observed interval is approximately $[0.338753,1.795128]$. Its endpoints are positive and it is generally asymmetric around the MLE. The coverage statement refers to repeated samples from this exponential model. Model misspecification, dependence, or selected observations would require a different sampling argument.

#### Problem PR-71: rate estimation from positive observations, case 3

The observations $2, 4, 1, 3, 2, 5, 1, 2$ are modeled as independent exponential variables with density $f(x;\lambda)=\lambda e^{-\lambda x}$ for $x>0$ and $\lambda>0$. Find the MLE, estimated asymptotic standard error, and an exact 95% interval.

Solution. Write $S=\sum_i x_i=20$ and $n=8$. Ignoring factors independent of the rate, the likelihood and log-likelihood are

$$
L(\lambda)=\lambda^n e^{-\lambda S},\qquad
\ell(\lambda)=n\log\lambda-\lambda S.
$$

The score equation is $n/\lambda-S=0$, giving $\hat\lambda=n/S=\frac{2}{5}$. Since $\ell''(\lambda)=-n/\lambda^2<0$, the log-likelihood is strictly concave. It tends to minus infinity at both ends of the positive parameter line for these data, proving that this is the unique global maximizer.

The Fisher information is $I_n(\lambda)=n/\lambda^2$. Substituting the estimate gives asymptotic standard error

$$
\widehat{\operatorname{se}}(\hat\lambda)
=\frac{\hat\lambda}{\sqrt n}\approx 0.141421.
$$

This is a large-sample approximation, not an exact small-sample sampling standard deviation.

For an exact interval, the sum of $n$ independent exponential variables has a gamma distribution. Equivalently, $2\lambda S\sim\chi^2_{2n}$. Thus

$$
P\left(
\frac{\chi^2_{2n,0.025}}{2S}
\le\lambda\le
\frac{\chi^2_{2n,0.975}}{2S}
\right)=0.95.
$$

The observed interval is approximately $[0.172692,0.721134]$. Its endpoints are positive and it is generally asymmetric around the MLE. The coverage statement refers to repeated samples from this exponential model. Model misspecification, dependence, or selected observations would require a different sampling argument.

#### Problem PR-72: rate estimation from positive observations, case 4

The observations $\frac{1}{4}, \frac{3}{4}, 1, \frac{1}{2}, \frac{3}{2}$ are modeled as independent exponential variables with density $f(x;\lambda)=\lambda e^{-\lambda x}$ for $x>0$ and $\lambda>0$. Find the MLE, estimated asymptotic standard error, and an exact 95% interval.

Solution. Write $S=\sum_i x_i=4$ and $n=5$. Ignoring factors independent of the rate, the likelihood and log-likelihood are

$$
L(\lambda)=\lambda^n e^{-\lambda S},\qquad
\ell(\lambda)=n\log\lambda-\lambda S.
$$

The score equation is $n/\lambda-S=0$, giving $\hat\lambda=n/S=\frac{5}{4}$. Since $\ell''(\lambda)=-n/\lambda^2<0$, the log-likelihood is strictly concave. It tends to minus infinity at both ends of the positive parameter line for these data, proving that this is the unique global maximizer.

The Fisher information is $I_n(\lambda)=n/\lambda^2$. Substituting the estimate gives asymptotic standard error

$$
\widehat{\operatorname{se}}(\hat\lambda)
=\frac{\hat\lambda}{\sqrt n}\approx 0.559017.
$$

This is a large-sample approximation, not an exact small-sample sampling standard deviation.

For an exact interval, the sum of $n$ independent exponential variables has a gamma distribution. Equivalently, $2\lambda S\sim\chi^2_{2n}$. Thus

$$
P\left(
\frac{\chi^2_{2n,0.025}}{2S}
\le\lambda\le
\frac{\chi^2_{2n,0.975}}{2S}
\right)=0.95.
$$

The observed interval is approximately $[0.405872,2.560397]$. Its endpoints are positive and it is generally asymmetric around the MLE. The coverage statement refers to repeated samples from this exponential model. Model misspecification, dependence, or selected observations would require a different sampling argument.

### 9.7 A Gaussian mean with unknown variance

Replacing an unknown standard deviation by a sample estimate changes the exact pivot from normal to Student's t under Gaussian sampling. The denominator uses $n-1$ in the sample variance. The t critical value and the standard error play separate roles.

#### Problem PR-73: Student interval from a small sample, case 1

For independent $N(\mu,\sigma^2)$ observations $2, 4, 5, 7, 7$ with both parameters unknown, form an exact 95% confidence interval for $\mu$.

Solution. The sample contains $n=5$ observations with mean $\bar x=5$. Centering every observation at this mean gives the sum of squared deviations

$$
\sum_i(x_i-\bar x)^2=18.
$$

The sample variance and standard error of the sample mean are

$$
s^2=\frac1{n-1}\sum_i(x_i-\bar x)^2=\frac{9}{2},\qquad
\operatorname{se}(\bar x)=\frac s{\sqrt n}\approx 0.948683.
$$

Under independent Gaussian sampling, $\sqrt n(\bar X-\mu)/\sigma$ is standard normal, $(n-1)s^2/\sigma^2$ is chi-squared with $n-1$ degrees of freedom, and the two quantities are independent. Dividing the first by the square root of the second divided by its degrees of freedom gives

$$
\frac{\bar X-\mu}{s/\sqrt n}\sim t_{n-1}.
$$

The 97.5th percentile for $4$ degrees of freedom is approximately $2.776445$. Solving the central probability inequality for $\mu$ gives the observed interval

$$
\bar x\pm t_{n-1,0.975}\frac s{\sqrt n}
\approx [2.366033,7.633967].
$$

The interval is wider than one using the normal critical value with the same standard error, reflecting uncertainty about the variance. The exact t distribution depends on Gaussian sampling. For independent non-Gaussian observations with finite variance, a similar interval may be justified asymptotically, but this small sample alone does not establish that approximation. Dependence would also change the standard error of the mean.

#### Problem PR-74: Student interval from a small sample, case 2

For independent $N(\mu,\sigma^2)$ observations $10, 11, 9, 12, 13, 11$ with both parameters unknown, form an exact 95% confidence interval for $\mu$.

Solution. The sample contains $n=6$ observations with mean $\bar x=11$. Centering every observation at this mean gives the sum of squared deviations

$$
\sum_i(x_i-\bar x)^2=10.
$$

The sample variance and standard error of the sample mean are

$$
s^2=\frac1{n-1}\sum_i(x_i-\bar x)^2=2,\qquad
\operatorname{se}(\bar x)=\frac s{\sqrt n}\approx 0.577350.
$$

Under independent Gaussian sampling, $\sqrt n(\bar X-\mu)/\sigma$ is standard normal, $(n-1)s^2/\sigma^2$ is chi-squared with $n-1$ degrees of freedom, and the two quantities are independent. Dividing the first by the square root of the second divided by its degrees of freedom gives

$$
\frac{\bar X-\mu}{s/\sqrt n}\sim t_{n-1}.
$$

The 97.5th percentile for $5$ degrees of freedom is approximately $2.570582$. Solving the central probability inequality for $\mu$ gives the observed interval

$$
\bar x\pm t_{n-1,0.975}\frac s{\sqrt n}
\approx [9.515874,12.484126].
$$

The interval is wider than one using the normal critical value with the same standard error, reflecting uncertainty about the variance. The exact t distribution depends on Gaussian sampling. For independent non-Gaussian observations with finite variance, a similar interval may be justified asymptotically, but this small sample alone does not establish that approximation. Dependence would also change the standard error of the mean.

#### Problem PR-75: Student interval from a small sample, case 3

For independent $N(\mu,\sigma^2)$ observations $1, 1, 2, 3, 5, 6, 3, 3$ with both parameters unknown, form an exact 95% confidence interval for $\mu$.

Solution. The sample contains $n=8$ observations with mean $\bar x=3$. Centering every observation at this mean gives the sum of squared deviations

$$
\sum_i(x_i-\bar x)^2=22.
$$

The sample variance and standard error of the sample mean are

$$
s^2=\frac1{n-1}\sum_i(x_i-\bar x)^2=\frac{22}{7},\qquad
\operatorname{se}(\bar x)=\frac s{\sqrt n}\approx 0.626783.
$$

Under independent Gaussian sampling, $\sqrt n(\bar X-\mu)/\sigma$ is standard normal, $(n-1)s^2/\sigma^2$ is chi-squared with $n-1$ degrees of freedom, and the two quantities are independent. Dividing the first by the square root of the second divided by its degrees of freedom gives

$$
\frac{\bar X-\mu}{s/\sqrt n}\sim t_{n-1}.
$$

The 97.5th percentile for $7$ degrees of freedom is approximately $2.364624$. Solving the central probability inequality for $\mu$ gives the observed interval

$$
\bar x\pm t_{n-1,0.975}\frac s{\sqrt n}
\approx [1.517893,4.482107].
$$

The interval is wider than one using the normal critical value with the same standard error, reflecting uncertainty about the variance. The exact t distribution depends on Gaussian sampling. For independent non-Gaussian observations with finite variance, a similar interval may be justified asymptotically, but this small sample alone does not establish that approximation. Dependence would also change the standard error of the mean.

#### Problem PR-76: Student interval from a small sample, case 4

For independent $N(\mu,\sigma^2)$ observations $-2, 0, 1, 2, 4$ with both parameters unknown, form an exact 95% confidence interval for $\mu$.

Solution. The sample contains $n=5$ observations with mean $\bar x=1$. Centering every observation at this mean gives the sum of squared deviations

$$
\sum_i(x_i-\bar x)^2=20.
$$

The sample variance and standard error of the sample mean are

$$
s^2=\frac1{n-1}\sum_i(x_i-\bar x)^2=5,\qquad
\operatorname{se}(\bar x)=\frac s{\sqrt n}\approx 1.000000.
$$

Under independent Gaussian sampling, $\sqrt n(\bar X-\mu)/\sigma$ is standard normal, $(n-1)s^2/\sigma^2$ is chi-squared with $n-1$ degrees of freedom, and the two quantities are independent. Dividing the first by the square root of the second divided by its degrees of freedom gives

$$
\frac{\bar X-\mu}{s/\sqrt n}\sim t_{n-1}.
$$

The 97.5th percentile for $4$ degrees of freedom is approximately $2.776445$. Solving the central probability inequality for $\mu$ gives the observed interval

$$
\bar x\pm t_{n-1,0.975}\frac s{\sqrt n}
\approx [-1.776445,3.776445].
$$

The interval is wider than one using the normal critical value with the same standard error, reflecting uncertainty about the variance. The exact t distribution depends on Gaussian sampling. For independent non-Gaussian observations with finite variance, a similar interval may be justified asymptotically, but this small sample alone does not establish that approximation. Dependence would also change the standard error of the mean.

### 9.8 Dependence and the variance of an average

For a stationary process, a sample mean contains all pairwise covariances, not only individual variances. A positive autocorrelation can make a long sequence contain much less information about the mean than an equally long independent sample. Negative autocorrelation can have the opposite effect.

#### Problem PR-77: AR covariance and effective sample size, case 1

A stationary process has $\operatorname{Cov}(X_t,X_{t-k})=\gamma_0\rho^{|k|}$ with $\rho=\frac{4}{5}$ and $\gamma_0=1$. For $n=20$, find the exact variance of its sample mean and compare it with the independent and long-run formulas.

Solution. Start from the variance of a sum:

$$
\operatorname{Var}(\bar X_n)
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
\operatorname{Cov}(X_i,X_j).
$$

There are $n$ diagonal terms of size $\gamma_0$. For each lag $k\ge1$, there are $n-k$ pairs above the diagonal and the same number below it. Therefore

$$
\operatorname{Var}(\bar X_n)
=\frac{\gamma_0}n
\left[1+2\sum_{k=1}^{n-1}
\left(1-\frac kn\right)\rho^k\right]
\approx 0.35115292.
$$

The bracket equals approximately $7.023058$. An independence calculation would give $\gamma_0/n=0.05000000$, omitting every off-diagonal covariance.

Since $|\rho|<1$, the covariance series is absolutely summable. Its sum is the long-run variance

$$
\Omega=\gamma_0+2\sum_{k=1}^\infty\gamma_0\rho^k
=\gamma_0\frac{1+\rho}{1-\rho}=9.
$$

The large-sample approximation to the variance of the mean is $\Omega/n=0.45000000$. Its difference from the exact finite-sample result comes from the weights $1-k/n$ and the truncation at lag $n-1$.

One descriptive effective sample size solves $\gamma_0/n_{\mathrm{eff}}=\operatorname{Var}(\bar X_n)$, giving $n_{\mathrm{eff}}\approx 2.8478$. This number concerns estimation of this mean under this covariance model. It is not a universal replacement for sample size in every statistic or test. A dependence-aware bootstrap or long-run variance estimator must preserve or estimate the covariance terms that the independent calculation discards.

#### Problem PR-78: AR covariance and effective sample size, case 2

A stationary process has $\operatorname{Cov}(X_t,X_{t-k})=\gamma_0\rho^{|k|}$ with $\rho=\frac{1}{2}$ and $\gamma_0=4$. For $n=50$, find the exact variance of its sample mean and compare it with the independent and long-run formulas.

Solution. Start from the variance of a sum:

$$
\operatorname{Var}(\bar X_n)
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
\operatorname{Cov}(X_i,X_j).
$$

There are $n$ diagonal terms of size $\gamma_0$. For each lag $k\ge1$, there are $n-k$ pairs above the diagonal and the same number below it. Therefore

$$
\operatorname{Var}(\bar X_n)
=\frac{\gamma_0}n
\left[1+2\sum_{k=1}^{n-1}
\left(1-\frac kn\right)\rho^k\right]
\approx 0.23360000.
$$

The bracket equals approximately $2.920000$. An independence calculation would give $\gamma_0/n=0.08000000$, omitting every off-diagonal covariance.

Since $|\rho|<1$, the covariance series is absolutely summable. Its sum is the long-run variance

$$
\Omega=\gamma_0+2\sum_{k=1}^\infty\gamma_0\rho^k
=\gamma_0\frac{1+\rho}{1-\rho}=12.
$$

The large-sample approximation to the variance of the mean is $\Omega/n=0.24000000$. Its difference from the exact finite-sample result comes from the weights $1-k/n$ and the truncation at lag $n-1$.

One descriptive effective sample size solves $\gamma_0/n_{\mathrm{eff}}=\operatorname{Var}(\bar X_n)$, giving $n_{\mathrm{eff}}\approx 17.1233$. This number concerns estimation of this mean under this covariance model. It is not a universal replacement for sample size in every statistic or test. A dependence-aware bootstrap or long-run variance estimator must preserve or estimate the covariance terms that the independent calculation discards.

#### Problem PR-79: AR covariance and effective sample size, case 3

A stationary process has $\operatorname{Cov}(X_t,X_{t-k})=\gamma_0\rho^{|k|}$ with $\rho=-\frac{1}{2}$ and $\gamma_0=1$. For $n=20$, find the exact variance of its sample mean and compare it with the independent and long-run formulas.

Solution. Start from the variance of a sum:

$$
\operatorname{Var}(\bar X_n)
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
\operatorname{Cov}(X_i,X_j).
$$

There are $n$ diagonal terms of size $\gamma_0$. For each lag $k\ge1$, there are $n-k$ pairs above the diagonal and the same number below it. Therefore

$$
\operatorname{Var}(\bar X_n)
=\frac{\gamma_0}n
\left[1+2\sum_{k=1}^{n-1}
\left(1-\frac kn\right)\rho^k\right]
\approx 0.01777778.
$$

The bracket equals approximately $0.355556$. An independence calculation would give $\gamma_0/n=0.05000000$, omitting every off-diagonal covariance.

Since $|\rho|<1$, the covariance series is absolutely summable. Its sum is the long-run variance

$$
\Omega=\gamma_0+2\sum_{k=1}^\infty\gamma_0\rho^k
=\gamma_0\frac{1+\rho}{1-\rho}=\frac{1}{3}.
$$

The large-sample approximation to the variance of the mean is $\Omega/n=0.01666667$. Its difference from the exact finite-sample result comes from the weights $1-k/n$ and the truncation at lag $n-1$.

One descriptive effective sample size solves $\gamma_0/n_{\mathrm{eff}}=\operatorname{Var}(\bar X_n)$, giving $n_{\mathrm{eff}}\approx 56.2500$. This number concerns estimation of this mean under this covariance model. It is not a universal replacement for sample size in every statistic or test. A dependence-aware bootstrap or long-run variance estimator must preserve or estimate the covariance terms that the independent calculation discards.

#### Problem PR-80: AR covariance and effective sample size, case 4

A stationary process has $\operatorname{Cov}(X_t,X_{t-k})=\gamma_0\rho^{|k|}$ with $\rho=\frac{9}{10}$ and $\gamma_0=2$. For $n=100$, find the exact variance of its sample mean and compare it with the independent and long-run formulas.

Solution. Start from the variance of a sum:

$$
\operatorname{Var}(\bar X_n)
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
\operatorname{Cov}(X_i,X_j).
$$

There are $n$ diagonal terms of size $\gamma_0$. For each lag $k\ge1$, there are $n-k$ pairs above the diagonal and the same number below it. Therefore

$$
\operatorname{Var}(\bar X_n)
=\frac{\gamma_0}n
\left[1+2\sum_{k=1}^{n-1}
\left(1-\frac kn\right)\rho^k\right]
\approx 0.34400096.
$$

The bracket equals approximately $17.200048$. An independence calculation would give $\gamma_0/n=0.02000000$, omitting every off-diagonal covariance.

Since $|\rho|<1$, the covariance series is absolutely summable. Its sum is the long-run variance

$$
\Omega=\gamma_0+2\sum_{k=1}^\infty\gamma_0\rho^k
=\gamma_0\frac{1+\rho}{1-\rho}=38.
$$

The large-sample approximation to the variance of the mean is $\Omega/n=0.38000000$. Its difference from the exact finite-sample result comes from the weights $1-k/n$ and the truncation at lag $n-1$.

One descriptive effective sample size solves $\gamma_0/n_{\mathrm{eff}}=\operatorname{Var}(\bar X_n)$, giving $n_{\mathrm{eff}}\approx 5.8139$. This number concerns estimation of this mean under this covariance model. It is not a universal replacement for sample size in every statistic or test. A dependence-aware bootstrap or long-run variance estimator must preserve or estimate the covariance terms that the independent calculation discards.

<!-- further-distinct-problems-PR -->

### 9.9 Further distribution and inference problems

These additional problems are constructed for these notes.

#### Problem PR-81: a Poisson model with a random intensity

Let $\Lambda$ have gamma density $b^k\lambda^{k-1}e^{-b\lambda}/\Gamma(k)$ for $\lambda>0$, with $k,b>0$. Conditional on $\Lambda$, let $N$ be Poisson with mean $\Lambda$. Find the marginal mass function, mean, and variance of $N$.

Solution. Multiply the conditional Poisson mass by the mixing density and integrate:

$$
P(N=n)=\frac{b^k}{n!\Gamma(k)}
\int_0^\infty\lambda^{n+k-1}e^{-(b+1)\lambda}d\lambda
=\frac{\Gamma(n+k)}{n!\Gamma(k)}
\frac{b^k}{(b+1)^{n+k}}.
$$

The gamma-integral identity follows from substituting $u=(b+1)\lambda$ into its defining integral. The expression is a negative-binomial mass function, with this parameterization valid for noninteger positive $k$ as well.

Conditional mean and variance both equal $\Lambda$. The tower property and total variance give

$$
E[N]=E[\Lambda]=\frac kb,\qquad
\operatorname{Var}(N)=E[\Lambda]+\operatorname{Var}(\Lambda)
=\frac kb+\frac k{b^2}.
$$

The variance exceeds the mean because intensity variation adds an extra component. For $k=2,b=1$, the mean is two and variance four, while the probability of zero is $1/4$.

This is a generative mixture calculation. It does not require interpreting the mixing distribution as a subjective prior. The extra dispersion comes from heterogeneity in the conditional intensity, even though each conditional count is exactly Poisson.

#### Problem PR-82: multinomial covariance and a log ratio

For $n$ independent categorical observations with positive category probabilities $p_1,\ldots,p_d$, derive the asymptotic covariance of sample proportions and the delta-method variance of $\log(\hat p_1/\hat p_2)$.

Solution. Represent one observation by the indicator vector $I$ with exactly one entry equal to one. Then $E[I]=p$, each diagonal second moment is $p_i$, and off-diagonal products vanish. Therefore

$$
\Sigma=\operatorname{Cov}(I)=\operatorname{diag}(p)-pp^\top.
$$

The multivariate CLT gives $\sqrt n(\hat p-p)\Rightarrow N(0,\Sigma)$. The covariance is singular because the proportions always sum to one. For any vector $a$, its quadratic form is the variance of a discrete variable taking values $a_i$ with positive probabilities $p_i$. It vanishes exactly when all $a_i$ are equal, so the rank is $d-1$.

For $g(p)=\log p_1-\log p_2$, the gradient has entries $1/p_1,-1/p_2$ and zeros elsewhere. Substituting into the covariance form gives

$$
\nabla g^\top\Sigma\nabla g=\frac1{p_1}+\frac1{p_2}.
$$

The negative covariance between distinct proportions contributes to this result and must not be dropped.

For counts $(20,30,50)$, the estimated log ratio is $\log(2/3)$ and its asymptotic standard error is $\sqrt{(5+10/3)/100}=1/\sqrt{12}$. A log-scale interval is $\log(2/3)\pm1.96/\sqrt{12}$; exponentiating gives a positive interval for the probability ratio.

Zero observed counts require separate treatment because the log and plug-in gradient are undefined. The approximation is a fixed-interior-probability limit, rather than a uniform guarantee for categories becoming arbitrarily rare.

#### Problem PR-83: a variance-stabilizing transformation

For a Bernoulli sample mean with fixed $0<p<1$, find a differentiable increasing transformation whose square-root-sample-size asymptotic variance is one. Construct the resulting approximate interval and state its boundary limitation.

Solution. The sample proportion has asymptotic variance $p(1-p)$. The delta method gives transformed variance $g'(p)^2p(1-p)$. Requiring it to equal one and choosing the increasing sign gives $g'(p)=1/\sqrt{p(1-p)}$.

Integration gives $g(p)=2\arcsin\sqrt p$, up to an irrelevant additive constant. Thus

$$
\sqrt n\,[g(\hat p)-g(p)]\Rightarrow N(0,1).
$$

An approximate interval for $g(p)$ is $g(\hat p)\pm z_{1-\alpha/2}/\sqrt n$. Its parameter range is $[0,\pi]$, so intersect with that range before applying the increasing inverse $p=\sin^2(g/2)$.

For $n=100,\hat p=0.1$, the center is $2\arcsin\sqrt{0.1}\approx0.643501$. At 95%, the untruncated transformed endpoints are approximately $0.447501$ and $0.839501$. Applying the inverse gives a bounded probability interval rather than directly adding a symmetric probability-scale error.

The transformation derivative diverges at zero and one. Its fixed-interior delta-method derivation does not establish exact coverage at those boundaries or for very rare successes in a small sample. A bounded interval is not automatically an accurately calibrated interval. The variance-stabilizing property is an asymptotic calculation under the stated sampling model.

#### Problem PR-84: a normal likelihood ratio and Student's t

For an iid normal sample with unknown variance and nonzero residual sum of squares, derive the likelihood-ratio statistic for $H_0:\mu=\mu_0$. Relate its exact rejection rule to the two-sided t test.

Solution. Under the unrestricted model, the MLEs are $\hat\mu=\bar X$ and $\hat v=n^{-1}\sum_i(X_i-\bar X)^2$. Under the null, the variance MLE is $v_0=n^{-1}\sum_i(X_i-\mu_0)^2$.

The sum-of-squares decomposition gives $v_0=\hat v+(\bar X-\mu_0)^2$. At either fitted variance, the maximized log likelihood differs only through $-n\log v/2$. Therefore

$$
-2\log\Lambda=n\log\left(\frac{v_0}{\hat v}\right)
=n\log\left(1+\frac{t^2}{n-1}\right),
\qquad
t=\frac{\sqrt n(\bar X-\mu_0)}s.
$$

The statistic is strictly increasing in $|t|$. Since $t\sim t_{n-1}$ exactly under the Gaussian null, an exact size-$\alpha$ LR test rejects when

$$
-2\log\Lambda>
n\log\left(1+\frac{t_{n-1,1-\alpha/2}^2}{n-1}\right).
$$

This is exactly the usual two-sided t rejection region. Using the asymptotic chi-square critical value for the LR statistic is a different finite-sample calibration, though the two become close in large samples.

For $n=5,\bar x=5,s^2=4,\mu_0=4$, $t=\sqrt5/2$ and the LR statistic is $5\log(1+5/16)$. The t magnitude is below the 5% critical value $2.776$, so the exact test does not reject.

The relation was derived for positive fitted variance. An all-equal sample would require handling the boundary likelihood rather than inserting zero into these ratios.

#### Problem PR-85: an exact interval for a uniform endpoint

For independent uniform observations on $(0,\theta)$, let $M$ be the maximum. Construct a finite-sample interval with coverage $1-\alpha$ using only $M$ and derive its endpoint factor.

Solution. The pivot $U=M/\theta$ has CDF $u^n$ on $[0,1]$. Also $M\le\theta$ almost surely. Seek an interval $[M,cM]$ with $c>1$.

Its coverage is

$$
P_\theta(M\le\theta\le cM)
=P_\theta(M/\theta\ge1/c)
=1-c^{-n}.
$$

Setting this equal to $1-\alpha$ gives $c=\alpha^{-1/n}$. Thus

$$
\left[M,\frac{M}{\alpha^{1/n}}\right]
$$

has exact coverage for every positive endpoint under the iid uniform model.

For $n=10$ and $\alpha=0.05$, the factor is approximately $1.34928$. The lower endpoint is the observed maximum because every smaller parameter value assigns zero likelihood to the sample. The upper endpoint reflects the chance that all observations happen to remain below a fraction of the true endpoint.

The confidence statement holds for each fixed $\theta$ over repeated samples. It does not assign probability $1-\alpha$ to a fixed parameter after observing the interval. It also uses the endpoint model's exact pivot, not the regular Fisher-information or normal approximation that fails for this moving-support family.

#### Problem PR-86: familywise error without independent tests

For $m$ tests with valid null p-values, prove Bonferroni familywise error control and the Holm step-down rule. Apply both at level $0.05$ to ordered p-values $0.001,0.0054,0.02,0.08,0.1,0.2,0.3,0.4,0.6,0.8$.

Solution. Validity means $P(p_i\le t)\le t$ under its true null. Bonferroni rejects when $p_i\le\alpha/m$. If $m_0$ nulls are true, the union bound gives probability of any false rejection at most $m_0\alpha/m\le\alpha$. No independence was used.

Holm orders the p-values and compares the $j$th to $\alpha/(m-j+1)$, rejecting consecutively until the first failure. To prove control, consider the first true null in the ordered list, at position $j$. If any true null is rejected, that first one must pass its threshold.

There are at least $m_0$ remaining hypotheses at that position, so $\alpha/(m-j+1)\le\alpha/m_0$. A false rejection therefore implies that the minimum of the true-null p-values is at most $\alpha/m_0$. The union bound over the $m_0$ true nulls gives probability at most $\alpha$.

For the supplied ten tests, Bonferroni uses threshold $0.005$ and rejects only the first. Holm rejects the first at $0.005$, then the second at $0.05/9\approx0.005556$. The third fails threshold $0.05/8=0.00625$, so the procedure stops and rejects exactly two.

Familywise control concerns the probability of at least one false rejection. It is not the same target as controlling the expected proportion of false discoveries. The proof depends on valid individual p-values; adjustment cannot repair invalid sampling assumptions in the underlying tests.
