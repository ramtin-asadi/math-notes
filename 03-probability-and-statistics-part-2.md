# Probability and Statistics — Part 2

Probability describes a model before observations arrive. Statistics uses observations to learn about an unknown feature of a model. These notes develop both, starting with events and ending with inference for dependent observations.

Sources:

- Joe Blitzstein and Jessica Hwang, [Harvard Stat 110 and its probability text](https://stat110.hsites.harvard.edu/). The [strategic practice collection](https://stat110.hsites.harvard.edu/strategic-practice-problems) has its own section numbering within each practice sheet.
- Philippe Rigollet, [MIT 18.650, Statistics for Applications, Fall 2016](https://ocw.mit.edu/courses/18-650-statistics-for-applications-fall-2016/), especially Problem Sets 1–6 and 10.

The source links identify the selected problems; their statements are paraphrased. Additional problems and worked examples are constructed for these notes. Throughout, iid means independent and identically distributed. The notation $`N(\mu,\sigma^2)`$ uses variance as its second argument, and $`\mathop{\mathrm{Exp}} (\lambda)`$ uses a rate, with mean $`1/\lambda`$.

[Contents](README.md) · [Part 1](03-probability-and-statistics-part-1.md)

Contents:

- [8. Resampling and dependence](#8-resampling-and-dependence)
- [9. Further worked problems](#9-further-worked-problems)

## 8. Resampling and dependence

### 8.1 Theory

The nonparametric bootstrap treats the empirical distribution $`F_n=n^{-1}\sum_i\delta_{X_i}`$ as an estimate of the population distribution. Conditional on the data, draw a new sample of size $`n`$ independently from $`F_n`$ and recompute the estimator. Repeating this approximates a sampling distribution when the bootstrap is valid for the particular statistic and data-generating assumptions.

For a regular iid mean, the conditional bootstrap distribution of $`\sqrt n(\bar X^*-\bar X)`$ approximates that of $`\sqrt n(\bar X-\mu)`$. The centering is essential: bootstrap randomness is generated around the observed empirical mean.

A percentile interval uses bootstrap estimator quantiles directly. A basic interval reflects centered-error quantiles around the observed estimate. A studentized interval resamples a standardized statistic, requiring a standard-error estimate in each resample. These constructions need not give the same endpoints, and none is universally valid.

Under serial dependence, independent resampling destroys the sequence of lag relationships. Moving-block bootstrap samples consecutive blocks and joins them until the desired sample length is reached. Under appropriate stationarity, mixing, moment, and block-length conditions, this preserves enough local dependence for many smooth statistics. Typically the block length $`\ell`$ must grow while $`\ell/n\to0`$.

### 8.2 Derivations and proofs

#### Exact conditional variance of a bootstrap mean

Conditional on the observed data, each bootstrap draw has mean $`\bar X`$ and variance $`v_n=n^{-1}\sum_i(X_i-\bar X)^2`$. The draws are conditionally independent, so

```math
E_*[\bar X^*]=\bar X,\qquad
\mathop{\mathrm{Var}} _*(\bar X^*)=\frac{v_n}{n}.
```

For iid finite-variance observations, $`v_n\to_p\sigma^2`$, giving the correct first-order variance scale. For dependent observations, $`v_n`$ estimates only lag-zero variance, while the mean depends on covariances at all lags.

#### Variance of a stationary sample mean

For a weakly stationary series with autocovariance $`\gamma_k`$,

```math
\mathop{\mathrm{Var}} (\bar X_n)
=\frac1{n^2}\sum_{t=1}^n\sum_{s=1}^n\gamma_{t-s}
=\frac1n\left[\gamma_0+
2\sum_{k=1}^{n-1}\left(1-\frac kn\right)\gamma_k\right].
```

There are $`n-k`$ pairs at each positive lag. If autocovariances are absolutely summable, multiplying by $`n`$ and taking a limit gives the long-run variance $`\Omega=\gamma_0+2\sum_{k\ge1}\gamma_k`$. An iid bootstrap instead targets $`\gamma_0`$. Positive serial correlation can therefore cause substantial understatement of uncertainty.

#### The bootstrap probability is conditional on the data

Given observations $`x_1,\ldots,x_n`$, an ordinary bootstrap draw selects each observed value with probability $`1/n`$, independently across bootstrap positions. The starred expectation $`E_*`$ refers to this conditional resampling distribution, with the original data held fixed.

For the resampled mean,

```math
E_*[\bar X^*]=\bar x,\qquad
\mathop{\mathrm{Var}} _*(\bar X^*)=
\frac1{n^2}\sum_{i=1}^n(x_i-\bar x)^2.
```

The second formula is the empirical-distribution variance divided by $`n`$. Its denominator differs slightly from using the unbiased sample variance divided by $`n`$. The bootstrap imitates sampling from the empirical distribution; it does not automatically apply every finite-sample bias correction.

A vector of bootstrap counts has a multinomial distribution. Counts $`k_1,\ldots,k_n`$ summing to $`n`$ have probability $`n!/(k_1!\cdots k_n!)(1/n)^n`$. This permits exact enumeration for very small samples, without simulation noise.

#### Example PR-E28: a small sample with degenerate resamples

Take data $`(1,2,6)`$. The mean is three and the sum of squared centered observations is fourteen. The conditional bootstrap variance of the mean is $`14/9`$. There are twenty-seven equally likely ordered resamples.

Three resamples contain only one repeated value: $`(1,1,1)`$, $`(2,2,2)`$, and $`(6,6,6)`$. Their total probability is $`1/9`$. Each has zero resampled standard deviation, so a studentized bootstrap statistic dividing by that standard deviation is undefined on a nonnegligible part of this tiny resampling space.

The bootstrap principle does not guarantee that every chosen statistic behaves regularly in every small sample. A percentile interval, a studentized interval, and an exact confidence procedure are different constructions with different assumptions.

#### Blocks retain only a finite stretch of dependence

For a moving-block bootstrap of a time series, form overlapping length-$`\ell`$ segments, draw whole segments, concatenate, and truncate to the desired sample length. Dependence within a selected segment is retained. Dependence across newly joined segments is broken.

Choosing $`\ell=1`$ reduces to independent resampling. Choosing one block as long as the sample supplies very little resampling variation. Standard consistency results for weakly dependent stationary sequences typically require $`\ell\to\infty`$ and $`\ell/n\to0`$, together with moment and dependence conditions. These rate conditions alone do not cover every process.

Circular blocks wrap the end of the series to the beginning to reduce unequal endpoint representation, while stationary-bootstrap variants use random block lengths. These are design choices, not universal repairs. Structural breaks, long memory, and nonstationary levels can require other methods.

Increasing the number of bootstrap replications reduces Monte Carlo noise in a quantile estimate. It does not repair an invalid resampling model. The sample size, dependence assumptions, statistic regularity, and block design determine whether the resampling distribution approximates the intended sampling distribution.

### 8.3 Worked examples

#### Example PR-E18: enumerating every resample

The observed sample is $`(0,2)`$, and the statistic is its mean.

Solution. The four equally likely ordered bootstrap samples are $`(0,0),(0,2),(2,0),(2,2)`$, with means $`0,1,1,2`$. Thus the bootstrap mean has probabilities $`1/4,1/2,1/4`$ at those values, conditional mean one, and variance $`1/2`$. The formula above gives $`v_n=[(0-1)^2+(2-1)^2]/2=1`$ and $`v_n/n=1/2`$.

#### Example PR-E19: serial correlation changes the uncertainty

For a stationary AR(1) with coefficient $`\rho=0.8`$, $`\gamma_k=\gamma_0\rho^{\vert k\vert }`$. The long-run variance is

```math
\Omega=\gamma_0\left(1+\frac{2\rho}{1-\rho}\right)
=9\gamma_0.
```

For a long sample, the true sample-mean standard error is about three times the iid formula. Resampling individual observations would reproduce the marginal distribution and still miss this factor.

#### Example PR-E20: overlapping blocks

For data $`(x_1,\ldots,x_6)`$ and block length three, the moving blocks are $`(x_1,x_2,x_3)`$, $`(x_2,x_3,x_4)`$, $`(x_3,x_4,x_5)`$, $`(x_4,x_5,x_6)`$. Draw two blocks independently and concatenate them to obtain a six-observation resample. Within each block the original order is preserved; joins introduce artificial boundaries. Longer blocks preserve longer lags but leave fewer effectively distinct blocks.

### 8.4 Additional solved problems

#### Problem PR-41: percentile and basic intervals

Constructed problem. An estimate is $`\widehat\theta=4`$, and bootstrap 2.5% and 97.5% quantiles are $`3.2`$ and $`5.1`$. Compute two interval constructions.

Solution. The percentile interval is $`[3.2,5.1]`$. The basic interval reflects the estimator quantiles about the observed estimate:

```math
[2\widehat\theta-5.1,\ 2\widehat\theta-3.2]=[2.9,4.8].
```

The difference arises because the bootstrap distribution is asymmetric around four. These endpoints alone do not prove 95% coverage; that requires validity of the resampling approximation for the estimator under the sampling model.

#### Problem PR-42: studentization

Constructed problem. If bootstrap quantiles of $`T^*=(\widehat\theta^*-\widehat\theta)/\widehat{\mathop{\mathrm{se}} }^*`$ are $`q_L,q_U`$, derive the interval.

Solution. Approximate the event $`q_L\le(\widehat\theta-\theta)/\widehat{\mathop{\mathrm{se}} }\le q_U`$. Solving for $`\theta`$ gives

```math
[\widehat\theta-q_U\widehat{\mathop{\mathrm{se}} },
\widehat\theta-q_L\widehat{\mathop{\mathrm{se}} }].
```

The order reverses when solving the two inequalities because the target is subtracted in the numerator. Every resample needs its own usable standard error. With the two-point sample from Example PR-E18, some resamples have zero sample variance, so ordinary studentization can be undefined.

#### Problem PR-43: an effective sample size for a persistent mean

Constructed problem. For a long stationary AR(1) sample, find $`n_{\mathrm{eff}}`$ such that $`\mathop{\mathrm{Var}} (\bar X_n)\approx\gamma_0/n_{\mathrm{eff}}`$.

Solution. The long-run variance calculation gives $`\mathop{\mathrm{Var}} (\bar X_n)\approx\gamma_0(1+\rho)/[n(1-\rho)]`$. Therefore $`n_{\mathrm{eff}}\approx n(1-\rho)/(1+\rho)`$. At $`\rho=0.8`$ this is $`n/9`$. Negative correlation can give an effective size above $`n`$ for estimating a mean; the quantity summarizes a variance comparison and is not an actual count of independent observations.

#### Problem PR-44: a sample maximum and bootstrap failure

Constructed problem. Explain a difficulty with the ordinary bootstrap for the endpoint estimator from a uniform sample.

Solution. A bootstrap sample drawn from observed values can never exceed the observed maximum $`M_n`$. It also repeats that exact endpoint with conditional probability $`1-(1-1/n)^n\to1-e^{-1}`$. The true continuous-sample maximum has no atom at the old observed maximum. Thus the usual scaled maximum bootstrap distribution has a persistent atom absent from its target limit. Smooth-mean bootstrap intuition does not automatically transfer to boundary estimators.

#### Problem PR-45: pairs versus residual resampling

Constructed problem. In regression, distinguish resampling observation pairs from resampling residuals.

Solution. Pairs resampling draws $`(X_i,Y_i)`$ together and preserves their empirical joint relationship, appropriate to an iid random-design sampling model under suitable regularity. Residual resampling keeps regressors fixed, draws centered fitted residuals, and constructs $`Y_i^*=X_i^\top\widehat\beta+u_i^*`$. Its simplest form assumes identically distributed errors independent of the design. Under heteroskedasticity, pooling residuals across all regressor values can impose a false constant error distribution; a suitably designed wild bootstrap is a common alternative.

#### Problem PR-46: block length extremes

Constructed problem. Explain the consequences of choosing block length one or the full sample length.

Solution. Length one reduces moving-block resampling to independent observation resampling, losing all nonzero-lag dependence. Length $`n`$ leaves only one ordinary overlapping block, the original full sample, so repeated selection produces no useful sampling variability. Intermediate lengths trade dependence preservation against the number of distinct blocks. Conditions such as $`\ell\to\infty`$ and $`\ell/n\to0`$ exclude both extremes asymptotically.

#### Problem PR-47: nonstationarity is not repaired by blocking

Constructed problem. Suppose $`X_t=t/n+\varepsilon_t`$ with iid centered errors. Explain the problem with treating ordinary blocks as exchangeable.

Solution. Blocks near the beginning have systematically smaller means than blocks near the end. Randomly rearranging them changes the deterministic time pattern, so the resampling distribution no longer represents the original fixed trend. A procedure might first model and remove the trend, then resample a justified residual process, or use a resampling method designed for the relevant local nonstationarity. The series should not be declared stationary merely because each block is short.

#### Problem PR-48: simulation error in a bootstrap quantile

Constructed problem. Distinguish increasing the original sample size $`n`$ from increasing the number of bootstrap draws $`B`$.

Solution. Increasing $`B`$ approximates the conditional bootstrap distribution more accurately. It does not add new original information or remove a mismatch between that distribution and the true sampling law. For a fixed tail probability $`p`$, the Monte Carlo proportion of bootstrap statistics below a threshold has conditional standard deviation $`\sqrt{p(1-p)/B}`$. Quantile error also depends on the density near the target quantile. Increasing $`n`$ changes the statistical information; increasing $`B`$ reduces only simulation noise conditional on that information.

<!-- extended-practice-start -->

## 9. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 9.1 Bayes calculations with rare events

A positive result can be common among true cases and still have a modest posterior probability when true cases are rare. The denominator must include positive results from both populations. Computing their separate probability masses makes the calculation easier to inspect.

#### Problem PR-49: posterior after an imperfect signal, case 1

Let $`D`$ be an event with prior probability $`\frac{1}{100}`$. A signal satisfies $`P(+\mid D)=\frac{9}{10}`$ and $`P(+\mid D^c)=\frac{1}{20}`$. Find the posterior after a positive and after a negative signal.

Solution. Separate the positive-signal event according to whether $`D`$ occurred:

```math
P(D\cap +)=P(D)P(+\mid D)=\frac{9}{1000},\qquad
P(D^c\cap +)=P(D^c)P(+\mid D^c)=\frac{99}{2000}.
```

The two events are disjoint and exhaust all positive signals. Therefore $`P(+)=\frac{117}{2000}`$. The conditional probability is

```math
P(D\mid +)=\frac{P(D\cap +)}{P(+)}=\frac{2}{13}
\approx 0.153846.
```

The sensitivity $`P(+\mid D)`$ and this posterior condition in opposite directions. They answer different questions and need not be close.

For a negative signal, replace sensitivity and false-positive probability by their complements. The total negative probability is $`1-P(+)=\frac{1883}{2000}`$, while $`P(D\cap -)=\frac{1}{1000}`$. Hence

```math
P(D\mid -)=\frac{2}{1883}\approx 0.001062.
```

An odds calculation gives an independent check of the positive posterior:

```math
\frac{P(D\mid +)}{1-P(D\mid +)}
=\frac{P(D)}{1-P(D)}\,
\frac{P(+\mid D)}{P(+\mid D^c)}.
```

The first factor is prior odds and the second is the likelihood ratio of the observed signal. Multiplying and converting odds back to probability recovers $`\frac{2}{13}`$. Neither calculation assumes that successive signals would be conditionally independent. If several observations are combined, that additional assumption must be stated or their joint likelihood must be supplied.

#### Problem PR-50: posterior after an imperfect signal, case 2

Let $`D`$ be an event with prior probability $`\frac{1}{1000}`$. A signal satisfies $`P(+\mid D)=\frac{99}{100}`$ and $`P(+\mid D^c)=\frac{1}{100}`$. Find the posterior after a positive and after a negative signal.

Solution. Separate the positive-signal event according to whether $`D`$ occurred:

```math
P(D\cap +)=P(D)P(+\mid D)=\frac{99}{100000},\qquad
P(D^c\cap +)=P(D^c)P(+\mid D^c)=\frac{999}{100000}.
```

The two events are disjoint and exhaust all positive signals. Therefore $`P(+)=\frac{549}{50000}`$. The conditional probability is

```math
P(D\mid +)=\frac{P(D\cap +)}{P(+)}=\frac{11}{122}
\approx 0.090164.
```

The sensitivity $`P(+\mid D)`$ and this posterior condition in opposite directions. They answer different questions and need not be close.

For a negative signal, replace sensitivity and false-positive probability by their complements. The total negative probability is $`1-P(+)=\frac{49451}{50000}`$, while $`P(D\cap -)=\frac{1}{100000}`$. Hence

```math
P(D\mid -)=\frac{1}{98902}\approx 0.000010.
```

An odds calculation gives an independent check of the positive posterior:

```math
\frac{P(D\mid +)}{1-P(D\mid +)}
=\frac{P(D)}{1-P(D)}\,
\frac{P(+\mid D)}{P(+\mid D^c)}.
```

The first factor is prior odds and the second is the likelihood ratio of the observed signal. Multiplying and converting odds back to probability recovers $`\frac{11}{122}`$. Neither calculation assumes that successive signals would be conditionally independent. If several observations are combined, that additional assumption must be stated or their joint likelihood must be supplied.

#### Problem PR-51: posterior after an imperfect signal, case 3

Let $`D`$ be an event with prior probability $`\frac{1}{5}`$. A signal satisfies $`P(+\mid D)=\frac{4}{5}`$ and $`P(+\mid D^c)=\frac{1}{10}`$. Find the posterior after a positive and after a negative signal.

Solution. Separate the positive-signal event according to whether $`D`$ occurred:

```math
P(D\cap +)=P(D)P(+\mid D)=\frac{4}{25},\qquad
P(D^c\cap +)=P(D^c)P(+\mid D^c)=\frac{2}{25}.
```

The two events are disjoint and exhaust all positive signals. Therefore $`P(+)=\frac{6}{25}`$. The conditional probability is

```math
P(D\mid +)=\frac{P(D\cap +)}{P(+)}=\frac{2}{3}
\approx 0.666667.
```

The sensitivity $`P(+\mid D)`$ and this posterior condition in opposite directions. They answer different questions and need not be close.

For a negative signal, replace sensitivity and false-positive probability by their complements. The total negative probability is $`1-P(+)=\frac{19}{25}`$, while $`P(D\cap -)=\frac{1}{25}`$. Hence

```math
P(D\mid -)=\frac{1}{19}\approx 0.052632.
```

An odds calculation gives an independent check of the positive posterior:

```math
\frac{P(D\mid +)}{1-P(D\mid +)}
=\frac{P(D)}{1-P(D)}\,
\frac{P(+\mid D)}{P(+\mid D^c)}.
```

The first factor is prior odds and the second is the likelihood ratio of the observed signal. Multiplying and converting odds back to probability recovers $`\frac{2}{3}`$. Neither calculation assumes that successive signals would be conditionally independent. If several observations are combined, that additional assumption must be stated or their joint likelihood must be supplied.

#### Problem PR-52: posterior after an imperfect signal, case 4

Let $`D`$ be an event with prior probability $`\frac{1}{50}`$. A signal satisfies $`P(+\mid D)=\frac{19}{20}`$ and $`P(+\mid D^c)=\frac{1}{100}`$. Find the posterior after a positive and after a negative signal.

Solution. Separate the positive-signal event according to whether $`D`$ occurred:

```math
P(D\cap +)=P(D)P(+\mid D)=\frac{19}{1000},\qquad
P(D^c\cap +)=P(D^c)P(+\mid D^c)=\frac{49}{5000}.
```

The two events are disjoint and exhaust all positive signals. Therefore $`P(+)=\frac{18}{625}`$. The conditional probability is

```math
P(D\mid +)=\frac{P(D\cap +)}{P(+)}=\frac{95}{144}
\approx 0.659722.
```

The sensitivity $`P(+\mid D)`$ and this posterior condition in opposite directions. They answer different questions and need not be close.

For a negative signal, replace sensitivity and false-positive probability by their complements. The total negative probability is $`1-P(+)=\frac{607}{625}`$, while $`P(D\cap -)=\frac{1}{1000}`$. Hence

```math
P(D\mid -)=\frac{5}{4856}\approx 0.001030.
```

An odds calculation gives an independent check of the positive posterior:

```math
\frac{P(D\mid +)}{1-P(D\mid +)}
=\frac{P(D)}{1-P(D)}\,
\frac{P(+\mid D)}{P(+\mid D^c)}.
```

The first factor is prior odds and the second is the likelihood ratio of the observed signal. Multiplying and converting odds back to probability recovers $`\frac{95}{144}`$. Neither calculation assumes that successive signals would be conditionally independent. If several observations are combined, that additional assumption must be stated or their joint likelihood must be supplied.

### 9.2 Deriving a transformed distribution

A monotone transformation can be handled by the CDF without guessing a density. The support is determined first. Differentiating the CDF then gives the density, and moments can be checked independently by integrating with respect to the original uniform variable.

#### Problem PR-53: power of a uniform variable, case 1

Let $`U\sim\mathop{\mathrm{Uniform}} (0,1)`$ and $`Y=U^p`$ with $`p=2`$. Find the support, CDF, density, mean, variance, and median.

Solution. Because $`p\gt 0`$, the transformation is strictly increasing and maps the unit interval into itself. Thus $`Y`$ has support $`[0,1]`$. For $`0\lt y\lt 1`$,

```math
P(Y\le y)=P(U^p\le y)=P(U\le y^{1/p})=y^{1/p}.
```

The complete CDF is zero for $`y\le0`$, equals $`y^{1/p}`$ inside the interval, and is one for $`y\ge1`$. There are no point masses at the endpoints because $`U`$ has none.

Differentiating on the interior gives

```math
f_Y(y)=\frac1p y^{1/p-1},\qquad 0\lt y\lt 1.
```

Its integral is $`[y^{1/p}]_0^1=1`$. A density can grow without bound near zero and still integrate to one; boundedness of a density is not a requirement.

For any real $`k`$ with $`pk\gt -1`$, direct integration over $`U`$ gives

```math
E[Y^k]=E[U^{pk}]=\int_0^1u^{pk}\,du=\frac1{pk+1}.
```

In particular, $`E[Y]=\frac{1}{3}`$ and $`E[Y^2]=\frac{1}{5}`$. Subtracting the square of the mean gives $`\mathop{\mathrm{Var}} (Y)=\frac{4}{45}`$. The same values follow by integrating $`y f_Y(y)`$ and $`y^2f_Y(y)`$, which checks the Jacobian factor in the density.

The median $`m`$ satisfies $`m^{1/p}=1/2`$, so $`m=2^{-p}=2^{-2}`$. It generally differs from the mean. The transformation preserves the ordering of observations, but changes the distances between them and therefore changes averages and variances.

#### Problem PR-54: power of a uniform variable, case 2

Let $`U\sim\mathop{\mathrm{Uniform}} (0,1)`$ and $`Y=U^p`$ with $`p=3`$. Find the support, CDF, density, mean, variance, and median.

Solution. Because $`p\gt 0`$, the transformation is strictly increasing and maps the unit interval into itself. Thus $`Y`$ has support $`[0,1]`$. For $`0\lt y\lt 1`$,

```math
P(Y\le y)=P(U^p\le y)=P(U\le y^{1/p})=y^{1/p}.
```

The complete CDF is zero for $`y\le0`$, equals $`y^{1/p}`$ inside the interval, and is one for $`y\ge1`$. There are no point masses at the endpoints because $`U`$ has none.

Differentiating on the interior gives

```math
f_Y(y)=\frac1p y^{1/p-1},\qquad 0\lt y\lt 1.
```

Its integral is $`[y^{1/p}]_0^1=1`$. A density can grow without bound near zero and still integrate to one; boundedness of a density is not a requirement.

For any real $`k`$ with $`pk\gt -1`$, direct integration over $`U`$ gives

```math
E[Y^k]=E[U^{pk}]=\int_0^1u^{pk}\,du=\frac1{pk+1}.
```

In particular, $`E[Y]=\frac{1}{4}`$ and $`E[Y^2]=\frac{1}{7}`$. Subtracting the square of the mean gives $`\mathop{\mathrm{Var}} (Y)=\frac{9}{112}`$. The same values follow by integrating $`y f_Y(y)`$ and $`y^2f_Y(y)`$, which checks the Jacobian factor in the density.

The median $`m`$ satisfies $`m^{1/p}=1/2`$, so $`m=2^{-p}=2^{-3}`$. It generally differs from the mean. The transformation preserves the ordering of observations, but changes the distances between them and therefore changes averages and variances.

#### Problem PR-55: power of a uniform variable, case 3

Let $`U\sim\mathop{\mathrm{Uniform}} (0,1)`$ and $`Y=U^p`$ with $`p=\frac{1}{2}`$. Find the support, CDF, density, mean, variance, and median.

Solution. Because $`p\gt 0`$, the transformation is strictly increasing and maps the unit interval into itself. Thus $`Y`$ has support $`[0,1]`$. For $`0\lt y\lt 1`$,

```math
P(Y\le y)=P(U^p\le y)=P(U\le y^{1/p})=y^{1/p}.
```

The complete CDF is zero for $`y\le0`$, equals $`y^{1/p}`$ inside the interval, and is one for $`y\ge1`$. There are no point masses at the endpoints because $`U`$ has none.

Differentiating on the interior gives

```math
f_Y(y)=\frac1p y^{1/p-1},\qquad 0\lt y\lt 1.
```

Its integral is $`[y^{1/p}]_0^1=1`$. A density can grow without bound near zero and still integrate to one; boundedness of a density is not a requirement.

For any real $`k`$ with $`pk\gt -1`$, direct integration over $`U`$ gives

```math
E[Y^k]=E[U^{pk}]=\int_0^1u^{pk}\,du=\frac1{pk+1}.
```

In particular, $`E[Y]=\frac{2}{3}`$ and $`E[Y^2]=\frac{1}{2}`$. Subtracting the square of the mean gives $`\mathop{\mathrm{Var}} (Y)=\frac{1}{18}`$. The same values follow by integrating $`y f_Y(y)`$ and $`y^2f_Y(y)`$, which checks the Jacobian factor in the density.

The median $`m`$ satisfies $`m^{1/p}=1/2`$, so $`m=2^{-p}=2^{-\frac{1}{2}}`$. It generally differs from the mean. The transformation preserves the ordering of observations, but changes the distances between them and therefore changes averages and variances.

#### Problem PR-56: power of a uniform variable, case 4

Let $`U\sim\mathop{\mathrm{Uniform}} (0,1)`$ and $`Y=U^p`$ with $`p=\frac{1}{3}`$. Find the support, CDF, density, mean, variance, and median.

Solution. Because $`p\gt 0`$, the transformation is strictly increasing and maps the unit interval into itself. Thus $`Y`$ has support $`[0,1]`$. For $`0\lt y\lt 1`$,

```math
P(Y\le y)=P(U^p\le y)=P(U\le y^{1/p})=y^{1/p}.
```

The complete CDF is zero for $`y\le0`$, equals $`y^{1/p}`$ inside the interval, and is one for $`y\ge1`$. There are no point masses at the endpoints because $`U`$ has none.

Differentiating on the interior gives

```math
f_Y(y)=\frac1p y^{1/p-1},\qquad 0\lt y\lt 1.
```

Its integral is $`[y^{1/p}]_0^1=1`$. A density can grow without bound near zero and still integrate to one; boundedness of a density is not a requirement.

For any real $`k`$ with $`pk\gt -1`$, direct integration over $`U`$ gives

```math
E[Y^k]=E[U^{pk}]=\int_0^1u^{pk}\,du=\frac1{pk+1}.
```

In particular, $`E[Y]=\frac{3}{4}`$ and $`E[Y^2]=\frac{3}{5}`$. Subtracting the square of the mean gives $`\mathop{\mathrm{Var}} (Y)=\frac{3}{80}`$. The same values follow by integrating $`y f_Y(y)`$ and $`y^2f_Y(y)`$, which checks the Jacobian factor in the density.

The median $`m`$ satisfies $`m^{1/p}=1/2`$, so $`m=2^{-p}=2^{-\frac{1}{3}}`$. It generally differs from the mean. The transformation preserves the ordering of observations, but changes the distances between them and therefore changes averages and variances.

### 9.3 Conditional expectation from a joint table

For binary variables, a four-cell table determines every marginal and conditional quantity. The conditional mean $`E[Y\mid X]`$ is itself a random variable with one value for each possible $`X`$. The law of total variance separates average within-group uncertainty from variation between those conditional means.

#### Problem PR-57: binary table and total variance, case 1

The rows $`X=0,1`$ and columns $`Y=0,1`$ have joint probabilities

```math
\begin{pmatrix}\frac{2}{5}&\frac{1}{10}\\[0pt] \frac{1}{5}&\frac{3}{10}\end{pmatrix}
```

. Find $`E[Y\mid X]`$, covariance, and both parts of the law of total variance.

Solution. Sum each row to obtain $`P(X=0)=\frac{1}{2}`$ and $`P(X=1)=\frac{1}{2}`$. Sum the second column to obtain $`P(Y=1)=\frac{2}{5}`$. All row probabilities are positive, so both conditional distributions are defined.

Since $`Y`$ is binary, its conditional expectation is its conditional probability of one:

```math
E[Y\mid X=0]=\frac{P(X=0,Y=1)}{P(X=0)}=\frac{1}{5},\qquad
E[Y\mid X=1]=\frac{3}{5}.
```

Thus the random variable $`E[Y\mid X]`$ equals the first number when $`X=0`$ and the second when $`X=1`$. Averaging with row probabilities gives $`E[E[Y\mid X]]=\frac{2}{5}=E[Y]`$, which checks the conditioning calculation.

Also $`E[X]=\frac{1}{2}`$ and $`E[XY]=P(X=1,Y=1)=\frac{3}{10}`$. Therefore

```math
\mathop{\mathrm{Cov}} (X,Y)=E[XY]-E[X]E[Y]=\frac{1}{10}.
```

Within each row, $`Y`$ is Bernoulli with success probability $`m_x`$. Its conditional variance is $`m_x(1-m_x)`$. Averaging these two variances gives

```math
E[\mathop{\mathrm{Var}} (Y\mid X)]=\frac{1}{5}.
```

The variance of the conditional mean, computed around its mean $`E[Y]`$, is

```math
\mathop{\mathrm{Var}} (E[Y\mid X])=\frac{1}{25}.
```

Their sum is $`\frac{6}{25}`$, exactly $`P(Y=1)(1-P(Y=1))`$. The equality explains the decomposition: learning $`X`$ removes the between-group part of uncertainty on average, while the within-group part remains. Individual conditional variances need not all be smaller than the unconditional variance.

#### Problem PR-58: binary table and total variance, case 2

The rows $`X=0,1`$ and columns $`Y=0,1`$ have joint probabilities

```math
\begin{pmatrix}\frac{1}{10}&\frac{3}{10}\\[0pt] \frac{1}{5}&\frac{2}{5}\end{pmatrix}
```

. Find $`E[Y\mid X]`$, covariance, and both parts of the law of total variance.

Solution. Sum each row to obtain $`P(X=0)=\frac{2}{5}`$ and $`P(X=1)=\frac{3}{5}`$. Sum the second column to obtain $`P(Y=1)=\frac{7}{10}`$. All row probabilities are positive, so both conditional distributions are defined.

Since $`Y`$ is binary, its conditional expectation is its conditional probability of one:

```math
E[Y\mid X=0]=\frac{P(X=0,Y=1)}{P(X=0)}=\frac{3}{4},\qquad
E[Y\mid X=1]=\frac{2}{3}.
```

Thus the random variable $`E[Y\mid X]`$ equals the first number when $`X=0`$ and the second when $`X=1`$. Averaging with row probabilities gives $`E[E[Y\mid X]]=\frac{7}{10}=E[Y]`$, which checks the conditioning calculation.

Also $`E[X]=\frac{3}{5}`$ and $`E[XY]=P(X=1,Y=1)=\frac{2}{5}`$. Therefore

```math
\mathop{\mathrm{Cov}} (X,Y)=E[XY]-E[X]E[Y]=-\frac{1}{50}.
```

Within each row, $`Y`$ is Bernoulli with success probability $`m_x`$. Its conditional variance is $`m_x(1-m_x)`$. Averaging these two variances gives

```math
E[\mathop{\mathrm{Var}} (Y\mid X)]=\frac{5}{24}.
```

The variance of the conditional mean, computed around its mean $`E[Y]`$, is

```math
\mathop{\mathrm{Var}} (E[Y\mid X])=\frac{1}{600}.
```

Their sum is $`\frac{21}{100}`$, exactly $`P(Y=1)(1-P(Y=1))`$. The equality explains the decomposition: learning $`X`$ removes the between-group part of uncertainty on average, while the within-group part remains. Individual conditional variances need not all be smaller than the unconditional variance.

#### Problem PR-59: binary table and total variance, case 3

The rows $`X=0,1`$ and columns $`Y=0,1`$ have joint probabilities

```math
\begin{pmatrix}\frac{1}{2}&\frac{1}{5}\\[0pt] \frac{1}{10}&\frac{1}{5}\end{pmatrix}
```

. Find $`E[Y\mid X]`$, covariance, and both parts of the law of total variance.

Solution. Sum each row to obtain $`P(X=0)=\frac{7}{10}`$ and $`P(X=1)=\frac{3}{10}`$. Sum the second column to obtain $`P(Y=1)=\frac{2}{5}`$. All row probabilities are positive, so both conditional distributions are defined.

Since $`Y`$ is binary, its conditional expectation is its conditional probability of one:

```math
E[Y\mid X=0]=\frac{P(X=0,Y=1)}{P(X=0)}=\frac{2}{7},\qquad
E[Y\mid X=1]=\frac{2}{3}.
```

Thus the random variable $`E[Y\mid X]`$ equals the first number when $`X=0`$ and the second when $`X=1`$. Averaging with row probabilities gives $`E[E[Y\mid X]]=\frac{2}{5}=E[Y]`$, which checks the conditioning calculation.

Also $`E[X]=\frac{3}{10}`$ and $`E[XY]=P(X=1,Y=1)=\frac{1}{5}`$. Therefore

```math
\mathop{\mathrm{Cov}} (X,Y)=E[XY]-E[X]E[Y]=\frac{2}{25}.
```

Within each row, $`Y`$ is Bernoulli with success probability $`m_x`$. Its conditional variance is $`m_x(1-m_x)`$. Averaging these two variances gives

```math
E[\mathop{\mathrm{Var}} (Y\mid X)]=\frac{22}{105}.
```

The variance of the conditional mean, computed around its mean $`E[Y]`$, is

```math
\mathop{\mathrm{Var}} (E[Y\mid X])=\frac{16}{525}.
```

Their sum is $`\frac{6}{25}`$, exactly $`P(Y=1)(1-P(Y=1))`$. The equality explains the decomposition: learning $`X`$ removes the between-group part of uncertainty on average, while the within-group part remains. Individual conditional variances need not all be smaller than the unconditional variance.

#### Problem PR-60: binary table and total variance, case 4

The rows $`X=0,1`$ and columns $`Y=0,1`$ have joint probabilities

```math
\begin{pmatrix}\frac{1}{5}&\frac{1}{5}\\[0pt] \frac{1}{2}&\frac{1}{10}\end{pmatrix}
```

. Find $`E[Y\mid X]`$, covariance, and both parts of the law of total variance.

Solution. Sum each row to obtain $`P(X=0)=\frac{2}{5}`$ and $`P(X=1)=\frac{3}{5}`$. Sum the second column to obtain $`P(Y=1)=\frac{3}{10}`$. All row probabilities are positive, so both conditional distributions are defined.

Since $`Y`$ is binary, its conditional expectation is its conditional probability of one:

```math
E[Y\mid X=0]=\frac{P(X=0,Y=1)}{P(X=0)}=\frac{1}{2},\qquad
E[Y\mid X=1]=\frac{1}{6}.
```

Thus the random variable $`E[Y\mid X]`$ equals the first number when $`X=0`$ and the second when $`X=1`$. Averaging with row probabilities gives $`E[E[Y\mid X]]=\frac{3}{10}=E[Y]`$, which checks the conditioning calculation.

Also $`E[X]=\frac{3}{5}`$ and $`E[XY]=P(X=1,Y=1)=\frac{1}{10}`$. Therefore

```math
\mathop{\mathrm{Cov}} (X,Y)=E[XY]-E[X]E[Y]=-\frac{2}{25}.
```

Within each row, $`Y`$ is Bernoulli with success probability $`m_x`$. Its conditional variance is $`m_x(1-m_x)`$. Averaging these two variances gives

```math
E[\mathop{\mathrm{Var}} (Y\mid X)]=\frac{11}{60}.
```

The variance of the conditional mean, computed around its mean $`E[Y]`$, is

```math
\mathop{\mathrm{Var}} (E[Y\mid X])=\frac{2}{75}.
```

Their sum is $`\frac{21}{100}`$, exactly $`P(Y=1)(1-P(Y=1))`$. The equality explains the decomposition: learning $`X`$ removes the between-group part of uncertainty on average, while the within-group part remains. Individual conditional variances need not all be smaller than the unconditional variance.

### 9.4 Gaussian conditioning with correlation

For a jointly Gaussian pair, subtracting the linear prediction produces a Gaussian residual uncorrelated with the conditioning variable. Joint Gaussianity then upgrades zero covariance to independence. Without that distributional assumption, the same calculation gives a linear projection rather than necessarily the conditional mean.

#### Problem PR-61: conditional Gaussian prediction, case 1

Let $`(X,Y)`$ be jointly Gaussian with mean $`(1,2)`$ and covariance

```math
\begin{pmatrix}4&3\\[0pt] 3&9\end{pmatrix}
```

. Find the law of $`Y\mid X=3`$ and the conditional probability that $`Y\gt 2`$.

Solution. The covariance matrix is valid and nonsingular: its first diagonal entry is $`4\gt 0`$ and its determinant is $`27\gt 0`$. The correlation is $`0.500000`$, which lies strictly between minus one and one.

Choose $`b=\mathop{\mathrm{Cov}} (X,Y)/\mathop{\mathrm{Var}} (X)=\frac{3}{4}`$ and define

```math
R=Y-\mu_Y-b(X-\mu_X).
```

Then $`E[R]=0`$ and $`\mathop{\mathrm{Cov}} (R,X)=\mathop{\mathrm{Cov}} (Y,X)-b\mathop{\mathrm{Var}} (X)=0`$. Because $`(R,X)`$ is jointly Gaussian, this zero covariance proves independence. Expanding its variance gives

```math
\mathop{\mathrm{Var}} (R)=\sigma_Y^2-\frac{\sigma_{XY}^2}{\sigma_X^2}
=\frac{27}{4}\gt 0.
```

Conditional on $`X=3`$, the linear term becomes a constant while $`R`$ keeps its unconditional Gaussian law. Therefore

```math
Y\mid X=3\sim
N\left(\frac{7}{2},\frac{27}{4}\right).
```

Standardizing the requested threshold gives

```math
P(Y\gt 2\mid X=3)
=1-\Phi\left(\frac{2-(\frac{7}{2})}{\sqrt{\frac{27}{4}}}\right)
\approx 0.718149.
```

The conditional variance does not depend on the observed value of $`X`$ in this Gaussian model. It equals the unconditional variance minus an explained component. A singular covariance matrix would require separate treatment: the residual variance could be zero and the conditional law could be a point mass, making division by its standard deviation invalid.

#### Problem PR-62: conditional Gaussian prediction, case 2

Let $`(X,Y)`$ be jointly Gaussian with mean $`(0,-1)`$ and covariance

```math
\begin{pmatrix}1&-1\\[0pt] -1&4\end{pmatrix}
```

. Find the law of $`Y\mid X=2`$ and the conditional probability that $`Y\gt -1`$.

Solution. The covariance matrix is valid and nonsingular: its first diagonal entry is $`1\gt 0`$ and its determinant is $`3\gt 0`$. The correlation is $`-0.500000`$, which lies strictly between minus one and one.

Choose $`b=\mathop{\mathrm{Cov}} (X,Y)/\mathop{\mathrm{Var}} (X)=-1`$ and define

```math
R=Y-\mu_Y-b(X-\mu_X).
```

Then $`E[R]=0`$ and $`\mathop{\mathrm{Cov}} (R,X)=\mathop{\mathrm{Cov}} (Y,X)-b\mathop{\mathrm{Var}} (X)=0`$. Because $`(R,X)`$ is jointly Gaussian, this zero covariance proves independence. Expanding its variance gives

```math
\mathop{\mathrm{Var}} (R)=\sigma_Y^2-\frac{\sigma_{XY}^2}{\sigma_X^2}
=3\gt 0.
```

Conditional on $`X=2`$, the linear term becomes a constant while $`R`$ keeps its unconditional Gaussian law. Therefore

```math
Y\mid X=2\sim
N\left(-3,3\right).
```

Standardizing the requested threshold gives

```math
P(Y\gt -1\mid X=2)
=1-\Phi\left(\frac{-1-(-3)}{\sqrt{3}}\right)
\approx 0.124107.
```

The conditional variance does not depend on the observed value of $`X`$ in this Gaussian model. It equals the unconditional variance minus an explained component. A singular covariance matrix would require separate treatment: the residual variance could be zero and the conditional law could be a point mass, making division by its standard deviation invalid.

#### Problem PR-63: conditional Gaussian prediction, case 3

Let $`(X,Y)`$ be jointly Gaussian with mean $`(-2,1)`$ and covariance

```math
\begin{pmatrix}9&3\\[0pt] 3&4\end{pmatrix}
```

. Find the law of $`Y\mid X=1`$ and the conditional probability that $`Y\gt 1`$.

Solution. The covariance matrix is valid and nonsingular: its first diagonal entry is $`9\gt 0`$ and its determinant is $`27\gt 0`$. The correlation is $`0.500000`$, which lies strictly between minus one and one.

Choose $`b=\mathop{\mathrm{Cov}} (X,Y)/\mathop{\mathrm{Var}} (X)=\frac{1}{3}`$ and define

```math
R=Y-\mu_Y-b(X-\mu_X).
```

Then $`E[R]=0`$ and $`\mathop{\mathrm{Cov}} (R,X)=\mathop{\mathrm{Cov}} (Y,X)-b\mathop{\mathrm{Var}} (X)=0`$. Because $`(R,X)`$ is jointly Gaussian, this zero covariance proves independence. Expanding its variance gives

```math
\mathop{\mathrm{Var}} (R)=\sigma_Y^2-\frac{\sigma_{XY}^2}{\sigma_X^2}
=3\gt 0.
```

Conditional on $`X=1`$, the linear term becomes a constant while $`R`$ keeps its unconditional Gaussian law. Therefore

```math
Y\mid X=1\sim
N\left(2,3\right).
```

Standardizing the requested threshold gives

```math
P(Y\gt 1\mid X=1)
=1-\Phi\left(\frac{1-(2)}{\sqrt{3}}\right)
\approx 0.718149.
```

The conditional variance does not depend on the observed value of $`X`$ in this Gaussian model. It equals the unconditional variance minus an explained component. A singular covariance matrix would require separate treatment: the residual variance could be zero and the conditional law could be a point mass, making division by its standard deviation invalid.

#### Problem PR-64: conditional Gaussian prediction, case 4

Let $`(X,Y)`$ be jointly Gaussian with mean $`(3,0)`$ and covariance

```math
\begin{pmatrix}4&1\\[0pt] 1&1\end{pmatrix}
```

. Find the law of $`Y\mid X=1`$ and the conditional probability that $`Y\gt 0`$.

Solution. The covariance matrix is valid and nonsingular: its first diagonal entry is $`4\gt 0`$ and its determinant is $`3\gt 0`$. The correlation is $`0.500000`$, which lies strictly between minus one and one.

Choose $`b=\mathop{\mathrm{Cov}} (X,Y)/\mathop{\mathrm{Var}} (X)=\frac{1}{4}`$ and define

```math
R=Y-\mu_Y-b(X-\mu_X).
```

Then $`E[R]=0`$ and $`\mathop{\mathrm{Cov}} (R,X)=\mathop{\mathrm{Cov}} (Y,X)-b\mathop{\mathrm{Var}} (X)=0`$. Because $`(R,X)`$ is jointly Gaussian, this zero covariance proves independence. Expanding its variance gives

```math
\mathop{\mathrm{Var}} (R)=\sigma_Y^2-\frac{\sigma_{XY}^2}{\sigma_X^2}
=\frac{3}{4}\gt 0.
```

Conditional on $`X=1`$, the linear term becomes a constant while $`R`$ keeps its unconditional Gaussian law. Therefore

```math
Y\mid X=1\sim
N\left(-\frac{1}{2},\frac{3}{4}\right).
```

Standardizing the requested threshold gives

```math
P(Y\gt 0\mid X=1)
=1-\Phi\left(\frac{0-(-\frac{1}{2})}{\sqrt{\frac{3}{4}}}\right)
\approx 0.281851.
```

The conditional variance does not depend on the observed value of $`X`$ in this Gaussian model. It equals the unconditional variance minus an explained component. A singular covariance matrix would require separate treatment: the residual variance could be zero and the conditional law could be a point mass, making division by its standard deviation invalid.

### 9.5 Finite-sample bounds and normal approximations

Chebyshev's inequality uses only a variance and gives a guaranteed upper bound. The central limit theorem gives an approximation whose finite-sample accuracy depends on the distribution and threshold. For Bernoulli samples, the binomial distribution allows an exact comparison with both.

#### Problem PR-65: Bernoulli mean error probability, case 1

Let $`X_1,\ldots,X_n`$ be independent Bernoulli variables with $`p=\frac{1}{2}`$ and $`n=100`$. Compare a variance bound, a normal approximation, and the exact probability of $`\vert \bar X-p\vert \ge \frac{1}{10}`$.

Solution. Independence removes cross-covariances, so

```math
E[\bar X]=p,\qquad
\mathop{\mathrm{Var}} (\bar X)=\frac{p(1-p)}n=\frac{1}{400}.
```

Chebyshev's inequality therefore gives

```math
P(\vert \bar X-p\vert \ge\epsilon)
\le \min\left(1,\frac{p(1-p)}{n\epsilon^2}\right)
=\frac{1}{4}.
```

The minimum with one records the trivial probability bound when the variance inequality is uninformative. This is a guaranteed finite-sample statement under the given independence and variance assumptions.

The normal approximation standardizes by $`\sqrt{p(1-p)/n}`$. The threshold is $`z=2.000000`$ standard deviations, giving the uncorrected approximation $`2[1-\Phi(z)]\approx 0.045500`$.

For an exact answer, let $`S=n\bar X`$. Then $`S\sim\mathop{\mathrm{Binomial}} (n,p)`$, and the event is

```math
S\le 40\quad\text{or}\quad S\ge 60.
```

The floor and ceiling follow from the discrete support and the non-strict inequality in the event. Summing the binomial masses in these two tails gives approximately $`0.056888`$.

The exact value is below the Chebyshev bound, as required. It need not equal the normal approximation. A continuity correction can improve the latter by moving the integer cutoffs to half-integer boundaries. The central limit theorem itself is a limit statement for fixed standardized thresholds; it does not guarantee accurate probabilities for every far-tail event in a small sample.

#### Problem PR-66: Bernoulli mean error probability, case 2

Let $`X_1,\ldots,X_n`$ be independent Bernoulli variables with $`p=\frac{1}{10}`$ and $`n=200`$. Compare a variance bound, a normal approximation, and the exact probability of $`\vert \bar X-p\vert \ge \frac{1}{20}`$.

Solution. Independence removes cross-covariances, so

```math
E[\bar X]=p,\qquad
\mathop{\mathrm{Var}} (\bar X)=\frac{p(1-p)}n=\frac{9}{20000}.
```

Chebyshev's inequality therefore gives

```math
P(\vert \bar X-p\vert \ge\epsilon)
\le \min\left(1,\frac{p(1-p)}{n\epsilon^2}\right)
=\frac{9}{50}.
```

The minimum with one records the trivial probability bound when the variance inequality is uninformative. This is a guaranteed finite-sample statement under the given independence and variance assumptions.

The normal approximation standardizes by $`\sqrt{p(1-p)/n}`$. The threshold is $`z=2.357023`$ standard deviations, giving the uncorrected approximation $`2[1-\Phi(z)]\approx 0.018422`$.

For an exact answer, let $`S=n\bar X`$. Then $`S\sim\mathop{\mathrm{Binomial}} (n,p)`$, and the event is

```math
S\le 10\quad\text{or}\quad S\ge 30.
```

The floor and ceiling follow from the discrete support and the non-strict inequality in the event. Summing the binomial masses in these two tails gives approximately $`0.024398`$.

The exact value is below the Chebyshev bound, as required. It need not equal the normal approximation. A continuity correction can improve the latter by moving the integer cutoffs to half-integer boundaries. The central limit theorem itself is a limit statement for fixed standardized thresholds; it does not guarantee accurate probabilities for every far-tail event in a small sample.

#### Problem PR-67: Bernoulli mean error probability, case 3

Let $`X_1,\ldots,X_n`$ be independent Bernoulli variables with $`p=\frac{1}{4}`$ and $`n=80`$. Compare a variance bound, a normal approximation, and the exact probability of $`\vert \bar X-p\vert \ge \frac{1}{8}`$.

Solution. Independence removes cross-covariances, so

```math
E[\bar X]=p,\qquad
\mathop{\mathrm{Var}} (\bar X)=\frac{p(1-p)}n=\frac{3}{1280}.
```

Chebyshev's inequality therefore gives

```math
P(\vert \bar X-p\vert \ge\epsilon)
\le \min\left(1,\frac{p(1-p)}{n\epsilon^2}\right)
=\frac{3}{20}.
```

The minimum with one records the trivial probability bound when the variance inequality is uninformative. This is a guaranteed finite-sample statement under the given independence and variance assumptions.

The normal approximation standardizes by $`\sqrt{p(1-p)/n}`$. The threshold is $`z=2.581989`$ standard deviations, giving the uncorrected approximation $`2[1-\Phi(z)]\approx 0.009823`$.

For an exact answer, let $`S=n\bar X`$. Then $`S\sim\mathop{\mathrm{Binomial}} (n,p)`$, and the event is

```math
S\le 10\quad\text{or}\quad S\ge 30.
```

The floor and ceiling follow from the discrete support and the non-strict inequality in the event. Summing the binomial masses in these two tails gives approximately $`0.013599`$.

The exact value is below the Chebyshev bound, as required. It need not equal the normal approximation. A continuity correction can improve the latter by moving the integer cutoffs to half-integer boundaries. The central limit theorem itself is a limit statement for fixed standardized thresholds; it does not guarantee accurate probabilities for every far-tail event in a small sample.

#### Problem PR-68: Bernoulli mean error probability, case 4

Let $`X_1,\ldots,X_n`$ be independent Bernoulli variables with $`p=\frac{4}{5}`$ and $`n=150`$. Compare a variance bound, a normal approximation, and the exact probability of $`\vert \bar X-p\vert \ge \frac{1}{10}`$.

Solution. Independence removes cross-covariances, so

```math
E[\bar X]=p,\qquad
\mathop{\mathrm{Var}} (\bar X)=\frac{p(1-p)}n=\frac{2}{1875}.
```

Chebyshev's inequality therefore gives

```math
P(\vert \bar X-p\vert \ge\epsilon)
\le \min\left(1,\frac{p(1-p)}{n\epsilon^2}\right)
=\frac{8}{75}.
```

The minimum with one records the trivial probability bound when the variance inequality is uninformative. This is a guaranteed finite-sample statement under the given independence and variance assumptions.

The normal approximation standardizes by $`\sqrt{p(1-p)/n}`$. The threshold is $`z=3.061862`$ standard deviations, giving the uncorrected approximation $`2[1-\Phi(z)]\approx 0.002200`$.

For an exact answer, let $`S=n\bar X`$. Then $`S\sim\mathop{\mathrm{Binomial}} (n,p)`$, and the event is

```math
S\le 105\quad\text{or}\quad S\ge 135.
```

The floor and ceiling follow from the discrete support and the non-strict inequality in the event. Summing the binomial masses in these two tails gives approximately $`0.003070`$.

The exact value is below the Chebyshev bound, as required. It need not equal the normal approximation. A continuity correction can improve the latter by moving the integer cutoffs to half-integer boundaries. The central limit theorem itself is a limit statement for fixed standardized thresholds; it does not guarantee accurate probabilities for every far-tail event in a small sample.

### 9.6 Exponential likelihood and an exact interval

An exponential rate is the reciprocal of a population mean. The maximum-likelihood estimate consequently involves the reciprocal of the sample mean. The sum has a gamma law, which also gives an exact confidence interval and avoids relying entirely on a large-sample approximation.

#### Problem PR-69: rate estimation from positive observations, case 1

The observations $`1, 2, 3, 2, 2`$ are modeled as independent exponential variables with density $`f(x;\lambda)=\lambda e^{-\lambda x}`$ for $`x\gt 0`$ and $`\lambda\gt 0`$. Find the MLE, estimated asymptotic standard error, and an exact 95% interval.

Solution. Write $`S=\sum_i x_i=10`$ and $`n=5`$. Ignoring factors independent of the rate, the likelihood and log-likelihood are

```math
L(\lambda)=\lambda^n e^{-\lambda S},\qquad
\ell(\lambda)=n\log\lambda-\lambda S.
```

The score equation is $`n/\lambda-S=0`$, giving $`\hat\lambda=n/S=\frac{1}{2}`$. Since $`\ell''(\lambda)=-n/\lambda^2\lt 0`$, the log-likelihood is strictly concave. It tends to minus infinity at both ends of the positive parameter line for these data, proving that this is the unique global maximizer.

The Fisher information is $`I_n(\lambda)=n/\lambda^2`$. Substituting the estimate gives asymptotic standard error

```math
\widehat{\mathop{\mathrm{se}} }(\hat\lambda)
=\frac{\hat\lambda}{\sqrt n}\approx 0.223607.
```

This is a large-sample approximation, not an exact small-sample sampling standard deviation.

For an exact interval, the sum of $`n`$ independent exponential variables has a gamma distribution. Equivalently, $`2\lambda S\sim\chi^2_{2n}`$. Thus

```math
P\left(
\frac{\chi^2_{2n,0.025}}{2S}
\le\lambda\le
\frac{\chi^2_{2n,0.975}}{2S}
\right)=0.95.
```

The observed interval is approximately $`[0.162349,1.024159]`$. Its endpoints are positive and it is generally asymmetric around the MLE. The coverage statement refers to repeated samples from this exponential model. Model misspecification, dependence, or selected observations would require a different sampling argument.

#### Problem PR-70: rate estimation from positive observations, case 2

The observations $`\frac{1}{2}, 1, \frac{3}{2}, 2, 1, \frac{1}{2}`$ are modeled as independent exponential variables with density $`f(x;\lambda)=\lambda e^{-\lambda x}`$ for $`x\gt 0`$ and $`\lambda\gt 0`$. Find the MLE, estimated asymptotic standard error, and an exact 95% interval.

Solution. Write $`S=\sum_i x_i=\frac{13}{2}`$ and $`n=6`$. Ignoring factors independent of the rate, the likelihood and log-likelihood are

```math
L(\lambda)=\lambda^n e^{-\lambda S},\qquad
\ell(\lambda)=n\log\lambda-\lambda S.
```

The score equation is $`n/\lambda-S=0`$, giving $`\hat\lambda=n/S=\frac{12}{13}`$. Since $`\ell''(\lambda)=-n/\lambda^2\lt 0`$, the log-likelihood is strictly concave. It tends to minus infinity at both ends of the positive parameter line for these data, proving that this is the unique global maximizer.

The Fisher information is $`I_n(\lambda)=n/\lambda^2`$. Substituting the estimate gives asymptotic standard error

```math
\widehat{\mathop{\mathrm{se}} }(\hat\lambda)
=\frac{\hat\lambda}{\sqrt n}\approx 0.376845.
```

This is a large-sample approximation, not an exact small-sample sampling standard deviation.

For an exact interval, the sum of $`n`$ independent exponential variables has a gamma distribution. Equivalently, $`2\lambda S\sim\chi^2_{2n}`$. Thus

```math
P\left(
\frac{\chi^2_{2n,0.025}}{2S}
\le\lambda\le
\frac{\chi^2_{2n,0.975}}{2S}
\right)=0.95.
```

The observed interval is approximately $`[0.338753,1.795128]`$. Its endpoints are positive and it is generally asymmetric around the MLE. The coverage statement refers to repeated samples from this exponential model. Model misspecification, dependence, or selected observations would require a different sampling argument.

#### Problem PR-71: rate estimation from positive observations, case 3

The observations $`2, 4, 1, 3, 2, 5, 1, 2`$ are modeled as independent exponential variables with density $`f(x;\lambda)=\lambda e^{-\lambda x}`$ for $`x\gt 0`$ and $`\lambda\gt 0`$. Find the MLE, estimated asymptotic standard error, and an exact 95% interval.

Solution. Write $`S=\sum_i x_i=20`$ and $`n=8`$. Ignoring factors independent of the rate, the likelihood and log-likelihood are

```math
L(\lambda)=\lambda^n e^{-\lambda S},\qquad
\ell(\lambda)=n\log\lambda-\lambda S.
```

The score equation is $`n/\lambda-S=0`$, giving $`\hat\lambda=n/S=\frac{2}{5}`$. Since $`\ell''(\lambda)=-n/\lambda^2\lt 0`$, the log-likelihood is strictly concave. It tends to minus infinity at both ends of the positive parameter line for these data, proving that this is the unique global maximizer.

The Fisher information is $`I_n(\lambda)=n/\lambda^2`$. Substituting the estimate gives asymptotic standard error

```math
\widehat{\mathop{\mathrm{se}} }(\hat\lambda)
=\frac{\hat\lambda}{\sqrt n}\approx 0.141421.
```

This is a large-sample approximation, not an exact small-sample sampling standard deviation.

For an exact interval, the sum of $`n`$ independent exponential variables has a gamma distribution. Equivalently, $`2\lambda S\sim\chi^2_{2n}`$. Thus

```math
P\left(
\frac{\chi^2_{2n,0.025}}{2S}
\le\lambda\le
\frac{\chi^2_{2n,0.975}}{2S}
\right)=0.95.
```

The observed interval is approximately $`[0.172692,0.721134]`$. Its endpoints are positive and it is generally asymmetric around the MLE. The coverage statement refers to repeated samples from this exponential model. Model misspecification, dependence, or selected observations would require a different sampling argument.

#### Problem PR-72: rate estimation from positive observations, case 4

The observations $`\frac{1}{4}, \frac{3}{4}, 1, \frac{1}{2}, \frac{3}{2}`$ are modeled as independent exponential variables with density $`f(x;\lambda)=\lambda e^{-\lambda x}`$ for $`x\gt 0`$ and $`\lambda\gt 0`$. Find the MLE, estimated asymptotic standard error, and an exact 95% interval.

Solution. Write $`S=\sum_i x_i=4`$ and $`n=5`$. Ignoring factors independent of the rate, the likelihood and log-likelihood are

```math
L(\lambda)=\lambda^n e^{-\lambda S},\qquad
\ell(\lambda)=n\log\lambda-\lambda S.
```

The score equation is $`n/\lambda-S=0`$, giving $`\hat\lambda=n/S=\frac{5}{4}`$. Since $`\ell''(\lambda)=-n/\lambda^2\lt 0`$, the log-likelihood is strictly concave. It tends to minus infinity at both ends of the positive parameter line for these data, proving that this is the unique global maximizer.

The Fisher information is $`I_n(\lambda)=n/\lambda^2`$. Substituting the estimate gives asymptotic standard error

```math
\widehat{\mathop{\mathrm{se}} }(\hat\lambda)
=\frac{\hat\lambda}{\sqrt n}\approx 0.559017.
```

This is a large-sample approximation, not an exact small-sample sampling standard deviation.

For an exact interval, the sum of $`n`$ independent exponential variables has a gamma distribution. Equivalently, $`2\lambda S\sim\chi^2_{2n}`$. Thus

```math
P\left(
\frac{\chi^2_{2n,0.025}}{2S}
\le\lambda\le
\frac{\chi^2_{2n,0.975}}{2S}
\right)=0.95.
```

The observed interval is approximately $`[0.405872,2.560397]`$. Its endpoints are positive and it is generally asymmetric around the MLE. The coverage statement refers to repeated samples from this exponential model. Model misspecification, dependence, or selected observations would require a different sampling argument.

### 9.7 A Gaussian mean with unknown variance

Replacing an unknown standard deviation by a sample estimate changes the exact pivot from normal to Student's t under Gaussian sampling. The denominator uses $`n-1`$ in the sample variance. The t critical value and the standard error play separate roles.

#### Problem PR-73: Student interval from a small sample, case 1

For independent $`N(\mu,\sigma^2)`$ observations $`2, 4, 5, 7, 7`$ with both parameters unknown, form an exact 95% confidence interval for $`\mu`$.

Solution. The sample contains $`n=5`$ observations with mean $`\bar x=5`$. Centering every observation at this mean gives the sum of squared deviations

```math
\sum_i(x_i-\bar x)^2=18.
```

The sample variance and standard error of the sample mean are

```math
s^2=\frac1{n-1}\sum_i(x_i-\bar x)^2=\frac{9}{2},\qquad
\mathop{\mathrm{se}} (\bar x)=\frac s{\sqrt n}\approx 0.948683.
```

Under independent Gaussian sampling, $`\sqrt n(\bar X-\mu)/\sigma`$ is standard normal, $`(n-1)s^2/\sigma^2`$ is chi-squared with $`n-1`$ degrees of freedom, and the two quantities are independent. Dividing the first by the square root of the second divided by its degrees of freedom gives

```math
\frac{\bar X-\mu}{s/\sqrt n}\sim t_{n-1}.
```

The 97.5th percentile for $`4`$ degrees of freedom is approximately $`2.776445`$. Solving the central probability inequality for $`\mu`$ gives the observed interval

```math
\bar x\pm t_{n-1,0.975}\frac s{\sqrt n}
\approx [2.366033,7.633967].
```

The interval is wider than one using the normal critical value with the same standard error, reflecting uncertainty about the variance. The exact t distribution depends on Gaussian sampling. For independent non-Gaussian observations with finite variance, a similar interval may be justified asymptotically, but this small sample alone does not establish that approximation. Dependence would also change the standard error of the mean.

#### Problem PR-74: Student interval from a small sample, case 2

For independent $`N(\mu,\sigma^2)`$ observations $`10, 11, 9, 12, 13, 11`$ with both parameters unknown, form an exact 95% confidence interval for $`\mu`$.

Solution. The sample contains $`n=6`$ observations with mean $`\bar x=11`$. Centering every observation at this mean gives the sum of squared deviations

```math
\sum_i(x_i-\bar x)^2=10.
```

The sample variance and standard error of the sample mean are

```math
s^2=\frac1{n-1}\sum_i(x_i-\bar x)^2=2,\qquad
\mathop{\mathrm{se}} (\bar x)=\frac s{\sqrt n}\approx 0.577350.
```

Under independent Gaussian sampling, $`\sqrt n(\bar X-\mu)/\sigma`$ is standard normal, $`(n-1)s^2/\sigma^2`$ is chi-squared with $`n-1`$ degrees of freedom, and the two quantities are independent. Dividing the first by the square root of the second divided by its degrees of freedom gives

```math
\frac{\bar X-\mu}{s/\sqrt n}\sim t_{n-1}.
```

The 97.5th percentile for $`5`$ degrees of freedom is approximately $`2.570582`$. Solving the central probability inequality for $`\mu`$ gives the observed interval

```math
\bar x\pm t_{n-1,0.975}\frac s{\sqrt n}
\approx [9.515874,12.484126].
```

The interval is wider than one using the normal critical value with the same standard error, reflecting uncertainty about the variance. The exact t distribution depends on Gaussian sampling. For independent non-Gaussian observations with finite variance, a similar interval may be justified asymptotically, but this small sample alone does not establish that approximation. Dependence would also change the standard error of the mean.

#### Problem PR-75: Student interval from a small sample, case 3

For independent $`N(\mu,\sigma^2)`$ observations $`1, 1, 2, 3, 5, 6, 3, 3`$ with both parameters unknown, form an exact 95% confidence interval for $`\mu`$.

Solution. The sample contains $`n=8`$ observations with mean $`\bar x=3`$. Centering every observation at this mean gives the sum of squared deviations

```math
\sum_i(x_i-\bar x)^2=22.
```

The sample variance and standard error of the sample mean are

```math
s^2=\frac1{n-1}\sum_i(x_i-\bar x)^2=\frac{22}{7},\qquad
\mathop{\mathrm{se}} (\bar x)=\frac s{\sqrt n}\approx 0.626783.
```

Under independent Gaussian sampling, $`\sqrt n(\bar X-\mu)/\sigma`$ is standard normal, $`(n-1)s^2/\sigma^2`$ is chi-squared with $`n-1`$ degrees of freedom, and the two quantities are independent. Dividing the first by the square root of the second divided by its degrees of freedom gives

```math
\frac{\bar X-\mu}{s/\sqrt n}\sim t_{n-1}.
```

The 97.5th percentile for $`7`$ degrees of freedom is approximately $`2.364624`$. Solving the central probability inequality for $`\mu`$ gives the observed interval

```math
\bar x\pm t_{n-1,0.975}\frac s{\sqrt n}
\approx [1.517893,4.482107].
```

The interval is wider than one using the normal critical value with the same standard error, reflecting uncertainty about the variance. The exact t distribution depends on Gaussian sampling. For independent non-Gaussian observations with finite variance, a similar interval may be justified asymptotically, but this small sample alone does not establish that approximation. Dependence would also change the standard error of the mean.

#### Problem PR-76: Student interval from a small sample, case 4

For independent $`N(\mu,\sigma^2)`$ observations $`-2, 0, 1, 2, 4`$ with both parameters unknown, form an exact 95% confidence interval for $`\mu`$.

Solution. The sample contains $`n=5`$ observations with mean $`\bar x=1`$. Centering every observation at this mean gives the sum of squared deviations

```math
\sum_i(x_i-\bar x)^2=20.
```

The sample variance and standard error of the sample mean are

```math
s^2=\frac1{n-1}\sum_i(x_i-\bar x)^2=5,\qquad
\mathop{\mathrm{se}} (\bar x)=\frac s{\sqrt n}\approx 1.000000.
```

Under independent Gaussian sampling, $`\sqrt n(\bar X-\mu)/\sigma`$ is standard normal, $`(n-1)s^2/\sigma^2`$ is chi-squared with $`n-1`$ degrees of freedom, and the two quantities are independent. Dividing the first by the square root of the second divided by its degrees of freedom gives

```math
\frac{\bar X-\mu}{s/\sqrt n}\sim t_{n-1}.
```

The 97.5th percentile for $`4`$ degrees of freedom is approximately $`2.776445`$. Solving the central probability inequality for $`\mu`$ gives the observed interval

```math
\bar x\pm t_{n-1,0.975}\frac s{\sqrt n}
\approx [-1.776445,3.776445].
```

The interval is wider than one using the normal critical value with the same standard error, reflecting uncertainty about the variance. The exact t distribution depends on Gaussian sampling. For independent non-Gaussian observations with finite variance, a similar interval may be justified asymptotically, but this small sample alone does not establish that approximation. Dependence would also change the standard error of the mean.

### 9.8 Dependence and the variance of an average

For a stationary process, a sample mean contains all pairwise covariances, not only individual variances. A positive autocorrelation can make a long sequence contain much less information about the mean than an equally long independent sample. Negative autocorrelation can have the opposite effect.

#### Problem PR-77: AR covariance and effective sample size, case 1

A stationary process has $`\mathop{\mathrm{Cov}} (X_t,X_{t-k})=\gamma_0\rho^{\vert k\vert }`$ with $`\rho=\frac{4}{5}`$ and $`\gamma_0=1`$. For $`n=20`$, find the exact variance of its sample mean and compare it with the independent and long-run formulas.

Solution. Start from the variance of a sum:

```math
\mathop{\mathrm{Var}} (\bar X_n)
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
\mathop{\mathrm{Cov}} (X_i,X_j).
```

There are $`n`$ diagonal terms of size $`\gamma_0`$. For each lag $`k\ge1`$, there are $`n-k`$ pairs above the diagonal and the same number below it. Therefore

```math
\mathop{\mathrm{Var}} (\bar X_n)
=\frac{\gamma_0}n
\left[1+2\sum_{k=1}^{n-1}
\left(1-\frac kn\right)\rho^k\right]
\approx 0.35115292.
```

The bracket equals approximately $`7.023058`$. An independence calculation would give $`\gamma_0/n=0.05000000`$, omitting every off-diagonal covariance.

Since $`\vert \rho\vert \lt 1`$, the covariance series is absolutely summable. Its sum is the long-run variance

```math
\Omega=\gamma_0+2\sum_{k=1}^\infty\gamma_0\rho^k
=\gamma_0\frac{1+\rho}{1-\rho}=9.
```

The large-sample approximation to the variance of the mean is $`\Omega/n=0.45000000`$. Its difference from the exact finite-sample result comes from the weights $`1-k/n`$ and the truncation at lag $`n-1`$.

One descriptive effective sample size solves $`\gamma_0/n_{\mathrm{eff}}=\mathop{\mathrm{Var}} (\bar X_n)`$, giving $`n_{\mathrm{eff}}\approx 2.8478`$. This number concerns estimation of this mean under this covariance model. It is not a universal replacement for sample size in every statistic or test. A dependence-aware bootstrap or long-run variance estimator must preserve or estimate the covariance terms that the independent calculation discards.

#### Problem PR-78: AR covariance and effective sample size, case 2

A stationary process has $`\mathop{\mathrm{Cov}} (X_t,X_{t-k})=\gamma_0\rho^{\vert k\vert }`$ with $`\rho=\frac{1}{2}`$ and $`\gamma_0=4`$. For $`n=50`$, find the exact variance of its sample mean and compare it with the independent and long-run formulas.

Solution. Start from the variance of a sum:

```math
\mathop{\mathrm{Var}} (\bar X_n)
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
\mathop{\mathrm{Cov}} (X_i,X_j).
```

There are $`n`$ diagonal terms of size $`\gamma_0`$. For each lag $`k\ge1`$, there are $`n-k`$ pairs above the diagonal and the same number below it. Therefore

```math
\mathop{\mathrm{Var}} (\bar X_n)
=\frac{\gamma_0}n
\left[1+2\sum_{k=1}^{n-1}
\left(1-\frac kn\right)\rho^k\right]
\approx 0.23360000.
```

The bracket equals approximately $`2.920000`$. An independence calculation would give $`\gamma_0/n=0.08000000`$, omitting every off-diagonal covariance.

Since $`\vert \rho\vert \lt 1`$, the covariance series is absolutely summable. Its sum is the long-run variance

```math
\Omega=\gamma_0+2\sum_{k=1}^\infty\gamma_0\rho^k
=\gamma_0\frac{1+\rho}{1-\rho}=12.
```

The large-sample approximation to the variance of the mean is $`\Omega/n=0.24000000`$. Its difference from the exact finite-sample result comes from the weights $`1-k/n`$ and the truncation at lag $`n-1`$.

One descriptive effective sample size solves $`\gamma_0/n_{\mathrm{eff}}=\mathop{\mathrm{Var}} (\bar X_n)`$, giving $`n_{\mathrm{eff}}\approx 17.1233`$. This number concerns estimation of this mean under this covariance model. It is not a universal replacement for sample size in every statistic or test. A dependence-aware bootstrap or long-run variance estimator must preserve or estimate the covariance terms that the independent calculation discards.

#### Problem PR-79: AR covariance and effective sample size, case 3

A stationary process has $`\mathop{\mathrm{Cov}} (X_t,X_{t-k})=\gamma_0\rho^{\vert k\vert }`$ with $`\rho=-\frac{1}{2}`$ and $`\gamma_0=1`$. For $`n=20`$, find the exact variance of its sample mean and compare it with the independent and long-run formulas.

Solution. Start from the variance of a sum:

```math
\mathop{\mathrm{Var}} (\bar X_n)
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
\mathop{\mathrm{Cov}} (X_i,X_j).
```

There are $`n`$ diagonal terms of size $`\gamma_0`$. For each lag $`k\ge1`$, there are $`n-k`$ pairs above the diagonal and the same number below it. Therefore

```math
\mathop{\mathrm{Var}} (\bar X_n)
=\frac{\gamma_0}n
\left[1+2\sum_{k=1}^{n-1}
\left(1-\frac kn\right)\rho^k\right]
\approx 0.01777778.
```

The bracket equals approximately $`0.355556`$. An independence calculation would give $`\gamma_0/n=0.05000000`$, omitting every off-diagonal covariance.

Since $`\vert \rho\vert \lt 1`$, the covariance series is absolutely summable. Its sum is the long-run variance

```math
\Omega=\gamma_0+2\sum_{k=1}^\infty\gamma_0\rho^k
=\gamma_0\frac{1+\rho}{1-\rho}=\frac{1}{3}.
```

The large-sample approximation to the variance of the mean is $`\Omega/n=0.01666667`$. Its difference from the exact finite-sample result comes from the weights $`1-k/n`$ and the truncation at lag $`n-1`$.

One descriptive effective sample size solves $`\gamma_0/n_{\mathrm{eff}}=\mathop{\mathrm{Var}} (\bar X_n)`$, giving $`n_{\mathrm{eff}}\approx 56.2500`$. This number concerns estimation of this mean under this covariance model. It is not a universal replacement for sample size in every statistic or test. A dependence-aware bootstrap or long-run variance estimator must preserve or estimate the covariance terms that the independent calculation discards.

#### Problem PR-80: AR covariance and effective sample size, case 4

A stationary process has $`\mathop{\mathrm{Cov}} (X_t,X_{t-k})=\gamma_0\rho^{\vert k\vert }`$ with $`\rho=\frac{9}{10}`$ and $`\gamma_0=2`$. For $`n=100`$, find the exact variance of its sample mean and compare it with the independent and long-run formulas.

Solution. Start from the variance of a sum:

```math
\mathop{\mathrm{Var}} (\bar X_n)
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
\mathop{\mathrm{Cov}} (X_i,X_j).
```

There are $`n`$ diagonal terms of size $`\gamma_0`$. For each lag $`k\ge1`$, there are $`n-k`$ pairs above the diagonal and the same number below it. Therefore

```math
\mathop{\mathrm{Var}} (\bar X_n)
=\frac{\gamma_0}n
\left[1+2\sum_{k=1}^{n-1}
\left(1-\frac kn\right)\rho^k\right]
\approx 0.34400096.
```

The bracket equals approximately $`17.200048`$. An independence calculation would give $`\gamma_0/n=0.02000000`$, omitting every off-diagonal covariance.

Since $`\vert \rho\vert \lt 1`$, the covariance series is absolutely summable. Its sum is the long-run variance

```math
\Omega=\gamma_0+2\sum_{k=1}^\infty\gamma_0\rho^k
=\gamma_0\frac{1+\rho}{1-\rho}=38.
```

The large-sample approximation to the variance of the mean is $`\Omega/n=0.38000000`$. Its difference from the exact finite-sample result comes from the weights $`1-k/n`$ and the truncation at lag $`n-1`$.

One descriptive effective sample size solves $`\gamma_0/n_{\mathrm{eff}}=\mathop{\mathrm{Var}} (\bar X_n)`$, giving $`n_{\mathrm{eff}}\approx 5.8139`$. This number concerns estimation of this mean under this covariance model. It is not a universal replacement for sample size in every statistic or test. A dependence-aware bootstrap or long-run variance estimator must preserve or estimate the covariance terms that the independent calculation discards.

<!-- further-distinct-problems-PR -->

### 9.9 Further distribution and inference problems

These additional problems are constructed for these notes.

#### Problem PR-81: a Poisson model with a random intensity

Let $`\Lambda`$ have gamma density $`b^k\lambda^{k-1}e^{-b\lambda}/\Gamma(k)`$ for $`\lambda\gt 0`$, with $`k,b\gt 0`$. Conditional on $`\Lambda`$, let $`N`$ be Poisson with mean $`\Lambda`$. Find the marginal mass function, mean, and variance of $`N`$.

Solution. Multiply the conditional Poisson mass by the mixing density and integrate:

```math
P(N=n)=\frac{b^k}{n!\Gamma(k)}
\int_0^\infty\lambda^{n+k-1}e^{-(b+1)\lambda}d\lambda
=\frac{\Gamma(n+k)}{n!\Gamma(k)}
\frac{b^k}{(b+1)^{n+k}}.
```

The gamma-integral identity follows from substituting $`u=(b+1)\lambda`$ into its defining integral. The expression is a negative-binomial mass function, with this parameterization valid for noninteger positive $`k`$ as well.

Conditional mean and variance both equal $`\Lambda`$. The tower property and total variance give

```math
E[N]=E[\Lambda]=\frac kb,\qquad
\mathop{\mathrm{Var}} (N)=E[\Lambda]+\mathop{\mathrm{Var}} (\Lambda)
=\frac kb+\frac k{b^2}.
```

The variance exceeds the mean because intensity variation adds an extra component. For $`k=2,b=1`$, the mean is two and variance four, while the probability of zero is $`1/4`$.

This is a generative mixture calculation. It does not require interpreting the mixing distribution as a subjective prior. The extra dispersion comes from heterogeneity in the conditional intensity, even though each conditional count is exactly Poisson.

#### Problem PR-82: multinomial covariance and a log ratio

For $`n`$ independent categorical observations with positive category probabilities $`p_1,\ldots,p_d`$, derive the asymptotic covariance of sample proportions and the delta-method variance of $`\log(\hat p_1/\hat p_2)`$.

Solution. Represent one observation by the indicator vector $`I`$ with exactly one entry equal to one. Then $`E[I]=p`$, each diagonal second moment is $`p_i`$, and off-diagonal products vanish. Therefore

```math
\Sigma=\mathop{\mathrm{Cov}} (I)=\mathop{\mathrm{diag}} (p)-pp^\top.
```

The multivariate CLT gives $`\sqrt n(\hat p-p)\Rightarrow N(0,\Sigma)`$. The covariance is singular because the proportions always sum to one. For any vector $`a`$, its quadratic form is the variance of a discrete variable taking values $`a_i`$ with positive probabilities $`p_i`$. It vanishes exactly when all $`a_i`$ are equal, so the rank is $`d-1`$.

For $`g(p)=\log p_1-\log p_2`$, the gradient has entries $`1/p_1,-1/p_2`$ and zeros elsewhere. Substituting into the covariance form gives

```math
\nabla g^\top\Sigma\nabla g=\frac1{p_1}+\frac1{p_2}.
```

The negative covariance between distinct proportions contributes to this result and must not be dropped.

For counts $`(20,30,50)`$, the estimated log ratio is $`\log(2/3)`$ and its asymptotic standard error is $`\sqrt{(5+10/3)/100}=1/\sqrt{12}`$. A log-scale interval is $`\log(2/3)\pm1.96/\sqrt{12}`$; exponentiating gives a positive interval for the probability ratio.

Zero observed counts require separate treatment because the log and plug-in gradient are undefined. The approximation is a fixed-interior-probability limit, rather than a uniform guarantee for categories becoming arbitrarily rare.

#### Problem PR-83: a variance-stabilizing transformation

For a Bernoulli sample mean with fixed $`0\lt p\lt 1`$, find a differentiable increasing transformation whose square-root-sample-size asymptotic variance is one. Construct the resulting approximate interval and state its boundary limitation.

Solution. The sample proportion has asymptotic variance $`p(1-p)`$. The delta method gives transformed variance $`g'(p)^2p(1-p)`$. Requiring it to equal one and choosing the increasing sign gives $`g'(p)=1/\sqrt{p(1-p)}`$.

Integration gives $`g(p)=2\arcsin\sqrt p`$, up to an irrelevant additive constant. Thus

```math
\sqrt n\,[g(\hat p)-g(p)]\Rightarrow N(0,1).
```

An approximate interval for $`g(p)`$ is $`g(\hat p)\pm z_{1-\alpha/2}/\sqrt n`$. Its parameter range is $`[0,\pi]`$, so intersect with that range before applying the increasing inverse $`p=\sin^2(g/2)`$.

For $`n=100,\hat p=0.1`$, the center is $`2\arcsin\sqrt{0.1}\approx0.643501`$. At 95%, the untruncated transformed endpoints are approximately $`0.447501`$ and $`0.839501`$. Applying the inverse gives a bounded probability interval rather than directly adding a symmetric probability-scale error.

The transformation derivative diverges at zero and one. Its fixed-interior delta-method derivation does not establish exact coverage at those boundaries or for very rare successes in a small sample. A bounded interval is not automatically an accurately calibrated interval. The variance-stabilizing property is an asymptotic calculation under the stated sampling model.

#### Problem PR-84: a normal likelihood ratio and Student's t

For an iid normal sample with unknown variance and nonzero residual sum of squares, derive the likelihood-ratio statistic for $`H_0:\mu=\mu_0`$. Relate its exact rejection rule to the two-sided t test.

Solution. Under the unrestricted model, the MLEs are $`\hat\mu=\bar X`$ and $`\hat v=n^{-1}\sum_i(X_i-\bar X)^2`$. Under the null, the variance MLE is $`v_0=n^{-1}\sum_i(X_i-\mu_0)^2`$.

The sum-of-squares decomposition gives $`v_0=\hat v+(\bar X-\mu_0)^2`$. At either fitted variance, the maximized log likelihood differs only through $`-n\log v/2`$. Therefore

```math
-2\log\Lambda=n\log\left(\frac{v_0}{\hat v}\right)
=n\log\left(1+\frac{t^2}{n-1}\right),
\qquad
t=\frac{\sqrt n(\bar X-\mu_0)}s.
```

The statistic is strictly increasing in $`\vert t\vert`$. Since $`t\sim t_{n-1}`$ exactly under the Gaussian null, an exact size-$`\alpha`$ LR test rejects when

```math
-2\log\Lambda\gt
n\log\left(1+\frac{t_{n-1,1-\alpha/2}^2}{n-1}\right).
```

This is exactly the usual two-sided t rejection region. Using the asymptotic chi-square critical value for the LR statistic is a different finite-sample calibration, though the two become close in large samples.

For $`n=5,\bar x=5,s^2=4,\mu_0=4`$, $`t=\sqrt5/2`$ and the LR statistic is $`5\log(1+5/16)`$. The t magnitude is below the 5% critical value $`2.776`$, so the exact test does not reject.

The relation was derived for positive fitted variance. An all-equal sample would require handling the boundary likelihood rather than inserting zero into these ratios.

#### Problem PR-85: an exact interval for a uniform endpoint

For independent uniform observations on $`(0,\theta)`$, let $`M`$ be the maximum. Construct a finite-sample interval with coverage $`1-\alpha`$ using only $`M`$ and derive its endpoint factor.

Solution. The pivot $`U=M/\theta`$ has CDF $`u^n`$ on $`[0,1]`$. Also $`M\le\theta`$ almost surely. Seek an interval $`[M,cM]`$ with $`c\gt 1`$.

Its coverage is

```math
P_\theta(M\le\theta\le cM)
=P_\theta(M/\theta\ge1/c)
=1-c^{-n}.
```

Setting this equal to $`1-\alpha`$ gives $`c=\alpha^{-1/n}`$. Thus

```math
\left[M,\frac{M}{\alpha^{1/n}}\right]
```

has exact coverage for every positive endpoint under the iid uniform model.

For $`n=10`$ and $`\alpha=0.05`$, the factor is approximately $`1.34928`$. The lower endpoint is the observed maximum because every smaller parameter value assigns zero likelihood to the sample. The upper endpoint reflects the chance that all observations happen to remain below a fraction of the true endpoint.

The confidence statement holds for each fixed $`\theta`$ over repeated samples. It does not assign probability $`1-\alpha`$ to a fixed parameter after observing the interval. It also uses the endpoint model's exact pivot, not the regular Fisher-information or normal approximation that fails for this moving-support family.

#### Problem PR-86: familywise error without independent tests

For $`m`$ tests with valid null p-values, prove Bonferroni familywise error control and the Holm step-down rule. Apply both at level $`0.05`$ to ordered p-values $`0.001,0.0054,0.02,0.08,0.1,0.2,0.3,0.4,0.6,0.8`$.

Solution. Validity means $`P(p_i\le t)\le t`$ under its true null. Bonferroni rejects when $`p_i\le\alpha/m`$. If $`m_0`$ nulls are true, the union bound gives probability of any false rejection at most $`m_0\alpha/m\le\alpha`$. No independence was used.

Holm orders the p-values and compares the $`j`$th to $`\alpha/(m-j+1)`$, rejecting consecutively until the first failure. To prove control, consider the first true null in the ordered list, at position $`j`$. If any true null is rejected, that first one must pass its threshold.

There are at least $`m_0`$ remaining hypotheses at that position, so $`\alpha/(m-j+1)\le\alpha/m_0`$. A false rejection therefore implies that the minimum of the true-null p-values is at most $`\alpha/m_0`$. The union bound over the $`m_0`$ true nulls gives probability at most $`\alpha`$.

For the supplied ten tests, Bonferroni uses threshold $`0.005`$ and rejects only the first. Holm rejects the first at $`0.005`$, then the second at $`0.05/9\approx0.005556`$. The third fails threshold $`0.05/8=0.00625`$, so the procedure stops and rejects exactly two.

Familywise control concerns the probability of at least one false rejection. It is not the same target as controlling the expected proportion of false discoveries. The proof depends on valid individual p-values; adjustment cannot repair invalid sampling assumptions in the underlying tests.

[Contents](README.md) · [Part 1](03-probability-and-statistics-part-1.md)
