# Econometrics and Time Series — Part 1

The same regression formula can describe a projection, estimate a statistical model, or identify a causal effect. The algebra may be identical, but the assumptions supporting those interpretations are different. These notes state the interpretation before using an estimator.

Sources:

- Joshua Angrist, [MIT 14.32, Econometrics, Spring 2007](https://ocw.mit.edu/courses/14-32-econometrics-spring-2007/).
- Victor Chernozhukov and Iván Fernández-Val, [MIT 14.382, Econometrics, Spring 2017](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/), lecture notes, their numbered problems, and the archived exams included with the course.
- Anna Mikusheva, [MIT 14.384, Time Series Analysis, Fall 2013](https://ocw.mit.edu/courses/14-384-time-series-analysis-fall-2013/), lecture notes and the explicitly labeled bank of sample problems.

Source problems below select theoretical questions or identified subparts. Numerical examples use small constructed data, so they can be checked without software or an external dataset. The time-series problem bank is a sample collection, not a complete set of the original assignments.

[Contents](README.md) · [Part 2](04-econometrics-and-time-series-part-2.md)

Contents:

- [1. Linear regression and projection](#1-linear-regression-and-projection)
- [2. Identification and misspecification](#2-identification-and-misspecification)
- [3. Variance estimation and dependence](#3-variance-estimation-and-dependence)
- [4. Instrumental variables](#4-instrumental-variables)
- [5. Generalized method of moments](#5-generalized-method-of-moments)
- [6. Panel data and treatment comparisons](#6-panel-data-and-treatment-comparisons)
- [7. Stationary and nonstationary time series](#7-stationary-and-nonstationary-time-series)

## 1. Linear regression and projection

### 1.1 Theory

Let $`Y`$ be a scalar outcome and $`X`$ a $`k`$-dimensional regressor vector, often including a constant. Assume finite second moments and $`Q=E[XX^\top]`$ invertible. The population linear projection coefficient solves

```math
\beta=\arg\min_b E[(Y-X^\top b)^2].
```

The objective is a strictly convex quadratic, so $`\beta=Q^{-1}E[XY]`$. Its residual $`u=Y-X^\top\beta`$ satisfies $`E[Xu]=0`$ by construction. This orthogonality does not imply the stronger condition $`E[u\mid X]=0`$. A conditional mean can be nonlinear while its best linear approximation remains well-defined.

In a sample, let $`X`$ also denote the $`n\times k`$ design matrix, with row $`X_i^\top`$, and let $`y`$ be the outcome column. If its columns are independent,

```math
\widehat\beta=(X^\top X)^{-1}X^\top y,\quad
\widehat y=P_Xy,\quad
\widehat u=M_Xy,
```

where $`P_X=X(X^\top X)^{-1}X^\top`$ and $`M_X=I-P_X`$. These are sample algebraic identities. They do not require random sampling, exogeneity, or Gaussian errors.

If an intercept is included, sample residuals sum to zero because the constant column is one of the columns to which they are orthogonal. Without an intercept, this conclusion need not hold. A zero sample residual-regressor inner product is guaranteed by least squares even when the regressor is endogenous; it cannot validate exogeneity.

### 1.2 Derivations and proofs

#### OLS and its residual orthogonality

Expand $`S(b)=(y-Xb)^\top(y-Xb)`$. Its gradient is $`-2X^\top(y-Xb)`$, so stationarity gives $`X^\top X\widehat\beta=X^\top y`$. The Hessian $`2X^\top X`$ is PD under full column rank, proving a unique minimum. Substituting back gives $`X^\top\widehat u=0`$.

Under a structural model $`y=X\beta_0+u`$, subtract $`\beta_0`$:

```math
\widehat\beta-\beta_0=(X^\top X)^{-1}X^\top u.
```

This identity exposes exactly where statistical assumptions enter. Conditional mean zero gives unbiasedness. A law of large numbers for $`n^{-1}X^\top u`$ gives consistency if $`n^{-1}X^\top X`$ approaches an invertible limit. A CLT for the score sum gives an asymptotic distribution.

#### Frisch-Waugh-Lovell

Partition the design into a regressor of interest $`d`$ and controls $`W`$, and minimize $`\Vert y-d\alpha-W\gamma\Vert ^2`$. For fixed $`\alpha`$, the best controls coefficient is $`(W^\top W)^{-1}W^\top(y-d\alpha)`$. Substitution leaves

```math
\min_\alpha\Vert M_Wy-M_Wd\,\alpha\Vert ^2.
```

Hence, provided $`d^\top M_Wd\gt 0`$,

```math
\widehat\alpha=(d^\top M_Wd)^{-1}d^\top M_Wy.
```

Residualize both the outcome and the target regressor against the same controls, then regress one residual on the other. Residualizing only the outcome while leaving the denominator $`d^\top d`$ unchanged generally gives a different answer. The condition $`d^\top M_Wd\gt 0`$ says that the target has variation left after the controls are removed.

#### A population projection is defined by its approximation class

Let $`X`$ contain a constant and have nonsingular second-moment matrix $`Q=E[XX^\top]`$. The population coefficient $`\beta=Q^{-1}E[XY]`$ minimizes expected squared error among linear functions of $`X`$. Its residual $`u=Y-X^\top\beta`$ satisfies $`E[Xu]=0`$.

This does not imply $`E[u\mid X]=0`$. A nonlinear conditional mean can leave a residual with zero correlation against every included regressor and nonzero conditional expectation. Enlarging the regressor list with squares or interactions changes the projection space and can change its coefficients.

In sample, $`\hat u=(I-P_X)y`$ is orthogonal to the columns of $`X`$ whenever OLS is well defined. That identity holds even under endogeneity, misspecification, and dependent errors. It describes how the coefficients were fitted; it is not a diagnostic proving the underlying model assumptions.

Partialling out controls has a geometric explanation. If $`W`$ is the control matrix, every residualized variable $`M_Wz`$ is perpendicular to the control space. Regressing $`M_Wy`$ on $`M_WX`$ therefore uses only variation not representable by those controls. The inverse in this residualized regression requires $`X^\top M_WX`$ nonsingular. A regressor perfectly explained by the controls cannot have a separate coefficient identified from the sample.

#### Leverage and deleting an observation

The diagonal $`h_i=x_i^\top(X^\top X)^{-1}x_i`$ of the hat matrix measures how strongly observation $`i`$ contributes to its own fitted value. For a full-rank design, $`0\le h_i\le1`$ and $`\sum_i h_i=p`$, where $`p`$ counts fitted coefficients including the intercept.

When deleting observation $`i`$ leaves full rank, the deleted-observation prediction error is

```math
y_i-x_i^\top\hat\beta_{(-i)}
=\frac{\hat u_i}{1-h_i}.
```

One derivation applies the rank-one inverse update to $`X^\top X-x_ix_i^\top`$ and substitutes $`X^\top y-x_iy_i`$. The denominator shows why a small training residual at high leverage can coexist with a much larger deletion error.

#### Example EC-E19: leverage at endpoints

Use an intercept and $`x=(0,1,2,3)`$. The mean regressor is $`1.5`$ and its centered sum of squares is five. Thus $`h_i=1/4+(x_i-1.5)^2/5`$, giving $`(0.7,0.3,0.3,0.7)`$. They sum to two, the number of coefficients.

An endpoint residual $`0.3`$ corresponds to deletion error $`0.3/(1-0.7)=1`$. An interior residual of the same size gives $`0.3/0.7`$. The design position changes how much an observation can pull the fitted line toward itself. Leverage alone does not determine influence; the residual size also matters.

### 1.3 Worked examples

#### Example EC-E01: a slope by centered sums

Use observations $`(x,y)=(0,1),(1,2),(2,2)`$. Their means are $`1`$ and $`5/3`$. The centered cross-product sum is $`(-1)(-2/3)+0(1/3)+1(1/3)=1`$, and the centered $`x`$ square sum is two. Thus the slope is $`1/2`$ and intercept $`7/6`$. This is FWL with the constant as the only control.

#### Example EC-E02: removing a group control

Four observations have $`d=(0,1,0,1)`$, group indicator $`g=(0,0,1,1)`$, and outcomes $`(1,3,4,8)`$. Include a constant and $`g`$.

Solution. Subtract group means: the target residual is $`(-1/2,1/2,-1/2,1/2)`$ and outcome residual is $`(-1,1,-2,2)`$. Their cross-product is three and the target residual sum of squares is one. The common adjusted slope is three, averaging the within-group differences two and four because the two groups have equal sizes and treatment variation.

### 1.4 Selected exercises

#### Problem EC-01: a binary regressor reproduces two conditional means

Source: [MIT 14.32, Problem Set II, B.1](https://ocw.mit.edu/courses/14-32-econometrics-spring-2007/9865af44aed8e5f2486ed3b2f533fac6_psy07n02.pdf).

Problem. Show that the population linear regression of $`Y`$ on an intercept and a nondegenerate binary variable reproduces the two conditional means.

Solution. Let $`P(D=1)=p\in(0,1)`$ and $`m_j=E[Y\mid D=j]`$. Then $`EY=(1-p)m_0+pm_1`$, $`E[DY]=pm_1`$, and $`\mathop{\mathrm{Var}} (D)=p(1-p)`$. Consequently

```math
\mathop{\mathrm{Cov}} (D,Y)=pm_1-p[(1-p)m_0+pm_1]
=p(1-p)(m_1-m_0).
```

The projection slope is $`m_1-m_0`$ and the intercept is $`EY-p(m_1-m_0)=m_0`$. The fitted function equals the conditional mean at both possible regressor values. This follows because a function on a two-point domain has exactly two free values.

#### Problem EC-02: linearity of population residualization

Source: [MIT 14.382, Lecture 1, Problems, 5](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/c62d33e015c910b0d126bcc9344cf2c5_MIT14_382S17_lec1.pdf).

Problem. For a fixed finite-second-moment control vector with nonsingular second-moment matrix, prove that population residualization is a linear operator.

Solution. For fixed controls $`W`$ with invertible $`E[WW^\top]`$, define $`R(V)=V-W^\top E[WW^\top]^{-1}E[WV]`$. Linearity of expectation gives

```math
R(aV+bU)=aV+bU-W^\top Q_W^{-1}
(aE[WV]+bE[WU])=aR(V)+bR(U).
```

All terms exist for finite-second-moment variables by Cauchy-Schwarz. The operator is fixed by the control space; changing controls between the two residualizations would invalidate the identity.

### 1.5 Additional solved problems

#### Problem EC-03: orthogonality without a correct linear conditional mean

Constructed problem. Let $`X`$ be uniform on $`(-1,1)`$ and $`Y=X^2`$. Find the projection on a constant and $`X`$.

Solution. Symmetry gives $`EX=EX^3=0`$, so the slope is zero and the intercept is $`EY=1/3`$. The residual $`u=X^2-1/3`$ has mean zero and $`E[Xu]=0`$. Yet $`E[u\mid X]=X^2-1/3`$, generally nonzero. The linear projection conditions hold even though the conditional mean is quadratic.

#### Problem EC-04: perfect collinearity

Constructed problem. Suppose a regression contains a constant and both group indicators $`D`$ and $`1-D`$.

Solution. The columns satisfy $`\mathbf1=D+(1-D)`$, so the design is rank deficient. If coefficients are $`(a,b,c)`$, adding $`(t,-t,-t)`$ leaves every fitted value unchanged. Dropping one indicator gives a unique parameterization of the same fitted-value space, provided both groups occur. The nonuniqueness comes from representation, rather than an inability to estimate the two group means.

#### Problem EC-05: sample orthogonality cannot test exogeneity

Constructed problem. A researcher finds $`X^\top\widehat u=0`$ and concludes $`E[Xu]=0`$.

Solution. The sample equality is enforced by the OLS normal equations for every outcome vector. If $`Y=\beta_0X+u`$ with $`E[Xu]\ne0`$, OLS can absorb the correlated part of $`u`$ into its slope, leaving a fitted residual orthogonal to $`X`$. Its population projection target is then $`\beta_0+E[Xu]/E[X^2]`$, rather than the structural $`\beta_0`$. A mechanical property of the fitting rule supplies no independent evidence for the identifying moment.

## 2. Identification and misspecification

### 2.1 Theory

A parameter is identified if the observable distribution and the stated restrictions determine it uniquely. More observations reduce sampling uncertainty around an identified quantity. They do not choose between observationally equivalent structural explanations.

In a structural regression, strict conditional exogeneity $`E[u_i\mid X]=0`$ is a finite-sample condition on the full design. The weaker moment condition $`E[X_iu_i]=0`$ can suffice for consistency under suitable sampling laws. A regressor is endogenous relative to the specified structural error if that identifying orthogonality fails.

Omitted variables, simultaneous determination, measurement error, and misspecified functional forms can all produce endogeneity. The direction of bias must be derived from their covariances; it should not be guessed from a label.

### 2.2 Derivations and proofs

#### Omitted-variable bias

Suppose $`Y=\alpha+\beta X+\gamma Z+u`$ with $`E[u\mid X,Z]=0`$, and regress $`Y`$ on a constant and $`X`$ only. Its population slope is

```math
\widetilde\beta
=\frac{\mathop{\mathrm{Cov}} (X,Y)}{\mathop{\mathrm{Var}} (X)}
=\beta+\gamma\frac{\mathop{\mathrm{Cov}} (X,Z)}{\mathop{\mathrm{Var}} (X)}.
```

The omitted term shifts the slope only when it both affects $`Y`$ and is correlated with $`X`$. If $`\gamma\gt 0`$ and $`\mathop{\mathrm{Cov}} (X,Z)\gt 0`$, the bias is positive. Reversing either sign reverses the bias. With controls $`W`$, apply the same formula to variables residualized against $`W`$.

#### Classical measurement error

Let the true model be $`Y=\alpha+\beta X^*+u`$, but observe $`X=X^*+v`$. Assume $`v`$ is mean zero and uncorrelated with $`X^*,u`$, and $`u`$ is uncorrelated with $`X^*`$. Then

```math
\mathop{\mathrm{plim}} \widehat\beta
=\beta\frac{\mathop{\mathrm{Var}} (X^*)}
{\mathop{\mathrm{Var}} (X^*)+\mathop{\mathrm{Var}} (v)}.
```

The numerator covariance stays $`\beta\mathop{\mathrm{Var}} (X^*)`$, while noise increases the denominator variance. This attenuation result depends on the classical error assumptions; nonclassical measurement error need not attenuate toward zero.

#### Simultaneous determination creates regressor-error covariance

Suppose $`Y=\beta D+u`$ and $`D=\gamma Y+v`$, where the two disturbances are centered and uncorrelated. If $`1-\gamma\beta\ne0`$, solve the system:

```math
D=\frac{\gamma u+v}{1-\gamma\beta}.
```

Even though $`u,v`$ are uncorrelated, $`D`$ contains $`u`$. Its covariance with the outcome error is $`\gamma\sigma_u^2/(1-\gamma\beta)`$. The population OLS slope in the outcome equation is

```math
\beta+
\frac{\gamma\sigma_u^2(1-\gamma\beta)}
{\gamma^2\sigma_u^2+\sigma_v^2}.
```

With $`\beta=1`$, $`\gamma=1/2`$, and both disturbance variances one, the OLS slope is $`1.2`$. Increasing the sample size estimates this projection more precisely; it does not remove its difference from the structural coefficient.

An instrument can enter the equation determining $`D`$ while remaining excluded from the outcome error. Both restrictions require an economic or experimental argument. The reduced-form algebra shows what follows from those assumptions but cannot establish them.

#### Two different measurement-error models

In classical measurement error, observe $`X=X^*+\eta`$ while $`Y=\beta X^*+u`$. If the centered measurement error is uncorrelated with $`X^*,u`$ and $`E[X^*u]=0`$, then

```math
\frac{\mathop{\mathrm{Cov}} (X,Y)}{\mathop{\mathrm{Var}} (X)}
=\beta\frac{\mathop{\mathrm{Var}} (X^*)}
{\mathop{\mathrm{Var}} (X^*)+\mathop{\mathrm{Var}} (\eta)}.
```

For a single regressor this shrinks the slope toward zero. The assumptions are doing the work: correlated or systematic measurement error can produce other directions of bias, and multiple-regressor effects need not all be simple attenuation.

In a Berkson model, instead write $`X^*=X+\eta`$ with $`\eta`$ independent of the assigned or observed $`X`$. Substitution gives $`Y=\beta X+(\beta\eta+u)`$. Under the corresponding orthogonality assumptions, the slope on $`X`$ is not attenuated, though the composite error variance increases. The two models reverse which variable is treated as the underlying predictor.

#### Example EC-E20: equal noise ratios, different conclusions

Let $`\beta=2`$ and both the true-regressor variance and classical error variance be one. The classical observed-regressor slope is one. In the Berkson formulation with observed-regressor variance one and independent assignment error variance one, the population slope remains two.

Both descriptions contain “measurement noise of variance one.” Their conditional structure differs. Stating the measurement equation is necessary before using an attenuation formula.

### 2.3 Worked examples

#### Example EC-E03: a signed omitted-variable calculation

Let $`\beta=2`$, $`\gamma=3`$, $`\mathop{\mathrm{Var}} (X)=4`$, and $`\mathop{\mathrm{Cov}} (X,Z)=-2`$. The short regression slope is $`2+3(-2)/4=1/2`$. The downward shift is $`3/2`$. A positive effect of the omitted variable is compatible with a negative bias when the omitted variable covaries negatively with the included one.

#### Example EC-E04: simultaneous determination

Suppose centered variables satisfy $`Y=\beta X+u`$ and $`X=\delta Y+v`$, with $`1-\delta\beta\ne0`$ and independent centered $`u,v`$.

Solution. Substitution gives $`X=(\delta u+v)/(1-\delta\beta)`$, so $`\mathop{\mathrm{Cov}} (X,u)=\delta\mathop{\mathrm{Var}} (u)/(1-\delta\beta)`$. Unless $`\delta=0`$ or the structural shock is degenerate, OLS's exogeneity moment fails. The error enters the regressor through the feedback equation.

### 2.4 Selected exercises

For Problems EC-06 and EC-07, use $`Y=X\beta+\varepsilon`$ with $`\varepsilon_i=\sigma(e_i-1)`$, iid standard exponential $`e_i`$, and errors independent of the full-rank design.

#### Problem EC-06: Gauss-Markov without normal errors

Source: [MIT 14.382, archived Midterm 2006, Question 1](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/cdc4509ac2df5df8f6f0a018607f92f6_MIT_14_382S17_midterm06.pdf).

Problem. In the linear model with independent errors $`u_i=\sigma(e_i-1)`$, $`e_i`$ exponential with rate one and independent of the regressors, check the Gauss-Markov conditions.

Solution. Since $`Ee_i=1`$ and $`\mathop{\mathrm{Var}} (e_i)=1`$, the errors have conditional mean zero, conditional variance $`\sigma^2`$, and no cross-observation covariance. With a linear model and full-column-rank $`X`$, these satisfy the usual Gauss-Markov assumptions. The errors are skewed, but Gaussianity is not required for that theorem. Normal-error finite-sample $`t`$ and $`F`$ reference laws are a separate issue.

#### Problem EC-07: conditional mean, covariance, and distribution

Source: [MIT 14.382, archived Midterm 2006, Question 2](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/cdc4509ac2df5df8f6f0a018607f92f6_MIT_14_382S17_midterm06.pdf).

Problem. For the shifted-exponential error model in EC-06, derive the conditional mean and covariance of OLS and determine whether its finite-sample distribution is Gaussian.

Solution. Apply the OLS error identity:

```math
E[\widehat\beta\mid X]=\beta,\qquad
\mathop{\mathrm{Var}} (\widehat\beta\mid X)=\sigma^2(X^\top X)^{-1}.
```

The estimator is a linear combination of centered exponential errors, so it is generally not exactly Gaussian. For an intercept-only design, $`\widehat\beta-\beta=\sigma(\bar e-1)`$ has support $`[-\sigma,\infty)`$ and is visibly nonnormal. A large-sample normal approximation can still hold after centering and scaling under suitable moment and design conditions.

### 2.5 Additional solved problems

#### Problem EC-08: an invalid functional form

Constructed problem. If $`Y=X^2+u`$, $`E[u\mid X]=0`$, and $`X`$ is uniform on $`(0,1)`$, find the population regression on $`1,X`$.

Solution. The moments are $`EX=1/2`$, $`EX^2=1/3`$, $`EX^3=1/4`$, and $`\mathop{\mathrm{Var}} (X)=1/12`$. Thus the slope is $`(1/4-(1/2)(1/3))/(1/12)=1`$ and intercept $`1/3-1/2=-1/6`$. The residual conditional mean is $`X^2-X+1/6`$, which is not zero. The slope describes a linear projection of a curved mean rather than a constant marginal effect.

#### Problem EC-09: an observational equivalence

Constructed problem. Explain why any proposed scalar slope $`b`$ can be fitted to the identity $`Y=bX+u_b`$ without an error restriction.

Solution. Define $`u_b=Y-bX`$. This reproduces the observed $`Y`$ for every $`b`$. The observable joint law of $`(X,Y)`$ alone does not rule out any of these decompositions. A condition such as $`E[Xu_b]=0`$ selects the projection slope if $`E[X^2]\gt 0`$; a causal interpretation requires an error restriction supported by the causal model. Merely naming $`u_b`$ a structural shock adds no identifying information.

#### Problem EC-10: measuring reliability

Constructed problem. A classical-error regressor has signal variance nine and noise variance three. Compute the probability limit for true slope four.

Solution. The reliability ratio is $`9/(9+3)=3/4`$, so the short regression slope tends to three. If independent repeated measurements $`X_1=X^*+v_1`$ and $`X_2=X^*+v_2`$ are available, their covariance identifies the signal variance under the additional independence assumptions. Their correlation is not automatically a correction factor unless their noise variances and scaling are handled correctly.

## 3. Variance estimation and dependence

### 3.1 Theory

Under $`E[u\mid X]=0`$, the exact conditional covariance is

```math
\mathop{\mathrm{Var}} (\widehat\beta\mid X)
=(X^\top X)^{-1}X^\top\Omega_X X(X^\top X)^{-1},
```

where $`\Omega_X=\mathop{\mathrm{Var}} (u\mid X)`$. If $`\Omega_X=\sigma^2I`$, this reduces to $`\sigma^2(X^\top X)^{-1}`$. Heteroskedasticity changes the diagonal variances; serial correlation introduces off-diagonal covariances.

For iid observations with finite appropriate moments, define $`Q=E[X_iX_i^\top]`$ and $`S=E[X_iX_i^\top u_i^2]`$. Then $`\sqrt n(\widehat\beta-\beta)\Rightarrow N(0,Q^{-1}SQ^{-1})`$. The covariance of the unscaled estimator is approximately this expression divided by $`n`$.

An HC0 sample covariance is

```math
\widehat V_{\mathrm{HC0}}
=(X^\top X)^{-1}
\left(\sum_iX_iX_i^\top\widehat u_i^2\right)
(X^\top X)^{-1}.
```

This allows unequal conditional variances but assumes cross-observation dependence is absent or otherwise negligible. It is not automatically appropriate for a time series.

For time-series scores $`g_t=X_tu_t`$, the asymptotic middle matrix is the long-run covariance $`\sum_{k=-\infty}^\infty E[g_tg_{t-k}^\top]`$. A Bartlett/Newey-West estimator uses weighted sample lag covariances with weights $`1-k/(L+1)`$. Consistency requires conditions on dependence, moments, bandwidth growth, and residual estimation.

Clustering groups observations whose errors may be correlated within the same cluster, while treating different clusters as independent or sufficiently weakly dependent. The effective asymptotic number of independent units is the number of clusters. Many observations in a few clusters do not justify a many-cluster normal approximation.

### 3.2 Derivations and proofs

#### The sandwich from the estimation error

Write

```math
\sqrt n(\widehat\beta-\beta)
=\left(\frac1n\sum_iX_iX_i^\top\right)^{-1}
\frac1{\sqrt n}\sum_iX_iu_i.
```

A law of large numbers sends the first factor to $`Q^{-1}`$. An iid vector CLT sends the second to $`N(0,S)`$. Slutsky produces the sandwich covariance. With time dependence, the first limit can remain the same while the score CLT has a different covariance, namely its long-run covariance.

#### Gauss-Markov proof

Consider any estimator linear in $`y`$, $`Cy`$, that is conditionally unbiased for every $`\beta`$. It must satisfy $`CX=I`$. Write $`C=(X^\top X)^{-1}X^\top+D`$, so $`DX=0`$. Under conditional covariance $`\sigma^2I`$,

```math
\mathop{\mathrm{Var}} (Cy\mid X)
=\sigma^2(X^\top X)^{-1}+\sigma^2DD^\top.
```

The cross terms vanish because $`X^\top D^\top=0`$. Since $`DD^\top`$ is PSD, OLS has the smallest covariance in the linear unbiased class. This says nothing about all nonlinear unbiased estimators under a specified nonnormal distribution.

#### Covariance of the score determines OLS uncertainty

Let $`z_t=x_tu_t`$ and $`Q=E[x_tx_t^\top]`$. Under consistency, a suitable law of large numbers, and a central limit theorem,

```math
\sqrt n(\hat\beta-\beta)
=Q^{-1}\frac1{\sqrt n}\sum_t z_t+o_p(1).
```

For serially dependent scores, their limiting covariance is the long-run sum $`\Omega=\sum_{k=-\infty}^{\infty}\mathop{\mathrm{Cov}} (z_t,z_{t-k})`$. The limiting coefficient covariance is $`Q^{-1}\Omega Q^{-1}`$. The sampling covariance of $`\hat\beta`$ is approximately this matrix divided by $`n`$.

For independent observations, only lag zero remains. A heteroskedasticity-robust estimator replaces it by squared-residual score products. HC1 rescales HC0 by $`n/(n-p)`$; HC2 divides squared residuals by $`1-h_i`$; HC3 divides them by $`(1-h_i)^2`$. These are finite-sample leverage corrections with different properties, rather than guarantees of exact inference under arbitrary heteroskedasticity.

For time dependence, a Bartlett-kernel HAC estimate uses

```math
\hat\Omega=\hat\Gamma_0+
\sum_{k=1}^L\left(1-\frac{k}{L+1}\right)
(\hat\Gamma_k+\hat\Gamma_k^\top),
\qquad
\hat\Gamma_k=\frac1n\sum_{t=k+1}^n\hat z_t\hat z_{t-k}^\top.
```

Bandwidth, moment conditions, and dependence assumptions determine consistency. Letting the bandwidth grow while remaining small relative to the sample is part of common results, but does not cover every nonstationary or long-memory process.

#### Clusters change the number of independent pieces

For independent clusters $`g`$, allow arbitrary error covariance within each cluster. The cluster sandwich sums the outer products of cluster scores $`X_g^\top\hat u_g`$. Its large-sample justification depends on many sufficiently independent clusters and restrictions preventing one cluster from dominating. A large number of rows inside a few clusters does not provide the same asymptotic information.

#### Example EC-E21: shared noise within a group

Let $`Y_{gi}=\mu+\eta_g+\epsilon_{gi}`$ for $`G`$ independent groups of equal size $`m`$. Shared effects have variance $`\tau^2`$, individual errors variance $`\sigma^2`$, and all distinct components are independent.

The grand mean has variance

```math
\frac{\tau^2}{G}+\frac{\sigma^2}{Gm}.
```

Treating all rows as independent would give $`(\tau^2+\sigma^2)/(Gm)`$. For $`G=20`$, $`m=10`$, and both variances one, the correct variance is $`0.055`$ while the independent-row calculation gives $`0.01`$. The standard-error ratio is $`\sqrt{5.5}\approx2.345`$. Adding more people within a group reduces the individual component but cannot eliminate the shared group component.

### 3.3 Worked examples

#### Example EC-E05: a heteroskedastic average

Estimate a common mean from two independent observations with variances one and nine.

Solution. The ordinary average is unbiased with variance $`(1+9)/4=5/2`$. An inverse-variance weighted mean gives weights $`9/10,1/10`$ and variance $`(9/10)^2(1)+(1/10)^2(9)=9/10`$. Homoskedastic OLS efficiency does not extend unchanged to this known unequal-variance setting.

#### Example EC-E06: autocorrelated errors

In an intercept-only model with stationary AR(1) errors and coefficient $`\rho=1/2`$, the long-run variance is $`\gamma_0(1+\rho)/(1-\rho)=3\gamma_0`$. The iid standard error targets $`\sqrt{\gamma_0/n}`$, while the long-sample dependent-data standard error is approximately $`\sqrt{3\gamma_0/n}`$. The coefficient estimate can remain consistent while its naive uncertainty estimate fails.

### 3.4 Selected exercises

#### Problem EC-11: the meaning of BLUE

Source: [MIT 14.382, archived Midterm 2006, Question 3](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/cdc4509ac2df5df8f6f0a018607f92f6_MIT_14_382S17_midterm06.pdf).

Problem. Determine whether OLS is BLUE in the nonnormal error model of EC-06, stating the class of estimators compared.

Solution. BLUE means best linear unbiased estimator, with best interpreted by the PSD ordering of covariance matrices. Under the independent shifted-exponential errors in Problems EC-06 and EC-07, conditional homoskedasticity and exogeneity hold. The Gauss-Markov proof above applies and makes OLS BLUE despite the nonnormal sampling law. The class restriction to estimators linear in the outcome is part of the theorem.

#### Problem EC-12: asymptotic normality with primitive conditions

Source: [MIT 14.382, archived Midterm 2006, Question 6](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/cdc4509ac2df5df8f6f0a018607f92f6_MIT_14_382S17_midterm06.pdf).

Problem. Give primitive iid moment and rank conditions yielding asymptotic normality of OLS in the shifted-exponential model, and derive its limiting covariance.

Solution. Assume iid rows, fixed regressor dimension, $`E\Vert X_i\Vert ^2\lt \infty`$, PD $`Q=E[X_iX_i^\top]`$, and independence of the exponential errors from the regressors. Then $`E[X_iu_i]=0`$ and $`E[\Vert X_iu_i\Vert ^2]=\sigma^2E\Vert X_i\Vert ^2\lt \infty`$. The LLN for $`X_iX_i^\top`$ and vector CLT for $`X_iu_i`$ yield

```math
\sqrt n(\widehat\beta-\beta)\Rightarrow N(0,\sigma^2Q^{-1}).
```

The exponentially skewed finite-sample errors have finite variance, so they satisfy the CLT assumptions. Gaussian errors are unnecessary for this limit.

#### Problem EC-13: estimating the limiting covariance

Source: [MIT 14.382, archived Midterm 2006, Question 7](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/cdc4509ac2df5df8f6f0a018607f92f6_MIT_14_382S17_midterm06.pdf).

Problem. Construct consistent estimators of the disturbance variance, the limiting covariance of $`\sqrt n(\hat\beta-\beta)`$, and the sampling covariance of $`\hat\beta`$ under the conditions of EC-12.

Solution. Put $`\widehat Q=X^\top X/n`$ and $`s^2=\widehat u^\top\widehat u/(n-k)`$. Under the conditions of EC-12, $`\widehat Q\to_pQ`$. For $`d=\widehat\beta-\beta=o_p(1)`$,

```math
\frac1n\sum_i\widehat u_i^2
=\frac1n\sum_i u_i^2
-2d^\top\frac1n\sum_iX_iu_i+d^\top\widehat Qd\to_p\sigma^2.
```

The first term converges by the LLN; the next two vanish. Since $`n/(n-k)\to1`$, $`s^2\to_p\sigma^2`$. Therefore $`s^2\widehat Q^{-1}`$ consistently estimates the covariance of the $`\sqrt n`$ limit. Divide by $`n`$ to estimate the covariance of $`\widehat\beta`$ itself.

#### Problem EC-14: a one-sided asymptotic test

Source: [MIT 14.382, archived Midterm 2006, Question 8](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/cdc4509ac2df5df8f6f0a018607f92f6_MIT_14_382S17_midterm06.pdf).

Problem. Using the asymptotic normal approximation in the shifted-exponential model, test $`H_0:\beta_j=0`$ against $`\beta_j\lt 0`$ at asymptotic level 5%, with a justified standard error.

Solution. For $`H_0:\beta_j=0`$ against $`\beta_j\lt 0`$, use $`t_j=\widehat\beta_j/\sqrt{s^2[(X^\top X)^{-1}]_{jj}}`$. The numerator CLT and consistent denominator give $`t_j\Rightarrow N(0,1)`$ under the null. Reject when $`t_j\lt -1.64485`$ for asymptotic 5% level. An exact $`t_{n-k}`$ law has not been established because the errors are not normal.

### 3.5 Additional solved problems

#### Problem EC-15: an HC sandwich by hand

Constructed problem. In an intercept-only regression, simplify HC0 and compare it with the conventional sample-mean variance estimate.

Solution. Here $`X^\top X=n`$ and $`\widehat u_i=y_i-\bar y`$. Thus HC0 is $`\sum_i\widehat u_i^2/n^2`$. The usual estimate $`s^2/n`$, using denominator $`n-1`$, is $`\sum_i\widehat u_i^2/[n(n-1)]`$. Multiplying HC0 by $`n/(n-1)`$ makes them equal. Their first-order limits agree; the finite-sample degrees-of-freedom convention differs.

#### Problem EC-16: cluster covariance

Constructed problem. Write a covariance estimator for independent clusters $`g=1,\ldots,G`$ with arbitrary within-cluster dependence.

Solution. Let $`s_g=\sum_{i\in g}X_i\widehat u_i`$. An unadjusted cluster sandwich is

```math
(X^\top X)^{-1}\left(\sum_{g=1}^Gs_gs_g^\top\right)(X^\top X)^{-1}.
```

Expanding each outer product includes all within-cluster cross-observation products. HC0 keeps only the diagonal individual terms and would omit these. Consistency requires an asymptotic framework with enough independent clusters and control of cluster sizes and moments. Small-cluster corrections and reference distributions require additional choices.

#### Problem EC-17: a Bartlett HAC estimate

Constructed problem. For a scalar score series, sample lag covariances are $`\widehat\gamma_0=4`$, $`\widehat\gamma_1=1`$, $`\widehat\gamma_2=1/2`$. Compute Bartlett HAC with bandwidth two.

Solution. The weights are $`2/3`$ and $`1/3`$, so

```math
\widehat\Omega=4+2[(2/3)(1)+(1/3)(1/2)]=\frac{17}{3}.
```

The lag-zero-only estimate is four. The positive sample lag covariances increase the estimate of long-run uncertainty. The weights taper high lags, whose sample estimates use fewer pairs and can be noisier. The bandwidth must be specified when reporting the estimate.

## 4. Instrumental variables

### 4.1 Theory

An instrument $`Z`$ supplies variation that is related to an endogenous regressor $`D`$ and satisfies a structural error restriction. In a centered scalar model $`Y=\beta D+u`$, the identifying conditions are

```math
E[Zu]=0,\qquad E[ZD]\ne0.
```

The first is validity and the second relevance. An exclusion restriction gives a model-based reason why $`Z`$ affects the outcome only through the included endogenous channel, after controls. Statistical relevance alone does not establish exclusion or exogeneity.

With a constant, use covariance versions of these expressions. With controls, first residualize $`Y,D,Z`$ against them. With several instruments, identification requires a full-column-rank population cross-moment between instruments and endogenous regressors.

Two-stage least squares uses the part of the endogenous regressors projected onto the instrument space:

```math
\widehat\beta_{\mathrm{2SLS}}=(X^\top P_ZX)^{-1}X^\top P_Zy.
```

Here $`X`$ includes all regressors, and $`Z`$ includes instruments plus exogenous controls. The original structural outcome is used in the second stage. Standard errors from treating the first-stage fitted values as ordinary fixed regressors generally miss the structural error calculation.

### 4.2 Derivations and proofs

#### Just-identified IV

The restriction $`E[Z(Y-\beta D)]=0`$ gives $`\beta=E[ZY]/E[ZD]`$ when the denominator is nonzero. Replace expectations by averages:

```math
\widehat\beta_{\mathrm{IV}}=\frac{\sum_iZ_iY_i}{\sum_iZ_iD_i}.
```

Its error is $`\sum_iZ_iu_i/\sum_iZ_iD_i`$. Under iid finite moments and strong relevance, the denominator divided by $`n`$ approaches a nonzero constant and the numerator divided by $`n`$ approaches zero, proving consistency. A CLT yields asymptotic variance $`E[Z^2u^2]/E[ZD]^2`$ for the $`\sqrt n`$ error.

If relevance becomes weak, the denominator is no longer reliably separated from zero. The ratio can become unstable and skewed; standard normal approximations based on a fixed nonzero denominator need not describe it.

#### Partialling out controls in IV

For $`Y=D\alpha+W\beta+u`$ and instruments $`(Z,W)`$, the control moment is $`W^\top(Y-D\alpha-W\beta)=0`$. Solve it for $`\beta`$ and substitute into the instrument moment:

```math
Z^\top M_W(Y-D\alpha)=0.
```

Thus scalar just-identified IV is $`\widehat\alpha=(Z^\top M_WD)^{-1}Z^\top M_WY`$. With multiple excluded instruments, apply 2SLS to the residualized variables. This is a matrix identity; structural validity of the residualized moment is still an assumption.

#### Rank conditions with several instruments

For $`Y=X^\top\beta+u`$ and instrument vector $`Z\in\mathbb R^q`$, the population restrictions are $`E[Z(Y-X^\top\beta)]=0`$. Identification of $`p`$ coefficients requires the $`q\times p`$ matrix $`E[ZX^\top]`$ to have column rank $`p`$. Having at least $`p`$ instruments is necessary but not sufficient: several instruments can shift exactly the same regressor combination.

For a full-rank sample instrument matrix, $`P_Z=Z(Z^\top Z)^{-1}Z^\top`$. Two-stage least squares gives

```math
\hat\beta=(X^\top P_ZX)^{-1}X^\top P_Zy.
```

The second-stage residual is orthogonal to $`P_ZX`$. In an overidentified model it need not be orthogonal to every column of $`Z`$. There may be more instrument moment equations than coefficients, so generally they cannot all be solved exactly.

Fitted first-stage regressors are generated quantities. Computing a naive ordinary second-stage covariance as though they were the original regressors and the fitted residual were the structural error does not yield the usual correct 2SLS inference formula.

#### Weak relevance is a limiting-distribution problem

For a scalar first stage $`D=\pi Z+v`$ with normalized instrument variance, fixed nonzero $`\pi`$ supports the usual ratio expansion. If $`\pi_n=c/\sqrt n`$, the denominator's deterministic signal and random sampling error have the same order. The IV estimator can then behave like a ratio of correlated normals, rather than a concentrated normal around the target.

The first-stage concentration measure is proportional to $`n\pi^2/\sigma_v^2`$ under this normalization. Replacing $`Z`$ by a large multiple changes its numerical coefficient but not its information content. Relevance is about the regressor variation explained relative to noise.

#### Example EC-E22: invalid instruments that agree with one another

Let $`D=\Pi^\top Z+v`$ with $`E[Zv]=0`$, and suppose $`Y=\beta D+c\Pi^\top Z+u`$ with $`E[Zu]=0`$. Then

```math
E[Z(Y-bD)]
=(\beta+c-b)E[ZD].
```

All instrument moments vanish at $`b=\beta+c`$, although the intended outcome error at $`\beta`$ contains the excluded direct effect. Multiple instruments can therefore agree on the same shifted coefficient.

An overidentification test can detect certain incompatibilities among moments, but passing it does not prove exclusion. It cannot distinguish this observationally aligned violation from a different structural coefficient with valid moments. Identification assumptions require support beyond a failure to reject a specification test.

### 4.3 Worked examples

#### Example EC-E07: a scalar IV calculation

Suppose sample centered cross-products are $`\sum Z_iD_i=20`$ and $`\sum Z_iY_i=50`$. The IV slope is $`2.5`$. If $`D=\pi Z+v`$ and $`Y=2D+u`$, the difference from two comes from the sample score $`\sum Z_iu_i=10`$. With a much smaller first-stage cross-product, the same score error produces a larger slope error.

#### Example EC-E08: a relevant but invalid instrument

Let $`D=Z+v`$ and $`Y=2D+Z+e`$, with independent centered $`Z,v,e`$ and $`\mathop{\mathrm{Var}} (Z)\gt 0`$.

Solution. The instrument strongly predicts $`D`$, but the structural error relative to slope two is $`u=Z+e`$. Then $`E[Zu]=\mathop{\mathrm{Var}} (Z)\gt 0`$. IV converges to $`E[ZY]/E[ZD]=3`$, not two. The direct $`Z`$ effect violates exclusion.

### 4.4 Selected exercises

#### Problem EC-18: OLS, reduced forms, and IV identification

Source: [MIT 14.382, Lecture 2, Appendix A, Problem 1](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/dbf8d6d02238c7ad02d665d7202b0aff_MIT14_382S17_lec2.pdf).

Problem. Analyze the scalar IV model using population orthogonality and relevance: obtain the OLS probability limit, reduced forms, and the identifying IV moment equation.

Solution. For centered $`Y=\alpha D+u`$ with $`E[Du]\ne0`$, OLS approaches $`\alpha+E[Du]/E[D^2]`$ and is inconsistent for $`\alpha`$. Suppose $`D=\pi Z+v`$ with $`E[Zv]=0`$, $`E[Zu]=0`$, and $`\pi\ne0`$. The reduced form is $`Y=\alpha\pi Z+(\alpha v+u)`$. Its slope on $`Z`$ is $`\alpha\pi`$, while the first-stage slope is $`\pi`$. Taking their ratio identifies $`\alpha`$. Replacing those two slopes by consistent sample estimates gives indirect least squares, algebraically the scalar IV estimator. With controls, use residualized variables and corresponding conditional or partialled moments.

#### Problem EC-19: partialling out in a linear IV system

Source: [MIT 14.382, archived Final 2016, Question 1](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/dc4003ef9928b631c627562e6c11edef_MIT_14_382S17_ExamF16.pdf).

Problem. In a linear IV system with included controls, derive the instrumented coefficient after partialling the controls out of outcome, endogenous regressors, and instruments.

Solution. In the stated outcome and first-stage equations, apply $`M_W`$ to remove their $`W`$ terms: $`\widetilde Y=\widetilde D\alpha+\widetilde u`$ and $`\widetilde D=\widetilde Z\delta+\widetilde v`$. The IV moments are then expressed only in residualized variables. Provided the residualized instrument-regressor cross-moment has full rank, estimation of $`\alpha`$ can proceed in this smaller system. Afterward recover the controls coefficient by regressing $`Y-D\widehat\alpha`$ on $`W`$. The algebra removes nuisance linear coefficients; it does not make invalid excluded instruments valid.

### 4.5 Additional solved problems

#### Problem EC-20: weak relevance and the denominator

Constructed problem. In $`D_i=\pi Z_i+v_i`$, compare fixed $`\pi\ne0`$ with $`\pi=c/\sqrt n`$.

Solution. With $`EZ_i^2=q\gt 0`$, the denominator $`\sum Z_iD_i`$ has deterministic leading term $`n\pi q`$ and random variation typically of order $`\sqrt n`$. At fixed $`\pi`$, the leading term dominates. With $`\pi=c/\sqrt n`$, both have order $`\sqrt n`$. The IV estimation error then involves a ratio whose denominator remains random at first order. Dividing by an estimated conventional standard error does not restore the fixed-relevance normal proof.

#### Problem EC-21: two instruments and one regressor

Constructed problem. Explain why two relevant instruments do not automatically identify two different structural effects.

Solution. If the structural model contains one endogenous regressor, two instruments supply two moment restrictions for the same coefficient. They can improve precision or provide an overidentifying restriction, conditional on validity. To identify two endogenous coefficients, the matrix $`E[ZX^\top]`$ must have rank two. If both instruments predict only the same single combination of regressors, the rank remains one and the separate effects are not identified.

#### Problem EC-22: an Anderson-Rubin idea

Constructed problem. For a hypothesized $`\beta_0`$, explain a test based on $`Z_i(Y_i-\beta_0D_i)`$.

Solution. Under the null and valid instruments, this moment has mean zero regardless of how strongly $`Z`$ predicts $`D`$. A CLT for its sample average, with an appropriate covariance estimate, yields a test without dividing by the first-stage coefficient. Invert the tests across $`\beta_0`$ to form a confidence set. Weak relevance can make the set wide, unbounded, or disconnected; that reflects limited identifying information rather than a computational failure.

## 5. Generalized method of moments

### 5.1 Theory

GMM begins with restrictions $`E[g_i(\theta_0)]=0`$, where $`g_i`$ is an $`m`$-dimensional vector and the parameter has dimension $`k\le m`$. These restrictions come from a model: exogeneity, an Euler equation, a distributional assumption, or another justified relationship. Writing a moment does not make it valid.

Let $`\bar g_n(\theta)=n^{-1}\sum_i g_i(\theta)`$. GMM chooses a parameter minimizing

```math
Q_n(\theta)=\bar g_n(\theta)^\top W_n\bar g_n(\theta),
```

where $`W_n`$ is positive definite or appropriately positive semidefinite. In a just-identified problem, a solution setting every sample moment to zero is unaffected by the weighting matrix. With more moments than parameters, they generally cannot all be zero, so the weights determine how their discrepancies are balanced.

Identification concerns the population criterion. If $`E[g_i(\theta)]=0`$ at several parameter values, minimizing sample moments need not select a unique structural target. Local identification is supported by a full-column-rank Jacobian $`G=E[\partial g_i(\theta_0)/\partial\theta^\top]`$. It is not a global uniqueness proof.

For independent observations, the moment covariance is $`S=E[g_i(\theta_0)g_i(\theta_0)^\top]`$. For time series it is the long-run covariance. Under regularity and strong identification,

```math
\sqrt n(\widehat\theta-\theta_0)\Rightarrow
N(0,(G^\top WG)^{-1}G^\top WSWG(G^\top WG)^{-1}).
```

Efficient GMM uses $`W=S^{-1}`$ and has covariance $`(G^\top S^{-1}G)^{-1}`$. This efficiency claim compares estimators using the same valid moments and a regular local asymptotic framework.

### 5.2 Derivations and proofs

#### Linear IV as GMM

For $`g_i(\beta)=Z_i(Y_i-X_i^\top\beta)`$, define $`a=Z^\top y/n`$ and $`B=Z^\top X/n`$. Then $`\bar g_n(\beta)=a-B\beta`$. Differentiation gives

```math
B^\top W_nB\,\widehat\beta=B^\top W_na.
```

If the coefficient matrix is invertible, $`\widehat\beta=(B^\top W_nB)^{-1}B^\top W_na`$. Choosing $`W_n=(Z^\top Z/n)^{-1}`$ gives 2SLS after cancellation of the $`n`$ factors. The asymptotic optimal weighting instead uses the estimated covariance of the instrument-error products, which need not be proportional to $`E[ZZ^\top]`$ under heteroskedasticity.

#### Linearization and the role of weights

The first-order condition is approximately $`G^\top W\bar g_n(\widehat\theta)=0`$. Expand $`\bar g_n(\widehat\theta)=\bar g_n(\theta_0)+G(\widehat\theta-\theta_0)+o_p(n^{-1/2})`$:

```math
\sqrt n(\widehat\theta-\theta_0)
=-(G^\top WG)^{-1}G^\top W\sqrt n\bar g_n(\theta_0)+o_p(1).
```

Applying the moment CLT gives the stated covariance. Smoothness, consistency, a uniform remainder bound, and full rank justify this expansion. In a nonlinear problem, solving a first-order condition alone does not establish the global minimizer or consistency.

#### The GMM linearization and efficient weight

Let $`g_n(\theta)`$ be a $`q`$-vector of sample moments, with $`g(\theta_0)=0`$ and derivative $`D=\partial g(\theta_0)/\partial\theta^\top`$ of full column rank $`p`$. Suppose $`\sqrt n g_n(\theta_0)\Rightarrow N(0,\Omega)`$. For a positive-definite limiting weight $`W`$, the first-order condition of the quadratic moment criterion gives

```math
\sqrt n(\hat\theta-\theta_0)
=-(D^\top WD)^{-1}D^\top W\sqrt n g_n(\theta_0)+o_p(1).
```

Its covariance is

```math
(D^\top WD)^{-1}
D^\top W\Omega WD
(D^\top WD)^{-1}.
```

Taking $`W=\Omega^{-1}`$ simplifies this to $`(D^\top\Omega^{-1}D)^{-1}`$. The inverse covariance downweights noisy combinations and accounts for correlation among moments. It does not simply assign inverse marginal variances when the off-diagonal entries matter.

A two-step estimator first obtains a consistent coefficient with a convenient weight, estimates $`\Omega`$ using that coefficient's moments, and minimizes again with the estimated inverse. Consistency of the first step, the covariance estimate, and the relevant moment laws are assumptions behind this procedure.

#### Why the J statistic has fewer degrees of freedom

At the true parameter, the whitened moment vector is asymptotically standard normal in $`q`$ dimensions. Estimating $`p`$ parameters removes its projection onto the $`p`$-dimensional column space of $`\Omega^{-1/2}D`$. The remaining vector lies in the orthogonal complement, of dimension $`q-p`$.

The squared length of that residual Gaussian vector gives

```math
J=n\,g_n(\hat\theta)^\top\hat\Omega^{-1}g_n(\hat\theta)
\Rightarrow\chi^2_{q-p}.
```

This argument requires correct moments, rank, an appropriate CLT, and consistent weighting. With exactly as many moments as coefficients, a locally solved moment system has no remaining overidentifying directions and the statistic supplies no specification test.

#### Example EC-E23: two mean restrictions for one parameter

Suppose two sample means $`a,b`$ estimate the same $`\theta`$ and $`\sqrt n((a,b)-\theta\mathbf1)`$ has identity covariance. Efficient GMM minimizes $`(a-\theta)^2+(b-\theta)^2`$, so $`\hat\theta=(a+b)/2`$.

The two residual moments are $`((a-b)/2,(b-a)/2)`$. Thus $`J=n(a-b)^2/2`$ and has a one-degree-of-freedom asymptotic reference under the common-mean hypothesis.

For $`a=1`$, $`b=1.2`$, and $`n=100`$, $`\hat\theta=1.1`$ and $`J=2`$. The approximate p-value is $`0.1573`$. The conclusion concerns compatibility of these two moment estimates under the specified covariance; a small difference need not imply that both are unbiased for the intended structural parameter.

### 5.3 Worked examples

#### Example EC-E09: estimating one mean from two noisy measurements

Two unbiased moment measurements of $`\theta`$ have independent asymptotic errors with variances one and four. Efficient weighting is proportional to the inverse variances, giving normalized weights $`4/5`$ and $`1/5`$. The combined asymptotic variance is $`(4/5)^2+(1/5)^2(4)=4/5`$. Equal weighting would give $`5/4`$.

#### Example EC-E10: redundant moments

Suppose $`g_2(\theta)=2g_1(\theta)`$ for every observation and parameter. There are two written moment equations but only one independent restriction. Their covariance matrix is singular. Counting the two rows as two independent overidentifying restrictions would be wrong; remove the redundancy or use a rank-aware formulation.

### 5.4 Selected exercises

#### Problem EC-23: the overidentification statistic

Source: [MIT 14.382, Lecture 4, Problems, 3](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/e1aecf6a41da164603168343fcce76bd_MIT14_382S17_lec4.pdf).

Problem. Derive the limiting distribution of the optimally weighted overidentification statistic under correctly specified moments and a full-rank derivative.

Solution. Assume $`m`$ valid moments, $`k`$ parameters, PD $`S`$, full-column-rank $`G`$, consistent efficient weighting, and the regular expansion above. Let $`Z=S^{-1/2}\sqrt n\bar g_n(\theta_0)\Rightarrow N_m(0,I)`$. Define $`A=S^{-1/2}G`$ and its orthogonal projector $`P_A=A(A^\top A)^{-1}A^\top`$. Substituting the efficient estimator expansion gives

```math
S^{-1/2}\sqrt n\bar g_n(\widehat\theta)
=(I-P_A)Z+o_p(1).
```

Therefore $`J=n\bar g_n(\widehat\theta)^\top\widehat S^{-1}\bar g_n(\widehat\theta)\Rightarrow Z^\top(I-P_A)Z\sim\chi_{m-k}^2`$. The last matrix is an idempotent projector of rank $`m-k`$. The statistic tests joint compatibility of the moment restrictions, not which instrument is invalid if it rejects.

#### Problem EC-24: stacking OLS and IV moments

Source: [MIT 14.382, Lecture 4, Problems, 1](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/e1aecf6a41da164603168343fcce76bd_MIT14_382S17_lec4.pdf).

Problem. Stack OLS and IV moments to study the difference between their coefficient estimates under a common valid null, retaining their cross covariance.

Solution. In the scalar centered setting, let $`\gamma`$ be the projection slope and $`\alpha`$ the IV slope. Stack

```math
g_i(\gamma,\alpha)=
\binom{D_i(Y_i-D_i\gamma)}
{Z_i(Y_i-D_i\alpha)}.
```

The Jacobian is diagonal with entries $`-E[D^2]`$ and $`-E[ZD]`$. Under iid finite moments, nonzero identifying entries, and a moment CLT, the stacked estimator is jointly asymptotically normal with covariance $`V=G^{-1}SG^{-\top}`$. For $`c=(-1,1)^\top`$, $`\sqrt n[(\widehat\alpha-\widehat\gamma)-(\alpha-\gamma)]\Rightarrow N(0,c^\top Vc)`$. A Wald statistic for equality is $`n(\widehat\alpha-\widehat\gamma)^2/(c^\top\widehat Vc)\Rightarrow\chi_1^2`$ when the variance is positive. The cross-covariance of the estimators is included automatically by stacking.

#### Problem EC-25: indirect inference for an MA(2)

Source: [MIT 14.384, sample problem bank, Problem 7](https://ocw.mit.edu/courses/14-384-time-series-analysis-fall-2013/cbfa217670dde974ad60704bcdcc23e6_MIT14_384F13_problems.pdf).

Problem. For an invertible MA(2), derive the mapping from structural parameters to the selected autocovariances and explain identification by those auxiliary moments.

Solution. For $`Y_t=\mu+e_t+\theta_1e_{t-1}+\theta_2e_{t-2}`$ with iid centered Gaussian innovations of variance $`\sigma^2`$, the population autocovariances are

```math
\gamma_0=\sigma^2(1+\theta_1^2+\theta_2^2),\quad
\gamma_1=\sigma^2\theta_1(1+\theta_2),\quad
\gamma_2=\sigma^2\theta_2,\quad\gamma_j=0\ (j\gt 2).
```

As auxiliary statistics use the sample mean and several sample autocovariances. For each trial structural parameter, simulate a long series, compute the same statistics, and minimize their weighted distance from the observed statistics. Equivalently use auxiliary autoregressive coefficients if their binding function identifies the structural parameters.

Restrict to an invertible MA parameterization to avoid observationally equivalent moving-average representations. Use more independent auxiliary statistics than parameters for a specification test. Its covariance must account for sampling variation and, for finite simulation length, simulation error as well. Under identification, differentiability, suitable stationary CLTs, and consistent covariance estimation, a minimized-distance statistic has a chi-square limit with degrees of freedom equal to the number of independent auxiliary restrictions minus the parameter dimension.

### 5.5 Additional solved problems

#### Problem EC-26: conditional moments imply unconditional moments

Constructed problem. If $`E[u\mid Z]=0`$, show $`E[h(Z)u]=0`$ for integrable products.

Solution. The tower property gives $`E[h(Z)u]=E[h(Z)E[u\mid Z]]=0`$. Thus each chosen instrument function $`h`$ supplies an unconditional moment. The converse for a finite selected list of functions need not hold: satisfying a few unconditional correlations does not enforce a zero conditional mean at every $`Z`$. Richer moment sets can carry more information but require appropriate estimation and identification conditions.

#### Problem EC-27: a moment that identifies two roots

Constructed problem. Let $`g(\theta)=E[Y]-\theta^2`$ and suppose $`EY=4`$. Discuss identification.

Solution. The zero set is $`\lbrace -2,2\rbrace`$. At either root the derivative $`-2\theta`$ is nonzero, giving local identification, but globally the parameter is ambiguous. Restricting $`\theta\gt 0`$ selects two. The sample objective $`(\bar Y-\theta^2)^2`$ can likewise have two minima. A local solver returning one minimum does not prove that the model uniquely identifies it.

#### Problem EC-28: a J-test cannot validate an exactly identified model

Constructed problem. Explain the overidentification degrees of freedom when $`m=k`$.

Solution. Under exact identification, the sample moments can often be set exactly to zero, so the minimized criterion is zero regardless of whether the structural interpretation is credible. The asymptotic residual projector has rank $`m-k=0`$. There is no overidentifying restriction left to test. Instrument validity must then rest on external model reasoning or other evidence.

#### Problem EC-29: misspecified GMM

Constructed problem. If no $`\theta`$ makes $`E[g(\theta)]=0`$, describe the target.

Solution. With a fixed limiting weight $`W`$, the estimator can converge to the pseudo-true minimizer of $`E[g(\theta)]^\top W E[g(\theta)]`$, assuming a uniform LLN and a uniquely separated minimum. Changing $`W`$ can change that target. The correctly specified covariance and chi-square formulas need not carry over, because the moment vector at the target is not zero. A small numerical objective is not itself a proof that misspecification is absent.

#### Problem EC-30: scaling a moment

Constructed problem. Replace moment vector $`g`$ by $`Ag`$ for an invertible matrix $`A`$. Determine the corresponding efficient weight.

Solution. The transformed covariance is $`ASA^\top`$, whose inverse is $`A^{-\top}S^{-1}A^{-1}`$. The efficient quadratic becomes $`(Ag)^\top A^{-\top}S^{-1}A^{-1}(Ag)=g^\top S^{-1}g`$. Thus properly transformed efficient GMM is invariant to an invertible rescaling or recombination of the moments. Keeping the old numerical weight after rescaling generally changes the criterion.

## 6. Panel data and treatment comparisons

### 6.1 Theory

Panel data observe the same units across time. In $`Y_{it}=\alpha_i+X_{it}^\top\beta+u_{it}`$, the unit effect $`\alpha_i`$ captures time-invariant unobserved differences. Fixed effects allow $`\alpha_i`$ to be correlated with regressors. They remove it by subtracting within-unit time averages:

```math
Y_{it}-\bar Y_i=(X_{it}-\bar X_i)^\top\beta+(u_{it}-\bar u_i).
```

This transformation does not remove time-varying omitted variables. A standard sufficient exogeneity condition is $`E[u_{it}\mid X_{i1},\ldots,X_{iT},\alpha_i]=0`$ for every date. Predetermined regressors or lagged outcomes can require a different estimator.

Random effects models impose restrictions such as $`E[\alpha_i\mid X_{i1},\ldots,X_{iT}]=E\alpha_i`$ along with an error covariance model. Their efficiency can improve if these restrictions hold, but the fixed-effects model does not grant them automatically.

Difference-in-differences compares outcome changes across treated and untreated groups. With two periods and no anticipation, its causal interpretation requires parallel trends in untreated potential outcomes for the two groups. Similar observed pretrends are evidence about this assumption, not a proof of the unobserved post-treatment counterfactual.

### 6.2 Derivations and proofs

#### The difference-in-differences estimand

Let $`D=1`$ mark the group treated only after period zero. Write $`Y_1(1),Y_1(0)`$ for treated and untreated potential outcomes. The observed DID is

```math
\Delta=[E(Y_1\mid D=1)-E(Y_0\mid D=1)]
-[E(Y_1\mid D=0)-E(Y_0\mid D=0)].
```

Add and subtract $`E[Y_1(0)\mid D=1]`$. Under no treatment in the control group and no anticipation,

```math
\Delta=E[Y_1(1)-Y_1(0)\mid D=1]
+\lbrace E[Y_1(0)-Y_0(0)\mid D=1]
-E[Y_1(0)-Y_0(0)\mid D=0]\rbrace .
```

The first term is the average treatment effect on treated units. Parallel trends sets the brace to zero. A violation of size $`b`$ adds exactly $`b`$ to the DID estimand.

#### Two-period fixed effects and differencing

With $`T=2`$, demeaning gives deviations $`\pm(Y_{i2}-Y_{i1})/2`$, and the same for regressors. The two demeaned observations carry the same slope information as the single first-difference equation $`\Delta Y_i=\Delta X_i^\top\beta+\Delta u_i`$. With more periods, ordinary first differences and within demeaning generally produce different weighting and efficiency.

#### Within transformation needs strict exogeneity

In $`y_{it}=x_{it}^\top\beta+\alpha_i+u_{it}`$, within transformation removes $`\alpha_i`$ exactly. It also replaces the error by $`u_{it}-\bar u_i`$ and the regressor by $`x_{it}-\bar x_i`$. The transformed regressor depends on the entire observed regressor path.

A sufficient exogeneity condition is

```math
E[u_{it}\mid x_{i1},\ldots,x_{iT},\alpha_i]=0.
```

Contemporaneous orthogonality $`E[x_{it}u_{it}]=0`$ alone is weaker. Future regressors may respond to today's error, making the within-transformed regressor correlated with the transformed error.

#### Example EC-E24: feedback biases a two-period within regression

Let $`x_{i1}=0`$ and $`x_{i2}=u_{i1}`$, with independent mean-zero errors of common variance $`\sigma^2`$ across dates and individuals. At each date the regressor is uncorrelated with that date's error.

For two dates, within estimation equals first-difference estimation. The differenced regressor is $`\Delta x_i=u_{i1}`$ and the differenced error is $`u_{i2}-u_{i1}`$. Their covariance is $`-\sigma^2`$, while the regressor variance is $`\sigma^2`$. Therefore the slope converges to $`\beta-1`$.

Removing the fixed effect did not remove feedback endogeneity. The failure can be read directly from the cross-date restriction $`E[x_{i2}u_{i1}]\ne0`$. This is why a full regressor-path condition appears in standard fixed-effects consistency results.

#### Random effects and correlated group heterogeneity

Random-effects GLS assumes a relationship between group heterogeneity and the regressors, not merely that the intercept is random. With $`E[\alpha_i\mid X_i]=0`$ and suitable idiosyncratic errors, the within-group error covariance is $`\sigma_u^2I+\sigma_\alpha^2\mathbf1\mathbf1^\top`$.

For a balanced panel, GLS can be expressed through quasi-demeaning with factor

```math
\theta=1-\sqrt{\frac{\sigma_u^2}{\sigma_u^2+T\sigma_\alpha^2}}.
```

When the group variance is zero, the factor is zero and the transformation reduces to pooled OLS. As the group component dominates, the factor approaches one and the transformation approaches full demeaning.

Including group means of time-varying regressors can model correlation between heterogeneity and the observed regressor path under a specified conditional-mean structure. It is not an automatic substitute for every possible relationship between $`\alpha_i`$ and $`X_i`$.

For DID, the corresponding causal restriction concerns untreated trends, together with no anticipation and comparable group composition. Differencing removes constant levels but not differential trends. A saturated regression reproduces the four observed means; those means do not reveal the missing treated-group counterfactual.

### 6.3 Worked examples

#### Example EC-E11: a within estimate

Two units have regressor histories $`(1,2)`$ and $`(2,4)`$ and outcomes $`(3,5)`$ and $`(7,11)`$. Their differences are $`\Delta X=(1,2)`$ and $`\Delta Y=(2,4)`$. The common differenced slope is $`(1\cdot2+2\cdot4)/(1^2+2^2)=2`$. Separate level intercepts are removed without needing to estimate them first.

#### Example EC-E12: a DID calculation and a trend violation

Treated-group means move from ten to fifteen; control means move from eight to ten. DID equals $`(15-10)-(10-8)=3`$. If the treated group would have risen by four without treatment while the control group would have risen by two, the trend-bias term is two, so the treatment effect is only one. The arithmetic alone cannot distinguish these counterfactual stories.

### 6.4 Selected exercises

#### Problem EC-31: fixed effects, first differences, and clustering

Source: [MIT 14.382, archived Final 2016, Question 4](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/dc4003ef9928b631c627562e6c11edef_MIT_14_382S17_ExamF16.pdf), the FD/FE construction and standard-error interpretation.

Problem. For the panel model, derive fixed-effects and first-difference transformations and state appropriate covariance treatment. This solution covers the theoretical transformation and clustering parts, rather than the source's empirical data tasks.

Solution. FD subtracts the previous observation for the same school, removing a school-specific level effect. FE subtracts each school's time mean, likewise removing its time-invariant effect. Time effects in the source table absorb changes shared across schools. Each reported coefficient measures a partial association in the corresponding transformed regression, conditional on the included current, lagged, and other covariates.

School-level clustering allows residual products from different dates of the same school to enter the covariance estimate. It does not treat each school-year as an independent observation. The standard errors in parentheses use this cluster construction, while bracketed entries use the stated bootstrap procedure. Neither transformation alone establishes a causal spending effect; time-varying confounding and the exogeneity assumptions still need justification.

#### Problem EC-32: treatment allocation under a budget

Source: [MIT 14.32, Problem Set I, B.2(a–c)](https://ocw.mit.edu/courses/14-32-econometrics-spring-2007/300124445f4e4e08cd6a4a9fb610abf2_psy07n01.pdf).

Problem. Under equal outcome variances and a fixed budget with per-person cost $`c_0`$ plus treatment cost $`c_T`$, find the treatment share minimizing the variance of a difference in means and explain the zero-additional-cost case.

Solution. Under the stated no-effect planning assumption, let both groups' outcome variance be $`v`$. With treated fraction $`p`$, the variance of the difference in means is $`v/[np]+v/[n(1-p)]=v/[np(1-p)]`$. For fixed $`n`$, maximize $`p(1-p)`$, giving $`p=1/2`$.

Now let baseline collection cost be $`c_0\gt 0`$ per person and additional treatment cost be $`c_T\ge0`$. With budget $`R`$, the continuous relaxation gives $`n=R/(c_0+c_Tp)`$. The variance becomes $`v(c_0+c_Tp)/[Rp(1-p)]`$. Differentiating the fraction gives numerator $`c_Tp^2+2c_0p-c_0`$. Its interior zero is

```math
p^*=\frac{\sqrt{c_0}}{\sqrt{c_0}+\sqrt{c_0+c_T}}.
```

It equals one half if treatment has no extra cost and is smaller otherwise. The variance diverges at either endpoint, so the unique stationary point is the global minimum. Actual sample sizes require integer rounding within the budget. Planning at no effect is useful for designing a test near its null and avoids pretending the unknown treatment effect is already known.

### 6.5 Additional solved problems

#### Problem EC-33: a time-invariant regressor

Constructed problem. Explain why a unit characteristic constant across all dates cannot be separately estimated in an unrestricted unit fixed-effects regression.

Solution. If $`Z_{it}=Z_i`$, then $`Z_{it}-\bar Z_i=0`$ for every date. Its column disappears after within transformation. Equivalently, its level effect can be absorbed into $`\alpha_i`$. Identification requires within-unit variation or restrictions linking unit effects to that characteristic. A larger number of dates does not create variation in a variable that is constant by definition.

#### Problem EC-34: serially correlated differenced errors

Constructed problem. If $`u_{it}`$ are independent over time with variance $`\sigma^2`$, compute covariance of neighboring differences.

Solution. Let $`\Delta u_t=u_t-u_{t-1}`$. Its variance is $`2\sigma^2`$. The neighboring covariance is

```math
E[(u_t-u_{t-1})(u_{t-1}-u_{t-2})]=-\sigma^2.
```

All other products vanish by independence and mean zero. Thus differencing independent level errors creates MA(1)-type serial dependence. Using independent-error standard errors on the differenced regression is not automatically justified.

#### Problem EC-35: a lagged dependent variable after demeaning

Constructed problem. Explain the source of bias in a short-panel fixed-effects autoregression.

Solution. In $`Y_{it}=\alpha_i+\rho Y_{i,t-1}+u_{it}`$, the lagged outcome includes earlier innovations. Demeaning the regressor subtracts a time average that contains terms affected by current and other errors, while demeaning the error also uses all dates. These transformed terms become correlated even if the original innovation is independent of the past. With fixed $`T`$ and many units, that correlation need not vanish. Removing unit effects is therefore insufficient to justify ordinary within OLS for this dynamic model.

#### Problem EC-36: a treatment indicator before and after

Constructed problem. Show that the interaction coefficient in a saturated two-group, two-period regression equals DID.

Solution. Write the four fitted cell means as $`a`$, $`a+b`$, $`a+c`$, and $`a+b+c+\tau`$ for control-before, treated-before, control-after, and treated-after. Their within-group changes are $`c`$ and $`c+\tau`$. Subtracting the changes gives $`\tau`$. In a saturated regression with group, period, and their interaction, the fitted cell means equal the observed cell means, so the interaction coefficient is the sample DID.

#### Problem EC-37: heterogeneous effects in a before-after comparison

Constructed problem. Explain why a treated group's raw before-after change need not equal its treatment effect.

Solution. It equals the treatment effect plus the group's untreated outcome trend. A control group is used to estimate that missing trend under a parallel-trends assumption. If treatment effects differ across units, the two-period DID targets the average effect for the treated group under the stated assumptions; it need not equal an effect for the entire population. The target population must be named before interpreting the coefficient.

#### Problem EC-38: random effects and correlated unit effects

Constructed problem. Let $`X_{it}=\alpha_i+v_{it}`$ with independent centered $`v_{it}`$ and $`\alpha_i`$. Show the problem with imposing exogeneity of unit effects.

Solution. $`\mathop{\mathrm{Cov}} (X_{it},\alpha_i)=\mathop{\mathrm{Var}} (\alpha_i)`$, positive when unit effects vary. A random-effects assumption requiring this covariance to be zero is violated. Within demeaning removes $`\alpha_i`$ from both sides of the slope equation, while a model that combines within and between variation under the false independence restriction can be inconsistent. The preferred estimator depends on which assumptions can be defended.

## 7. Stationary and nonstationary time series

### 7.1 Theory

A weakly stationary process has constant finite mean and autocovariances depending only on the lag: $`\gamma_k=\mathop{\mathrm{Cov}} (Y_t,Y_{t-k})`$. Strict stationarity requires invariance of every finite-dimensional joint distribution to a time shift. Finite-second-moment strict stationarity implies weak stationarity, but weak stationarity alone does not control all higher-order distributions.

White noise has mean zero, constant variance, and zero nonzero-lag covariance. It need not be independent unless that stronger condition is stated. An innovation relative to a filtration has conditional mean zero given past information, which is useful for forecasting.

For a causal AR(1),

```math
Y_t=c+\phi Y_{t-1}+\varepsilon_t,\qquad \vert \phi\vert \lt 1,
```

the stationary mean is $`\mu=c/(1-\phi)`$. The centered process has an infinite moving-average representation. Persistence is high when $`\phi`$ is close to one; a fixed coefficient close to one and an exact unit root have different asymptotic behavior.

An ARMA model combines autoregressive feedback and finite moving-average innovations. Its causal stationary representation requires autoregressive roots outside the unit circle under the convention $`1-\phi_1z-\cdots-\phi_pz^p`$. Invertibility of the MA part concerns recovery of innovations from observed history and requires the corresponding MA roots outside the circle.

### 7.2 Derivations and proofs

#### Stationary AR(1)

Iterate the centered recursion:

```math
Y_t-\mu=\sum_{j=0}^{k-1}\phi^j\varepsilon_{t-j}
+\phi^k(Y_{t-k}-\mu).
```

For a causal stationary process with finite variance, the last term tends to zero in $`L^2`$ when $`\vert \phi\vert \lt 1`$. The innovation series converges in $`L^2`$ because $`\sum_j\phi^{2j}\lt \infty`$. If innovations are uncorrelated with variance $`\sigma_\varepsilon^2`$,

```math
\gamma_0=\frac{\sigma_\varepsilon^2}{1-\phi^2},\qquad
\gamma_h=\phi^{\vert h\vert }\gamma_0.
```

The causal restriction matters. With $`\vert \phi\vert \gt 1`$, one can construct an anticipative stationary solution from future shocks, but it does not have the usual innovation and forecasting interpretation used here.

#### Multistep forecast

If future innovations have conditional mean zero given current information $`\mathcal F_t`$, then

```math
E[Y_{t+h}\mid\mathcal F_t]=\mu+\phi^h(Y_t-\mu).
```

The forecast error is $`\sum_{j=0}^{h-1}\phi^j\varepsilon_{t+h-j}`$. With conditionally uncorrelated constant-variance innovations, its conditional variance is $`\sigma_\varepsilon^2\sum_{j=0}^{h-1}\phi^{2j}`$. Gaussian innovations make the forecast distribution Gaussian; the mean and variance calculations themselves do not require Gaussianity.

#### Unit roots and stochastic trends

At $`\phi=1`$, $`Y_t=Y_0+ct+\sum_{j=1}^t\varepsilon_j`$. With deterministic $`Y_0`$ and independent variance-$`\sigma^2`$ increments, the variance is $`t\sigma^2`$, which grows over time. Differencing removes this unit root: $`\Delta Y_t=c+\varepsilon_t`$.

Regressing unrelated random walks in levels can produce large fitted associations whose usual stationary-regression $`t`$ reference laws fail. If a particular linear combination of integrated variables is stationary, they are cointegrated. That is a separate restriction, not a reason to treat every levels regression as stationary.

#### MA invertibility identifies the innovation representation

For $`X_t=\epsilon_t+\theta\epsilon_{t-1}`$ with white-noise variance $`\sigma^2`$,

```math
\gamma_0=\sigma^2(1+\theta^2),\qquad
\gamma_1=\theta\sigma^2,\qquad \gamma_k=0\ (\vert k\vert \gt 1).
```

The lag-one correlation is $`\theta/(1+\theta^2)`$. Replacing $`\theta`$ by $`1/\theta`$ and the innovation variance by $`\theta^2\sigma^2`$ leaves the covariance function unchanged. Under Gaussian innovations, it leaves the whole process law unchanged.

The invertibility convention $`\vert \theta\vert \lt 1`$ chooses the representation whose innovations can be recovered causally:

```math
\epsilon_t=\sum_{j=0}^\infty(-\theta)^jX_{t-j}.
```

Absolute geometric decay ensures the appropriate mean-square convergence. Stationarity of this finite MA does not require invertibility; invertibility addresses recovery and uniqueness of the chosen innovation representation.

#### A unit root changes the regression scaling

For a random walk $`X_t=X_{t-1}+\epsilon_t`$ with $`X_0=0`$, regress $`X_t`$ on $`X_{t-1}`$ without an intercept. Then

```math
\hat\rho-1=
\frac{\sum_tX_{t-1}\epsilon_t}{\sum_tX_{t-1}^2},\qquad
2\sum_tX_{t-1}\epsilon_t=X_n^2-\sum_t\epsilon_t^2.
```

Under an iid finite-variance functional central limit theorem, $`X_{\lfloor ns\rfloor}/(\sigma\sqrt n)`$ converges in process law to Brownian motion $`W_s`$. The denominator divided by $`n^2`$ converges to $`\sigma^2\int_0^1W_s^2ds`$, while the numerator divided by $`n`$ converges to $`\tfrac12\sigma^2(W_1^2-1)`$. Consequently

```math
n(\hat\rho-1)\Rightarrow
\frac{W_1^2-1}{2\int_0^1W_s^2ds}.
```

The rate is $`n`$ rather than $`\sqrt n`$, and the limit is nonnormal. The stated no-intercept specification matters; adding deterministic terms changes the limit. The functional limit theorem is an advanced input from the time-series source, while the displayed regression identity explains why its limit appears.

#### Example EC-E25: a common stochastic trend

Let $`X_t=W_t+u_t`$ and $`Y_t=2W_t+v_t`$, where $`W_t`$ is a random walk and $`u_t,v_t`$ form a jointly stationary process. Both levels contain a stochastic trend, but

```math
Y_t-2X_t=v_t-2u_t
```

is stationary. The variables are cointegrated with vector $`(-2,1)`$. A regression in differences discards this stationary long-run relation unless it is included through an error-correction term. In a vector error-correction representation, a reduced-rank matrix $`\alpha\beta^\top`$ links stationary combinations $`\beta^\top X_{t-1}`$ to subsequent changes. Neither nonstationarity alone nor a high level correlation proves cointegration.

### 7.3 Worked examples

#### Example EC-E13: a stable AR forecast

Let $`c=1`$, $`\phi=1/2`$, $`\sigma_\varepsilon^2=3`$, and $`Y_t=6`$. The stationary mean is two. The two-step forecast is $`2+(1/2)^2(6-2)=3`$, and its error variance is $`3(1+1/4)=15/4`$. The unconditional variance is $`3/(1-1/4)=4`$, the limiting forecast variance as the horizon grows.

#### Example EC-E14: a random-walk forecast

For $`Y_t=Y_{t-1}+\varepsilon_t`$, current value five, and innovation variance two, the $`h`$-step forecast is five and its error variance is $`2h`$. The forecast does not return to a fixed mean, and its uncertainty grows without a stationary ceiling.

#### Example EC-E15: a mean-reversion half-life

For $`0\lt \phi\lt 1`$, a conditional deviation shrinks by $`\phi^h`$. Setting $`\phi^h=1/2`$ gives $`h=\log(1/2)/\log\phi`$. At $`\phi=0.9`$ the half-life is about $`6.58`$ periods. This calculation refers to the conditional mean deviation; realized paths continue receiving new shocks.

### 7.4 Selected exercises

#### Problem EC-39: AR coefficients and an MA representation

Source: [MIT 14.384, sample problem bank, Problem 1](https://ocw.mit.edu/courses/14-384-time-series-analysis-fall-2013/cbfa217670dde974ad60704bcdcc23e6_MIT14_384F13_problems.pdf).

Problem. For a stable AR($`p`$) process, derive the recurrence obeyed by the coefficients in its causal moving-average representation.

Solution. Use $`\phi(z)=1-\sum_{j=1}^p a_jz^j`$ and $`\psi(z)=\sum_{k\ge0}\psi_kz^k`$, with $`\psi_0=1`$ and $`\psi_k=0`$ for negative indices. The representation implies $`\phi(z)\psi(z)=1`$. Equating each coefficient of $`z^k`$ gives

```math
\psi_k-\sum_{j=1}^{\min(p,k)}a_j\psi_{k-j}=0,\qquad k\ge1.
```

For $`k\ge p`$, this is the stated order-$`p`$ homogeneous recursion using all $`p`$ lags. For AR(1), $`\psi_k=a_1^k`$. Causality and convergence of the representation require the stable-root conditions; a formal coefficient identity alone does not establish stationarity.

#### Problem EC-40: subsampling at a fixed stable AR coefficient

Source: [MIT 14.384, sample problem bank, Problem 5(a)](https://ocw.mit.edu/courses/14-384-time-series-analysis-fall-2013/cbfa217670dde974ad60704bcdcc23e6_MIT14_384F13_problems.pdf).

Problem. For a fixed stable AR(1) coefficient, establish the asymptotic validity of a subsampling confidence interval under the source's subsample-quantile consistency assumption, with block size tending to infinity and block-to-sample ratio tending to zero.

Solution. For a fixed $`0\lt \rho_0\lt 1`$, assume the usual stationary AR estimator and standard-error limits, and the source's stated consistency of subsample quantiles. The full-sample $`t`$ statistic tends to $`N(0,1)`$. In a block of length $`b_T\to\infty`$, the statistic centered at $`\rho_0`$ has the same limit. Replacing that center by the full-sample estimate adds a term of order $`\sqrt{b_T}(\widehat\rho_T-\rho_0)=O_p(\sqrt{b_T/T})=o_p(1)`$. Thus its subsampling quantiles approach normal quantiles and coverage approaches $`1-\alpha`$. This fixed-stable-parameter argument does not establish uniform validity as $`\rho_0`$ approaches one.

### 7.5 Additional solved problems

#### Problem EC-41: an MA(1) covariance

Constructed problem. Let $`Y_t=\varepsilon_t+\theta\varepsilon_{t-1}`$ with iid centered innovations of variance $`\sigma^2`$.

Solution. The variance is $`\sigma^2(1+\theta^2)`$. At lag one, the only shared innovation is $`\varepsilon_{t-1}`$, with product coefficient $`\theta`$, so $`\gamma_1=\theta\sigma^2`$. For $`\vert h\vert \gt 1`$ there are no shared innovations and $`\gamma_h=0`$. The process is stationary for every finite $`\theta`$. Invertibility requires $`\vert \theta\vert \lt 1`$ under the usual chosen representation, a distinct issue from stationarity.

#### Problem EC-42: centering an AR model

Constructed problem. Explain why replacing $`Y_t`$ by $`Y_t-c`$ does not generally remove the intercept in an AR(1).

Solution. The stationary mean solves $`\mu=c+\phi\mu`$, so it is $`c/(1-\phi)`$, not $`c`$ unless $`\phi=0`$. Substituting $`Z_t=Y_t-\mu`$ gives $`Z_t=\phi Z_{t-1}+\varepsilon_t`$. Subtracting $`c`$ instead leaves intercept $`\phi c`$, so the claimed centered recursion would be incorrect.

#### Problem EC-43: a forecast with estimated parameters

Constructed problem. Distinguish innovation uncertainty from parameter uncertainty in $`\widehat c+\widehat\phi Y_t`$.

Solution. With known parameters, the one-step forecast error is just the next innovation. With estimated parameters, the error also contains $`(c-\widehat c)+(\phi-\widehat\phi)Y_t`$. Its variance depends on the estimator covariance and the current state, and may include cross terms if the estimation and forecast errors are not suitably separated. A plug-in innovation variance omits this parameter-estimation contribution.

#### Problem EC-44: a deterministic trend

Constructed problem. Compare $`Y_t=a+bt+\varepsilon_t`$ with a random walk with drift.

Solution. Subtracting $`a+bt`$ from the first model leaves the stationary error process when $`\varepsilon_t`$ is stationary. A random walk retains accumulated past innovations after its deterministic drift is removed, so its variance still grows with time. Detrending and differencing address different structures. Choosing between them is a modeling and inference question.

#### Problem EC-45: a simple cointegrated pair

Constructed problem. Let $`X_t`$ be a random walk and $`Y_t=2X_t+v_t`$, where $`v_t`$ is stationary.

Solution. Both level series inherit the stochastic trend from $`X_t`$, but $`Y_t-2X_t=v_t`$ is stationary. The vector $`(1,-2)`$ applied to $`(Y_t,X_t)`$ is a cointegrating vector. A regression of these levels has a different statistical theory from unrelated random walks. The existence of one stationary combination does not make the two-dimensional level process stationary.

#### Problem EC-46: negative autoregressive persistence

Constructed problem. For $`\phi=-1/2`$, describe the conditional mean forecast path around its stationary mean.

Solution. A current deviation $`d`$ has forecast deviations $`-d/2,d/4,-d/8,\ldots`$. Magnitudes shrink geometrically while signs alternate. The stationary variance formula still uses $`\phi^2`$ and remains positive. A positive-coefficient half-life formula based directly on $`\log\phi`$ is not real-valued; for the envelope of magnitudes use $`\log(1/2)/\log\vert \phi\vert`$ instead.

[Contents](README.md) · [Part 2](04-econometrics-and-time-series-part-2.md)
