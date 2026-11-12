# Econometrics and Time Series — Part 2

The same regression formula can describe a projection, estimate a statistical model, or identify a causal effect. The algebra may be identical, but the assumptions supporting those interpretations are different. These notes state the interpretation before using an estimator.

Sources:

- Joshua Angrist, [MIT 14.32, Econometrics, Spring 2007](https://ocw.mit.edu/courses/14-32-econometrics-spring-2007/).
- Victor Chernozhukov and Iván Fernández-Val, [MIT 14.382, Econometrics, Spring 2017](https://ocw.mit.edu/courses/14-382-econometrics-spring-2017/), lecture notes, their numbered problems, and the archived exams included with the course.
- Anna Mikusheva, [MIT 14.384, Time Series Analysis, Fall 2013](https://ocw.mit.edu/courses/14-384-time-series-analysis-fall-2013/), lecture notes and the explicitly labeled bank of sample problems.

Source problems below select theoretical questions or identified subparts. Numerical examples use small constructed data, so they can be checked without software or an external dataset. The time-series problem bank is a sample collection, not a complete set of the original assignments.

[Contents](README.md) · [Part 1](04-econometrics-and-time-series-part-1.md)

Contents:

- [8. Multivariate forecasting and state-space models](#8-multivariate-forecasting-and-state-space-models)
- [9. Further worked problems](#9-further-worked-problems)

## 8. Multivariate forecasting and state-space models

### 8.1 Theory

A VAR(1) is $`Y_t=c+AY_{t-1}+\varepsilon_t`$ for a vector-valued process. If the spectral radius of $`A`$ is below one, powers $`A^h`$ vanish and a causal stationary solution exists under finite-variance innovations. Its mean is $`\mu=(I-A)^{-1}c`$.

The $`h`$-step conditional forecast is $`\mu+A^h(Y_t-\mu)`$, and its innovation forecast-error covariance is

```math
V_h=\sum_{j=0}^{h-1}A^j\Sigma_\varepsilon(A^j)^\top.
```

These reduced-form equations describe predictive dynamics. A structural impulse response requires a model for contemporaneous shocks, such as a justified decomposition $`\varepsilon_t=Bu_t`$. Different choices of $`B`$ can yield different structural interpretations for the same reduced-form covariance.

Out-of-sample evaluation fits models using only information available before each forecast. An expanding window keeps all prior observations; a rolling window keeps only a fixed number of recent observations. Feature construction, preprocessing, lag selection, and model choice must also obey the information cutoff.

A linear Gaussian state-space model separates a hidden state from its noisy observations. The Kalman filter alternates a prediction step with an update using the new observation's surprise. It gives the conditional mean and covariance under the correctly specified linear Gaussian model; more general settings may retain a best-linear-prediction interpretation under weaker assumptions.

### 8.2 Derivations and proofs

#### VAR forecasts and innovations

Iterating $`Y_{t+1}-\mu=A(Y_t-\mu)+\varepsilon_{t+1}`$ gives

```math
Y_{t+h}-\mu=A^h(Y_t-\mu)
+\sum_{j=0}^{h-1}A^j\varepsilon_{t+h-j}.
```

Future conditional means zero remove the sum from the forecast. Conditional independence or appropriate uncorrelatedness of future innovations makes the covariance of the sum equal to the stated sum of covariance matrices. If the spectral radius is below one, the infinite version converges and solves $`\Gamma=A\Gamma A^\top+\Sigma_\varepsilon`$.

#### Scalar Kalman update

Suppose a prior state estimate is $`m`$ with variance $`P`$, and observe $`y=\alpha+v`$ with independent $`v\sim N(0,R)`$. The joint covariance of $`(\alpha,y)`$ has entries $`P,P,P+R`$. The conditional Gaussian formula gives

```math
K=\frac P{P+R},\qquad
m^+=m+K(y-m),\qquad
P^+=(1-K)P.
```

The gain is larger when the prior is uncertain and smaller when the measurement is noisy. For state dynamics $`\alpha_{t+1}=F\alpha_t+\eta_{t+1}`$, predict next mean $`Fm^+`$ and variance $`F^2P^++Q`$ before the next measurement.

#### A Kalman update is Gaussian conditioning

Suppose the state has prior $`x\mid\mathcal F_{t-1}\sim N(a,P)`$ and the observation is $`y=Hx+\eta`$, with independent $`\eta\sim N(0,R)`$. The innovation is $`e=y-Ha`$. Its covariance is $`S=HPH^\top+R`$ and its covariance with the state is $`PH^\top`$.

The joint Gaussian conditioning formula therefore gives

```math
K=PH^\top S^{-1},\qquad
m=a+Ke,\qquad
P_{\mathrm{post}}=P-KSK^\top.
```

An equivalent covariance expression is

```math
P_{\mathrm{post}}=(I-KH)P(I-KH)^\top+KRK^\top.
```

This Joseph form displays positive semidefiniteness and is often useful in rounded arithmetic. Substitution of $`KS=PH^\top`$ proves its equivalence to the simpler expression.

With state transition $`x_{t+1}=Fx_t+\xi_{t+1}`$ and independent noise covariance $`Q`$, the next prior is $`a_{t+1}=Fm`$ and $`P_{t+1}=FP_{\mathrm{post}}F^\top+Q`$. Filtering uses information through the current observation; smoothing uses later observations as well. A smoothed state estimate is therefore not automatically an admissible real-time forecasting input.

#### Forecast comparison needs an available information set

At origin $`t`$, a forecast for $`Y_{t+h}`$ must be measurable with respect to information actually available then. This includes fitted coefficients, transformations, variable selection, and data-release vintages. A feature calculated from the full future sample can violate timing even when the regression was fitted only on earlier rows.

An expanding window uses all observations available by each origin. A rolling window uses only a recent fixed or varying span. The choice trades estimation precision against adaptation to changing relationships. Neither choice guarantees better forecasts without an evaluation under a specified process or dataset.

For squared-error loss, define $`d_t=e_{A,t}^2-e_{B,t}^2`$. A positive mean favors model B. Under a stationary weak-dependence CLT with positive long-run variance $`\Omega_d`$,

```math
\frac{\sqrt n\,\bar d}{\sqrt{\hat\Omega_d}}
\Rightarrow N(0,1)
```

under equal expected loss. Multi-step forecast errors may overlap, so an independent-observation variance can be inappropriate. Nested-model or degeneracy cases can require a different reference distribution; this formula is not a universal forecast-comparison test.

#### Example EC-E26: a loss difference with known limiting variance

Suppose a valid evaluation yields $`n=100`$, average loss difference $`0.2`$, and long-run variance estimate one. The statistic is $`\sqrt{100}(0.2)=2`$ and the two-sided normal p-value is about $`0.0455`$.

The direction says model B had lower average squared loss. The inference uses a long-run variance of the loss differences, not a variance of raw outcomes or a training-fit standard error. It applies to the specified forecast procedures and evaluation population, with the stated time-series assumptions.

### 8.3 Worked examples

#### Example EC-E16: a two-variable VAR

Let

```math
A=\left(\begin{smallmatrix}1/2&1/4\\[0pt] 0&1/2\end{smallmatrix}\right)
```

, $`c=0`$, and $`Y_t=(2,4)^\top`$. The one-step forecast is $`(2,2)^\top`$. Since

```math
A^2=\left(\begin{smallmatrix}1/4&1/4\\[0pt] 0&1/4\end{smallmatrix}\right)
```

, the two-step forecast is $`(3/2,1)^\top`$. The second variable predicts the first through the off-diagonal term.

#### Example EC-E17: a measurement update

A state prior has mean ten and variance four. A measurement eleven has noise variance one. Then $`K=4/5`$, updated mean $`10+(4/5)(1)=10.8`$, and variance $`4/5`$. The estimate moves most of the way toward the more precise new measurement without treating it as exact.

#### Example EC-E18: forecast losses

Two models' squared forecast errors over four genuinely held-out dates are $`(1,4,1,4)`$ and $`(2,2,2,2)`$. Their mean losses are $`2.5`$ and two. The average loss difference is $`0.5`$. This is a descriptive sample comparison. Testing equal predictive ability needs the long-run variance of the loss difference, especially for overlapping multistep forecasts.

### 8.4 Selected exercises

#### Problem EC-47: principal-component factor estimation

Source: [MIT 14.384, sample problem bank, Problem 4(a–d)](https://ocw.mit.edu/courses/14-384-time-series-analysis-fall-2013/cbfa217670dde974ad60704bcdcc23e6_MIT14_384F13_problems.pdf).

Problem. Formulate principal-component factor estimation as a constrained least-squares problem and derive the eigenvector solution under the stated normalization.

Solution. For a $`T\times N`$ matrix $`X`$, approximate $`X=F\Lambda^\top+E`$, imposing $`F^\top F/T=I_k`$. At fixed $`F`$, differentiate $`\Vert X-F\Lambda^\top\Vert _F^2`$ to obtain $`\widehat\Lambda=X^\top F/T`$. Substitution gives

```math
\Vert X-F\widehat\Lambda^\top\Vert _F^2
=\Vert X\Vert _F^2-\frac1T\mathop{\mathrm{tr}} (F^\top XX^\top F).
```

Thus choose $`F`$ to maximize the trace under its normalization. If $`u_1,\ldots,u_k`$ are leading unit eigenvectors of $`XX^\top`$, take $`\widehat F=\sqrt T(u_1,\ldots,u_k)`$. Then $`\widehat\Lambda=X^\top\widehat F/T`$, and the common component is $`\widehat F\widehat\Lambda^\top=U_kU_k^\top X`$. Orthogonal rotations of the factors and loadings give the same common component, so the latent coordinates require a normalization and are not uniquely labeled by the observed data alone.

#### Problem EC-48: a constant hidden level

Source: [MIT 14.384, sample problem bank, Problem 8(a–c)](https://ocw.mit.edu/courses/14-384-time-series-analysis-fall-2013/cbfa217670dde974ad60704bcdcc23e6_MIT14_384F13_problems.pdf).

Problem. In a constant latent-level model with Gaussian observation noise and a Gaussian prior, derive recursive and closed-form posterior means and variances, including their long-sample limits.

Solution. The state is constant, $`\alpha_t=\alpha`$, and $`Y_t=\alpha+v_t`$ with iid variance-$`\sigma^2`$ Gaussian noise. Start with normal prior mean $`a`$ and variance $`p_0\gt 0`$. Repeated Gaussian updates add precisions:

```math
P_T^{-1}=p_0^{-1}+T/\sigma^2,\qquad
P_T=\frac{p_0}{1+Tp_0/\sigma^2},
```

```math
m_T=\frac{a/p_0+\sum_{t=1}^TY_t/\sigma^2}
{1/p_0+T/\sigma^2}.
```

The coefficient on each new observation is $`1/(T+\sigma^2/p_0)`$, tending to zero. The prior mean weight tends to zero, while the total data weight tends to one. Conditional on the fixed true state, the LLN gives $`m_T\to\alpha`$, independently of $`a`$ and positive finite $`p_0`$. With diffuse prior $`p_0\to\infty`$, $`m_T=\bar Y_T`$ and $`P_T=\sigma^2/T`$. A zero prior variance would instead encode certainty and prevent learning, so it is excluded from the prior-forgetting statement.

### 8.5 Additional solved problems

#### Problem EC-49: a forecast-comparison standard error

Constructed problem. Let $`d_t`$ be a stationary squared-error loss difference with mean $`\delta`$ and long-run variance $`\Omega_d\gt 0`$.

Solution. Under a suitable dependent-data CLT, $`\sqrt T(\bar d-\delta)\Rightarrow N(0,\Omega_d)`$. Estimate $`\Omega_d`$ by an appropriate HAC estimator and test $`\delta=0`$ using $`\sqrt T\bar d/\sqrt{\widehat\Omega_d}`$. A positive estimate favors the model whose loss was subtracted. Nested-model comparisons, estimation-window choices, and very short samples can change the reference theory; the generic CLT needs to be justified for the chosen evaluation design.

#### Problem EC-50: a timing violation

Constructed problem. A predictor is standardized using the mean and variance of the entire sample before splitting into training and test periods.

Solution. The scale parameters use future observations unavailable at forecast time. Even if no future outcome is inserted directly, the transformed historical predictor depends on later data. At each forecast origin, compute preprocessing quantities from the available training window, apply those fixed quantities to the forecast input, then evaluate the future outcome. Revised historical datasets require an additional distinction between release-date data and final revised data.

#### Problem EC-51: a diagonal VAR covariance

Constructed problem. Let $`A=\mathop{\mathrm{diag}} (a,b)`$ with $`\vert a\vert ,\vert b\vert \lt 1`$ and innovation covariance

```math
\left(\begin{smallmatrix}q&r\\[0pt] r&s\end{smallmatrix}\right)
```

. Solve the stationary covariance.

Solution. The equation $`\Gamma=A\Gamma A^\top+\Sigma_\varepsilon`$ gives $`\Gamma_{11}=a^2\Gamma_{11}+q`$, $`\Gamma_{22}=b^2\Gamma_{22}+s`$, and $`\Gamma_{12}=ab\Gamma_{12}+r`$. Thus the entries are $`q/(1-a^2)`$, $`s/(1-b^2)`$, and $`r/(1-ab)`$. Even with diagonal transition dynamics, contemporaneously correlated innovations can produce correlated state coordinates.

#### Problem EC-52: an impulse response requires a shock definition

Constructed problem. Compare a unit reduced-form innovation in coordinate one with a structural unit shock.

Solution. The first changes the horizon-$`h`$ conditional mean by $`A^he_1`$. If $`\varepsilon=Bu`$ and a structural shock sets $`u=e_j`$, the change is $`A^hBe_j`$. Unless $`B=I`$, these are different perturbations. A Cholesky factor is one possible $`B`$ consistent with the reduced-form covariance, but its ordering restrictions require interpretation.

#### Problem EC-53: an expanding and rolling window

Constructed problem. Explain the tradeoff in fitting a stable regression to a process whose coefficient changes.

Solution. An expanding window reduces sampling noise when a common coefficient applies throughout, but old observations can bias a fit toward an obsolete regime after a change. A rolling window discards old regimes at the cost of using fewer observations. Neither is uniformly better. Compare them at the same forecast origins and with window lengths selected without seeing their future evaluation losses.

#### Problem EC-54: a zero-noise measurement

Constructed problem. In the scalar Kalman update, consider $`R\downarrow0`$ with $`P\gt 0`$, then $`P\downarrow0`$ with $`R\gt 0`$.

Solution. In the first limit $`K=P/(P+R)\to1`$, so the updated state equals the essentially exact measurement and its variance tends to zero. In the second limit $`K\to0`$, so a perfectly known prior state barely moves in response to a noisy observation. If both variances are exactly zero and the values disagree, the assumed probability model is internally inconsistent; the formula's denominator vanishes and should not be used.

<!-- extended-practice-start -->

## 9. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 9.1 Regression arithmetic and the sampling model

A fitted line is an algebraic summary of observed data. A standard error adds a sampling model. A causal interpretation adds further assumptions about how the regressor was generated. These examples compute the first two objects while keeping their assumptions visible.

#### Problem EC-55: four-observation OLS calculation, case 1

Regress

```math
y=\begin{pmatrix}1\\[0pt] 2\\[0pt] 2\\[0pt] 5\end{pmatrix}
```

 on an intercept and $`x=(0,1,2,3)^\top`$. Compute the coefficients, residual variance estimate, and homoskedastic slope standard error.

Solution. The sufficient sums are $`n=4`$, $`\sum x_i=6`$, $`\sum x_i^2=14`$, $`\sum y_i=10`$, and $`\sum x_iy_i=21`$. Thus

```math
X^\top X=\begin{pmatrix}4&6\\[0pt] 6&14\end{pmatrix},
\qquad
(X^\top X)^{-1}=\frac1{20}
\begin{pmatrix}14&-6\\[0pt] -6&4\end{pmatrix}.
```

Multiplying by $`X^\top y`$ gives $`\hat\beta_0=\frac{7}{10}`$ and $`\hat\beta_1=\frac{6}{5}`$. The fitted values are

```math
\begin{pmatrix}\frac{7}{10}\\[0pt] \frac{19}{10}\\[0pt] \frac{31}{10}\\[0pt] \frac{43}{10}\end{pmatrix}
```

 and the residuals are

```math
\begin{pmatrix}\frac{3}{10}\\[0pt] \frac{1}{10}\\[0pt] -\frac{11}{10}\\[0pt] \frac{7}{10}\end{pmatrix}
```

. Their sum and their inner product with $`x`$ are both zero, checking the normal equations.

The residual sum of squares is $`\frac{9}{5}`$. With two fitted coefficients, the residual degrees of freedom are $`4-2=2`$, so

```math
s^2=\frac{\mathop{\mathrm{SSE}} }{n-2}=\frac{9}{10}.
```

Under a correctly specified conditional mean, fixed design or conditioning on $`X`$, and homoskedastic uncorrelated errors, $`\mathop{\mathrm{Var}} (\hat\beta\mid X)=\sigma^2(X^\top X)^{-1}`$. The slope diagonal entry is $`1/5`$, giving estimated slope standard error $`\sqrt{s^2/5}\approx 0.424264`$.

The coefficient calculation itself did not require normal errors. Exact small-sample t inference would additionally use independent Gaussian errors conditional on the design. Heteroskedasticity would change the covariance formula, and endogeneity would undermine estimation of a structural coefficient even if the fitted residuals remain orthogonal to the regressors in sample. Sample orthogonality is produced by minimization; it does not prove the population exogeneity condition.

#### Problem EC-56: four-observation OLS calculation, case 2

Regress

```math
y=\begin{pmatrix}4\\[0pt] 3\\[0pt] 3\\[0pt] 1\end{pmatrix}
```

 on an intercept and $`x=(0,1,2,3)^\top`$. Compute the coefficients, residual variance estimate, and homoskedastic slope standard error.

Solution. The sufficient sums are $`n=4`$, $`\sum x_i=6`$, $`\sum x_i^2=14`$, $`\sum y_i=11`$, and $`\sum x_iy_i=12`$. Thus

```math
X^\top X=\begin{pmatrix}4&6\\[0pt] 6&14\end{pmatrix},
\qquad
(X^\top X)^{-1}=\frac1{20}
\begin{pmatrix}14&-6\\[0pt] -6&4\end{pmatrix}.
```

Multiplying by $`X^\top y`$ gives $`\hat\beta_0=\frac{41}{10}`$ and $`\hat\beta_1=-\frac{9}{10}`$. The fitted values are

```math
\begin{pmatrix}\frac{41}{10}\\[0pt] \frac{16}{5}\\[0pt] \frac{23}{10}\\[0pt] \frac{7}{5}\end{pmatrix}
```

 and the residuals are

```math
\begin{pmatrix}-\frac{1}{10}\\[0pt] -\frac{1}{5}\\[0pt] \frac{7}{10}\\[0pt] -\frac{2}{5}\end{pmatrix}
```

. Their sum and their inner product with $`x`$ are both zero, checking the normal equations.

The residual sum of squares is $`\frac{7}{10}`$. With two fitted coefficients, the residual degrees of freedom are $`4-2=2`$, so

```math
s^2=\frac{\mathop{\mathrm{SSE}} }{n-2}=\frac{7}{20}.
```

Under a correctly specified conditional mean, fixed design or conditioning on $`X`$, and homoskedastic uncorrelated errors, $`\mathop{\mathrm{Var}} (\hat\beta\mid X)=\sigma^2(X^\top X)^{-1}`$. The slope diagonal entry is $`1/5`$, giving estimated slope standard error $`\sqrt{s^2/5}\approx 0.264575`$.

The coefficient calculation itself did not require normal errors. Exact small-sample t inference would additionally use independent Gaussian errors conditional on the design. Heteroskedasticity would change the covariance formula, and endogeneity would undermine estimation of a structural coefficient even if the fitted residuals remain orthogonal to the regressors in sample. Sample orthogonality is produced by minimization; it does not prove the population exogeneity condition.

#### Problem EC-57: four-observation OLS calculation, case 3

Regress

```math
y=\begin{pmatrix}-1\\[0pt] 1\\[0pt] 4\\[0pt] 4\end{pmatrix}
```

 on an intercept and $`x=(0,1,2,3)^\top`$. Compute the coefficients, residual variance estimate, and homoskedastic slope standard error.

Solution. The sufficient sums are $`n=4`$, $`\sum x_i=6`$, $`\sum x_i^2=14`$, $`\sum y_i=8`$, and $`\sum x_iy_i=21`$. Thus

```math
X^\top X=\begin{pmatrix}4&6\\[0pt] 6&14\end{pmatrix},
\qquad
(X^\top X)^{-1}=\frac1{20}
\begin{pmatrix}14&-6\\[0pt] -6&4\end{pmatrix}.
```

Multiplying by $`X^\top y`$ gives $`\hat\beta_0=-\frac{7}{10}`$ and $`\hat\beta_1=\frac{9}{5}`$. The fitted values are

```math
\begin{pmatrix}-\frac{7}{10}\\[0pt] \frac{11}{10}\\[0pt] \frac{29}{10}\\[0pt] \frac{47}{10}\end{pmatrix}
```

 and the residuals are

```math
\begin{pmatrix}-\frac{3}{10}\\[0pt] -\frac{1}{10}\\[0pt] \frac{11}{10}\\[0pt] -\frac{7}{10}\end{pmatrix}
```

. Their sum and their inner product with $`x`$ are both zero, checking the normal equations.

The residual sum of squares is $`\frac{9}{5}`$. With two fitted coefficients, the residual degrees of freedom are $`4-2=2`$, so

```math
s^2=\frac{\mathop{\mathrm{SSE}} }{n-2}=\frac{9}{10}.
```

Under a correctly specified conditional mean, fixed design or conditioning on $`X`$, and homoskedastic uncorrelated errors, $`\mathop{\mathrm{Var}} (\hat\beta\mid X)=\sigma^2(X^\top X)^{-1}`$. The slope diagonal entry is $`1/5`$, giving estimated slope standard error $`\sqrt{s^2/5}\approx 0.424264`$.

The coefficient calculation itself did not require normal errors. Exact small-sample t inference would additionally use independent Gaussian errors conditional on the design. Heteroskedasticity would change the covariance formula, and endogeneity would undermine estimation of a structural coefficient even if the fitted residuals remain orthogonal to the regressors in sample. Sample orthogonality is produced by minimization; it does not prove the population exogeneity condition.

#### Problem EC-58: four-observation OLS calculation, case 4

Regress

```math
y=\begin{pmatrix}2\\[0pt] 1\\[0pt] 3\\[0pt] 6\end{pmatrix}
```

 on an intercept and $`x=(0,1,2,3)^\top`$. Compute the coefficients, residual variance estimate, and homoskedastic slope standard error.

Solution. The sufficient sums are $`n=4`$, $`\sum x_i=6`$, $`\sum x_i^2=14`$, $`\sum y_i=12`$, and $`\sum x_iy_i=25`$. Thus

```math
X^\top X=\begin{pmatrix}4&6\\[0pt] 6&14\end{pmatrix},
\qquad
(X^\top X)^{-1}=\frac1{20}
\begin{pmatrix}14&-6\\[0pt] -6&4\end{pmatrix}.
```

Multiplying by $`X^\top y`$ gives $`\hat\beta_0=\frac{9}{10}`$ and $`\hat\beta_1=\frac{7}{5}`$. The fitted values are

```math
\begin{pmatrix}\frac{9}{10}\\[0pt] \frac{23}{10}\\[0pt] \frac{37}{10}\\[0pt] \frac{51}{10}\end{pmatrix}
```

 and the residuals are

```math
\begin{pmatrix}\frac{11}{10}\\[0pt] -\frac{13}{10}\\[0pt] -\frac{7}{10}\\[0pt] \frac{9}{10}\end{pmatrix}
```

. Their sum and their inner product with $`x`$ are both zero, checking the normal equations.

The residual sum of squares is $`\frac{21}{5}`$. With two fitted coefficients, the residual degrees of freedom are $`4-2=2`$, so

```math
s^2=\frac{\mathop{\mathrm{SSE}} }{n-2}=\frac{21}{10}.
```

Under a correctly specified conditional mean, fixed design or conditioning on $`X`$, and homoskedastic uncorrelated errors, $`\mathop{\mathrm{Var}} (\hat\beta\mid X)=\sigma^2(X^\top X)^{-1}`$. The slope diagonal entry is $`1/5`$, giving estimated slope standard error $`\sqrt{s^2/5}\approx 0.648074`$.

The coefficient calculation itself did not require normal errors. Exact small-sample t inference would additionally use independent Gaussian errors conditional on the design. Heteroskedasticity would change the covariance formula, and endogeneity would undermine estimation of a structural coefficient even if the fitted residuals remain orthogonal to the regressors in sample. Sample orthogonality is produced by minimization; it does not prove the population exogeneity condition.

### 9.2 Omitted-variable bias from population moments

An omitted regressor affects a population regression coefficient through two quantities: its effect on the outcome and its linear association with the retained regressor. The same omitted variable can create upward or downward bias. A covariance matrix must be internally valid before this formula is used.

#### Problem EC-59: short and long population regressions, case 1

Let centered variables satisfy $`Y=2X+(3)Z+u`$, with $`E[Xu]=E[Zu]=0`$ and covariance matrix of $`(X,Z)`$ equal to

```math
\begin{pmatrix}4&2\\[0pt] 2&4\end{pmatrix}
```

. Find the slope when $`Z`$ is omitted and verify the full regression coefficients.

Solution. The regressor covariance matrix is positive definite because $`\mathop{\mathrm{Var}} (X)=4\gt 0`$ and its determinant is $`12\gt 0`$. Thus both the short regression denominator and the full regression inverse are well defined.

The population short-regression slope is

```math
b_{\mathrm{short}}
=\frac{E[XY]}{E[X^2]}
=\beta+\gamma\frac{E[XZ]}{E[X^2]}
=\frac{7}{2}.
```

The omitted-variable term is $`\gamma\mathop{\mathrm{Cov}} (X,Z)/\mathop{\mathrm{Var}} (X)=\frac{3}{2}`$. It vanishes if either the omitted effect is zero or the two regressors are uncorrelated. Its sign depends on both factors.

For the full regression, the covariance vector with the outcome is

```math
E\left[\begin{pmatrix}X\\[0pt] Z\end{pmatrix}Y\right]
=\begin{pmatrix}14\\[0pt] 16\end{pmatrix}.
```

Solving the two population normal equations gives

```math
\begin{pmatrix}b_X\\[0pt] b_Z\end{pmatrix}
=
\begin{pmatrix}4&2\\[0pt] 2&4\end{pmatrix}^{-1}
\begin{pmatrix}14\\[0pt] 16\end{pmatrix}
=\begin{pmatrix}2\\[0pt] 3\end{pmatrix}.
```

The coefficients recover the specified values because the error is orthogonal to both regressors and the covariance matrix is nonsingular.

An iid law of large numbers with suitable moments makes sample OLS converge to these population projections. Calling the short-regression difference a causal bias additionally requires that the full equation have a defensible structural interpretation. The algebra by itself describes population linear prediction. It does not establish which variables should be controlled for in a causal study.

#### Problem EC-60: short and long population regressions, case 2

Let centered variables satisfy $`Y=1X+(-2)Z+u`$, with $`E[Xu]=E[Zu]=0`$ and covariance matrix of $`(X,Z)`$ equal to

```math
\begin{pmatrix}2&1\\[0pt] 1&3\end{pmatrix}
```

. Find the slope when $`Z`$ is omitted and verify the full regression coefficients.

Solution. The regressor covariance matrix is positive definite because $`\mathop{\mathrm{Var}} (X)=2\gt 0`$ and its determinant is $`5\gt 0`$. Thus both the short regression denominator and the full regression inverse are well defined.

The population short-regression slope is

```math
b_{\mathrm{short}}
=\frac{E[XY]}{E[X^2]}
=\beta+\gamma\frac{E[XZ]}{E[X^2]}
=0.
```

The omitted-variable term is $`\gamma\mathop{\mathrm{Cov}} (X,Z)/\mathop{\mathrm{Var}} (X)=-1`$. It vanishes if either the omitted effect is zero or the two regressors are uncorrelated. Its sign depends on both factors.

For the full regression, the covariance vector with the outcome is

```math
E\left[\begin{pmatrix}X\\[0pt] Z\end{pmatrix}Y\right]
=\begin{pmatrix}0\\[0pt] -5\end{pmatrix}.
```

Solving the two population normal equations gives

```math
\begin{pmatrix}b_X\\[0pt] b_Z\end{pmatrix}
=
\begin{pmatrix}2&1\\[0pt] 1&3\end{pmatrix}^{-1}
\begin{pmatrix}0\\[0pt] -5\end{pmatrix}
=\begin{pmatrix}1\\[0pt] -2\end{pmatrix}.
```

The coefficients recover the specified values because the error is orthogonal to both regressors and the covariance matrix is nonsingular.

An iid law of large numbers with suitable moments makes sample OLS converge to these population projections. Calling the short-regression difference a causal bias additionally requires that the full equation have a defensible structural interpretation. The algebra by itself describes population linear prediction. It does not establish which variables should be controlled for in a causal study.

#### Problem EC-61: short and long population regressions, case 3

Let centered variables satisfy $`Y=-1X+(4)Z+u`$, with $`E[Xu]=E[Zu]=0`$ and covariance matrix of $`(X,Z)`$ equal to

```math
\begin{pmatrix}3&-1\\[0pt] -1&2\end{pmatrix}
```

. Find the slope when $`Z`$ is omitted and verify the full regression coefficients.

Solution. The regressor covariance matrix is positive definite because $`\mathop{\mathrm{Var}} (X)=3\gt 0`$ and its determinant is $`5\gt 0`$. Thus both the short regression denominator and the full regression inverse are well defined.

The population short-regression slope is

```math
b_{\mathrm{short}}
=\frac{E[XY]}{E[X^2]}
=\beta+\gamma\frac{E[XZ]}{E[X^2]}
=-\frac{7}{3}.
```

The omitted-variable term is $`\gamma\mathop{\mathrm{Cov}} (X,Z)/\mathop{\mathrm{Var}} (X)=-\frac{4}{3}`$. It vanishes if either the omitted effect is zero or the two regressors are uncorrelated. Its sign depends on both factors.

For the full regression, the covariance vector with the outcome is

```math
E\left[\begin{pmatrix}X\\[0pt] Z\end{pmatrix}Y\right]
=\begin{pmatrix}-7\\[0pt] 9\end{pmatrix}.
```

Solving the two population normal equations gives

```math
\begin{pmatrix}b_X\\[0pt] b_Z\end{pmatrix}
=
\begin{pmatrix}3&-1\\[0pt] -1&2\end{pmatrix}^{-1}
\begin{pmatrix}-7\\[0pt] 9\end{pmatrix}
=\begin{pmatrix}-1\\[0pt] 4\end{pmatrix}.
```

The coefficients recover the specified values because the error is orthogonal to both regressors and the covariance matrix is nonsingular.

An iid law of large numbers with suitable moments makes sample OLS converge to these population projections. Calling the short-regression difference a causal bias additionally requires that the full equation have a defensible structural interpretation. The algebra by itself describes population linear prediction. It does not establish which variables should be controlled for in a causal study.

#### Problem EC-62: short and long population regressions, case 4

Let centered variables satisfy $`Y=3X+(2)Z+u`$, with $`E[Xu]=E[Zu]=0`$ and covariance matrix of $`(X,Z)`$ equal to

```math
\begin{pmatrix}5&0\\[0pt] 0&1\end{pmatrix}
```

. Find the slope when $`Z`$ is omitted and verify the full regression coefficients.

Solution. The regressor covariance matrix is positive definite because $`\mathop{\mathrm{Var}} (X)=5\gt 0`$ and its determinant is $`5\gt 0`$. Thus both the short regression denominator and the full regression inverse are well defined.

The population short-regression slope is

```math
b_{\mathrm{short}}
=\frac{E[XY]}{E[X^2]}
=\beta+\gamma\frac{E[XZ]}{E[X^2]}
=3.
```

The omitted-variable term is $`\gamma\mathop{\mathrm{Cov}} (X,Z)/\mathop{\mathrm{Var}} (X)=0`$. It vanishes if either the omitted effect is zero or the two regressors are uncorrelated. Its sign depends on both factors.

For the full regression, the covariance vector with the outcome is

```math
E\left[\begin{pmatrix}X\\[0pt] Z\end{pmatrix}Y\right]
=\begin{pmatrix}15\\[0pt] 2\end{pmatrix}.
```

Solving the two population normal equations gives

```math
\begin{pmatrix}b_X\\[0pt] b_Z\end{pmatrix}
=
\begin{pmatrix}5&0\\[0pt] 0&1\end{pmatrix}^{-1}
\begin{pmatrix}15\\[0pt] 2\end{pmatrix}
=\begin{pmatrix}3\\[0pt] 2\end{pmatrix}.
```

The coefficients recover the specified values because the error is orthogonal to both regressors and the covariance matrix is nonsingular.

An iid law of large numbers with suitable moments makes sample OLS converge to these population projections. Calling the short-regression difference a causal bias additionally requires that the full equation have a defensible structural interpretation. The algebra by itself describes population linear prediction. It does not establish which variables should be controlled for in a causal study.

### 9.3 Heteroskedastic OLS and known-variance GLS

With fixed design, heteroskedasticity does not bias OLS when errors have conditional mean zero. It does change precision. When the error variances are known, GLS reweights observations by inverse variance and can reduce the estimator's variance. A robust covariance estimate instead changes the estimated uncertainty while retaining the OLS point estimate.

#### Problem EC-63: precision under unequal error variances, case 1

In the through-origin model $`Y_i=\beta x_i+u_i`$, take fixed $`x=(1,2,3)`$, independent mean-zero errors, and $`\mathop{\mathrm{Var}} (u_i)=(1, 4, 9)`$. Derive the OLS and GLS variances.

Solution. The OLS estimator is $`\hat\beta=\sum_i x_iY_i/\sum_i x_i^2`$. Substituting the model yields

```math
\hat\beta-\beta=\frac{\sum_i x_iu_i}{14}.
```

Conditional mean zero proves unbiasedness. Independence makes the variance of the numerator the sum of $`x_i^2\sigma_i^2`$, so

```math
\mathop{\mathrm{Var}} (\hat\beta\mid X)
=\frac{\sum_i x_i^2\sigma_i^2}{(\sum_i x_i^2)^2}
=\frac{1}{2}\approx 0.50000000.
```

This is the scalar sandwich formula. It cannot generally be replaced by one common variance times $`1/14`$.

If the variances are known, minimize $`\sum_i(Y_i-\beta x_i)^2/\sigma_i^2`$. Differentiating gives

```math
\hat\beta_{\mathrm{GLS}}
=\frac{\sum_i x_iY_i/\sigma_i^2}{\sum_i x_i^2/\sigma_i^2},
\qquad
\hat\beta_{\mathrm{GLS}}-\beta
=\frac{\sum_i x_iu_i/\sigma_i^2}{\sum_i x_i^2/\sigma_i^2}.
```

The variance of the numerator is $`\sum_i x_i^2/\sigma_i^2=3`$. Therefore $`\mathop{\mathrm{Var}} (\hat\beta_{\mathrm{GLS}}\mid X)=\frac{1}{3}\approx 0.33333333`$.

Cauchy–Schwarz gives $`(\sum_i x_i^2)^2\le(\sum_i x_i^2\sigma_i^2)(\sum_i x_i^2/\sigma_i^2)`$, proving that the GLS variance is no larger. A heteroskedasticity-robust standard error for OLS estimates the first variance using squared residuals; it does not turn the OLS coefficients into GLS coefficients. Estimating the variance weights creates feasible GLS and introduces an additional estimation step.

#### Problem EC-64: precision under unequal error variances, case 2

In the through-origin model $`Y_i=\beta x_i+u_i`$, take fixed $`x=(1,2,3)`$, independent mean-zero errors, and $`\mathop{\mathrm{Var}} (u_i)=(9, 4, 1)`$. Derive the OLS and GLS variances.

Solution. The OLS estimator is $`\hat\beta=\sum_i x_iY_i/\sum_i x_i^2`$. Substituting the model yields

```math
\hat\beta-\beta=\frac{\sum_i x_iu_i}{14}.
```

Conditional mean zero proves unbiasedness. Independence makes the variance of the numerator the sum of $`x_i^2\sigma_i^2`$, so

```math
\mathop{\mathrm{Var}} (\hat\beta\mid X)
=\frac{\sum_i x_i^2\sigma_i^2}{(\sum_i x_i^2)^2}
=\frac{17}{98}\approx 0.17346939.
```

This is the scalar sandwich formula. It cannot generally be replaced by one common variance times $`1/14`$.

If the variances are known, minimize $`\sum_i(Y_i-\beta x_i)^2/\sigma_i^2`$. Differentiating gives

```math
\hat\beta_{\mathrm{GLS}}
=\frac{\sum_i x_iY_i/\sigma_i^2}{\sum_i x_i^2/\sigma_i^2},
\qquad
\hat\beta_{\mathrm{GLS}}-\beta
=\frac{\sum_i x_iu_i/\sigma_i^2}{\sum_i x_i^2/\sigma_i^2}.
```

The variance of the numerator is $`\sum_i x_i^2/\sigma_i^2=\frac{91}{9}`$. Therefore $`\mathop{\mathrm{Var}} (\hat\beta_{\mathrm{GLS}}\mid X)=\frac{9}{91}\approx 0.09890110`$.

Cauchy–Schwarz gives $`(\sum_i x_i^2)^2\le(\sum_i x_i^2\sigma_i^2)(\sum_i x_i^2/\sigma_i^2)`$, proving that the GLS variance is no larger. A heteroskedasticity-robust standard error for OLS estimates the first variance using squared residuals; it does not turn the OLS coefficients into GLS coefficients. Estimating the variance weights creates feasible GLS and introduces an additional estimation step.

#### Problem EC-65: precision under unequal error variances, case 3

In the through-origin model $`Y_i=\beta x_i+u_i`$, take fixed $`x=(1,2,3)`$, independent mean-zero errors, and $`\mathop{\mathrm{Var}} (u_i)=(1, 1, 16)`$. Derive the OLS and GLS variances.

Solution. The OLS estimator is $`\hat\beta=\sum_i x_iY_i/\sum_i x_i^2`$. Substituting the model yields

```math
\hat\beta-\beta=\frac{\sum_i x_iu_i}{14}.
```

Conditional mean zero proves unbiasedness. Independence makes the variance of the numerator the sum of $`x_i^2\sigma_i^2`$, so

```math
\mathop{\mathrm{Var}} (\hat\beta\mid X)
=\frac{\sum_i x_i^2\sigma_i^2}{(\sum_i x_i^2)^2}
=\frac{149}{196}\approx 0.76020408.
```

This is the scalar sandwich formula. It cannot generally be replaced by one common variance times $`1/14`$.

If the variances are known, minimize $`\sum_i(Y_i-\beta x_i)^2/\sigma_i^2`$. Differentiating gives

```math
\hat\beta_{\mathrm{GLS}}
=\frac{\sum_i x_iY_i/\sigma_i^2}{\sum_i x_i^2/\sigma_i^2},
\qquad
\hat\beta_{\mathrm{GLS}}-\beta
=\frac{\sum_i x_iu_i/\sigma_i^2}{\sum_i x_i^2/\sigma_i^2}.
```

The variance of the numerator is $`\sum_i x_i^2/\sigma_i^2=\frac{89}{16}`$. Therefore $`\mathop{\mathrm{Var}} (\hat\beta_{\mathrm{GLS}}\mid X)=\frac{16}{89}\approx 0.17977528`$.

Cauchy–Schwarz gives $`(\sum_i x_i^2)^2\le(\sum_i x_i^2\sigma_i^2)(\sum_i x_i^2/\sigma_i^2)`$, proving that the GLS variance is no larger. A heteroskedasticity-robust standard error for OLS estimates the first variance using squared residuals; it does not turn the OLS coefficients into GLS coefficients. Estimating the variance weights creates feasible GLS and introduces an additional estimation step.

#### Problem EC-66: precision under unequal error variances, case 4

In the through-origin model $`Y_i=\beta x_i+u_i`$, take fixed $`x=(1,2,3)`$, independent mean-zero errors, and $`\mathop{\mathrm{Var}} (u_i)=(4, 2, 3)`$. Derive the OLS and GLS variances.

Solution. The OLS estimator is $`\hat\beta=\sum_i x_iY_i/\sum_i x_i^2`$. Substituting the model yields

```math
\hat\beta-\beta=\frac{\sum_i x_iu_i}{14}.
```

Conditional mean zero proves unbiasedness. Independence makes the variance of the numerator the sum of $`x_i^2\sigma_i^2`$, so

```math
\mathop{\mathrm{Var}} (\hat\beta\mid X)
=\frac{\sum_i x_i^2\sigma_i^2}{(\sum_i x_i^2)^2}
=\frac{39}{196}\approx 0.19897959.
```

This is the scalar sandwich formula. It cannot generally be replaced by one common variance times $`1/14`$.

If the variances are known, minimize $`\sum_i(Y_i-\beta x_i)^2/\sigma_i^2`$. Differentiating gives

```math
\hat\beta_{\mathrm{GLS}}
=\frac{\sum_i x_iY_i/\sigma_i^2}{\sum_i x_i^2/\sigma_i^2},
\qquad
\hat\beta_{\mathrm{GLS}}-\beta
=\frac{\sum_i x_iu_i/\sigma_i^2}{\sum_i x_i^2/\sigma_i^2}.
```

The variance of the numerator is $`\sum_i x_i^2/\sigma_i^2=\frac{21}{4}`$. Therefore $`\mathop{\mathrm{Var}} (\hat\beta_{\mathrm{GLS}}\mid X)=\frac{4}{21}\approx 0.19047619`$.

Cauchy–Schwarz gives $`(\sum_i x_i^2)^2\le(\sum_i x_i^2\sigma_i^2)(\sum_i x_i^2/\sigma_i^2)`$, proving that the GLS variance is no larger. A heteroskedasticity-robust standard error for OLS estimates the first variance using squared residuals; it does not turn the OLS coefficients into GLS coefficients. Estimating the variance weights creates feasible GLS and introduces an additional estimation step.

### 9.4 Instrument relevance and exclusion

An instrument must shift the endogenous regressor and be excluded from the outcome error. A large first-stage effect addresses relevance but cannot repair an excluded-variable violation. The following moment calculation displays these two requirements separately.

#### Problem EC-67: IV estimand with a direct instrument effect, case 1

Centered variables satisfy $`D=\pi Z+v`$ and $`Y=\beta D+\delta Z+u`$, with $`E[Z^2]=1`$ and $`E[Zv]=E[Zu]=0`$. Take $`\beta=2`$, $`\pi=\frac{1}{2}`$, $`\delta=0`$. Find the population IV ratio using $`Z`$.

Solution. The first-stage covariance is

```math
E[ZD]=\pi E[Z^2]+E[Zv]=\pi=\frac{1}{2}\ne0.
```

The instrument is relevant in the population. The reduced-form covariance is

```math
E[ZY]=\beta E[ZD]+\delta E[Z^2]+E[Zu]
=\beta\pi+\delta=1.
```

Dividing the reduced form by the first stage gives

```math
\beta_{\mathrm{IV}}
=\frac{E[ZY]}{E[ZD]}
=\beta+\frac\delta\pi=2.
```

If the intended structural equation omits $`Z`$, its composite error is $`e=\delta Z+u`$. Then $`E[Ze]=\delta`$, so exclusion and orthogonality hold only when this direct effect is zero under the supplied assumptions.

The deviation from the target coefficient is $`\delta/\pi=0`$. A small first stage can amplify even a small exclusion violation. This is a population identification problem, separate from the finite-sample distribution of the IV estimator.

Multiplying the instrument by a nonzero constant multiplies both numerator and denominator by the same constant and leaves the ratio unchanged. It also leaves the substantive relevance and exclusion conditions unchanged. Rescaling cannot strengthen an instrument.

In a sample, the ratio of sample covariances converges to this population value under an appropriate law of large numbers when $`\pi`$ is fixed away from zero. If relevance shrinks with sample size, ordinary strong-instrument normal approximations may fail. A nonzero reported first-stage coefficient alone is not enough to justify those approximations.

#### Problem EC-68: IV estimand with a direct instrument effect, case 2

Centered variables satisfy $`D=\pi Z+v`$ and $`Y=\beta D+\delta Z+u`$, with $`E[Z^2]=1`$ and $`E[Zv]=E[Zu]=0`$. Take $`\beta=2`$, $`\pi=\frac{1}{2}`$, $`\delta=\frac{1}{10}`$. Find the population IV ratio using $`Z`$.

Solution. The first-stage covariance is

```math
E[ZD]=\pi E[Z^2]+E[Zv]=\pi=\frac{1}{2}\ne0.
```

The instrument is relevant in the population. The reduced-form covariance is

```math
E[ZY]=\beta E[ZD]+\delta E[Z^2]+E[Zu]
=\beta\pi+\delta=\frac{11}{10}.
```

Dividing the reduced form by the first stage gives

```math
\beta_{\mathrm{IV}}
=\frac{E[ZY]}{E[ZD]}
=\beta+\frac\delta\pi=\frac{11}{5}.
```

If the intended structural equation omits $`Z`$, its composite error is $`e=\delta Z+u`$. Then $`E[Ze]=\delta`$, so exclusion and orthogonality hold only when this direct effect is zero under the supplied assumptions.

The deviation from the target coefficient is $`\delta/\pi=\frac{1}{5}`$. A small first stage can amplify even a small exclusion violation. This is a population identification problem, separate from the finite-sample distribution of the IV estimator.

Multiplying the instrument by a nonzero constant multiplies both numerator and denominator by the same constant and leaves the ratio unchanged. It also leaves the substantive relevance and exclusion conditions unchanged. Rescaling cannot strengthen an instrument.

In a sample, the ratio of sample covariances converges to this population value under an appropriate law of large numbers when $`\pi`$ is fixed away from zero. If relevance shrinks with sample size, ordinary strong-instrument normal approximations may fail. A nonzero reported first-stage coefficient alone is not enough to justify those approximations.

#### Problem EC-69: IV estimand with a direct instrument effect, case 3

Centered variables satisfy $`D=\pi Z+v`$ and $`Y=\beta D+\delta Z+u`$, with $`E[Z^2]=1`$ and $`E[Zv]=E[Zu]=0`$. Take $`\beta=1`$, $`\pi=\frac{1}{20}`$, $`\delta=\frac{1}{10}`$. Find the population IV ratio using $`Z`$.

Solution. The first-stage covariance is

```math
E[ZD]=\pi E[Z^2]+E[Zv]=\pi=\frac{1}{20}\ne0.
```

The instrument is relevant in the population. The reduced-form covariance is

```math
E[ZY]=\beta E[ZD]+\delta E[Z^2]+E[Zu]
=\beta\pi+\delta=\frac{3}{20}.
```

Dividing the reduced form by the first stage gives

```math
\beta_{\mathrm{IV}}
=\frac{E[ZY]}{E[ZD]}
=\beta+\frac\delta\pi=3.
```

If the intended structural equation omits $`Z`$, its composite error is $`e=\delta Z+u`$. Then $`E[Ze]=\delta`$, so exclusion and orthogonality hold only when this direct effect is zero under the supplied assumptions.

The deviation from the target coefficient is $`\delta/\pi=2`$. A small first stage can amplify even a small exclusion violation. This is a population identification problem, separate from the finite-sample distribution of the IV estimator.

Multiplying the instrument by a nonzero constant multiplies both numerator and denominator by the same constant and leaves the ratio unchanged. It also leaves the substantive relevance and exclusion conditions unchanged. Rescaling cannot strengthen an instrument.

In a sample, the ratio of sample covariances converges to this population value under an appropriate law of large numbers when $`\pi`$ is fixed away from zero. If relevance shrinks with sample size, ordinary strong-instrument normal approximations may fail. A nonzero reported first-stage coefficient alone is not enough to justify those approximations.

#### Problem EC-70: IV estimand with a direct instrument effect, case 4

Centered variables satisfy $`D=\pi Z+v`$ and $`Y=\beta D+\delta Z+u`$, with $`E[Z^2]=1`$ and $`E[Zv]=E[Zu]=0`$. Take $`\beta=-1`$, $`\pi=-2`$, $`\delta=\frac{1}{5}`$. Find the population IV ratio using $`Z`$.

Solution. The first-stage covariance is

```math
E[ZD]=\pi E[Z^2]+E[Zv]=\pi=-2\ne0.
```

The instrument is relevant in the population. The reduced-form covariance is

```math
E[ZY]=\beta E[ZD]+\delta E[Z^2]+E[Zu]
=\beta\pi+\delta=\frac{11}{5}.
```

Dividing the reduced form by the first stage gives

```math
\beta_{\mathrm{IV}}
=\frac{E[ZY]}{E[ZD]}
=\beta+\frac\delta\pi=-\frac{11}{10}.
```

If the intended structural equation omits $`Z`$, its composite error is $`e=\delta Z+u`$. Then $`E[Ze]=\delta`$, so exclusion and orthogonality hold only when this direct effect is zero under the supplied assumptions.

The deviation from the target coefficient is $`\delta/\pi=-\frac{1}{10}`$. A small first stage can amplify even a small exclusion violation. This is a population identification problem, separate from the finite-sample distribution of the IV estimator.

Multiplying the instrument by a nonzero constant multiplies both numerator and denominator by the same constant and leaves the ratio unchanged. It also leaves the substantive relevance and exclusion conditions unchanged. Rescaling cannot strengthen an instrument.

In a sample, the ratio of sample covariances converges to this population value under an appropriate law of large numbers when $`\pi`$ is fixed away from zero. If relevance shrinks with sample size, ordinary strong-instrument normal approximations may fail. A nonzero reported first-stage coefficient alone is not enough to justify those approximations.

### 9.5 Efficiently combining two moment estimates

Two unbiased estimates of the same scalar parameter can contain correlated noise. An efficient linear combination uses their covariance, not only their individual variances. With strong positive correlation, an optimal weight may be negative because one estimate helps remove a shared noise component.

#### Problem EC-71: two-moment efficient weighting, case 1

Suppose $`\sqrt n((\hat\mu_1,\hat\mu_2)^\top-\mu\mathbf1)`$ has asymptotic covariance

```math
\Omega=\begin{pmatrix}1&0\\[0pt] 0&4\end{pmatrix}
```

, with $`n=100`$. Find the minimum-variance linear combination whose weights sum to one.

Solution. Write the combination as $`\hat\mu(w)=w\hat\mu_1+(1-w)\hat\mu_2`$. Its asymptotic variance before division by $`n`$ is

```math
V(w)=w^2v_1+(1-w)^2v_2+2w(1-w)c.
```

The coefficient of $`w^2`$ is $`v_1+v_2-2c=5\gt 0`$. Differentiating and setting the derivative to zero gives

```math
w^*=\frac{v_2-c}{v_1+v_2-2c}=\frac{4}{5},\qquad
1-w^*=\frac{1}{5}.
```

The second derivative is $`2(v_1+v_2-2c)\gt 0`$, so this solution is the unique variance minimum over all real weights summing to one. Positivity of individual weights was not imposed.

Substitution simplifies the minimum asymptotic variance to

```math
V(w^*)=\frac{v_1v_2-c^2}{v_1+v_2-2c}
=\frac{4}{5}.
```

The estimated sampling variance for sample size $`n`$ is consequently $`V(w^*)/n=\frac{1}{125}`$, with standard error approximately $`0.089443`$.

The covariance matrix is positive definite because its determinant is $`4\gt 0`$ and its first diagonal entry is positive. Thus the efficient variance is strictly positive. The vector formula is $`w=\Omega^{-1}\mathbf1/(\mathbf1^\top\Omega^{-1}\mathbf1)`$, and the minimized variance is $`(\mathbf1^\top\Omega^{-1}\mathbf1)^{-1}`$.

These are the same inverse-covariance weights that arise in efficient GMM when the two moment estimates identify one common parameter. If either moment is systematically biased for that parameter, optimal variance weighting does not remove the bias. Efficiency is a precision statement conditional on valid moment restrictions.

#### Problem EC-72: two-moment efficient weighting, case 2

Suppose $`\sqrt n((\hat\mu_1,\hat\mu_2)^\top-\mu\mathbf1)`$ has asymptotic covariance

```math
\Omega=\begin{pmatrix}2&1\\[0pt] 1&3\end{pmatrix}
```

, with $`n=200`$. Find the minimum-variance linear combination whose weights sum to one.

Solution. Write the combination as $`\hat\mu(w)=w\hat\mu_1+(1-w)\hat\mu_2`$. Its asymptotic variance before division by $`n`$ is

```math
V(w)=w^2v_1+(1-w)^2v_2+2w(1-w)c.
```

The coefficient of $`w^2`$ is $`v_1+v_2-2c=3\gt 0`$. Differentiating and setting the derivative to zero gives

```math
w^*=\frac{v_2-c}{v_1+v_2-2c}=\frac{2}{3},\qquad
1-w^*=\frac{1}{3}.
```

The second derivative is $`2(v_1+v_2-2c)\gt 0`$, so this solution is the unique variance minimum over all real weights summing to one. Positivity of individual weights was not imposed.

Substitution simplifies the minimum asymptotic variance to

```math
V(w^*)=\frac{v_1v_2-c^2}{v_1+v_2-2c}
=\frac{5}{3}.
```

The estimated sampling variance for sample size $`n`$ is consequently $`V(w^*)/n=\frac{1}{120}`$, with standard error approximately $`0.091287`$.

The covariance matrix is positive definite because its determinant is $`5\gt 0`$ and its first diagonal entry is positive. Thus the efficient variance is strictly positive. The vector formula is $`w=\Omega^{-1}\mathbf1/(\mathbf1^\top\Omega^{-1}\mathbf1)`$, and the minimized variance is $`(\mathbf1^\top\Omega^{-1}\mathbf1)^{-1}`$.

These are the same inverse-covariance weights that arise in efficient GMM when the two moment estimates identify one common parameter. If either moment is systematically biased for that parameter, optimal variance weighting does not remove the bias. Efficiency is a precision statement conditional on valid moment restrictions.

#### Problem EC-73: two-moment efficient weighting, case 3

Suppose $`\sqrt n((\hat\mu_1,\hat\mu_2)^\top-\mu\mathbf1)`$ has asymptotic covariance

```math
\Omega=\begin{pmatrix}1&\frac{3}{2}\\[0pt] \frac{3}{2}&4\end{pmatrix}
```

, with $`n=100`$. Find the minimum-variance linear combination whose weights sum to one.

Solution. Write the combination as $`\hat\mu(w)=w\hat\mu_1+(1-w)\hat\mu_2`$. Its asymptotic variance before division by $`n`$ is

```math
V(w)=w^2v_1+(1-w)^2v_2+2w(1-w)c.
```

The coefficient of $`w^2`$ is $`v_1+v_2-2c=2\gt 0`$. Differentiating and setting the derivative to zero gives

```math
w^*=\frac{v_2-c}{v_1+v_2-2c}=\frac{5}{4},\qquad
1-w^*=-\frac{1}{4}.
```

The second derivative is $`2(v_1+v_2-2c)\gt 0`$, so this solution is the unique variance minimum over all real weights summing to one. Positivity of individual weights was not imposed.

Substitution simplifies the minimum asymptotic variance to

```math
V(w^*)=\frac{v_1v_2-c^2}{v_1+v_2-2c}
=\frac{7}{8}.
```

The estimated sampling variance for sample size $`n`$ is consequently $`V(w^*)/n=\frac{7}{800}`$, with standard error approximately $`0.093541`$.

The covariance matrix is positive definite because its determinant is $`\frac{7}{4}\gt 0`$ and its first diagonal entry is positive. Thus the efficient variance is strictly positive. The vector formula is $`w=\Omega^{-1}\mathbf1/(\mathbf1^\top\Omega^{-1}\mathbf1)`$, and the minimized variance is $`(\mathbf1^\top\Omega^{-1}\mathbf1)^{-1}`$.

These are the same inverse-covariance weights that arise in efficient GMM when the two moment estimates identify one common parameter. If either moment is systematically biased for that parameter, optimal variance weighting does not remove the bias. Efficiency is a precision statement conditional on valid moment restrictions.

#### Problem EC-74: two-moment efficient weighting, case 4

Suppose $`\sqrt n((\hat\mu_1,\hat\mu_2)^\top-\mu\mathbf1)`$ has asymptotic covariance

```math
\Omega=\begin{pmatrix}4&\frac{1}{2}\\[0pt] \frac{1}{2}&1\end{pmatrix}
```

, with $`n=50`$. Find the minimum-variance linear combination whose weights sum to one.

Solution. Write the combination as $`\hat\mu(w)=w\hat\mu_1+(1-w)\hat\mu_2`$. Its asymptotic variance before division by $`n`$ is

```math
V(w)=w^2v_1+(1-w)^2v_2+2w(1-w)c.
```

The coefficient of $`w^2`$ is $`v_1+v_2-2c=4\gt 0`$. Differentiating and setting the derivative to zero gives

```math
w^*=\frac{v_2-c}{v_1+v_2-2c}=\frac{1}{8},\qquad
1-w^*=\frac{7}{8}.
```

The second derivative is $`2(v_1+v_2-2c)\gt 0`$, so this solution is the unique variance minimum over all real weights summing to one. Positivity of individual weights was not imposed.

Substitution simplifies the minimum asymptotic variance to

```math
V(w^*)=\frac{v_1v_2-c^2}{v_1+v_2-2c}
=\frac{15}{16}.
```

The estimated sampling variance for sample size $`n`$ is consequently $`V(w^*)/n=\frac{3}{160}`$, with standard error approximately $`0.136931`$.

The covariance matrix is positive definite because its determinant is $`\frac{15}{4}\gt 0`$ and its first diagonal entry is positive. Thus the efficient variance is strictly positive. The vector formula is $`w=\Omega^{-1}\mathbf1/(\mathbf1^\top\Omega^{-1}\mathbf1)`$, and the minimized variance is $`(\mathbf1^\top\Omega^{-1}\mathbf1)^{-1}`$.

These are the same inverse-covariance weights that arise in efficient GMM when the two moment estimates identify one common parameter. If either moment is systematically biased for that parameter, optimal variance weighting does not remove the bias. Efficiency is a precision statement conditional on valid moment restrictions.

### 9.6 Difference in differences and a trend violation

The difference-in-differences calculation removes group levels and a common time change. Its causal interpretation requires a statement about the treated group's untreated change. A post-treatment outcome table alone cannot verify that counterfactual trend.

#### Problem EC-75: two-group two-period treatment comparison, case 1

Treated-group means before and after are $`(10,16)`$; control-group means are $`(8,11)`$. Compute DID. Then suppose the treated group's untreated trend would exceed the control trend by $`0`$; derive the treatment effect consistent with that assumption.

Solution. The treated change is $`16-10=6`$ and the control change is $`11-8=3`$. Therefore the observed difference in differences is

```math
\widehat{\mathrm{DID}}=(16-10)-(11-8)=3.
```

Under parallel untreated trends, the treated group's counterfactual post-treatment mean would be its baseline plus the control change, namely $`13`$. Subtracting that counterfactual from the observed post-treatment mean gives the DID effect $`3`$.

Now replace parallel trends by the supplied deviation. The untreated change for the treated group is $`3+(0)=3`$. Its counterfactual post-treatment level becomes $`13`$. The corresponding treatment effect is

```math
\tau=16-(13)=3
=\mathrm{DID}-(0).
```

This is an exact bias decomposition: the observed DID equals the treatment effect plus the difference in untreated trends. The observed cell means cannot separately determine both quantities.

The same DID number is the interaction coefficient in a saturated regression on an intercept, a treated-group indicator, a post-period indicator, and their interaction. The intercept is $`8`$, the group coefficient is $`2`$, the time coefficient is $`3`$, and the interaction is $`3`$. These four coefficients reproduce all four cell means.

That regression identity does not prove parallel trends. In applications, no anticipation, a stable comparison group, and an appropriate sampling and clustering scheme also need attention. With only these four means, there is no information about the sampling variance or the number of independent treatment-assignment units.

#### Problem EC-76: two-group two-period treatment comparison, case 2

Treated-group means before and after are $`(5,9)`$; control-group means are $`(7,8)`$. Compute DID. Then suppose the treated group's untreated trend would exceed the control trend by $`1`$; derive the treatment effect consistent with that assumption.

Solution. The treated change is $`9-5=4`$ and the control change is $`8-7=1`$. Therefore the observed difference in differences is

```math
\widehat{\mathrm{DID}}=(9-5)-(8-7)=3.
```

Under parallel untreated trends, the treated group's counterfactual post-treatment mean would be its baseline plus the control change, namely $`6`$. Subtracting that counterfactual from the observed post-treatment mean gives the DID effect $`3`$.

Now replace parallel trends by the supplied deviation. The untreated change for the treated group is $`1+(1)=2`$. Its counterfactual post-treatment level becomes $`7`$. The corresponding treatment effect is

```math
\tau=9-(7)=2
=\mathrm{DID}-(1).
```

This is an exact bias decomposition: the observed DID equals the treatment effect plus the difference in untreated trends. The observed cell means cannot separately determine both quantities.

The same DID number is the interaction coefficient in a saturated regression on an intercept, a treated-group indicator, a post-period indicator, and their interaction. The intercept is $`7`$, the group coefficient is $`-2`$, the time coefficient is $`1`$, and the interaction is $`3`$. These four coefficients reproduce all four cell means.

That regression identity does not prove parallel trends. In applications, no anticipation, a stable comparison group, and an appropriate sampling and clustering scheme also need attention. With only these four means, there is no information about the sampling variance or the number of independent treatment-assignment units.

#### Problem EC-77: two-group two-period treatment comparison, case 3

Treated-group means before and after are $`(20,18)`$; control-group means are $`(12,11)`$. Compute DID. Then suppose the treated group's untreated trend would exceed the control trend by $`-2`$; derive the treatment effect consistent with that assumption.

Solution. The treated change is $`18-20=-2`$ and the control change is $`11-12=-1`$. Therefore the observed difference in differences is

```math
\widehat{\mathrm{DID}}=(18-20)-(11-12)=-1.
```

Under parallel untreated trends, the treated group's counterfactual post-treatment mean would be its baseline plus the control change, namely $`19`$. Subtracting that counterfactual from the observed post-treatment mean gives the DID effect $`-1`$.

Now replace parallel trends by the supplied deviation. The untreated change for the treated group is $`-1+(-2)=-3`$. Its counterfactual post-treatment level becomes $`17`$. The corresponding treatment effect is

```math
\tau=18-(17)=1
=\mathrm{DID}-(-2).
```

This is an exact bias decomposition: the observed DID equals the treatment effect plus the difference in untreated trends. The observed cell means cannot separately determine both quantities.

The same DID number is the interaction coefficient in a saturated regression on an intercept, a treated-group indicator, a post-period indicator, and their interaction. The intercept is $`12`$, the group coefficient is $`8`$, the time coefficient is $`-1`$, and the interaction is $`-1`$. These four coefficients reproduce all four cell means.

That regression identity does not prove parallel trends. In applications, no anticipation, a stable comparison group, and an appropriate sampling and clustering scheme also need attention. With only these four means, there is no information about the sampling variance or the number of independent treatment-assignment units.

#### Problem EC-78: two-group two-period treatment comparison, case 4

Treated-group means before and after are $`(3,10)`$; control-group means are $`(4,6)`$. Compute DID. Then suppose the treated group's untreated trend would exceed the control trend by $`2`$; derive the treatment effect consistent with that assumption.

Solution. The treated change is $`10-3=7`$ and the control change is $`6-4=2`$. Therefore the observed difference in differences is

```math
\widehat{\mathrm{DID}}=(10-3)-(6-4)=5.
```

Under parallel untreated trends, the treated group's counterfactual post-treatment mean would be its baseline plus the control change, namely $`5`$. Subtracting that counterfactual from the observed post-treatment mean gives the DID effect $`5`$.

Now replace parallel trends by the supplied deviation. The untreated change for the treated group is $`2+(2)=4`$. Its counterfactual post-treatment level becomes $`7`$. The corresponding treatment effect is

```math
\tau=10-(7)=3
=\mathrm{DID}-(2).
```

This is an exact bias decomposition: the observed DID equals the treatment effect plus the difference in untreated trends. The observed cell means cannot separately determine both quantities.

The same DID number is the interaction coefficient in a saturated regression on an intercept, a treated-group indicator, a post-period indicator, and their interaction. The intercept is $`4`$, the group coefficient is $`-1`$, the time coefficient is $`2`$, and the interaction is $`5`$. These four coefficients reproduce all four cell means.

That regression identity does not prove parallel trends. In applications, no anticipation, a stable comparison group, and an appropriate sampling and clustering scheme also need attention. With only these four means, there is no information about the sampling variance or the number of independent treatment-assignment units.

### 9.7 AR forecasts and forecast-error uncertainty

For a stable AR(1), iterating the recursion produces a deterministic conditional mean and a sum of future innovations. The observed current value affects the conditional mean, while the innovation variance controls forecast uncertainty. A nonzero unconditional mean must be centered correctly.

#### Problem EC-79: centered AR forecast, case 1

Let $`X_{t+1}-\mu=\rho(X_t-\mu)+\varepsilon_{t+1}`$ with independent mean-zero innovations of variance $`1`$. Take $`\rho=\frac{1}{2}`$, $`\mu=2`$, $`X_t=5`$. Find the $`h=3`$ forecast and its mean squared error.

Solution. One substitution gives $`X_{t+2}-\mu=\rho^2(X_t-\mu)+\rho\varepsilon_{t+1}+\varepsilon_{t+2}`$. Repeating the same step yields

```math
X_{t+h}-\mu
=\rho^h(X_t-\mu)+\sum_{j=1}^h\rho^{h-j}\varepsilon_{t+j}.
```

Future innovations have conditional mean zero, so

```math
\hat X_{t+h\mid t}=\mu+\rho^h(X_t-\mu)=\frac{19}{8}.
```

Subtracting this forecast leaves exactly the innovation sum. Its cross terms have zero expectation by independence, giving

```math
\mathop{\mathrm{MSE}} _h
=\sigma_\varepsilon^2\sum_{k=0}^{h-1}\rho^{2k}
=\sigma_\varepsilon^2\frac{1-\rho^{2h}}{1-\rho^2}
=\frac{21}{16}\approx 1.312500.
```

The geometric-series formula is valid because $`\vert \rho\vert \lt 1`$ here. The stationary variance is $`\sigma_\varepsilon^2/(1-\rho^2)=\frac{4}{3}`$, and the forecast error variance approaches it as the horizon grows. The conditional mean approaches $`\mu`$. For negative $`\rho`$, the forecast deviation alternates sign across horizons even though the forecast-error variance still increases.

If innovations are Gaussian, a conditional 95% prediction interval is the forecast plus or minus $`1.96\sqrt{\mathop{\mathrm{MSE}} _h}`$ when parameters are known. Without Gaussian innovations, the conditional mean and MSE remain correct, but that normal interval need not have exact coverage. Estimating parameters adds uncertainty not included in this known-parameter calculation.

#### Problem EC-80: centered AR forecast, case 2

Let $`X_{t+1}-\mu=\rho(X_t-\mu)+\varepsilon_{t+1}`$ with independent mean-zero innovations of variance $`2`$. Take $`\rho=\frac{4}{5}`$, $`\mu=0`$, $`X_t=-2`$. Find the $`h=4`$ forecast and its mean squared error.

Solution. One substitution gives $`X_{t+2}-\mu=\rho^2(X_t-\mu)+\rho\varepsilon_{t+1}+\varepsilon_{t+2}`$. Repeating the same step yields

```math
X_{t+h}-\mu
=\rho^h(X_t-\mu)+\sum_{j=1}^h\rho^{h-j}\varepsilon_{t+j}.
```

Future innovations have conditional mean zero, so

```math
\hat X_{t+h\mid t}=\mu+\rho^h(X_t-\mu)=-\frac{512}{625}.
```

Subtracting this forecast leaves exactly the innovation sum. Its cross terms have zero expectation by independence, giving

```math
\mathop{\mathrm{MSE}} _h
=\sigma_\varepsilon^2\sum_{k=0}^{h-1}\rho^{2k}
=\sigma_\varepsilon^2\frac{1-\rho^{2h}}{1-\rho^2}
=\frac{72242}{15625}\approx 4.623488.
```

The geometric-series formula is valid because $`\vert \rho\vert \lt 1`$ here. The stationary variance is $`\sigma_\varepsilon^2/(1-\rho^2)=\frac{50}{9}`$, and the forecast error variance approaches it as the horizon grows. The conditional mean approaches $`\mu`$. For negative $`\rho`$, the forecast deviation alternates sign across horizons even though the forecast-error variance still increases.

If innovations are Gaussian, a conditional 95% prediction interval is the forecast plus or minus $`1.96\sqrt{\mathop{\mathrm{MSE}} _h}`$ when parameters are known. Without Gaussian innovations, the conditional mean and MSE remain correct, but that normal interval need not have exact coverage. Estimating parameters adds uncertainty not included in this known-parameter calculation.

#### Problem EC-81: centered AR forecast, case 3

Let $`X_{t+1}-\mu=\rho(X_t-\mu)+\varepsilon_{t+1}`$ with independent mean-zero innovations of variance $`1`$. Take $`\rho=-\frac{1}{2}`$, $`\mu=1`$, $`X_t=4`$. Find the $`h=3`$ forecast and its mean squared error.

Solution. One substitution gives $`X_{t+2}-\mu=\rho^2(X_t-\mu)+\rho\varepsilon_{t+1}+\varepsilon_{t+2}`$. Repeating the same step yields

```math
X_{t+h}-\mu
=\rho^h(X_t-\mu)+\sum_{j=1}^h\rho^{h-j}\varepsilon_{t+j}.
```

Future innovations have conditional mean zero, so

```math
\hat X_{t+h\mid t}=\mu+\rho^h(X_t-\mu)=\frac{5}{8}.
```

Subtracting this forecast leaves exactly the innovation sum. Its cross terms have zero expectation by independence, giving

```math
\mathop{\mathrm{MSE}} _h
=\sigma_\varepsilon^2\sum_{k=0}^{h-1}\rho^{2k}
=\sigma_\varepsilon^2\frac{1-\rho^{2h}}{1-\rho^2}
=\frac{21}{16}\approx 1.312500.
```

The geometric-series formula is valid because $`\vert \rho\vert \lt 1`$ here. The stationary variance is $`\sigma_\varepsilon^2/(1-\rho^2)=\frac{4}{3}`$, and the forecast error variance approaches it as the horizon grows. The conditional mean approaches $`\mu`$. For negative $`\rho`$, the forecast deviation alternates sign across horizons even though the forecast-error variance still increases.

If innovations are Gaussian, a conditional 95% prediction interval is the forecast plus or minus $`1.96\sqrt{\mathop{\mathrm{MSE}} _h}`$ when parameters are known. Without Gaussian innovations, the conditional mean and MSE remain correct, but that normal interval need not have exact coverage. Estimating parameters adds uncertainty not included in this known-parameter calculation.

#### Problem EC-82: centered AR forecast, case 4

Let $`X_{t+1}-\mu=\rho(X_t-\mu)+\varepsilon_{t+1}`$ with independent mean-zero innovations of variance $`\frac{1}{4}`$. Take $`\rho=\frac{9}{10}`$, $`\mu=3`$, $`X_t=1`$. Find the $`h=5`$ forecast and its mean squared error.

Solution. One substitution gives $`X_{t+2}-\mu=\rho^2(X_t-\mu)+\rho\varepsilon_{t+1}+\varepsilon_{t+2}`$. Repeating the same step yields

```math
X_{t+h}-\mu
=\rho^h(X_t-\mu)+\sum_{j=1}^h\rho^{h-j}\varepsilon_{t+j}.
```

Future innovations have conditional mean zero, so

```math
\hat X_{t+h\mid t}=\mu+\rho^h(X_t-\mu)=\frac{90951}{50000}.
```

Subtracting this forecast leaves exactly the innovation sum. Its cross terms have zero expectation by independence, giving

```math
\mathop{\mathrm{MSE}} _h
=\sigma_\varepsilon^2\sum_{k=0}^{h-1}\rho^{2k}
=\sigma_\varepsilon^2\frac{1-\rho^{2h}}{1-\rho^2}
=\frac{342800821}{400000000}\approx 0.857002.
```

The geometric-series formula is valid because $`\vert \rho\vert \lt 1`$ here. The stationary variance is $`\sigma_\varepsilon^2/(1-\rho^2)=\frac{25}{19}`$, and the forecast error variance approaches it as the horizon grows. The conditional mean approaches $`\mu`$. For negative $`\rho`$, the forecast deviation alternates sign across horizons even though the forecast-error variance still increases.

If innovations are Gaussian, a conditional 95% prediction interval is the forecast plus or minus $`1.96\sqrt{\mathop{\mathrm{MSE}} _h}`$ when parameters are known. Without Gaussian innovations, the conditional mean and MSE remain correct, but that normal interval need not have exact coverage. Estimating parameters adds uncertainty not included in this known-parameter calculation.

### 9.8 A triangular VAR with cross-variable propagation

In an upper-triangular VAR, the second variable can affect future values of the first without feedback in the opposite direction. A multi-period forecast error contains propagated shocks from every intermediate date. It is not generally just the one-step innovation covariance multiplied by the horizon.

#### Problem EC-83: two-variable VAR propagation, case 1

Let $`X_{t+1}=AX_t+\varepsilon_{t+1}`$,

```math
A=\begin{pmatrix}\frac{1}{2}&1\\[0pt] 0&\frac{1}{4}\end{pmatrix}
```

, with serially independent mean-zero innovations having covariance $`I_2`$. From $`X_t=(1,2)^\top`$, find the $`h=2`$ conditional mean and forecast-error covariance.

Solution. The eigenvalues are the diagonal entries $`a=\frac{1}{2}`$ and $`d=\frac{1}{4}`$, both inside the unit circle. Thus this VAR has a causal stationary solution under finite innovation variances. The off-diagonal entry affects propagation even though it does not change these eigenvalues.

For an upper-triangular matrix with $`a\ne d`$,

```math
A^h=
\begin{pmatrix}
a^h&b\sum_{j=0}^{h-1}a^{h-1-j}d^j\\[0pt]
0&d^h
\end{pmatrix},
\qquad
\sum_{j=0}^{h-1}a^{h-1-j}d^j=\frac{a^h-d^h}{a-d}.
```

At the requested horizon, direct multiplication gives

```math
A^h=\begin{pmatrix}\frac{1}{4}&\frac{3}{4}\\[0pt] 0&\frac{1}{16}\end{pmatrix}
```

. The conditional mean is therefore

```math
A^hX_t=\begin{pmatrix}\frac{7}{4}\\[0pt] \frac{1}{8}\end{pmatrix}
```

.

Iterating the stochastic recursion gives forecast error $`\sum_{j=1}^h A^{h-j}\varepsilon_{t+j}`$. Innovations at different dates are uncorrelated, so their propagated covariance matrices add:

```math
\Sigma_h=\sum_{k=0}^{h-1}A^k I_2(A^k)^\top
=\begin{pmatrix}\frac{9}{4}&\frac{1}{4}\\[0pt] \frac{1}{4}&\frac{17}{16}\end{pmatrix}.
```

The $`k=0`$ term is $`I_2`$. Each additional term is positive semidefinite, so forecast uncertainty increases in the positive-semidefinite ordering as the horizon expands. Off-diagonal covariance can appear even though the original innovations have diagonal covariance.

The second coordinate's shock affects its own future and, through $`b`$, the first coordinate's future. A shock to the first coordinate does not enter the second equation in this triangular system. These are propagation statements for the specified reduced-form innovations. Giving them a structural economic interpretation requires additional identification assumptions.

#### Problem EC-84: two-variable VAR propagation, case 2

Let $`X_{t+1}=AX_t+\varepsilon_{t+1}`$,

```math
A=\begin{pmatrix}\frac{4}{5}&\frac{1}{5}\\[0pt] 0&\frac{1}{2}\end{pmatrix}
```

, with serially independent mean-zero innovations having covariance $`I_2`$. From $`X_t=(1,2)^\top`$, find the $`h=3`$ conditional mean and forecast-error covariance.

Solution. The eigenvalues are the diagonal entries $`a=\frac{4}{5}`$ and $`d=\frac{1}{2}`$, both inside the unit circle. Thus this VAR has a causal stationary solution under finite innovation variances. The off-diagonal entry affects propagation even though it does not change these eigenvalues.

For an upper-triangular matrix with $`a\ne d`$,

```math
A^h=
\begin{pmatrix}
a^h&b\sum_{j=0}^{h-1}a^{h-1-j}d^j\\[0pt]
0&d^h
\end{pmatrix},
\qquad
\sum_{j=0}^{h-1}a^{h-1-j}d^j=\frac{a^h-d^h}{a-d}.
```

At the requested horizon, direct multiplication gives

```math
A^h=\begin{pmatrix}\frac{64}{125}&\frac{129}{500}\\[0pt] 0&\frac{1}{8}\end{pmatrix}
```

. The conditional mean is therefore

```math
A^hX_t=\begin{pmatrix}\frac{257}{250}\\[0pt] \frac{1}{4}\end{pmatrix}
```

.

Iterating the stochastic recursion gives forecast error $`\sum_{j=1}^h A^{h-j}\varepsilon_{t+j}`$. Innovations at different dates are uncorrelated, so their propagated covariance matrices add:

```math
\Sigma_h=\sum_{k=0}^{h-1}A^k I_2(A^k)^\top
=\begin{pmatrix}\frac{5393}{2500}&\frac{33}{200}\\[0pt] \frac{33}{200}&\frac{21}{16}\end{pmatrix}.
```

The $`k=0`$ term is $`I_2`$. Each additional term is positive semidefinite, so forecast uncertainty increases in the positive-semidefinite ordering as the horizon expands. Off-diagonal covariance can appear even though the original innovations have diagonal covariance.

The second coordinate's shock affects its own future and, through $`b`$, the first coordinate's future. A shock to the first coordinate does not enter the second equation in this triangular system. These are propagation statements for the specified reduced-form innovations. Giving them a structural economic interpretation requires additional identification assumptions.

#### Problem EC-85: two-variable VAR propagation, case 3

Let $`X_{t+1}=AX_t+\varepsilon_{t+1}`$,

```math
A=\begin{pmatrix}-\frac{1}{2}&\frac{1}{2}\\[0pt] 0&\frac{1}{3}\end{pmatrix}
```

, with serially independent mean-zero innovations having covariance $`I_2`$. From $`X_t=(1,2)^\top`$, find the $`h=3`$ conditional mean and forecast-error covariance.

Solution. The eigenvalues are the diagonal entries $`a=-\frac{1}{2}`$ and $`d=\frac{1}{3}`$, both inside the unit circle. Thus this VAR has a causal stationary solution under finite innovation variances. The off-diagonal entry affects propagation even though it does not change these eigenvalues.

For an upper-triangular matrix with $`a\ne d`$,

```math
A^h=
\begin{pmatrix}
a^h&b\sum_{j=0}^{h-1}a^{h-1-j}d^j\\[0pt]
0&d^h
\end{pmatrix},
\qquad
\sum_{j=0}^{h-1}a^{h-1-j}d^j=\frac{a^h-d^h}{a-d}.
```

At the requested horizon, direct multiplication gives

```math
A^h=\begin{pmatrix}-\frac{1}{8}&\frac{7}{72}\\[0pt] 0&\frac{1}{27}\end{pmatrix}
```

. The conditional mean is therefore

```math
A^hX_t=\begin{pmatrix}\frac{5}{72}\\[0pt] \frac{2}{27}\end{pmatrix}
```

.

Iterating the stochastic recursion gives forecast error $`\sum_{j=1}^h A^{h-j}\varepsilon_{t+j}`$. Innovations at different dates are uncorrelated, so their propagated covariance matrices add:

```math
\Sigma_h=\sum_{k=0}^{h-1}A^k I_2(A^k)^\top
=\begin{pmatrix}\frac{113}{72}&\frac{17}{108}\\[0pt] \frac{17}{108}&\frac{91}{81}\end{pmatrix}.
```

The $`k=0`$ term is $`I_2`$. Each additional term is positive semidefinite, so forecast uncertainty increases in the positive-semidefinite ordering as the horizon expands. Off-diagonal covariance can appear even though the original innovations have diagonal covariance.

The second coordinate's shock affects its own future and, through $`b`$, the first coordinate's future. A shock to the first coordinate does not enter the second equation in this triangular system. These are propagation statements for the specified reduced-form innovations. Giving them a structural economic interpretation requires additional identification assumptions.

#### Problem EC-86: two-variable VAR propagation, case 4

Let $`X_{t+1}=AX_t+\varepsilon_{t+1}`$,

```math
A=\begin{pmatrix}\frac{2}{3}&-1\\[0pt] 0&\frac{1}{3}\end{pmatrix}
```

, with serially independent mean-zero innovations having covariance $`I_2`$. From $`X_t=(1,2)^\top`$, find the $`h=2`$ conditional mean and forecast-error covariance.

Solution. The eigenvalues are the diagonal entries $`a=\frac{2}{3}`$ and $`d=\frac{1}{3}`$, both inside the unit circle. Thus this VAR has a causal stationary solution under finite innovation variances. The off-diagonal entry affects propagation even though it does not change these eigenvalues.

For an upper-triangular matrix with $`a\ne d`$,

```math
A^h=
\begin{pmatrix}
a^h&b\sum_{j=0}^{h-1}a^{h-1-j}d^j\\[0pt]
0&d^h
\end{pmatrix},
\qquad
\sum_{j=0}^{h-1}a^{h-1-j}d^j=\frac{a^h-d^h}{a-d}.
```

At the requested horizon, direct multiplication gives

```math
A^h=\begin{pmatrix}\frac{4}{9}&-1\\[0pt] 0&\frac{1}{9}\end{pmatrix}
```

. The conditional mean is therefore

```math
A^hX_t=\begin{pmatrix}-\frac{14}{9}\\[0pt] \frac{2}{9}\end{pmatrix}
```

.

Iterating the stochastic recursion gives forecast error $`\sum_{j=1}^h A^{h-j}\varepsilon_{t+j}`$. Innovations at different dates are uncorrelated, so their propagated covariance matrices add:

```math
\Sigma_h=\sum_{k=0}^{h-1}A^k I_2(A^k)^\top
=\begin{pmatrix}\frac{22}{9}&-\frac{1}{3}\\[0pt] -\frac{1}{3}&\frac{10}{9}\end{pmatrix}.
```

The $`k=0`$ term is $`I_2`$. Each additional term is positive semidefinite, so forecast uncertainty increases in the positive-semidefinite ordering as the horizon expands. Off-diagonal covariance can appear even though the original innovations have diagonal covariance.

The second coordinate's shock affects its own future and, through $`b`$, the first coordinate's future. A shock to the first coordinate does not enter the second equation in this triangular system. These are propagation statements for the specified reduced-form innovations. Giving them a structural economic interpretation requires additional identification assumptions.

<!-- further-distinct-problems-EC -->

### 9.9 Further estimation and dependence problems

These additional problems are constructed for these notes.

#### Problem EC-87: a robust slope variance by hand

For an intercept regression with $`x=(0,1,2,3)`$ and $`y=(1,2,2,5)`$, compute the OLS slope and its homoskedastic, HC0, HC1, and HC2 variance estimates. Compare the numbers without assuming an exact test distribution.

Solution. The centered regressor sum of squares is five. The OLS intercept is $`0.7`$ and slope $`1.2`$, giving residuals $`(0.3,0.1,-1.1,0.7)`$. Their sum and regressor inner product are zero.

The residual sum of squares is $`1.8`$. With two residual degrees of freedom, $`s^2=0.9`$ and the homoskedastic slope variance estimate is $`s^2/5=0.18`$.

The slope influence weights in the linear OLS formula are $`(x_i-\bar x)/5=(-0.3,-0.1,0.1,0.3)`$. HC0 therefore gives

```math
\widehat{\mathop{\mathrm{Var}} }_{\mathrm{HC0}}(\hat\beta_1)
=0.09(0.09)+0.01(0.01)+0.01(1.21)+0.09(0.49)
=0.0644.
```

HC1 multiplies by $`n/(n-p)=2`$, producing $`0.1288`$. The leverages are $`(0.7,0.3,0.3,0.7)`$. HC2 divides each squared residual contribution by one minus its leverage:

```math
\frac{0.0081}{0.3}+\frac{0.0001}{0.7}
+\frac{0.0121}{0.7}+\frac{0.0441}{0.3}
=\frac{67}{350}\approx0.191429.
```

Robust variance estimates need not exceed the homoskedastic estimate. They estimate a different covariance expression and include different finite-sample corrections. With only four observations, their asymptotic justification does not provide an exact small-sample t law. The calculation displays the estimators and their leverage adjustment, rather than claiming that the largest or smallest number is automatically the correct uncertainty.

#### Problem EC-88: a full-rank IV system with a weak direction

Let

```math
E[ZD^\top]=G=\begin{pmatrix}1&1\\[0pt] 1&1+\epsilon\end{pmatrix}
```

 with $`\epsilon\gt 0`$ and $`E[ZY]=G(1,2)^\top`$. Solve the two identifying moment equations. Then perturb only the second entry of $`E[ZY]`$ by $`\delta`$ and analyze sensitivity.

Solution. The determinant is $`\epsilon`$, so the instrument-regressor matrix has full rank for every stated positive value. Its inverse is

```math
G^{-1}=\frac1\epsilon
\begin{pmatrix}1+\epsilon&-1\\[0pt] -1&1\end{pmatrix}.
```

The unperturbed outcome cross-moment is $`(3,3+2\epsilon)^\top`$, and multiplication recovers $`\beta=(1,2)^\top`$.

A perturbation $`(0,\delta)^\top`$ changes the coefficient by

```math
\Delta\beta=G^{-1}\binom0\delta
=\binom{-\delta/\epsilon}{\delta/\epsilon}.
```

Thus exact population rank and numerical strength are separate issues. As the rows approach equality, the data identify the sum of coefficients more strongly than their difference. A small moment perturbation can move the coefficients substantially in opposite directions while preserving their sum.

For $`\epsilon=10^{-3}`$ and $`\delta=10^{-4}`$, the changed coefficients are $`(0.9,2.1)`$. At $`\epsilon=0`$ the matrix loses rank and the two coefficients cannot be separately identified.

These cross-moments can arise from valid instruments: take $`E[ZZ^\top]=I`$ and $`D=G^\top Z+v`$ with $`E[Zv]=0`$, and an outcome error orthogonal to $`Z`$. The example does not rely on an exclusion violation. It shows the ill-conditioning of identification even within a valid moment model. Sampling inference must account for how close the relevant rank direction is to zero.

#### Problem EC-89: a complete two-period fixed-effects calculation

For three individuals at two dates, let regressor pairs be $`(1,2),(2,4),(3,3)`$ and outcome pairs $`(2,5),(1,5),(4,3)`$. Fit one common slope and individual intercepts by fixed effects. Compute the slope, intercepts, and residuals.

Solution. First differences are $`\Delta x=(1,2,0)`$ and $`\Delta y=(3,4,-1)`$. With two dates, within OLS equals a through-origin regression of differences, giving

```math
\hat\beta=\frac{1(3)+2(4)+0(-1)}{1^2+2^2+0^2}=\frac{11}5.
```

To check the equivalence, each individual's demeaned regressors are $`(-\Delta x_i/2,\Delta x_i/2)`$ and similarly for outcomes. Summing products gives half the difference-regression numerator and denominator, leaving the same ratio.

Recover $`\hat\alpha_i=\bar y_i-\hat\beta\bar x_i`$. The three values are $`1/5,-18/5,-31/10`$. The residual pairs are

```math
(-2/5,2/5),\qquad(1/5,-1/5),\qquad(1/2,-1/2).
```

Each pair sums to zero. The first two individuals' contributions to the within regressor-residual inner product are $`2/5`$ and $`-2/5`$, while the third contributes zero. Thus the slope normal equation holds. The total residual sum of squares is $`9/10`$.

The third individual has no within-regressor variation and supplies no slope denominator. Its outcomes still determine its intercept and contribute residual information. With all three individual intercepts included, a separate unrestricted common intercept would be redundant; a normalization or omitted dummy would be required.

The arithmetic estimates the within projection. A structural interpretation needs the appropriate strict-exogeneity and sampling assumptions, which cannot be established from these six observations alone.

#### Problem EC-90: uncorrelated returns with dependent squares

Let $`z_t`$ be iid standard normal and consider $`r_t=\sigma_tz_t`$ with $`\sigma_t^2=\omega+\alpha r_{t-1}^2`$, $`\omega\gt 0`$, $`0\le\alpha\lt 1`$. Construct a stationary causal variance process, derive return variance, and determine when the fourth moment exists. Relate return and squared-return correlations.

Solution. Substituting $`r_{t-1}^2=\sigma_{t-1}^2z_{t-1}^2`$ repeatedly gives

```math
\sigma_t^2=\omega\left[
1+\sum_{j=1}^\infty\alpha^j\prod_{k=1}^jz_{t-k}^2
\right].
```

All terms are nonnegative. Independence and $`E[z^2]=1`$ give finite expected sum $`\omega/(1-\alpha)`$, so the series is finite almost surely. It is a time-shifted function of iid past shocks and therefore stationary and causal. Substitution verifies the variance recursion.

Conditional on the past, $`E[r_t\mid\mathcal F_{t-1}]=0`$ and $`E[r_t^2\mid\mathcal F_{t-1}]=\sigma_t^2`$. Thus return variance is $`v=\omega/(1-\alpha)`$. Every nonzero-lag return covariance is zero by conditioning on the earlier information.

For the fourth moment, $`E[z^4]=3`$. If finite, $`m_4=E[r_t^4]`$ satisfies

```math
m_4=3(\omega^2+2\omega\alpha v+\alpha^2m_4).
```

Consequently

```math
m_4=\frac{3\omega^2(1+\alpha)}
{(1-\alpha)(1-3\alpha^2)}
```

when $`3\alpha^2\lt 1`$. This condition is also sufficient: the $`L^2`$ norms of terms in the variance series are bounded by a geometric series with ratio $`\alpha\sqrt3\lt 1`$. If $`3\alpha^2\ge1`$, the moment equation with its strictly positive constant term rules out a finite fourth moment.

When finite, kurtosis is $`3(1-\alpha^2)/(1-3\alpha^2)`$. Also $`r_t^2=\omega+\alpha r_{t-1}^2+\eta_t`$ with $`\eta_t=\sigma_t^2(z_t^2-1)`$ a martingale difference. The centered squares therefore have autocorrelation $`\alpha^k`$ at lag $`k`$. At $`\alpha=1/2`$, returns have zero linear autocorrelation while squared-return correlations decay geometrically and kurtosis equals nine.

Zero return autocorrelation has not established independence. The conditional variance carries dependence that is visible in the squares.

[Contents](README.md) · [Part 1](04-econometrics-and-time-series-part-1.md)
