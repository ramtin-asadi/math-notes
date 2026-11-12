# Stochastic Calculus — Part 2

Brownian motion has continuous paths but does not behave like a differentiable function. Its increments are small enough to have zero mean and large enough that their squares accumulate. Stochastic calculus keeps track of that accumulated second-order contribution.

Sources:

- Steven Lalley, [Brownian Motion and Stochastic Calculus, University of Chicago](https://galton.uchicago.edu/~lalley/Courses/385/2014/index.html), including the linked lecture notes and the assignment PDFs currently dated 2016. Citations identify the actual document and problem rather than assuming its date from the older index.
- [MIT 18.642, Topics in Mathematics with Applications in Finance, Fall 2024](https://ocw.mit.edu/courses/18-642-topics-in-mathematics-with-applications-in-finance-fall-2024/), especially the stochastic-calculus and pricing lectures.
- The prerequisite treatment of [conditional expectation and Gaussian variables](README.md).

Unless otherwise stated, $`W`$ is a standard Brownian motion relative to a filtration $`(\mathcal F_t)`$ on a probability space $`(\Omega,\mathcal F,P)`$. The filtration contains the information already observed and satisfies the usual completeness and right-continuity conditions. Equalities between stochastic processes are understood up to the appropriate almost-sure versions.

[Contents](README.md) · [Part 1](06-stochastic-calculus-part-1.md)

Contents:

- [8. Further worked problems](#8-further-worked-problems)

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

[Contents](README.md) · [Part 1](06-stochastic-calculus-part-1.md)
