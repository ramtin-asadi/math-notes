# Dynamic Optimization — Part 2

A dynamic problem links today's choice to tomorrow's available choices. The state records the information needed for that link. The action is what can be chosen now. A value function records the best achievable future objective from a specified state, including the current period under the convention used here.

Sources:

- Guido Lorenzoni, [MIT 14.451, Dynamic Optimization Methods with Applications, Fall 2009](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/), lecture notes and Problem Sets 1–4.
- Thomas Sargent and John Stachurski, [Dynamic Programming](https://dp.quantecon.org/), especially the volume on finite states, fixed points, Markov chains, and decision processes.

Discount factors satisfy $`0\lt \beta\lt 1`$ unless a finite-horizon exception is stated. A stochastic transition is conditional on the current state and chosen action. Examples are constructed here; source problems and subparts are identified individually.

[Contents](README.md) · [Part 1](05-dynamic-optimization-part-1.md)

Contents:

- [7. Further worked problems](#7-further-worked-problems)

## 7. Further worked problems

The problems in this section are constructed for these notes. Each statement supplies its data and assumptions. Closely related cases are placed together so that the changes in the calculation can be followed directly.

### 7.1 Finite-horizon consumption and the budget

An Euler equation compares neighboring dates, but the total budget determines the scale of consumption. With log utility and a known interest factor, the full allocation can be obtained directly from a finite constrained problem. This also checks terminal wealth and feasibility.

#### Problem DO-45: log consumption over 3 dates, case 1

Maximize $`\sum_{t=0}^{2}\beta^t\log c_t`$ with positive consumption and present-value budget $`\sum_{t=0}^{2}c_t/R^t\le w`$. Let $`\beta=\frac{1}{2}`$, $`R=2`$, and $`w=12`$. Find the entire allocation.

Solution. The budget binds because utility is strictly increasing in every consumption coordinate. Introduce a positive multiplier $`\lambda`$ on present-value spending. The first-order conditions are

```math
\frac{\beta^t}{c_t}=\frac\lambda{R^t},
\qquad c_t=\frac{(\beta R)^t}\lambda.
```

Consequently $`c_{t+1}/c_t=\beta R`$. Substituting the full sequence into the budget determines its initial level:

```math
w=c_0\sum_{t=0}^{m-1}\beta^t,\qquad
c_0=\frac w{\sum_{t=0}^{m-1}\beta^t}=\frac{48}{7}.
```

The consumption sequence, from date zero to the last date, is

```math
(c_0,\ldots,c_{m-1})=(\frac{48}{7}, \frac{48}{7}, \frac{48}{7}).
```

Every entry is positive. Its discounted cost is exactly $`w=12`$, and its utility is approximately $`3.369259`$.

For a sequential feasibility check, set $`a_0=w`$ and $`a_{t+1}=R(a_t-c_t)`$. The resulting asset balances, including the terminal one, are

```math
(a_0,\ldots,a_m)=(12, \frac{72}{7}, \frac{48}{7}, 0).
```

The balances before the last date fund the remaining planned consumption and are nonnegative. The terminal balance is zero: leaving positive wealth after the final date would waste utility.

Finally, the objective is strictly concave on the positive orthant and the budget set is convex. The first-order conditions therefore identify the unique global optimum. As a consumption coordinate approaches zero, its log utility tends to minus infinity, excluding a zero-consumption boundary optimum. The Euler equation alone would have left the multiplier, and hence every consumption level, undetermined.

#### Problem DO-46: log consumption over 4 dates, case 2

Maximize $`\sum_{t=0}^{3}\beta^t\log c_t`$ with positive consumption and present-value budget $`\sum_{t=0}^{3}c_t/R^t\le w`$. Let $`\beta=\frac{3}{4}`$, $`R=1`$, and $`w=10`$. Find the entire allocation.

Solution. The budget binds because utility is strictly increasing in every consumption coordinate. Introduce a positive multiplier $`\lambda`$ on present-value spending. The first-order conditions are

```math
\frac{\beta^t}{c_t}=\frac\lambda{R^t},
\qquad c_t=\frac{(\beta R)^t}\lambda.
```

Consequently $`c_{t+1}/c_t=\beta R`$. Substituting the full sequence into the budget determines its initial level:

```math
w=c_0\sum_{t=0}^{m-1}\beta^t,\qquad
c_0=\frac w{\sum_{t=0}^{m-1}\beta^t}=\frac{128}{35}.
```

The consumption sequence, from date zero to the last date, is

```math
(c_0,\ldots,c_{m-1})=(\frac{128}{35}, \frac{96}{35}, \frac{72}{35}, \frac{54}{35}).
```

Every entry is positive. Its discounted cost is exactly $`w=10`$, and its utility is approximately $`2.642114`$.

For a sequential feasibility check, set $`a_0=w`$ and $`a_{t+1}=R(a_t-c_t)`$. The resulting asset balances, including the terminal one, are

```math
(a_0,\ldots,a_m)=(10, \frac{222}{35}, \frac{18}{5}, \frac{54}{35}, 0).
```

The balances before the last date fund the remaining planned consumption and are nonnegative. The terminal balance is zero: leaving positive wealth after the final date would waste utility.

Finally, the objective is strictly concave on the positive orthant and the budget set is convex. The first-order conditions therefore identify the unique global optimum. As a consumption coordinate approaches zero, its log utility tends to minus infinity, excluding a zero-consumption boundary optimum. The Euler equation alone would have left the multiplier, and hence every consumption level, undetermined.

#### Problem DO-47: log consumption over 3 dates, case 3

Maximize $`\sum_{t=0}^{2}\beta^t\log c_t`$ with positive consumption and present-value budget $`\sum_{t=0}^{2}c_t/R^t\le w`$. Let $`\beta=\frac{9}{10}`$, $`R=\frac{10}{9}`$, and $`w=9`$. Find the entire allocation.

Solution. The budget binds because utility is strictly increasing in every consumption coordinate. Introduce a positive multiplier $`\lambda`$ on present-value spending. The first-order conditions are

```math
\frac{\beta^t}{c_t}=\frac\lambda{R^t},
\qquad c_t=\frac{(\beta R)^t}\lambda.
```

Consequently $`c_{t+1}/c_t=\beta R`$. Substituting the full sequence into the budget determines its initial level:

```math
w=c_0\sum_{t=0}^{m-1}\beta^t,\qquad
c_0=\frac w{\sum_{t=0}^{m-1}\beta^t}=\frac{900}{271}.
```

The consumption sequence, from date zero to the last date, is

```math
(c_0,\ldots,c_{m-1})=(\frac{900}{271}, \frac{900}{271}, \frac{900}{271}).
```

Every entry is positive. Its discounted cost is exactly $`w=9`$, and its utility is approximately $`3.252748`$.

For a sequential feasibility check, set $`a_0=w`$ and $`a_{t+1}=R(a_t-c_t)`$. The resulting asset balances, including the terminal one, are

```math
(a_0,\ldots,a_m)=(9, \frac{1710}{271}, \frac{900}{271}, 0).
```

The balances before the last date fund the remaining planned consumption and are nonnegative. The terminal balance is zero: leaving positive wealth after the final date would waste utility.

Finally, the objective is strictly concave on the positive orthant and the budget set is convex. The first-order conditions therefore identify the unique global optimum. As a consumption coordinate approaches zero, its log utility tends to minus infinity, excluding a zero-consumption boundary optimum. The Euler equation alone would have left the multiplier, and hence every consumption level, undetermined.

#### Problem DO-48: log consumption over 4 dates, case 4

Maximize $`\sum_{t=0}^{3}\beta^t\log c_t`$ with positive consumption and present-value budget $`\sum_{t=0}^{3}c_t/R^t\le w`$. Let $`\beta=\frac{2}{3}`$, $`R=\frac{3}{2}`$, and $`w=16`$. Find the entire allocation.

Solution. The budget binds because utility is strictly increasing in every consumption coordinate. Introduce a positive multiplier $`\lambda`$ on present-value spending. The first-order conditions are

```math
\frac{\beta^t}{c_t}=\frac\lambda{R^t},
\qquad c_t=\frac{(\beta R)^t}\lambda.
```

Consequently $`c_{t+1}/c_t=\beta R`$. Substituting the full sequence into the budget determines its initial level:

```math
w=c_0\sum_{t=0}^{m-1}\beta^t,\qquad
c_0=\frac w{\sum_{t=0}^{m-1}\beta^t}=\frac{432}{65}.
```

The consumption sequence, from date zero to the last date, is

```math
(c_0,\ldots,c_{m-1})=(\frac{432}{65}, \frac{432}{65}, \frac{432}{65}, \frac{432}{65}).
```

Every entry is positive. Its discounted cost is exactly $`w=16`$, and its utility is approximately $`4.559722`$.

For a sequential feasibility check, set $`a_0=w`$ and $`a_{t+1}=R(a_t-c_t)`$. The resulting asset balances, including the terminal one, are

```math
(a_0,\ldots,a_m)=(16, \frac{912}{65}, \frac{144}{13}, \frac{432}{65}, 0).
```

The balances before the last date fund the remaining planned consumption and are nonnegative. The terminal balance is zero: leaving positive wealth after the final date would waste utility.

Finally, the objective is strictly concave on the positive orthant and the budget set is convex. The first-order conditions therefore identify the unique global optimum. As a consumption coordinate approaches zero, its log utility tends to minus infinity, excluding a zero-consumption boundary optimum. The Euler equation alone would have left the multiplier, and hence every consumption level, undetermined.

#### Problem DO-49: log consumption over 5 dates, case 5

Maximize $`\sum_{t=0}^{4}\beta^t\log c_t`$ with positive consumption and present-value budget $`\sum_{t=0}^{4}c_t/R^t\le w`$. Let $`\beta=\frac{4}{5}`$, $`R=\frac{5}{4}`$, and $`w=20`$. Find the entire allocation.

Solution. The budget binds because utility is strictly increasing in every consumption coordinate. Introduce a positive multiplier $`\lambda`$ on present-value spending. The first-order conditions are

```math
\frac{\beta^t}{c_t}=\frac\lambda{R^t},
\qquad c_t=\frac{(\beta R)^t}\lambda.
```

Consequently $`c_{t+1}/c_t=\beta R`$. Substituting the full sequence into the budget determines its initial level:

```math
w=c_0\sum_{t=0}^{m-1}\beta^t,\qquad
c_0=\frac w{\sum_{t=0}^{m-1}\beta^t}=\frac{12500}{2101}.
```

The consumption sequence, from date zero to the last date, is

```math
(c_0,\ldots,c_{m-1})=(\frac{12500}{2101}, \frac{12500}{2101}, \frac{12500}{2101}, \frac{12500}{2101}, \frac{12500}{2101}).
```

Every entry is positive. Its discounted cost is exactly $`w=20`$, and its utility is approximately $`5.994792`$.

For a sequential feasibility check, set $`a_0=w`$ and $`a_{t+1}=R(a_t-c_t)`$. The resulting asset balances, including the terminal one, are

```math
(a_0,\ldots,a_m)=(20, \frac{36900}{2101}, \frac{30500}{2101}, \frac{22500}{2101}, \frac{12500}{2101}, 0).
```

The balances before the last date fund the remaining planned consumption and are nonnegative. The terminal balance is zero: leaving positive wealth after the final date would waste utility.

Finally, the objective is strictly concave on the positive orthant and the budget set is convex. The first-order conditions therefore identify the unique global optimum. As a consumption coordinate approaches zero, its log utility tends to minus infinity, excluding a zero-consumption boundary optimum. The Euler equation alone would have left the multiplier, and hence every consumption level, undetermined.

### 7.2 Solving every stationary policy in a small MDP

For two states and two actions there are only four deterministic stationary policies. Each policy has its own linear value equation. Computing all four and then checking the Bellman inequalities gives a complete finite-state solution, including the policy that actually implements the value.

#### Problem DO-50: two-state discounted control, case 1

Use discount $`\beta=\frac{1}{2}`$. Action 0 has rewards

```math
\begin{pmatrix}2\\[0pt] 0\end{pmatrix}
```

 and transition

```math
P_0=\begin{pmatrix}\frac{3}{4}&\frac{1}{4}\\[0pt] \frac{1}{2}&\frac{1}{2}\end{pmatrix}
```

; action 1 has rewards

```math
\begin{pmatrix}0\\[0pt] 3\end{pmatrix}
```

 and transition

```math
P_1=\begin{pmatrix}\frac{1}{4}&\frac{3}{4}\\[0pt] \frac{1}{4}&\frac{3}{4}\end{pmatrix}
```

. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta\lt 1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`\frac{24}{7}`$ | $`\frac{8}{7}`$ |
| $`(0,1)`$ | $`\frac{13}{3}`$ | $`\frac{17}{3}`$ |
| $`(1,0)`$ | $`0`$ | $`0`$ |
| $`(1,1)`$ | $`\frac{9}{4}`$ | $`\frac{21}{4}`$ |

The policy $`\pi^*=(0,1)`$ has value

```math
v^*=\begin{pmatrix}\frac{13}{3}\\[0pt] \frac{17}{3}\end{pmatrix}
```

. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=\frac{13}{3},\quad Q(1,1)=\frac{8}{3},\qquad
Q(2,0)=\frac{5}{2},\quad Q(2,1)=\frac{17}{3}.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

#### Problem DO-51: two-state discounted control, case 2

Use discount $`\beta=\frac{4}{5}`$. Action 0 has rewards

```math
\begin{pmatrix}1\\[0pt] 2\end{pmatrix}
```

 and transition

```math
P_0=\begin{pmatrix}\frac{4}{5}&\frac{1}{5}\\[0pt] \frac{1}{5}&\frac{4}{5}\end{pmatrix}
```

; action 1 has rewards

```math
\begin{pmatrix}3\\[0pt] 0\end{pmatrix}
```

 and transition

```math
P_1=\begin{pmatrix}\frac{1}{2}&\frac{1}{2}\\[0pt] \frac{3}{4}&\frac{1}{4}\end{pmatrix}
```

. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta\lt 1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`\frac{85}{13}`$ | $`\frac{110}{13}`$ |
| $`(0,1)`$ | $`\frac{25}{6}`$ | $`\frac{25}{8}`$ |
| $`(1,0)`$ | $`\frac{235}{19}`$ | $`\frac{210}{19}`$ |
| $`(1,1)`$ | $`10`$ | $`\frac{15}{2}`$ |

The policy $`\pi^*=(1,0)`$ has value

```math
v^*=\begin{pmatrix}\frac{235}{19}\\[0pt] \frac{210}{19}\end{pmatrix}
```

. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=\frac{203}{19},\quad Q(1,1)=\frac{235}{19},\qquad
Q(2,0)=\frac{210}{19},\quad Q(2,1)=\frac{183}{19}.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

#### Problem DO-52: two-state discounted control, case 3

Use discount $`\beta=\frac{2}{3}`$. Action 0 has rewards

```math
\begin{pmatrix}0\\[0pt] 4\end{pmatrix}
```

 and transition

```math
P_0=\begin{pmatrix}\frac{1}{2}&\frac{1}{2}\\[0pt] \frac{1}{4}&\frac{3}{4}\end{pmatrix}
```

; action 1 has rewards

```math
\begin{pmatrix}2\\[0pt] 1\end{pmatrix}
```

 and transition

```math
P_1=\begin{pmatrix}\frac{3}{4}&\frac{1}{4}\\[0pt] \frac{1}{2}&\frac{1}{2}\end{pmatrix}
```

. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta\lt 1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`\frac{24}{5}`$ | $`\frac{48}{5}`$ |
| $`(0,1)`$ | $`1`$ | $`2`$ |
| $`(1,0)`$ | $`\frac{15}{2}`$ | $`\frac{21}{2}`$ |
| $`(1,1)`$ | $`\frac{27}{5}`$ | $`\frac{21}{5}`$ |

The policy $`\pi^*=(1,0)`$ has value

```math
v^*=\begin{pmatrix}\frac{15}{2}\\[0pt] \frac{21}{2}\end{pmatrix}
```

. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=6,\quad Q(1,1)=\frac{15}{2},\qquad
Q(2,0)=\frac{21}{2},\quad Q(2,1)=7.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

#### Problem DO-53: two-state discounted control, case 4

Use discount $`\beta=\frac{3}{4}`$. Action 0 has rewards

```math
\begin{pmatrix}-1\\[0pt] 2\end{pmatrix}
```

 and transition

```math
P_0=\begin{pmatrix}\frac{9}{10}&\frac{1}{10}\\[0pt] \frac{1}{2}&\frac{1}{2}\end{pmatrix}
```

; action 1 has rewards

```math
\begin{pmatrix}1\\[0pt] 1\end{pmatrix}
```

 and transition

```math
P_1=\begin{pmatrix}\frac{1}{5}&\frac{4}{5}\\[0pt] \frac{1}{10}&\frac{9}{10}\end{pmatrix}
```

. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta\lt 1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`-\frac{19}{7}`$ | $`\frac{11}{7}`$ |
| $`(0,1)`$ | $`-\frac{5}{2}`$ | $`\frac{5}{2}`$ |
| $`(1,0)`$ | $`\frac{292}{49}`$ | $`\frac{332}{49}`$ |
| $`(1,1)`$ | $`4`$ | $`4`$ |

The policy $`\pi^*=(1,0)`$ has value

```math
v^*=\begin{pmatrix}\frac{292}{49}\\[0pt] \frac{332}{49}\end{pmatrix}
```

. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=\frac{173}{49},\quad Q(1,1)=\frac{292}{49},\qquad
Q(2,0)=\frac{332}{49},\quad Q(2,1)=\frac{295}{49}.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

#### Problem DO-54: two-state discounted control, case 5

Use discount $`\beta=\frac{9}{10}`$. Action 0 has rewards

```math
\begin{pmatrix}3\\[0pt] 1\end{pmatrix}
```

 and transition

```math
P_0=\begin{pmatrix}\frac{2}{3}&\frac{1}{3}\\[0pt] \frac{1}{3}&\frac{2}{3}\end{pmatrix}
```

; action 1 has rewards

```math
\begin{pmatrix}1\\[0pt] 4\end{pmatrix}
```

 and transition

```math
P_1=\begin{pmatrix}\frac{1}{3}&\frac{2}{3}\\[0pt] \frac{2}{3}&\frac{1}{3}\end{pmatrix}
```

. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta\lt 1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`\frac{150}{7}`$ | $`\frac{130}{7}`$ |
| $`(0,1)`$ | $`33`$ | $`34`$ |
| $`(1,0)`$ | $`10`$ | $`10`$ |
| $`(1,1)`$ | $`\frac{310}{13}`$ | $`\frac{340}{13}`$ |

The policy $`\pi^*=(0,1)`$ has value

```math
v^*=\begin{pmatrix}33\\[0pt] 34\end{pmatrix}
```

. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=33,\quad Q(1,1)=\frac{313}{10},\qquad
Q(2,0)=\frac{313}{10},\quad Q(2,1)=34.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

### 7.3 Stationary distributions and actual convergence

A stationary distribution remains unchanged after one transition. Convergence from a given initial state is a separate assertion. A two-state chain permits an exact formula showing why periodicity can prevent convergence even when the stationary distribution is unique.

#### Problem DO-55: two-state Markov law, case 1

Let

```math
P=\begin{pmatrix}\frac{4}{5}&\frac{1}{5}\\[0pt] \frac{2}{5}&\frac{3}{5}\end{pmatrix}
```

 on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q\gt 0`$, the unique solution is

```math
\pi_0=\frac q{p+q}=\frac{2}{3},\qquad
\pi_1=\frac p{p+q}=\frac{1}{3}.
```

For the time-dependent law, set $`z_n=P(X_n=1\mid X_0=0)`$. Conditioning on the current state gives

```math
z_{n+1}=p(1-z_n)+(1-q)z_n=p+(1-p-q)z_n.
```

Subtract its fixed point $`\pi_1`$ to obtain $`z_{n+1}-\pi_1=(1-p-q)(z_n-\pi_1)`$. Since $`z_0=0`$,

```math
z_n=\pi_1[1-(1-p-q)^n]
=\frac{1}{3}[1-(\frac{2}{5})^n].
```

The date-$`n`$ distribution is $`(1-z_n,z_n)`$. Substituting $`n=0`$ recovers the initial state, and substituting $`n=1`$ gives $`z_1=p`$, checking the recursion.

Here $`\vert 1-p-q\vert =\frac{2}{5}\lt 1`$, so the power term vanishes and the state distribution converges to $`(\frac{2}{3},\frac{1}{3})`$. A negative second eigenvalue gives alternating deviations; a positive one gives deviations with a constant sign.

The eigenvalues of the transition matrix are one and $`1-p-q=\frac{2}{5}`$. The first eigenvalue preserves probability mass; the second controls how deviations from stationarity evolve. This calculation shows why a proof of an invariant distribution alone does not prove convergence to it. Irreducibility ensures uniqueness here, while aperiodicity rules out the alternating exceptional case.

#### Problem DO-56: two-state Markov law, case 2

Let

```math
P=\begin{pmatrix}\frac{1}{2}&\frac{1}{2}\\[0pt] \frac{1}{2}&\frac{1}{2}\end{pmatrix}
```

 on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q\gt 0`$, the unique solution is

```math
\pi_0=\frac q{p+q}=\frac{1}{2},\qquad
\pi_1=\frac p{p+q}=\frac{1}{2}.
```

For the time-dependent law, set $`z_n=P(X_n=1\mid X_0=0)`$. Conditioning on the current state gives

```math
z_{n+1}=p(1-z_n)+(1-q)z_n=p+(1-p-q)z_n.
```

Subtract its fixed point $`\pi_1`$ to obtain $`z_{n+1}-\pi_1=(1-p-q)(z_n-\pi_1)`$. Since $`z_0=0`$,

```math
z_n=\pi_1[1-(1-p-q)^n]
=\frac{1}{2}[1-(0)^n].
```

The date-$`n`$ distribution is $`(1-z_n,z_n)`$. Substituting $`n=0`$ recovers the initial state, and substituting $`n=1`$ gives $`z_1=p`$, checking the recursion.

Here $`\vert 1-p-q\vert =0\lt 1`$, so the power term vanishes and the state distribution converges to $`(\frac{1}{2},\frac{1}{2})`$. A negative second eigenvalue gives alternating deviations; a positive one gives deviations with a constant sign.

The eigenvalues of the transition matrix are one and $`1-p-q=0`$. The first eigenvalue preserves probability mass; the second controls how deviations from stationarity evolve. This calculation shows why a proof of an invariant distribution alone does not prove convergence to it. Irreducibility ensures uniqueness here, while aperiodicity rules out the alternating exceptional case.

#### Problem DO-57: two-state Markov law, case 3

Let

```math
P=\begin{pmatrix}\frac{9}{10}&\frac{1}{10}\\[0pt] \frac{3}{10}&\frac{7}{10}\end{pmatrix}
```

 on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q\gt 0`$, the unique solution is

```math
\pi_0=\frac q{p+q}=\frac{3}{4},\qquad
\pi_1=\frac p{p+q}=\frac{1}{4}.
```

For the time-dependent law, set $`z_n=P(X_n=1\mid X_0=0)`$. Conditioning on the current state gives

```math
z_{n+1}=p(1-z_n)+(1-q)z_n=p+(1-p-q)z_n.
```

Subtract its fixed point $`\pi_1`$ to obtain $`z_{n+1}-\pi_1=(1-p-q)(z_n-\pi_1)`$. Since $`z_0=0`$,

```math
z_n=\pi_1[1-(1-p-q)^n]
=\frac{1}{4}[1-(\frac{3}{5})^n].
```

The date-$`n`$ distribution is $`(1-z_n,z_n)`$. Substituting $`n=0`$ recovers the initial state, and substituting $`n=1`$ gives $`z_1=p`$, checking the recursion.

Here $`\vert 1-p-q\vert =\frac{3}{5}\lt 1`$, so the power term vanishes and the state distribution converges to $`(\frac{3}{4},\frac{1}{4})`$. A negative second eigenvalue gives alternating deviations; a positive one gives deviations with a constant sign.

The eigenvalues of the transition matrix are one and $`1-p-q=\frac{3}{5}`$. The first eigenvalue preserves probability mass; the second controls how deviations from stationarity evolve. This calculation shows why a proof of an invariant distribution alone does not prove convergence to it. Irreducibility ensures uniqueness here, while aperiodicity rules out the alternating exceptional case.

#### Problem DO-58: two-state Markov law, case 4

Let

```math
P=\begin{pmatrix}0&1\\[0pt] 1&0\end{pmatrix}
```

 on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q\gt 0`$, the unique solution is

```math
\pi_0=\frac q{p+q}=\frac{1}{2},\qquad
\pi_1=\frac p{p+q}=\frac{1}{2}.
```

For the time-dependent law, set $`z_n=P(X_n=1\mid X_0=0)`$. Conditioning on the current state gives

```math
z_{n+1}=p(1-z_n)+(1-q)z_n=p+(1-p-q)z_n.
```

Subtract its fixed point $`\pi_1`$ to obtain $`z_{n+1}-\pi_1=(1-p-q)(z_n-\pi_1)`$. Since $`z_0=0`$,

```math
z_n=\pi_1[1-(1-p-q)^n]
=\frac{1}{2}[1-(-1)^n].
```

The date-$`n`$ distribution is $`(1-z_n,z_n)`$. Substituting $`n=0`$ recovers the initial state, and substituting $`n=1`$ gives $`z_1=p`$, checking the recursion.

Here the second eigenvalue is minus one. Starting in state zero, the chain alternates deterministically, so the probability of state one alternates between zero and one. It does not converge. Its time averages still approach the stationary proportions.

The eigenvalues of the transition matrix are one and $`1-p-q=-1`$. The first eigenvalue preserves probability mass; the second controls how deviations from stationarity evolve. This calculation shows why a proof of an invariant distribution alone does not prove convergence to it. Irreducibility ensures uniqueness here, while aperiodicity rules out the alternating exceptional case.

#### Problem DO-59: two-state Markov law, case 5

Let

```math
P=\begin{pmatrix}0&1\\[0pt] \frac{1}{2}&\frac{1}{2}\end{pmatrix}
```

 on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q\gt 0`$, the unique solution is

```math
\pi_0=\frac q{p+q}=\frac{1}{3},\qquad
\pi_1=\frac p{p+q}=\frac{2}{3}.
```

For the time-dependent law, set $`z_n=P(X_n=1\mid X_0=0)`$. Conditioning on the current state gives

```math
z_{n+1}=p(1-z_n)+(1-q)z_n=p+(1-p-q)z_n.
```

Subtract its fixed point $`\pi_1`$ to obtain $`z_{n+1}-\pi_1=(1-p-q)(z_n-\pi_1)`$. Since $`z_0=0`$,

```math
z_n=\pi_1[1-(1-p-q)^n]
=\frac{2}{3}[1-(-\frac{1}{2})^n].
```

The date-$`n`$ distribution is $`(1-z_n,z_n)`$. Substituting $`n=0`$ recovers the initial state, and substituting $`n=1`$ gives $`z_1=p`$, checking the recursion.

Here $`\vert 1-p-q\vert =\frac{1}{2}\lt 1`$, so the power term vanishes and the state distribution converges to $`(\frac{1}{3},\frac{2}{3})`$. A negative second eigenvalue gives alternating deviations; a positive one gives deviations with a constant sign.

The eigenvalues of the transition matrix are one and $`1-p-q=-\frac{1}{2}`$. The first eigenvalue preserves probability mass; the second controls how deviations from stationarity evolve. This calculation shows why a proof of an invariant distribution alone does not prove convergence to it. Irreducibility ensures uniqueness here, while aperiodicity rules out the alternating exceptional case.

### 7.4 Stopping with two possible offers

At each date an offer is observed. Accepting ends the problem and pays the offer immediately. Rejecting costs $`c`$ now and gives the discounted value of a fresh independent offer next date. Writing the timing explicitly determines the continuation value.

#### Problem DO-60: two-offer stopping rule, case 1

The offer is $`H=10`$ with probability $`p=\frac{1}{2}`$ and $`L=2`$ otherwise. Let $`\beta=\frac{9}{10}`$ and rejection cost $`c=1`$. Find the value before seeing the offer and an optimal acceptance rule.

Solution. Let $`v`$ be the value before the next offer is observed. After an offer $`x`$ is seen, acceptance gives $`x`$ and rejection gives the common continuation value $`C=-c+\beta v`$. Therefore

```math
v=p\max(H,C)+(1-p)\max(L,C).
```

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H\lt H`$.

There are two remaining candidate policies. If both offers are accepted, $`v=pH+(1-p)L`$, and consistency requires $`-c+\beta v\le L`$. If only the high offer is accepted, then

```math
v=pH+(1-p)(-c+\beta v),
\qquad
v=\frac{pH-(1-p)c}{1-(1-p)\beta}.
```

This second candidate is consistent when its continuation value is at least $`L`$. Checking the inequalities for the supplied data gives

```math
v=\frac{90}{11},\qquad C=-c+\beta v=\frac{70}{11}.
```

The optimal rule is to reject the low offer and accept the high offer. At equality with the low offer, either accepting it or rejecting it gives the same value.

The right-hand side of the scalar Bellman equation is a contraction with modulus no greater than $`\beta`$, so a consistent fixed point is unique. The cost is paid only after rejection in this timing convention. Charging it before every observation would alter both candidate equations and the resulting threshold.

#### Problem DO-61: two-offer stopping rule, case 2

The offer is $`H=8`$ with probability $`p=\frac{1}{5}`$ and $`L=5`$ otherwise. Let $`\beta=\frac{4}{5}`$ and rejection cost $`c=1`$. Find the value before seeing the offer and an optimal acceptance rule.

Solution. Let $`v`$ be the value before the next offer is observed. After an offer $`x`$ is seen, acceptance gives $`x`$ and rejection gives the common continuation value $`C=-c+\beta v`$. Therefore

```math
v=p\max(H,C)+(1-p)\max(L,C).
```

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H\lt H`$.

There are two remaining candidate policies. If both offers are accepted, $`v=pH+(1-p)L`$, and consistency requires $`-c+\beta v\le L`$. If only the high offer is accepted, then

```math
v=pH+(1-p)(-c+\beta v),
\qquad
v=\frac{pH-(1-p)c}{1-(1-p)\beta}.
```

This second candidate is consistent when its continuation value is at least $`L`$. Checking the inequalities for the supplied data gives

```math
v=\frac{28}{5},\qquad C=-c+\beta v=\frac{87}{25}.
```

The optimal rule is to accept both offers. At equality with the low offer, either accepting it or rejecting it gives the same value.

The right-hand side of the scalar Bellman equation is a contraction with modulus no greater than $`\beta`$, so a consistent fixed point is unique. The cost is paid only after rejection in this timing convention. Charging it before every observation would alter both candidate equations and the resulting threshold.

#### Problem DO-62: two-offer stopping rule, case 3

The offer is $`H=12`$ with probability $`p=\frac{1}{4}`$ and $`L=1`$ otherwise. Let $`\beta=\frac{9}{10}`$ and rejection cost $`c=0`$. Find the value before seeing the offer and an optimal acceptance rule.

Solution. Let $`v`$ be the value before the next offer is observed. After an offer $`x`$ is seen, acceptance gives $`x`$ and rejection gives the common continuation value $`C=-c+\beta v`$. Therefore

```math
v=p\max(H,C)+(1-p)\max(L,C).
```

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H\lt H`$.

There are two remaining candidate policies. If both offers are accepted, $`v=pH+(1-p)L`$, and consistency requires $`-c+\beta v\le L`$. If only the high offer is accepted, then

```math
v=pH+(1-p)(-c+\beta v),
\qquad
v=\frac{pH-(1-p)c}{1-(1-p)\beta}.
```

This second candidate is consistent when its continuation value is at least $`L`$. Checking the inequalities for the supplied data gives

```math
v=\frac{120}{13},\qquad C=-c+\beta v=\frac{108}{13}.
```

The optimal rule is to reject the low offer and accept the high offer. At equality with the low offer, either accepting it or rejecting it gives the same value.

The right-hand side of the scalar Bellman equation is a contraction with modulus no greater than $`\beta`$, so a consistent fixed point is unique. The cost is paid only after rejection in this timing convention. Charging it before every observation would alter both candidate equations and the resulting threshold.

#### Problem DO-63: two-offer stopping rule, case 4

The offer is $`H=6`$ with probability $`p=\frac{3}{4}`$ and $`L=4`$ otherwise. Let $`\beta=\frac{1}{2}`$ and rejection cost $`c=1`$. Find the value before seeing the offer and an optimal acceptance rule.

Solution. Let $`v`$ be the value before the next offer is observed. After an offer $`x`$ is seen, acceptance gives $`x`$ and rejection gives the common continuation value $`C=-c+\beta v`$. Therefore

```math
v=p\max(H,C)+(1-p)\max(L,C).
```

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H\lt H`$.

There are two remaining candidate policies. If both offers are accepted, $`v=pH+(1-p)L`$, and consistency requires $`-c+\beta v\le L`$. If only the high offer is accepted, then

```math
v=pH+(1-p)(-c+\beta v),
\qquad
v=\frac{pH-(1-p)c}{1-(1-p)\beta}.
```

This second candidate is consistent when its continuation value is at least $`L`$. Checking the inequalities for the supplied data gives

```math
v=\frac{11}{2},\qquad C=-c+\beta v=\frac{7}{4}.
```

The optimal rule is to accept both offers. At equality with the low offer, either accepting it or rejecting it gives the same value.

The right-hand side of the scalar Bellman equation is a contraction with modulus no greater than $`\beta`$, so a consistent fixed point is unique. The cost is paid only after rejection in this timing convention. Charging it before every observation would alter both candidate equations and the resulting threshold.

#### Problem DO-64: two-offer stopping rule, case 5

The offer is $`H=20`$ with probability $`p=\frac{1}{10}`$ and $`L=0`$ otherwise. Let $`\beta=\frac{19}{20}`$ and rejection cost $`c=2`$. Find the value before seeing the offer and an optimal acceptance rule.

Solution. Let $`v`$ be the value before the next offer is observed. After an offer $`x`$ is seen, acceptance gives $`x`$ and rejection gives the common continuation value $`C=-c+\beta v`$. Therefore

```math
v=p\max(H,C)+(1-p)\max(L,C).
```

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H\lt H`$.

There are two remaining candidate policies. If both offers are accepted, $`v=pH+(1-p)L`$, and consistency requires $`-c+\beta v\le L`$. If only the high offer is accepted, then

```math
v=pH+(1-p)(-c+\beta v),
\qquad
v=\frac{pH-(1-p)c}{1-(1-p)\beta}.
```

This second candidate is consistent when its continuation value is at least $`L`$. Checking the inequalities for the supplied data gives

```math
v=2,\qquad C=-c+\beta v=-\frac{1}{10}.
```

The optimal rule is to accept both offers. At equality with the low offer, either accepting it or rejecting it gives the same value.

The right-hand side of the scalar Bellman equation is a contraction with modulus no greater than $`\beta`$, so a consistent fixed point is unique. The cost is paid only after rejection in this timing convention. Charging it before every observation would alter both candidate equations and the resulting threshold.

### 7.5 Log growth with full depreciation

The log-growth model is one of the few nonlinear dynamic programs with an explicit value and policy. Solving it for several parameter sets allows the policy, steady state, Euler equation, and transversality condition to be checked together.

#### Problem DO-65: explicit growth solution, case 1

Maximize $`\sum_{t\ge0}\beta^t\log c_t`$ subject to $`c_t+k_{t+1}=Ak_t^\alpha`$, positive consumption and capital, with $`\alpha=\frac{1}{2}`$, $`\beta=\frac{4}{5}`$, and $`A=2`$. Derive the candidate value, policy, and steady state.

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0\lt k'\lt Ak^\alpha`$ and tends to minus infinity at both endpoints.

The first-order condition is

```math
-\frac1{Ak^\alpha-k'}+\frac{\beta B}{k'}=0.
```

Matching the coefficient of $`\log k`$ after substitution gives $`B=\alpha(1+\beta B)`$, so $`B=\alpha/(1-\alpha\beta)=\frac{5}{6}`$. Consequently

```math
k'=\alpha\beta Ak^\alpha=\frac{4}{5}k^{\frac{1}{2}},
\qquad
c=(1-\alpha\beta)Ak^\alpha=\frac{6}{5}k^{\frac{1}{2}}.
```

The constant equation yields

```math
a=\frac{\log(1-\alpha\beta)+[\log A+\alpha\beta\log(\alpha\beta)]/(1-\alpha\beta)}{1-\beta}
\approx 0.167796.
```

The positive steady state solves $`k^*=\alpha\beta A(k^*)^\alpha`$:

```math
k^*=(\alpha\beta A)^{1/(1-\alpha)}=\frac{16}{25},
\qquad c^*\approx 0.960000.
```

For a direct Euler check, the resource policy implies $`k_{t+1}/c_t=\alpha\beta/(1-\alpha\beta)`$. Substituting the next-period consumption rule into $`1/c_t=\beta A\alpha k_{t+1}^{\alpha-1}/c_{t+1}`$ verifies equality.

The discounted terminal marginal-value expression is $`\beta^t k_{t+1}/c_t=\beta^t\alpha\beta/(1-\alpha\beta)\to0`$. This is the relevant transversality check. Since log utility is unbounded below and the state set is not a compact bounded-utility model, a generic sup-norm contraction claim is insufficient. The earlier verification argument, using the Bellman inequality and the terminal term, justifies the candidate value under the stated feasible paths.

#### Problem DO-66: explicit growth solution, case 2

Maximize $`\sum_{t\ge0}\beta^t\log c_t`$ subject to $`c_t+k_{t+1}=Ak_t^\alpha`$, positive consumption and capital, with $`\alpha=\frac{2}{3}`$, $`\beta=\frac{3}{4}`$, and $`A=3`$. Derive the candidate value, policy, and steady state.

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0\lt k'\lt Ak^\alpha`$ and tends to minus infinity at both endpoints.

The first-order condition is

```math
-\frac1{Ak^\alpha-k'}+\frac{\beta B}{k'}=0.
```

Matching the coefficient of $`\log k`$ after substitution gives $`B=\alpha(1+\beta B)`$, so $`B=\alpha/(1-\alpha\beta)=\frac{4}{3}`$. Consequently

```math
k'=\alpha\beta Ak^\alpha=\frac{3}{2}k^{\frac{2}{3}},
\qquad
c=(1-\alpha\beta)Ak^\alpha=\frac{3}{2}k^{\frac{2}{3}}.
```

The constant equation yields

```math
a=\frac{\log(1-\alpha\beta)+[\log A+\alpha\beta\log(\alpha\beta)]/(1-\alpha\beta)}{1-\beta}
\approx 3.243721.
```

The positive steady state solves $`k^*=\alpha\beta A(k^*)^\alpha`$:

```math
k^*=(\alpha\beta A)^{1/(1-\alpha)}=\frac{27}{8},
\qquad c^*\approx 3.375000.
```

For a direct Euler check, the resource policy implies $`k_{t+1}/c_t=\alpha\beta/(1-\alpha\beta)`$. Substituting the next-period consumption rule into $`1/c_t=\beta A\alpha k_{t+1}^{\alpha-1}/c_{t+1}`$ verifies equality.

The discounted terminal marginal-value expression is $`\beta^t k_{t+1}/c_t=\beta^t\alpha\beta/(1-\alpha\beta)\to0`$. This is the relevant transversality check. Since log utility is unbounded below and the state set is not a compact bounded-utility model, a generic sup-norm contraction claim is insufficient. The earlier verification argument, using the Bellman inequality and the terminal term, justifies the candidate value under the stated feasible paths.

#### Problem DO-67: explicit growth solution, case 3

Maximize $`\sum_{t\ge0}\beta^t\log c_t`$ subject to $`c_t+k_{t+1}=Ak_t^\alpha`$, positive consumption and capital, with $`\alpha=\frac{3}{4}`$, $`\beta=\frac{2}{3}`$, and $`A=2`$. Derive the candidate value, policy, and steady state.

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0\lt k'\lt Ak^\alpha`$ and tends to minus infinity at both endpoints.

The first-order condition is

```math
-\frac1{Ak^\alpha-k'}+\frac{\beta B}{k'}=0.
```

Matching the coefficient of $`\log k`$ after substitution gives $`B=\alpha(1+\beta B)`$, so $`B=\alpha/(1-\alpha\beta)=\frac{3}{2}`$. Consequently

```math
k'=\alpha\beta Ak^\alpha=1k^{\frac{3}{4}},
\qquad
c=(1-\alpha\beta)Ak^\alpha=1k^{\frac{3}{4}}.
```

The constant equation yields

```math
a=\frac{\log(1-\alpha\beta)+[\log A+\alpha\beta\log(\alpha\beta)]/(1-\alpha\beta)}{1-\beta}
\approx 0.000000.
```

The positive steady state solves $`k^*=\alpha\beta A(k^*)^\alpha`$:

```math
k^*=(\alpha\beta A)^{1/(1-\alpha)}=1,
\qquad c^*\approx 1.000000.
```

For a direct Euler check, the resource policy implies $`k_{t+1}/c_t=\alpha\beta/(1-\alpha\beta)`$. Substituting the next-period consumption rule into $`1/c_t=\beta A\alpha k_{t+1}^{\alpha-1}/c_{t+1}`$ verifies equality.

The discounted terminal marginal-value expression is $`\beta^t k_{t+1}/c_t=\beta^t\alpha\beta/(1-\alpha\beta)\to0`$. This is the relevant transversality check. Since log utility is unbounded below and the state set is not a compact bounded-utility model, a generic sup-norm contraction claim is insufficient. The earlier verification argument, using the Bellman inequality and the terminal term, justifies the candidate value under the stated feasible paths.

#### Problem DO-68: explicit growth solution, case 4

Maximize $`\sum_{t\ge0}\beta^t\log c_t`$ subject to $`c_t+k_{t+1}=Ak_t^\alpha`$, positive consumption and capital, with $`\alpha=\frac{1}{2}`$, $`\beta=\frac{9}{10}`$, and $`A=1`$. Derive the candidate value, policy, and steady state.

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0\lt k'\lt Ak^\alpha`$ and tends to minus infinity at both endpoints.

The first-order condition is

```math
-\frac1{Ak^\alpha-k'}+\frac{\beta B}{k'}=0.
```

Matching the coefficient of $`\log k`$ after substitution gives $`B=\alpha(1+\beta B)`$, so $`B=\alpha/(1-\alpha\beta)=\frac{10}{11}`$. Consequently

```math
k'=\alpha\beta Ak^\alpha=\frac{9}{20}k^{\frac{1}{2}},
\qquad
c=(1-\alpha\beta)Ak^\alpha=\frac{11}{20}k^{\frac{1}{2}}.
```

The constant equation yields

```math
a=\frac{\log(1-\alpha\beta)+[\log A+\alpha\beta\log(\alpha\beta)]/(1-\alpha\beta)}{1-\beta}
\approx -12.511615.
```

The positive steady state solves $`k^*=\alpha\beta A(k^*)^\alpha`$:

```math
k^*=(\alpha\beta A)^{1/(1-\alpha)}=\frac{81}{400},
\qquad c^*\approx 0.247500.
```

For a direct Euler check, the resource policy implies $`k_{t+1}/c_t=\alpha\beta/(1-\alpha\beta)`$. Substituting the next-period consumption rule into $`1/c_t=\beta A\alpha k_{t+1}^{\alpha-1}/c_{t+1}`$ verifies equality.

The discounted terminal marginal-value expression is $`\beta^t k_{t+1}/c_t=\beta^t\alpha\beta/(1-\alpha\beta)\to0`$. This is the relevant transversality check. Since log utility is unbounded below and the state set is not a compact bounded-utility model, a generic sup-norm contraction claim is insufficient. The earlier verification argument, using the Bellman inequality and the terminal term, justifies the candidate value under the stated feasible paths.

#### Problem DO-69: explicit growth solution, case 5

Maximize $`\sum_{t\ge0}\beta^t\log c_t`$ subject to $`c_t+k_{t+1}=Ak_t^\alpha`$, positive consumption and capital, with $`\alpha=\frac{2}{3}`$, $`\beta=\frac{9}{10}`$, and $`A=2`$. Derive the candidate value, policy, and steady state.

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0\lt k'\lt Ak^\alpha`$ and tends to minus infinity at both endpoints.

The first-order condition is

```math
-\frac1{Ak^\alpha-k'}+\frac{\beta B}{k'}=0.
```

Matching the coefficient of $`\log k`$ after substitution gives $`B=\alpha(1+\beta B)`$, so $`B=\alpha/(1-\alpha\beta)=\frac{5}{3}`$. Consequently

```math
k'=\alpha\beta Ak^\alpha=\frac{6}{5}k^{\frac{2}{3}},
\qquad
c=(1-\alpha\beta)Ak^\alpha=\frac{4}{5}k^{\frac{2}{3}}.
```

The constant equation yields

```math
a=\frac{\log(1-\alpha\beta)+[\log A+\alpha\beta\log(\alpha\beta)]/(1-\alpha\beta)}{1-\beta}
\approx 0.503388.
```

The positive steady state solves $`k^*=\alpha\beta A(k^*)^\alpha`$:

```math
k^*=(\alpha\beta A)^{1/(1-\alpha)}=\frac{216}{125},
\qquad c^*\approx 1.152000.
```

For a direct Euler check, the resource policy implies $`k_{t+1}/c_t=\alpha\beta/(1-\alpha\beta)`$. Substituting the next-period consumption rule into $`1/c_t=\beta A\alpha k_{t+1}^{\alpha-1}/c_{t+1}`$ verifies equality.

The discounted terminal marginal-value expression is $`\beta^t k_{t+1}/c_t=\beta^t\alpha\beta/(1-\alpha\beta)\to0`$. This is the relevant transversality check. Since log utility is unbounded below and the state set is not a compact bounded-utility model, a generic sup-norm contraction claim is insufficient. The earlier verification argument, using the Bellman inequality and the terminal term, justifies the candidate value under the stated feasible paths.

### 7.6 A verified continuous-time quadratic controller

For scalar linear dynamics and a quadratic objective, the HJB equation reduces to an algebraic equation for one coefficient. The correct solution must also satisfy a discounted terminal condition. Solving only the algebraic equation can leave an inadmissible branch.

#### Problem DO-70: scalar discounted LQ control, case 1

Minimize $`\int_0^\infty e^{-\rho t}(qx_t^2+ru_t^2)\,dt`$ with $`\dot x=ax+bu`$, $`x_0=x`$. Use $`\rho=1`$, $`a=0`$, $`b=1`$, $`r=1`$, $`q=2`$.

Solution. For a candidate $`V(x)=Px^2`$, the HJB equation is

```math
\rho Px^2=\min_u\left(qx^2+ru^2+2Px(ax+bu)\right).
```

The expression is strictly convex in $`u`$ because $`r\gt 0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=1`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-1x_t,\qquad
\dot x_t=(-1)x_t,\qquad
x_t=e^{(-1)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=3\gt 0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

The achieved cost can be integrated directly:

```math
J(u^*)=
\frac{q+r(bP/r)^2}{\rho-2(a-b^2P/r)}x^2
=1x^2=V(x).
```

For any admissible control satisfying the corresponding vanishing terminal term, completing the square in the HJB expression yields

```math
J(u)=Px^2+\int_0^\infty e^{-\rho t}
r\left(u_t+\frac{bP}r x_t\right)^2\,dt.
```

The added integral is nonnegative and vanishes for the proposed feedback. This verifies optimality and explains why the chosen Riccati solution is the appropriate one. Admissibility and the terminal term are parts of that verification, not optional additions to an algebraic root calculation.

#### Problem DO-71: scalar discounted LQ control, case 2

Minimize $`\int_0^\infty e^{-\rho t}(qx_t^2+ru_t^2)\,dt`$ with $`\dot x=ax+bu`$, $`x_0=x`$. Use $`\rho=2`$, $`a=1`$, $`b=1`$, $`r=2`$, $`q=2`$.

Solution. For a candidate $`V(x)=Px^2`$, the HJB equation is

```math
\rho Px^2=\min_u\left(qx^2+ru^2+2Px(ax+bu)\right).
```

The expression is strictly convex in $`u`$ because $`r\gt 0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=2`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-1x_t,\qquad
\dot x_t=(0)x_t,\qquad
x_t=e^{(0)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=2\gt 0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

The achieved cost can be integrated directly:

```math
J(u^*)=
\frac{q+r(bP/r)^2}{\rho-2(a-b^2P/r)}x^2
=2x^2=V(x).
```

For any admissible control satisfying the corresponding vanishing terminal term, completing the square in the HJB expression yields

```math
J(u)=Px^2+\int_0^\infty e^{-\rho t}
r\left(u_t+\frac{bP}r x_t\right)^2\,dt.
```

The added integral is nonnegative and vanishes for the proposed feedback. This verifies optimality and explains why the chosen Riccati solution is the appropriate one. Admissibility and the terminal term are parts of that verification, not optional additions to an algebraic root calculation.

#### Problem DO-72: scalar discounted LQ control, case 3

Minimize $`\int_0^\infty e^{-\rho t}(qx_t^2+ru_t^2)\,dt`$ with $`\dot x=ax+bu`$, $`x_0=x`$. Use $`\rho=1`$, $`a=-1`$, $`b=2`$, $`r=1`$, $`q=\frac{5}{2}`$.

Solution. For a candidate $`V(x)=Px^2`$, the HJB equation is

```math
\rho Px^2=\min_u\left(qx^2+ru^2+2Px(ax+bu)\right).
```

The expression is strictly convex in $`u`$ because $`r\gt 0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=\frac{1}{2}`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-1x_t,\qquad
\dot x_t=(-3)x_t,\qquad
x_t=e^{(-3)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=7\gt 0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

The achieved cost can be integrated directly:

```math
J(u^*)=
\frac{q+r(bP/r)^2}{\rho-2(a-b^2P/r)}x^2
=\frac{1}{2}x^2=V(x).
```

For any admissible control satisfying the corresponding vanishing terminal term, completing the square in the HJB expression yields

```math
J(u)=Px^2+\int_0^\infty e^{-\rho t}
r\left(u_t+\frac{bP}r x_t\right)^2\,dt.
```

The added integral is nonnegative and vanishes for the proposed feedback. This verifies optimality and explains why the chosen Riccati solution is the appropriate one. Admissibility and the terminal term are parts of that verification, not optional additions to an algebraic root calculation.

#### Problem DO-73: scalar discounted LQ control, case 4

Minimize $`\int_0^\infty e^{-\rho t}(qx_t^2+ru_t^2)\,dt`$ with $`\dot x=ax+bu`$, $`x_0=x`$. Use $`\rho=3`$, $`a=1`$, $`b=1`$, $`r=1`$, $`q=6`$.

Solution. For a candidate $`V(x)=Px^2`$, the HJB equation is

```math
\rho Px^2=\min_u\left(qx^2+ru^2+2Px(ax+bu)\right).
```

The expression is strictly convex in $`u`$ because $`r\gt 0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=2`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-2x_t,\qquad
\dot x_t=(-1)x_t,\qquad
x_t=e^{(-1)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=5\gt 0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

The achieved cost can be integrated directly:

```math
J(u^*)=
\frac{q+r(bP/r)^2}{\rho-2(a-b^2P/r)}x^2
=2x^2=V(x).
```

For any admissible control satisfying the corresponding vanishing terminal term, completing the square in the HJB expression yields

```math
J(u)=Px^2+\int_0^\infty e^{-\rho t}
r\left(u_t+\frac{bP}r x_t\right)^2\,dt.
```

The added integral is nonnegative and vanishes for the proposed feedback. This verifies optimality and explains why the chosen Riccati solution is the appropriate one. Admissibility and the terminal term are parts of that verification, not optional additions to an algebraic root calculation.

#### Problem DO-74: scalar discounted LQ control, case 5

Minimize $`\int_0^\infty e^{-\rho t}(qx_t^2+ru_t^2)\,dt`$ with $`\dot x=ax+bu`$, $`x_0=x`$. Use $`\rho=2`$, $`a=0`$, $`b=2`$, $`r=2`$, $`q=4`$.

Solution. For a candidate $`V(x)=Px^2`$, the HJB equation is

```math
\rho Px^2=\min_u\left(qx^2+ru^2+2Px(ax+bu)\right).
```

The expression is strictly convex in $`u`$ because $`r\gt 0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=1`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-1x_t,\qquad
\dot x_t=(-2)x_t,\qquad
x_t=e^{(-2)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=6\gt 0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

The achieved cost can be integrated directly:

```math
J(u^*)=
\frac{q+r(bP/r)^2}{\rho-2(a-b^2P/r)}x^2
=1x^2=V(x).
```

For any admissible control satisfying the corresponding vanishing terminal term, completing the square in the HJB expression yields

```math
J(u)=Px^2+\int_0^\infty e^{-\rho t}
r\left(u_t+\frac{bP}r x_t\right)^2\,dt.
```

The added integral is nonnegative and vanishes for the proposed feedback. This verifies optimality and explains why the chosen Riccati solution is the appropriate one. Admissibility and the terminal term are parts of that verification, not optional additions to an algebraic root calculation.

<!-- further-distinct-problems-DO -->

### 7.7 Further decision and verification problems

These additional problems are constructed for these notes.

#### Problem DO-75: terminal rewards change the policy

There are states zero and one and two decision dates, followed by a terminal payment. Staying pays the current state and leaves it unchanged. Switching pays $`-0.2`$ and changes the state. The discount factor is $`0.9`$ and terminal payments are $`g(0)=5,g(1)=0`$. Solve backward and identify any time-dependent action.

Solution. Let $`V_2=(5,0)`$. With one decision remaining, state zero gives stay value $`0+0.9(5)=4.5`$ and switch value $`-0.2+0.9(0)=-0.2`$. State one gives stay value one and switch value $`-0.2+0.9(5)=4.3`$. Thus $`V_1=(4.5,4.3)`$, with staying in state zero and switching in state one.

At the initial date, state zero gives stay value $`0.9(4.5)=4.05`$ and switch value $`-0.2+0.9(4.3)=3.67`$. State one gives stay value $`1+0.9(4.3)=4.87`$ and switch value $`-0.2+0.9(4.5)=3.85`$. Hence

```math
V_0=(4.05,4.87).
```

Initially the policy stays in both states. At the next date it switches from state one. The reward and transition rules did not change; the time remaining until terminal payment did.

Starting in state one, the optimal path receives one immediately, switches next date at discounted cost $`0.18`$, and receives discounted terminal payment $`0.9^2(5)=4.05`$. The total is $`1-0.18+4.05=4.87`$, checking the recursion by direct enumeration of the chosen path.

A stationary policy on the original two-state space cannot describe this finite-horizon solution. Including the date in the state restores a stationary description on the enlarged state space.

#### Problem DO-76: sensitivity of a discounted value to model errors

For two finite-state discounted control models with the same actions and discount $`\beta`$, suppose rewards differ by at most $`\epsilon`$ and transition rows differ in summed absolute probabilities by at most $`\delta`$. Rewards of the second model are bounded in magnitude by $`C`$. Bound the difference of optimal values.

Solution. Write their Bellman operators $`T,\widetilde T`$ and values $`V,\widetilde V`$. For any bounded $`v`$, a single action's reward-plus-continuation difference is at most $`\epsilon+\beta\delta\Vert v\Vert _\infty`$. The maximum inequality gives the same bound for the operators.

Since $`\Vert \widetilde V\Vert _\infty\le C/(1-\beta)`$,

```math
\Vert V-\widetilde V\Vert _\infty
\le\Vert TV-T\widetilde V\Vert _\infty
+\Vert T\widetilde V-\widetilde T\widetilde V\Vert _\infty
\le\beta\Vert V-\widetilde V\Vert _\infty+
\epsilon+\frac{\beta\delta C}{1-\beta}.
```

Rearranging yields

```math
\Vert V-\widetilde V\Vert _\infty
\le\frac{\epsilon}{1-\beta}
+\frac{\beta\delta C}{(1-\beta)^2}.
```

The transition discrepancy was defined as the full sum of absolute row differences. If total variation is defined with a factor one half, its numerical bound must be converted before substitution.

For identical transitions, every fixed policy's value changes by at most $`\epsilon/(1-\beta)`$. Let a policy be optimal in the second model. Comparing its true value, its second-model value, and the true optimal value gives a true policy loss at most $`2\epsilon/(1-\beta)`$.

For $`\epsilon=0.01,\beta=0.9`$ this last bound is $`0.2`$. These bounds show why near-unit discounting magnifies model errors as well as iteration errors. They do not identify the direction in which a particular optimal action will change.

#### Problem DO-77: a value supremum without an optimal policy

There is one state, actions $`a\in(0,1)`$, reward $`a`$, and a deterministic return to the same state. Let $`0\lt \beta\lt 1`$. Find the value supremum and determine whether any admissible policy attains it.

Solution. Every reward is below one, so every discounted payoff is at most $`1/(1-\beta)`$. A stationary action $`a=1-\eta`$ gives payoff $`(1-\eta)/(1-\beta)`$, which approaches that bound as $`\eta\downarrow0`$. The value supremum is therefore $`V^*=1/(1-\beta)`$.

No deterministic policy attains it. For any action sequence,

```math
\frac1{1-\beta}-\sum_{t=0}^\infty\beta^ta_t
=\sum_{t=0}^\infty\beta^t(1-a_t)\gt 0.
```

Already the first term is strictly positive. Allowing the actions to approach one later does not recover that initial shortfall. Randomization over actions strictly below one also leaves positive expected first-period shortfall.

The Bellman operator defined using a supremum is $`Tv=1+\beta v`$ and has the fixed point above. It correctly identifies the value supremum but has no maximizing action. A contraction proves a unique value in its function space; it does not replace an action-attainment assumption.

For any requested payoff gap $`\epsilon\gt 0`$, choose $`0\lt \eta\lt \epsilon(1-\beta)`$. The corresponding constant policy is epsilon optimal. Closing the action set to $`[0,1]`$ creates an exact optimal policy, always choosing one. The missing endpoint is the entire source of nonattainment in this example.

#### Problem DO-78: multiple invariant laws and a transient state

For states one, two, and three, take

```math
P=\begin{pmatrix}1&0&0\\[0pt] 0&1&0\\[0pt] p&1-p&0\end{pmatrix}
```

 with $`0\lt p\lt 1`$. Find every invariant distribution and the limiting law from an arbitrary initial distribution.

Solution. States one and two are separate absorbing classes. State three is transient and is left after one step. If $`\pi=\pi P`$, its third coordinate must be zero because the third column is zero.

Conversely, any distribution $`(\lambda,1-\lambda,0)`$ with $`0\le\lambda\le1`$ remains unchanged. These are all invariant distributions.

Starting with $`q=(q_1,q_2,q_3)`$, one transition gives

```math
qP=(q_1+pq_3,\ q_2+(1-p)q_3,\ 0).
```

This distribution is already invariant, so it is the limiting law. The mixture weight depends on the initial law because there is no communication between the two absorbing classes.

If the chain starts in state three, a single realized path enters state one with probability $`p`$ and state two otherwise, then remains there forever. Its time-average state frequencies are either $`(1,0,0)`$ or $`(0,1,0)`$. The ensemble limiting distribution is their probability mixture $`(p,1-p,0)`$.

An invariant law can therefore be a mixture without describing the long-run frequencies of every individual path as the same deterministic vector. Irreducibility, or another condition joining the recurrent structure, is needed for the usual unique ergodic conclusion.

#### Problem DO-79: discounted values and average rewards

For a two-state chain with transition probabilities $`p`$ from zero to one and $`q`$ from one to zero, both positive, derive the discounted value of state rewards $`r_0,r_1`$. Find the limit of $`(1-\beta)V_\beta`$ as $`\beta\uparrow1`$.

Solution. The stationary probabilities are $`\pi_0=q/(p+q)`$ and $`\pi_1=p/(p+q)`$. Let $`\bar r=\pi_0r_0+\pi_1r_1`$ and $`\eta=1-p-q`$.

The two-state transition recursion gives

```math
E_i[r(X_t)]=\bar r+\eta^t(r_i-\bar r).
```

Discounting and summing the two geometric series yields

```math
V_\beta(i)=\frac{\bar r}{1-\beta}
+\frac{r_i-\bar r}{1-\beta\eta}.
```

For every $`\beta\lt 1`$ the series converge, including a periodic chain with $`\eta=-1`$. Since positive $`p,q`$ imply $`\eta\lt 1`$, multiplying by $`1-\beta`$ and taking the limit gives $`\bar r`$ in both states.

For $`p=1/4,q=1/2,r_0=0,r_1=3`$, the stationary average reward is one and $`\eta=1/4`$. At $`\beta=0.8`$, values are $`5-1/0.8=3.75`$ and $`5+2/0.8=7.5`$.

The diverging common term $`\bar r/(1-\beta)`$ represents the growing effective horizon. The state-dependent difference remains bounded as discounting vanishes. This limit concerns discounted reward averages. Convergence of the state distribution at every date can still fail in a periodic chain, while its discounted average has this well-defined limit.

#### Problem DO-80: finite-horizon growth saving shares

In fully depreciating log-utility growth with production $`Ak^\alpha`$, find the coefficient of log capital and the saving share with $`m`$ consumption dates remaining and no terminal capital value.

Solution. With one date remaining, all output is consumed: $`V_1(k)=\log A+\alpha\log k`$. Thus $`b_1=\alpha`$ and the final saving share is zero. This final case is handled directly; it does not require evaluating a future log value at zero capital.

For $`m\ge2`$, suppose $`V_{m-1}=a_{m-1}+b_{m-1}\log k`$. Maximizing the Bellman objective gives saving share $`s_m=\beta b_{m-1}/(1+\beta b_{m-1})`$ and coefficient recursion $`b_m=\alpha(1+\beta b_{m-1})`$.

Solving the finite geometric recursion gives

```math
b_m=\alpha\sum_{j=0}^{m-1}(\alpha\beta)^j
=\frac{\alpha[1-(\alpha\beta)^m]}{1-\alpha\beta}.
```

Substitution simplifies the share to

```math
s_m=\frac{\alpha\beta[1-(\alpha\beta)^{m-1}]}
{1-(\alpha\beta)^m},\qquad m\ge2.
```

These shares increase with the remaining horizon because $`b_{m-1}`$ increases and the map $`x\mapsto\beta x/(1+\beta x)`$ is increasing. They approach $`\alpha\beta`$, the infinite-horizon share.

For $`\alpha=1/2,\beta=4/5`$, the one-, two-, and three-date shares are $`0,2/7,14/39`$. The limiting share is $`2/5`$. A finite-horizon agent saves less because capital left after the final consumption date has no value. The policy difference comes from the terminal condition and persists even though preferences and production are unchanged.

#### Problem DO-81: log growth with bounded iid productivity

Let current output be $`Az k^\alpha`$ with observed productivity $`z`$. Future productivities are iid and bounded above and away from zero. With $`0\lt \alpha,\beta\lt 1`$ and log consumption utility, derive a value of the form $`a+B\log k+C\log z`$ and its saving policy.

Solution. The Bellman equation is

```math
V(k,z)=\max_{0\lt k'\lt Azk^\alpha}
\lbrace \log(Azk^\alpha-k')+\beta E[V(k',Z')]\rbrace .
```

Put $`m=E[\log Z]`$, which is finite under the stated bounds. The candidate gives expected continuation $`a+B\log k'+Cm`$. Stationarity implies a saving fraction $`s=\beta B/(1+\beta B)`$.

After substitution, matching the log-capital and log-productivity coefficients yields

```math
B=\alpha(1+\beta B)=\frac{\alpha}{1-\alpha\beta},
\qquad C=1+\beta B=\frac1{1-\alpha\beta}.
```

Thus $`s=\alpha\beta`$ and the policy is $`k'=\alpha\beta Azk^\alpha`$. Matching constants gives

```math
a=\frac1{1-\beta}\left[
\log(1-\alpha\beta)+
\frac{\log A+\alpha\beta\log(\alpha\beta)+\beta m}{1-\alpha\beta}
\right].
```

The maximand is strictly concave in positive next capital and diverges negatively at both boundaries. The candidate action is therefore its unique optimum.

For verification, bounded productivity supplies a uniform upper bound on feasible capital from each initial state. Finite-horizon log-value recursions provide the upper bound as in DO-24. Under the candidate policy, log capital satisfies a stable linear recursion with bounded innovations, so it stays bounded along each path and the discounted terminal candidate tends to zero. Telescoping the Bellman equalities then attains the candidate.

The log-capital recursion is $`\ell_{t+1}=\log(\alpha\beta A)+\alpha\ell_t+\log Z_t`$. Its stationary mean is $`[\log(\alpha\beta A)+m]/(1-\alpha)`$ and variance is $`\mathop{\mathrm{Var}} (\log Z)/(1-\alpha^2)`$. Exponentiating the mean log capital does not give mean capital in general.

#### Problem DO-82: a finite-horizon stochastic quadratic control

Minimize $`E[\int_t^T u_s^2/2\,ds+gX_T^2/2]`$ for $`dX_s=u_sds+\sigma dW_s`$, $`X_t=x`$, and $`g\gt 0`$. Derive the value and optimal feedback, including the term caused by noise.

Solution. The HJB equation is $`V_t+\min_u(u^2/2+uV_x+\sigma^2V_{xx}/2)=0`$ with terminal value $`gx^2/2`$. Take $`V(t,x)=P(t)x^2/2+C(t)`$.

Stationarity in the control gives $`u^*=-P(t)x`$. Substitution and coefficient matching give

```math
P'=P^2,\quad P(T)=g,\qquad
C'=-\frac{\sigma^2}{2}P,\quad C(T)=0.
```

Solving backward yields

```math
P(t)=\frac g{1+g(T-t)},\qquad
C(t)=\frac{\sigma^2}{2}\log[1+g(T-t)].
```

Thus the feedback is $`u_s=-gX_s/[1+g(T-s)]`$. The value includes a positive noise term even when the initial state is zero.

For a verification, complete the square in the HJB expression. Itô's formula, expectation, and the terminal condition give for admissible square-integrable controls

```math
E\left[\int_t^T\frac{u_s^2}{2}ds+\frac g2X_T^2\right]
=V(t,x)+\frac12E\int_t^T(u_s+P(s)X_s)^2ds.
```

Finite-horizon moment bounds justify the stochastic-integral expectation. The final term is nonnegative and vanishes for the proposed feedback, proving optimality.

As $`g`$ grows, the noise contribution grows logarithmically. Enforcing an increasingly tight terminal target in the presence of additive Brownian noise has a growing control cost. Taking a singular hard-terminal limit requires its own admissibility analysis; it is not supplied by the finite positive-$`g`$ formula alone.

[Contents](README.md) · [Part 1](05-dynamic-optimization-part-1.md)
