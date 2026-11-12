# Dynamic Optimization

A dynamic problem links today's choice to tomorrow's available choices. The state records the information needed for that link. The action is what can be chosen now. A value function records the best achievable future objective from a specified state, including the current period under the convention used here.

Sources:

- Guido Lorenzoni, [MIT 14.451, Dynamic Optimization Methods with Applications, Fall 2009](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/), lecture notes and Problem Sets 1–4.
- Thomas Sargent and John Stachurski, [Dynamic Programming](https://dp.quantecon.org/), especially the volume on finite states, fixed points, Markov chains, and decision processes.

Discount factors satisfy $`0<\beta<1`$ unless a finite-horizon exception is stated. A stochastic transition is conditional on the current state and chosen action. Examples are constructed here; source problems and subparts are identified individually.

Contents:

1. [Sequential decision problems](#1-sequential-decision-problems)
2. [Bellman equations and contractions](#2-bellman-equations-and-contractions)
3. [Policy evaluation and improvement](#3-policy-evaluation-and-improvement)
4. [Optimal stopping](#4-optimal-stopping)
5. [Consumption, saving, and growth](#5-consumption-saving-and-growth)
6. [Stochastic decisions and continuous time](#6-stochastic-decisions-and-continuous-time)
7. [Further worked problems](#7-further-worked-problems)

## 1. Sequential decision problems

### 1.1 Theory

At date $`t`$, observe state $`s_t`$, choose a feasible action $`a_t\in A_t(s_t)`$, receive reward $`r_t(s_t,a_t)`$, and move to the next state. In a deterministic model the transition is $`s_{t+1}=f_t(s_t,a_t)`$. In a stochastic model it is a conditional probability law $`P_t(ds'\mid s_t,a_t)`$.

A finite-horizon objective with terminal reward $`g`$ is

```math
\max E\left[\sum_{t=0}^{T-1}\beta^t r_t(s_t,a_t)+\beta^Tg(s_T)\right].
```

A policy describes a decision rule. A history-dependent policy may use every past observation. If the stated state summarizes everything affecting rewards and transitions, an optimal Markov policy can be chosen in standard finite-state finite-horizon problems. A state that omits relevant information can make a seemingly Markov decision rule inadequate.

Backward induction starts from the terminal condition and solves one date at a time in reverse order. It does not choose today's action by today's reward alone. Each candidate action is evaluated using both its immediate reward and the optimized continuation from the state it creates.

For continuously varying actions, existence of a best action needs conditions such as a nonempty compact feasible set and continuous objective. If only a supremum exists, a maximizing action may be absent. A Bellman equation written with a maximum carries an attainment claim.

### 1.2 Derivations and proofs

#### Bellman recursion from the sequence problem

Define $`V_T(s)=g(s)`$. At date $`t`$, choose the first action and then an optimal continuation. Conditional on the resulting next state,

```math
V_t(s)=\max_{a\in A_t(s)}
\left\{r_t(s,a)+\beta E[V_{t+1}(s_{t+1})\mid s,a]\right\}.
```

For the upper bound, any full feasible plan has a first action whose later reward cannot exceed $`V_{t+1}`$. Its total is bounded by the displayed maximum. For the lower bound, choose an action attaining that maximum and use an optimal continuation at each next state. In finite state and action spaces this construction requires only finitely many comparisons. In general spaces, measurable selection and existence of continuation policies need separate assumptions.

#### One-period deviations and finite horizons

If a candidate policy chooses a maximizing Bellman action at every date and state, induction proves it optimal. At the last decision date, the continuation is already fixed by $`g`$. Assuming optimality from date $`t+1`$, the chosen action at date $`t`$ maximizes immediate reward plus exactly that best future value. Proceed backward to date zero.

#### The state summarizes the relevant history

A sequential problem begins with a history of observations and actions. A state is a smaller record sufficient to determine current feasible actions, current rewards, and the conditional law of the next state. Two histories assigned the same state must lead to the same continuation problem.

For saving, current wealth alone may be insufficient when labor income has persistent shocks; the current income state may also be needed. For a work-rest problem, the length of the current work streak matters because it changes the reward from another day of work. Omitting that variable creates a recursion that forgets relevant information.

Timing must be stated as carefully as the variables. A standard convention is: observe $`s_t`$, choose $`a_t`$, receive $`r_t(s_t,a_t)`$, then draw $`s_{t+1}`$ from $`P_t(\cdot\mid s_t,a_t)`$. If the next shock is observed before choosing the action, the expectation and maximization appear in a different order.

For a finite state and action problem, backward induction gives

```math
V_T(s)=g(s),\qquad
V_t(s)=\max_{a\in A_t(s)}
\left[r_t(s,a)+\beta\sum_{s'}P_t(s'\mid s,a)V_{t+1}(s')\right].
```

To prove the recursion, any admissible continuation after choosing $`a`$ has value no larger than $`V_{t+1}`$ in each possible next state. This supplies the upper bound. Conversely, choose a maximizing current action and then a maximizing continuation in each next state. Finite sets ensure these choices exist and achieve the displayed bound. This proves the equality by induction.

For noncompact action sets, a supremum may fail to be attained. Writing a maximizing policy then requires additional assumptions, rather than a change of notation from supremum to maximum.

#### Example DO-E13: integer resources and a finite horizon

There are two consumption dates, no interest, and initial inventory two. Consumption must be an integer between zero and remaining inventory. Utility is $`\sqrt c`$ and $`\beta=1/2`$. There is no terminal salvage value.

Solution. At the final date all remaining inventory is consumed because utility is increasing. The one-date value is $`V_1(x)=\sqrt x`$ for $`x=0,1,2`$. At the initial date,

```math
V_0(2)=\max_{c\in\{0,1,2\}}
\left(\sqrt c+\tfrac12\sqrt{2-c}\right).
```

The three values are $`\sqrt2/2`$, $`3/2`$, and $`\sqrt2`$. The middle action is optimal, so consumption is one at each date and total discounted utility is $`3/2`$.

The terminal value is zero after the last consumption date. Giving terminal inventory a positive salvage value would change the final-date decision. The integer action set also matters: a continuous resource allocation would generally allow a different split.

### 1.3 Worked examples

#### Example DO-E01: consume over two dates

An agent has two units of a nonstorable-income-free resource, chooses $`c_0,c_1>0`$ with $`c_0+c_1=2`$, and maximizes $`\log c_0+\tfrac12\log c_1`$.

Solution. At the final date, $`V_1(s)=\log s`$, since all remaining resources are consumed. The first Bellman problem is $`\max_{0<c<2}\{\log c+\tfrac12\log(2-c)\}`$. Its derivative $`1/c-1/[2(2-c)]`$ vanishes at $`c_0=4/3`$, leaving $`c_1=2/3`$. The second derivative is negative throughout, and the objective tends to $`-\infty`$ at either endpoint, proving the unique global optimum.

#### Example DO-E02: a replacement decision

A machine is good or worn. At the final date, operating yields rewards three and one respectively. One date earlier, a good machine can operate for three and become worn; replacement costs two and leaves it good next date. Let $`\beta=0.9`$.

Solution. The final values are $`V_1(G)=3`$, $`V_1(W)=1`$. For a good machine, operating gives $`3+0.9(1)=3.9`$, while replacing without current operating revenue gives $`-2+0.9(3)=0.7`$. Operating is optimal. If replacement also allowed current operation, its current reward must be changed in the model; silently combining those timing conventions would solve a different problem.

### 1.4 Selected exercises

For the next three problems, use the source's finite saving model. At the start of a date, assets are $`a`$, gross return is $`R=1+r=1/\beta`$, income is $`y`$, and next assets obey $`a'=Ra+y-c`$. There are $`m`$ consumption dates remaining, terminal assets are zero, utility is strictly increasing and strictly concave, and borrowing is allowed subject to this terminal feasibility and utility's domain.

#### Problem DO-01: the last consumption date

Source: [MIT 14.451, Problem Set 1, Section 2, Question 1](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/6830cccb461b78379fca41f342bb6f0d_MIT14_451F09_pset1.pdf).

Problem. In the saving model of this section with $`R=1/\beta`$ and no terminal wealth value, solve the problem with one consumption date remaining.

Solution. With $`m=1`$, terminal feasibility sets $`a'=0`$, so $`c=Ra+y`$ and $`V_1(a)=u(Ra+y)`$. Consuming less leaves unused terminal resources and lowers utility; consuming more violates the terminal budget. The proposed general formula $`c_m=y+(1-\beta)Ra/(1-\beta^m)`$ reduces to the same result at $`m=1`$, and its value multiplier $`(1-\beta^m)/(1-\beta)`$ equals one.

#### Problem DO-02: induction for finite-horizon consumption

Source: [MIT 14.451, Problem Set 1, Section 2, Question 2](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/6830cccb461b78379fca41f342bb6f0d_MIT14_451F09_pset1.pdf).

Problem. For the same strictly concave saving problem, use induction to derive the value and common consumption level with $`m`$ dates remaining.

Solution. Suppose the formula holds with $`m-1`$ dates left. In the Bellman objective $`u(c)+\beta V_{m-1}(a')`$, differentiate using $`a'=Ra+y-c`$. The induction formula gives $`V_{m-1}'(a')=R\,u'(c_{m-1}(a'))`$ because its annuity factor cancels the derivative of its consumption rule. Since $`\beta R=1`$, stationarity requires $`u'(c)=u'(c_{m-1}(a'))`$. Strict concavity makes $`u'`$ injective, so current and future optimal consumption agree.

Enforce the intertemporal budget with common consumption $`c`$:

```math
\sum_{j=0}^{m-1}\beta^jc
=Ra+y\sum_{j=0}^{m-1}\beta^j.
```

Thus $`c=y+(1-\beta)Ra/(1-\beta^m)`$. The value is $`V_m(a)=[(1-\beta^m)/(1-\beta)]u(c)`$. Strict concavity of the feasible consumption objective proves global optimality and uniqueness; feasibility and the utility domain must permit the stated common level.

#### Problem DO-03: the infinite-horizon limit

Source: [MIT 14.451, Problem Set 1, Section 2, Question 3](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/6830cccb461b78379fca41f342bb6f0d_MIT14_451F09_pset1.pdf).

Problem. Take the finite-horizon saving solution to the infinite-horizon limit and identify the consumption and asset rules under a solvency restriction.

Solution. As $`m\to\infty`$, $`\beta^m\to0`$, so consumption tends to $`c_\infty=y+(1-\beta)Ra=y+(R-1)a`$. The value tends to $`u(c_\infty)/(1-\beta)`$ whenever this finite utility is defined. The asset rule becomes $`a'=Ra+y-c_\infty=a`$: the agent consumes income plus interest and preserves principal. This infinite-horizon policy satisfies discounted terminal-wealth decay $`\beta^Ta\to0`$. Without a no-Ponzi or comparable solvency restriction, the infinite-horizon budget problem could allow artificial borrowing plans.

### 1.5 Additional solved problems

#### Problem DO-04: changing the discount factor

Constructed problem. In Example DO-E01, replace $`1/2`$ by arbitrary $`\beta>0`$.

Solution. Stationarity is $`1/c_0=\beta/c_1`$, so $`c_1=\beta c_0`$. The resource constraint gives $`c_0=2/(1+\beta)`$ and $`c_1=2\beta/(1+\beta)`$. Greater future weight shifts resources toward the second date. Even $`\beta>1`$ is mathematically admissible here because the horizon is finite; the infinite-horizon contraction argument later requires a discount below one.

#### Problem DO-05: a nonattained Bellman supremum

Constructed problem. A final decision chooses $`a\in(0,1)`$ and earns reward $`a`$.

Solution. The value is $`\sup_{0<a<1}a=1`$. No feasible action attains it. For every $`\varepsilon>0`$, choosing $`a>1-\varepsilon`$ is $`\varepsilon`$-optimal. Writing a maximum and choosing $`a=1`$ would falsely enlarge the feasible set. Compactness of the closed interval would restore attainment.

#### Problem DO-06: an insufficient state

Constructed problem. A reward depends on whether the previous action differed from the current one. Explain why the current external shock alone is insufficient.

Solution. Two histories with the same current shock but different previous actions have different switching costs for the same current choice. Their optimization problems are different. Include the previous action in the state, for example $`s_t=(z_t,a_{t-1})`$. Then the immediate reward and transition depend only on the augmented state and current action, restoring the intended Markov structure.

## 2. Bellman equations and contractions

### 2.1 Theory

For a bounded reward function and discounted infinite horizon, the Bellman operator acts on bounded value functions:

```math
(Tv)(s)=\sup_{a\in A(s)}
\{r(s,a)+\beta E[v(s')\mid s,a]\}.
```

Let $`\mathcal B(S)`$ be the space of bounded real functions on the state set, with norm $`\Vert v\Vert _\infty=\sup_s\vert v(s)\vert `$. On finite or countable states with all functions measurable, this is a complete normed space. On general measurable spaces restrict to bounded measurable functions and require that the operator preserves measurability.

If rewards are uniformly bounded by $`M`$, every discounted policy value has magnitude at most $`M/(1-\beta)`$. The operator maps bounded functions to bounded functions because the continuation term is bounded by $`\beta\Vert v\Vert _\infty`$. With finite actions the supremum is a maximum.

Two properties recur in proofs. Monotonicity means $`v\le w`$ pointwise implies $`Tv\le Tw`$. Discounting means $`T(v+c)=Tv+\beta c`$ for constants $`c`$. Both follow directly from the definition when transition probabilities sum to one.

### 2.2 Derivations and proofs

#### The supremum inequality

For any bounded families $`f_a,g_a`$, $`f_a\le g_a+\sup_b\vert f_b-g_b\vert \le\sup_bg_b+\sup_b\vert f_b-g_b\vert `$. Taking the supremum over $`a`$ and reversing their roles proves

```math
\vert \sup_af_a-\sup_ag_a\vert \le\sup_a\vert f_a-g_a\vert .
```

This does not require the same action to maximize both expressions.

#### Bellman contraction

The immediate rewards cancel when comparing the same action under $`v,w`$. Also $`\vert E[v(s')-w(s')\mid s,a]\vert \le\Vert v-w\Vert _\infty`$. Applying the supremum inequality gives

```math
\Vert Tv-Tw\Vert _\infty\le\beta\Vert v-w\Vert _\infty.
```

The factor below one makes iteration reduce discrepancies globally. A large number of states does not change that factor; a discount near one does make convergence slow.

#### Banach fixed-point theorem with its proof

Start from any bounded $`v_0`$ and set $`v_{n+1}=Tv_n`$. Consecutive differences satisfy $`\Vert v_{n+1}-v_n\Vert \le\beta^n\Vert v_1-v_0\Vert `$. For $`m>n`$, sum them:

```math
\Vert v_m-v_n\Vert \le\frac{\beta^n}{1-\beta}\Vert v_1-v_0\Vert .
```

The sequence is Cauchy, so completeness gives a limit $`v^*`$. The contraction inequality makes $`T`$ continuous, yielding $`Tv^*=\lim Tv_n=\lim v_{n+1}=v^*`$. If $`w^*`$ is another fixed point, $`\Vert v^*-w^*\Vert \le\beta\Vert v^*-w^*\Vert `$, forcing equality of the functions.

The error bounds are

```math
\Vert v_n-v^*\Vert \le\beta^n\Vert v_0-v^*\Vert ,\qquad
\Vert v-v^*\Vert \le\frac{\Vert Tv-v\Vert }{1-\beta}.
```

The second follows from $`\Vert v-v^*\Vert \le\Vert v-Tv\Vert +\beta\Vert v-v^*\Vert `$. It turns a computable Bellman residual into a value-error certificate.

#### The space on which the contraction acts

Let $`B(S)`$ be the bounded real functions on a state set, with $`\Vert v\Vert _\infty=\sup_s\vert v(s)\vert `$. This space is complete. If $`v_n`$ is Cauchy in sup norm, then for each state the numbers $`v_n(s)`$ are Cauchy and have a real limit $`v(s)`$. Fixing one late function bounds all the others uniformly, so the limit is bounded. Passing to the pointwise limit in $`\vert v_n(s)-v_m(s)\vert <\epsilon`$ and then taking the supremum proves uniform convergence to $`v`$.

With $`\vert r(s,a)\vert \le C`$ and $`0<\beta<1`$, the Bellman operator maps bounded functions to bounded functions because $`\vert Tv(s)\vert \le C+\beta\Vert v\Vert _\infty`$. For any two functions,

```math
\left\vert \sup_a F(a)-\sup_a G(a)\right\vert
\le\sup_a\vert F(a)-G(a)\vert .
```

Apply this with the two continuation values. A probability-weighted average of their difference has magnitude at most their sup-norm difference, so $`\Vert Tv-Tw\Vert _\infty\le\beta\Vert v-w\Vert _\infty`$.

The fixed point is consequently unique within this bounded-function space, and $`\Vert V^*\Vert _\infty\le C/(1-\beta)`$. A model with unbounded rewards or an unbounded value needs a different argument or a weighted norm. Naming the contraction theorem without checking its space does not establish its conclusion.

#### Residuals, value error, and policy error

Suppose an approximate value $`v`$ has Bellman residual $`\epsilon=\Vert Tv-v\Vert _\infty`$. Adding and subtracting $`Tv`$ gives

```math
\Vert v-V^*\Vert _\infty
\le\epsilon+\beta\Vert v-V^*\Vert _\infty,
\qquad
\Vert v-V^*\Vert _\infty\le\frac\epsilon{1-\beta}.
```

This is an a posteriori bound: it can be checked from the current approximation without knowing the exact value.

Let $`\pi`$ be greedy for $`v`$, so $`T_\pi v=Tv`$. Put $`\delta=\Vert v-V^*\Vert _\infty`$. Both $`TV^*`$ and $`T_\pi V^*`$ differ from their corresponding evaluations at $`v`$ by at most $`\beta\delta`$. Thus $`\Vert TV^*-T_\pi V^*\Vert _\infty\le2\beta\delta`$. Comparing the fixed point of $`T_\pi`$ with $`V^*`$ gives

```math
\Vert V^*-V_\pi\Vert _\infty
\le\frac{2\beta}{1-\beta}\delta
\le\frac{2\beta}{(1-\beta)^2}\epsilon.
```

These bounds can be conservative, but their assumptions and units are explicit.

#### Example DO-E14: interpreting a stopping tolerance

If $`\beta=0.95`$ and the sup-norm Bellman residual is at most $`0.005`$, the value error bound is $`0.005/0.05=0.1`$. The displayed greedy-policy bound is $`2(0.95)(0.005)/(0.05)^2=3.8`$. A small residual need not imply an equally small policy-loss certificate when discounting is weak.

If the task requires a value error of at most $`0.01`$, a sufficient residual tolerance is $`(1-\beta)0.01=0.0005`$. The requested accuracy should determine the stopping tolerance, rather than treating the last iterate difference as an error in the value itself.

### 2.3 Worked examples

#### Example DO-E03: one state and a constant reward

If the only available action pays three each date and $`\beta=0.8`$, then $`Tv=3+0.8v`$. The fixed point is fifteen. Starting from zero gives $`v_n=15(1-0.8^n)`$. The exact error is $`15(0.8^n)`$, showing directly how the discount controls iteration speed.

#### Example DO-E04: a residual stopping rule

At $`\beta=0.95`$, a candidate has residual at most $`10^{-4}`$ in every state. The value-error bound is $`10^{-4}/0.05=0.002`$. To guarantee error at most $`10^{-4}`$, require residual at most $`5\times10^{-6}`$. A numerical residual and a value error are different quantities.

### 2.4 Selected exercises

#### Problem DO-07: working and resting as a recursive problem

Source: [MIT 14.451, Problem Set 2, Section 2, Question 1](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/47c8930ae9beae1093f1c49fd89c2886_MIT14_451F09_pset2.pdf).

Problem. Formulate the work-rest problem recursively using the number of consecutive work dates as the state.

Solution. Let $`n`$ be the number of consecutive work dates since the last rest. Working gives $`r(n)`$ and moves to $`n+1`$; resting gives zero and resets to zero. Hence

```math
V(n)=\max\{r(n)+\beta V(n+1),\ \beta V(0)\}.
```

The previous work streak is sufficient to determine both current reward and next state. The initial state is zero. A sequence description chooses work/rest decisions over all dates; the recursive form compresses their effect into the current streak.

#### Problem DO-08: bounded dynamic programming on countable states

Source: [MIT 14.451, Problem Set 2, Section 2, Question 2](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/47c8930ae9beae1093f1c49fd89c2886_MIT14_451F09_pset2.pdf).

Problem. Under the stated bounded-reward assumption, prove existence and uniqueness of a bounded work-rest value and existence of an optimal stationary policy.

Solution. Assume the source's additional boundedness condition on $`r(n)`$. The space of bounded sequences indexed by nonnegative integers is complete: a sup-norm Cauchy sequence has coordinatewise limits, the same uniform bound passes to the limit, and the Cauchy inequality passes to uniform convergence. The two-action operator maps this space to itself and is a contraction with factor $`\beta`$ by the preceding proof.

Its unique fixed point is the value function. For any policy, Bellman dominance bounds its finite-horizon reward by the fixed point plus a bounded discounted terminal term, which vanishes. A maximizing action exists in each state because there are two choices; following such actions makes the Bellman inequalities equalities and attains the fixed point. Continuity in a Euclidean state variable is unnecessary here.

### 2.5 Additional solved problems

#### Problem DO-09: an a priori iteration count

Constructed problem. Rewards satisfy $`\vert r\vert \le2`$, $`\beta=0.9`$, and iteration starts from zero. Bound the iterations needed for value error at most $`0.01`$.

Solution. The true value norm is at most $`2/(1-0.9)=20`$. Thus the error after $`n`$ steps is at most $`20(0.9)^n`$. Require $`(0.9)^n\le0.0005`$, or $`n\ge\log(0.0005)/\log(0.9)\approx72.14`$. Seventy-three iterations suffice under the bound. Actual convergence can be faster, but the guarantee uses only the supplied reward bound and discount.

#### Problem DO-10: reward perturbation

Constructed problem. Two problems have the same feasible actions and transitions but rewards differing uniformly by at most $`\varepsilon`$.

Solution. Their operators satisfy $`\Vert Tv-\widetilde Tv\Vert \le\varepsilon`$. For their fixed points, $`\Vert v^*-\widetilde v^*\Vert \le\beta\Vert v^*-\widetilde v^*\Vert +\varepsilon`$, so the difference is at most $`\varepsilon/(1-\beta)`$. This also follows by bounding each date's reward change and summing the discounted series. The bound controls values; a tiny reward change can still switch a nearly tied optimal action.

#### Problem DO-11: no discount contraction at one

Constructed problem. Let $`Tv=r+v`$ in a one-state problem with $`\beta=1`$.

Solution. If $`r>0`$, the infinite total reward is unbounded and no finite fixed point exists. If $`r=0`$, every scalar is a fixed point, so uniqueness fails. The Lipschitz factor is one rather than below one. Undiscounted stochastic shortest-path or average-reward problems require other assumptions and different arguments.

## 3. Policy evaluation and improvement

### 3.1 Theory

For a finite-state Markov decision process, a stationary policy $`\pi`$ chooses one action in each state. It induces reward vector $`r_\pi`$ and row-stochastic transition matrix $`P_\pi`$. Its value satisfies $`v_\pi=r_\pi+\beta P_\pi v_\pi`$, so

```math
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\sum_{j\ge0}\beta^jP_\pi^j`$ converges in the sup operator norm. Each power has norm one, so the geometric discount makes the series summable.

Policy improvement chooses actions maximizing $`r(s,a)+\beta P_a(s,\cdot)v_\pi`$. Evaluation tells us the value of following one policy forever. Improvement asks whether one action now, followed by that policy, gives a better result. The improvement theorem then connects this one-step comparison to the value of following the improved policy forever.

### 3.2 Derivations and proofs

#### Policy improvement

For a greedy improved policy $`\pi'`$, $`T_{\pi'}v_\pi=Tv_\pi\ge T_\pi v_\pi=v_\pi`$. Monotonicity gives $`T_{\pi'}^nv_\pi\ge v_\pi`$ for every $`n`$. Contraction of $`T_{\pi'}`$ makes these iterates converge to $`v_{\pi'}`$. Therefore $`v_{\pi'}\ge v_\pi`$ pointwise.

If $`Tv_\pi=v_\pi`$, uniqueness of the Bellman fixed point implies $`v_\pi=v^*`$, so the old policy is optimal. With finitely many policies, strict improvements and a consistent rule for ties lead to finite termination of policy iteration.

#### A greedy policy from an approximate value

Let $`\pi`$ be greedy for a bounded approximation $`v`$ and suppose $`\Vert v-v^*\Vert \le\varepsilon`$. Then

```math
\Vert v^*-v_\pi\Vert
\le\Vert Tv^*-Tv\Vert +\Vert T_\pi v-T_\pi v_\pi\Vert
\le\beta\varepsilon+\beta(\varepsilon+\Vert v^*-v_\pi\Vert ).
```

Thus $`\Vert v^*-v_\pi\Vert \le2\beta\varepsilon/(1-\beta)`$. This is a worst-case guarantee. If action value gaps are large enough, an approximate value can already identify the exact optimal policy.

#### Policy evaluation as discounted occupancy

For a fixed finite-state policy, $`V_\pi=(I-\beta P_\pi)^{-1}r_\pi`$. The inverse is

```math
(I-\beta P_\pi)^{-1}=\sum_{t=0}^\infty\beta^tP_\pi^t.
```

Its $`(i,j)`$ entry is the expected discounted number of visits to state $`j`$ from state $`i`$. All entries are nonnegative, and each row sums to $`1/(1-\beta)`$. Multiplying by the reward vector therefore adds rewards over the discounted state visits.

This interpretation proves monotonicity of policy values in rewards and explains why the inverse can propagate a local improvement to other states. It is not necessary to compute the inverse explicitly in a numerical implementation; solving the linear system avoids that extra step.

#### An exact policy-improvement identity

Let $`\pi'`$ be greedy for $`V_\pi`$. Define its one-step advantage $`d=T_{\pi'}V_\pi-V_\pi`$. Since retaining the old action is available, $`d\ge0`$ componentwise. Subtract the policy equations to obtain

```math
(I-\beta P_{\pi'})(V_{\pi'}-V_\pi)=d.
```

The nonnegative inverse gives $`V_{\pi'}-V_\pi\ge0`$. If one component of $`d`$ is strictly positive, value improves at that state and at any state that can reach it with positive probability under $`\pi'`$. Not every state must improve strictly.

In a finite discounted problem, keeping the old action whenever it ties for the maximum prevents pointless switches among equivalent policies. Every genuine change with a strict advantage improves a value component; with finitely many policies, policy iteration eventually reaches a policy satisfying all Bellman inequalities.

#### Example DO-E15: one improvement solves the model

Use $`\beta=1/2`$. Action zero stays in the current state and pays rewards $`(1,0)`$. Action one switches to the other state and pays rewards $`(0,2)`$.

Solution. Starting with action zero in both states gives $`V=(2,0)^\top`$. In state one, staying gives $`1+(1/2)2=2`$ while switching gives zero. In state two, staying gives zero while switching gives $`2+(1/2)2=3`$. The improved policy stays in state one and switches in state two.

Its equations are $`v_1=1+(1/2)v_1`$ and $`v_2=2+(1/2)v_1`$, giving $`v=(2,3)^\top`$. Check greediness again using this new vector. In state one, switching gives $`(1/2)3=1.5<2`$. In state two, staying gives $`(1/2)3=1.5<3`$. The policy remains greedy and is optimal.

Evaluation and improvement used different operations. Evaluation held a policy fixed and solved its linear system. Improvement compared actions using the evaluated vector. Mixing value components from unrelated policies during an evaluation would not solve any one policy equation.

### 3.3 Worked examples

#### Example DO-E05: policy evaluation in two states

Let $`r=(1,2)^\top`$, $`P=\left(\begin{smallmatrix}1/2&1/2\\0&1\end{smallmatrix}\right)`$, and $`\beta=1/2`$.

Solution. The second equation is $`v_2=2+v_2/2`$, giving four. The first is $`v_1=1+v_1/4+v_2/4`$, hence $`3v_1/4=2`$ and $`v_1=8/3`$. The better state is absorbing, so starting in the first state includes the possibility of reaching it later.

#### Example DO-E06: a policy improvement

At state one in the preceding example, a new action pays zero now and moves to state two with certainty. Its one-step value is $`0+(1/2)(4)=2`$, below $`8/3`$, so switching is not an improvement. If it instead pays one now, its one-step value is three, above $`8/3`$, so it improves the policy. Evaluate the changed transition matrix to obtain its full value rather than stopping at the one-step comparison.

### 3.4 Selected exercises

#### Problem DO-12: positive and decreasing work-streak values

Source: [MIT 14.451, Problem Set 2, Section 2, Question 3](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/47c8930ae9beae1093f1c49fd89c2886_MIT14_451F09_pset2.pdf).

Problem. With a decreasing work reward and $`r(0)>0`$, prove that the work-streak value is nonincreasing in the streak and strictly positive.

Solution. Start value iteration from zero. If $`v(n)`$ is nonincreasing and $`r(n)`$ is decreasing, both $`r(n)+\beta v(n+1)`$ and the constant $`\beta v(0)`$ are nonincreasing in $`n`$, so their maximum is nonincreasing. Induction and uniform convergence show $`V(n)`$ is nonincreasing.

The value is strictly positive in every state: rest once to reset, then work one date for $`r(0)>0`$, then use the feasible always-rest policy. That plan gives $`\beta r(0)>0`$ from any state; at state zero, working immediately gives at least $`r(0)>0`$. An all-rest plan also shows nonnegativity before this strict improvement.

#### Problem DO-13: cycles of work and rest

Source: [MIT 14.451, Problem Set 2, Section 2, Question 4](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/47c8930ae9beae1093f1c49fd89c2886_MIT14_451F09_pset2.pdf).

Problem. When sufficiently long work streaks yield negative rewards, show that an optimal work-rest policy has a finite threshold and describe its repeating cycle.

Solution. The advantage of working over resting is $`D(n)=r(n)+\beta V(n+1)-\beta V(0)`$. It decreases with $`n`$, so an optimal stationary policy has a work region followed by a rest region, allowing ties at a threshold. At zero, always resting would give zero, while working gives a positive feasible reward, so at least one work date occurs.

If $`r(n)<0`$ for $`n\ge\bar n`$, then $`V(n+1)\le V(0)`$ gives $`D(n)<0`$ there; resting is strictly better. Thus a threshold $`n^*`$ satisfies $`1\le n^*\le\bar n`$. Starting just after rest, the state visits $`0,1,\ldots,n^*-1`$ while working, then rests and resets. The optimal path repeats a cycle of $`n^*`$ work dates and one rest date.

#### Problem DO-14: patience and the work cycle

Source: [MIT 14.451, Problem Set 2, Section 2, Question 5](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/47c8930ae9beae1093f1c49fd89c2886_MIT14_451F09_pset2.pdf).

Problem. For $`r(n)=3-2n`$, compare the possible optimal work-rest cycles as the discount factor changes.

Solution. With $`r(n)=3-2n`$, only streak rewards three and one are positive. The threshold result leaves cycles of length one or two work dates. Their values from state zero are

```math
V^{(1)}=\frac3{1-\beta^2},\qquad
V^{(2)}=\frac{3+\beta}{1-\beta^3}.
```

The denominators account for the zero-reward rest date. Comparing after multiplying by positive denominators gives a difference numerator $`\beta(1-\beta)(1-2\beta)`$ in favor of the two-work cycle. Hence two work dates are better for $`\beta<1/2`$, one is better for $`\beta>1/2`$, and both are optimal at equality. A patient worker values resetting sooner to regain the high reward three.

### 3.5 Additional solved problems

#### Problem DO-15: a stationary Markov distribution

Constructed problem. For $`P=\left(\begin{smallmatrix}0.8&0.2\\0.4&0.6\end{smallmatrix}\right)`$, solve $`\pi=\pi P`$.

Solution. Write $`\pi=(p,1-p)`$. The first equation is $`p=0.8p+0.4(1-p)`$, so $`0.6p=0.4`$ and $`p=2/3`$. Thus $`\pi=(2/3,1/3)`$. Both states communicate and each has a self-transition, so the finite chain is irreducible and aperiodic: the stationary distribution is unique and time distributions converge to it. Irreducibility supplies uniqueness; aperiodicity supplies convergence from arbitrary initial states.

#### Problem DO-16: an invariant distribution without convergence

Constructed problem. Analyze $`P=\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)`$.

Solution. The unique stationary distribution is $`(1/2,1/2)`$, since invariance requires equal masses. Starting in state one, the distribution alternates forever between $`(1,0)`$ and $`(0,1)`$, so it does not converge. The chain is irreducible but periodic with period two. Its time-average occupation frequencies still approach one half for each state.

#### Problem DO-17: absorbing classes and nonunique invariance

Constructed problem. For $`P=I_2`$, characterize stationary distributions and discounted values.

Solution. Every probability vector is invariant because $`\pi I_2=\pi`$. There are two closed absorbing classes. Nevertheless for any reward vector and $`\beta<1`$, the discounted policy value is uniquely $`r/(1-\beta)`$. Nonuniqueness of an invariant distribution does not contradict uniqueness of a discounted Bellman fixed point; the two equations solve different problems.

## 4. Optimal stopping

### 4.1 Theory

A stopping problem gives a choice between taking a terminal payoff now and continuing. If $`g(s)`$ is the stopping payoff and continuing gives current reward $`r(s)`$, the Bellman equation is

```math
V(s)=\max\{g(s),r(s)+\beta E[V(s')\mid s]\}.
```

The stopping region contains states where stopping is at least as good as continuation. A threshold policy requires additional order or monotonicity structure; the existence of two actions by itself does not guarantee one.

A stopping time is a date whose occurrence can be determined from information available by that date. A rule depending on a future shock is not an admissible stopping decision at the present date.

### 4.2 Derivations and proofs

#### Independent offers and a reservation threshold

Suppose each date offers $`W\sim U(0,1)`$ independently. Accepting pays $`W`$ once and ends the problem. Rejecting costs $`c\ge0`$ now and leads to another offer after discounting. The continuation value $`C=-c+\beta E[V(W)]`$ does not depend on the current offer. Hence $`V(w)=\max(w,C)`$.

If $`0\le C\le1`$, then $`E[\max(W,C)]=\int_0^CC\,dw+\int_C^1w\,dw=(1+C^2)/2`$. Thus

```math
C=-c+\frac\beta2(1+C^2).
```

For $`c=0`$, the root in $`[0,1]`$ is $`C=[1-\sqrt{1-\beta^2}]/\beta`$. The other quadratic root exceeds one and does not satisfy the assumed region. If the cost is sufficiently high that $`C<0`$, every nonnegative offer is accepted, and the correct continuation calculation is instead $`C=-c+\beta/2`$. Always verify the region used in an integral formula.

#### A finite-horizon stopping envelope

Let $`G_t`$ be an integrable adapted payoff already expressed in date-zero discounted units, for dates $`0,\ldots,T`$. Set

```math
Y_T=G_T,\qquad
Y_t=\max(G_t,E[Y_{t+1}\mid\mathcal F_t]).
```

The process $`Y`$ dominates the payoff and is a supermartingale because $`Y_t\ge E[Y_{t+1}\mid\mathcal F_t]`$. It is the smallest supermartingale with this domination property: any other such process dominates $`Y_T`$, and backward induction then makes it dominate every earlier $`Y_t`$.

Define $`\tau^*`$ as the first date at which $`Y_t=G_t`$. It exists no later than $`T`$. Before that date the conditional continuation equality holds, so the stopped envelope is a martingale. Bounded stopping and the terminal equality give $`E[G_{\tau^*}]=E[Y_0]`$. For any other stopping time $`\tau\le T`$, domination and the supermartingale property give $`E[G_\tau]\le E[Y_0]`$. This proves optimality.

If payoffs are expressed in current-date units instead, the continuation term includes the discount factor. The two conventions are equivalent when applied consistently. They should not be mixed within one recursion.

#### When a reservation threshold is justified

For independent identically distributed offers and no other state variable, the value after observing an offer compares that offer with one continuation constant. Acceptance is therefore a threshold rule. If offers are persistent, the continuation value can depend on the current offer or another state. A threshold then needs a monotonicity argument; it is not automatic from the word “stopping.”

Finite-horizon thresholds can change with the number of opportunities remaining. A long remaining horizon makes rejection less costly than rejection on the final date. Terminal conditions are therefore directly reflected in the acceptance rule.

#### Example DO-E16: several remaining uniform offers

Offers are independent uniform draws on $`[0,1]`$. There is no discounting or search cost, and the last offer must be accepted. Find the value before observing an offer when one, two, or three offers remain.

Solution. Finite horizon permits $`\beta=1`$. With one offer left, $`v_1=E[U]=1/2`$. With two left, reject only when the current offer is below $`v_1`$:

```math
v_2=E[\max(U,1/2)]
=\int_0^{1/2}\frac12du+\int_{1/2}^1u\,du
=\frac58.
```

For any continuation constant $`a\in[0,1]`$, $`E[\max(U,a)]=(1+a^2)/2`$. Thus $`v_3=(1+(5/8)^2)/2=89/128`$. The initial acceptance threshold with three offers left is $`5/8`$; with two left it is $`1/2`$; on the final date every offer is accepted.

The increasing values reflect additional options. The infinite-horizon undiscounted limit is a different problem and is not justified by the bounded discounted contraction argument.

### 4.3 Worked examples

#### Example DO-E07: a reservation offer

At $`\beta=0.8`$ and zero rejection cost, $`C=(1-0.6)/0.8=1/2`$. Accept offers at least one half. The acceptance probability each date is one half, so the number of observed offers until acceptance has mean two. The expected accepted offer conditional on acceptance is $`3/4`$, but the time-zero value before seeing an offer is $`E[V(W)]=(1+C^2)/2=5/8`$ because delayed acceptance is discounted.

#### Example DO-E08: finite-horizon stopping

There are two offers left, discount is one, and offers are iid $`U(0,1)`$. The last offer must be accepted, so its pre-observation value is $`1/2`$. At the preceding date, accept if the current offer is at least $`1/2`$. The value before that first offer is $`E\max(W,1/2)=5/8`$. A threshold can depend on remaining horizon; it need not equal an infinite-horizon threshold.

### 4.4 Additional solved problems

#### Problem DO-18: a two-point offer distribution

Constructed problem. Offers are zero or one with equal probabilities, there is no rejection cost, and $`\beta<1`$.

Solution. Reject zero and accept one is a candidate. Before observing an offer its value $`v`$ satisfies $`v=(1/2)(1)+(1/2)\beta v`$, so $`v=1/(2-\beta)`$. At zero, continuation $`\beta/(2-\beta)>0`$ exceeds stopping. At one, stopping one exceeds continuation because $`\beta/(2-\beta)<1`$. Both comparisons verify the candidate policy globally.

#### Problem DO-19: an inadmissible future rule

Constructed problem. A proposed rule stops today if tomorrow's offer will be below today's offer.

Solution. The decision event depends on tomorrow's unseen random variable and cannot generally be determined from today's information. It is not a stopping-time rule under the offer filtration. Its expected payoff can exceed the feasible optimum because it grants foresight. A comparison with admissible Bellman policies would be invalid.

#### Problem DO-20: sensitivity of a reservation threshold

Constructed problem. In the uniform-offer model with an interior threshold, differentiate $`C=-c+\beta(1+C^2)/2`$.

Solution. Holding $`\beta`$ fixed gives $`(1-\beta C)dC/dc=-1`$, so the threshold decreases as rejection becomes more costly. Holding $`c`$ fixed gives $`(1-\beta C)dC/d\beta=(1+C^2)/2>0`$, so patience raises the threshold. The denominator is positive for $`C\le1`$ and $`\beta<1`$, justifying both signs.

#### Problem DO-21: a kink in a value function

Source: [MIT 14.451, Problem Set 2, Section 3, Question 1](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/47c8930ae9beae1093f1c49fd89c2886_MIT14_451F09_pset2.pdf).

Problem. For $`f(x,y)=\max\{-(y+1)^2,\ x-(y-1)^2\}`$ and $`y\in[-2,2]`$, find the maximized value and all optimizing choices and examine continuity and differentiability.

Solution. For $`f(x,y)=\max\{-(y+1)^2,\ x-(y-1)^2\}`$ and $`y\in[-2,2]`$, each branch is bounded above by zero or $`x`$, respectively. These bounds are attained at $`y=-1`$ and $`y=1`$. Hence $`V(x)=\max(0,x)`$. For $`x<0`$, the unique optimizer is $`-1`$; for $`x>0`$, it is $`1`$; at zero both are optimal. The value is not differentiable at zero and no single-valued optimizing selection is continuous there, even though both optimizers are interior.

#### Problem DO-22: unique plans and many steady states

Source: [MIT 14.451, Problem Set 3, Section 1, Questions 2–3](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/3cc76fb7d321d16eb12dfae8a085d4f0_MIT14_451F09_pset3.pdf).

Problem. Use the return $`F(x,y)=-(x-y)^2`$ with nonnegative states to investigate whether a unique optimal plan from each initial state forces a unique steady state.

Solution. With $`F(x,y)=-(x-y)^2`$ and every nonnegative next state feasible, each reward is nonpositive. Keeping $`x_{t+1}=x_t`$ forever achieves zero, the largest possible total. Any path that changes state at a finite date incurs a strictly negative discounted reward and can never offset it with positive rewards, so it is not optimal. The unique optimal path from each $`x_0`$ is constant, and every nonnegative state is a steady state. The objective is not jointly strictly concave: it is flat along the diagonal direction $`(1,1)`$. Unique optimal continuation from a fixed initial condition does not force a unique steady state across all initial conditions.

## 5. Consumption, saving, and growth

### 5.1 Theory

In a saving problem, the budget links marginal utility across dates. At an interior optimum with deterministic gross return $`R_{t+1}`$,

```math
u'(c_t)=\beta R_{t+1}u'(c_{t+1}).
```

This Euler equation says that a small unit of forgone current consumption must have the same utility value as the extra feasible consumption it produces next date. It is a necessary local condition under the stated interiority and differentiability assumptions. It does not by itself impose terminal feasibility or rule out implausible borrowing plans.

A borrowing constraint can change equality to an inequality. If $`a_{t+1}\ge0`$ and a nonnegative multiplier $`\lambda_t`$ is attached to this constraint in the appropriate sign convention, the condition is $`u'(c_t)=\beta E_t[R_{t+1}u'(c_{t+1})]+\lambda_t`$. When the constraint binds, the agent would prefer more current consumption financed by borrowing but cannot do so.

In a growth problem, current output is allocated between consumption and next capital. With full depreciation, $`c_t+k_{t+1}=f(k_t)`$. A steady state is a state preserved by the optimal policy: $`k^*=g(k^*)`$. A steady state can exist without every initial condition converging to it. Stability is a property of the policy dynamics, established separately.

### 5.2 Derivations and proofs

#### The Euler equation as a feasible perturbation

Suppose the budget is $`a_{t+1}=R_ta_t+y_t-c_t`$. Reduce $`c_t`$ by $`\varepsilon`$ and increase $`c_{t+1}`$ by $`R_{t+1}\varepsilon`$, leaving later assets unchanged. If both perturbed consumptions and asset choices remain feasible, the first-order utility change is

```math
-\beta^tu'(c_t)\varepsilon+
\beta^{t+1}u'(c_{t+1})R_{t+1}\varepsilon.
```

At an interior optimum either sign of a small perturbation is possible, so its coefficient must be zero. This proves the Euler equality. With uncertainty, the current perturbation is chosen using current information and the future utility change is conditionally averaged.

#### A transversality condition

For a standard nonnegative-asset infinite-horizon saving problem, a common sufficient terminal condition accompanying concavity, feasibility, and Euler/KKT conditions is

```math
\lim_{T\to\infty}E[\beta^T u'(c_T)a_{T+1}]=0.
```

Its form depends on asset timing and the model. To see its role, use concavity to bound the objective difference between any feasible comparison plan and the candidate by the sum of marginal-utility-weighted consumption differences. Substitute budget differences and use the Euler conditions to telescope that sum. Interior terms cancel, leaving an initial term, zero for a shared initial state, and a discounted terminal asset term. A suitable transversality limit makes the upper bound nonpositive.

This is why Euler equations alone are insufficient: telescoping a finite horizon leaves a terminal term, and infinite-horizon optimality needs control of that term.

#### Borrowing constraints turn Euler equalities into inequalities

Write current resources as $`a+y`$ and let next-period assets satisfy

```math
c+\frac{a'}R=a+y,\qquad a'\ge0.
```

For an interior choice of $`a'`$, differentiating $`u(a+y-a'/R)+\beta E[V(a',z')]`$ gives $`u'(c)=\beta R E[V_a(a',z')]`$. An envelope condition identifies $`V_a`$ with next period's marginal utility under the stated wealth convention.

At the lower bound $`a'=0`$, increasing saving is the only feasible local direction. The derivative of the maximand in that direction must be nonpositive:

```math
u'(c)\ge\beta R E[u'(c')].
```

With multiplier $`\lambda\ge0`$, this can be written $`u'(c)=\beta R E[u'(c')]+R\lambda`$ and $`\lambda a'=0`$. A binding borrowing constraint can make present marginal utility strictly larger than the discounted future marginal utility. Forcing equality would ask the agent to borrow when borrowing is prohibited.

Under uncertainty, the expectation applies to marginal utility, not to consumption before applying utility's derivative. For log utility, $`E[1/c']`$ is generally not $`1/E[c']`$. Jensen's inequality gives the former at least as large when consumption is positive and the expectations exist. This identity is one source of effects from consumption risk, but a general saving comparison also depends on the rest of the model.

#### Terminal wealth and transversality

Finite-horizon resource constraints telescope into a present-value budget plus terminal wealth. When terminal wealth has no utility or bequest value, leaving it unused cannot be optimal for an increasing utility function. In an infinite horizon there is no final date at which to spend it, so an Euler equation alone can admit paths that save too much forever.

A transversality condition rules out a positive limiting discounted marginal value of unused assets. In a deterministic interior consumption-saving convention, a common form is $`\lim_{t\to\infty}\beta^t u'(c_t)a_{t+1}=0`$. The precise timing factor depends on whether assets are measured before or after interest; a fixed positive factor does not change a zero limit.

This is distinct from a no-Ponzi restriction, which prevents financing consumption by debt that is never repaid. One controls excess unconsumed value; the other restricts borrowing paths. A verification proof should state the conditions it actually uses.

#### Example DO-E17: a binding borrowing limit

There are two dates. Initial resources are one, next income is two, interest is zero, utility is $`\log c_0+\tfrac12\log c_1`$, and saving $`a`$ must be nonnegative.

Solution. Feasibility gives $`c_0=1-a`$ and $`c_1=2+a`$, with $`0\le a<1`$. The objective derivative is $`-1/(1-a)+(1/2)/(2+a)`$, which is negative throughout the interval. Thus $`a^*=0`$ and consumption is $`(1,2)`$.

The Euler inequality reads $`1\ge(1/2)(1/2)=1/4`$. Its strictness corresponds to the binding borrowing limit. Without that limit, total resources three would be allocated as $`c_0=2`$ and $`c_1=1`$, requiring $`a=-1`$. The unconstrained Euler solution is therefore infeasible in the stated problem.

### 5.3 Worked examples

#### Example DO-E09: log utility and consumption growth

For $`u(c)=\log c`$ and constant $`R`$, Euler gives $`1/c_t=\beta R/c_{t+1}`$, hence $`c_{t+1}=\beta R c_t`$. If $`\beta R=1`$, consumption is flat; if it exceeds one, interior consumption grows. The budget and solvency conditions determine its initial level.

#### Example DO-E10: a binding borrowing constraint

In a two-date model with $`R=1`$, $`\beta=1`$, incomes $`(1,3)`$, and no initial assets, the unconstrained log-utility solution would consume two each date, requiring borrowing one initially. Under $`a_1\ge0`$, choose $`a_1=0`$, so consumption is $`(1,3)`$. The marginal utility inequality is $`1>1/3`$: the agent would like to shift future resources toward today but the constraint prevents it.

### 5.4 Selected exercises

For the next four problems, consider the source's full-depreciation growth model with $`f(k)=Ak^\alpha`$, $`A>0`$, $`0<\alpha<1`$, log utility, and $`0<\beta<1`$.

#### Problem DO-23: solving the growth Bellman equation

Source: [MIT 14.451, Problem Set 3, Section 2, Question 1](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/3cc76fb7d321d16eb12dfae8a085d4f0_MIT14_451F09_pset3.pdf).

Problem. For fully depreciating Cobb-Douglas growth and log utility, solve the Bellman equation using an affine function of log capital.

Solution. The Bellman equation is $`V(k)=\max_{0<k'<Ak^\alpha}\{\log(Ak^\alpha-k')+\beta V(k')\}`$. Guess $`V(k)=a+b\log k`$. At fixed $`k`$, stationarity gives $`-1/(Ak^\alpha-k')+\beta b/k'=0`$, so $`k'=\beta bAk^\alpha/(1+\beta b)`$.

Matching the coefficient of $`\log k`$ gives $`b=\alpha(1+\beta b)`$, hence $`b=\alpha/(1-\alpha\beta)`$. Thus the saving share is $`\alpha\beta`$:

```math
k'=\alpha\beta Ak^\alpha,\qquad
c=(1-\alpha\beta)Ak^\alpha.
```

Matching constants gives

```math
a=\frac1{1-\beta}\left[
\log(1-\alpha\beta)
+\frac{\log A+\alpha\beta\log(\alpha\beta)}{1-\alpha\beta}
\right].
```

The scalar objective is strictly concave in $`k'`$ and tends to $`-\infty`$ at both endpoints, so the candidate choice uniquely maximizes it.

#### Problem DO-24: verifying the infinite-horizon value

Source: [MIT 14.451, Problem Set 3, Section 2, Question 2](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/3cc76fb7d321d16eb12dfae8a085d4f0_MIT14_451F09_pset3.pdf).

Problem. Verify the log-growth candidate as the actual infinite-horizon value, accounting for its unbounded-below utility rather than invoking a bounded-function contraction.

Solution. The unbounded-below log value is outside the earlier bounded-function contraction theorem, so that theorem cannot verify the guess. First bound all feasible capital from above by $`K=\max(k_0,A^{1/(1-\alpha)})`$: if $`k_t\le K`$, then $`k_{t+1}\le Ak_t^\alpha\le K`$. Therefore every period's utility is at most the finite number $`M=\log(AK^\alpha)`$.

Let $`V_N`$ be the finite-horizon value with $`N+1`$ consumption dates and no terminal capital value. Backward induction gives $`V_N(k)=a_N+b_N\log k`$, with $`b_0=\alpha`$ and $`b_N=\alpha(1+\beta b_{N-1})`$. Hence $`b_N\to b`$. The constant recursion has form $`a_N=\beta a_{N-1}+h(b_{N-1})`$, where $`h`$ is continuous near the positive limit; the discounted recursion implies $`a_N\to a`$ by splitting the weighted sum into a finite old portion and a uniformly small recent coefficient error.

Any infinite feasible plan has total utility at most $`V_N(k_0)+\beta^{N+1}M/(1-\beta)`$, with the convention allowing value $`-\infty`$. Taking limits gives an upper bound $`V(k_0)`$. Under the proposed policy, $`\log k_{t+1}=\log(\alpha\beta A)+\alpha\log k_t`$ converges to a finite constant, so its Bellman terminal term $`\beta^TV(k_T)`$ vanishes. Telescoping Bellman equalities then gives total utility exactly $`V(k_0)`$. This proves optimality and the value formula.

#### Problem DO-25: the growth steady state

Source: [MIT 14.451, Problem Set 3, Section 2, Question 3](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/3cc76fb7d321d16eb12dfae8a085d4f0_MIT14_451F09_pset3.pdf).

Problem. Find the positive steady state of the log-growth policy and determine its stability from every positive initial capital.

Solution. Solve $`k^*=\alpha\beta A(k^*)^\alpha`$ to obtain $`k^*=(\alpha\beta A)^{1/(1-\alpha)}>0`$. Let $`z_t=\log k_t-\log k^*`$. The policy gives $`z_{t+1}=\alpha z_t`$, so $`z_t=\alpha^tz_0\to0`$ for every positive initial capital. Thus the positive steady state is globally stable on $`k>0`$. Also $`g'(k)=\alpha^2\beta Ak^{\alpha-1}`$ and $`g'(k^*)=\alpha`$, using the steady-state identity. Steady consumption is $`(1-\alpha\beta)A(k^*)^\alpha`$.

#### Problem DO-26: the same local policy derivative from Euler

Source: [MIT 14.451, Problem Set 3, Section 2, Question 4](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/3cc76fb7d321d16eb12dfae8a085d4f0_MIT14_451F09_pset3.pdf).

Problem. Differentiate the growth Euler equation around the steady state to recover the local policy derivative and distinguish its stable and unstable roots.

Solution. The Euler equation is $`c_{t+1}=\beta f'(k_{t+1})c_t`$. If $`k_{t+1}=g(k_t)`$ and $`c_t=f(k_t)-g(k_t)`$, substitute to get

```math
f(g(k))-g(g(k))
=\beta f'(g(k))[f(k)-g(k)].
```

Differentiate at $`k^*`$ and write $`d=g'(k^*)`$. Since $`f'(k^*)=1/\beta`$, the result is

```math
(1/\beta-d)d=\beta f''(k^*)d\,c^*+(1/\beta-d).
```

Here $`\beta f''(k^*)c^*=(\alpha-1)(1-\alpha\beta)/(\alpha\beta)`$. Substitution reduces the quadratic to roots $`d=\alpha`$ and $`d=1/(\alpha\beta)`$. The second exceeds one and belongs to the unstable Euler branch. The stable optimal policy has derivative $`\alpha`$, agreeing with the explicit solution. Selecting a root requires stability or the appropriate terminal condition; Euler alone presents both branches.

### 5.5 Additional solved problems

#### Problem DO-27: a numerical growth policy

Constructed problem. Take $`A=2`$, $`\alpha=1/2`$, and $`\beta=0.9`$.

Solution. The saving share is $`0.45`$, policy is $`k'=0.9\sqrt k`$, and steady capital is $`0.9^2=0.81`$. At $`k=1`$, output is two, next capital is $`0.9`$, and consumption is $`1.1`$. At the steady state output is $`1.8`$ and consumption is $`0.99`$. The derivative of the policy at the steady state is one half, so small log deviations halve each date.

#### Problem DO-28: retaining depreciated capital

Constructed problem. With $`k'=Ak^\alpha+(1-\delta)k-c`$, derive the interior steady-state Euler condition.

Solution. A unit saved becomes capital whose next gross marginal return is $`\alpha A(k')^{\alpha-1}+1-\delta`$. At a stationary interior consumption level, Euler gives $`1=\beta[\alpha A(k^*)^{\alpha-1}+1-\delta]`$. Thus $`k^*=[\alpha A/(1/\beta-1+\delta)]^{1/(1-\alpha)}`$ when the denominator is positive. The full-depreciation saving-share policy need not remain valid after changing the technology's resource equation.

#### Problem DO-29: terminal wealth changes the last decision

Constructed problem. At the final date, maximize $`\log c+\eta\log a'`$ with $`c+a'=w`$ and $`\eta>0`$.

Solution. Stationarity gives $`1/c=\eta/a'`$, so $`c=w/(1+\eta)`$ and $`a'=\eta w/(1+\eta)`$. The terminal asset is not zero because it carries a positive terminal payoff. Backward induction must begin from the specified terminal objective, rather than automatically imposing consumption of all resources.

#### Problem DO-30: comparative statics of the steady state

Constructed problem. In the full-depreciation log-growth model, differentiate $`\log k^*`$ with respect to $`\log A`$ and $`\log\beta`$.

Solution. Since $`\log k^*=[\log\alpha+\log\beta+\log A]/(1-\alpha)`$, each elasticity is $`1/(1-\alpha)`$. Productivity and patience both increase steady capital. These are comparisons of different stationary solutions; an unexpected change in a parameter also creates a transition path from the existing capital stock.

## 6. Stochastic decisions and continuous time

### 6.1 Theory

With Markov uncertainty $`z`$, a stochastic Bellman equation includes the conditional distribution of the next shock:

```math
V(x,z)=\max_{a\in A(x,z)}
\left\{r(x,z,a)+\beta\sum_{z'}P(z,z')V(f(x,z,a),z')\right\}.
```

An expectation over unconditional shock probabilities is generally wrong if shocks are persistent. The transition law must also reflect whether today's action affects tomorrow's state or today's reward.

For a finite Markov chain, a stationary distribution solves $`\pi=\pi P`$ and $`\sum_i\pi_i=1`$. Existence follows by compactness of the probability simplex, for example by taking limits of averages of iterated distributions. Irreducibility ensures uniqueness. Aperiodicity, together with irreducibility, gives convergence of time distributions. Discounted value-function uniqueness does not require any of these chain properties.

In continuous time, a discounted control problem can lead to a Hamilton-Jacobi-Bellman equation. For controlled diffusion $`dX=b(X,a)dt+\sigma(X,a)dW`$ and reward $`r`$, its formal stationary equation is

```math
\rho V(x)=\sup_a\{r(x,a)+b(x,a)V'(x)
+\tfrac12\sigma(x,a)^2V''(x)\}.
```

A classical verification theorem needs a sufficiently smooth candidate, an admissible maximizing control, integrability of stochastic terms, and an appropriate discounted terminal condition. Writing the HJB equation supplies a candidate, rather than completing verification.

### 6.2 Derivations and proofs

#### From a short decision interval to HJB

Over a small interval $`h`$, dynamic programming gives approximately $`V(x)=\sup_a\{r(x,a)h+e^{-\rho h}E[V(X_h)]\}`$. For smooth $`V`$, Itô's formula gives $`E[V(X_h)]=V(x)+h[bV'+\tfrac12\sigma^2V'']+o(h)`$ under local integrability. Expand $`e^{-\rho h}=1-\rho h+o(h)`$, subtract $`V(x)`$, divide by $`h`$, and let $`h\downarrow0`$. This yields the HJB expression. Interchanging the optimization and limiting operations requires assumptions; the calculation explains the equation's terms.

#### A discounted linear-quadratic control

Minimize $`\int_0^\infty e^{-\rho t}(qX_t^2+ru_t^2)dt`$ subject to $`\dot X=aX+bu`$, with $`q,r,\rho>0`$ and $`b\ne0`$. Guess $`V(x)=Px^2`$. The HJB minimization gives $`u^*=-bPx/r`$ and

```math
\rho P=q+2aP-\frac{b^2P^2}{r}.
```

The positive solution is

```math
P=\frac{(2a-\rho)+\sqrt{(2a-\rho)^2+4(b^2/r)q}}{2b^2/r}.
```

Under this feedback, $`X_t=x_0e^{(a-b^2P/r)t}`$. The discounted terminal value has exponential rate $`2a-2b^2P/r-\rho=-\sqrt{(2a-\rho)^2+4(b^2/r)q}<0`$, so it vanishes. Completing the HJB square yields

```math
q x^2+ru^2+V'(x)(ax+bu)-\rho V(x)
=r\left(u+\frac{bP}{r}x\right)^2\ge0.
```

Integrate this identity along an admissible comparison path. With the appropriate terminal limit, it proves no control has lower cost, while the feedback makes the square zero. The Riccati root is therefore verified, rather than chosen solely because it is algebraically positive.

#### Expectation and action timing

The stochastic Bellman expression $`\max_a E[r(s,a,Z)+\beta V(s')]`$ assumes the action is chosen before the unresolved shock $`Z`$. If the shock is observed first, the value becomes $`E[\max_a(r(s,a,Z)+\beta V(s'))]`$. The second value is at least the first because the decision can respond to information.

For a fair binary shock and two actions, suppose the reward is one if the action matches the shock and zero otherwise. Choosing before observation gives expected reward one half. Choosing after observation gives reward one. This simple difference persists in more complicated models: a policy depending on tomorrow's realized shock is not admissible when tomorrow's shock is outside today's information.

A transition matrix describes conditional probabilities after the action. When the action changes transitions, the continuation value must use that action's row. Averaging actions first or reusing a transition from another policy changes the decision problem.

#### Deriving a smooth HJB equation

For a finite-horizon diffusion control problem, let the flow reward be $`f(t,x,a)`$, discount rate $`\rho`$, and dynamics $`dX=b(t,X,a)dt+\sigma(t,X,a)dW`$. Over a short interval $`h`$, the dynamic programming principle compares

```math
f(t,x,a)h+e^{-\rho h}E[V(t+h,X_{t+h})]
```

with the current value. Assume the candidate is $`C^{1,2}`$ and the short-time expansion is justified uniformly over the actions being compared. Itô's formula gives

```math
E[V(t+h,X_{t+h})]
=V+h(V_t+bV_x+\tfrac12\sigma^2V_{xx})+o(h).
```

Subtract $`V`$, divide by $`h`$, and let $`h`$ decrease to zero. The resulting equation is

```math
0=V_t+\sup_a\left(f+bV_x+\tfrac12\sigma^2V_{xx}-\rho V\right),
\qquad V(T,x)=g(x).
```

For a time-homogeneous infinite-horizon problem, a stationary value removes $`V_t`$. This is a derivation under smoothness assumptions, not a proof that every control problem has a smooth solution. A verification theorem additionally checks admissibility, integrability, and boundary or terminal conditions.

#### Example DO-E18: log saving with a random return

An agent has one unit now, saves fraction $`a\in(0,1)`$, and consumes the invested proceeds next date. The gross return is $`1/2`$ or $`3/2`$ with equal probability. Maximize $`\log(1-a)+\beta E[\log(Ra)]`$ for $`\beta>0`$.

Solution. Positivity gives $`\log(Ra)=\log a+\log R`$. The expectation is

```math
\log(1-a)+\beta\log a+\frac\beta2\log(3/4).
```

The last term is independent of $`a`$. Differentiating gives $`-1/(1-a)+\beta/a=0`$, so $`a^*=\beta/(1+\beta)`$. The second derivative is negative and utility tends to minus infinity at either endpoint, proving uniqueness.

Replacing the return by its mean one leaves the saving fraction unchanged in this specific log, two-date model, but raises value by $`-(\beta/2)\log(3/4)>0`$. The equality of saving fractions relies on the logarithmic separation; it is not a general claim that return risk has no effect on saving.

### 6.3 Worked examples

#### Example DO-E11: a stochastic continuation value

If current state one transitions to values ten and four with probabilities $`0.7,0.3`$, its expected next value is $`8.2`$. With immediate reward two and $`\beta=0.9`$, the action value is $`2+0.9(8.2)=9.38`$. Averaging next values with equal weights would give a different result and ignore the stated transition law.

#### Example DO-E12: a scalar continuous-time control

Set $`a=0`$, $`b=q=r=\rho=1`$. The Riccati equation is $`P^2+P-1=0`$, so $`P=(\sqrt5-1)/2`$. The optimal control is $`u=-Px`$, and the state decays as $`e^{-Pt}`$. The other root is negative and cannot represent a nonnegative minimum of a nonnegative cost functional.

### 6.4 Selected exercises

For Problems DO-31 through DO-35, use the source's switching problem. The inherited action is $`x\in\{a,b\}`$, the current shock is $`z\in\{A,B\}`$, shocks persist with probability $`p`$, matching today's action to today's shock earns one, and changing the inherited action costs $`\kappa\ge0`$. Assume $`0<p<1`$.

#### Problem DO-31: the switching threshold

Source: [MIT 14.451, Problem Set 4, Section 1, Question 1](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/7c7cfdadc358dee770358d67c2744525_MIT14_451F09_pset4.pdf).

Problem. In the two-action switching-cost model of this section, derive the cost threshold separating inaction and active matching and verify both candidate policies.

Solution. By symmetry let $`H=V(a,A)=V(b,B)`$ and $`L=V(a,B)=V(b,A)`$. In a matched state, keeping the action gives $`1+\beta[pH+(1-p)L]`$. In a mismatched state, keeping gives $`\beta[(1-p)H+pL]`$, while switching gives $`1-\kappa+\beta[pH+(1-p)L]`$.

If mismatch is left unchanged, subtracting the two value equations gives $`D=H-L=1+\beta(2p-1)D`$, so $`D=1/[1-\beta(2p-1)]`$. Switching would improve the mismatched state by $`D-\kappa`$. Thus inaction is optimal for

```math
\kappa\ge\kappa^*(p)=\frac1{1-\beta(2p-1)}.
```

If mismatches are corrected, the two equations have identical continuation and $`D=\kappa`$. The switching advantage is $`1-\kappa[1-\beta(2p-1)]`$, nonnegative exactly for $`\kappa\le\kappa^*(p)`$. In matched states keeping dominates switching by $`1+\kappa+\beta(2p-1)D`$, positive in either regime. These checks verify both policies, with ties on the boundary.

#### Problem DO-32: persistence and responsiveness

Source: [MIT 14.451, Problem Set 4, Section 1, Question 2](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/7c7cfdadc358dee770358d67c2744525_MIT14_451F09_pset4.pdf).

Problem. Determine how the switching threshold responds to shock persistence, holding the discount factor and switching cost fixed.

Solution. Differentiate $`\kappa^*(p)`$ to get $`2\beta/[1-\beta(2p-1)]^2>0`$. At a fixed switching cost, increasing persistence can move the parameters from inaction to active matching. A switch then creates a match that is more likely to remain useful next date. This comparison holds the discount and cost fixed.

#### Problem DO-33: the controlled transition matrices

Source: [MIT 14.451, Problem Set 4, Section 1, Question 3](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/7c7cfdadc358dee770358d67c2744525_MIT14_451F09_pset4.pdf).

Problem. Write the state transition matrices under inaction and active matching, using the state ordering $`(a,A),(a,B),(b,A),(b,B)`$.

Solution. Order states as $`(a,A),(a,B),(b,A),(b,B)`$. Under inaction,

```math
P_I=\begin{pmatrix}
p&1-p&0&0\\1-p&p&0&0\\[0pt]
0&0&p&1-p\\0&0&1-p&p
\end{pmatrix}.
```

Under active matching, the chosen next inherited action matches the current shock, while the next shock still changes randomly:

```math
P_M=\begin{pmatrix}
p&1-p&0&0\\0&0&1-p&p\\[0pt]
p&1-p&0&0\\0&0&1-p&p
\end{pmatrix}.
```

Every row sums to one. The next pair need not be matched because tomorrow's shock can differ from today's.

#### Problem DO-34: recurrent classes under the policies

Source: [MIT 14.451, Problem Set 4, Section 1, Question 4](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/7c7cfdadc358dee770358d67c2744525_MIT14_451F09_pset4.pdf).

Problem. For interior shock transition probabilities, determine the closed recurrent classes under each strict policy regime and discuss the policy-tie boundary.

Solution. In the strict inaction region, the inherited action never changes. The two closed irreducible classes are $`\{(a,A),(a,B)\}`$ and $`\{(b,A),(b,B)\}`$. In the strict active-matching region, paths can move between all four states because both shock transitions have positive probabilities. There is one closed irreducible class. Self-transitions at matched states make this chain aperiodic. On the policy-tie boundary, different optimal selections can generate different chain structures.

#### Problem DO-35: the matching regime's invariant distribution

Source: [MIT 14.451, Problem Set 4, Section 1, Question 5](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/7c7cfdadc358dee770358d67c2744525_MIT14_451F09_pset4.pdf).

Problem. Find and verify the invariant distribution in the active-matching regime.

Solution. In stationary operation the shock has probabilities one half each, and the inherited action matches the previous shock. Therefore

```math
\pi_M=\left(\frac p2,\frac{1-p}2,\frac{1-p}2,\frac p2\right).
```

For a direct check, the first coordinate of $`\pi_MP_M`$ is $`p(\pi_1+\pi_3)=p/2`$, and the second is $`(1-p)(\pi_1+\pi_3)=(1-p)/2`$. The other two are analogous. Irreducibility ensures this invariant distribution is unique.

#### Problem DO-36: distributions supported on a closed class

Source: [MIT 14.451, Problem Set 4, Section 2, Question 1](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/7c7cfdadc358dee770358d67c2744525_MIT14_451F09_pset4.pdf).

Problem. Prove that a probability distribution supported on a closed class remains supported there after one transition.

Solution. If $`E`$ is a closed class, then $`\sum_{j\in E}P_{ij}=1`$ for each $`i\in E`$. A distribution $`q`$ supported on $`E`$ satisfies

```math
\sum_{j\in E}(qP)_j
=\sum_{i\in E}q_i\sum_{j\in E}P_{ij}
=\sum_{i\in E}q_i=1.
```

Thus one step leaves its entire probability mass in the class. Closure, rather than mere communication between some states, is the needed condition.

#### Problem DO-37: uniqueness from a common positive column

Source: [MIT 14.451, Problem Set 4, Section 2, Question 2](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/7c7cfdadc358dee770358d67c2744525_MIT14_451F09_pset4.pdf).

Problem. For a finite transition matrix with a common strictly positive column, prove uniqueness of an invariant distribution and determine its support relative to a closed class.

Solution. In the finite chain, suppose $`P_{ij}>0`$ for a common $`j`$ and every $`i`$. Let $`\varepsilon=\min_iP_{ij}>0`$. If $`\varepsilon<1`$, each row decomposes as $`\varepsilon\delta_j+(1-\varepsilon)R_i`$ for a stochastic matrix $`R`$. Hence the map $`q\mapsto qP`$ contracts total-variation distance by at most $`1-\varepsilon`$. If $`\varepsilon=1`$, every row is already $`\delta_j`$ and uniqueness is immediate. Completeness of the finite simplex gives a unique invariant distribution.

A closed class $`E`$ supports an invariant distribution by applying the same finite-state existence argument to its restriction. Extending that distribution by zero outside $`E`$ yields an invariant distribution for the full chain. Uniqueness forces the full invariant distribution to give zero mass to states outside $`E`$.

#### Problem DO-38: positivity inside the recurrent class

Source: [MIT 14.451, Problem Set 4, Section 2, Question 3](https://ocw.mit.edu/courses/14-451-dynamic-optimization-methods-with-applications-fall-2009/7c7cfdadc358dee770358d67c2744525_MIT14_451F09_pset4.pdf).

Problem. Prove that an invariant distribution supported on a closed irreducible class is strictly positive at every state of that class.

Solution. Interpret an ergodic class here as a closed irreducible class. Since its invariant mass sums to one, at least one state $`i`$ has $`\pi_i>0`$. For any other state $`j`$ in the class, irreducibility gives an integer $`n`$ with $`(P^n)_{ij}>0`$. Stationarity then gives $`\pi_j=\sum_k\pi_k(P^n)_{kj}\ge\pi_i(P^n)_{ij}>0`$. Together with the zero mass outside from DO-37, the class is exactly $`\{i:\pi_i>0\}`$.

### 6.5 Additional solved problems

#### Problem DO-39: values in the switching regime

Constructed problem. Compute $`H,L`$ explicitly in both regimes of DO-31.

Solution. Under inaction, $`H+L=1+\beta(H+L)`$, so the sum is $`1/(1-\beta)`$. With $`D=1/[1-\beta(2p-1)]`$, $`H=\tfrac12[(1-\beta)^{-1}+D]`$ and $`L=\tfrac12[(1-\beta)^{-1}-D]`$.

Under matching, $`H-L=\kappa`$ and $`H=1+\beta[H-(1-p)\kappa]`$. Thus $`H=[1-\beta(1-p)\kappa]/(1-\beta)`$ and $`L=H-\kappa`$. The two formulas agree on the threshold boundary.

#### Problem DO-40: invariant distributions in the inaction regime

Constructed problem. Characterize all invariant distributions for $`P_I`$.

Solution. Each two-state shock block has invariant probabilities $`(1/2,1/2)`$. Since there is no flow between the two action classes, any mixture

```math
\pi_\lambda=(\lambda/2,\lambda/2,(1-\lambda)/2,(1-\lambda)/2),
\qquad0\le\lambda\le1,
```

is invariant. Conversely each invariant distribution has some fixed mass $`\lambda`$ in the first class and must be stationary within each occupied block, giving exactly these mixtures.

#### Problem DO-41: a continuous-time noise term

Constructed problem. Add independent state noise $`\sigma\,dW_t`$ to the linear-quadratic problem and seek $`V(x)=Px^2+C`$.

Solution. The HJB contains $`\tfrac12\sigma^2V''=\sigma^2P`$. The coefficient of $`x^2`$ gives the same Riccati equation for $`P`$, and the constant equation is $`\rho C=\sigma^2P`$, so $`C=\sigma^2P/\rho`$. The optimal feedback remains $`u=-bPx/r`$ because the additive noise does not depend on the action. Verification also requires finite discounted second moments, satisfied for the stabilizing Riccati root under this discounted cost.

#### Problem DO-42: a policy can change while value changes little

Constructed problem. In one state, two actions have the same transition and rewards zero and $`\varepsilon>0`$.

Solution. The second action is uniquely optimal for every positive $`\varepsilon`$, while at $`\varepsilon=0`$ both are optimal. The value difference from the all-zero problem is only $`\varepsilon/(1-\beta)`$ and tends to zero. Thus uniform continuity of values under reward perturbations does not imply continuity of a selected optimal action at a tie.

#### Problem DO-43: a maximum-theorem check

Constructed problem. State conditions supporting continuity of a one-step Bellman value.

Solution. If the state and action spaces are metric, the feasible-action correspondence is nonempty compact-valued and continuous, and the current-plus-continuation objective is continuous in state and action, the maximum theorem gives a continuous value and an upper-hemicontinuous nonempty compact argmax correspondence. A unique maximizing action then gives a continuous policy under these conditions. The theorem is stated here; applying it still requires checking continuity of the expected continuation term and of the feasible correspondence.

#### Problem DO-44: an error bound for a chosen policy

Constructed problem. Suppose $`\Vert Tv-v\Vert \le\delta`$ and $`\pi`$ is greedy for $`v`$. Bound its suboptimality.

Solution. The residual bound gives $`\Vert v-v^*\Vert \le\delta/(1-\beta)`$. Insert this into the greedy-policy estimate from Section 3:

```math
\Vert v^*-v_\pi\Vert \le\frac{2\beta\delta}{(1-\beta)^2}.
```

This bound can be conservative but needs only the residual and discount. A small residual has to be judged relative to $`(1-\beta)`$, especially when the decision maker is very patient.

<!-- extended-practice-start -->

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

Use discount $`\beta=\frac{1}{2}`$. Action 0 has rewards $`\begin{pmatrix}2\\0\end{pmatrix}`$ and transition $`P_0=\begin{pmatrix}\frac{3}{4}&\frac{1}{4}\\\frac{1}{2}&\frac{1}{2}\end{pmatrix}`$; action 1 has rewards $`\begin{pmatrix}0\\3\end{pmatrix}`$ and transition $`P_1=\begin{pmatrix}\frac{1}{4}&\frac{3}{4}\\\frac{1}{4}&\frac{3}{4}\end{pmatrix}`$. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta<1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`\frac{24}{7}`$ | $`\frac{8}{7}`$ |
| $`(0,1)`$ | $`\frac{13}{3}`$ | $`\frac{17}{3}`$ |
| $`(1,0)`$ | $`0`$ | $`0`$ |
| $`(1,1)`$ | $`\frac{9}{4}`$ | $`\frac{21}{4}`$ |

The policy $`\pi^*=(0,1)`$ has value $`v^*=\begin{pmatrix}\frac{13}{3}\\\frac{17}{3}\end{pmatrix}`$. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=\frac{13}{3},\quad Q(1,1)=\frac{8}{3},\qquad
Q(2,0)=\frac{5}{2},\quad Q(2,1)=\frac{17}{3}.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

#### Problem DO-51: two-state discounted control, case 2

Use discount $`\beta=\frac{4}{5}`$. Action 0 has rewards $`\begin{pmatrix}1\\2\end{pmatrix}`$ and transition $`P_0=\begin{pmatrix}\frac{4}{5}&\frac{1}{5}\\\frac{1}{5}&\frac{4}{5}\end{pmatrix}`$; action 1 has rewards $`\begin{pmatrix}3\\0\end{pmatrix}`$ and transition $`P_1=\begin{pmatrix}\frac{1}{2}&\frac{1}{2}\\\frac{3}{4}&\frac{1}{4}\end{pmatrix}`$. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta<1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`\frac{85}{13}`$ | $`\frac{110}{13}`$ |
| $`(0,1)`$ | $`\frac{25}{6}`$ | $`\frac{25}{8}`$ |
| $`(1,0)`$ | $`\frac{235}{19}`$ | $`\frac{210}{19}`$ |
| $`(1,1)`$ | $`10`$ | $`\frac{15}{2}`$ |

The policy $`\pi^*=(1,0)`$ has value $`v^*=\begin{pmatrix}\frac{235}{19}\\\frac{210}{19}\end{pmatrix}`$. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=\frac{203}{19},\quad Q(1,1)=\frac{235}{19},\qquad
Q(2,0)=\frac{210}{19},\quad Q(2,1)=\frac{183}{19}.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

#### Problem DO-52: two-state discounted control, case 3

Use discount $`\beta=\frac{2}{3}`$. Action 0 has rewards $`\begin{pmatrix}0\\4\end{pmatrix}`$ and transition $`P_0=\begin{pmatrix}\frac{1}{2}&\frac{1}{2}\\\frac{1}{4}&\frac{3}{4}\end{pmatrix}`$; action 1 has rewards $`\begin{pmatrix}2\\1\end{pmatrix}`$ and transition $`P_1=\begin{pmatrix}\frac{3}{4}&\frac{1}{4}\\\frac{1}{2}&\frac{1}{2}\end{pmatrix}`$. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta<1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`\frac{24}{5}`$ | $`\frac{48}{5}`$ |
| $`(0,1)`$ | $`1`$ | $`2`$ |
| $`(1,0)`$ | $`\frac{15}{2}`$ | $`\frac{21}{2}`$ |
| $`(1,1)`$ | $`\frac{27}{5}`$ | $`\frac{21}{5}`$ |

The policy $`\pi^*=(1,0)`$ has value $`v^*=\begin{pmatrix}\frac{15}{2}\\\frac{21}{2}\end{pmatrix}`$. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=6,\quad Q(1,1)=\frac{15}{2},\qquad
Q(2,0)=\frac{21}{2},\quad Q(2,1)=7.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

#### Problem DO-53: two-state discounted control, case 4

Use discount $`\beta=\frac{3}{4}`$. Action 0 has rewards $`\begin{pmatrix}-1\\2\end{pmatrix}`$ and transition $`P_0=\begin{pmatrix}\frac{9}{10}&\frac{1}{10}\\\frac{1}{2}&\frac{1}{2}\end{pmatrix}`$; action 1 has rewards $`\begin{pmatrix}1\\1\end{pmatrix}`$ and transition $`P_1=\begin{pmatrix}\frac{1}{5}&\frac{4}{5}\\\frac{1}{10}&\frac{9}{10}\end{pmatrix}`$. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta<1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`-\frac{19}{7}`$ | $`\frac{11}{7}`$ |
| $`(0,1)`$ | $`-\frac{5}{2}`$ | $`\frac{5}{2}`$ |
| $`(1,0)`$ | $`\frac{292}{49}`$ | $`\frac{332}{49}`$ |
| $`(1,1)`$ | $`4`$ | $`4`$ |

The policy $`\pi^*=(1,0)`$ has value $`v^*=\begin{pmatrix}\frac{292}{49}\\\frac{332}{49}\end{pmatrix}`$. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=\frac{173}{49},\quad Q(1,1)=\frac{292}{49},\qquad
Q(2,0)=\frac{332}{49},\quad Q(2,1)=\frac{295}{49}.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

#### Problem DO-54: two-state discounted control, case 5

Use discount $`\beta=\frac{9}{10}`$. Action 0 has rewards $`\begin{pmatrix}3\\1\end{pmatrix}`$ and transition $`P_0=\begin{pmatrix}\frac{2}{3}&\frac{1}{3}\\\frac{1}{3}&\frac{2}{3}\end{pmatrix}`$; action 1 has rewards $`\begin{pmatrix}1\\4\end{pmatrix}`$ and transition $`P_1=\begin{pmatrix}\frac{1}{3}&\frac{2}{3}\\\frac{2}{3}&\frac{1}{3}\end{pmatrix}`$. Rewards are received before the next state is drawn. Find an optimal stationary policy.

Solution. A policy $`\pi=(a,b)`$ selects action $`a`$ in state one and action $`b`$ in state two. Form $`r_\pi`$ and $`P_\pi`$ by selecting the corresponding reward and transition row at each state. Policy evaluation solves

```math
v_\pi=r_\pi+\beta P_\pi v_\pi,\qquad
v_\pi=(I-\beta P_\pi)^{-1}r_\pi.
```

The inverse exists because $`\beta<1`$ and a stochastic matrix has sup-norm operator bound one. Equivalently, the inverse is the convergent series $`\sum_{k\ge0}\beta^kP_\pi^k`$. Exact evaluation of all policies gives:

| Policy | Value in state one | Value in state two |
|---|---:|---:|
| $`(0,0)`$ | $`\frac{150}{7}`$ | $`\frac{130}{7}`$ |
| $`(0,1)`$ | $`33`$ | $`34`$ |
| $`(1,0)`$ | $`10`$ | $`10`$ |
| $`(1,1)`$ | $`\frac{310}{13}`$ | $`\frac{340}{13}`$ |

The policy $`\pi^*=(0,1)`$ has value $`v^*=\begin{pmatrix}33\\34\end{pmatrix}`$. To verify that this value is optimal, calculate the reward-plus-continuation values of both actions using this same vector:

```math
Q(1,0)=33,\quad Q(1,1)=\frac{313}{10},\qquad
Q(2,0)=\frac{313}{10},\quad Q(2,1)=34.
```

At each state, the selected action attains the larger entry and that entry equals the corresponding component of $`v^*`$. Thus $`v^*=Tv^*`$ for the Bellman operator $`(Tv)(s)=\max_u[r_u(s)+\beta P_u(s,\cdot)v]`$.

This operator is a contraction with modulus $`\beta`$ on the finite-dimensional sup norm. Its fixed point is unique, so the verified vector is the optimal value. A tie would allow either maximizing action at that state. Merely choosing the largest immediate reward would ignore the transition effect on future rewards; the displayed $`Q`$ values include both.

### 7.3 Stationary distributions and actual convergence

A stationary distribution remains unchanged after one transition. Convergence from a given initial state is a separate assertion. A two-state chain permits an exact formula showing why periodicity can prevent convergence even when the stationary distribution is unique.

#### Problem DO-55: two-state Markov law, case 1

Let $`P=\begin{pmatrix}\frac{4}{5}&\frac{1}{5}\\\frac{2}{5}&\frac{3}{5}\end{pmatrix}`$ on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q>0`$, the unique solution is

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

Here $`\vert 1-p-q\vert =\frac{2}{5}<1`$, so the power term vanishes and the state distribution converges to $`(\frac{2}{3},\frac{1}{3})`$. A negative second eigenvalue gives alternating deviations; a positive one gives deviations with a constant sign.

The eigenvalues of the transition matrix are one and $`1-p-q=\frac{2}{5}`$. The first eigenvalue preserves probability mass; the second controls how deviations from stationarity evolve. This calculation shows why a proof of an invariant distribution alone does not prove convergence to it. Irreducibility ensures uniqueness here, while aperiodicity rules out the alternating exceptional case.

#### Problem DO-56: two-state Markov law, case 2

Let $`P=\begin{pmatrix}\frac{1}{2}&\frac{1}{2}\\\frac{1}{2}&\frac{1}{2}\end{pmatrix}`$ on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q>0`$, the unique solution is

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

Here $`\vert 1-p-q\vert =0<1`$, so the power term vanishes and the state distribution converges to $`(\frac{1}{2},\frac{1}{2})`$. A negative second eigenvalue gives alternating deviations; a positive one gives deviations with a constant sign.

The eigenvalues of the transition matrix are one and $`1-p-q=0`$. The first eigenvalue preserves probability mass; the second controls how deviations from stationarity evolve. This calculation shows why a proof of an invariant distribution alone does not prove convergence to it. Irreducibility ensures uniqueness here, while aperiodicity rules out the alternating exceptional case.

#### Problem DO-57: two-state Markov law, case 3

Let $`P=\begin{pmatrix}\frac{9}{10}&\frac{1}{10}\\\frac{3}{10}&\frac{7}{10}\end{pmatrix}`$ on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q>0`$, the unique solution is

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

Here $`\vert 1-p-q\vert =\frac{3}{5}<1`$, so the power term vanishes and the state distribution converges to $`(\frac{3}{4},\frac{1}{4})`$. A negative second eigenvalue gives alternating deviations; a positive one gives deviations with a constant sign.

The eigenvalues of the transition matrix are one and $`1-p-q=\frac{3}{5}`$. The first eigenvalue preserves probability mass; the second controls how deviations from stationarity evolve. This calculation shows why a proof of an invariant distribution alone does not prove convergence to it. Irreducibility ensures uniqueness here, while aperiodicity rules out the alternating exceptional case.

#### Problem DO-58: two-state Markov law, case 4

Let $`P=\begin{pmatrix}0&1\\1&0\end{pmatrix}`$ on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q>0`$, the unique solution is

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

Let $`P=\begin{pmatrix}0&1\\\frac{1}{2}&\frac{1}{2}\end{pmatrix}`$ on states zero and one. Find its stationary distribution and the distribution at date $`n`$ starting from state zero.

Solution. Write $`\pi=(\pi_0,\pi_1)`$. The stationarity equation reduces to flow balance $`\pi_0p=\pi_1q`$, together with $`\pi_0+\pi_1=1`$. Since $`p,q>0`$, the unique solution is

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

Here $`\vert 1-p-q\vert =\frac{1}{2}<1`$, so the power term vanishes and the state distribution converges to $`(\frac{1}{3},\frac{2}{3})`$. A negative second eigenvalue gives alternating deviations; a positive one gives deviations with a constant sign.

The eigenvalues of the transition matrix are one and $`1-p-q=-\frac{1}{2}`$. The first eigenvalue preserves probability mass; the second controls how deviations from stationarity evolve. This calculation shows why a proof of an invariant distribution alone does not prove convergence to it. Irreducibility ensures uniqueness here, while aperiodicity rules out the alternating exceptional case.

### 7.4 Stopping with two possible offers

At each date an offer is observed. Accepting ends the problem and pays the offer immediately. Rejecting costs $`c`$ now and gives the discounted value of a fresh independent offer next date. Writing the timing explicitly determines the continuation value.

#### Problem DO-60: two-offer stopping rule, case 1

The offer is $`H=10`$ with probability $`p=\frac{1}{2}`$ and $`L=2`$ otherwise. Let $`\beta=\frac{9}{10}`$ and rejection cost $`c=1`$. Find the value before seeing the offer and an optimal acceptance rule.

Solution. Let $`v`$ be the value before the next offer is observed. After an offer $`x`$ is seen, acceptance gives $`x`$ and rejection gives the common continuation value $`C=-c+\beta v`$. Therefore

```math
v=p\max(H,C)+(1-p)\max(L,C).
```

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H<H`$.

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

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H<H`$.

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

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H<H`$.

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

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H<H`$.

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

The high offer is accepted. To see why, the strategy of rejecting forever has nonpositive discounted reward, while an immediate nonnegative offer is available; also $`v\le H`$ because no accepted payoff exceeds $`H`$ and rejection costs are nonnegative. Thus $`C\le\beta H<H`$.

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

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0<k'<Ak^\alpha`$ and tends to minus infinity at both endpoints.

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

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0<k'<Ak^\alpha`$ and tends to minus infinity at both endpoints.

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

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0<k'<Ak^\alpha`$ and tends to minus infinity at both endpoints.

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

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0<k'<Ak^\alpha`$ and tends to minus infinity at both endpoints.

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

Solution. Use the candidate $`V(k)=a+B\log k`$. The Bellman maximand for next capital $`k'`$ is $`\log(Ak^\alpha-k')+\beta a+\beta B\log k'`$. It is strictly concave on $`0<k'<Ak^\alpha`$ and tends to minus infinity at both endpoints.

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

The expression is strictly convex in $`u`$ because $`r>0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=1`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-1x_t,\qquad
\dot x_t=(-1)x_t,\qquad
x_t=e^{(-1)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=3>0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

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

The expression is strictly convex in $`u`$ because $`r>0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=2`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-1x_t,\qquad
\dot x_t=(0)x_t,\qquad
x_t=e^{(0)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=2>0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

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

The expression is strictly convex in $`u`$ because $`r>0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=\frac{1}{2}`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-1x_t,\qquad
\dot x_t=(-3)x_t,\qquad
x_t=e^{(-3)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=7>0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

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

The expression is strictly convex in $`u`$ because $`r>0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=2`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-2x_t,\qquad
\dot x_t=(-1)x_t,\qquad
x_t=e^{(-1)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=5>0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

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

The expression is strictly convex in $`u`$ because $`r>0`$. Its minimizer is $`u^*(x)=-(bP/r)x`$. Substituting gives the Riccati equation

```math
\rho P=q+2aP-\frac{b^2P^2}r.
```

For these parameters, $`P=1`$ satisfies this equation. The feedback rule and its closed-loop state are

```math
u_t=-1x_t,\qquad
\dot x_t=(-2)x_t,\qquad
x_t=e^{(-2)t}x.
```

The discounted squared state decays at rate $`\rho-2(a-b^2P/r)=6>0`$. Hence $`e^{-\rho T}V(x_T)\to0`$, including a case where a state path might grow slowly enough to remain square-integrable after discounting.

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

There is one state, actions $`a\in(0,1)`$, reward $`a`$, and a deterministic return to the same state. Let $`0<\beta<1`$. Find the value supremum and determine whether any admissible policy attains it.

Solution. Every reward is below one, so every discounted payoff is at most $`1/(1-\beta)`$. A stationary action $`a=1-\eta`$ gives payoff $`(1-\eta)/(1-\beta)`$, which approaches that bound as $`\eta\downarrow0`$. The value supremum is therefore $`V^*=1/(1-\beta)`$.

No deterministic policy attains it. For any action sequence,

```math
\frac1{1-\beta}-\sum_{t=0}^\infty\beta^ta_t
=\sum_{t=0}^\infty\beta^t(1-a_t)>0.
```

Already the first term is strictly positive. Allowing the actions to approach one later does not recover that initial shortfall. Randomization over actions strictly below one also leaves positive expected first-period shortfall.

The Bellman operator defined using a supremum is $`Tv=1+\beta v`$ and has the fixed point above. It correctly identifies the value supremum but has no maximizing action. A contraction proves a unique value in its function space; it does not replace an action-attainment assumption.

For any requested payoff gap $`\epsilon>0`$, choose $`0<\eta<\epsilon(1-\beta)`$. The corresponding constant policy is epsilon optimal. Closing the action set to $`[0,1]`$ creates an exact optimal policy, always choosing one. The missing endpoint is the entire source of nonattainment in this example.

#### Problem DO-78: multiple invariant laws and a transient state

For states one, two, and three, take $`P=\begin{pmatrix}1&0&0\\0&1&0\\p&1-p&0\end{pmatrix}`$ with $`0<p<1`$. Find every invariant distribution and the limiting law from an arbitrary initial distribution.

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

For every $`\beta<1`$ the series converge, including a periodic chain with $`\eta=-1`$. Since positive $`p,q`$ imply $`\eta<1`$, multiplying by $`1-\beta`$ and taking the limit gives $`\bar r`$ in both states.

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

Let current output be $`Az k^\alpha`$ with observed productivity $`z`$. Future productivities are iid and bounded above and away from zero. With $`0<\alpha,\beta<1`$ and log consumption utility, derive a value of the form $`a+B\log k+C\log z`$ and its saving policy.

Solution. The Bellman equation is

```math
V(k,z)=\max_{0<k'<Azk^\alpha}
\{\log(Azk^\alpha-k')+\beta E[V(k',Z')]\}.
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

The log-capital recursion is $`\ell_{t+1}=\log(\alpha\beta A)+\alpha\ell_t+\log Z_t`$. Its stationary mean is $`[\log(\alpha\beta A)+m]/(1-\alpha)`$ and variance is $`\mathop{\mathrm{Var}}\nolimits(\log Z)/(1-\alpha^2)`$. Exponentiating the mean log capital does not give mean capital in general.

#### Problem DO-82: a finite-horizon stochastic quadratic control

Minimize $`E[\int_t^T u_s^2/2\,ds+gX_T^2/2]`$ for $`dX_s=u_sds+\sigma dW_s`$, $`X_t=x`$, and $`g>0`$. Derive the value and optimal feedback, including the term caused by noise.

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
