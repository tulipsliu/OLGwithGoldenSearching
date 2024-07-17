# Exercise of Overlapping Generation Model with policy reform analysis

I establish a OLG model and code it in Fortran language. It is a plausible practice for my learning of economic dynamics.

It is also a gift of my birthday in 2024. Happy birthday to me. (May 31, 2024.)



# The Bellman Equation


Households' problem. The dynamic programming problem of households can be written as
$$
\begin{aligned}
V_t(z) &= \max\limits_{c,l,a^+} u(c,1-l) + \beta E \left[ V_{t+1}(z^+) \middle \vert \eta \right]	\\
\text{s.t.}\quad a^+ + p_t c &= (1+r_t^n) a + w_t^n h l + pen, \quad a^+ \geq0, l\geq0  \\
\text{and}\quad  \eta^+ &= \rho\eta + \epsilon^+ \quad \text{with}\quad
\varepsilon^+ \sim N(0,\sigma_\epsilon^2),
\end{aligned}	
$$
where $z=(j,a,\theta,\eta)$ again is the vector of individual state variables. Note that we put a time index on the value function and on prices. This will be necessary as soon as we compute transitional dynamic of the model. The terminal condition for the value function is
$$
	V_t(z) = 0 \quad\text{for}\quad z= (J+1,a,\theta,\eta),
$$
which means we assume that the household doesn't value what is happening after death.




## Enjoy


											Sincerely   yours
                                               Daniel Tulpen Liu.   