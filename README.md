# HW-Computational-Methods-in-Finance
# HW 1
## Lec 1 Random Number Generation, Simulation of Random Variates
Generate random numbers follow Binomial distribution, Exponentially distributed, Normally (Box- Muller, Polar-Marsaglia) and calculate the probability. 
## Lec 2 Monte Carlo Simulations
* Use monte carlo to compute the expected value $A(t) = E(W_t^2+sin(W_t))$ .
* Simulate price path and compute the option price.
## Lec 3 Low Discrepancy Sequences and Simulation of Stochastic Processes
* Compute call option price (Euler’s discretization scheme, Milshtein’s discretization scheme) and greeks
* Simulate paths that stock price follows 2-factor model for stock prices with stochastic volatility, use the Full Truncation, Partial Truncation, and Reflection method
* Use 2-dimensional Halton sequence to estimate the integral.

# HW 2
## Lec 4 Binomial-tree & Trinomial-tree
Use Binomial-tree, Trinomial-tree model to price American option and calculate greeks
## Lec 5 LSMC
Use the LSMC (Least Square Monte Carlo method) to compute the American option price, using the first 𝑘 of the Laguerre Polynomials, Hermite Polynomials, Simple Monomials. 
## Lec 6 Numerical PDE method
Use Explicit Finite-Difference method, Implicit Finite-Difference method, Crank-Nicolson Finite-Difference method to price American option

# HW 3
## Lec 7 Exotic Options, Variance Swaps, Jump-Diffusions
* Assume that the value of a collateral follows a jump-diffusion process: $\frac{dV_t}{V_t^{-}} = \mu dt + \sigma dW_t + \gamma dJ_t$. The stopping time $\tau = \min \{ t\ge0, V_t\le q_t L_t \}, L_t= a - bc^{12t}$ is the first time when the relative value of the collateral (with respect to the outstanding loan balance) crosses a threshold, which will be viewed as the “optimal exercise boundary” of the option to default. We assume the embedded default option will be exercised at time $\tau$, if and only if $\tau\le T$. If the default option is exercised at time 𝜏 then the “payoff” to the borrower is: $(L_{\tau}-\epsilon V_{\tau})^+$. We need to compute the value of the default option, the default probability and the expected option exercise time.
* The stock price process follows 2-factor model. Estimate the Price of a Down-and-Out Put option with time-dependent barrier.
## Lec 8 Fixed Income Securities
* The short-term interest rate follows CIR model. Use Monte Carlo Simulation to find the price of a coupon-paying bond, the price of a European Call option a Pure Discount Bond and use Implicit Finite-Difference Method to find the price of a European Call option a Pure Discount Bond.
* The short-term interest rate follows G2++ model. Use Monte Carlo Simulation to find the price of a European Call option a Pure Discount Bond.

## Lec 9,10 MBS
Use Numerix Prepayment Model to compute the price of MBS, Option-Adjusted-Spread (OAS), IO and PO tranches. 
