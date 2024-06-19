//
//  main.cpp
//  hw3
//
//  Created by Yujie Yi on 2024/5/19.
//

#include <iostream>
using namespace std;
#include <cmath>
#include <complex>
#include <string>
#include <tuple>
#include <utility>
#include <chrono>
#include <algorithm>
#include <vector>
#include <list>

#include <Eigen/Sparse>
#include <Eigen/Dense>
using namespace Eigen;
typedef Eigen::SparseMatrix<double> SparseMat;
typedef Eigen::VectorXd Vec;

using namespace Eigen;

#include <random>
random_device rd;
mt19937 gen(rd());
normal_distribution<> dis(0.0, 1.0);

template<typename T1, typename T2>
auto max(T1 a, T2 b) -> typename std::common_type<T1, T2>::type {
    return (a > b) ? a : b;
}
// Question 1
// Assume that the value of a collateral follows a jump-diffusion process: 𝑑𝑉𝑡/𝑉𝑡− = 𝜇𝑑𝑡 + 𝜎𝑑𝑊𝑡 + 𝛾𝑑𝐽𝑡 where 𝜇, 𝜎, 𝛾 < 0, and 𝑉0 are given, 𝐽 is a Poisson process, with intensity 𝜆1, independent of the Brownian Motion process 𝑊. 𝑉𝑡− is the value process before jump occurs at time t (if any). Consider a collateralized loan, with a contract rate per period r and maturity T on the above-collateral, and assume the outstanding balance of that loan follows this process: 𝐿𝑡 = 𝑎 − 𝑏𝑐^{12𝑡} where 𝑎 > 0, 𝑏 > 0 , 𝑐 > 1, 𝑎𝑛𝑑 𝐿0 are given. We have that 𝐿_𝑇 = 0. Define the following stopping time: 𝜏 = 𝑚𝑖𝑛{𝑡 ≥ 0: 𝑉𝑡 ≤ 𝑞_𝑡𝐿_𝑡} Here, 𝑞𝑡 is a known function of time. This stopping time is the first time when the relative value of the collateral (with respect to the outstanding loan balance) crosses a threshold, which will be viewed as the “optimal exercise boundary” of the option to default. We assume the embedded default option will be exercised at time 𝜏, if and only if 𝜏 ≤ 𝑇. If the default option is exercised at time 𝜏 then the “payoff” to the borrower is: (𝐿𝜏 − 𝜖𝑉𝜏)+.

// Write the code as a function Proj3_2func.* that takes 𝜆1 and T as parameters, setting defaults if these parameters are not supplied, and outputs the default option value, the default probability, and the expected default option exercise time. Function specification: function [D, Prob, Et] = Proj3_2func(lambda1, T)
// Define a function to model and simulate a project finance scenario using a Monte Carlo approach
VectorXd Proj3_2(float lambda1 = 0.2, float T = 5, float V0 = 20000, float L0 = 22000, float mu = -0.1, float sigma = 0.2, float gamma = -0.4, float r0 = 0.055, float delta = 0.25, float lambda2 = 0.4, float alpha = 0.7, float epsilon = 0.95, int N1 = 10000, int N2 = 10000)
{
    // Calculating terms for loan payment, risk assessment and default threshold
    float R = r0 + delta * lambda2, r = R / 12.0, n = T * 12, PMT = (L0*r)/(1-1/pow(1+r,n));
    float a = PMT/r, b = PMT/(r*pow(1+r,n)), c = 1+r;
    float beta = (epsilon - alpha)/T;
    float D = 0, Prob = 0, Et = 0;
    float Vt, qt, Lt, dt = T/N2, t;
    exponential_distribution<> exp_dis(lambda1);
    default_random_engine gen;
    normal_distribution<float> dis(0, 1);
    
    // Monte Carlo simulation for N1 trials
    for(int i = 0; i < N1; i++)
    {
        t = 0;
        Vt = V0;
        int j = 1;
        // Loop over the simulation time period
        while(t < T) {
            double jumpTime = exp_dis(gen); // Time until the next jump
            // Process the asset value between jumps
            while (t + jumpTime < T && t + jumpTime <= dt * j) {
                Vt *= exp(mu * jumpTime + sigma * sqrt(jumpTime) * dis(gen));
                Vt *= (1 + gamma);
                t += jumpTime;
                jumpTime = exp_dis(gen);
            }
            double interval = min(dt, T - t);
            Vt *= exp(mu * interval + sigma * sqrt(interval) * dis(gen));
            t += interval;
            Lt = a - b * pow(c,12*t);
            qt = alpha + beta * t;
            // Check for default condition
            if(Vt <= qt * Lt) {
                Prob += 1;
                Et += t;
                D += exp(-r * t) * max(Lt - epsilon * Vt, 0.0);
                break;
            }
            j++;
        }
    }
    // Calculate expected time to default, loss given default, and default probability
    Et /= Prob;
    D /= N1;
    Prob /= N1;
    VectorXd results(3);
    results << D, Prob, Et;
    return results;
}

// Helper function to print Eigen vectors as Python-style lists
void printAsPythonList(const VectorXd& vec) {
    cout << "[";
    for (int i = 0; i < vec.size(); ++i) {
        cout << vec[i];
        if (i < vec.size() - 1) {
            cout << ", ";
        }
    }
    cout << "]" << endl;
}

// (a) Estimate the value of the default option for the following ranges of parameters: 𝜆1 from 0.05 to 0.4 in increments of 0.05; T from 3 to 8 in increments of 1;


// (b) Estimate the default probability for the following ranges of parameters: 𝜆1 from 0.05 to 0.4 in increments of 0.05; T from 3 to 8 in increments of 1;


// (c) Find the Expected option Exercise Time of the default option, conditional on 𝜏 < 𝑇. That is, estimate 𝐸(𝜏 | 𝜏 < 𝑇) for the following ranges of parameters: 𝜆1 from 0.05 to 0.4 in increments of 0.05; T from 3 to 8 in increments of 1


// Question 2
// Consider the following 2-factor model for a stock price process, under the risk-neutral measure: {𝑑𝑆𝑡 = 𝑟𝑆𝑡𝑑𝑡 + √𝑣𝑡 𝑆𝑡 𝑑𝑊𝑡, 𝑑𝑣𝑡 = (𝛼 + 𝛽𝑣𝑡)𝑑𝑡 + 𝛾√𝑣𝑡 𝑑𝐵𝑡} where 𝑊𝑡 and 𝐵𝑡 are correlated Brownian Motion processes with 𝑑𝑊𝑡𝑑𝐵𝑡 = 𝜌𝑑𝑡. Default parameter values: 𝑣0 = 0.1, 𝛼 = 0.45, 𝛽 = −5.105, 𝛾 = 0.25, 𝑆0 = $100, 𝑟 = 0.05, 𝜌 = −0.75, 𝐾 = $100, 𝑇 = 1.

// (a) Estimate the Price (𝑃1) of a Down-and-Out Put option with the barrier at 𝑆_𝑏^1(𝑡) = 94.
// Define a constant barrier function for the Down-and-Out Put option.
float Sb1(float t, float T) {
    // This function returns a constant barrier of 94 regardless of time.
    return 94;
}

// Function to estimate the price of a Down-and-Out Put option using Monte Carlo simulation.
float DOP(float(*Sb)(float,float), float V0 = 0.1, float alpha = 0.45, float beta = -5.105, float gamma = 0.25, float S0 = 100, float r = 0.05, float rho = -0.75, float K = 100, float T = 1, int N1 = 1000, int N2 = 1000) {
    float St, Vt, dt = T / N2, dW, dB, t, flag, P_do = 0;
    default_random_engine gen;
    normal_distribution<float> dis(0, 1);

    // Run N1 simulations for the option pricing.
    for (int i = 0; i < N1; i++) {
        St = S0; // Start stock price
        Vt = V0; // Start volatility
        t = 0;
        flag = 1; // Flag to check if the option is still active (not knocked out)

        // Time stepping through the option's life
        for (int j = 0; j <= N2; j++) {
            t += dt;
            dW = sqrt(dt) * dis(gen); // Random walk component for stock
            dB = rho * dW + sqrt(1 - pow(rho, 2.0)) * sqrt(dt) * dis(gen); // Random walk for volatility

            // Update stock and volatility using stochastic differential equations
            St += r * St * dt + sqrt(max(Vt, 0)) * St * dW;
            Vt += (alpha + beta * max(Vt, 0)) * dt + gamma * sqrt(max(Vt, 0)) * dB;

            // Check if the barrier is breached
            if (St <= Sb(t, T)) {
                flag = 0;
                break;
            }
        }

        // If not knocked out, calculate payoff
        if (flag == 1) {
            P_do += max(K - St, 0.0f);
        }
    }

    // Average the payoff and discount back to present value
    P_do *= exp(-r * T) / N1;
    return P_do;
}

// (b) Estimate the Price (𝑃2) of a Down-and-Out Put option with time-dependent barrier 𝑆_𝑏^2(𝑡) = 6t/𝑇 + 91.
float Sb2(float t, float T)
{
    return 6*t/T + 91;
}

// (c) Estimate the Price (𝑃3) of a Down-and-Out Put option with time-dependent barrier 𝑆_𝑏^3(𝑡) = − 6t/𝑇 + 97.
float Sb3(float t, float T)
{
    return -6*t/T + 97;
}

// Notes: • All options in parts a, b, and c have payoffs similar to the European Put option; however, these options become void (the contract is cancelled), if the underlying security’s price crosses it and goes below the barrier 𝑆_𝑏^𝑖(𝑡) at any time during the life of the option.
// • You may use any method to price the securities – Monte Carlo Simulations, Binomial or Trinomial Tree Methods, the PDE approach, etc.
// • If you use Monte Carlo simulations in this problem, use the “full truncation” method to simulate the volatility-process, 𝑣𝑘+1. The Euler discretization scheme, in that case, will be as follows: 𝑆𝑘+1 = 𝑆𝑘 + 𝑟𝑆𝑘∆𝑡 + 𝑆𝑘√𝑣_𝑘^+ (∆𝑊_{𝑘+1}), 𝑣𝑘+1 = 𝑣𝑘 + (𝛼 + 𝛽𝑣_𝑘^+)∆𝑡 + 𝛾√𝑣_𝑘^+ (∆𝐵_{𝑘+1}) where 𝑣_𝑘^+ = max(𝑣𝑘 ,0).


// Question 3
// Assume the dynamics of the short-term interest rate, under the risk-neutral measure, are given by the following SDE (CIR model): 𝑑𝑟𝑡 = 𝜅(𝑟̅ − 𝑟𝑡)𝑑𝑡 + 𝜎√𝑟𝑡𝑑𝑊𝑡 with 𝑟0 = 5%, 𝜎 = 12%, 𝜅 = 0.92, 𝑟̅ = 5.5%.
// (a) Use Monte Carlo Simulation to find the price of a coupon-paying bond, with Face Value of $1,000, paying semiannual coupons of $30, maturing in 𝑇 = 4 years: 𝑃(0, 𝐶, 𝑇) = 𝔼_0^∗[∑_𝑖=1^8 𝐶𝑖 ∗ 𝑒𝑥𝑝 (− ∫_0^Ti 𝑟(𝑠)𝑑𝑠 )] where 𝐶 = {𝐶𝑖 = $30 for 𝑖 = 1,2, ... ,7; and 𝐶8 = $1,030}, and 𝑇 = {𝑇1, 𝑇2, 𝑇3, 𝑇4, 𝑇5, 𝑇6, 𝑇7, 𝑇8} = {0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4}.

// Simulate interest rates using the Cox-Ingersoll-Ross (CIR) model
MatrixXd CIR_simulate(float r0, float sigma, float k, float r_bar, float T, float dt, int N) {
    int n = int(T / dt);
    MatrixXd r(N, n);
    default_random_engine gen;
    normal_distribution<float> dis(0, 1);
    float dW;
    r.col(0) = VectorXd::Constant(N, r0);  // Initialize the first column with the initial rate r0
    
    // Simulate the rate paths
    for (int i = 0; i < N; i++) {
        for (int j = 1; j < n; j++) {
            dW = sqrt(dt) * dis(gen);  // Generate a random shock
            // Update the rate using the CIR model differential equation
            r(i, j) = r(i, j-1) + k * (r_bar - r(i, j-1)) * dt + sigma * sqrt(max(r(i, j-1), 0.0f)) * dW;
        }
    }
    return r;
}

// Calculate the price of a coupon bond under the CIR model
float CIR_Coupon_Bond(float face_value, float coupon, float r0, float sigma, float k, float r_bar, float T, float dt, list<float> T_list, int N = 1000) {
    MatrixXd r = CIR_simulate(r0, sigma, k, r_bar, T, dt, N);
    float price = 0;

    // Compute the price by discounting the coupon payments and the face value
    for (float t : T_list) {
        MatrixXd rt = r.block(0, 0, N, int(t / dt));
        price += coupon * (-rt.rowwise().sum() * dt).array().exp().mean();  // Discounting coupon payments
        if (t == T_list.back()) {
            price += face_value * (-rt.rowwise().sum() * dt).array().exp().mean();  // Discounting the face value
        }
    }
    return price;
}




// (b) Use Monte Carlo Simulation to find at time 𝑡 = 0 the price 𝑐𝑀𝐶(𝑡, 𝑇, 𝑆) of a European Call option, with strike price of 𝐾 = $980 and expiration in 𝑇 = 0.5 years on a Pure Discount Bond that has Face Value of $1,000 and matures in 𝑆 = 1 year: 𝑐_{𝑀𝐶}(𝑡, 𝑇, 𝑆) = 𝔼_𝑡^∗[𝑒𝑥𝑝 (− ∫_t^T 𝑟(𝑢)𝑑𝑢) ∗ max(𝑃(𝑇, 𝑆) − 𝐾, 0)]

// Monte Carlo Simulation to find the price of a European Call option on a Pure Discount Bond
float CMC_Bond(float face_value, float K, float T, float S, float r0, float sigma, float k, float r_bar, float dt, int N = 1000, int M = 1000) {
    MatrixXd r1 = CIR_simulate(r0, sigma, k, r_bar, T, dt, N);
    float price = 0;

    // Simulate the option price
    for (int i = 0; i < N; i++) {
        MatrixXd r2 = CIR_simulate(r1(i, int(T / dt) - 1), sigma, k, r_bar, S - T, dt, M);
        price += exp(-r1.row(i).sum() * dt) * max(face_value * (-r2.rowwise().sum() * dt).array().exp().mean() - K, 0.0);
    }
    return price / N;
}


// (c) Use the Implicit Finite-Difference Method to find at time 𝑡 = 0 the price 𝑐𝑃𝐷𝐸(𝑡, 𝑇, 𝑆) of a European Call option, with strike price of 𝐾 = $980 and expiration in 𝑇 = 0.5 years on a Pure Discount Bond that has Face Value of $1,000 and matures in 𝑆 = 1 year. The PDE is given as follows 𝜕𝑐/𝜕𝑡 + 1/2 𝜎2𝑟 𝜕2𝑐/𝜕𝑟2 + 𝜅(𝑟̅ − 𝑟) 𝜕𝑐/𝜕𝑟 − 𝑟𝑐 = 0 with 𝑐(𝑇, 𝑇, 𝑆) = max(𝑃(𝑇, 𝑆) − 𝐾, 0), and 𝑃(𝑇, 𝑆) is computed explicitly.
// Function to calculate the price of a pure discount bond using the CIR model
float pure_bond(float k, float sigma, float r_bar, float t, float S, float rt) {
    float h1, h2, h3, A, B;
    // Calculate parameters for the analytic solution
    h1 = sqrt(pow(k, 2) + 2 * pow(sigma, 2));  // h1 is a helper coefficient in the formula, dependent on the mean reversion speed k and volatility sigma
    h2 = (k + h1) / 2;  // h2 modifies the mean reversion speed by incorporating the effect of h1
    h3 = 2 * k * r_bar / pow(sigma, 2);  // h3 represents the long-term mean level adjusted by the volatility

    // A and B are components of the formula to compute the bond price analytically
    A = pow((h1 * exp(h2 * (S - t))) / (h2 * (exp(h1 * (S - t)) - 1) + h1), h3);  // A calculates the discounting factor
    B = (exp(h1 * (S - t)) - 1) / (h2 * (exp(h1 * (S - t)) - 1) + h1);  // B is a discounting factor multiplier that adjusts for the change in rates over time

    // The price of the bond is then calculated using the exponential of the negative product of B and the current interest rate (rt)
    return A * exp(-B * rt);
}

// Implicit Finite-Difference Method to price a European Call option
VectorXd IFD_Bond(float face_value, float K, float T, float S, float r0, float sigma, float k, float r_bar, int N = 200, int M = 1000) {
    float dr = 3 * r0 / M, dt = T / N;
    MatrixXd Pumd(N+1, 3), A(N-1, N-1);
    VectorXd C(N-1), B(N-1);

    // Setup the coefficients for the finite difference matrix
    for (int i = 1; i <= N; i++) {
        Pumd.row(i) << dt * (-pow(sigma/dr, 2) * i / 2 + k * (i * dr - r_bar) / (2 * dr)),
                       1 + pow(sigma, 2) * i * dt / dr + i * dr * dt,
                       dt * (-pow(sigma/dr, 2) * i / 2 - k * (i * dr - r_bar) / (2 * dr));

        if (i == 1) {
            A.row(N-2).segment(N-3, 2) = Pumd.row(1).segment(0, 2);
        } else if (i == N-1) {
            A.row(0).segment(0, 2) = Pumd.row(N-1).segment(1, 2);
        } else if (i == N) {
            break;
        } else {
            A.row(N-i-1).segment(N-i-2, 3) = Pumd.row(i);
        }
    }

    // Set boundary conditions for the option at maturity
    for (int i = 1; i <= N-1; i++) {
        C(i-1) = max(face_value * pure_bond(k, sigma, r_bar, T, S, (N-i) * dr) - K, 0.0);
    }

    // Solve the system using LU decomposition to find the option price at t=0
    for (int i = M-1; i >= 0; i--) {
        B(N-2) = Pumd(N, 2) * (face_value * pure_bond(k, sigma, r_bar, i * dt, S, (N-1) * dr) - K);
        C = A.lu().solve(C - B);
    }
    return C;
}

/*
VectorXd IFD_Bond(float face_value, float K, float T, float S, float r0, float sigma, float k, float r_bar, int N = 200, int M = 1000)
{
    float dr = 2 * r0 / M, dt = T / N;
    MatrixXd Pumd(N+1,3);
    MatrixXd A(N+1,N+1);
    VectorXd C(N+1), B(N+1);
    A(0,0) = 1;
    A(0,1) = -1;
    A(N,N-1) = 1;
    A(N,N) = -1;
    for(int i = 1; i < N; i++)
    {
        Pumd.row(i) << dt*(-pow(sigma,2)*i/(2*dr)-k*r_bar/(2*dr)+k*i/2), 1+pow(sigma,2)*i*dt/dr+i*dr*dt, dt*(-pow(sigma,2)*i/(2*dr)+k*r_bar/(2*dr)-k*i/2);
        A.row(N-i).segment(N-i-1,3) = Pumd.row(i);
    }
    //initiate the boudary value
    for(int i = 0; i <= N; i++){C(i) = max(face_value * cal_P(k,sigma,r_bar,T,S,i*dr) - K,0);}
    
    // iterate to get Put price(t=0)
    B(N) = face_value * (cal_P(k,sigma,r_bar,T,S,0) - cal_P(k,sigma,r_bar,T,S,dr));
    cout << B(N) << endl;
    for(int i = M-1; i >= 0; i--)
    {
        B.segment(1,N-1) = C.segment(1,N-1);
        C = A.lu().solve(B);
    }
    return C;
}
*/




// Question 4
// Assume the dynamics of the short-term interest rate, under the risk-neutral measure, are given by the following system of SDEs (G2++ model): 𝑑𝑥𝑡 = −𝑎𝑥𝑡𝑑𝑡 + 𝜎𝑑𝑊𝑡, 𝑑𝑦𝑡 = −𝑏𝑦𝑡𝑑𝑡 + 𝜂𝑑𝑊𝑡2, 𝑟𝑡 = 𝑥𝑡 + 𝑦𝑡 + 𝜙𝑡, Default parameter values: 𝑥0 = 𝑦0 = 0, 𝜙0 = 𝑟0 = 5.5%, 𝑑𝑊𝑡1𝑑𝑊𝑡2 = 𝜌𝑑𝑡, 𝜌 = 0.7, 𝑎 = 0.1, 𝑏 = 0.3, 𝜎 =5%, 𝜂 = 9%. Assume 𝜙𝑡 = 𝑐𝑜𝑛𝑠𝑡 = 5.5% for any 𝑡 ≥ 0. Use Monte Carlo Simulation to find at time 𝑡 = 0 the price 𝑝(𝑡, 𝑇, 𝑆, 𝐾, 𝜌) of a European Put option, with strike price of 𝐾 = $950, expiration in 𝑇 = 0.5 years on a Pure Discount Bond with Face value of $1,000 that matures in 𝑆 = 1 year. Compare it with the price found by the explicit formula and comment on it.

// Function to simulate the G2++ model for interest rates over a time period T
MatrixXd simulate_G2(float a, float b, float r0, float rho, float sigma, float eta, float phi, float T, float dt, int N, VectorXd& x0, VectorXd& y0) {
    float x, y, w1, w2;
    int M = int(T / dt);  // Number of time steps
    MatrixXd r(N, M + 1);  // Matrix to store simulated interest rates
    default_random_engine gen;
    normal_distribution<float> dis(0, 1);

    // Simulate each path
    for (int i = 0; i < N; i++) {
        x = x0(i);
        y = y0(i);
        r(i, 0) = r0;  // Set initial rate
        
        // Time evolution of the rate
        for (int j = 1; j <= M; j++) {
            w1 = sqrt(dt) * dis(gen);
            w2 = rho * w1 + sqrt(1 - pow(rho, 2)) * sqrt(dt) * dis(gen);
            x += -a * x * dt + sigma * w1;
            y += -b * y * dt + eta * w2;
            r(i, j) = x + y + phi;  // Calculate the rate at each step
        }
        x0(i) = x;  // Update last values for continuation
        y0(i) = y;
    }
    return r;
}

// Monte Carlo simulation to price a European call option on a bond
float Monte_put_G2(float a, float b, float r0, float rho, float sigma, float eta, float phi, float T, float S, float dt, float face_value, float K, int N = 1000, int M = 1000) {
    VectorXd x0 = VectorXd::Zero(N), y0 = VectorXd::Zero(N);
    MatrixXd r = simulate_G2(a, b, r0, rho, sigma, eta, phi, T, dt, N, x0, y0);
    float price = 0;

    // Price the option by averaging the discounted payoff from all simulated paths
    for (int i = 0; i < N; i++) {
        VectorXd x2 = VectorXd::Constant(M, x0(i)), y2 = VectorXd::Constant(M, y0(i));
        MatrixXd r2 = simulate_G2(a, b, r(i, int(T/dt)), rho, sigma, eta, phi, S - T, dt, M, x2, y2);
        price += exp(-r.row(i).sum() * dt) * max(face_value * (-r2.rowwise().sum() * dt).array().exp().mean() - K, 0.0);
    }
    return price / N;
}

// Analytical pricing of a bond under the two-factor Gaussian model
float two_factor_price(float t, float T, float phi, float a, float b, float sigma, float eta, float rho, float xt, float yt) {
    float A, B, C, D, E, F;
    A = -phi * (T - t);  // Discount factor from the model's long-term mean adjustment
    B = -(1 - exp(-a * (T - t))) / a * xt;  // Contribution from factor x
    C = -(1 - exp(-b * (T - t))) / b * yt;  // Contribution from factor y
    D = pow(sigma, 2) / (2 * pow(a, 2)) * (T - t + 2/a * exp(-a * (T-t)) - 1/(2*a) * exp(-2 * a * (T-t)) - 3/(2*a));
    E = pow(eta, 2) / (2 * pow(b, 2)) * (T - t + 2/b * exp(-b * (T-t)) - 1/(2*b) * exp(-2 * b * (T-t)) - 3/(2*b));
    F = sigma * eta * rho / (a * b) * (T - t + (exp(-a * (T-t)) - 1) / a + (exp(-b * (T-t)) - 1) / b - (exp(-(a+b) * (T-t)) - 1) / (a+b));
    return exp(A + B + C + D + E + F);  // Total bond price
}

// Calculate the cumulative distribution function for a normal distribution
double norm_cdf(double value) {
    return 0.5 * erfc(-value * M_SQRT1_2);
}

// Explicit formula for pricing a European call option on a pure discount bond
float explicit_put_G2(float T, float S, float face_value, float K, float r0, float sigma, float eta, float a, float b, float rho, float phi) {
    float P_S = two_factor_price(0, S, phi, a, b, sigma, eta, rho, 0, 0), P_T = two_factor_price(0, T, phi, a, b, sigma, eta, rho, 0, 0);
    float Sigma = pow(sigma, 2) / (2 * pow(a, 3)) * pow((1 - exp(-a * (S-T))), 2) * (1 - exp(-2 * a * T)) + pow(eta, 2) / (2 * pow(b, 3)) * pow((1 - exp(-b * (S-T))), 2) * (1 - exp(-2 * b * T)) + 2 * sigma * eta * rho / (a * b * (a + b)) * (1 - exp(-a * (S-T))) * (1 - exp(-b * (S-T))) * (1 - exp(-(a+b) * T));
    float d1 = log(P_T * K / P_S / face_value) / sqrt(Sigma) - 0.5 * sqrt(Sigma), d2 = d1 + sqrt(Sigma);
    return -P_S * face_value * norm_cdf(d1) + P_T * K * norm_cdf(d2);  // Black-Scholes type formula
}




// Question 5
// Consider a 30-year MBS with a fixed weighted-average-coupon, 𝑊𝐴𝐶 = 8%. Monthly cash flows are starting in January of this year. The Notional Amount of the Pool is $100,000. Use the CIR model of interest rates, 𝑑𝑟𝑡 = 𝜅(𝑟̅ − 𝑟𝑡)𝑑𝑡 + 𝜎√𝑟𝑡𝑑𝑊𝑡, with the following default parameters: 𝑟0 = 0.078, 𝑘 = 0.6, 𝑟̅ = 0.08, 𝜎 = 0.12. Consider the Numerix Prepayment Model in all problems below.

// (a) Compute the price of the MBS. The code should be generic: the user is prompted for inputs and the program runs and gives the output.

// Function to simulate and price Mortgage Backed Securities (MBS)
tuple<float, float> MBS(float WAC, float r0, float k, float r_bar, float sigma, float notion, int T, bool tranches = false, int N = 1000) {
    float r = WAC / 12;  // Monthly weighted average coupon
    MatrixXd rt(N, (T+10)*12+1), PV(N, T*12+1);  // rt stores rate paths, PV stores present values
    float price = 0, interest = 0, principal = 0, discount;
    VectorXd SY(12);  // Seasonality factors
    SY << 0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22, 1.23, 0.98;
    
    rt.col(0).setConstant(r0);  // Initialize rates
    default_random_engine gen;
    normal_distribution<float> dis(0, 1);

    // Simulate interest rate paths using a modified CIR model
    for (int t = 1; t <= (T+10)*12; t++) {
        for (int i = 0; i < N; i++) {
            rt(i, t) = rt(i, t-1) + k * (r_bar - rt(i, t-1)) / 12 + sigma * sqrt(rt(i, t-1)/12) * dis(gen);
        }
    }

    PV.col(0).setConstant(notion);  // Initial principal
    // Calculate cash flows for MBS
    for (int t = 1; t <= T*12; t++) {
        for (int i = 0; i < N; i++) {
            float bond = pure_bond(k, sigma, r_bar, 0, 10, rt(i, t-1));
            float rT = -log(bond) / 10;  // Recalculate interest rate based on bond price
            float CPR = (0.28 + 0.14 * atan(-8.57 + 430 * (WAC - rT))) * (0.3 + 0.7 * PV(i, t-1) / notion) * min(1, t/30) * SY((t-1)%12);  // Calculate prepayment rate
            float ct = (PV(i, t-1) * r) / (1 - pow(1+r, -T*12+t-1)) + PV(i, t-1) * (1-r*(1/(1-pow(1+r, -T*12+t-1))-1)) * (1-pow(1-CPR, 1/12.0));  // Total cash flow
            discount = exp(-rt.row(i).segment(0, t).sum() / 12);  // Discount factor
            interest += discount * r * PV(i, t-1);  // Accumulate interest
            principal += discount * (ct - r * PV(i, t-1));  // Accumulate principal repayment
            price += discount * ct;  // Accumulate total price
            PV(i, t) = PV(i, t-1) - (ct - r * PV(i, t-1));  // Update remaining principal
        }
    }

    if (!tranches)
        return make_tuple(price / N, interest / N);  // Return average price and interest for standard MBS
    else
        return make_tuple(interest / N, principal / N);  // Return interest and principal for tranches
}



// (b) Compute the Option-Adjusted-Spread (OAS) if the Market Price of MBS is 𝑃̂ =$98,000.
// Function to compute Option-Adjusted-Spread (OAS) given the Market Price of the MBS
float cal_OAS(float marketPrice, float lower, float upper, float WAC, float r0, float k, float r_bar, float sigma, float notion, int T) {
    auto priceFunction = [&](float oas) {
        auto [price, _] = MBS(WAC, r0 + oas, k, r_bar, sigma, notion, T);
        return marketPrice - price;
    };

    float low = lower, high = upper, oas = 0;
    // Bisection method to find the OAS
    while (high - low > 0.0001) {
        oas = (low + high) / 2;
        float priceDiff = priceFunction(oas);
        if (fabs(priceDiff) < 0.01)  // Convergence criterion
            return oas;
        else if (priceDiff > 0)
            low = oas;
        else
            high = oas;
    }
    return oas;
}

// (c) Consider the MBS described above and the IO and PO tranches. Price the IO and PO tranches
// Function to price the Interest Only (IO) and Principal Only (PO) tranches of an MBS
pair<float, float> cal_IO_PO(float WAC, float r0, float k, float r_bar, float sigma, float notion, int T) {
    auto [interest, principal] = MBS(WAC, r0, k, r_bar, sigma, notion, T, true);
    return {interest, principal};  // Return the calculated interest and principal specifically for tranches
}


int main()
{
    // Q1, Q3(c)
    // Q1
    // printAsPythonList(Proj3_2());
    /* cout << "[";
    for(float lambda = 0.05; lambda <= 0.45; lambda = lambda + 0.05)
    {
        for(float T = 3; T <= 8; T = T + 1)
        {
            // cout << lambda << ", " << T << endl;
            printAsPythonList(Proj3_2(lambda, T));
            cout << ", ";
        }
    }
    cout << "]";*/
    
    // Q2
    // cout << DOP(Sb1) << endl;
    // cout << DOP(Sb2) << endl;
    // cout << DOP(Sb3) << endl;
    
    // Q3
    // list<float> T_list = {0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4};
    // cout << "Q3(a):" << CIR_Coupon_Bond(1000, 30, 0.05, 0.12, 0.92, 0.055, 4, 0.01, T_list) << endl;
    // cout << "Q3(b):" << CMC_Bond(1000, 980, 0.5, 1, 0.05, 0.12, 0.92, 0.055,0.01) << endl;
    // cout << "Q3(c):" << endl;
    // VectorXd C = IFD_Bond(1000, 980, 0.5, 1, 0.05, 0.12, 0.92, 0.055);
    // printAsPythonList(C);
    // cout << C[50] << endl;
    
    // Q4
    /*
    cout << "[";
    for(float rho = -0.7; rho <= 0.8; rho = rho + 0.1)
    {
        cout << Monte_put_G2(0.1, 0.3, 0.055, rho, 0.05, 0.09, 0.055, 0.5, 1, 0.01, 1000, 950) << ", ";
    }
    cout << "]" << endl;
    
    cout << "[";
    for(float rho = -0.7; rho <= 0.8; rho = rho + 0.1)
    {
        cout << explicit_put_G2(0.5, 1, 1000, 950, 0.055, 0.05, 0.09, 0.1, 0.3, rho, 0.055) << ", ";
    }
    cout << "]" << endl;*/
    
    // cout << "Q4 Monte Carlo:" << Monte_put_G2(0.1, 0.3, 0.055, 0.7, 0.05, 0.09, 0.055, 0.5, 1, 0.5/100, 1000, 950) << endl;
    // cout << "Q4 Explicit:" << explicit_put_G2(0.5, 1, 1000, 950, 0.055, 0.05, 0.09, 0.1, 0.3, 0.7, 0.055) << endl;
    
    // Q5
    float WAC, r0, k, r_bar, sigma, notion, T;
    // cout << "Input WAC, r0, k, r_bar, sigma, notion, T: ";  // 0.08 0.078 0.6 0.08 0.12 100000 30
    // cin >> WAC >> r0 >> k >> r_bar >> sigma >> notion >> T;
    WAC = 0.08;
    r0 = 0.078;
    k = 0.6;
    r_bar = 0.08;
    sigma = 0.12;
    notion = 100000;
    T = 30;
    /*
    auto [price,interest] = MBS(WAC, r0, k, r_bar, sigma, notion, T);
    cout << "Q5(a)" << price << endl;
    cout << "Q5(b)" << cal_OAS(98000, 0.003, 0.005, 0.08, 0.078, 0.6, 0.08, 0.12, 100000, 30) << endl;
    auto [IO, PO] = cal_IO_PO(0.08, 0.078, 0.6, 0.08, 0.12, 100000, 30);
    cout << "Q5(a): IO is " << IO << ", PO is " << PO << endl;*/
    cout << "[";
    for(r_bar = 0.04; r_bar <= 0.11; r_bar = r_bar + 0.01)
    {
        auto [price,interest] = MBS(WAC, r0, k, r_bar, sigma, notion, T);
        cout << price << ", " << endl;
    }
    cout << "]" << endl;
    return 0;
}
