//
//  main.cpp
//  computational hw
//
//  Created by Yujie Yi on 2024/4/15.
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

#include <random>
random_device rd;
mt19937 gen(rd());
normal_distribution<> dis(0.0, 1.0);
// Question 1
// Using the LGM method, generate Uniformly distributed random numbers on [0,1] to do the following:
// (a) Generate 1,000 random numbers with Binomial distribution with 𝑛 = 44 and 𝑝 = 0.64. Compute the probability that the random variable X, that has Binomial (44, 0.64) distribution, is at least 40: 𝑃(𝑋 ≥ 40).

unsigned long long current_random = 566;
unsigned long long a = pow(2,19)-1, b = 0, m = pow(2,61)-1;
int N_simulate = 10000;
double uniform()
{
    double ans = (a * current_random + b) % m;
    current_random = ans;
    return ans/m;
}
float prob(float *arr, float num, int size)
{
    // calculate the P(arr>=num)
    float ans=0;
    for (int i = 0; i < size; ++i)
    {
        if(arr[i] >= num)
            ans += 1;
    }
    return ans/size;
}
void binomial(int N, int n, float p)
{
    // simulate binomial random numbers
    float arr[N];
    for (int i = 0; i < N; ++i)
    {
        arr[i] = 0;
        for (int j = 0; j < n; ++j)
        {
            if(uniform() < p)
                arr[i] += 1;
        }
    }
    // calculate 𝑃(𝑋 ≥ 40)
    cout << "P(X >= 40) = " << prob(arr, 35, N) << endl;
    // print 1000 binomial random numbers
    cout << "The 1000 binomial random numbers are ";
    for (int i = 0; i < N; ++i)
    {
        cout << arr[i] << ", ";
    }
    cout << endl;
}

// (b) Generate 10,000 Exponentially distributed random numbers with parameter 𝜆 = 1.5. Estimate 𝑃(𝑋 ≥ 1); 𝑃(𝑋 ≥ 4); and compute the empirical mean and the standard deviation of the sequence of 10,000 numbers. Draw the histogram by using the 10,000 numbers you have generated.
template<typename T>
T arr_mean(T array[], int N)
{
    // int N = sizeof(array) / sizeof(array[0]);
    T sum = 0;
    for(int i = 0; i < N; i++)
    {
        sum += array[i];
    }
    return sum/float(N);
}

void exponential(float lambda, float X1, float X2, int N = 10000)
{
    // simulate exponentially distributed random numbers
    float arr[N], std=0, mean;
    for(int i = 0; i < N; i++)
    {
        arr[i] = - lambda * log(1-uniform());
    }
    mean = arr_mean(arr,N);
    for(int i = 0; i < N; i++)
    {
        std += pow(arr[i]-mean,2);
    }
    cout << "The empirical mean is " << mean << " and the standard deviation is " << sqrt(std/N) << endl;
    cout << "P(X >= 1) = " << prob(arr, X1, N) << ", P(X >= 4) = " << prob(arr, X2, N) << endl;
    // print 10,000 Exponentially random numbers
    cout << "The 10,000 Exponentially random numbers: " ;
    for (int i = 0; i < N; ++i)
    {
        cout << arr[i] << ", ";
    }
    cout << endl;
}


// (c) Generate 5,000 Normally distributed random numbers with mean 0 and variance 1, by using the Box- Muller Method.
float box_muller(int N)
{
    float arr[2*N], u1, u2;
    for(int i = 0; i < N; i++)
    {
        u1 = uniform();
        u2 = uniform();
        arr[2*i] = sqrt(-2*log(u1))*cos(2*M_PI*u2);
        arr[2*i+1] = sqrt(-2*log(u1))*sin(2*M_PI*u2);
    }
    if(N>1)
    {
        cout << "The 5,000 box_muller random numbers: " ;
        for (int i = 0; i < 2*N; ++i)
        {
            cout << arr[i] << ", ";
        }
        cout << endl;
    }
    return arr[0];
}

// (d) Generate 5,000 Normally distributed random numbers with mean 0 and variance 1, by using the Polar-Marsaglia method.
float polar_marsaglia(int N)
{
    float arr[2*N];
    float u1, u2, v1, v2, w;
    int i = 0;
    while(i < 2*N)
    {
        u1 = uniform();
        u2 = uniform();
        v1 = 2*u1-1;
        v2 = 2*u2-1;
        w = pow(v1,2)+pow(v2,2);
        if (w <= 1) {
            float multiplier = sqrt((-2 * log(w)) / w);
            arr[i] = v1 * multiplier; // First normal variable
            arr[i+1] = v2 * multiplier; // Second normal variable
            i += 2;
        }
    }
    cout << "The 5,000 polar_marsaglia random numbers: " ;
    for (int i = 0; i < 2*N; ++i)
    {
        cout << arr[i] << ", ";
    }
    cout << endl;
    return arr[0];
}


/* (e) Now compare the efficiencies of the two above-algorithms, by comparing the execution times to generate 5,000 normally distributed random numbers by the two methods. Which one is more efficient? If you do not see a clear difference, you need to increase the number of generated realizations of random variables to 10,000, 20,000, etc */





/* --Question 2-- */
/* (a) Estimate the following expected values by simulation: A(t)=E(W_t^2+sin(W_t)), B(t)=E(e^(t/2)cos(W_t)), for t = 1,3,5. Here, W_t is a Standard Wiener Process. */
void calculate_AB(int t)
{
    double w, A[N_simulate], B[N_simulate];
    for(int i = 0; i < N_simulate; i++)
    {
        w = sqrt(t) * dis(gen);
        //cout << w << endl;
        A[i] = pow(w,2.0) + sin(w);
        B[i] = exp(t/2.0) * cos(w);
    }
    cout << "At time " << t << " - A(t):" << arr_mean(A,N_simulate) << ", B(t):" << arr_mean(B,N_simulate) << endl;
}

/* (b) How are the values of 𝐵(𝑡) (for the cases t = 1,3,5) related? */


/* (c) Now use a variance reduction technique (whichever you want) to compute the expected value 𝐵(5). Do you see any improvements? Comment on your findings. */
float B_anti_vars(int t)
{
    float w, u, B[N_simulate];
    for(int i = 0; i < N_simulate; i++)
    {
        w = sqrt(t) * dis(gen);
        u = - 0.3 * w + sqrt(1-pow(0.3,2.0)) * t * dis(gen);
        B[i] = (exp(t/2.0)*cos(w) + exp(t/2.0)*cos(w))/2.0;
    }
    return arr_mean(B, N_simulate);
}


// Antithetic Variates
// A(t) = sum((W_t^2+sin(W_t) + (-W_t)^2+sin(-W_t))/2)/N = sum(W_t^2)/N
// B(t) = sum((exp(t/2)*cos(W_t) + exp(t/2)*cos(-W_t))/2)/N = sum(exp(t/2)*cos(W_t))/N




/* --Question 3-- */
/* Let S_t be a Geometric Brownian Motion process: S_t=S_0e^(sigma W_t+(r-sigma^2/2)t) , where r=0.055,sigma=0.2,𝑆_0=$100; W_t is a Standard Brownian Motion process (Standard Wiener process).

 (a) Estimate the price c of a European Call option on the stock with 𝑇=5,𝑋=$100 by using Monte Carlo simulation. */
template<typename T1, typename T2>
auto max(T1 a, T2 b) -> typename std::common_type<T1, T2>::type {
    return (a > b) ? a : b;
}
float call(float S0, float r, float sigma, float T, float X)
{
    float arr[N_simulate];
    for(int i = 0; i < N_simulate; i++)
    {
        arr[i] = max(exp(-r*T)*(S0*exp(sigma*sqrt(T)*dis(gen)+(r-pow(sigma,2.0)/2.0)*T) - X), 0);
    }
    return arr_mean(arr, N_simulate);
}

/* (b) Compute the exact value of the option c using the Black-Scholes formula. */
double norm_cdf(double value)
{
    return 0.5 * erfc(-value * M_SQRT1_2);
}

float bs_call(float S0, float r, float sigma, float t, float X)
{
    double d1 = (log(S0 / X) + (r + 0.5 * sigma * sigma) * t) / (sigma * sqrt(t));
    double d2 = d1 - sigma * sqrt(t);
    return S0 * norm_cdf(d1) - X * exp(-r * t) * norm_cdf(d2);
}

/* (c) Use variance reduction techniques (whichever one(s) you want) to estimate the price in part (a) again using the same number of simulations. Did the accuracy improve? Compare your findings and comment. */

float anti_vars(float S0, float r, float sigma, float t, float X)
{
    float arr[N_simulate], w;
    for(int i = 0; i < N_simulate; i++)
    {
        w = sqrt(t) * dis(gen);
        arr[i] = (max(exp(-r*t)*(S0*exp(sigma*w + (r-pow(sigma,2.0)/2.0)*t) - X), 0) + max(exp(-r*t)*(S0*exp(-sigma*w + (r-pow(sigma,2.0)/2.0)*t) - X), 0))/2.0;
    }
    return arr_mean(arr, N_simulate);
}







/* --Question 4-- */

/* (a) For each integer number 𝑛 from 1 to 10, use 1,000 simulations of 𝑆𝑛 to estimate 𝐸(𝑆𝑛), where 𝑆𝑡 is a Geometric Brownian Motion process: 𝑆𝑡=𝑆0𝑒(𝜎𝑊𝑡+(𝑟−𝜎^2/2)𝑡), where 𝑟=0.055,𝜎=0.20,𝑆0=$88. Plot all of the above 𝐸(𝑆𝑛), for 𝑛 ranging from 1 to 10, in one graph. */
void ES_n(float r, float sigma, float S0)
{
    float ES[10];
    for(int n = 1; n <= 10; n++)
    {
        ES[n-1] = 0;
        for(int i = 0; i < N_simulate; i++)
        {
            ES[n-1] += S0*exp(sigma*dis(gen)*sqrt(n) + (r-pow(sigma,2.0)/2.0)*n);
        }
        ES[n-1] /= float(N_simulate);
    }
    for(int n = 1; n <= 10; n++)
    {
        cout << ES[n-1] << ", ";
    }
    cout << endl;
}


/* (b) Now simulate 3 paths of 𝑆𝑡 for 0≤𝑡≤10 (defined in part (a)) by dividing up the interval [0, 10] into 1,000 equal parts  */
void simulate_paths(float r, float sigma, float S0, float n)
{
    float S[int(n)][1000], w, dt = 0.01;
    for(int i = 0; i < n; i++)
    {
        S[i][0] = S0;
        cout << "[ " << S0;
        for(int j = 1; j < 1000; j++)
        {
            w = sqrt(dt) * dis(gen);
            S[i][j] = S[i][j-1] + r*S[i][j-1]*dt + sigma*S[i][j-1]*w;
            cout << ", " << S[i][j];
        }
        cout << "]" << endl;
    }
}


/* (c) Plot your data from parts (a) and (b) in one graph. */



/* (d) What would happen to the 𝐸(𝑆𝑛) graph if you increased 𝜎 from 20% to 30%? What would happen to the 3 plots of 𝑆𝑡 for 0≤𝑡≤10, if you increased 𝜎 from 20% to 30%? */




/* --Question 5-- */

/* (a) Write a code to compute prices of European Call options via Monte Carlo simulation of paths of the stock price process. Use Euler’s discretization scheme to discretize the SDE for the stock price process. The code should be generic: for any input of the 5 model parameters - 𝑆0 , 𝑇, 𝑋, 𝑟, 𝜎 − the output is the corresponding price of the European call option and the standard error of the estimate. */
void euler_call(float S0, float T, float X, float r, float sigma, float dt)
{
    float S, call[N_simulate], std=0;
    //dS=rSdt+sigmaSdW
    for(int i = 0; i < N_simulate; i++)
    {
        S = S0;
        for(float t = 0; t < T; t = t + dt)
        {
            S += r*S*dt + sigma*S*sqrt(dt)*dis(gen);
        }
        call[i] = exp(-r*T) * max(S-X,0);
    }
    float mean = arr_mean(call, N_simulate);
    cout << mean << ", ";
    for(int i = 0; i < N_simulate; i++)
    {
        std += pow(call[i]-mean,2.0);
    }
    cout << sqrt(std) << endl;
}
 
/* (b) Write a code to compute prices of European Call options via Monte Carlo simulation of paths of the stock price process. Use Milshtein’s discretization scheme to discretize the SDE for the stock price process. The code should be generic: for any input of the 5 model parameters - 𝑆0 , 𝑇, 𝑋, 𝑟, 𝜎 − the output is the corresponding price of the European call option and the standard error of the estimate. */
void milshtein_call(float S0, float T, float X, float r, float sigma, float dt)
{
    float S, call[N_simulate], std=0, w;
    for(int i = 0; i < N_simulate; i++)
    {
        S = S0;
        for(float t = 0; t < T; t = t + dt)
        {
            w = sqrt(dt)*dis(gen);
            S += r*S*dt + sigma*S*w + pow(sigma,2.0)*S*(pow(w,2.0)-dt)*0.5;
        }
        call[i] = exp(-r*T)*max(S-X,0);
    }
    double mean = arr_mean(call, N_simulate);
    cout << mean << ", ";
    for(int i = 0; i < N_simulate; i++)
    {
        std += pow(call[i]-mean,2.0);
    }
    cout << sqrt(std) << endl;
}
 
 /* (c) Write code to compute the prices of European Call options by using the Black-Scholes formula. Use the approximation of 𝑁(∙) described in Chapter 3. The code should be generic: for any input values of the 5 parameters - 𝑆0, 𝑇, 𝑋, 𝑟, 𝜎 - the output is the corresponding price of the European call option. */
double cdf_N(double x)
{
    double d_1 = 0.0498673470;
    double d_2 = 0.0211410061;
    double d_3 = 0.0032776263;
    double d_4 = 0.0000380036;
    double d_5 = 0.0000488906;
    double d_6 = 0.0000053830;
    if(x >= 0)
        return 1 - 0.5 * pow(1 + d_1*x + d_2*pow(x,2.0) + d_3*pow(x,3.0) + d_4*pow(x,4.0) + d_5*pow(x,5.0) + d_6*pow(x,6.0),-16.0);
    else
        return 1 - cdf_N(-x);
}
double bs_call_cdf(float S0, float r, float sigma, float t, float X)
{
    double d1 = (log(S0 / X) + (r + 0.5 * pow(sigma,2.0)) * t) / (sigma * sqrt(t));
    double d2 = d1 - sigma * sqrt(t);
    return S0 * cdf_N(d1) - X * exp(-r * t) * cdf_N(d2);
}


/* (d) Use the results of (a) to (c) to compare the two schemes of parts (a) and (b) by computing the following European call option prices: 𝑋 = 100, 𝜎 = 0.25, 𝑟 = 0.055, 𝑇 = 0.5, and the range [95, 104] for 𝑆0, with a step size of 1. */


 
 /* (e) Estimate the European call option’s greeks – delta (Δ), gamma (Γ), theta (𝜃), and vega (𝜈) - and graph them as functions of the initial stock price 𝑆0. Use 𝑋 = 100, 𝜎 = 0.25, 𝑟 = 0.055 and 𝑇 = 0.5 in your estimations. Use the range [95, 105] for 𝑆0 , with a step size of 1. You will have 4 different graphs for each of the 4 greeks. In all cases, 𝑑𝑡 (time-step) should be user-defined. Use 𝑑𝑡 = 0.05 as a default value. */
void greeks(float r, float sigma, float T, float X, float dt = 0.05)
{
    float delta[11], gamma[11], theta[11], vega[11];
    float S0_list[11] = {95,96,97,98,99,100,101,102,103,104,105};
    cout << "[" ;
    for(int i=0;i<=10;i++)
    {
        delta[i] = (bs_call_cdf(S0_list[i]+1,r,sigma,T,X) - bs_call_cdf(S0_list[i]-1,r,sigma,T,X))/2.0;
        gamma[i] = ((bs_call_cdf(S0_list[i]+2,r,sigma,T,X) - bs_call_cdf(S0_list[i],r,sigma,T,X))/2.0-(bs_call_cdf(S0_list[i],r,sigma,T,X) - bs_call_cdf(S0_list[i]-2,r,sigma,T,X))/2.0)/2.0;
        theta[i] = (bs_call_cdf(S0_list[i],r,sigma,T-0.1,X) - bs_call_cdf(S0_list[i],r,sigma,T+0.1,X))/0.2;
        vega[i] = (bs_call_cdf(S0_list[i],r,sigma+0.01,T,X) - bs_call_cdf(S0_list[i],r,sigma-0.01,T,X))/0.02;
        cout << "[ " << delta[i] << ", " << gamma[i]  << ", " << theta[i]  << ", " << vega[i] << " ]" << endl;
    }
    cout << "]" ;
}


/* --Question 6-- */
/* Consider the following 2-factor model for stock prices with stochastic volatility:𝑑𝑆𝑡 = 𝑟𝑆𝑡 𝑑𝑡 + √𝑉𝑡 𝑆𝑡 𝑑𝑊𝑡^1, 𝑑𝑉𝑡 = 𝛼(𝛽 − 𝑉𝑡)𝑑𝑡 + 𝜎√𝑉𝑡 𝑑𝑊𝑡^2 where the Brownian Motion processes above are correlated: 𝑑𝑊𝑡^1𝑑𝑊𝑡^2 = 𝜌𝑑𝑡, where the correlation 𝜌 is a constant in [−1,1]. Estimate the price of a European Call option (via Monte Carlo simulation) that has a strike price of 𝑋 and matures in 𝑇 years. Use the following default parameters of the model: 𝜌 = −0.6, 𝑟 = 0.055, 𝑆0 = $100, 𝑋 = $100, 𝑉0 = 0.05, 𝜎 = 0.42, 𝛼 = 5.8, 𝛽 = 0.0625, 𝑑𝑡 = 0.05, 𝑁 = 10,000. Use the Full Truncation, Partial Truncation, and Reflection methods, and provide 3 price estimates by using the tree methods.*/
void sto_vol(float rho, float r, float S0, float X, float V0, float sigma, float alpha, float beta,float dt, int N)
{
    double S_ft, V_ft, c_ft[N_simulate], S_pt, V_pt, c_pt[N_simulate], S_r=S0, V_r=V0, c_r[N_simulate], z1, z2;
    for(int j=0;j<N_simulate;j++)
    {
        //initiate S equals S0, V equals V0
        S_ft = S_pt = S_r = S0;
        V_ft = V_pt = V_r = V0;
        for(int i = 1; i < N; i++)
        {
            z1 = dis(gen);
            z2 = rho * z1 + sqrt(1-pow(rho,2.0))*dis(gen);
            // full truncation
            S_ft += r*S_ft*dt + sqrt(max(V_ft,0))*S_ft*sqrt(dt)*z1;
            V_ft += alpha*(beta-max(V_ft,0))*dt + sigma * sqrt(max(V_ft,0)*dt)*z2;
            // partial truncation
            S_pt += r*S_pt*dt + sqrt(max(V_pt,0))*S_pt*sqrt(dt)*z1;
            V_pt += alpha*(beta-V_pt)*dt + sigma * sqrt(max(V_pt,0)*dt)*z2;
            // reflection methods
            S_r += r*S_r*dt + sqrt(abs(V_r))*S_r*sqrt(dt)*z1;
            V_r = abs(V_r) + alpha*(beta-abs(V_r))*dt + sigma * sqrt(abs(V_r)*dt)*z2;
        }
        // save call option price max(S-X,0)
        c_ft[j] = max(S_ft-X,0);
        c_pt[j] = max(S_pt-X,0);
        c_r[j] = max(S_r-X,0);
    }
    cout << "Full Truncation: " << exp(-r*dt*N)*arr_mean(c_ft, N_simulate) << endl;
    cout << "Partial Truncation: " << exp(-r*dt*N)*arr_mean(c_pt, N_simulate) << endl;
    cout << "Reflection method: " << exp(-r*dt*N)*arr_mean(c_r, N_simulate) << endl;
}





/* --Question 7-- */
/* The objective of this exercise is to compare a sample of Pseudo-Random numbers with a sample of Quasi-Monte Carlo numbers of 𝑈𝑛𝑖𝑓𝑜𝑟𝑚[0,1]𝑥[0,1]: Use 2-dimensional Halton sequences to estimate the following integral: 𝐼=∫∫𝑒^{−𝑥𝑦}(𝑠𝑖𝑛6𝜋𝑥+𝑐𝑜𝑠^{1/3}2𝜋𝑦)𝑑𝑥𝑑𝑦 Default parameter values: N=10,000; (2,3) for bases. */
class HaltonSequence {
public:
    HaltonSequence(int base) : base(base) {}
    double next()
    {
        double f = 1, r = 0;
        int i = num;
        while (i > 0)
        {
            f = f / base;
            r = r + f * (i % base);
            i = i / base;
        }
        num++;
        return r;
    }
private:
    int base;
    int num = 0;
};

double func(double x, double y)
{
    if(cos(2 * M_PI * y) > 0)
        return exp(-x*y)*(sin(6 * M_PI * x) + pow(cos(2 * M_PI * y),1.0/3.0));
    else
        return exp(-x*y)*(sin(6 * M_PI * x) - pow(-cos(2 * M_PI * y),1.0/3.0));
}

double integral_7()
{
    HaltonSequence x(2), y(3);
    double ans = 0;
    x.next();
    y.next();
    for(int i = 0; i < N_simulate; i++)
    {
        ans += func(x.next(),y.next());
    }
    return ans/N_simulate;
}



int main()
{
    // Question 1
    cout << "Question 1: " << endl;
    cout << "1(a): " << endl;
    binomial(1000, 44, 0.64);
    cout << "1(b): " << endl;
    exponential(1/1.5, 1, 4);
    auto time_1 = chrono::high_resolution_clock::now();
    cout << "1(c):" << endl;
    box_muller(2500);
    auto time_2 = chrono::high_resolution_clock::now();
    cout << "1(d):" << endl;
    polar_marsaglia(2500);
    auto time_3 = chrono::high_resolution_clock::now();
    cout << "1(e)" << endl;
    chrono::duration<double> delta_1 = time_2 - time_1;
    chrono::duration<double> delta_2 = time_3 - time_2;
    cout << "The time of box_muller is  " << delta_1.count() << " s\n";
    cout << "The time of polar_marsaglia is  " << delta_2.count() << " s\n";
    
    
    // Question 2
    // 2(a)
    cout << "Question 2: " << endl;
    cout << "2(a):" << endl;
    calculate_AB(1);
    calculate_AB(3);
    calculate_AB(5);
    cout << B_anti_vars(5) << endl;
    
    
    // Question 3
    float S0=100, r=0.055, sigma=0.2;
    float T=5, X=100, dt=0.05;
    cout << "Question 3: " << endl;
    cout << "3(a): " << call(S0, r, sigma, T, X) << endl;
    cout << "3(b): " << bs_call(S0, r, sigma, T, X) << endl;
    cout << "3(c): " << anti_vars(S0, r, sigma, T, X) << endl;
    
    
    // Question 4
    cout << "Question 4: " << endl;
    cout << "4(a): ";
    sigma = 0.2;
    S0 = 88;
    ES_n(r, sigma, S0);
    cout << "4(b): " << endl;
    simulate_paths(r, sigma, S0, 3);
    cout << "4(c): " << endl;
    cout << "4(d): " << endl;
    sigma = 0.3;
    ES_n(r, sigma, S0);
    simulate_paths(r, sigma, S0, 3);
    
    
    // Question 5
    T = 0.5;
    sigma = 0.25;
    r = 0.055;
    dt = T/100;
    X = 100;
    cout << "Question 5: " << endl;
    cout << "The estimate " << " standard error " << endl;
    cout << "euler_call:" << endl;
    for(S0 = 95; S0 < 105; S0++)
    {
        euler_call(S0, T, X, r, sigma,dt);
    }
    cout << "milshtein_call:" << endl;
    for(S0 = 95; S0 < 105; S0++)
    {
        milshtein_call(S0, T, X, r, sigma,dt);
    }
    cout << "bs_call_cdf:" << endl;
    for(S0 = 95; S0 < 105; S0++)
    {
        cout << bs_call_cdf(S0, r, sigma, T, X) << ", ";
    }
    greeks(r, sigma, T, X, dt);
    
    
    // Question 6
    float rho=-0.6, V0=0.05, alpha=5.8, beta=0.0625;
    S0 = 100;
    r = 0.055;
    X = 100;
    sigma = 0.42;
    dt = 0.05;
    cout << "Question 6: " << endl;
    sto_vol(rho, r, S0, X, V0, sigma, alpha, beta, dt, N_simulate);
    
    
    //Question 7
    cout << "Question 7: " << endl;
    cout << integral_7();
    
    
    return 0;
}
