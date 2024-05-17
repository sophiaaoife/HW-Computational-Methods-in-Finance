//
//  main.cpp
//  HW2_2
//
//  Created by Yujie Yi on 2024/5/7.
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
// Compare the convergence rates of the four methods below by doing the following: Use the Binomial Method to price a 6-month American Put option with the following information: the risk-free interest rate is 5.5% per annum; the volatility is 25% per annum; the current stock price is $180; and the strike price is $170. Divide the time interval into 𝑛 parts to estimate the price of this option. Use 𝑛 = 20, 40, 80, 100, 200, 500, to estimate the price and draw all resulting prices in one graph, where the horizontal axis measures 𝑛, and the vertical one the price of the option. (a) Use the binomial method in which u=1/d, d=c− √c2 − 1, c = 1/2 (𝑒−𝑟∆ + 𝑒(𝑟+𝜎2)∆), 𝑝 = (𝑒𝑟∆ − 𝑑)/(𝑢 − 𝑑)   (b) Use the binomial method in which 𝑢 = 𝑒[(𝑟−𝜎2/2)∆+𝜎√∆], 𝑑 = 𝑒[(𝑟−𝜎2/2)∆−𝜎√∆], 𝑝 = 1/2

void binomial_1(int n, float T, float r, float sigma, float S0, float K)
{
    float S[n+1], C[n+1];
    float c = (exp(-r*T/n)+exp((r+pow(sigma,2))*T/n))/2;
    float d = c-sqrt(pow(c,2)-1), u = 1/d;
    float p = (exp(r*T/n)-d)/(u-d);
    for(int i = 0; i <= n; i++)
    {
        for(int j = 0; j <= n-i; j++)
        {
            S[j] = S0 * pow(u,n-i-j) * pow(d,j);
            if(i == 0)
                C[j] = max(K-S[j],0);
            else
                C[j] = max(max(K-S[j],0),exp(-r*T/n)*(p*C[j]+(1-p)*C[j+1]));
        }
    }
    cout << C[0] << ", ";
}

void binomial_2(int n, float T, float r, float sigma, float S0, float K)
{
    float S[n+1],C[n+1];
    float u = exp((r-pow(sigma,2)/2)*T/n+sigma*sqrt(T/n)), d = exp((r-pow(sigma,2)/2)*T/n-sigma*sqrt(T/n));
    float p = 0.5; //cannot write 1/2
    for(int i = 0; i <= n; i++)
    {
        for(int j = 0; j <= n-i; j++)
        {
            S[j] = S0 * pow(u,n-i-j) * pow(d,j);
            if(i == 0)
                C[j] = max(K-S[j],0);
            else
                C[j] = max(max(K-S[j],0),exp(-r*T/n)*(p*C[j]+(1-p)*C[j+1]));
        }
    }
    cout << C[0] << ", ";
}

// Question 2
// Consider the following information on the stock of a company and American put options on it: 𝑆0 = $180, 𝑋 = $170, 𝑟 = 0.055, 𝜎 = 0.25, 𝑇 = 6 months, 𝜇 = 0.15. Using the CRR Binomial tree method, estimate the following and draw their graphs:

// (i) Delta of the put option as a function of 𝑆0, for 𝑆0 ranging from $170 to $190, in increments of $2.

float CRR(float T, float r, float sigma, float S0, float K, int n = 100)
{
    // use this function to calculate the CRR binomial tree method put option price
    float S[n+1],C[n+1];
    float u = exp(sigma*sqrt(T/n)), d = exp(-sigma*sqrt(T/n));
    float p = (exp(r*T/n)-d)/(u-d);
    for(int i = 0; i <= n; i++)
    {
        for(int j = 0; j <= n-i; j++)
        {
            S[j] = S0 * pow(u,n-i-j) * pow(d,j);
            if(i == 0)
                C[j] = max(K-S[j],0);
            else
                C[j] = max(max(K-S[j],0),exp(-r*T/n)*(p*C[j]+(1-p)*C[j+1]));
        }
    }
    return C[0];
}

void delta_S(float T, float r, float sigma, float S1, float S2, float K, float increments)
{
    int N = (S2-S1)/increments;
    float delta[N+1];
    for(int i = 0; i <= N; i++)
    {
        delta[i] = (CRR(T,r,sigma,S1+increments*i+1,K) - CRR(T,r,sigma,S1+increments*i-1,K))/2;
        cout << delta[i] << ", ";
    }
    cout << endl;
}


// (ii) Delta of the put option, as a function of T (time to expiration), T ranging from 0 to 0.18 in increments of 0.003.
void delta_T(float T, float r, float sigma, float S0, float K, float increments)
{
    int N = T/increments;
    float delta[N+1];
    for(int i = 0; i <= N; i++)
    {
        delta[i] = (CRR(T+increments*i,r,sigma,S0+1,K) - CRR(T+increments*i,r,sigma,S0-1,K))/2;
        cout << delta[i] << ", ";
    }
    cout << endl;
}

// (iii) Theta of the put option, as a function of T (time to expiration), T ranging from 0 to 0.18 in increments of 0.003.
void theta_T(float T, float r, float sigma, float S0, float K, float increments)
{
    int N = T/increments;
    float theta[N+1];
    for(int i = 0; i <= N; i++)
    {
        // if i == 0 we can only calculate theta on one side
        if(i == 0)
            theta[i] = (CRR(T+increments*i,r,sigma,S0,K) - CRR(T+increments*i+0.001,r,sigma,S0,K))/0.001;
        else
            theta[i] = (CRR(T+increments*i-0.001,r,sigma,S0,K) - CRR(T+increments*i+0.001,r,sigma,S0,K))/0.002;
        cout << theta[i] << ", ";
    }
    cout << endl;
}


// (iv) Vega of the put option, as a function of 𝑆0, for 𝑆0 ranging from $170 to $190, in increments of $2.
void vega_S(float T, float r, float sigma, float S1, float S2, float K, float increments)
{
    int N = (S2-S1)/increments;
    float vega[N+1];
    for(int i = 0; i <= N; i++)
    {
        vega[i] = (CRR(T,r,sigma+0.01,S1+increments*i,K) - CRR(T,r,sigma-0.01,S1+increments*i,K))/0.02;
        cout << vega[i] << ", ";
    }
    cout << endl;
}



// Question 3
// Compare the convergence rates of the two methods, (a) and (b), described below, by doing the following: Use the Trinomial-tree method to price a 6-month American put option with the following information: the risk-free interest rate is 5.5% per annum, the volatility is 25% per annum, the current stock price is $180, and the strike price is $170. Divide the time interval into 𝑛 equal parts to estimate the option price. Use n = 20, 40, 70, 80, 100, 200, 500; to estimate option prices and draw them all in one graph, where the horizontal axis measures 𝑛, and the vertical one measures option price. The two methods are in (a) and (b) below

// (a) Use the Trinomial-tree method applied to the stock price-process (𝑆𝑡) in which 𝑢 = 1/𝑑 , 𝑑 = 𝑒−𝜎√3∆, 𝑝𝑑 = (𝑟∆(1−𝑢)+(𝑟∆)2+𝜎2∆)/(𝑢−𝑑)(1−𝑑) , 𝑝𝑢 = (𝑟∆(1−𝑑)+(𝑟∆)2+𝜎2∆)(𝑢−𝑑)(𝑢−1) , 𝑝𝑚 = 1 − 𝑝𝑢 − 𝑝𝑑
void trinomial_1(int n, float T, float r, float sigma, float S0, float K)
{
    float S[2*n+1], P[2*n+1];
    float d = exp(-sigma*sqrt(3*T/n)), u = 1/d;
    float pd = (r*T/n*(1-u)+pow(r*T/n,2)+pow(sigma,2)*T/n)/(u-d)/(1-d);
    float pu = (r*T/n*(1-d)+pow(r*T/n,2)+pow(sigma,2)*T/n)/(u-d)/(u-1), pm = 1-pu-pd;
    for(int i = 0; i <= n; i++)
    {
        for(int j = 0; j <= 2*(n-i); j++)
        {
            S[j] = S0 * pow(u,n-i-j);
            if(i == 0)
                P[j] = max(K-S[j],0);
            else
                P[j] = max(max(K-S[j],0),exp(-r*T/n)*(pu*P[j]+pm*P[j+1]+pd*P[j+2]));
        }
    }
    cout << P[0] << ", ";
}


// (b) Use the Trinomial-tree method applied to the Log-stock price-process (𝑋𝑡) in which ∆𝑋𝑢 = 𝜎√3∆, ∆𝑋𝑑 = −𝜎√3∆, 𝑝𝑑 = 1/2 ((𝜎2∆+(𝑟−𝜎2/2)^2∆^2)/∆𝑋𝑢2 − (𝑟−𝜎2/2)∆/∆𝑋𝑢), 𝑝𝑢 = ((𝜎2∆+(𝑟−𝜎2/2)^2∆^2)/∆𝑋𝑢2 + (𝑟−𝜎2/2)∆/∆𝑋𝑢), 𝑝𝑚 = 1 − 𝑝𝑢 − 𝑝𝑑
void trinomial_2(int n, float T, float r, float sigma, float S0, float K)
{
    float X[2*n+1], P[2*n+1];
    float u = sigma*sqrt(3*T/n);
    float pd = ((pow(sigma,2)*T/n+pow(r-pow(sigma,2)/2,2)*pow(T/n,2))/pow(u,2)-(r-pow(sigma,2)/2)*T/n/u)/2;
    float pu = ((pow(sigma,2)*T/n+pow(r-pow(sigma,2)/2,2)*pow(T/n,2))/pow(u,2)+(r-pow(sigma,2)/2)*T/n/u)/2;
    float pm = 1-pu-pd;
    //cout << u << ", "<< d << ", "<< pu << ", "<< pd << ", "<< pm << ", " << endl;
    for(int i = 0; i <= n; i++)
    {
        //cout << i << ": ";
        for(int j = 0; j <= 2*(n-i); j++)
        {
            X[j] = log(S0) + (n-i-j) * u;
            if(i == 0)
                P[j] = max(K-exp(X[j]),0);
            else
                P[j] = max(max(K-exp(X[j]),0),exp(-r*T/n)*(pu*P[j]+pm*P[j+1]+pd*P[j+2]));
            //cout << P[j] << ", ";
        }
        //cout << endl;
    }
    cout << P[0] << ", ";
}


// Question 4
// Consider the following information on the stock of company XYZ: The current stock price is $180, and the volatility of the stock price is 𝜎 = 25% per annum. Assume the prevailing risk-free rate is 𝑟 = 5.5% per annum. Use the following method to price the specified option:

// (a) Use the LSMC method with N=100,000 paths simulations (50,000 plus 50,000 antithetic variates) and a time step of ∆= 1/√𝑁 to price an American Put option with strike price of 𝑋 = $170 and maturity of 0.5-years and 1.5-years. Use the first 𝑘 of the Laguerre Polynomials for 𝑘 = 2, 3, 4, 5. (That is, you will compute 8 prices here). Compare the prices for the 4 cases, 𝑘 = 2, 3, 4, 5 and comment on the choice of k.
double Laguerre(int k, double x)
{
    // Compute the Laguerre polynomial value for a given x based on order k
    switch(k)
    {
        case 0: return exp(-x/2);
        case 1: return exp(-x/2)*(1-x);
        case 2: return exp(-x/2)*(1-2*x+pow(x,2)/2.0);
        case 3: return exp(-x/2)*(1-3*x+3*pow(x,2)/2.0-pow(x,3)/6.0);
        case 4: return exp(-x/2)*(1-4*x+3*pow(x,2)-2*pow(x,3)/3.0+pow(x,4)/24.0);
    }
    return 0;
}

SparseMat simulate_path(vector<vector<double>>& S, float S0, int N, int n, float r, float sigma, float dt, float K)
{
    // antithetic variates
    float w;
    for(int i = 0; i < N/2; i++){
        S[0][2*i] = S0;
        S[0][2*i+1] = S0;
        for(int j = 1; j < n; j++){
            w = sqrt(dt) * dis(gen);
            S[j][2*i] = S[j-1][2*i] + r*S[j-1][2*i]*dt + sigma*S[j-1][2*i]*w;
            S[j][2*i+1] = S[j-1][2*i+1] + r*S[j-1][2*i+1]*dt - sigma*S[j-1][2*i+1]*w;
        }
    }
    SparseMat EV(N,n);
    EV.setZero();
    for(int i = 0; i < N; i++){
        for(int j = 0; j < n; j++){
            if(K-S[j][i]>0){EV.coeffRef(i,j) = K-S[j][i];}
        }
    }
    return EV;
}

float LSMC(float T, float r, float sigma, double S0, double K, int k, double(*base_func)(int,double), int N = 10000)
{
    // Compute LSMC for American options using regression
    float dt = 1/sqrt(N);
    int n = int(T/dt)+1; // calculate how many timestamp
    
    // use Sparse Matrix to save space
    SparseMat Index(n, N);
    SparseMat P(n, N);
    MatrixXd A(k, k);
    // initiate as 0
    Index.setZero();
    P.setZero();
    A.setZero();
    vector<vector<double>> S(n, vector<double>(N, 0.0));
    VectorXd b = VectorXd::Zero(k), a = VectorXd::Zero(k);
    
    // 1. simulate N paths
    SparseMat EV(N, n);
    EV = simulate_path(S, S0, N, n, r, sigma, dt, K);
    EV.makeCompressed();
    
    // set the initial index
    int count = 0;
    for (SparseMatrix<double>::InnerIterator it(EV, n-1); it; ++it) {
        if(EV.coeffRef(it.row(),n-1)<=0){cout << "error!!!" << endl;}
        Index.coeffRef(n-1,it.row()) = 1;
    }
    
    // 2. Calculate continuation values and make exercise decisions
    int j = 0;
    for(int i = n-1; i >= 1; i--)
    {
        A.setZero();
        b.setZero();
        // calculate parameter a s.t. A'a=b
        for (SparseMatrix<double>::InnerIterator it(EV, i-1); it; ++it) {
            j = int(it.row());
            for(int m = 0; m < k; m++){
                for(int q = 0; q < k; q++){A(m,q) += base_func(m,S[i-1][j]/K) * base_func(q,S[i-1][j]/K);}
                for(int q = i; q < n; q++){b(m) += Index.coeffRef(q,j)*exp(-r*(q-i+1)*dt)*max(1-S[q][j]/K,0) * base_func(m,S[i-1][j]/K);}
            }
        }
        a = A.lu().solve(b);
        
        // compare EV and CV
        for (SparseMatrix<double>::InnerIterator it(EV, i-1); it; ++it) {
            j = int(it.row());
            float CV = 0;
            for(int m = 0; m < k; m++){CV += a[m]*base_func(m,S[i-1][j]/K);}
            if(EV.coeffRef(j,i-1)/K > CV){
                Index.coeffRef(i-1,j) = 1;
                for(int q = i; q < n; q++){Index.coeffRef(q,j) = 0;}
                P.coeffRef(i-1,j) = exp(-r*(i-1)*dt) * EV.coeffRef(j,i-1);
                count += 1;
            }
        }
    }

    // Last Calculate and return the average payoff
    float mean = 0;
    Index.makeCompressed();
    for (int k = 0; k < Index.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(Index, k); it; ++it) {
            j = int(it.row());
            if(Index.coeffRef(j,k) == 1)
            {
                mean += P.coeffRef(j,k);
            }
        }
    }
    return mean/N;
}

// (b) Use the LSMC method with N=100,000 paths simulations (50,000 plus 50,000 antithetic variates) and a time step of ∆= 1/√𝑁 to price an American Put option with strike price of 𝑋 = $170 and maturity of 0.5-years and 1.5-years. Use the first 𝑘 of the Hermite Polynomials for 𝑘 = 2, 3, 4, 5. (That is, you will compute 8 prices here). Compare the prices for the 4 cases, 𝑘 = 2, 3, 4, 5 and comment on the choice of k.
double Hermite(int k, double x)
{
    switch(k)
    {
        case 0: return 1;
        case 1: return 2*x;
        case 2: return 4*pow(x,2)-2;
        case 3: return 8*pow(x,3)-12*x;
        case 4: return 16*pow(x,4)-56*pow(x,2)+16;
    }
    return 0;
}

// (c) Use the LSMC method with N=100,000 paths simulations (50,000 plus 50,000 antithetic variates) and a time step of ∆= 1/√𝑁 to price an American Put option with strike price of 𝑋 = $170 and maturity of 0.5-years and 1.5-years. Use the first 𝑘 of the Simple Monomials for 𝑘 = 2, 3, 4, 5. (That is, you will compute 8 prices here). Compare the prices for the 4 cases, 𝑘 = 2, 3, 4, 5 and comment on the choice of k.
double Monomials(int k, double x)
{
    switch(k)
    {
        case 0: return 1;
        case 1: return x;
        case 2: return pow(x,2);
        case 3: return pow(x,3);
        case 4: return pow(x,4);
    }
    return 0;
}

// (d) Compare all your findings above and commen



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

// Question 5
// Consider the following information on the stock of company XYZ: The volatility of the stock price is 𝜎 = 25% per annum. Assume the prevailing risk-free rate is 𝑟 = 5.5% per annum. Use the 𝑋 = 𝑙𝑛(𝑆) transformation of the Black-Scholes PDE, and 𝛥𝑡 = 0.002, with 𝛥𝑋 = 𝜎√𝛥𝑡; then with 𝛥𝑋 = 𝜎√3𝛥𝑡; then with 𝛥𝑋 = 𝜎√4𝛥𝑡, and a uniform grid (on X) to price an American Put option with strike price of 𝑋 = $170, expiration of 6 months and current stock prices ranging from $170 to $190; using the specified methods below:

// (a) Explicit Finite-Difference method,
MatrixXd EFD(float dt, float sigma, float r, float dx, float S1, float S2, float K, float T)
{
    int n = T/dt, N = int(log(S2/S1)/dx+1);
    float Pu = dt/2*(pow(sigma/dx,2)+(r-pow(sigma,2)/2)/dx);
    float Pd = dt/2*(pow(sigma/dx,2)-(r-pow(sigma,2)/2)/dx);
    float Pm = 1-dt*pow(sigma/dx,2)-r*dt;
    MatrixXd A(N+1,N+1);
    VectorXd P(N+1), B(N+1);
    
    // initiate A
    VectorXd p(3);
    p << Pu,Pm,Pd;
    A.row(0).segment(0,3) = p;
    A.row(N).segment(N-2,3) = p;
    for(int i = 1; i <= N-1; i++){A.row(i).segment(i-1,3) = p;}
    
    //initiate the boudary value
    for(int i = 0; i <= N; i++){P(i) = max(K - S1*exp((N-i)*dx),0);}
    
    // iterate to get Put price(t=0)
    B(N) = S1*(exp(dx)-1);
    for(int i = n-1; i >= 0; i--)
    {
        P = A * P + B;
        // American condition
        for(int j = 0; j <= N; j++){P(j) = max(P(j), max(K - S1*exp((N-j)*dx),0));}
    }
    return P;
}

// (b) Implicit Finite-Difference method,
MatrixXd IFD(float dt, float sigma, float r, float dx, float S1, float S2, float K, float T)
{
    int n = T/dt, N = int(log(S2/S1)/dx+1);
    float Pu = -dt/2*(pow(sigma/dx,2)+(r-pow(sigma,2)/2)/dx);
    float Pd = -dt/2*(pow(sigma/dx,2)-(r-pow(sigma,2)/2)/dx);
    float Pm = 1+dt*pow(sigma/dx,2)+r*dt;
    MatrixXd A(N+1,N+1);
    VectorXd P(N+1), B(N+1);
    
    // initiate A
    VectorXd p(3);
    p << Pu,Pm,Pd;
    A(0,0) = 1;
    A(0,1) = -1;
    A(N,N-1) = 1;
    A(N,N) = -1;
    for(int i = 1; i <= N-1; i++){A.row(i).segment(i-1,3) = p;}
    
    //initiate the boudary value
    for(int i = 0; i <= N; i++){P(i) = max(K - S1*exp((N-i)*dx),0);}
    
    // iterate to get Put price(t=0)
    B(N) = - S1*(exp(dx)-1);
    for(int i = n-1; i >= 0; i--)
    {
        B.segment(1,N-1) = P.segment(1,N-1);
        P = A.colPivHouseholderQr().solve(B);
        // American condition
        for(int j = 0; j <= N; j++){P(j) = max(P(j), max(K - S1*exp((N-j)*dx),0));}
    }
    return P;
}

// (c) Crank-Nicolson Finite-Difference method.
MatrixXd CN_FD(float dt, float sigma, float r, float dx, float S1, float S2, float K, float T)
{
    int n = T/dt, N = int(log(S2/S1)/dx+1);
    float Pu = -dt*(pow(sigma/dx,2)+(r-pow(sigma,2)/2)/dx)/4;
    float Pd = -dt*(pow(sigma/dx,2)-(r-pow(sigma,2)/2)/dx)/4;
    float Pm = 1+dt*pow(sigma/dx,2)/2+r*dt/2;
    MatrixXd A(N+1,N+1);
    VectorXd P(N+1), B(N+1);
    
    // initiate A
    VectorXd p(3);
    p << Pu,Pm,Pd;
    A(0,0) = 1;
    A(0,1) = -1;
    A(N,N-1) = 1;
    A(N,N) = -1;
    for(int i = 1; i <= N-1; i++){A.row(i).segment(i-1,3) = p;}
    
    //initiate the boudary value
    for(int i = 0; i <= N; i++){P(i) = max(K - S1*exp((N-i)*dx),0);}
    
    // iterate to get Put price(t=0)
    B(N) = - S1*(exp(dx)-1);
    for(int i = n-1; i >= 0; i--)
    {
        B.segment(1,N-1) = - Pu * P.segment(0,N-1);
        B.segment(1,N-1) += - (Pm-2) * P.segment(1,N-1);
        B.segment(1,N-1) += - Pd * P.segment(2,N-1);
        P = A.colPivHouseholderQr().solve(B);
        // American condition
        for(int j = 0; j <= N; j++){P(j) = max(P(j), max(K - S1*exp((N-j)*dx),0));}
    }
    return P;
}





vector<double> elementWiseMultiply(const vector<double>& v1, const vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw invalid_argument("Vectors must be of the same size.");
    }

    vector<double> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] * v2[i];
    }
    return result;
}


// Question 6
// Consider the following information on the stock of company XYZ: The volatility of the stock price is 𝜎 = 25% per annum. Assume the prevailing risk-free rate is 𝑟 = 5.5% per annum. Use the Black-Scholes PDE (for S) to price American Put options with strike prices of 𝐾 = $170, expiration of 6 months and current stock prices for a range from $170 to $190; using the specified methods below: Choose 𝛥𝑡 = 0.002, with 𝛥𝑆 = 0.5, or with 𝛥𝑆 = 1.

// (a) Explicit Finite-Difference method,
MatrixXd EFD_S(float dt, float sigma, float r, float dS, float S1, float S2, float K, float T)
{
    int n = T/dt, N = int((S2-S1)/dS/2);
    float S0 = (S1+S2)/2;
    MatrixXd Pumd(2*N+1,3);
    VectorXd P(2*N+1), B(2*N+1);
    Pumd.setZero();
    P.setZero();
    B.setZero();
    
    // initiate A
    MatrixXd A = MatrixXd::Identity(2*N+1, 2*N+1);
    
    // initiate Pumd
    for(int i = N; i >= -N; i--)
    {
        Pumd.row(N-i) << dt/2*(pow(sigma*(S0+i*dS)/dS,2)+r*(S0+i*dS)/dS), 1-dt*pow(sigma*(S0+i*dS)/dS,2)+r*dt, dt/2*(pow(sigma*(S0+i*dS)/dS,2)-r*(S0+i*dS)/dS);
        if(i!=N and i!=-N){A.row(N-i).segment(N-i-1,3) = Pumd.row(N-i);}
    }
    A.row(0).segment(0,3) = Pumd.row(1);
    A.row(2*N-1).segment(2*N-2,3) = Pumd.row(2*N-1);
    
    //initiate the boudary value
    for(int i = 0; i <= 2*N; i++){P(i) = max(K - (S0+(N-i)*dS),0);}
    
    // iterate to get Put price(t=0)
    B(2*N) = dS;
    //Pumd.col(1).array() -= 2;
    for(int i = n-1; i >= 0; i--)
    {
        P = A.lu().solve(B);
        // American condition
        for(int j = 0; j <= 2*N; j++){P(j) = max(P(j), max(K - (S0+(N-j)*dS),0));}
    }
    return P;
}


// (b) Implicit Finite-Difference method,
MatrixXd IFD_S(float dt, float sigma, float r, float dS, float S1, float S2, float K, float T)
{
    int n = T/dt, N = int((S2-S1)/dS/2);
    float S0 = (S1+S2)/2;
    MatrixXd Pumd(2*N+1,3);
    MatrixXd A(2*N+1,2*N+1);
    VectorXd P(2*N+1), B(2*N+1);
    Pumd.setZero();
    A.setZero();
    P.setZero();
    B.setZero();
    
    // initiate A
    VectorXd p(3);
    A(0,0) = 1;
    A(0,1) = -1;
    A(2*N,2*N-1) = 1;
    A(2*N,2*N) = -1;
    for(int i = N-1; i >= -N+1; i--)
    {
        Pumd.row(N-i) << -dt/2*(pow(sigma*(S0+i*dS)/dS,2)+r*(S0+i*dS)/dS), 1+dt*pow(sigma*(S0+i*dS)/dS,2)+r*dt, -dt/2*(pow(sigma*(S0+i*dS)/dS,2)-r*(S0+i*dS)/dS);
        A.row(N-i).segment(N-i-1,3) = Pumd.row(N-i);
    }
    
    //initiate the boudary value
    for(int i = 0; i <= 2*N; i++){P(i) = max(K - (S0+(N-i)*dS),0);}
    //printAsPythonList(P);
    // iterate to get Put price(t=0)
    B(2*N) = - dS;
    for(int i = n-1; i >= 0; i--)
    {
        B.segment(1,2*N-1) = P.segment(1,2*N-1);
        P = A.lu().solve(B);
        // American condition
        for(int j = 0; j <= 2*N; j++){P(j) = max(P(j), max(K - (S0+(N-j)*dS),0));}
    }
    return P;
}


// (c) Crank-Nicolson Finite-Difference method.
MatrixXd CN_FD_S(float dt, float sigma, float r, float dS, float S1, float S2, float K, float T)
{
    int n = T/dt, N = int((S2-S1)/dS/2);
    float S0 = (S1+S2)/2;
    
    MatrixXd Pumd(2*N+1,3);
    MatrixXd A(2*N+1,2*N+1);
    VectorXd B(2*N+1), P(2*N+1);
    Pumd.setZero();
    A.setZero();
    P.setZero();
    B.setZero();
    
    // initiate A
    VectorXd p(3);
    A(0,0) = 1;
    A(0,1) = -1;
    A(2*N,2*N-1) = 1;
    A(2*N,2*N) = -1;
    
    for(int i = N; i >= -N; i--)
    {
        Pumd.row(N-i) << -dt/4*(pow(sigma*(S0+i*dS)/dS,2)+r*(S0+i*dS)/dS), 1+dt/2*(pow(sigma*(S0+i*dS)/dS,2)+r), -dt/4*(pow(sigma*(S0+i*dS)/dS,2)-r*(S0+i*dS)/dS);
        if(i!=N and i!=-N){A.row(N-i).segment(N-i-1,3) = Pumd.row(N-i);}
    }
    
    //initiate the boudary value
    for(int i = 0; i <= 2*N; i++){P(i) = max(K - (S0+(N-i)*dS),0);}
    
    // iterate to get Put price(t=0)
    B(2*N) = - dS;
    Pumd.col(1).array() -= 2;
    for(int i = n-1; i >= 0; i--)
    {
        B.segment(1,2*N-1) = - (Pumd.col(0).segment(1,2*N-1).array() * P.segment(0,2*N-1).array()).matrix();
        B.segment(1,2*N-1) += - (Pumd.col(1).segment(1,2*N-1).array() * P.segment(1,2*N-1).array()).matrix();
        B.segment(1,2*N-1) += - (Pumd.col(2).segment(1,2*N-1).array() * P.segment(2,2*N-1).array()).matrix();
        P = A.lu().solve(B);
        // American condition
        for(int j = 0; j <= 2*N; j++){P(j) = max(P(j), max(K - (S0+(N-j)*dS),0));}
    }
    return P;
}



int main()
{
    
    // Q1
    int n_list[] = {20, 40, 80, 100, 200, 500};
    float T = 0.5, r = 0.055, sigma = 0.25, S0 = 180, K = 170;
     /*cout << "1(a): ";
    for(int i = 0; i < 6; i++)
    {
        binomial_1(n_list[i], T, r, sigma, S0, K);
    }
    cout << endl << "1(b): ";
    for(int i = 0; i < 6; i++)
    {
        binomial_2(n_list[i], T, r, sigma, S0, K);
    }
    
    // Q2
    float S1 = 170, S2 = 190;
    cout << endl << "2(a): ";
    delta_S(T, r, sigma, S1, S2, K, 2);
    cout << "2(b): ";
    delta_T(0.18, r, sigma, S0, K, 0.003);
    cout << "2(c): ";
    theta_T(0.18, r, sigma, S0, K, 0.003);
    cout << "2(d): ";
    vega_S(T, r, sigma, S1, S2, K, 2);
    
    //Q3
    cout << "3(a): ";
    for(int i = 0; i < 6; i++)
    {
        trinomial_1(n_list[i], T, r, sigma, S0, K);
    }
    cout << endl << "3(b): ";
    for(int i = 0; i < 6; i++)
    {
        trinomial_2(n_list[i], T, r, sigma, S0, K);
    }
    
    // Q4
    
    cout << "4(a): ";
    for(float T = 1.5; T < 2; T+=1)
    {
        for(int n = 2; n <= 5; n++)
        {
            cout << T  << ", " << n << ": " << LSMC(T, r, sigma, S0, K, n, Laguerre)<<endl;
        }
    }
    cout << "4(b): ";
    for(float T = 0.5; T < 2; T+=1)
    {
        for(int n = 2; n <= 5; n++)
        {
            cout << T  << ", " << n << ": " << LSMC(T, r, sigma, S0, K, n, Hermite)<<endl;
        }
    }
    cout << "4(c): ";
    for(float T = 0.5; T < 2; T+=1)
    {
        for(int n = 2; n <= 5; n++)
        {
            cout << T  << ", " << n << ": " << LSMC(T, r, sigma, S0, K, n, Monomials)<<endl;
        }
    }*/
    //Q5
    cout << "5(a): ";
    float dt = 0.002, S1 = 170, S2 = 190;
    /*printAsPythonList(EFD(dt, sigma, r, sigma*sqrt(dt), S1, S2, K, T));
    printAsPythonList(EFD(dt, sigma, r, sigma*sqrt(3*dt), S1, S2, K, T));
    printAsPythonList(EFD(dt, sigma, r, sigma*sqrt(4*dt), S1, S2, K, T));
    cout << "5(b): ";
    printAsPythonList(IFD(dt, sigma, r, sigma*sqrt(dt), S1, S2, K, T));
    printAsPythonList(IFD(dt, sigma, r, sigma*sqrt(3*dt), S1, S2, K, T));
    printAsPythonList(IFD(dt, sigma, r, sigma*sqrt(4*dt), S1, S2, K, T));
    cout << "5(c): ";
    printAsPythonList(CN_FD(dt, sigma, r, sigma*sqrt(dt), S1, S2, K, T));
    printAsPythonList(CN_FD(dt, sigma, r, sigma*sqrt(3*dt), S1, S2, K, T));
    printAsPythonList(CN_FD(dt, sigma, r, sigma*sqrt(4*dt), S1, S2, K, T));*/
    
    //Q6
    cout << "6(a): ";
    printAsPythonList(EFD_S(dt, sigma, r, 0.5, S1, S2, K, T));
    printAsPythonList(EFD_S(dt, sigma, r, 1, S1, S2, K, T));
    cout << "6(b): ";
    //printAsPythonList(IFD_S(dt, sigma, r, 0.5, S1, S2, K, T));
    //printAsPythonList(IFD_S(dt, sigma, r, 1, S1, S2, K, T));
    cout << "6(c): ";
    //printAsPythonList(CN_FD_S(dt, sigma, r, 0.5, S1, S2, K, T));
    //printAsPythonList(CN_FD_S(dt, sigma, r, 1, S1, S2, K, T));
    return 0;
}
