//
//  main.cpp
//  endterm_part2
//
//  Created by Yitian Fu on 12/8/17.
//  Copyright © 2017 Yitian Fu. All rights reserved.
//


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "normdist.h"
#define UNIT_STEP(x) ((x)>0?(1):(0))
using namespace std;

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility, barrier_price;
int no_of_trials, no_of_discrete_barriers;


double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time){
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time) {  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

double down_and_in_call_option_theory(const double& S,       // underlying price
                                      const double& K,       // strike (exercise) price,
                                      const double& r,       // interest rate
                                      const double& sigma,   // volatility
                                      const double& time,  //time to maturity
                                      const double& H)     //barrier price
{
    
    double lamda=(r+(pow(sigma, 2))/2)/pow(sigma, 2);
    double y=(log((H*H)/(S*K)))/(sigma*sqrt(time))+lamda*sigma*sqrt(time);
    return S*pow(H/S,2*lamda)*N(y)-K*exp(-r*time)*pow(H/S,2*lamda-2)*N(y-sigma*sqrt(time));
}

double down_and_out_call_option_theory(const double& S,       // underlying price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time,  //time to maturity
                                       const double& H)     //barrier price
{
    double lamda=(r+(pow(sigma, 2))/2)/pow(sigma, 2);
    double x1= log(S/H)/(sigma*sqrt(time))+lamda*sigma*sqrt(time);
    double y1= log(H/S)/(sigma*sqrt(time))+lamda*sigma*sqrt(time);
    return S*N(x1)-K*exp(-r*time)*N(x1-sigma*sqrt(time))-S*pow(H/S,2*lamda)*N(y1)+K*exp(-r*time)*pow(H/S,2*lamda-2)*N(y1-sigma*sqrt(time));
}

double down_and_in_put_option_theory(const double& S,       // underlying price
                                     const double& K,       // strike (exercise) price,
                                     const double& r,       // interest rate
                                     const double& sigma,   // volatility
                                     const double& time,  //time to maturity
                                     const double& H)     //barrier price
{
    double lamda=(r+(pow(sigma, 2))/2)/pow(sigma, 2);
    double y=(log((H*H)/(S*K)))/(sigma*sqrt(time))+lamda*sigma*sqrt(time);
    double x1= log(S/H)/(sigma*sqrt(time))+lamda*sigma*sqrt(time);
    double y1= log(H/S)/(sigma*sqrt(time))+lamda*sigma*sqrt(time);
    return -S*N(-x1)+K*exp(-r*time)*N(-x1+sigma*sqrt(time))+S*pow(H/S,2*lamda)*(N(y)-N(y1))-K*exp(-r*time)*pow(H/S,2*lamda-2)*(N(y-sigma*sqrt(time))-N(y1-sigma*sqrt(time)));
}

double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
};



double max(double a, double b) {
    return (b < a )? a:b;
}

double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}


double discrete_probability_stock_path_not_hit_barrier(const double &current_stock_price){
    
    if (current_stock_price <= barrier_price)
        
    {
        return 1.0;
        
    }
    else
    {
    double probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant = 1.0;
    for ( int j=1; j<= no_of_discrete_barriers; j++)
    {
        double mean_at_sampling_instant_j =initial_stock_price+
        ( ((double) j)/((double) no_of_discrete_barriers)*(current_stock_price-initial_stock_price) );
        double variance_at_sampling_instant_j =( ((double) j)/((double) no_of_discrete_barriers) )*expiration_time*
        (1.0 - ((double) j)/((double) no_of_discrete_barriers));
        probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *=
        (1.0 - N((barrier_price - mean_at_sampling_instant_j)/sqrt(variance_at_sampling_instant_j)));
    
    }
        
      return probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant;
    }
    
}

int main (int argc, char* argv[])
{
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%lf", &risk_free_rate);
    sscanf (argv[3], "%lf", &volatility);
    sscanf (argv[4], "%lf", &initial_stock_price);
    sscanf (argv[5], "%lf", &strike_price);
    sscanf (argv[6], "%d", &no_of_trials);
    sscanf (argv[7], "%d", &no_of_discrete_barriers);
    sscanf (argv[8], "%lf", &barrier_price);
    
    
    cout << "--------------------------------" << endl;
    cout << "European Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of Discrete Barriers = " << no_of_discrete_barriers<<endl;
    cout << "Barrier Price =  " << barrier_price <<endl;
    cout << "--------------------------------" << endl;
    
    double discrete_barrier_call_option_price_via_simulation = 0.0;
    double discrete_barrier_put_option_price_via_simulation = 0.0;
    double discrete_barrier_call_option_price_adjust =0.0;
    double discrete_barrier_put_option_price_adjust =0.0;
    
    
    double deltaT = expiration_time/(double)no_of_discrete_barriers;
    double R = (risk_free_rate - 0.5*pow(volatility,2))*deltaT;
    double SD = volatility*sqrt(deltaT);
    
    
    
    for ( int i =0 ;i<no_of_trials; i++)
    {
        double S1=initial_stock_price;
        double S2=initial_stock_price;
        double S3=initial_stock_price;
        double S4=initial_stock_price;
        double index1 =1.0;
        double index2 =1.0;
        double index3 =1.0;
        double index4 =1.0;
     
        
        
        
        for ( int j=0; j<no_of_discrete_barriers; j++)
        {
            double x = get_uniform();
            double y = get_uniform();
            double a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            double b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
            
            
            S1=S1*exp(R + SD*a);
            S2=S2*exp(R - SD*a);
            S3=S3*exp(R + SD*b);
            S4=S4*exp(R - SD*b);
            
            if (S1 <= barrier_price)
            {
                index1 = 0.0;
            }
            if (S2 <= barrier_price)
            {
                index2 = 0.0;
            }
            if (S3 <= barrier_price)
            {
                index3 = 0.0;
            }
            if (S4 <= barrier_price)
            {
                index4 = 0.0;
            }
            
           
            
            
        }
        
        
        
        discrete_barrier_call_option_price_via_simulation += (max(0.0,S1-strike_price)*index1+max(0.0,S2-strike_price)*index2+max(0.0,S3-strike_price)*index3+max(0.0,S4-strike_price)*index4)/4.0;
        discrete_barrier_put_option_price_via_simulation += (max(0.0,strike_price-S1)*index1+max(0.0,strike_price-S2)*index2+max(0.0,strike_price-S3)*index3+max(0.0,strike_price-S4)*index4)/4.0;
        
        //
        
        
        
    
        
        discrete_barrier_call_option_price_adjust +=((discrete_probability_stock_path_not_hit_barrier(S1))*max(0.0,S1-strike_price)+
                                            (discrete_probability_stock_path_not_hit_barrier(S2))*max(0.0,S2-strike_price)+
                                            (discrete_probability_stock_path_not_hit_barrier(S3))*max(0.0,S3-strike_price)+
                                            (discrete_probability_stock_path_not_hit_barrier(S4))*max(0.0,S4-strike_price))/4.0;
        discrete_barrier_put_option_price_adjust +=((discrete_probability_stock_path_not_hit_barrier(S1))*max(0.0,strike_price-S1)+
                                           (discrete_probability_stock_path_not_hit_barrier(S2))*max(0.0,strike_price-S2)+
                                           (discrete_probability_stock_path_not_hit_barrier(S3))*max(0.0,strike_price-S3)+
                                           (discrete_probability_stock_path_not_hit_barrier(S4))*max(0.0,strike_price-S4))/4.0;
        
    }
    
    
    
    discrete_barrier_call_option_price_via_simulation = exp(-risk_free_rate*expiration_time)*( discrete_barrier_call_option_price_via_simulation /((double) no_of_trials));
    discrete_barrier_put_option_price_via_simulation= exp(-risk_free_rate*expiration_time)*( discrete_barrier_put_option_price_via_simulation /((double) no_of_trials));
    discrete_barrier_call_option_price_adjust =exp(-risk_free_rate*expiration_time)*(discrete_barrier_call_option_price_adjust/((double) no_of_trials));
    discrete_barrier_put_option_price_adjust =exp(-risk_free_rate*expiration_time)*(discrete_barrier_put_option_price_adjust/((double) no_of_trials));
  
    
    cout << "The average Call Price ivia explicit simulation of price paths = " << discrete_barrier_call_option_price_via_simulation << endl;
    cout <<"The average call price with Brownian-Bridge correction on the final price  = " <<discrete_barrier_call_option_price_adjust<<endl;
   
    cout << "--------------------------------" << endl;
    cout << "The average Put Price ivia explicit simulation of price paths =  " << discrete_barrier_put_option_price_via_simulation<<endl;
    cout <<"The average put price with Brownian-Bridge correction on the final price = " <<discrete_barrier_put_option_price_adjust<<endl;
   
    
    cout << "--------------------------------" << endl;
}





