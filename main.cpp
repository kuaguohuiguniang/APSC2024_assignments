#include <iostream>
#include <cmath>
#include <functional>

typedef std::function<double(std::vector<double> &)> func;
typedef std::function<std::vector<double>(std::vector<double> &)> dfunc;

// structure containing all coefficients
struct problem
{
    int max_it = 10000;
    double tol_r = 1e-9;    // here I changed the tolrance from 1e-6 to 1e-9 
    double tol_s = 1e-9;    // to obtain better convergence result
    double alpha_0 = 0.1;
    double mu = 0.2;
    std::vector<double> x0{0., 0.};
};

// the function f to de solved
double f(std::vector<double> & x){
    double result = x[0]*x[1]+4*pow(x[0],4)+pow(x[1],2)+3*x[0];
    return result;
}

// derivative of f, df
std::vector<double> df(std::vector<double> & x){
    std::vector<double> result;
    result.push_back(x[1]+16*pow(x[0],3)+3);
    result.push_back(x[0]+2*x[1]);
    return result;
}

// compute vector norm
double norm(std::vector<double> x){
    double result = 0.;
    for (unsigned i=0; i < x.size(); i++){
        result += pow(x[i],2);
    }
    result = pow(result, 0.5);
    return result;
}

// enable vector minus vector
std::vector<double> operator-(const std::vector<double>& lhs, const std::vector<double>& rhs){
    std::vector<double> result;
    for (unsigned i=0; i < lhs.size(); i++){
        result.push_back(lhs[i]-rhs[i]);
    }
    return result;
}

// enable scalar times vector
std::vector<double> operator*(const double& lhs, const std::vector<double>& rhs){
    std::vector<double> result;
    for (unsigned i=0; i < rhs.size(); i++){
        result.push_back(lhs*rhs[i]);
    }
    return result;
}

// the solver, used inverse decay method
std::vector<double> solver(func f, dfunc df, struct problem p){
    double alpha = p.alpha_0;
    std::vector<double> x1 = p.x0;
    std::vector<double> x2 = x1-alpha*df(x1);
    double res_r = fabs(f(x1)-f(x2));
    double res_s = norm(x1-x2);
    int it = 1;
    while (res_r > p.tol_r && res_s > p.tol_s && it < p.max_it){
        //alpha = p.alpha_0*exp(-p.mu*it);
        alpha = p.alpha_0/(1+p.mu*it);
        x1 = x2;
        x2 = x1-alpha*df(x1);
        res_r = fabs(f(x1)-f(x2));
        res_s = norm(x1-x2);
        it++;
    }
    return x2;
}

int main(){
    struct problem p;   // Using default coefficients
    std::vector<double> x = solver(f, df, p);
    std::cout<<"x="<<x[0]<<' '<<x[1]<<std::endl;
    std::cout<<"f(x)="<<f(x)<<std::endl;
    return 0;
}

