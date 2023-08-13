#include <iostream>
#include <cmath>
#include "Utils.hpp"
using namespace std;

void displayV(double *v, int M){

    cout << "............" << endl;
    for (int i=0;i<M;i++){
        cout << v[i] << endl;
    }
    cout << "............" << endl;

}

double normalCDF(double x){
    return std::erfc(-x / std::sqrt(2)) / 2;
}

double utility(double c){
    //if (c <= 0){
    //    return -10000000;
    //} else {
    return log(c);
    //}
}

double h(double y, double y_bar){
    if (y <= y_bar){
        return y;
    } else {
        return y_bar;
    }
}
