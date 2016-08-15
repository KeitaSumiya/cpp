//
//  main.cpp
//  integration
//
//  Created by 角谷啓太 on 2016/08/13.
//  Copyright © 2016年 角谷啓太. All rights reserved.
//

#include <iostream>
#include <math.h>

void flow(double t, double x, double y, double u, double v, double m, double *dxdt, double *dydt, double *dudt, double *dvdt) {
    *dxdt = u / m;
    *dydt = v / m;
    *dudt = - x / pow(x * x + y * y, 3/2);
    *dvdt = - y / pow(x * x + y * y, 3/2);
}

void runge(double dt, double m, double *t, double *x, double *y, double *u, double *v) {
    double dxdt1, dxdt2, dxdt3, dxdt4;
    double dydt1, dydt2, dydt3, dydt4;
    double dudt1, dudt2, dudt3, dudt4;
    double dvdt1, dvdt2, dvdt3, dvdt4;
    flow(*t,      *x,            *y,            *u,            *v,            m, &dxdt1, &dydt1, &dudt1, &dvdt1);
    flow(*t+dt/2, *x+dxdt1*dt/2, *y+dydt1*dt/2, *u+dudt1*dt/2, *v+dvdt1*dt/2, m, &dxdt2, &dydt2, &dudt2, &dvdt2);
    flow(*t+dt/2, *x+dxdt2*dt/2, *y+dydt2*dt/2, *u+dudt2*dt/2, *v+dvdt2*dt/2, m, &dxdt3, &dydt3, &dudt3, &dvdt3);
    flow(*t+dt,   *x+dxdt3*dt,   *y+dydt3*dt,   *u+dudt3*dt,   *v+dvdt3*dt,   m, &dxdt4, &dydt4, &dudt4, &dvdt4);
   
    *x = *x + dt * (dxdt1 + 2 * dxdt2 + 2 * dxdt3 + dxdt4)/6;
    *y = *y + dt * (dydt1 + 2 * dydt2 + 2 * dydt3 + dydt4)/6;
    *u = *u + dt * (dudt1 + 2 * dudt2 + 2 * dudt3 + dudt4)/6;
    *v = *v + dt * (dvdt1 + 2 * dvdt2 + 2 * dvdt3 + dvdt4)/6;
    *t = *t + dt;
}

int main(void) {
    double t = 0;
    double t_end = 10;
    double dt = 0.01;
    double x = 1;
    double y = 0;
    double u = 0;
    double En = -0.5;
    double r = sqrt(x*x + y*y);
    double v = sqrt(2*(En - u*u/2 + 1/r));
    //    double v = 1;
    double m = 1;
    int i = 0;
    int ycut = 0;
    double yold = y;

    printf("%lf %lf %lf %lf %lf %3d \n", t, x, y, u, v, ycut);
    while (t < t_end) {
        runge(dt, m, &t, &x, &y, &u, &v);
        if (y*yold<0) {
            ycut++;
        }
        printf("%lf %lf %lf %lf %lf %3d \n", t, x, y, u, v, ycut);
        yold = y;
        i++;
    }
}

