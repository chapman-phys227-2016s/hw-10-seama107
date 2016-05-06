#! /usr/bin/env python

"""
File: Runge_Kutta4.py
Copyright (c) 2016 Michael Seaman
License: MIT
Description: Implements the fourth order
Runge Kutta for solving a set of 2 coupled ODEs
"""

from ODE_solver import ODE_solver
import numpy as np
import matplotlib.pyplot as plt

class Runge_Kutta4(ODE_solver):
    """
    Implements the solve_n method using Heun's method
    """
    def solve_n(self):
        output_u = np.zeros(self.n+1) + self.u0
        output_v = np.zeros(self.n+1) + self.v0
        t_list = np.linspace(self.a, self.b, self.n+1)
        delta_t = (self.b - self.a)/ float(self.n)
        for i in xrange(len(t_list) - 1):
            k1_u = delta_t * self.dudt(output_u[i], output_v[i], t_list[i])
            k1_v = delta_t * self.dvdt(output_u[i], output_v[i], t_list[i])
            k2_u = delta_t * self.dudt(output_u[i] + .5 * k1_u, output_v[i] + .5 * k1_v, t_list[i] + delta_t * .5)
            k2_v = delta_t * self.dvdt(output_u[i] + .5 * k1_u, output_v[i] + .5 * k1_v, t_list[i] + delta_t * .5)
            k3_u = delta_t * self.dudt(output_u[i] + .5 * k2_u, output_v[i] + .5 * k2_v, t_list[i] + delta_t * .5)
            k3_v = delta_t * self.dvdt(output_u[i] + .5 * k2_u, output_v[i] + .5 * k2_v, t_list[i] + delta_t * .5)
            k4_u = delta_t * self.dudt(output_u[i] + k1_u, output_v[i] + k1_v, t_list[i] + delta_t)
            k4_v = delta_t * self.dvdt(output_u[i] + k1_u, output_v[i] + k1_v, t_list[i] + delta_t)
            output_u[i+1] = output_u[i] + (k1_u + 2*k2_u + 2*k3_u + k4_u)/6.0
            output_v[i+1] = output_v[i] + (k1_v + 2*k2_v + 2*k3_v + k4_v)/6.0
        return output_u, output_v