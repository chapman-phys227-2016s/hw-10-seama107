#! /usr/bin/env python

"""
File: Heun.py
Copyright (c) 2016 Michael Seaman
License: MIT
Description: Implements Heun's method for solving 
a set of 2 coupled ODEs
"""

from ODE_solver import ODE_solver
import numpy as np
import matplotlib.pyplot as plt

class Heun(ODE_solver):
    """
    Implements the solve_n method using Heun's method
    """
    def solve_n(self):
        output_u = np.zeros(self.n+1) + self.u0
        output_v = np.zeros(self.n+1) + self.v0
        t_list = np.linspace(self.a, self.b, self.n+1)
        delta_t = (self.b - self.a)/ float(self.n)
        for i in xrange(len(t_list) - 1):
            geussed_u_next = output_u[i] + delta_t * self.dudt(output_u[i], output_v[i], t_list[i])
            geussed_v_next = output_v[i] + delta_t * self.dvdt(output_u[i], output_v[i], t_list[i])
            output_u[i+1] = output_u[i] + .5 * delta_t * (self.dudt(output_u[i], output_v[i], t_list[i]) + self.dudt(geussed_u_next, geussed_v_next, t_list[i+1]) )
            output_v[i+1] = output_v[i] + .5 * delta_t * (self.dvdt(output_u[i], output_v[i], t_list[i]) + self.dvdt(geussed_u_next, geussed_v_next, t_list[i+1]) )
        return output_u, output_v

