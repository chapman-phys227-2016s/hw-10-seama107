#! /usr/bin/env python

"""
File: ODE_solver.py
Copyright (c) 2016 Michael Seaman
License: MIT
Description: The base class for all ODE_solvers in homework 10
that will be implemented by Euler, Heuns, RK-2, and RK-4
"""

import numpy as np
from np import pi
import matplotlib.pyplot as plt

class ODE_solver():
    """
    Super class that other ODE solving classes implement. Is limited 
    to two coupled equations: u and v.
    
    """
    def __init__(self, dudt = -1, dvdt = -1, u0 = 0, v0 = 0, a = 0, b = 10 * pi,  n = 10):
        self.u0 = u0
        self.v0 = v0
        self.a = a
        self.b = b
        self.n = n
        if dudt == -1:
            self.dudt = default_dudt
        else:
            self.dudt = dudt
        if dvdt == -1:
            self.dvdt = default_dvdt
        else:
            self.dvdt = dudt

    def default_dudt(u, v, t):
        return v
    
    def default_dvdt(u, v, t):
        return -1 * u
    

    def solve_n(self):
        """
        Abstract - Implemented by it's children
        """
        raise NotImplementedError("no rule in class %s" % self.__class__.__name__)

    def generate_graph_with_exact(self, exact_u, exact_v):
        """
        Generates images ('.png's) of the plots of the ODE_solver's
        solutions as well as the exact solutions.
        """
        
        approx_values = solve_n()
        u_list_approx = approx_values[0]
        v_list_approx = approx_values[1]
        
        u_e = np.vectorize(exact_u)
        v_e = np.vectorize(exact_v)
        t_list = np.linspace(self.a, self.b, self.n)
        u_list_exact = u_e(t_list)
        v_list_exact = v_e(t_list)
        
        fig, ax = plt.subplots(nrows = 1, ncols = 1)
        ax.plot(u_list_approx, v_list_approx, 'r')
        ax.plot(u_list_exact, v_list_exact, 'g--')
        ax.xlabel("u(t)")
        ax.ylabel("v(t)")
        ax.title("Approx ODE solution with Exact solution for n = " + str(self.n))
        ax.grid(True)
        fig.savefig('coupledODE_%3d.png' % (self.n))
        plt.close(fig)
