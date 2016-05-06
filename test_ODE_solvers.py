#! /usr/bin/env python

"""
File: test_ODE_solvers.py
Copyright (c) 2016 Michael Seaman
License: MIT
Description: The base class for all ODE_solvers in homework 10
that will be implemented by Euler, Heuns, RK-2, and RK-4
"""

import numpy as np
from Euler import Euler
from Heun import Heun
from Runge_Kutta2 import Runge_Kutta2
from Runge_Kutta4 import Runge_Kutta4

def linear_dudt(u, v, t):
    return 2

def linear_dvdt(u, v, t):
    return -.5

def general_ODE_linear_tst(ODEsolver_class):
    """
    The general linear test case method - Takes a class name as parameter
    and tests that it approximates a line perfectly (as all of these methods
    should). They run from t = 0 to 10 with steps of .1. dvdt = -.5 and dudt= 2
    """
    a = 0
    b = 10
    n = 100
    u0 = 0
    v0 = 0
    t_list = np.linspace(a, b, n+1)
    delta_t = (b - a)/ float(n)
    analyitic_result_u = t_list * linear_dudt(u0, v0, a)
    analyitic_result_v = t_list * linear_dvdt(u0, v0, a)
    ODEsolver = ODEsolver_class(dudt = linear_dudt, dvdt = linear_dvdt, u0 = u0, v0 = v0, a = a, b = b,  n = n)
    appr_results = ODEsolver.solve_n()
    print appr_results
    approx_result_u = appr_results[0]
    approx_result_v = appr_results[1]
    print "exact u: " + str(analyitic_result_u[:10])
    print "appr u: " + str(approx_result_u[:10])
    print "exact v: " + str(analyitic_result_v[:10])
    print "appr v: " + str(approx_result_v[:10])
    print np.allclose(analyitic_result_u, approx_result_u)
    print np.allclose(analyitic_result_v, approx_result_v)
    return np.allclose(analyitic_result_u, approx_result_u) and np.allclose(analyitic_result_v, approx_result_v)

def test_Euler_linear():
    assert general_ODE_linear_tst(Euler)

def test_Heun_linear():
    assert general_ODE_linear_tst(Heun)

def test_Runge_Kutta2_linear():
    assert general_ODE_linear_tst(Runge_Kutta2)

def test_Runge_Kutta4_linear():
    assert general_ODE_linear_tst(Runge_Kutta4)
