#! /usr/bin/env python

"""
File: graph_printer.py
Copyright (c) 2016 Michael Seaman
License: MIT
Description: Runnable module that creates
graphs of the functions given by ODE solvers
along with their analyitical solutions
"""

from Euler import Euler
from Heun import Heun
from Runge_Kutta2 import Runge_Kutta2
from Runge_Kutta4 import Runge_Kutta4
from math import sin, cos



e = Euler(u0 = 1, v0 = 0, n = 8)
e.generate_graph_with_exact(sin, cos)

e = Euler(u0 = 1, v0 = 0, n = 100)
e.generate_graph_with_exact(sin, cos)

e = Euler(u0 = 1, v0 = 0, n = 1000)
e.generate_graph_with_exact(sin, cos)

e = Heun(u0 = 1, v0 = 0, n = 8)
e.generate_graph_with_exact(sin, cos)

e = Heun(u0 = 1, v0 = 0, n = 25)
e.generate_graph_with_exact(sin, cos)

e = Heun(u0 = 1, v0 = 0, n = 100)
e.generate_graph_with_exact(sin, cos)

e = Runge_Kutta2(u0 = 1, v0 = 0, n = 8)
e.generate_graph_with_exact(sin, cos)

e = Runge_Kutta2(u0 = 1, v0 = 0, n = 25)
e.generate_graph_with_exact(sin, cos)

e = Runge_Kutta2(u0 = 1, v0 = 0, n = 100)
e.generate_graph_with_exact(sin, cos)

e = Runge_Kutta4(u0 = 1, v0 = 0, n = 8)
e.generate_graph_with_exact(sin, cos)

e = Runge_Kutta4(u0 = 1, v0 = 0, n = 25)
e.generate_graph_with_exact(sin, cos)

e = Runge_Kutta4(u0 = 1, v0 = 0, n = 100)
e.generate_graph_with_exact(sin, cos)

