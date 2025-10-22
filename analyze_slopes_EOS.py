#!/usr/bin/env python 
from pylab import *
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline

"""
-----------------------------------------------------------------------------------------------|
 TESTING MIE-GRUNEISEN BY FITTING EOS TABLES                                                   |
                                                                                               |
This code fits our EOS data tables to infer the difference in thermal pressure predicted       |
by Mie Gruneisen model versus the actual thermal pressure.                                     |
                                                                                               |
Felipe Gonzalez                                                          Berkeley, 10/22/2025  |
-----------------------------------------------------------------------------------------------|
Last modified on:                                                                  10/22/2025
"""

