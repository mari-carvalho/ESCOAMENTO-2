# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 16:12:23 2023

@author: Usuário
"""

rs = ((( P / 112.727 ) + 12.34 )*(dg ** 0.8439)*(10 * a_rs)) ** 1.73184   # P em psia e T em F
bo = 1.0113 + (7.2046 * 10 ^ -5) * ( ((((( P / 112.727 ) + 12.34 )*(dg ^ 0.8439)*(10 * a_rs)) ^ 1.73184) ^ 0.3738) * (dg ^ 0.2914/ do ^ 0.6265) + 0.24626 * (T ^ 0.5371)) ^ 3.0936)   # T em F

(1.0113 + (7.2046 * 10 **-5 ) * ( ( ( ( ( (P / 112.727) + 12.34) ** (dg ** 0.8439) * (10 * a_rs) ) * 1.73184) * 0.3738) * (dg ** 0.2914 / 0.6265) + 0.24626 ** (T ** 0.5371) )