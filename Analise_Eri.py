# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 08:17:58 2023

@author: diogo
"""

import numpy as np
import matplotlib.pyplot as plt

plt.style.use('aa.mplstyle')
plt.rcParams['font.size'] = 16
arn = np.zeros(6000)

with open("EpsilonEridani.txt", "r") as arquivo:

    ind =0  

    for linha in arquivo:


        lista = eval(linha)
       

        arn[ind] = lista[0]
        ind += 1
#print(arn)
print("Desvio",np.std(arn[:ind]))
plt.hist(arn[:ind], bins=50)
plt.axvline(x=2.634235837475673, color="red", label="n_* = 2.63")
plt.axvline(x=2.634235837475673+np.std(arn[:ind]), color="orange", label="Desvio padrão = "+str(round(np.std(arn[:ind]),2)))
plt.axvline(x=2.634235837475673-np.std(arn[:ind]), color="orange")
plt.xlabel("n")
plt.title("Monte Carlo Epsilon")
plt.legend()
plt.savefig("Epsilon.png")
