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

with open("ThetaOrionr.txt", "r") as arquivo:

    ind =0  

    for linha in arquivo:


        lista = eval(linha)
       

        arn[ind] = lista[0]
        ind += 1
#print(arn)
print("Desvio",np.std(arn[:ind]))
plt.hist(arn[:ind], bins=50)
plt.axvline(x=3.018202851614777, color="red", label="n_* = 3.02")
plt.axvline(x=3.018202851614777+np.std(arn[:ind]), color="orange", label="Desvio padrão = "+str(round(np.std(arn[:ind]),1)))
plt.axvline(x=3.018202851614777-np.std(arn[:ind]), color="orange")
plt.xlabel("n")
plt.title("Monte Carlo Theta")
plt.legend()
plt.savefig("Theta Ori C.png")
