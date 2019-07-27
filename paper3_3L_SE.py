# -*- coding: utf-8 -*-
"""
Created on Sat Jul  7 10:49:37 2018

@author: luyuwei
"""

# matplotlib inline
import matplotlib.pyplot as plt
import numpy as np

from qutip import *

# parameter definition
h_bar = 6.626 * 10 ** (-34) * 2 * np.pi;
ev = 1.602 * 10 ** (-19);

g0 = 0.02901 * ev / h_bar;

wc1 = 2.787 * ev / h_bar / g0;  # cavity frequency
wc2 = 2.890 * ev / h_bar  / g0;  # cavity frequency
wc3 = 2.950 * ev / h_bar / g0;  # cavity frequency
wa = 2.7870 * ev / h_bar / g0;  # atom frequency
#g1=0;g2=0;
g1  = 0.02901 * ev / h_bar / g0;  # coupling strength
g2  = 0.04295 * ev / h_bar / g0;  # coupling strength
g3  = 0.10170 * ev / h_bar / g0;  # coupling strength
kappa1 = 0.04194 * ev / h_bar / g0;       # cavity dissipation rate
kappa2 = 0.05914 * ev / h_bar / g0;       # cavity dissipation rate
kappa3 = 0.06790 * ev / h_bar / g0;       # cavity dissipation rate
gamma = 0.0 * ev / h_bar / g0;       # atomic dissipation rate
N = 3;              # number of cavity fock states
use_rwa = True;


tlist = np.linspace(0, 4, 121);

# intial state
psi0 = tensor(basis(N,0), basis(N,0), basis(N,0), basis(2,1));    # start with an excited atom

# operators
a1  = tensor(destroy(N), qeye(N), qeye(N), qeye(2));
a2  = tensor(qeye(N), destroy(N), qeye(N), qeye(2));
a3  = tensor(qeye(N), qeye(N), destroy(N), qeye(2));
sm = tensor(qeye(N), qeye(N), qeye(N), destroy(2));

# Hamiltonian
if use_rwa:
    H = wc1 * a1.dag() * a1 + wc2 * a2.dag() * a2 + wc3 * a3.dag() * a3 + wa * sm.dag() * sm + g1 * (a1.dag() * sm + a1 * sm.dag()) + g2 * (a2.dag() * sm + a2 * sm.dag()) + g3 * (a3.dag() * sm + a3 * sm.dag());
else:
    H = wc3 * a3.dag() * a3 + wa * sm.dag() * sm + g3 * (a3.dag() + a3) * (sm + sm.dag());
    
    
c_ops = [];

# cavity relaxation
rate1 = kappa1;
rate2 = kappa2;
rate3 = kappa3;
rate4 = gamma;
c_ops.append(np.sqrt(rate1) * a1);
c_ops.append(np.sqrt(rate2) * a2);
c_ops.append(np.sqrt(rate3) * a3);
c_ops.append(np.sqrt(rate4) * sm);

output = mesolve(H, psi0, tlist, c_ops, [a1.dag() * a1, a2.dag() * a2, a3.dag() * a3, sm.dag() * sm]);

# final plot 1
n_c1 = output.expect[0]
n_c2 = output.expect[1]
n_c3 = output.expect[2]
n_c = output.expect[0] + output.expect[1] + output.expect[2]
n_a = output.expect[3]

fig, axes = plt.subplots(1, 1, figsize=(10,6))

axes.plot(tlist, n_c, label="Cavity")
axes.plot(tlist, n_a, label="Atom excited state")
axes.plot(tlist, n_c1, '--', label="Cavity 1")
axes.plot(tlist, n_c2, '--', label="Cavity 2")
axes.plot(tlist, n_c3, '--', label="Cavity 3")
axes.legend(loc=0)
axes.set_xlabel('Time')
axes.set_ylabel('Occupation probability')
axes.set_title('Vacuum Rabi oscillations')    
axes.set_xlim(0, 4);  
axes.set_ylim(0, 1.005);    
    
