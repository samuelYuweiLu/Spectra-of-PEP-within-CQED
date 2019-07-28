# matplotlib inline
import matplotlib.pyplot as plt
import numpy as np

from qutip import *

# parameter definition
h_bar = 6.626 * 10 ** (-34) / ( 2 * np.pi );
ev = 1.602 * 10 ** (-19);

g0 = 0.02901 * ev / h_bar;

wc1 = 2.787 * ev / h_bar / g0;  # cavity frequency
wa = 2.887 * ev / h_bar / g0;  # atom frequency
g1  = 0.02901 * ev / h_bar / g0;  # coupling strength
kappa1 = 0.04194 * ev / h_bar / g0;       # cavity dissipation rate
N = 2;              # number of cavity fock states



# intial state
psi0 = tensor(basis(N,0), basis(2,1));    # start with an excited atom

# operators
a1  = tensor(destroy(N), qeye(2));
sm = tensor(qeye(N), destroy(2));


# Hamiltonian
H = ( wc1 - 0 ) * a1.dag() * a1 + ( wa - 0 ) * sm.dag() * sm + g1 * (a1.dag() * sm + a1 * sm.dag());


c_ops = [];

# cavity relaxation
rate1 = kappa1;
c_ops.append(np.sqrt(rate1) * a1);
   
tlist = [ 0, 10 ];
output = mesolve(H, psi0, tlist, c_ops, []);

taulist = np.linspace(0, 15, 1501);
corr_vec = correlation_2op_1t(H, output.states[0], taulist, c_ops, (sm.dag()*sm - sm*sm.dag()), (sm.dag()*sm - sm*sm.dag()),solver='me')


# exact two-time correlation for the same atomic operator
# see e.g. Eq.(59) in PRA 97, 052101 (2018)

g0 = 0.02901;
l = 0.04194/2;
#w = 2.787;
y = 2 * g0 **2 / l;
sz = 1;
a = complex(0, np.sqrt( 2 * y / l - 1 ));
At = [ np.exp(-1/2*l*t/g0)*(np.cosh(l/2*a*t/g0)+1/a*np.sinh(l/2*a*t/g0)) for t in taulist ]
At = np.array(At);
A0 = np.exp(-1/2*l*0)*np.cosh(l/2*a*0)+1/a*np.sinh(l/2*a*0)

Czzt = 1 - (np.abs(At)**2+np.abs(A0)**2-2*np.conjugate(At)*At*A0)*(1+sz);


fig, ax = plt.subplots(1, 1)
ax.plot(taulist, np.real(Czzt), '-', label='Exact',linewidth=1.2)
ax.plot(taulist, np.real(corr_vec), '--', label='Qutip',linewidth=2.0)
ax.set_xlim(0, 15)
ax.legend(loc=0, fontsize=14)
ax.set_xlabel(r'$g_{1}t$', fontsize=16)
ax.set_ylabel(r'$Re[\sigma_z(t+\tau) \sigma_z(\tau)]$', fontsize=16)





