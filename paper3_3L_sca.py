# matplotlib inline
import matplotlib.pyplot as plt
import numpy as np

from qutip import *

# parameter definition
h_bar = 6.626 * 10 ** (-34) / ( 2 * np.pi );
ev = 1.602 * 10 ** (-19);

g0 = 0.02901 * ev / h_bar;

wc1 = 2.787 * ev / h_bar / g0;  # cavity frequency
wc2 = 2.890 * ev / h_bar  / g0;  # cavity frequency
wc3 = 2.950 * ev / h_bar / g0;  # cavity frequency
wa = wc1;  # atom frequency
g1  = 0.02901 * ev / h_bar / g0;  # coupling strength
g2  = 0.04295 * ev / h_bar / g0;  # coupling strength
g3  = 0.10170 * ev / h_bar / g0;  # coupling strength
kappa1 = 0.04194 * ev / h_bar / g0;       # cavity dissipation rate
kappa2 = 0.05914 * ev / h_bar / g0;       # cavity dissipation rate
kappa3 = 0.06790 * ev / h_bar / g0;       # cavity dissipation rate
gamma = 0.02 * ev / h_bar / g0;       # atomic dissipation rate
N = 3;              # number of cavity fock states
omega = 0.01 * kappa1;


# intial state
psi0 = tensor(basis(N,0), basis(N,0), basis(N,0), basis(2,0));    # start with an unexcited atom

# operators
a1  = tensor(destroy(N), qeye(N), qeye(N), qeye(2));
a2  = tensor(qeye(N), destroy(N), qeye(N), qeye(2));
a3  = tensor(qeye(N), qeye(N), destroy(N), qeye(2));
sm = tensor(qeye(N), qeye(N), qeye(N), destroy(2));


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

ue = 1;
uc = 13.6;

# Hamiltonian
def calulcate_avg_photons(wl):
    H = ( wc1 - wl ) * a1.dag() * a1 + ( wc2 - wl ) * a2.dag() * a2 + ( wc3 - wl ) * a3.dag() * a3 + ( wa - wl ) * sm.dag() * sm + g1 * (a1.dag() * sm + a1 * sm.dag()) + g2 * (a2.dag() * sm + a2 * sm.dag()) + g3 * (a3.dag() * sm + a3 * sm.dag()) + omega / 2 * ( uc * ( a1.dag() + a1 ) + ue * ( sm.dag() + sm ) );
   
    rho_ss = steadystate(H, c_ops)

    n_c1 = expect(a1.dag() * a1, rho_ss)
    n_ac1 = expect(a1.dag() * sm, rho_ss)
    n_a = expect(sm.dag() * sm, rho_ss)
    
    return n_c1, n_a, n_ac1
    

wl_vec = np.linspace(2.787-0.2, 2.787+0.2, 201) * ev / h_bar / g0

    
n_c1_vec = []
n_ac1_vec = []
n_a_vec = []


num = 1
for wl in wl_vec:
    n_c1_avg, n_a_avg, n_ac1_avg = calulcate_avg_photons(wl)
    n_c1_vec.append(n_c1_avg)
    n_ac1_vec.append(n_ac1_avg)
    n_a_vec.append(n_a_avg)
#    print(num)
    num = num + 1


# final plot 1
fig, axes = plt.subplots(1, 1, figsize=(10,6))


dev = 1000*(wl_vec-wc1) / (ev/h_bar/g0);

axes.plot(dev, 300 * ue**2 * np.array( n_a_vec ), label="Atomic")
axes.plot(dev, uc**2 * np.array( n_c1_vec ), label="dipolar mode")
axes.plot(dev, 10 * 2 * ue * uc * np.array( np.real( n_ac1_vec ) ), label="atomic & dipolar mode")
axes.legend(loc=0)
axes.set_xlabel('Detuning (meV)')
axes.set_ylabel('Insentity')
axes.set_title('Steady-state scattering spectra')    
axes.set_xlim(-150, 150);  
  





