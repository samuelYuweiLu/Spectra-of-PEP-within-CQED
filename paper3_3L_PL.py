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
wl = wc1;  # laser frequency
g1  = 0.02901 * ev / h_bar / g0;  # coupling strength
g2  = 0.04295 * ev / h_bar / g0;  # coupling strength
g3  = 0.10170 * ev / h_bar / g0;  # coupling strength
kappa1 = 0.04194 * ev / h_bar / g0;       # cavity dissipation rate
kappa2 = 0.05914 * ev / h_bar / g0;       # cavity dissipation rate
kappa3 = 0.06790 * ev / h_bar / g0;       # cavity dissipation rate
gamma = 0.02 * ev / h_bar / g0;       # atomic dissipation rate
N = 3;              # number of cavity fock states
omega = 0.1* kappa1;



# intial state
psi0 = tensor(basis(N,0), basis(N,0), basis(N,0), basis(2,0));    # start with an unexcited atom

# operators
a1  = tensor(destroy(N), qeye(N), qeye(N), qeye(2));
a2  = tensor(qeye(N), destroy(N), qeye(N), qeye(2));
a3  = tensor(qeye(N), qeye(N), destroy(N), qeye(2));
sm = tensor(qeye(N), qeye(N), qeye(N), destroy(2));


# Hamiltonian
H = ( wc1 - wl ) * a1.dag() * a1 + ( wc2 - wl ) * a2.dag() * a2 + ( wc3 - wl ) * a3.dag() * a3 + ( wa - wl ) * sm.dag() * sm + g1 * (a1.dag() * sm + a1 * sm.dag()) + g2 * (a2.dag() * sm + a2 * sm.dag()) + g3 * (a3.dag() * sm + a3 * sm.dag()) + omega / 2 * ( ( sm.dag() + sm ) );


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


taulist = np.linspace(0, 50, 5001);
h_exp0 = [ complex(0,wl*t) for t in taulist ];
h_exp1 = np.exp(np.array(h_exp0));
corr_vec = correlation_2op_1t(H, None, taulist, c_ops, sm, sm.dag(),solver='me')
wlist, spec = spectrum_correlation_fft(taulist, h_exp1 * (corr_vec-corr_vec[-1]))


fig, ax = plt.subplots(1, 1)
ax.plot(taulist, np.real(corr_vec-corr_vec[-1]))
ax.plot(taulist, np.imag(corr_vec-corr_vec[-1]))
ax.set_xlim(0, 15)



fig, ax = plt.subplots(1, 1)
ax.plot(1000*(wl-wlist) / (ev/h_bar/g0), spec)
ax.set_xlabel('Detuning (meV)')
ax.set_ylabel('Power spectrum')
ax.set_title('Vacuum Rabi splitting')
ax.set_xlim(-150, 300)

