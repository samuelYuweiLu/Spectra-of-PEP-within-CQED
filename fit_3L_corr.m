% clear all;clc

t = linspace( 0, 15, 1501 );

g0 = 0.02901;
l = 0.04194/2;
y = 2*g0^2/l;
w = 2.787;
d = 0.1;
sz = 1;
a = sqrt(1 - 2*y*l/(l + 1i*d)^2);
At= exp(-1i*w*t/g0 - 1/2*(l + 1i*d)*t/g0).*(cosh((l + 1i*d)/2*a*t/g0) + 1/a*sinh((l + 1i*d)/2*a*t/g0));
A0= exp(-1i*w*0 - 1/2*(l + 1i*d)*0).*(cosh((l + 1i*d)/2*a*0) + 1/a*sinh((l + 1i*d)/2*a*0));
czzt_detuning= 1 - (abs(At).^2 + abs(A0).^2 -  2*conj(At).*At.*A0)*(1 + sz);

plot( t, real(czzt), 'LineWidth', 2.0 );hold on;
plot( t, real(rc), '--', 'LineWidth', 2 );hold on;
plot( t, real(czzt_detuning), 'LineWidth', 2 );hold on;
plot( t, real(rc_detuning), '--', 'LineWidth', 2 );hold on;
xlabel( 'g_{1}t' );
ylabel( 'Re[\sigma_z(\tau)\sigma_z(0)]' );
legend( 'exact, \Delta\omega=0', 'QRT, \Delta\omega=0', 'exact, \Delta\omega=0.1eV', 'QRT, \Delta\omega=0.1eV' );
set( gca, 'Fontname', 'Arial' );
set( gca, 'Fontsize', 16 );
xlim( [ 0, 15 ] );
ylim( [ -1.05, 1.05 ] );