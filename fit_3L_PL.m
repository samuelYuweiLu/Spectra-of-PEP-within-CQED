% clear all;clc

units;
c = 299792458;
h = 6.626e-34;
hb = h / ( 2 * pi );
e0 = 8.85419e-12;
ev = 1.602e-19;
d = 40 * 3.33564e-30;

g0 = 29.01 / 1000;

g1 = 29.01 / 1000;
g2 = 42.95 / 1000;
g3 = 101.7 / 1000;

y1 = 41.94 / 2 / 1000;
y2 = 59.14 / 2 / 1000;
y3 = 67.90 / 2 / 1000;

w1 = 2.787;
w2 = 2.89;
w3 = 2.95;

w = linspace( w3 - 0.4, w3 + 0.3, 2801 );

lc = 13.6;
le = 1;



wr = w*ev/hb;


wa = w1 - 1i * 0.02 / 2;
E0 = 1;
lamb_wa = (g2*ev/hb)^2 / pi * ( w - w2 )*ev/hb ./ ( ( w*ev/hb - w2*ev/hb ).^2 + (y2*ev/hb)^2 ) ...
    + (g3*ev/hb)^2 / pi * ( w - w3 )*ev/hb ./ ( ( w*ev/hb - w3*ev/hb ).^2 + (y3*ev/hb)^2 );
gamma_wa = (g2*ev/hb)^2 / pi * (y2*ev/hb) ./ ( ( w*ev/hb - w2*ev/hb ).^2 + (y2*ev/hb)^2 ) ...
    + (g3*ev/hb)^2 / pi * (y3*ev/hb) ./ ( ( w*ev/hb - w3*ev/hb ).^2 + (y3*ev/hb)^2 );
m_wa = ( lamb_wa - 1i * gamma_wa ) * pi;

lamb1_wa = (g1*ev/hb)^2 / pi * ( w - w1 )*ev/hb ./ ( ( w*ev/hb - w1*ev/hb ).^2 + (y1*ev/hb)^2 );
gamma1_wa = (g1*ev/hb)^2 / pi * (y1*ev/hb) ./ ( ( w*ev/hb - w1*ev/hb ).^2 + (y1*ev/hb)^2 );
m1_wa = ( lamb1_wa - 1i * gamma1_wa ) * pi;



u_w = 1 ./ ( 1i * ( wr + wa*ev/hb ) + (g1*ev/hb)^2 ./ (1i *( wr + ( w1*ev/hb - 1i * y1*ev/hb ) )) + (g2*ev/hb)^2 ./ (1i *( wr + ( w2*ev/hb - 1i * y2*ev/hb ) )) + (g3*ev/hb)^2 ./ (1i *( wr + ( w3*ev/hb - 1i * y3*ev/hb ) )) );
flu_spec = real( 1 ./ ( 1i * wr ) .* ( 1 ./ ( 1i * ( wr - wa*ev/hb ) - 1i * ( m_wa + m1_wa ) ) ) .* u_w );



%  figure4(b)
figure(10);plot( ( w - wa )  * 1000, flu_spec ./ max( flu_spec ), 'LineWidth', 1 );hold on;
figure(10);plot( dev, spec ./ max( spec ), '--', 'LineWidth', 2 );hold on;
xlabel( 'Detuning (meV)' );
ylabel( 'Intensity' );
title( 'PL spectra' );
legend( 'analytical', 'QuTip' );
set( gca, 'Fontname', 'Arial' );
set( gca, 'Fontsize', 16 );
xlim( [ -150, 300 ] );
ylim( [ -0.01, 1.05 ] );
% close(figure(10))

