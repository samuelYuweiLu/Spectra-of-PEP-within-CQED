% clear all;clc

%% this program is for paper5
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

%% spectral density
s0w = (g1*ev/hb)^2 / pi * (y1*ev/hb) ./ ( ( w*ev/hb - w1*ev/hb ).^2 + (y1*ev/hb)^2 ) ...
    + (g2*ev/hb)^2 / pi * (y2*ev/hb) ./ ( ( w*ev/hb - w2*ev/hb ).^2 + (y2*ev/hb)^2 ) ...
    + (g3*ev/hb)^2 / pi * (y3*ev/hb) ./ ( ( w*ev/hb - w3*ev/hb ).^2 + (y3*ev/hb)^2 );
slw = (g1*ev/hb)^2 / pi * ( w - w1 )*ev/hb ./ ( ( w*ev/hb - w1*ev/hb ).^2 + (y1*ev/hb)^2 ) ...
    + (g2*ev/hb)^2 / pi * ( w - w2 )*ev/hb ./ ( ( w*ev/hb - w2*ev/hb ).^2 + (y2*ev/hb)^2 ) ...
    + (g3*ev/hb)^2 / pi * ( w - w3 )*ev/hb ./ ( ( w*ev/hb - w3*ev/hb ).^2 + (y3*ev/hb)^2 );

enei = eV2nm ./ w;
k = 2 * pi ./ ( enei * 1e-9 );
s0_w = s0w / ( d^2 / ( pi * e0 * hb ) );
sl_w = slw / ( d^2 / ( pi * e0 * hb ) );
G0 = k.^3 / ( 6 * pi );

wr = w*ev/hb;


%% aa spectrum
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


    
Iaa = abs( E0 ./ ( 1i * ( wr - ( w1*ev/hb - 1i * y1*ev/hb ) ) ) .* ( ( 1i * ( wr - wa*ev/hb ) - 1i * m_wa ) * lc + 1i * (g1*ev/hb) * le ) ./ ( 1i * ( wr - wa*ev/hb ) - 1i * ( m_wa + m1_wa ) ) ).^2;  
Iss = abs( E0 .* ( -1i*m1_wa/(g1*ev/hb) * lc - 1i * le ) ./ ( 1i * ( wr - wa*ev/hb ) - 1i * ( m_wa + m1_wa ) ) ).^2;
Ias = 2 * real( 1i * conj( E0 ./ ( 1i * ( wr - ( w1*ev/hb - 1i * y1*ev/hb ) ) ) .* ( ( 1i * ( wr - wa*ev/hb ) - 1i * m_wa ) * lc + 1i * (g1*ev/hb) * le ) ./ ( 1i * ( wr - wa*ev/hb ) - 1i * ( m_wa + m1_wa ) ) ) .* ...
    ( E0 .* ( -1i*m1_wa/(g1*ev/hb) * lc - 1i * le ) ./ ( 1i * ( wr - wa*ev/hb ) - 1i * ( m_wa + m1_wa ) ) ) );




%  figure4(a)
figure(11);plot( ( w - wa ) * 1000, le^2 * Iss * 300  / max( lc^2 * Iaa ), 'LineWidth', 1 );hold on;
figure(11);plot( ( w - wa ) * 1000, lc^2 * Iaa  / max( lc^2 * Iaa ), 'LineWidth', 1 );hold on;
figure(11);plot( ( w - wa ) * 1000, le*lc * Ias * 10  / max( lc^2 * Iaa ), 'LineWidth', 1 );hold on;
% figure(11);plot( ( w - wa ) * 1000, le^2 * Iss + lc^2 * Iaa + le*lc * Ias, 'k--', 'LineWidth', 2 );hold on;


figure(11);plot( nev, le^2 * nav * 300 / max( lc^2 * ncv ), '--', 'LineWidth', 2 );hold on;
figure(11);plot( nev, lc^2 * ncv / max( lc^2 * ncv ), '--', 'LineWidth', 2 );hold on;
figure(11);plot( nev, le*lc * nacv * 10 / max( lc^2 * ncv ), '--', 'LineWidth', 2 );hold on;
ylabel( 'Intensity' );
xlabel( 'Detuning (meV)' );
title( 'Steady-state scattering spectra' );
legend( 'Iss, anayltical', 'Iaa, anayltical', 'Ias, anayltical', 'Iss, QuTip', 'Iaa, QuTip', 'Ias, QuTip' )
set( gca, 'Fontname', 'Arial' );
set( gca, 'Fontsize', 16 );
xlim( [ -150, 150 ] );
% close(figure(11))


