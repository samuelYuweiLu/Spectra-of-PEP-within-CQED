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

w = linspace( w3 - 0.4, w3 + 0.3, 1001 );


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


% Lamb shift
iG_f = s0_w; 
rG_f = zeros( length( iG_f ) );
for n = 1:length( iG_f )
    rG_f( n ) = ( trapz( wr( 1, 1:( n - 1 ) ), iG_f( 1, 1:( n - 1 ) ) ./ ( - wr( 1, 1:( n - 1 ) ) + wr( 1, n ) ), 2 ) + ...
        trapz( wr( 1, ( n + 1 ):end ), iG_f( 1, ( n + 1 ):end ) ./ ( - wr( 1, ( n + 1 ):end ) + wr( 1, n ) ), 2 ) ) / pi;
end
sgammax = iG_f * 2 * d^2 / ( hb * e0 );
slambshiftx = rG_f( :, 1 )' * d^2 / ( hb * e0 );

%% dipole spectrum
wl = w1;
sCw = 1 / pi * sgammax / 2 ./ ( ( ( wr - wl * ev / hb ) - slambshiftx ).^2 + ( 1 / 2 * sgammax ).^2 );
figure(3);plot( ( w - wl ) / g3, abs(sCw) / max( abs(sCw) ), 'LineWidth', 2 );hold on;
% figure(3);plot( ( w - wl ) / g3, sgammax / max( sgammax ), 'k--', 'LineWidth', 1.5 );hold on;
xlim( [ -2, 4 ] );
ylabel( 'Intensity');
xlabel( '( \omega - \omega_{\it{l}} ) / {\it{g_3}}' );
set( gca, 'Fontname', 'Arial' );
set( gca, 'Fontsize', 16 );
title( 'Dipole spectrum' );
% close(figure(3))

%% time domain
wr = wr';
dws = wr(2:end) - wr(1:end-1);
dws = [ dws(1) ;dws ];
dt = 0.002;
tmax = 2;
ta_0 = 2.56 * 10^-9 / 10^4;
t = linspace( 0, tmax, tmax / dt );
t = t * ta_0;
tt = repmat( t', [ 1, size( wr, 1 ) ] );
wss = repmat( wr', [ tmax / dt, 1 ] );
dwss = repmat( dws', [ tmax / dt, 1 ] );
sgammaxx = repmat( sgammax, [ tmax / dt, 1 ] );
slambshiftxx = repmat( slambshiftx, [ tmax / dt, 1 ] );
Ctx = 1 / pi * sgammaxx / 2 ./ ( ( ( wss - wl * ev / hb ) - slambshiftxx ).^2 + ( 1 / 2 * sgammaxx ).^2 ) .* exp( -1i * wss .* tt ) .* dwss;
Cet = abs( sum( Ctx, 2 ) ).^2;

figure(5);plot( t*g0*ev/hb, Cet, 'LineWidth', 1 );hold on;
figure(5);plot( tlist, ta, '--', 'LineWidth', 2 );hold on;
xlim( [ 0, 4 ] );
ylim( [ 0, 1 ] );
ylabel( 'Population');
xlabel( 't / (eV/\hbar)' );
legend( 'Green''s tensor', 'QuTip' );
set( gca, 'Fontname', 'Arial' );
set( gca, 'Fontsize', 16 );
% close(figure(5));

