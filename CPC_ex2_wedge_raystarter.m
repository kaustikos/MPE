clear variables;
%close all;

clc;

set(0, 'DefaultAxesFontSize', 20, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 20, 'DefaultTextFontName', 'Arial');

% ~~~~~~~~~~ z=0                                -> sea surface
% -   -   -    -                }   water column
% (o) -   -    -       z=z_s=100 m                  -> source
%  -    -   -    -              }   c = 1500, \rho = 1
%-------------------- z=h=200 m                     -> sea bottom
% oooooooooooooo       }    c=1700, \rho = 1.5
% oooooooooooooo



% problem parameters

f = 25;                 % source frequency
omeg = 2*pi*f;     % angular frequency
c = 1500;             % sound speed in the water
zs = 100;             % source depth

% domain and grid parameters

xmax = 25000;
dx = 50;
x = 0:dx:xmax;
nx = length(x);

ymax = 4000;
dy = 2;
y = -ymax:dy:ymax;
ny = length(y);

% bottm relief parameters

h0 = 200;

%!!!!! FLAT BOTTOM TEST CASE, example 1
talph = 200/4000;
hy = h0 + y*talph;


PEpars.x = x;
PEpars.y = y;

% Pade approximation orders
PEpars.nP =9;


% PML parameters

PEpars.sigma = -5;
PEpars.lThickness = 1000;

% the array containing the solution 
% i.e., complex acoustical pressure

P(1:ny,1:nx) = 0;

% computed by the CAMBALA solver of the spectral problem
% see https://github.com/Nauchnik/Acoustics-at-home

% wavenumbers for various values of water depth
kh = dlmread('ASA_wedge_kj/kj_wedge_att.txt');

% mode eigenfunctions at z = z_s = 100 m for various values of water depth
phizsh = dlmread('ASA_wedge_kj/phizs_wedge.txt');

% mode eigenfunctions at z = z_r = 30 m for various values of water depth
phizrh = dlmread('ASA_wedge_kj/phizr_wedge.txt');

% number of modes taken into account

mmax = 8;

% some auxiliary variables

ih1 = find(hy>=kh(end,1),1,'first');
ih2 = find(hy<=kh(1,1),1,'last');

iy0 = find(y>=0,1,'first');


for ii = 1:mmax
    
    disp(['Solving MPEs: ' int2str(ii)  ' of ' int2str(mmax)  ' modes'] );
    
    
    %PEpars.nP = nP;
    
    ky = interp1(kh(:,1),kh(:,ii+1),hy);
    ky(1:ih1) = ky(ih1);
    ky(ih2:end) = ky(ih2);
    
    krs = real(ky(iy0));
    
    PEpars.k2 = (ky.^2).';
    PEpars.k02 = (krs).^2;
    
    wmode(1:ny,1) = interp1(phizsh(:,1),phizsh(:,ii+1),hy);
    
    wmode_zs = wmode(iy0);
    
    % ray starter
    
    x0 = 50;
    ix0 = find(x>=x0,1,'first');
    
    th_y = atan(y/x0);
    thmax = pi*80/180;
    ith1 = find(th_y<=-thmax,1,'last');
    ith2 = find(th_y>=thmax,1,'first');
    A_th = ones(size(y));
    th_width = 5*pi/180;
    
    A_th(1:ith1) = exp(  (-(th_y(1:ith1) + thmax ).^2)/(th_width^2) );
    A_th(ith2:end) = exp(  (-(th_y(ith2:end) - thmax ).^2)/(th_width^2) );
    
    
    PEpars.ic(1:ny) = wmode_zs*( exp( 1i*krs*sqrt( x0^2 + y.^2 ) )  ).*A_th./sqrt(  8*pi*krs*sqrt( x0^2 + y.^2 )  ) ;
    
    
    % --- SSP solution for WAMPE ---
    
    Aj = WAPEsolveSSP_PML(PEpars);
    
    % ---
    
    % ray starter for x<x0
    
    Aj(1:ny,ix0:end) = Aj(1:ny,1:end-ix0+1);
    
    for jj = 1:ix0-1
        Aj(1:ny,jj) = wmode_zs*( exp( 1i*krs*sqrt( (x(jj))^2 + y.^2 ) )  ).*A_th./sqrt(  8*pi*krs*sqrt( (x(jj))^2 + y.^2 )  ) ;
    end;
    
    % end ray starter
    
    
    wmode(1:ny,1) = interp1(phizrh(:,1),phizrh(:,ii+1),hy);
    wmode(1:ih1) = wmode(ih1);
    wmode(ih2:end) = wmode(ih2);
    
    % computing sound pressure via the modal expansion
    
    P(1:ny,1:nx) = P(1:ny,1:nx) + Aj(1:ny,1:nx).*repmat(wmode,1,nx);

    
    % showing amplitudes of all modes
    
    TLxy = 20*log10(abs(Aj));
    
    figure;
    imagesc(x/1000,y/1000,TLxy );
    caxis([-80 -40]);
    xlabel('x, km');
    ylabel('y, km');

    
  
end;

TLxy = 20*log10(abs(4*pi*P));


figure;
imagesc(x/1000,y/1000,TLxy );
caxis([-80 -50]);
xlabel('x, km');
ylabel('y, km');



iy0 = find(y>=0,1,'first');
TLtrack = TLxy(iy0,:);
figure;
hold all;
plot(x/1000,TLtrack);

dlmwrite(['WAMPE_pwedge_f=' int2str(f) '_Hz_SSP_raystarter.txt'],[x.'/1000 TLtrack.'],'delimiter','\t','precision',4);

xx = dlmread('ASA_wedge_kj/cross_wedge_ASA.TL');
plot(xx(:,1),xx(:,2));
ylim([-80 -20]);
xlabel('x, km');

% xx = dlmread('ASA_wedge_kj/TLx_AT.txt');
% plot(xx(:,1)/1000,xx(:,2));


figure;
hold all;
ix = find(x>=10000,1,'first');
plot(y/1000,TLxy(:,ix),'linewidth',1.5,'color','r');
xx = dlmread('ASA_wedge_kj/TLy_x10km');
plot((xx(:,1)/1000)-4,xx(:,2)-20*log10(4*pi),'linewidth',1.5,'color','k','linestyle','--');
xlabel('y, km');
ylabel('TL, dB re 1 m');
grid on;

figure;
hold all;
ix = find(x>=25000,1,'first');
plot(y/1000,TLxy(:,ix),'linewidth',1.5,'color','r');
xx = dlmread('ASA_wedge_kj/TLy_x25km');
plot((xx(:,1)/1000)-4,xx(:,2)-20*log10(4*pi),'linewidth',1.5,'color','k','linestyle','--');
xlabel('y, km');
ylabel('TL, dB re 1 m');
grid on;

