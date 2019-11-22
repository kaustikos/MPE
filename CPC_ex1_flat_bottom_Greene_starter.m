clear variables;
close all;

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

xmax = 5000;
dx = 10;
x = 0:dx:xmax;
nx = length(x);

ymax = 4000;
dy = 2;
y = -ymax:dy:ymax;
ny = length(y);

% bottm relief parameters

h0 = 200;

%!!!!! FLAT BOTTOM TEST CASE, example 1
talph = 0;
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

mmax = 3;

% some auxiliary variables

ih1 = find(hy>=kh(end,1),1,'first');
ih2 = find(hy<=kh(1,1),1,'last');

iy0 = find(y>=0,1,'first');

MsolX(1:nx) = 0;
MsolY(1:ny) = 0;
Panalytical(1:ny,1:nx) = 0;


[Xx Yy] = meshgrid(x, y);
Rr = sqrt(Xx.^2 + Yy.^2);


for ii = 1:mmax
    
    disp(['Solving MPEs: ' int2str(ii)  ' of ' int2str(mmax)  ' modes'] );
    
    
    %PEpars.nP = nP;
    
    ky = interp1(kh(:,1),kh(:,ii+1),hy);
    ky(1:ih1) = ky(ih1);
    ky(ih2:end) = ky(ih2);
    
    krs = ky(iy0);
    
    PEpars.k2 = (ky.^2).';
    PEpars.k02 = (krs).^2;
    
    wmode(1:ny,1) = interp1(phizsh(:,1),phizsh(:,ii+1),hy);
    
    wmode_zs = wmode(iy0);
    
    % Greene starter
    
    PEpars.ic(1:ny) = wmode_zs*(1.4467 - 0.8402 * krs^2 * y.^2 ).*exp(  -krs^2 * y.^2 / 1.5256 )/(2*sqrt(pi));
    
    
    % --- SSP solution for WAMPE ---
    
    Aj = WAPEsolveSSP_PML(PEpars);
    
    % ---
    

    
    
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

    
    % reference solution
    
    Panalytical(1:ny,1:nx) = Panalytical(1:ny,1:nx) + 0.25*1i*wmode_zs*repmat(wmode,1,nx).*besselh(0,1,krs*Rr);
    MsolX(1:nx) = MsolX(1:nx) + 0.25*1i*wmode_zs*wmode(iy0,1)*besselh(0,1,krs*x);
    MsolY(1:ny) = MsolY(1:ny) + 0.25*1i*wmode_zs*(wmode(1:ny,1).').*besselh(0,1,krs*sqrt(x(end)^2 + y.^2));
end;

TLxy = 20*log10(abs(4*pi*P));

TLanalytical = 20*log10(abs(4*pi*Panalytical));

figure;
imagesc(x/1000,y/1000,TLxy );
caxis([-80 -50]);
xlabel('x, km');
ylabel('y, km');


figure;
imagesc(x/1000,y/1000,TLanalytical );
caxis([-80 -50]);
xlabel('x, km');
ylabel('y, km');


[Xx Yy] = meshgrid(x,y);
Rc = 3000;
th_c = -pi/2:0.01:pi/2;

xc = Rc*cos(th_c);
yc = Rc*sin(th_c);

TLc = interp2(Xx,Yy,TLanalytical,xc,yc);
TLc_greene = interp2(Xx,Yy,TLxy,xc,yc);

figure; 
hold all; 
plot(180*th_c/pi,TLc,'linewidth',1.5); 
plot(180*th_c/pi,TLc_greene,'linewidth',1.5,'linestyle','--');
xlim([-90 90]);
ylim([-80 -50]);
legend('analytical solution','SSP + Greene starter');

dlmwrite('greene_analyt.txt',[180*th_c.'/pi TLc.' TLc_greene.'],'delimiter','\t','precision',6);



iy0 = find(y>=0,1,'first');
TLtrack = TLxy(iy0,:);
figure;
hold all;
plot(x/1000,TLtrack);
plot(x/1000,trlo(4*pi*MsolX(1:nx)));

ylim([-80 -20]);
xlabel('x, km');


figure;
hold all;
plot(y/1000,TLxy(:,end));
plot(y/1000,trlo(4*pi*MsolY(1:ny)));

grid on;
ylabel('y, km');






