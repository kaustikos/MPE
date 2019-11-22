close all;
clear variables;
clc;


set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial'); 

cw = 1500;
cb = 1800;

rhow = 1;
rhob = 2;

betab = 0.5;

h0 = 20;
zs = 10;


freq = 150;

omeg = 2*pi*freq;

dz0 = 0.25;
opts.Tgr = 3;
opts.Ngr = 3;
opts.nmod = -1;
opts.Hb = 300;
opts.BotBC = 'D';

dhs = 15:-1:0;

kh = [];

% c(z) = c_0 – (\Delta c) tanh( (z-z_0)/\sigma),

%c0 = 1487.8 m/s, z_0 = 32 m, \sigma = 13 m, \Delta c = 33.17 m/s.

% c0 = 1490;
% z0 = 25;
% dc = 30;
% sigma = 10;
% 
% zssp = 0:2:100;
% cssp = c0 - dc*tanh( (zssp-z0)/sigma );

MP.HydrologyData = [   [0       cw];
                                 [100    cw]];
% 
% figure;
% plot(cssp,zssp);
% set(gca,'Ydir','reverse');

for ii = 1:length(dhs)
    
    dh = dhs(ii);
    
    MP.LayersData = [[0    cw  cw  rhow    rhow       0        0];
                            [h0+dh   cw  cb  rhow    rhob   0   betab]
                             ];
    
    
    if ii == 1
        [krs, wmode, dwmode] = ac_modesr(dz0,MP,freq,opts);
        
        nmod = length(krs);
        opts.nmod = nmod;
        kh = zeros(length(dhs),nmod); 
        
        %mgv0 = ModesGroupVelocitiesPekeris(freq,krs,MP);
        
        
        phihs = zeros(size(kh));
        phihzs = zeros(size(kh));
        phizs = zeros(size(kh));
        z = dz0*(0:size(wmode,1)-1);
        
    else
        
        [krs, wmode] = ac_modesr(dz0,MP,freq,opts);
    end;
    
    nmod = length(krs);  
    
    wnum_im_part = ModesAttCoeffs(dz0,freq,krs,wmode,MP);
    
    
    kh(ii,1:nmod) = krs(1:nmod) + 1i*wnum_im_part;
    
    ih1 = find(z>=h0+dh,1,'first');
    
    izs = find(z>=zs,1,'first');
    
    phihs(ii,1:nmod) = wmode(ih1,1:nmod);
    phihzs(ii,1:nmod) = dwmode(ih1,1:nmod);
    
    phizs(ii,1:nmod) = wmode(izs,1:nmod);
    
    
end;


dlmwrite('kj_canyon_att_150Hz_nocutoff.txt',[h0+dhs.' kh],'delimiter','\t','precision',10);
dlmwrite('phizs_canyon_150Hz_nocutoff.txt',[h0+dhs.' phizs],'delimiter','\t','precision',10);

figure;
hold all;
for ii = 1:size(kh,2);
    plot(h0+dhs,kh(:,ii));
    
end;


figure;
hold all;
for ii = 1:size(kh,2);
    plot(h0+dhs,phizs(:,ii));
    
end;  

