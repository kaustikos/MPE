function U_out = WAPEsolveSSP_PML_SelfStarter(PEpars)

% SSP solver for wide-angle parabolic equation in Cartesian coordinates
% self-starter algorithm is used. This code can be used for initializing
% the wide-angle PE solution

% The WAPE is solved for the quantity U
% WAPE is one-way approximation for the 2D Helmholtz equation
% U_{xx} + U_{yy} + k^2 U = 0
% in this function it is assumed that k = k(y)

% Perfectly matching layers (PMLs) are used for truncating domain in y at
% the endpoints PEpars.y(1) and PEpars.y(end)

% the input is supplied via the cell array PEpars

% Godin's energy-conserving form of the SSP algorithm is used

x = PEpars.x;
y_o = PEpars.y;
dx = x(2)-x(1);
dy = y_o(2)-y_o(1);
nx = length(x);
ny = length(y_o);

%% PML
% PML is basically imaginary part of b!
% we replace b with b1 = b + Q, where Q is
% Q = 0 for z_l < z < z_r
% Q = Q(z) = -i\sigma (z - z_r)^3, z > z_r
% z_r is the right boundary of the domain
% and similarly for z < z_l


sigma = 10;
if isfield(PEpars,'sigma');
    sigma = PEpars.sigma;
end;

lThickness = 100;
if isfield(PEpars,'lThickness');
    lThickness = PEpars.lThickness;
end;

nPML = max( 1,fix(lThickness/dy) );

ny_o = ny;
ny = ny + 2*nPML;

y(1:ny,1) = 0;
y(nPML+1:nPML+ny_o,1) = y_o(1:ny_o);
y(1:nPML,1) = y(nPML+1) - (nPML:-1:1)*dy;

y(nPML+ny_o+1:ny,1) = y(nPML+ny_o) + (1:nPML)*dy;

%iy0 = find(y>=0,1,'first');

% the 'real' solution is within nPML+1:nPML+ny_o
% the PML is 1:nPML and ny_o+nPML+1:ny_o+2*PML

if nPML > 1
    lShape(1:nPML,1) = sigma*(0:nPML-1).^3/((nPML - 1)^3);
    PML(1:ny,1) = 0;
    PML(ny-nPML+1:ny,1) = lShape(1:nPML,1);
    PML(1:nPML,1) = flipdim(lShape(1:nPML,1), 1);
end;


%% 

k2 = PEpars.k2;
k02 = PEpars.k02;
k2c = zeros(ny,1);

dLower(1:ny,1) = 0;
dUpper(1:ny,1) = 0;
dMain(1:ny,1) = 0;

IdM =  spdiags(ones(ny,1),0,ny,ny);


%% compute pade coeff exp(i*k0*h) -- the propagator of SSP method

nP = PEpars.nP;

theta = 1;

if isfield(PEpars,'theta');
    theta = PEpars.theta;
end;

hk0 = sqrt(k02)*dx;

[p1,q1] = pade_expsqrt(hk0,nP,nP);

[p2,q2] = pade_expsqrt(hk0,nP-1,nP);

p1 = p1/exp(1i*hk0);
p2 = p2/exp(1i*hk0);

qP = theta*q1 + (1-theta)*q2;
pP = theta*p1 + (1-theta)*[0; p2];

cP = -1./roots(pP);
bP = -1./roots(qP);

aP = zeros(size(bP));

for ii = 1:nP
    aP(ii) = (cP(ii) - bP(ii))*prod(   (cP([1:ii-1 ii+1:nP]) - bP(ii))./(bP([1:ii-1  ii+1:nP]) - bP(ii))  );
end;

if nP== 1
    aP(1) = cP(1) - bP(1);
end

dP = - aP./bP;
d0 = 1 - sum(dP);


%% compute Pade coeffs for (1+X)^(1/4) and (1+X)^(-1/4)

% this is necessary for Godin's energy-conserving approach

nP1 = nP;

[d01,dP1,bP1] = pade_1by4(nP1);

nP2 = nP;

[d02,dP2,bP2] = pade_m1by4(nP2);


%% The starter

ICU(1:ny,1) = 0;

krs = sqrt(k02);

% SStarter IC (self starter for the homogeneous medium)

ic(1:ny) = exp(  1i*krs * abs(y) )/(4*krs);

% correcting the self-starter for the use with PMLs

yQ = y(ny-nPML+1);
Eps = y(ny) - y(ny-nPML+1);

PMLdecay(1:nPML) = exp( (  (y(ny-nPML+1:ny) - yQ).^4)*krs*sigma/(4*Eps^3)  );

ic(ny-nPML+1:ny) = ic(ny-nPML+1:ny).*PMLdecay(1:nPML);
ic(1:nPML) = ic(1:nPML).*flipdim(PMLdecay(1:nPML), 2);


ICU(1:ny,1) = ic;

% making a transformation U -> V 

% O. A. Godin, Reciprocity and energy conservation within the parabolic approximation, 
% Wave motion 29 (2) (1999) 175-194.

% X. Antoine, Y. Huang, Y. Y. Lu, Computing high-frequency scattered fields by beam propagation methods: A prospective study, 
% Journal of Algorithms & Computational Technology 4 (2) (2010) 147-166.

ICV = d01*ICU;

for jj = 1:nP1
    
    X = FormLMatrixTri(1);
    
    W =  (IdM + bP1(jj)*X )\ICU;
    
    
    ICV(1:ny,1) = ICV(1:ny,1) + dP1(jj)*W;
end;

ICV = (k02^(1/4))*ICV;




%% initializing variables for the SSP solution main loop

V(1:ny,1:nx) = 0;
V(1:ny,1) = ICV;
Vp = V(1:ny,1);

for ii = 2:nx
    
    if mod(ii,1000) == 0
        disp([int2str(ii)  ' of ' int2str(nx)  ' steps'] );
    end;
    
    Vc = d0*Vp;
    
    
    
    for jj = 1:nP
        
        X = FormLMatrixTri(ii-1);
        
        W =  (IdM + bP(jj)*X )\Vp;
        
        
        Vc(1:ny,1) = Vc(1:ny,1) + dP(jj)*W;
    end;
    
    Vp = Vc*exp(1i*hk0);
    
    V(1:ny,ii) = Vp;
    

 end;



%% back-transforming from V to U

U(1:ny,1:nx) = 0;

for ii = 1:nx
    
    U(1:ny,ii) = d02*V(1:ny,ii);
    
    for jj = 1:nP2
        
        X = FormLMatrixTri(ii);
        
        W =  (IdM + bP2(jj)*X )\V(1:ny,ii);
        
        
        U(1:ny,ii) = U(1:ny,ii) + dP2(jj)*W;
    end;
    
    
end;

U = (k02^(-1/4))*U;

% Correction due to Self starter: apply sqrt(1+L)


nP3 = nP;

[cP3,bP3] = pade_sqrt(nP3);


X = FormLMatrixTri(1);

for ii = 1:nx
    
    for jj = 1:nP3;
        W =  (IdM + bP3(jj)*X )\U(1:ny,ii);
        U(1:ny,ii) = (IdM + cP3(jj)*X)*W;
    end;
    
    U(1:ny,ii) = sqrt(k02)*U(1:ny,ii);
    
end;

% output: truncate the PMLs and return only the original domain

U_out(1:ny_o,1:nx) = U(nPML+1:nPML+ny_o,1:nx);


%% A subroutine used for the SSP 3-diagonal matrix of the discretized operator X = L/k_0^2


    function opX = FormLMatrixTri(iis)
    
        
        
        if size(k2,2) > 1
            k2c(nPML+1:nPML+ny_o,1) = k2(1:ny_o,iis);
        else
            k2c(nPML+1:nPML+ny_o,1) = k2(1:ny_o,1);
        end;
        k2c(1:nPML,1) = k2c(nPML+1);
        k2c(nPML+ny_o+1:ny,1) = k2c(nPML+ny_o);
        
        
        
        dLower(2:ny-1,1) = (1/k02)*(1./(1 - 1i*PML(2:ny-1)))./(  1 - 1i*(  PML(2:ny-1) + PML(1:ny-2)  )/2 );
        
        dUpper(2:ny-1,1) = (1/k02)*(1./(1 - 1i*PML(2:ny-1)))./(  1 - 1i*(  PML(2:ny-1) + PML(3:ny)  )/2 );
        
        dMain(2:ny-1,1) = - (1/k02)*(1./(1 - 1i*PML(2:ny-1))).*...
            (  1./(  1 - 1i*(  PML(2:ny-1) + PML(3:ny)  )/2 ) +...
            1./(  1 - 1i*(  PML(2:ny-1) + PML(1:ny-2)  )/2 )  ) + ...
            (dy^2)*(  k2c(2:ny-1)/k02 -1  );
        
        Mc = spdiags([dUpper dMain dLower], -1:1, ny, ny).'; % this is done on purpose, so that in the line j of the matrix Mc there are dL(j), dM(j), dU(j)
        
        Mc(1,1) = 1/k02;
        Mc(ny,ny) = 1/k02;
        
        opX = Mc/(dy^2);
        
    end


end

