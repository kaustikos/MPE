function [cP,bP] = pade_sqrt(n)

N = 2*n;

%cTaylor(1:2*n) = -((-1/4).^(1:2*n)).*cumprod(1:2:4*n+1);

cTaylor = ((-1/2).^(1:N-1)).*cumprod(1:2:2*N-3 );


cTaylor = [1 0.5 0.5*cTaylor];

cTaylor = cTaylor./factorial(0:N);

[p,q] = pade0fromT(cTaylor,0,n,n);

cP = -1./roots(p);
bP = -1./roots(q);



