function [d0,dP,bP] = pade_m1by4(n)

N = 2*n;

%cTaylor(1:2*n) = -((-1/4).^(1:2*n)).*cumprod(1:2:4*n+1);

cTaylor = ((-1/4).^(0:N)).*[1  cumprod(1:4:1+4*(N-1) )];


%cTaylor = [1 0.25 0.25*cTaylor];

cTaylor = cTaylor./factorial(0:N);

[p,q] = pade0fromT(cTaylor,0,n,n);

cP = -1./roots(p);
bP = -1./roots(q);

aP = zeros(size(bP));

for ii = 1:n
    %aP(ii) = (cP(ii) - bP(ii))*prod(   (cP([[1:ii-1] [ii+1:nP]]) - bP(ii))./(bP([[1:ii-1] [ii+1:nP]]) - bP(ii))  );
    aP(ii) = (cP(ii) - bP(ii))*prod(   (cP([1:ii-1 ii+1:n]) - bP(ii))./(bP([1:ii-1  ii+1:n]) - bP(ii))  );
end;

if n== 1
   aP(1) = cP(1) - bP(1); 
end

dP = - aP./bP;

d0 = 1 - sum(dP);

