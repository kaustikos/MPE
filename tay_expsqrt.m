function cTaylor = tay_expsqrt(hk0,N)

% derivatives for Taylor coeffs

cTaylor(1:N) = 0;
cTaylor(1) = 1;

vp(1:2*(N)) = 0;
vc(1:2*(N)) = 0;
vp(2) = 1i*hk0/2;


for ii = 2:N
    
    cTaylor(ii) = sum(vp)/(factorial(ii-1));
    vc(1:2*(N)) = 0;
    
    for jj = 3:2*ii
        
        vc(jj) = -vp(jj-2)*(jj-3)/2 + vp(jj-1)*1i*hk0/2;
        
    end;
   
    
    
    vp = vc;
end;


cTaylor = cTaylor*exp(1i*hk0);
