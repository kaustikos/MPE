function [p,q]=pade0fromT(cTaylor,xo,n,m)
% PADE computes the Padè approximant of the function F at the point XO.
%
% [P,Q]=PADE(cTaylor,XO,N,M) returns two polynomial forms, P of order N and Q of order M,
% representing the denominator and the numerator of the rational form 
% which approximates the function F around XO. 
% cTaylor contains taylor coefficients of the approximated function
% N and M must be positive integer values.
%

% Reference paper: Baker, G. A. Jr. Essentials of Padé Approximants in 
% Theoretical Physics. New York: Academic Press, pp. 27-38, 1975.
% see also 
% http://mathworld.wolfram.com/PadeApproximant.html

% This routine has been programmed by Luigi Sanguigno, Phd.
% Affiliation: Italian Institute of Technology
% Modified and slightly corrected by Pavel Petrov PhD
% Affiliation: Il'ichev Pacific Oceanological Inst

% Calculate n+m Taylor coeffs at xo

a=cTaylor.';

% Calculation of Padè coefficients.
% pq=cat(2,cat(1,speye(n+1),zeros(m,n+1)),...
%     spdiags(repmat(reshape(-cat(1,a(end:-1:1),0),1,[]),n+m+1,1),...
%     -(n+m+1):0,n+m+1,n))\a;

% PP: note that the last argument of spdiags should be m

pq=cat(2,cat(1,speye(n+1),zeros(m,n+1)),...
    spdiags(repmat(reshape(-cat(1,a(end:-1:1),0),1,[]),n+m+1,1),...
    -(n+m+1):0,n+m+1,m))\a;

% Rewrite the output as evaluable polynomial forms.
p=shiftpoly(pq(n+1:-1:1),xo); q=shiftpoly(cat(1,pq(end:-1:n+2),1),xo);

end


function ps=shiftpoly(p,xo)
% Displaces the origin of -xo

% Initialize values.
ps=zeros(size(p)); q=1; base=[1;-xo]; ps(end)=p(end);

% Substitute the base polynomial in the original polynomial form.
for n=1:(length(p)-1)
    q=conv(q,base);
    ps(end-n:end)=ps(end-n:end)+p(end-n)*q;
end

end
    
