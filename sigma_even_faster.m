function Sig=sigma_even_faster(k,Beta,d,n,Mp,Mn,psi_positive,psi_negative)
%Remember_Bernie=memoize(@bernoulli_poly2);

pe=sum(exp(2*1i*(n/2)*psi_positive(1:Mp))./(k*d*sin(psi_positive(1:Mp))));
ne=sum(exp(-2*1i*(n/2)*psi_negative(1:Mn))./(k*d*sin(psi_negative(1:Mn))));

f=0;
for m=0:n/2
    f=f+(-1)^(m)*2^(2*m)*factorial(n/2+m-1)/(factorial(2*m)*factorial(n/2-m))*(2*pi/(k*d))^(2*m)*bernoulli_poly2(2*m,Beta*d/(2*pi));
end

Sig=2*(-1)^(n/2)*(pe+ne)+(1i/pi)*f;

end