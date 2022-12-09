function Sig=sigma_odd_faster(k,Beta,d,n,Mp,Mn,psi_positive,psi_negative)
%Remember_Bernie=memoize(@bernoulli_poly2);

pe=sum(exp(1i*(2*(n+1)/2-1)*psi_positive(1:Mp))./(k*d*sin(psi_positive(1:Mp))));
ne=sum(exp(-1i*(2*(n+1)/2-1)*psi_negative(1:Mn))./(k*d*sin(psi_negative(1:Mn))));

f=0;
for m=0:(n+1)/2-1
    f=f+(-1)^(m)*2^(2*m)*factorial((n+1)/2+m-1)/(factorial(2*m+1)*factorial((n+1)/2-m-1))*(2*pi/(k*d))^(2*m+1)*bernoulli_poly2(2*m+1,Beta*d/(2*pi));
end

Sig=2*(-1)^((n+1)/2)*1i*(pe+ne)+(2/pi)*f;

end