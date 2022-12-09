function Sig=sigma_zero(k,d,Mp,Mn,psi_positive,psi_negative,eulergammaF)

psi_negative=fliplr(psi_negative);
pe=sum(2./(k*d*sin(psi_positive(2:Mp)))+1i./(pi*(1:Mp-1)));
ne=sum(2./(k*d*sin(psi_negative(1:Mn)))+1i./(pi*(1:Mn)));
Sig = -1 - (2*1i/pi) * (eulergammaF + log(k*d/(4*pi))) + 2/(k*d*sin(psi_positive(1))) + pe + ne;    
end
