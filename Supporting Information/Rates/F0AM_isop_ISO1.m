function K = F0AM_isop_ISO1(T,A0,B0,C0,D0,E0,F0,G0)
K0 = D0.*exp(E0/T).*exp(1e8/(T.^3));
K1 = F0.*exp(G0/T);
K2 = C0.*K0/(K0+K1);
K =  A0.*(exp(B0/T)).*(1-K2);
