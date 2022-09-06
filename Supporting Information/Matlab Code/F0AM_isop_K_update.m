function K = F0AM_isop_K(Met)
% function K = MCMv331_K(Met)
% calculates generic rate constants used in MCMv3.3.1 mechanism.
%
% INPUTS:
% Met: a structure containing the following variables.
%   T: Terature in K (can be scalar or column vector)
%   M: number density in molec/cm^3 (same size as T)
%   H2O: water concentration in molec/cm3
%
% OUTPUTS:
% K: structure of rate constants. Each is size length(T) x 1
%
% 20080916 GMW
% 20110521 GMW  Updated from 3.1 to 3.2
% 20150603 GMW  Renamed from MCMrateK in UWCMv2.2.
%               Fixed MCM typo in parameterization for KBPPN (FCPPN-->FPPN in last line)
% 20150616 JBK  Updated to MCMv3.3.1
% 20160304 GMW  Changed output from name/value pair to structure, and input to structure.
% 20160321 GMW  Fixed some hole-overs from MCMv3.2 (updated to match documentation now).

NO = 1;
HO2 = 1;
RO2 = 1;
NO3 = 1;
O2 = 210000000;

struct2var(Met)

nk = 18; %number of rate constants
krx = nan(length(T),nk);
Knames = cell(nk,1);
i=0;

%% %%%COMPLEX RATE CONSTANTS%%%%%

i=i+1;
Knames{i}   = 'K_OH_CO'; % From JPL 2006
T3I = 1./T;
KLO1=5.9E-33*(300.*T3I).^(1.4);
KHI1=1.1E-12*(300.*T3I).^(-1.3);
XYRAT1=(KLO1.*M)./KHI1;
BLOG1=log10(XYRAT1);
FEXP1=1.0./(1.0+BLOG1.*BLOG1);
KCO1=KLO1.*M.*0.6.^FEXP1./(1.0+XYRAT1);
KLO2=1.5E-13*(300.*T3I).^(-0.6);
KHI2=2.1E9 *(300.*T3I).^(-6.1);
XYRAT2=KLO2.*M./KHI2;
BLOG2=log10(XYRAT2);
FEXP2=1.0./(1.0+BLOG2.*BLOG2);
KCO2=KLO2.*0.6.^FEXP2./(1.d0+XYRAT2);
KCO=KCO1+KCO2;
krx(:,i) = KCO;

i=i+1;
Knames{i} = 'KMT01'; %o3p + no = no2
K10 = 1.0e-31.*M.*(T./300).^-1.6 ;
K1I = 5.0e-11.*(T./300).^-0.3 ;
KR1 = K10./K1I ;
FC1 = 0.85 ;
NC1 = 0.75-1.27.*(log10(FC1)) ;
F1 = 10.^(log10(FC1)./(1+(log10(KR1)./NC1).^2)) ;
krx(:,i) = (K10.*K1I).*F1./(K10+K1I) ;

i=i+1;
Knames{i} = 'KMT02'; %o3p + no2 = no3
K20 = 1.3e-31.*M.*(T./300).^-1.5 ;
K2I = 2.3e-11.*(T./300).^0.24 ;
KR2 = K20./K2I ;
FC2 = 0.6 ;
NC2 = 0.75-1.27.*(log10(FC2)) ;
F2 = 10.^(log10(FC2)./(1+(log10(KR2)./NC2).^2)) ;
krx(:,i) = (K20.*K2I).*F2./(K20+K2I) ;

i=i+1;
Knames{i} = 'KMT03'; %no2 + no3 = n2o5
K30 = 3.6e-30.*M.*(T./300).^-4.1 ;
K3I = 1.9e-12.*(T./300).^0.2 ;
KR3 = K30./K3I ;
FC3 = 0.35 ;
NC3 = 0.75-1.27.*(log10(FC3)) ;
F3 = 10.^(log10(FC3)./(1+(log10(KR3)./NC3).^2)) ;
krx(:,i) = (K30.*K3I).*F3./(K30+K3I) ;

i=i+1;
Knames{i} = 'KMT04'; %n2o5 = no2 + no3
K40 = 1.3e-3.*M.*(T./300).^-3.5.*exp(-11000./T) ;
K4I = 9.7e+14.*(T./300).^0.1.*exp(-11080./T) ;
KR4 = K40./K4I ;
FC4 = 0.35 ;
NC4 = 0.75-1.27.*(log10(FC4)) ;
F4 = 10.^(log10(FC4)./(1+(log10(KR4)./NC4).^2)) ;
krx(:,i) = (K40.*K4I).*F4./(K40+K4I) ;

i=i+1;
Knames{i} = 'KMT05';
krx(:,i) = 1.44e-13.*(1+(M./4.2e19));

i=i+1;
Knames{i} = 'KMT06';
krx(:,i) = 1 + (1.40e-21.*exp(2200./T).*H2O) ;

i=i+1;
Knames{i} = 'KMT07'; %oh + no = hono
K70 = 7.4e-31.*M.*(T./300).^-2.4 ;
K7I = 3.3e-11.*(T./300).^-0.3 ;
KR7 = K70./K7I ;
FC7 = 0.81 ;
NC7 = 0.75-1.27.*(log10(FC7)) ;
F7 = 10.^(log10(FC7)./(1+(log10(KR7)./NC7).^2)) ;
krx(:,i) = (K70.*K7I).*F7./(K70+K7I) ;

i=i+1;
Knames{i} = 'KMT08'; %oh + no2 = hno3
K80 = 3.2e-30.*M.*(T./300).^-4.5 ;
K8I = 3e-11 ;
KR8 = K80./K8I ;
FC8 = 0.41 ;
NC8 = 0.75-1.27.*(log10(FC8)) ;
F8 = 10.^(log10(FC8)./(1+(log10(KR8)./NC8).^2));
krx(:,i) = (K80.*K8I).*F8./(K80+K8I) ;

i=i+1;
Knames{i} = 'KMT09'; %ho2 + no2 = ho2no2
K90 = 1.4e-31.*M.*(T./300).^-3.1 ;
K9I = 4.0e-12 ;
KR9 = K90./K9I ;
FC9 = 0.4 ;
NC9 = 0.75-1.27.*(log10(FC9)) ;
F9 = 10.^(log10(FC9)./(1+(log10(KR9)./NC9).^2)) ;
krx(:,i) = (K90.*K9I).*F9./(K90+K9I) ;

i=i+1;
Knames{i} = 'KMT10'; %ho2no2 = ho2 + no2
K100 = 4.10e-05.*M.*exp(-10650./T) ;
K10I = 6.0e+15.*exp(-11170./T) ;
KR10 = K100./K10I ;
FC10 = 0.4 ;
NC10 = 0.75-1.27.*(log10(FC10)) ;
F10 = 10.^(log10(FC10)./(1+(log10(KR10)./NC10).^2)) ;
krx(:,i) = (K100.*K10I).*F10./(K100+K10I) ;

i=i+1;
Knames{i} = 'KMT11'; %oh + hno3 = no3
k1 = 2.4e-14.*exp(460./T);
k3 = 6.5e-34.*exp(1335./T);
k4 = 2.7e-17.*exp(2199./T);
k2 = k3.*M./(1 + k3.*M./k4);
krx(:,i) = k1 + k2;

i=i+1;
Knames{i} = 'KMT12';
K120 = 2.5e-31.*M.*(T./300).^-2.6 ;
K12I = 2e-12 ;
KR12 = K120./K12I ;
FC12 = 0.53 ;
NC12 = 0.75-1.27.*(log10(FC12)) ;
F12 = 10.^(log10(FC12)./(1.0+(log10(KR12)./NC12).^2)) ;
krx(:,i) = (K120.*K12I.*F12)./(K120+K12I) ;

i=i+1;
Knames{i} = 'KMT13';
K130 = 2.5e-30.*M.*(T./300).^-5.5 ;
K13I = 1.8e-11 ;
KR13 = K130./K13I ;
FC13 = 0.36 ;
NC13 = 0.75-1.27.*(log10(FC13)) ;
F13 = 10.^(log10(FC13)./(1+(log10(KR13)./NC13).^2)) ;
krx(:,i) = (K130.*K13I).*F13./(K130+K13I) ;

i=i+1;
Knames{i} = 'KMT14';
K140 = 9.0e-5.*exp(-9690./T).*M ;
K14I = 1.1e+16.*exp(-10560./T) ;
KR14 = K140./K14I ;
FC14 = 0.36 ;
NC14 = 0.75-1.27.*(log10(FC14)) ;
F14 = 10.^(log10(FC14)./(1+(log10(KR14)./NC14).^2)) ;
krx(:,i) = (K140.*K14I).*F14./(K140+K14I) ;

i=i+1;
Knames{i} = 'KIHOO1';
K0 = (5.05E15).*exp((-12200)./T).*exp(1E8./T.^3);
K1 = (1.79E14).*exp((-8830)./T);
K2 = (9.33E-2).*K0./(K0+K1);
krx(:,i) =  (1.7E-11).*exp((390)./T).*(1-K2);

i=i+1;
Knames{i} = 'KIHOO4';
K0 = (2.22E9).*exp((-7160)./T).*exp(1E8./T.^3);
K1 = (1.75E14).*exp((-9054)./T);
K2 = (2.26E-1).*K0./(K0+K1);
krx(:,i) = (1.0E-11) .* exp((390)./T).*(1.-K2) ;

i=i+1;
Knames{i} = 'KISO1';
K0 = (5.05E15).*exp((-12200)./T).*exp(1E8./T.^3);
K1 = (1.79E14).*exp((-8830)./T);
K2 = (9.33E-2).*K0./(K0+K1);
krx(:,i) = (1.7E-11) .* exp((390)./T).*K2 ;

i=i+1;
Knames{i} = 'KISO4';
K0 = (2.22E9).*exp((-7160)./T).*exp(1E8./T.^3);
K1 = (1.75E14).*exp((-9054)./T);
K2 = (2.26E-1).*K0./(K0+K1);
krx(:,i) = (1.0E-11) .* exp((390)./T).*K2 ;

i=i+1;
Knames{i} = 'KRO2HO2' ;
krx(:,i) = 2.91e-13.*exp(1300/T) ;

i=i+1;
Knames{i} = 'KRO2NO' ;
krx(:,i) = 2.7e-12.*exp(360/T) ;

i=i+1;
Knames{i} = 'KRO2NO3' ;
krx(:,i) = 2.3e-12 ;

i=i+1;
Knames{i} = 'KRO2' ;
krx(:,i) = 1.26e-12.*0.01 ;

i=i+1;
Knames{i} = 'KAPHO2' ;
krx(:,i) = 5.2e-13.*exp(980/T) ;

i=i+1;
Knames{i} = 'KAPNO' ;
krx(:,i) = 7.5e-12.*exp(290/T) ;


i=i+1;
Knames{i} = 'KNO3AL' ;
krx(:,i) = 1.44e-12.*exp(-1862/T) ;

i=i+1;
Knames{i} = 'KDEC' ;
krx(:,i) = 1.00e06 ;

i=i+1;
Knames{i} = 'KROPRIM' ;
krx(:,i) = 2.50e-14.*exp(-300/T) ;

i=i+1;
Knames{i} = 'KROSEC' ;
krx(:,i) = 2.50e-14.*exp(-300/T) ;

i=i+1;
Knames{i} = 'KCH3O2' ;
krx(:,i) = 1.03e-13.*exp(365/T) ;

i=i+1;
Knames{i} = 'K298CH3O2' ;
krx(:,i) = 3.5e-13 ;

i=i+1;
Knames{i} = 'K14ISOM1' ;
krx(:,i) = 3.00e7.*exp(-5300/T) ;

i=i+1;
Knames{i} = 'KBPAN';
KD0 = 1.10e-05.*M.*exp(-10100/T) ;
KDI = 1.90e17.*exp(-14100/T) ;
KRD = KD0/KDI ;
FCD = 0.30 ;
NCD = 0.75-1.27.*(log10(FCD)) ;
FD = 10.^(log10(FCD)/(1+(log10(KRD)/NCD).^2)) ;
krx(:,i) = (KD0.*KDI).*FD/(KD0+KDI) ;

i=i+1;
Knames{i} = 'KFPAN' ;
KC0 = 3.28e-28.*M.*(T/300).^-6.87 ;
KCI = 1.125e-11.*(T/300).^-1.105 ;
KRC = KC0/KCI ;
FCC = 0.30 ;
NC = 0.75-1.27.*(log10(FCC)) ;
FC = 10.^(log10(FCC)/(1+(log10(KRC)/NC).^2)) ;
krx(:,i) = (KC0.*KCI).*FC/(KC0+KCI) ;

i=i+1;
Knames{i} = 'KMT15' ;
K150 = 8.6e-29.*M.*(T/300).^-3.1 ;
K15I = 9.0e-12.*(T/300).^-0.85 ;
KR15 = K150/K15I ;
FC15 = 0.48 ;
NC15 = 0.75-1.27.*(log10(FC15)) ;
F15 = 10.^(log10(FC15)/(1+(log10(KR15)/NC15).^2)) ;
krx(:,i) = (K150.*K15I).*F15/(K150+K15I) ;

i=i+1;
Knames{i} = 'KMT16' ;
K160 = 8e-27.*M.*(T/300).^-3.5 ;
K16I = 3.0e-11.*(T/300).^-1 ;
KR16 = K160/K16I ;
FC16 = 0.5 ;
NC16 = 0.75-1.27.*(log10(FC16)) ;
F16 = 10.^(log10(FC16)/(1+(log10(KR16)/NC16).^2)) ;
krx(:,i) = (K160.*K16I).*F16/(K160+K16I) ;

i=i+1;
Knames{i} = 'KMT17' ;
K170 = 5.0e-30.*M.*(T/300).^-1.5 ;
K17I = 1.0e-12 ;
KR17 = K170/K17I ;
FC17 = 0.17.*exp(-51/T)+exp(-T/204) ;
NC17 = 0.75-1.27.*(log10(FC17)) ;
F17 = 10.^(log10(FC17)/(1.0+(log10(KR17)/NC17).^2)) ;
krx(:,i) = (K170.*K17I.*F17)/(K170+K17I) ;

i=i+1;
Knames{i} = 'KMT18' ;
krx(:,i) = 9.5e-39.*O2.*exp(5270/T)/(1+7.5e-29.*O2.*exp(5610/T)) ;

i=i+1;
Knames{i} = 'KBPPN';
KPPN0 = 1.7e-03.*exp(-11280/T).*M ;
KPPNI = 8.3e16.*exp(-13940/T) ;
KRPPN = KPPN0/KPPNI ;
FCPPN = 0.36 ;
NCPPN = 0.75-1.27.*(log10(FCPPN)) ;
FPPN = 10.^(log10(FCPPN)/(1+(log10(KRPPN)/NCPPN).^2)) ;
krx(:,i) = (KPPN0.*KPPNI).*FCPPN/(KPPN0+KPPNI) ;

i=i+1;
Knames{i} = 'K_HO2_HO2';
K1    =3.00D-13.*exp(460./T) ;
K2    =2.1D-33.*exp(920./T).*M;
krx(:,i) = K1+K2;

i=i+1;
Knames{i} = 'K_HO2_HO2_H2O';
K1    =4.20D-34.*exp(2660./T) ;
K2    =2.94D-54.*exp(3120./T).*M;
krx(:,i) = K1+K2;

i=i+1;
Knames{i} = 'K_O3P_NO';
LPL 	= 9.10D-32.*((T./300).^-1.5).*M ;
HPL 	= 3.00D-11.*((T./300).^-0.0) ;
krx(:,i)=(LPL./(1+LPL./HPL)).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_OH_NO' ;
LPL     =(7.10e-31.*(T./300).^-2.6).*M;
HPL     =(3.60e-11.*(T./300).^-0.1);
krx(:,i) =(LPL./(1+LPL./HPL)).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_O3P_NO2_NO3';
LPL=3.40D-31.*((T./300).^-1.6).*M ;
HPL=2.30D-11.*((T./300).^-0.2) ;
krx(:,i)=(LPL./(1+LPL./HPL)).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_OH_NO2';
LPL     =(1.8e-30.*(T./300).^-3.0).*M;
HPL     =(2.8e-11.*(T./300).^-0.0);
krx(:,i)=(LPL./(1+LPL./HPL)).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_NO2_NO3';
LPL     =(2.40e-30.*(T./300).^-3).*M;
HPL     =(1.60e-12.*(T./300).^0.1);
krx(:,i) =(LPL./(1+LPL./HPL)).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_HO2_NO2'; 
LPL=1.90D-31.*((T./300).^-3.4).*M ;
HPL=4.00D-12.*((T./300).^-0.3) ;
krx(:,i) =(LPL./(1+LPL./HPL)).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_OH_SO2';
LPL     =(2.90e-31.*(T./300).^-4.1).*M;
HPL     =(1.70e-12.*(T./300).^-0.2);
krx(:,i) =(LPL./(1+LPL./HPL)).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_ETHE_OH';
LPL		= (1.00e-28.*(T./300).^-4.5).*M;
HPL		= (8.80e-12.*(T./300).^-0.85);
krx(:,i) =(LPL./(1+LPL./HPL) ).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_ACYE_OH';
LPL		= (5.50e-30.*(T./300).^-0.0).*M;
HPL		= (8.30e-13.*(T./300).^2.0);
krx(:,i) = (LPL./(1+LPL./HPL) ).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_MECO3_NO2';
LPL		= (9.70e-29.*(T./300).^-5.6).*M;
HPL		= (9.30e-12.*(T./300).^-1.5);
krx(:,i) =(LPL./(1+LPL./HPL) ).*0.6.^(1./(1+(log10(LPL./HPL)).^2));

i=i+1;
Knames{i} = 'K_N2O5';
LPL		= (2.40e-30.*(T./300).^-3.0).*M;
HPL		= (1.60e-12.*(T./300).^0.1);
KF =(LPL./(1+LPL./HPL) ).*0.6.^(1./(1+(log10(LPL./HPL)).^2));
A=1.72e26;
B=-10840;
krx(:,i) = KF.*(A.*exp(B./T));

i=i+1;
Knames{i} = 'K_HNO4';
LPL		= (1.90e-31.*(T./300).^-3.4).*M;
HPL		= (4.00e-12.*(T./300).^-0.3);
KF =(LPL./(1+LPL./HPL) ).*0.6.^(1./(1+(log10(LPL./HPL)).^2));
A=4.76e26;
B=-10900;
krx(:,i) = KF.*(A.*exp(B./T));

i=i+1;
Knames{i} = 'K_PAN';
LPL		= (9.70e-29.*(T./300).^-5.6).*M;
HPL		= (9.30e-12.*(T./300).^-1.5);
KF =(LPL./(1+LPL./HPL) ).*0.6.^(1./(1+(log10(LPL./HPL)).^2));
A=1.11e28;
B=-14000;
krx(:,i) = KF.*(A.*exp(B./T));

i=i+1;
Knames{i} = 'K_HO2_NO_HNO3';
K1 = ((6.095e-14).*(exp(270./T)).*((T./300).^-1));
K2 = ((6.857e-34).*(exp(270./T)).*((T./300).^1)).*M;
K3 = ((-5.968e-14).*(exp(270./T)));
krx(:,i) = K1+K2+K3;

i=i+1;
Knames{i} = 'K_OH_HNO3';
K0    =2.40D-14.*exp(460./T) ;
K2    =2.70D-17.*exp(2199./T) ;
K3    =6.50D-34.*exp(1335./T).*M ;
K1    = K3./(1+(K3./K2)) ;
krx(:,i) =K0+K1 ;

i=i+1;
Knames{i} = 'K_OH_CO';
K1    =1.44D-13.*exp(0./T) ;
K2    =2.74D-33.*exp(0./T).*M;
krx(:,i) = K1+K2;

%% accumulate
K = struct;
for i=1:length(Knames)
    K.(Knames{i}) = krx(:,i);
end
