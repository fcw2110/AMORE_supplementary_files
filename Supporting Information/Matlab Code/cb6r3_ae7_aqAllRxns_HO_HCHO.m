%Foam Reactions File based on the mech.def file for the cb6r3_ae7_aq mechanism.
% # of species   =  147
% # of reactions =  343
% file created by Bryan Place

% Set constant species by scaling to air number density
N2  =  0.780800000.*M;
O2  =  0.209500000.*M;
H2  =  0.000000560.*M;
CH4 =  0.000001850.*M;


SpeciesToAdd = {...
'NO2'; 'NO'; 'O'; 'O3'; 'NO3';'HCHO'; ...
'O1D'; 'HO'; 'HO2'; 'H2O2'; 'N2O5'; ...
'HNO3'; 'HONO'; 'PNA'; 'SO2'; 'SULF'; ...
'SULRXN'; 'ACO3'; 'MO2'; 'RO2'; 'PAN'; ...
'PACD'; 'AACD'; 'CXO3'; 'ALD2'; 'XO2H'; ...
'PANX'; 'HCHO'; 'MEPX'; 'MEOH'; 'ROOH'; ...
'XO2'; 'XO2N'; 'NTR1'; 'NTR2'; 'FACD'; ...
'CO'; 'HCO3'; 'ALDX'; 'GLYD'; 'GLY'; ...
'MGLY'; 'ETHA'; 'ETOH'; 'KET'; 'PAR'; ...
'ACET'; 'PRPA'; 'XPRP'; 'XPAR'; 'ROR'; ...
'ETHY'; 'ETH'; 'OLE'; 'IOLE'; 'ISOP'; ...
'ISO2'; 'ISOPRXN'; 'ISPD'; 'INTR'; 'ISPX'; ...
'HPLD'; 'OPO3'; 'EPOX'; 'EPX2'; 'TERP'; ...
'TRPRXN'; 'TERPNRO2'; 'APIN'; 'BENZENE'; 'CRES'; ...
'BZO2'; 'OPEN'; 'BENZRO2'; 'TOL'; 'TO2'; ...
'TOLRO2'; 'XOPN'; 'XYLMN'; 'XLO2'; 'XYLRO2'; ...
'NAPH'; 'PAHRO2'; 'CRO'; 'CAT1'; 'CRON'; ...
'OPAN'; 'ECH4'; 'CL2'; 'CL'; 'HOCL'; ...
'CLO'; 'FMCL'; 'HCL'; 'CLNO2'; 'CLNO3'; ...
'SVAVB2'; 'SVAVB3'; 'SVAVB4'; 'SVAVB1'; 'SESQ'; ...
'SESQRXN'; 'SOAALK'; 'H2NO3PIJ'; 'H2NO3PK'; 'ACLI'; ...
'ACLJ'; 'ACLK'; 'IEPOXP'; 'ASO4J'; 'AISO3J'; ...
'AGLYJ'; 'MTNO3'; 'AMTNO3J'; 'AMTHYDJ'; 'AAVB2J'; ...
'AOLGAJ'; 'AAVB3J'; 'AAVB4J'; 'AISO1J'; 'AOLGBJ'; ...
'AISO2J'; 'ASQTJ'; 'APOCI'; 'APNCOMI'; 'APOCJ'; ...
'APNCOMJ'; 'PCVOC'; 'PCSOARXN'; 'VLVPO1'; 'VSVPO1'; ...
'VSVPO2'; 'VSVPO3'; 'VIVPO1'; 'VLVOO1'; 'VLVOO2'; ...
'VSVOO2'; 'VSVOO3'; 'VSVOO1'; 'HCHO_PRIMARY'; 'ALD2_PRIMARY'; ...
'BUTADIENE13'; 'ACROLEIN'; 'ACRO_PRIMARY'; 'TOLU'; 'HG'; ...
'HGIIAER'; 'HGIIGAS'; };


AddSpecies


%   1, <R1>
i=i+1;
Rnames{   1} = 'NO2 = NO + O ';
k(:,i) = (JNO2_IUPAC10 ); 
Gstr{i,   1}='NO2';
fNO2(i)=fNO2(i)-1.0;
fNO(i)=fNO(i)+  1.000;fO(i)=fO(i)+  1.000;

%   2, <R2>
i=i+1;
Rnames{   2} = 'O + O2 + M = O3 ';
k(:,i) = (  5.6800E-34.*(T./300).^( -2.6000E+00) ).*O2.*M; 
Gstr{i,   1}='O';
fO(i)=fO(i)-1.0;
fO3(i)=fO3(i)+  1.000;

%   3, <R3>
i=i+1;
Rnames{   3} = 'O3 + NO = NO2 ';
k(:,i) = (  1.4000E-12.*exp( -1.3100E+03./T) ); 
Gstr{i,   1}='O3';Gstr{i,   2}='NO';
fO3(i)=fO3(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;

%   4, <R4>
i=i+1;
Rnames{   4} = 'O + NO + M = NO2 ';
k(:,i) = (  1.0000E-31.*(T./300).^( -1.6000E+00) ).*M; 
Gstr{i,   1}='O';Gstr{i,   2}='NO';
fO(i)=fO(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;

%   5, <R5>
i=i+1;
Rnames{   5} = 'O + NO2 = NO ';
k(:,i) = (  5.5000E-12.*exp(  1.8800E+02./T) ); 
Gstr{i,   1}='O';Gstr{i,   2}='NO2';
fO(i)=fO(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fNO(i)=fNO(i)+  1.000;

%   6, <R6>
i=i+1;
Rnames{   6} = 'O + NO2 = NO3 ';
xko =   1.3000E-31.*M.*exp(  0.0000E+00./T).*(T./300).^ -1.5000E+00;
xkinf =   2.3000E-11.*exp(  0.0000E+00./T).*(T./300).^  2.4000E-01;
xn =   1.0000E+00;
F =   6.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='O';Gstr{i,   2}='NO2';
fO(i)=fO(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;

%   7, <R7>
i=i+1;
Rnames{   7} = 'O + O3 =';
k(:,i) = (  8.0000E-12.*exp( -2.0600E+03./T) ); 
Gstr{i,   1}='O';Gstr{i,   2}='O3';
fO(i)=fO(i)-1.0;fO3(i)=fO3(i)-1.0;


%   8, <R8>
i=i+1;
Rnames{   8} = 'O3 = O ';
k(:,i) = (JO3_O3P_IUPAC10 ); 
Gstr{i,   1}='O3';
fO3(i)=fO3(i)-1.0;
fO(i)=fO(i)+  1.000;

%   9, <R9>
i=i+1;
Rnames{   9} = 'O3 = O1D ';
k(:,i) = (JO3_O1D_IUPAC10 ); 
Gstr{i,   1}='O3';
fO3(i)=fO3(i)-1.0;
fO1D(i)=fO1D(i)+  1.000;

%  10, <R10>
i=i+1;
Rnames{  10} = 'O1D + M = O ';
k(:,i) = (  2.2300E-11.*exp(  1.1500E+02./T) ).*M; 
Gstr{i,   1}='O1D';
fO1D(i)=fO1D(i)-1.0;
fO(i)=fO(i)+  1.000;

%  11, <R11>
i=i+1;
Rnames{  11} = 'O1D + H2O = 2.00000* HO ';
k(:,i) = (  2.1400E-10 ).*H2O; 
Gstr{i,   1}='O1D';
fO1D(i)=fO1D(i)-1.0;
fHO(i)=fHO(i)+  2.000;

%  12, <R12>
i=i+1;
Rnames{  12} = 'O3 +  HO = HO2 ';
k(:,i) = (  1.7000E-12.*exp( -9.4000E+02./T) ); 
Gstr{i,   1}='O3';Gstr{i,   2}='HO';
fO3(i)=fO3(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;

%  13, <R13>
i=i+1;
Rnames{  13} = 'O3 + HO2 =  HO ';
k(:,i) = (  2.0300E-16.*exp(  6.9300E+02./T).*(T./300).^(  4.5700E+00 ) ); 
Gstr{i,   1}='O3';Gstr{i,   2}='HO2';
fO3(i)=fO3(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHO(i)=fHO(i)+  1.000;

%  14, <R14>
i=i+1;
Rnames{  14} = ' HO + O = HO2 ';
k(:,i) = (  2.4000E-11.*exp(  1.1000E+02./T) ); 
Gstr{i,   1}='HO';Gstr{i,   2}='O';
fHO(i)=fHO(i)-1.0;fO(i)=fO(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;

%  15, <R15>
i=i+1;
Rnames{  15} = 'HO2 + O =  HO ';
k(:,i) = (  2.7000E-11.*exp(  2.2400E+02./T) ); 
Gstr{i,   1}='HO2';Gstr{i,   2}='O';
fHO2(i)=fHO2(i)-1.0;fO(i)=fO(i)-1.0;
fHO(i)=fHO(i)+  1.000;

%  16, <R16>
i=i+1;
Rnames{  16} = ' HO +  HO = O ';
k(:,i) = (  6.2000E-14.*exp(  9.4500E+02./T).*(T./300).^(  2.6000E+00 ) ); 
Gstr{i,   1}='HO';Gstr{i,   2}='HO';
fHO(i)=fHO(i)-1.0;fHO(i)=fHO(i)-1.0;
fO(i)=fO(i)+  1.000;

%  17, <R17>
i=i+1;
Rnames{  17} = ' HO +  HO = H2O2 ';
xko =   6.9000E-31.*M.*exp(  0.0000E+00./T).*(T./300).^ -8.0000E-01;
xkinf =   2.6000E-11.*exp(  0.0000E+00./T).*(T./300).^  0.0000E+00;
xn =   1.1300E+00;
F =   5.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='HO';Gstr{i,   2}='HO';
fHO(i)=fHO(i)-1.0;fHO(i)=fHO(i)-1.0;
fH2O2(i)=fH2O2(i)+  1.000;

%  18, <R18>
i=i+1;
Rnames{  18} = ' HO + HO2 =';
k(:,i) = (  4.8000E-11.*exp(  2.5000E+02./T) ); 
Gstr{i,   1}='HO';Gstr{i,   2}='HO2';
fHO(i)=fHO(i)-1.0;fHO2(i)=fHO2(i)-1.0;


%  19, <R19>
i=i+1;
Rnames{  19} = 'HO2 + HO2 = H2O2 ';
xk0 =   2.2000E-13.*exp(  6.0000E+02./T);
xk1 =   1.9000E-33.*exp(  9.8000E+02./T);
k(:,i) = (xk0+xk1.*M ); 
Gstr{i,   1}='HO2';Gstr{i,   2}='HO2';
fHO2(i)=fHO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fH2O2(i)=fH2O2(i)+  1.000;

%  20, <R20>
i=i+1;
Rnames{  20} = 'HO2 + HO2 + H2O = H2O2 ';
xk0 =   3.0800E-34.*exp(  2.8000E+03./T);
xk1 =   2.6600E-54.*exp(  3.1800E+03./T);
k(:,i) = (xk0+xk1.*M ).*H2O; 
Gstr{i,   1}='HO2';Gstr{i,   2}='HO2';
fHO2(i)=fHO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fH2O2(i)=fH2O2(i)+  1.000;

%  21, <R21>
i=i+1;
Rnames{  21} = 'H2O2 = 2.00000* HO ';
k(:,i) = (JH2O2_IUPAC10 ); 
Gstr{i,   1}='H2O2';
fH2O2(i)=fH2O2(i)-1.0;
fHO(i)=fHO(i)+  2.000;

%  22, <R22>
i=i+1;
Rnames{  22} = 'H2O2 +  HO = HO2 ';
k(:,i) = (  2.9000E-12.*exp( -1.6000E+02./T) ); 
Gstr{i,   1}='H2O2';Gstr{i,   2}='HO';
fH2O2(i)=fH2O2(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;

%  23, <R23>
i=i+1;
Rnames{  23} = 'H2O2 + O =  HO + HO2 ';
k(:,i) = (  1.4000E-12.*exp( -2.0000E+03./T) ); 
Gstr{i,   1}='H2O2';Gstr{i,   2}='O';
fH2O2(i)=fH2O2(i)-1.0;fO(i)=fO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

%  24, <R24>
i=i+1;
Rnames{  24} = 'NO + NO + O2 = 2.00000*NO2 ';
k(:,i) = (  3.3000E-39.*exp(  5.3000E+02./T) ).*O2; 
Gstr{i,   1}='NO';Gstr{i,   2}='NO';
fNO(i)=fNO(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  2.000;

%  25, <R25>
i=i+1;
Rnames{  25} = 'HO2 + NO =  HO + NO2 ';
k(:,i) = (  3.4500E-12.*exp(  2.7000E+02./T) ); 
Gstr{i,   1}='HO2';Gstr{i,   2}='NO';
fHO2(i)=fHO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

%  26, <R26>
i=i+1;
Rnames{  26} = 'NO2 + O3 = NO3 ';
k(:,i) = (  1.4000E-13.*exp( -2.4700E+03./T) ); 
Gstr{i,   1}='NO2';Gstr{i,   2}='O3';
fNO2(i)=fNO2(i)-1.0;fO3(i)=fO3(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;

%  27, <R27>
i=i+1;
Rnames{  27} = 'NO3 = NO2 + O ';
k(:,i) = (JNO3NO2_06 ); 
Gstr{i,   1}='NO3';
fNO3(i)=fNO3(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;fO(i)=fO(i)+  1.000;

%  28, <R28>
i=i+1;
Rnames{  28} = 'NO3 = NO ';
k(:,i) = (JNO3NO_06 ); 
Gstr{i,   1}='NO3';
fNO3(i)=fNO3(i)-1.0;
fNO(i)=fNO(i)+  1.000;

%  29, <R29>
i=i+1;
Rnames{  29} = 'NO3 + NO = 2.00000*NO2 ';
k(:,i) = (  1.8000E-11.*exp(  1.1000E+02./T) ); 
Gstr{i,   1}='NO3';Gstr{i,   2}='NO';
fNO3(i)=fNO3(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  2.000;

%  30, <R30>
i=i+1;
Rnames{  30} = 'NO3 + NO2 = NO + NO2 ';
k(:,i) = (  4.5000E-14.*exp( -1.2600E+03./T) ); 
Gstr{i,   1}='NO3';Gstr{i,   2}='NO2';
fNO3(i)=fNO3(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fNO(i)=fNO(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

%  31, <R31>
i=i+1;
Rnames{  31} = 'NO3 + O = NO2 ';
k(:,i) = (  1.7000E-11 ); 
Gstr{i,   1}='NO3';Gstr{i,   2}='O';
fNO3(i)=fNO3(i)-1.0;fO(i)=fO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;

%  32, <R32>
i=i+1;
Rnames{  32} = 'NO3 +  HO = HO2 + NO2 ';
k(:,i) = (  2.0000E-11 ); 
Gstr{i,   1}='NO3';Gstr{i,   2}='HO';
fNO3(i)=fNO3(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

%  33, <R33>
i=i+1;
Rnames{  33} = 'NO3 + HO2 =  HO + NO2 ';
k(:,i) = (  4.0000E-12 ); 
Gstr{i,   1}='NO3';Gstr{i,   2}='HO2';
fNO3(i)=fNO3(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHO(i)=fHO(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

%  34, <R34>
i=i+1;
Rnames{  34} = 'NO3 + O3 = NO2 ';
k(:,i) = (  1.0000E-17 ); 
Gstr{i,   1}='NO3';Gstr{i,   2}='O3';
fNO3(i)=fNO3(i)-1.0;fO3(i)=fO3(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;

%  35, <R35>
i=i+1;
Rnames{  35} = 'NO3 + NO3 = 2.00000*NO2 ';
k(:,i) = (  8.5000E-13.*exp( -2.4500E+03./T) ); 
Gstr{i,   1}='NO3';Gstr{i,   2}='NO3';
fNO3(i)=fNO3(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO2(i)=fNO2(i)+  2.000;

%  36, <R36>
i=i+1;
Rnames{  36} = 'NO3 + NO2 = N2O5 ';
xko =   3.6000E-30.*M.*exp(  0.0000E+00./T).*(T./300).^ -4.1000E+00;
xkinf =   1.9000E-12.*exp(  0.0000E+00./T).*(T./300).^  2.0000E-01;
xn =   1.3300E+00;
F =   3.5000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='NO3';Gstr{i,   2}='NO2';
fNO3(i)=fNO3(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fN2O5(i)=fN2O5(i)+  1.000;

%  37, <R37>
i=i+1;
Rnames{  37} = 'N2O5 = NO3 + NO2 ';
xko =   1.3000E-03.*M.*exp( -1.1000E+04./T).*(T./300).^ -3.5000E+00;
xkinf =   9.7000E+14.*exp( -1.1080E+04./T).*(T./300).^  1.0000E-01;
xn =   1.3300E+00;
F =   3.5000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='N2O5';
fN2O5(i)=fN2O5(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

%  38, <R38>
i=i+1;
Rnames{  38} = 'N2O5 = NO2 + NO3 ';
k(:,i) = (JN2O5_IUPAC10 ); 
Gstr{i,   1}='N2O5';
fN2O5(i)=fN2O5(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;fNO3(i)=fNO3(i)+  1.000;

%  39, <R39>
i=i+1;
Rnames{  39} = 'N2O5 + H2O = 2.00000*HNO3 ';
k(:,i) = (  1.0000E-22 ).*H2O; 
Gstr{i,   1}='N2O5';
fN2O5(i)=fN2O5(i)-1.0;
fHNO3(i)=fHNO3(i)+  2.000;

%  40, <R40>
i=i+1;
Rnames{  40} = 'NO +  HO = HONO ';
xko =   7.4000E-31.*M.*exp(  0.0000E+00./T).*(T./300).^ -2.4000E+00;
xkinf =   3.3000E-11.*exp(  0.0000E+00./T).*(T./300).^ -3.0000E-01;
xn =   8.7000E-01;
F =   8.1000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='NO';Gstr{i,   2}='HO';
fNO(i)=fNO(i)-1.0;fHO(i)=fHO(i)-1.0;
fHONO(i)=fHONO(i)+  1.000;

%  41, <R41>
i=i+1;
Rnames{  41} = 'NO + NO2 + H2O = 2.00000*HONO ';
k(:,i) = (  5.0000E-40 ).*H2O; 
Gstr{i,   1}='NO';Gstr{i,   2}='NO2';
fNO(i)=fNO(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fHONO(i)=fHONO(i)+  2.000;

%  42, <R42>
i=i+1;
Rnames{  42} = 'HONO + HONO = NO + NO2 ';
k(:,i) = (  1.0000E-20 ); 
Gstr{i,   1}='HONO';Gstr{i,   2}='HONO';
fHONO(i)=fHONO(i)-1.0;fHONO(i)=fHONO(i)-1.0;
fNO(i)=fNO(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

%  43, <R43>
i=i+1;
Rnames{  43} = 'HONO = NO +  HO ';
k(:,i) = (JHONO_IUPAC10 ); 
Gstr{i,   1}='HONO';
fHONO(i)=fHONO(i)-1.0;
fNO(i)=fNO(i)+  1.000;fHO(i)=fHO(i)+  1.000;

%  44, <R44>
i=i+1;
Rnames{  44} = 'HONO +  HO = NO2 ';
k(:,i) = (  2.5000E-12.*exp(  2.6000E+02./T) ); 
Gstr{i,   1}='HONO';Gstr{i,   2}='HO';
fHONO(i)=fHONO(i)-1.0;fHO(i)=fHO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;

%  45, <R45>
i=i+1;
Rnames{  45} = 'NO2 +  HO = HNO3 ';
xko =   1.8000E-30.*M.*exp(  0.0000E+00./T).*(T./300).^ -3.0000E+00;
xkinf =   2.8000E-11.*exp(  0.0000E+00./T).*(T./300).^  0.0000E+00;
xn =   1.0000E+00;
F =   6.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='NO2';Gstr{i,   2}='HO';
fNO2(i)=fNO2(i)-1.0;fHO(i)=fHO(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;

%  46, <R46>
i=i+1;
Rnames{  46} = 'HNO3 +  HO = NO3 ';
 xk0 =   2.4000E-14.*exp(  4.6000E+02./T);
 xk2 =   2.7000E-17.*exp(  2.1990E+03./T);
 xk3 =   6.5000E-34.*exp(  1.3350E+03./T);
k(:,i) = (xk0+xk3.*M./(1.0+xk3.*M./xk2) ); 
Gstr{i,   1}='HNO3';Gstr{i,   2}='HO';
fHNO3(i)=fHNO3(i)-1.0;fHO(i)=fHO(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;

%  47, <R47>
i=i+1;
Rnames{  47} = 'HNO3 =  HO + NO2 ';
k(:,i) = (JHNO3_IUPAC10 ); 
Gstr{i,   1}='HNO3';
fHNO3(i)=fHNO3(i)-1.0;
fHO(i)=fHO(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

%  48, <R48>
i=i+1;
Rnames{  48} = 'HO2 + NO2 = PNA ';
xko =   1.8000E-31.*M.*exp(  0.0000E+00./T).*(T./300).^ -3.2000E+00;
xkinf =   4.7000E-12.*exp(  0.0000E+00./T).*(T./300).^  0.0000E+00;
xn =   1.0000E+00;
F =   6.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='HO2';Gstr{i,   2}='NO2';
fHO2(i)=fHO2(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fPNA(i)=fPNA(i)+  1.000;

%  49, <R49>
i=i+1;
Rnames{  49} = 'PNA = HO2 + NO2 ';
xko =   4.1000E-05.*M.*exp( -1.0650E+04./T).*(T./300).^  0.0000E+00;
xkinf =   4.8000E+15.*exp( -1.1170E+04./T).*(T./300).^  0.0000E+00;
xn =   1.0000E+00;
F =   6.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='PNA';
fPNA(i)=fPNA(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

%  50, <R50>
i=i+1;
Rnames{  50} = 'PNA = 0.59000*HO2 +  0.59000*NO2 +  0.41000* HO +  0.41000*NO3 ';
k(:,i) = (JPNA_IUPAC10 ); 
Gstr{i,   1}='PNA';
fPNA(i)=fPNA(i)-1.0;
fHO2(i)=fHO2(i)+  0.590;fNO2(i)=fNO2(i)+  0.590;fHO(i)=fHO(i)+  0.410;fNO3(i)=fNO3(i)+  0.410;

%  51, <R51>
i=i+1;
Rnames{  51} = 'PNA +  HO = NO2 ';
k(:,i) = (  3.2000E-13.*exp(  6.9000E+02./T) ); 
Gstr{i,   1}='PNA';Gstr{i,   2}='HO';
fPNA(i)=fPNA(i)-1.0;fHO(i)=fHO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;

%  52, <R52>
i=i+1;
Rnames{  52} = 'SO2 +  HO = SULF + HO2 + SULRXN ';
xko =   4.5000E-31.*M.*exp(  0.0000E+00./T).*(T./300).^ -3.9000E+00;
xkinf =   1.3000E-12.*exp(  0.0000E+00./T).*(T./300).^ -7.0000E-01;
xn =   1.1000E+00;
F =   5.3000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='SO2';Gstr{i,   2}='HO';
fSO2(i)=fSO2(i)-1.0;fHO(i)=fHO(i)-1.0;
fSULF(i)=fSULF(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fSULRXN(i)=fSULRXN(i)+  1.000;

%  53, <R53>
i=i+1;
Rnames{  53} = 'ACO3 + NO = NO2 + MO2 + RO2 ';
k(:,i) = (  7.5000E-12.*exp(  2.9000E+02./T) ); 
Gstr{i,   1}='ACO3';Gstr{i,   2}='NO';
fACO3(i)=fACO3(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

%  54, <R54>
i=i+1;
Rnames{  54} = 'ACO3 + NO2 = PAN ';
xko =   2.7000E-28.*M.*exp(  0.0000E+00./T).*(T./300).^ -7.1000E+00;
xkinf =   1.2000E-11.*exp(  0.0000E+00./T).*(T./300).^ -9.0000E-01;
xn =   1.4100E+00;
F =   3.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='ACO3';Gstr{i,   2}='NO2';
fACO3(i)=fACO3(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fPAN(i)=fPAN(i)+  1.000;

%  55, <R55>
i=i+1;
Rnames{  55} = 'PAN = NO2 + ACO3 ';
xko =   4.9000E-03.*M.*exp( -1.2100E+04./T).*(T./300).^  0.0000E+00;
xkinf =   5.4000E+16.*exp( -1.3830E+04./T).*(T./300).^  0.0000E+00;
xn =   1.4100E+00;
F =   3.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='PAN';
fPAN(i)=fPAN(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;fACO3(i)=fACO3(i)+  1.000;

%  56, <R56>
i=i+1;
Rnames{  56} = 'PAN = 0.60000*NO2 +  0.60000*ACO3 +  0.40000*NO3 +  0.40000*MO2 +  0.40000*RO2 ';
k(:,i) = (JPAN_IUPAC10 ); 
Gstr{i,   1}='PAN';
fPAN(i)=fPAN(i)-1.0;
fNO2(i)=fNO2(i)+  0.600;fACO3(i)=fACO3(i)+  0.600;fNO3(i)=fNO3(i)+  0.400;fMO2(i)=fMO2(i)+  0.400;fRO2(i)=fRO2(i)+  0.400;

%  57, <R57>
i=i+1;
Rnames{  57} = 'ACO3 + HO2 = 0.41000*PACD +  0.15000*AACD +  0.15000*O3 +  0.44000*MO2 +  0.44000*RO2 +  0.44000* HO ';
k(:,i) = (  5.2000E-13.*exp(  9.8000E+02./T) ); 
Gstr{i,   1}='ACO3';Gstr{i,   2}='HO2';
fACO3(i)=fACO3(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fPACD(i)=fPACD(i)+  0.410;fAACD(i)=fAACD(i)+  0.150;fO3(i)=fO3(i)+  0.150;fMO2(i)=fMO2(i)+  0.440;fRO2(i)=fRO2(i)+  0.440;fHO(i)=fHO(i)+  0.440;

%  58, <R58>
i=i+1;
Rnames{  58} = 'ACO3 + RO2 = ACO3 ';
k(:,i) = (  8.9000E-13.*exp(  8.0000E+02./T) ); 
Gstr{i,   1}='ACO3';Gstr{i,   2}='RO2';
fACO3(i)=fACO3(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fACO3(i)=fACO3(i)+  1.000;

%  59, <R59>
i=i+1;
Rnames{  59} = 'ACO3 + ACO3 = 2.00000*MO2 +  2.00000*RO2 ';
k(:,i) = (  2.9000E-12.*exp(  5.0000E+02./T) ); 
Gstr{i,   1}='ACO3';Gstr{i,   2}='ACO3';
fACO3(i)=fACO3(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fMO2(i)=fMO2(i)+  2.000;fRO2(i)=fRO2(i)+  2.000;

%  60, <R60>
i=i+1;
Rnames{  60} = 'ACO3 + CXO3 = MO2 + ALD2 + XO2H +  2.00000*RO2 ';
k(:,i) = (  2.9000E-12.*exp(  5.0000E+02./T) ); 
Gstr{i,   1}='ACO3';Gstr{i,   2}='CXO3';
fACO3(i)=fACO3(i)-1.0;fCXO3(i)=fCXO3(i)-1.0;
fMO2(i)=fMO2(i)+  1.000;fALD2(i)=fALD2(i)+  1.000;fXO2H(i)=fXO2H(i)+  1.000;fRO2(i)=fRO2(i)+  2.000;

%  61, <R61>
i=i+1;
Rnames{  61} = 'CXO3 + NO = NO2 + ALD2 + XO2H + RO2 ';
k(:,i) = (  6.7000E-12.*exp(  3.4000E+02./T) ); 
Gstr{i,   1}='CXO3';Gstr{i,   2}='NO';
fCXO3(i)=fCXO3(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;fALD2(i)=fALD2(i)+  1.000;fXO2H(i)=fXO2H(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

%  62, <R62>
i=i+1;
Rnames{  62} = 'CXO3 + NO2 = PANX ';
k(:,i) = (k(:,  54) ); 
Gstr{i,   1}='CXO3';Gstr{i,   2}='NO2';
fCXO3(i)=fCXO3(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fPANX(i)=fPANX(i)+  1.000;

%  63, <R63>
i=i+1;
Rnames{  63} = 'PANX = NO2 + CXO3 ';
k(:,i) = (k(:,  55) ); 
Gstr{i,   1}='PANX';
fPANX(i)=fPANX(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;fCXO3(i)=fCXO3(i)+  1.000;

%  64, <R64>
i=i+1;
Rnames{  64} = 'PANX = 0.60000*NO2 +  0.60000*CXO3 +  0.40000*NO3 +  0.40000*ALD2 +  0.40000*XO2H +  0.40000*RO2 ';
k(:,i) = (JPAN_IUPAC10 ); 
Gstr{i,   1}='PANX';
fPANX(i)=fPANX(i)-1.0;
fNO2(i)=fNO2(i)+  0.600;fCXO3(i)=fCXO3(i)+  0.600;fNO3(i)=fNO3(i)+  0.400;fALD2(i)=fALD2(i)+  0.400;fXO2H(i)=fXO2H(i)+  0.400;fRO2(i)=fRO2(i)+  0.400;

%  65, <R65>
i=i+1;
Rnames{  65} = 'CXO3 + HO2 = 0.41000*PACD +  0.15000*AACD +  0.15000*O3 +  0.44000*ALD2 +  0.44000*XO2H +  0.44000*RO2 +  0.44000* HO ';
k(:,i) = (  5.2000E-13.*exp(  9.8000E+02./T) ); 
Gstr{i,   1}='CXO3';Gstr{i,   2}='HO2';
fCXO3(i)=fCXO3(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fPACD(i)=fPACD(i)+  0.410;fAACD(i)=fAACD(i)+  0.150;fO3(i)=fO3(i)+  0.150;fALD2(i)=fALD2(i)+  0.440;fXO2H(i)=fXO2H(i)+  0.440;fRO2(i)=fRO2(i)+  0.440;fHO(i)=fHO(i)+  0.440;

%  66, <R66>
i=i+1;
Rnames{  66} = 'CXO3 + RO2 = 0.80000*ALD2 +  0.80000*XO2H +  0.80000*RO2 ';
k(:,i) = (  8.9000E-13.*exp(  8.0000E+02./T) ); 
Gstr{i,   1}='CXO3';Gstr{i,   2}='RO2';
fCXO3(i)=fCXO3(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fALD2(i)=fALD2(i)+  0.800;fXO2H(i)=fXO2H(i)+  0.800;fRO2(i)=fRO2(i)+  0.800;

%  67, <R67>
i=i+1;
Rnames{  67} = 'CXO3 + CXO3 = 2.00000*ALD2 +  2.00000*XO2H +  2.00000*RO2 ';
k(:,i) = (  3.2000E-12.*exp(  5.0000E+02./T) ); 
Gstr{i,   1}='CXO3';Gstr{i,   2}='CXO3';
fCXO3(i)=fCXO3(i)-1.0;fCXO3(i)=fCXO3(i)-1.0;
fALD2(i)=fALD2(i)+  2.000;fXO2H(i)=fXO2H(i)+  2.000;fRO2(i)=fRO2(i)+  2.000;

%  68, <R68>
i=i+1;
Rnames{  68} = 'RO2 + NO = NO ';
k(:,i) = (  2.4000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='RO2';Gstr{i,   2}='NO';
fRO2(i)=fRO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO(i)=fNO(i)+  1.000;

%  69, <R69>
i=i+1;
Rnames{  69} = 'RO2 + HO2 = HO2 ';
k(:,i) = (  4.8000E-13.*exp(  8.0000E+02./T) ); 
Gstr{i,   1}='RO2';Gstr{i,   2}='HO2';
fRO2(i)=fRO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;

%  70, <R70>
i=i+1;
Rnames{  70} = 'RO2 + RO2 =';
k(:,i) = (  6.5000E-14.*exp(  5.0000E+02./T) ); 
Gstr{i,   1}='RO2';Gstr{i,   2}='RO2';
fRO2(i)=fRO2(i)-1.0;fRO2(i)=fRO2(i)-1.0;


%  71, <R71>
i=i+1;
Rnames{  71} = 'MO2 + NO = HCHO + HO2 + NO2 ';
k(:,i) = (  2.3000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='MO2';Gstr{i,   2}='NO';
fMO2(i)=fMO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fHCHO(i)=fHCHO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

%  72, <R72>
i=i+1;
Rnames{  72} = 'MO2 + HO2 = 0.90000*MEPX +  0.10000*HCHO ';
k(:,i) = (  3.8000E-13.*exp(  7.8000E+02./T) ); 
Gstr{i,   1}='MO2';Gstr{i,   2}='HO2';
fMO2(i)=fMO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fMEPX(i)=fMEPX(i)+  0.900;fHCHO(i)=fHCHO(i)+  0.100;

%  73, <R73>
i=i+1;
Rnames{  73} = 'MO2 + ACO3 = HCHO +  0.90000*HO2 +  0.90000*MO2 +  0.10000*AACD +  0.90000*RO2 ';
k(:,i) = (  2.0000E-12.*exp(  5.0000E+02./T) ); 
Gstr{i,   1}='MO2';Gstr{i,   2}='ACO3';
fMO2(i)=fMO2(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fHCHO(i)=fHCHO(i)+  1.000;fHO2(i)=fHO2(i)+  0.900;fMO2(i)=fMO2(i)+  0.900;fAACD(i)=fAACD(i)+  0.100;fRO2(i)=fRO2(i)+  0.900;

%  74, <R74>
i=i+1;
Rnames{  74} = 'MO2 + RO2 = 0.68500*HCHO +  0.31500*ME HO +  0.37000*HO2 + RO2 ';
k(:,i) = (k(:,  70) ); 
Gstr{i,   1}='MO2';Gstr{i,   2}='RO2';
fMO2(i)=fMO2(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fHCHO(i)=fHCHO(i)+  0.685;fMEOH(i)=fMEOH(i)+  0.315;fHO2(i)=fHO2(i)+  0.370;fRO2(i)=fRO2(i)+  1.000;

%  75, <R75>
i=i+1;
Rnames{  75} = 'XO2H + NO = NO2 + HO2 ';
k(:,i) = (  2.7000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='XO2H';Gstr{i,   2}='NO';
fXO2H(i)=fXO2H(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

%  76, <R76>
i=i+1;
Rnames{  76} = 'XO2H + HO2 = RO HO ';
k(:,i) = (  6.8000E-13.*exp(  8.0000E+02./T) ); 
Gstr{i,   1}='XO2H';Gstr{i,   2}='HO2';
fXO2H(i)=fXO2H(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fROOH(i)=fROOH(i)+  1.000;

%  77, <R77>
i=i+1;
Rnames{  77} = 'XO2H + ACO3 = 0.80000*HO2 +  0.80000*MO2 +  0.20000*AACD +  0.80000*RO2 ';
k(:,i) = (k(:,  58) ); 
Gstr{i,   1}='XO2H';Gstr{i,   2}='ACO3';
fXO2H(i)=fXO2H(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fHO2(i)=fHO2(i)+  0.800;fMO2(i)=fMO2(i)+  0.800;fAACD(i)=fAACD(i)+  0.200;fRO2(i)=fRO2(i)+  0.800;

%  78, <R78>
i=i+1;
Rnames{  78} = 'XO2H + RO2 = 0.60000*HO2 + RO2 ';
k(:,i) = (k(:,  70) ); 
Gstr{i,   1}='XO2H';Gstr{i,   2}='RO2';
fXO2H(i)=fXO2H(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fHO2(i)=fHO2(i)+  0.600;fRO2(i)=fRO2(i)+  1.000;

%  79, <R79>
i=i+1;
Rnames{  79} = 'XO2 + NO = NO2 ';
k(:,i) = (k(:,  75) ); 
Gstr{i,   1}='XO2';Gstr{i,   2}='NO';
fXO2(i)=fXO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;

%  80, <R80>
i=i+1;
Rnames{  80} = 'XO2 + HO2 = RO HO ';
k(:,i) = (k(:,  76) ); 
Gstr{i,   1}='XO2';Gstr{i,   2}='HO2';
fXO2(i)=fXO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fROOH(i)=fROOH(i)+  1.000;

%  81, <R81>
i=i+1;
Rnames{  81} = 'XO2 + ACO3 = 0.80000*MO2 +  0.20000*AACD +  0.80000*RO2 ';
k(:,i) = (k(:,  58) ); 
Gstr{i,   1}='XO2';Gstr{i,   2}='ACO3';
fXO2(i)=fXO2(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fMO2(i)=fMO2(i)+  0.800;fAACD(i)=fAACD(i)+  0.200;fRO2(i)=fRO2(i)+  0.800;

%  82, <R82>
i=i+1;
Rnames{  82} = 'XO2 + RO2 = RO2 ';
k(:,i) = (k(:,  70) ); 
Gstr{i,   1}='XO2';Gstr{i,   2}='RO2';
fXO2(i)=fXO2(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fRO2(i)=fRO2(i)+  1.000;

%  83, <R83>
i=i+1;
Rnames{  83} = 'XO2N + NO = 0.50000*NTR1 +  0.50000*NTR2 ';
k(:,i) = (k(:,  75) ); 
Gstr{i,   1}='XO2N';Gstr{i,   2}='NO';
fXO2N(i)=fXO2N(i)-1.0;fNO(i)=fNO(i)-1.0;
fNTR1(i)=fNTR1(i)+  0.500;fNTR2(i)=fNTR2(i)+  0.500;

%  84, <R84>
i=i+1;
Rnames{  84} = 'XO2N + HO2 = RO HO ';
k(:,i) = (k(:,  76) ); 
Gstr{i,   1}='XO2N';Gstr{i,   2}='HO2';
fXO2N(i)=fXO2N(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fROOH(i)=fROOH(i)+  1.000;

%  85, <R85>
i=i+1;
Rnames{  85} = 'XO2N + ACO3 = 0.80000*HO2 +  0.80000*MO2 +  0.20000*AACD +  0.80000*RO2 ';
k(:,i) = (k(:,  58) ); 
Gstr{i,   1}='XO2N';Gstr{i,   2}='ACO3';
fXO2N(i)=fXO2N(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fHO2(i)=fHO2(i)+  0.800;fMO2(i)=fMO2(i)+  0.800;fAACD(i)=fAACD(i)+  0.200;fRO2(i)=fRO2(i)+  0.800;

%  86, <R86>
i=i+1;
Rnames{  86} = 'XO2N + RO2 = RO2 ';
k(:,i) = (k(:,  70) ); 
Gstr{i,   1}='XO2N';Gstr{i,   2}='RO2';
fXO2N(i)=fXO2N(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fRO2(i)=fRO2(i)+  1.000;

%  87, <R87>
i=i+1;
Rnames{  87} = 'MEPX +  HO = 0.60000*MO2 +  0.60000*RO2 +  0.40000*HCHO +  0.40000* HO ';
k(:,i) = (  5.3000E-12.*exp(  1.9000E+02./T) ); 
Gstr{i,   1}='MEPX';Gstr{i,   2}='HO';
fMEPX(i)=fMEPX(i)-1.0;fHO(i)=fHO(i)-1.0;
fMO2(i)=fMO2(i)+  0.600;fRO2(i)=fRO2(i)+  0.600;fHCHO(i)=fHCHO(i)+  0.400;fHO(i)=fHO(i)+  0.400;

%  88, <R88>
i=i+1;
Rnames{  88} = 'MEPX = MO2 + RO2 +  HO ';
k(:,i) = (JMEPX_IUPAC10 ); 
Gstr{i,   1}='MEPX';
fMEPX(i)=fMEPX(i)-1.0;
fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;fHO(i)=fHO(i)+  1.000;

%  89, <R89>
i=i+1;
Rnames{  89} = 'RO HO +  HO = 0.54000*XO2H +  0.06000*XO2N +  0.60000*RO2 +  0.40000* HO ';
k(:,i) = (  5.3000E-12.*exp(  1.9000E+02./T) ); 
Gstr{i,   1}='ROOH';Gstr{i,   2}='HO';
fROOH(i)=fROOH(i)-1.0;fHO(i)=fHO(i)-1.0;
fXO2H(i)=fXO2H(i)+  0.540;fXO2N(i)=fXO2N(i)+  0.060;fRO2(i)=fRO2(i)+  0.600;fHO(i)=fHO(i)+  0.400;

%  90, <R90>
i=i+1;
Rnames{  90} = 'RO HO = HO2 +  HO ';
k(:,i) = (JMEPX_IUPAC10 ); 
Gstr{i,   1}='ROOH';
fROOH(i)=fROOH(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fHO(i)=fHO(i)+  1.000;

%  91, <R91>
i=i+1;
Rnames{  91} = 'NTR1 +  HO = NTR2 ';
k(:,i) = (  2.0000E-12 ); 
Gstr{i,   1}='NTR1';Gstr{i,   2}='HO';
fNTR1(i)=fNTR1(i)-1.0;fHO(i)=fHO(i)-1.0;
fNTR2(i)=fNTR2(i)+  1.000;

%  92, <R92>
i=i+1;
Rnames{  92} = 'NTR1 = NO2 ';
k(:,i) = (JNTR_IUPAC10 ); 
Gstr{i,   1}='NTR1';
fNTR1(i)=fNTR1(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;

%  93, <R93>
i=i+1;
Rnames{  93} = 'FACD +  HO = HO2 ';
k(:,i) = (  4.5000E-13 ); 
Gstr{i,   1}='FACD';Gstr{i,   2}='HO';
fFACD(i)=fFACD(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;

%  94, <R94>
i=i+1;
Rnames{  94} = 'AACD +  HO = MO2 + RO2 ';
k(:,i) = (  4.0000E-14.*exp(  8.5000E+02./T) ); 
Gstr{i,   1}='AACD';Gstr{i,   2}='HO';
fAACD(i)=fAACD(i)-1.0;fHO(i)=fHO(i)-1.0;
fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

%  95, <R95>
i=i+1;
Rnames{  95} = 'PACD +  HO = ACO3 ';
k(:,i) = (  5.3000E-12.*exp(  1.9000E+02./T) ); 
Gstr{i,   1}='PACD';Gstr{i,   2}='HO';
fPACD(i)=fPACD(i)-1.0;fHO(i)=fHO(i)-1.0;
fACO3(i)=fACO3(i)+  1.000;

%  96, <R96>
i=i+1;
Rnames{  96} = 'HCHO +  HO = HO2 + CO ';
k(:,i) = (  5.4000E-12.*exp(  1.3500E+02./T) ); 
Gstr{i,   1}='HCHO';Gstr{i,   2}='HO';
fHCHO(i)=fHCHO(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fCO(i)=fCO(i)+  1.000;

%  97, <R97>
i=i+1;
Rnames{  97} = 'HCHO = 2.00000*HO2 + CO ';
k(:,i) = (JFORM_R_IUPAC10 ); 
Gstr{i,   1}='HCHO';
fHCHO(i)=fHCHO(i)-1.0;
fHO2(i)=fHO2(i)+  2.000;fCO(i)=fCO(i)+  1.000;

%  98, <R98>
i=i+1;
Rnames{  98} = 'HCHO = CO ';
k(:,i) = (JFORM_M_IUPAC10 ); 
Gstr{i,   1}='HCHO';
fHCHO(i)=fHCHO(i)-1.0;
fCO(i)=fCO(i)+  1.000;

%  99, <R99>
i=i+1;
Rnames{  99} = 'HCHO + O =  HO + HO2 + CO ';
k(:,i) = (  3.4000E-11.*exp( -1.6000E+03./T) ); 
Gstr{i,   1}='HCHO';Gstr{i,   2}='O';
fHCHO(i)=fHCHO(i)-1.0;fO(i)=fO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fCO(i)=fCO(i)+  1.000;

% 100, <R100>
i=i+1;
Rnames{ 100} = 'HCHO + NO3 = HNO3 + HO2 + CO ';
k(:,i) = (  5.5000E-16 ); 
Gstr{i,   1}='HCHO';Gstr{i,   2}='NO3';
fHCHO(i)=fHCHO(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fCO(i)=fCO(i)+  1.000;

% 101, <R101>
i=i+1;
Rnames{ 101} = 'HCHO + HO2 = HCO3 ';
k(:,i) = (  9.7000E-15.*exp(  6.2500E+02./T) ); 
Gstr{i,   1}='HCHO';Gstr{i,   2}='HO2';
fHCHO(i)=fHCHO(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHCO3(i)=fHCO3(i)+  1.000;

% 102, <R102>
i=i+1;
Rnames{ 102} = 'HCO3 = HCHO + HO2 ';
k(:,i) = (  2.4000E+12.*exp( -7.0000E+03./T) ); 
Gstr{i,   1}='HCO3';
fHCO3(i)=fHCO3(i)-1.0;
fHCHO(i)=fHCHO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

% 103, <R103>
i=i+1;
Rnames{ 103} = 'HCO3 + NO = FACD + NO2 + HO2 ';
k(:,i) = (  5.6000E-12 ); 
Gstr{i,   1}='HCO3';Gstr{i,   2}='NO';
fHCO3(i)=fHCO3(i)-1.0;fNO(i)=fNO(i)-1.0;
fFACD(i)=fFACD(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

% 104, <R104>
i=i+1;
Rnames{ 104} = 'HCO3 + HO2 = 0.50000*MEPX +  0.50000*FACD +  0.20000* HO +  0.20000*HO2 ';
k(:,i) = (  5.6000E-15.*exp(  2.3000E+03./T) ); 
Gstr{i,   1}='HCO3';Gstr{i,   2}='HO2';
fHCO3(i)=fHCO3(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fMEPX(i)=fMEPX(i)+  0.500;fFACD(i)=fFACD(i)+  0.500;fHO(i)=fHO(i)+  0.200;fHO2(i)=fHO2(i)+  0.200;

% 105, <R105>
i=i+1;
Rnames{ 105} = 'ALD2 + O = ACO3 +  HO ';
k(:,i) = (  1.8000E-11.*exp( -1.1000E+03./T) ); 
Gstr{i,   1}='ALD2';Gstr{i,   2}='O';
fALD2(i)=fALD2(i)-1.0;fO(i)=fO(i)-1.0;
fACO3(i)=fACO3(i)+  1.000;fHO(i)=fHO(i)+  1.000;

% 106, <R106>
i=i+1;
Rnames{ 106} = 'ALD2 +  HO = ACO3 ';
k(:,i) = (  4.7000E-12.*exp(  3.4500E+02./T) ); 
Gstr{i,   1}='ALD2';Gstr{i,   2}='HO';
fALD2(i)=fALD2(i)-1.0;fHO(i)=fHO(i)-1.0;
fACO3(i)=fACO3(i)+  1.000;

% 107, <R107>
i=i+1;
Rnames{ 107} = 'ALD2 + NO3 = ACO3 + HNO3 ';
k(:,i) = (  1.4000E-12.*exp( -1.8600E+03./T) ); 
Gstr{i,   1}='ALD2';Gstr{i,   2}='NO3';
fALD2(i)=fALD2(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fACO3(i)=fACO3(i)+  1.000;fHNO3(i)=fHNO3(i)+  1.000;

% 108, <R108>
i=i+1;
Rnames{ 108} = 'ALD2 = MO2 + RO2 + CO + HO2 ';
k(:,i) = (JALD2_R_IUPAC10 ); 
Gstr{i,   1}='ALD2';
fALD2(i)=fALD2(i)-1.0;
fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;fCO(i)=fCO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

% 109, <R109>
i=i+1;
Rnames{ 109} = 'ALDX + O = CXO3 +  HO ';
k(:,i) = (  1.3000E-11.*exp( -8.7000E+02./T) ); 
Gstr{i,   1}='ALDX';Gstr{i,   2}='O';
fALDX(i)=fALDX(i)-1.0;fO(i)=fO(i)-1.0;
fCXO3(i)=fCXO3(i)+  1.000;fHO(i)=fHO(i)+  1.000;

% 110, <R110>
i=i+1;
Rnames{ 110} = 'ALDX +  HO = CXO3 ';
k(:,i) = (  4.9000E-12.*exp(  4.0500E+02./T) ); 
Gstr{i,   1}='ALDX';Gstr{i,   2}='HO';
fALDX(i)=fALDX(i)-1.0;fHO(i)=fHO(i)-1.0;
fCXO3(i)=fCXO3(i)+  1.000;

% 111, <R111>
i=i+1;
Rnames{ 111} = 'ALDX + NO3 = CXO3 + HNO3 ';
k(:,i) = (  6.3000E-15 ); 
Gstr{i,   1}='ALDX';Gstr{i,   2}='NO3';
fALDX(i)=fALDX(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fCXO3(i)=fCXO3(i)+  1.000;fHNO3(i)=fHNO3(i)+  1.000;

% 112, <R112>
i=i+1;
Rnames{ 112} = 'ALDX = ALD2 + XO2H + RO2 + CO + HO2 ';
k(:,i) = (JALDX_R_IUPAC10 ); 
Gstr{i,   1}='ALDX';
fALDX(i)=fALDX(i)-1.0;
fALD2(i)=fALD2(i)+  1.000;fXO2H(i)=fXO2H(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;fCO(i)=fCO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

% 113, <R113>
i=i+1;
Rnames{ 113} = 'GLYD +  HO = 0.20000*GLY +  0.20000*HO2 +  0.80000*ACO3 ';
k(:,i) = (  8.0000E-12 ); 
Gstr{i,   1}='GLYD';Gstr{i,   2}='HO';
fGLYD(i)=fGLYD(i)-1.0;fHO(i)=fHO(i)-1.0;
fGLY(i)=fGLY(i)+  0.200;fHO2(i)=fHO2(i)+  0.200;fACO3(i)=fACO3(i)+  0.800;

% 114, <R114>
i=i+1;
Rnames{ 114} = 'GLYD = 0.74000*HCHO +  0.89000*CO +  1.40000*HO2 +  0.15000*ME HO +  0.19000* HO +  0.11000*GLY +  0.11000*XO2H +  0.11000*RO2 ';
k(:,i) = (JGLYD_IUPAC10 ); 
Gstr{i,   1}='GLYD';
fGLYD(i)=fGLYD(i)-1.0;
fHCHO(i)=fHCHO(i)+  0.740;fCO(i)=fCO(i)+  0.890;fHO2(i)=fHO2(i)+  1.400;fMEOH(i)=fMEOH(i)+  0.150;fHO(i)=fHO(i)+  0.190;fGLY(i)=fGLY(i)+  0.110;fXO2H(i)=fXO2H(i)+  0.110;fRO2(i)=fRO2(i)+  0.110;

% 115, <R115>
i=i+1;
Rnames{ 115} = 'GLYD + NO3 = HNO3 + ACO3 ';
k(:,i) = (  1.4000E-12.*exp( -1.8600E+03./T) ); 
Gstr{i,   1}='GLYD';Gstr{i,   2}='NO3';
fGLYD(i)=fGLYD(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;fACO3(i)=fACO3(i)+  1.000;

% 116, <R116>
i=i+1;
Rnames{ 116} = 'GLY +  HO = 1.80000*CO +  0.20000*XO2 +  0.20000*RO2 + HO2 ';
k(:,i) = (  3.1000E-12.*exp(  3.4000E+02./T) ); 
Gstr{i,   1}='GLY';Gstr{i,   2}='HO';
fGLY(i)=fGLY(i)-1.0;fHO(i)=fHO(i)-1.0;
fCO(i)=fCO(i)+  1.800;fXO2(i)=fXO2(i)+  0.200;fRO2(i)=fRO2(i)+  0.200;fHO2(i)=fHO2(i)+  1.000;

% 117, <R117>
i=i+1;
Rnames{ 117} = 'GLY = 2.00000*HO2 +  2.00000*CO ';
k(:,i) = (JGLY_R_IUPAC10 ); 
Gstr{i,   1}='GLY';
fGLY(i)=fGLY(i)-1.0;
fHO2(i)=fHO2(i)+  2.000;fCO(i)=fCO(i)+  2.000;

% 118, <R118>
i=i+1;
Rnames{ 118} = 'GLY + NO3 = HNO3 +  1.50000*CO +  0.50000*XO2 +  0.50000*RO2 + HO2 ';
k(:,i) = (  1.4000E-12.*exp( -1.8600E+03./T) ); 
Gstr{i,   1}='GLY';Gstr{i,   2}='NO3';
fGLY(i)=fGLY(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;fCO(i)=fCO(i)+  1.500;fXO2(i)=fXO2(i)+  0.500;fRO2(i)=fRO2(i)+  0.500;fHO2(i)=fHO2(i)+  1.000;

% 119, <R119>
i=i+1;
Rnames{ 119} = 'MGLY = ACO3 + HO2 + CO ';
k(:,i) = (JMGLY_IUPAC10 ); 
Gstr{i,   1}='MGLY';
fMGLY(i)=fMGLY(i)-1.0;
fACO3(i)=fACO3(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fCO(i)=fCO(i)+  1.000;

% 120, <R120>
i=i+1;
Rnames{ 120} = 'MGLY + NO3 = HNO3 + ACO3 + XO2 + RO2 ';
k(:,i) = (  1.4000E-12.*exp( -1.8600E+03./T) ); 
Gstr{i,   1}='MGLY';Gstr{i,   2}='NO3';
fMGLY(i)=fMGLY(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;fACO3(i)=fACO3(i)+  1.000;fXO2(i)=fXO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 121, <R121>
i=i+1;
Rnames{ 121} = 'MGLY +  HO = ACO3 + CO ';
k(:,i) = (  1.9000E-12.*exp(  5.7500E+02./T) ); 
Gstr{i,   1}='MGLY';Gstr{i,   2}='HO';
fMGLY(i)=fMGLY(i)-1.0;fHO(i)=fHO(i)-1.0;
fACO3(i)=fACO3(i)+  1.000;fCO(i)=fCO(i)+  1.000;

% 122, <R122>
i=i+1;
Rnames{ 122} = ' HO + H2 = HO2 ';
k(:,i) = (  7.7000E-12.*exp( -2.1000E+03./T) ).*H2; 
Gstr{i,   1}='HO';
fHO(i)=fHO(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;

% 123, <R123>
i=i+1;
Rnames{ 123} = 'CO +  HO = HO2 ';
xk0 =   1.4400E-13.*exp(  0.0000E+00./T);
xk1 =   3.4300E-33.*exp(  0.0000E+00./T);
k(:,i) = (xk0+xk1.*M ); 
Gstr{i,   1}='CO';Gstr{i,   2}='HO';
fCO(i)=fCO(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;

% 124, <R124>
i=i+1;
Rnames{ 124} = ' HO + CH4 = MO2 + RO2 ';
k(:,i) = (  1.8500E-12.*exp( -1.6900E+03./T) ).*CH4; 
Gstr{i,   1}='HO';
fHO(i)=fHO(i)-1.0;
fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 125, <R125>
i=i+1;
Rnames{ 125} = 'ETHA +  HO = 0.99100*ALD2 +  0.99100*XO2H +  0.00900*XO2N + RO2 ';
k(:,i) = (  6.9000E-12.*exp( -1.0000E+03./T) ); 
Gstr{i,   1}='ETHA';Gstr{i,   2}='HO';
fETHA(i)=fETHA(i)-1.0;fHO(i)=fHO(i)-1.0;
fALD2(i)=fALD2(i)+  0.991;fXO2H(i)=fXO2H(i)+  0.991;fXO2N(i)=fXO2N(i)+  0.009;fRO2(i)=fRO2(i)+  1.000;

% 126, <R126>
i=i+1;
Rnames{ 126} = 'ME HO +  HO = HCHO + HO2 ';
k(:,i) = (  2.8500E-12.*exp( -3.4500E+02./T) ); 
Gstr{i,   1}='MEOH';Gstr{i,   2}='HO';
fMEOH(i)=fMEOH(i)-1.0;fHO(i)=fHO(i)-1.0;
fHCHO(i)=fHCHO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

% 127, <R127>
i=i+1;
Rnames{ 127} = 'ET HO +  HO = 0.95000*ALD2 +  0.90000*HO2 +  0.10000*XO2H +  0.10000*RO2 +  0.07800*HCHO +  0.01100*GLYD ';
k(:,i) = (  3.0000E-12.*exp(  2.0000E+01./T) ); 
Gstr{i,   1}='ETOH';Gstr{i,   2}='HO';
fETOH(i)=fETOH(i)-1.0;fHO(i)=fHO(i)-1.0;
fALD2(i)=fALD2(i)+  0.950;fHO2(i)=fHO2(i)+  0.900;fXO2H(i)=fXO2H(i)+  0.100;fRO2(i)=fRO2(i)+  0.100;fHCHO(i)=fHCHO(i)+  0.078;fGLYD(i)=fGLYD(i)+  0.011;

% 128, <R128>
i=i+1;
Rnames{ 128} = 'KET = 0.50000*ALD2 +  0.50000*ACO3 +  0.50000*XO2H +  0.50000*CXO3 +  0.50000*MO2 + RO2 -  2.50000*PAR ';
k(:,i) = (JKET_IUPAC10 ); 
Gstr{i,   1}='KET';
fKET(i)=fKET(i)-1.0;
fALD2(i)=fALD2(i)+  0.500;fACO3(i)=fACO3(i)+  0.500;fXO2H(i)=fXO2H(i)+  0.500;fCXO3(i)=fCXO3(i)+  0.500;fMO2(i)=fMO2(i)+  0.500;fRO2(i)=fRO2(i)+  1.000;fPAR(i)=fPAR(i)-  2.500;

% 129, <R129>
i=i+1;
Rnames{ 129} = 'ACET = 0.38000*CO +  1.38000*MO2 +  1.38000*RO2 +  0.62000*ACO3 ';
k(:,i) = (JACET_IUPAC10 ); 
Gstr{i,   1}='ACET';
fACET(i)=fACET(i)-1.0;
fCO(i)=fCO(i)+  0.380;fMO2(i)=fMO2(i)+  1.380;fRO2(i)=fRO2(i)+  1.380;fACO3(i)=fACO3(i)+  0.620;

% 130, <R130>
i=i+1;
Rnames{ 130} = 'ACET +  HO = HCHO + ACO3 + XO2 + RO2 ';
k(:,i) = (  1.4100E-12.*exp( -6.2060E+02./T) ); 
Gstr{i,   1}='ACET';Gstr{i,   2}='HO';
fACET(i)=fACET(i)-1.0;fHO(i)=fHO(i)-1.0;
fHCHO(i)=fHCHO(i)+  1.000;fACO3(i)=fACO3(i)+  1.000;fXO2(i)=fXO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 131, <R131>
i=i+1;
Rnames{ 131} = 'PRPA +  HO = XPRP ';
k(:,i) = (  7.6000E-12.*exp( -5.8500E+02./T) ); 
Gstr{i,   1}='PRPA';Gstr{i,   2}='HO';
fPRPA(i)=fPRPA(i)-1.0;fHO(i)=fHO(i)-1.0;
fXPRP(i)=fXPRP(i)+  1.000;

% 132, <R132>
i=i+1;
Rnames{ 132} = 'PAR +  HO = XPAR ';
k(:,i) = (  8.1000E-13 ); 
Gstr{i,   1}='PAR';Gstr{i,   2}='HO';
fPAR(i)=fPAR(i)-1.0;fHO(i)=fHO(i)-1.0;
fXPAR(i)=fXPAR(i)+  1.000;

% 133, <R133>
i=i+1;
Rnames{ 133} = 'ROR = 0.20000*KET +  0.42000*ACET +  0.74000*ALD2 +  0.37000*ALDX +  0.04000*XO2N +  0.94000*XO2H +  0.98000*RO2 +  0.02000*ROR -  2.70000*PAR ';
k(:,i) = (  5.7000E+12.*exp( -5.7800E+03./T) ); 
Gstr{i,   1}='ROR';
fROR(i)=fROR(i)-1.0;
fKET(i)=fKET(i)+  0.200;fACET(i)=fACET(i)+  0.420;fALD2(i)=fALD2(i)+  0.740;fALDX(i)=fALDX(i)+  0.370;fXO2N(i)=fXO2N(i)+  0.040;fXO2H(i)=fXO2H(i)+  0.940;fRO2(i)=fRO2(i)+  0.980;fROR(i)=fROR(i)+  0.020;fPAR(i)=fPAR(i)-  2.700;

% 134, <R134>
i=i+1;
Rnames{ 134} = 'ROR + O2 = KET + HO2 ';
k(:,i) = (  1.5000E-14.*exp( -2.0000E+02./T) ).*O2; 
Gstr{i,   1}='ROR';
fROR(i)=fROR(i)-1.0;
fKET(i)=fKET(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

% 135, <R135>
i=i+1;
Rnames{ 135} = 'ROR + NO2 = NTR1 ';
k(:,i) = (  8.6000E-12.*exp(  4.0000E+02./T) ); 
Gstr{i,   1}='ROR';Gstr{i,   2}='NO2';
fROR(i)=fROR(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fNTR1(i)=fNTR1(i)+  1.000;

% 136, <R136>
i=i+1;
Rnames{ 136} = 'ETHY +  HO = 0.70000*GLY +  0.70000* HO +  0.30000*FACD +  0.30000*CO +  0.30000*HO2 ';
xko =   5.0000E-30.*M.*exp(  0.0000E+00./T).*(T./300).^ -1.5000E+00;
xkinf =   1.0000E-12.*exp(  0.0000E+00./T).*(T./300).^  0.0000E+00;
xn =   1.3000E+00;
F =   3.7000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='ETHY';Gstr{i,   2}='HO';
fETHY(i)=fETHY(i)-1.0;fHO(i)=fHO(i)-1.0;
fGLY(i)=fGLY(i)+  0.700;fHO(i)=fHO(i)+  0.700;fFACD(i)=fFACD(i)+  0.300;fCO(i)=fCO(i)+  0.300;fHO2(i)=fHO2(i)+  0.300;

% 137, <R137>
i=i+1;
Rnames{ 137} = 'ETH + O = HCHO + HO2 + CO +  0.70000*XO2H +  0.70000*RO2 +  0.30000* HO ';
k(:,i) = (  1.0400E-11.*exp( -7.9200E+02./T) ); 
Gstr{i,   1}='ETH';Gstr{i,   2}='O';
fETH(i)=fETH(i)-1.0;fO(i)=fO(i)-1.0;
fHCHO(i)=fHCHO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fCO(i)=fCO(i)+  1.000;fXO2H(i)=fXO2H(i)+  0.700;fRO2(i)=fRO2(i)+  0.700;fHO(i)=fHO(i)+  0.300;

% 138, <R138>
i=i+1;
Rnames{ 138} = 'ETH +  HO = XO2H + RO2 +  1.56000*HCHO +  0.22000*GLYD ';
xko =   8.6000E-29.*M.*exp(  0.0000E+00./T).*(T./300).^ -3.1000E+00;
xkinf =   9.0000E-12.*exp(  0.0000E+00./T).*(T./300).^ -8.5000E-01;
xn =   1.1500E+00;
F =   4.8000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='ETH';Gstr{i,   2}='HO';
fETH(i)=fETH(i)-1.0;fHO(i)=fHO(i)-1.0;
fXO2H(i)=fXO2H(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;fHCHO(i)=fHCHO(i)+  1.560;fGLYD(i)=fGLYD(i)+  0.220;

% 139, <R139>
i=i+1;
Rnames{ 139} = 'ETH + O3 = HCHO +  0.51000*CO +  0.16000*HO2 +  0.16000* HO +  0.37000*FACD ';
k(:,i) = (  9.1000E-15.*exp( -2.5800E+03./T) ); 
Gstr{i,   1}='ETH';Gstr{i,   2}='O3';
fETH(i)=fETH(i)-1.0;fO3(i)=fO3(i)-1.0;
fHCHO(i)=fHCHO(i)+  1.000;fCO(i)=fCO(i)+  0.510;fHO2(i)=fHO2(i)+  0.160;fHO(i)=fHO(i)+  0.160;fFACD(i)=fFACD(i)+  0.370;

% 140, <R140>
i=i+1;
Rnames{ 140} = 'ETH + NO3 = 0.50000*NO2 +  0.50000*NTR1 +  0.50000*XO2H +  0.50000*XO2 + RO2 +  1.12500*HCHO ';
k(:,i) = (  3.3000E-12.*exp( -2.8800E+03./T) ); 
Gstr{i,   1}='ETH';Gstr{i,   2}='NO3';
fETH(i)=fETH(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO2(i)=fNO2(i)+  0.500;fNTR1(i)=fNTR1(i)+  0.500;fXO2H(i)=fXO2H(i)+  0.500;fXO2(i)=fXO2(i)+  0.500;fRO2(i)=fRO2(i)+  1.000;fHCHO(i)=fHCHO(i)+  1.125;

% 141, <R141>
i=i+1;
Rnames{ 141} = 'OLE + O = 0.20000*ALD2 +  0.30000*ALDX +  0.10000*HO2 +  0.20000*XO2H +  0.20000*CO +  0.20000*HCHO +  0.01000*XO2N +  0.21000*RO2 +  0.20000*PAR +  0.10000* HO ';
k(:,i) = (  1.0000E-11.*exp( -2.8000E+02./T) ); 
Gstr{i,   1}='OLE';Gstr{i,   2}='O';
fOLE(i)=fOLE(i)-1.0;fO(i)=fO(i)-1.0;
fALD2(i)=fALD2(i)+  0.200;fALDX(i)=fALDX(i)+  0.300;fHO2(i)=fHO2(i)+  0.100;fXO2H(i)=fXO2H(i)+  0.200;fCO(i)=fCO(i)+  0.200;fHCHO(i)=fHCHO(i)+  0.200;fXO2N(i)=fXO2N(i)+  0.010;fRO2(i)=fRO2(i)+  0.210;fPAR(i)=fPAR(i)+  0.200;fHO(i)=fHO(i)+  0.100;

% 142, <R142>
i=i+1;
Rnames{ 142} = 'OLE +  HO = 0.78100*HCHO +  0.48800*ALD2 +  0.48800*ALDX +  0.97600*XO2H +  0.19500*XO2 +  0.02400*XO2N +  1.19500*RO2 -  0.73000*PAR ';
xko =   8.0000E-27.*M.*exp(  0.0000E+00./T).*(T./300).^ -3.5000E+00;
xkinf =   3.0000E-11.*exp(  0.0000E+00./T).*(T./300).^ -1.0000E+00;
xn =   1.1300E+00;
F =   5.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='OLE';Gstr{i,   2}='HO';
fOLE(i)=fOLE(i)-1.0;fHO(i)=fHO(i)-1.0;
fHCHO(i)=fHCHO(i)+  0.781;fALD2(i)=fALD2(i)+  0.488;fALDX(i)=fALDX(i)+  0.488;fXO2H(i)=fXO2H(i)+  0.976;fXO2(i)=fXO2(i)+  0.195;fXO2N(i)=fXO2N(i)+  0.024;fRO2(i)=fRO2(i)+  1.195;fPAR(i)=fPAR(i)-  0.730;

% 143, <R143>
i=i+1;
Rnames{ 143} = 'OLE + O3 = 0.29500*ALD2 +  0.55500*HCHO +  0.27000*ALDX +  0.15000*XO2H +  0.15000*RO2 +  0.33400* HO +  0.08000*HO2 +  0.37800*CO +  0.07500*GLY +  0.07500*MGLY +  0.09000*FACD +  0.13000*AACD +  0.04000*H2O2 -  0.79000*PAR ';
k(:,i) = (  5.5000E-15.*exp( -1.8800E+03./T) ); 
Gstr{i,   1}='OLE';Gstr{i,   2}='O3';
fOLE(i)=fOLE(i)-1.0;fO3(i)=fO3(i)-1.0;
fALD2(i)=fALD2(i)+  0.295;fHCHO(i)=fHCHO(i)+  0.555;fALDX(i)=fALDX(i)+  0.270;fXO2H(i)=fXO2H(i)+  0.150;fRO2(i)=fRO2(i)+  0.150;fHO(i)=fHO(i)+  0.334;fHO2(i)=fHO2(i)+  0.080;fCO(i)=fCO(i)+  0.378;fGLY(i)=fGLY(i)+  0.075;fMGLY(i)=fMGLY(i)+  0.075;fFACD(i)=fFACD(i)+  0.090;fAACD(i)=fAACD(i)+  0.130;fH2O2(i)=fH2O2(i)+  0.040;fPAR(i)=fPAR(i)-  0.790;

% 144, <R144>
i=i+1;
Rnames{ 144} = 'OLE + NO3 = 0.50000*NO2 +  0.50000*NTR1 +  0.48000*XO2 +  0.48000*XO2H +  0.04000*XO2N + RO2 +  0.50000*HCHO +  0.25000*ALD2 +  0.37500*ALDX - PAR ';
k(:,i) = (  4.6000E-13.*exp( -1.1550E+03./T) ); 
Gstr{i,   1}='OLE';Gstr{i,   2}='NO3';
fOLE(i)=fOLE(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO2(i)=fNO2(i)+  0.500;fNTR1(i)=fNTR1(i)+  0.500;fXO2(i)=fXO2(i)+  0.480;fXO2H(i)=fXO2H(i)+  0.480;fXO2N(i)=fXO2N(i)+  0.040;fRO2(i)=fRO2(i)+  1.000;fHCHO(i)=fHCHO(i)+  0.500;fALD2(i)=fALD2(i)+  0.250;fALDX(i)=fALDX(i)+  0.375;fPAR(i)=fPAR(i)-  1.000;

% 145, <R145>
i=i+1;
Rnames{ 145} = 'IOLE + O = 1.24000*ALD2 +  0.66000*ALDX +  0.10000*XO2H +  0.10000*RO2 +  0.10000*CO +  0.10000*PAR ';
k(:,i) = (  2.3000E-11 ); 
Gstr{i,   1}='IOLE';Gstr{i,   2}='O';
fIOLE(i)=fIOLE(i)-1.0;fO(i)=fO(i)-1.0;
fALD2(i)=fALD2(i)+  1.240;fALDX(i)=fALDX(i)+  0.660;fXO2H(i)=fXO2H(i)+  0.100;fRO2(i)=fRO2(i)+  0.100;fCO(i)=fCO(i)+  0.100;fPAR(i)=fPAR(i)+  0.100;

% 146, <R146>
i=i+1;
Rnames{ 146} = 'IOLE +  HO = 1.30000*ALD2 +  0.70000*ALDX + XO2H + RO2 ';
k(:,i) = (  1.0500E-11.*exp(  5.1900E+02./T) ); 
Gstr{i,   1}='IOLE';Gstr{i,   2}='HO';
fIOLE(i)=fIOLE(i)-1.0;fHO(i)=fHO(i)-1.0;
fALD2(i)=fALD2(i)+  1.300;fALDX(i)=fALDX(i)+  0.700;fXO2H(i)=fXO2H(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 147, <R147>
i=i+1;
Rnames{ 147} = 'IOLE + O3 = 0.73200*ALD2 +  0.44200*ALDX +  0.12800*HCHO +  0.24500*CO +  0.50000* HO +  0.30000*XO2H +  0.30000*RO2 +  0.24000*GLY +  0.06000*MGLY +  0.29000*PAR +  0.08000*AACD +  0.08000*H2O2 ';
k(:,i) = (  4.7000E-15.*exp( -1.0130E+03./T) ); 
Gstr{i,   1}='IOLE';Gstr{i,   2}='O3';
fIOLE(i)=fIOLE(i)-1.0;fO3(i)=fO3(i)-1.0;
fALD2(i)=fALD2(i)+  0.732;fALDX(i)=fALDX(i)+  0.442;fHCHO(i)=fHCHO(i)+  0.128;fCO(i)=fCO(i)+  0.245;fHO(i)=fHO(i)+  0.500;fXO2H(i)=fXO2H(i)+  0.300;fRO2(i)=fRO2(i)+  0.300;fGLY(i)=fGLY(i)+  0.240;fMGLY(i)=fMGLY(i)+  0.060;fPAR(i)=fPAR(i)+  0.290;fAACD(i)=fAACD(i)+  0.080;fH2O2(i)=fH2O2(i)+  0.080;

% 148, <R148>
i=i+1;
Rnames{ 148} = 'IOLE + NO3 = 0.50000*NO2 +  0.50000*NTR1 +  0.48000*XO2 +  0.48000*XO2H +  0.04000*XO2N + RO2 +  0.50000*ALD2 +  0.62500*ALDX + PAR ';
k(:,i) = (  3.7000E-13 ); 
Gstr{i,   1}='IOLE';Gstr{i,   2}='NO3';
fIOLE(i)=fIOLE(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO2(i)=fNO2(i)+  0.500;fNTR1(i)=fNTR1(i)+  0.500;fXO2(i)=fXO2(i)+  0.480;fXO2H(i)=fXO2H(i)+  0.480;fXO2N(i)=fXO2N(i)+  0.040;fRO2(i)=fRO2(i)+  1.000;fALD2(i)=fALD2(i)+  0.500;fALDX(i)=fALDX(i)+  0.625;fPAR(i)=fPAR(i)+  1.000;

% 149, <R149>
i=i+1;
Rnames{ 149} = 'ISOP +  HO = ISO2 + RO2 + ISOPRXN ';
k(:,i) = (  2.7000E-11.*exp(  3.9000E+02./T) ); 
Gstr{i,   1}='ISOP';Gstr{i,   2}='HO';
fISOP(i)=fISOP(i)-1.0;fHO(i)=fHO(i)-1.0;
fISO2(i)=fISO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;fISOPRXN(i)=fISOPRXN(i)+  1.000;

% 150, <R150>
i=i+1;
Rnames{ 150} = 'ISOP + O = 0.75000*ISPD +  0.50000*HCHO +  0.25000*XO2 +  0.25000*RO2 +  0.25000*HO2 +  0.25000*CXO3 +  0.25000*PAR ';
k(:,i) = (  3.0000E-11 ); 
Gstr{i,   1}='ISOP';Gstr{i,   2}='O';
fISOP(i)=fISOP(i)-1.0;fO(i)=fO(i)-1.0;
fISPD(i)=fISPD(i)+  0.750;fHCHO(i)=fHCHO(i)+  0.500;fXO2(i)=fXO2(i)+  0.250;fRO2(i)=fRO2(i)+  0.250;fHO2(i)=fHO2(i)+  0.250;fCXO3(i)=fCXO3(i)+  0.250;fPAR(i)=fPAR(i)+  0.250;

% 151, <R151>
i=i+1;
Rnames{ 151} = 'ISO2 + NO = 0.10000*INTR +  0.90000*NO2 +  0.67300*HCHO +  0.90000*ISPD +  0.81800*HO2 +  0.08200*XO2H +  0.08200*RO2 ';
k(:,i) = (  2.3900E-12.*exp(  3.6500E+02./T) ); 
Gstr{i,   1}='ISO2';Gstr{i,   2}='NO';
fISO2(i)=fISO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fINTR(i)=fINTR(i)+  0.100;fNO2(i)=fNO2(i)+  0.900;fHCHO(i)=fHCHO(i)+  0.673;fISPD(i)=fISPD(i)+  0.900;fHO2(i)=fHO2(i)+  0.818;fXO2H(i)=fXO2H(i)+  0.082;fRO2(i)=fRO2(i)+  0.082;

% 152, <R152>
i=i+1;
Rnames{ 152} = 'ISO2 + HO2 = 0.88000*ISPX +  0.12000* HO +  0.12000*HO2 +  0.12000*HCHO +  0.12000*ISPD ';
k(:,i) = (  7.4300E-13.*exp(  7.0000E+02./T) ); 
Gstr{i,   1}='ISO2';Gstr{i,   2}='HO2';
fISO2(i)=fISO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fISPX(i)=fISPX(i)+  0.880;fHO(i)=fHO(i)+  0.120;fHO2(i)=fHO2(i)+  0.120;fHCHO(i)=fHCHO(i)+  0.120;fISPD(i)=fISPD(i)+  0.120;

% 153, <R153>
i=i+1;
Rnames{ 153} = 'ISO2 + ACO3 = 0.59800*HCHO + ISPD +  0.72800*HO2 +  0.07200*XO2H +  0.80000*MO2 +  0.20000*AACD +  0.87200*RO2 ';
k(:,i) = (k(:,  58) ); 
Gstr{i,   1}='ISO2';Gstr{i,   2}='ACO3';
fISO2(i)=fISO2(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fHCHO(i)=fHCHO(i)+  0.598;fISPD(i)=fISPD(i)+  1.000;fHO2(i)=fHO2(i)+  0.728;fXO2H(i)=fXO2H(i)+  0.072;fMO2(i)=fMO2(i)+  0.800;fAACD(i)=fAACD(i)+  0.200;fRO2(i)=fRO2(i)+  0.872;

% 154, <R154>
i=i+1;
Rnames{ 154} = 'ISO2 + RO2 = 0.59800*HCHO + ISPD +  0.72800*HO2 +  0.07200*XO2H +  1.07200*RO2 ';
k(:,i) = (k(:,  70) ); 
Gstr{i,   1}='ISO2';Gstr{i,   2}='RO2';
fISO2(i)=fISO2(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fHCHO(i)=fHCHO(i)+  0.598;fISPD(i)=fISPD(i)+  1.000;fHO2(i)=fHO2(i)+  0.728;fXO2H(i)=fXO2H(i)+  0.072;fRO2(i)=fRO2(i)+  1.072;

% 155, <R155>
i=i+1;
Rnames{ 155} = 'ISO2 = HO2 + HPLD ';
k(:,i) = (  3.3000E+09.*exp( -8.3000E+03./T) ); 
Gstr{i,   1}='ISO2';
fISO2(i)=fISO2(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fHPLD(i)=fHPLD(i)+  1.000;

% 156, <R156>
i=i+1;
Rnames{ 156} = 'ISOP + O3 = 0.60000*HCHO +  0.65000*ISPD +  0.15000*ALDX +  0.20000*CXO3 +  0.35000*PAR +  0.26600* HO +  0.20000*XO2 +  0.20000*RO2 +  0.06600*HO2 +  0.06600*CO ';
k(:,i) = (  1.0300E-14.*exp( -1.9950E+03./T) ); 
Gstr{i,   1}='ISOP';Gstr{i,   2}='O3';
fISOP(i)=fISOP(i)-1.0;fO3(i)=fO3(i)-1.0;
fHCHO(i)=fHCHO(i)+  0.600;fISPD(i)=fISPD(i)+  0.650;fALDX(i)=fALDX(i)+  0.150;fCXO3(i)=fCXO3(i)+  0.200;fPAR(i)=fPAR(i)+  0.350;fHO(i)=fHO(i)+  0.266;fXO2(i)=fXO2(i)+  0.200;fRO2(i)=fRO2(i)+  0.200;fHO2(i)=fHO2(i)+  0.066;fCO(i)=fCO(i)+  0.066;

% 157, <R157>
i=i+1;
Rnames{ 157} = 'ISOP + NO3 = 0.35000*NO2 +  0.65000*NTR2 +  0.64000*XO2H +  0.33000*XO2 +  0.03000*XO2N + RO2 +  0.35000*HCHO +  0.35000*ISPD + ISOPRXN ';
k(:,i) = (  3.0300E-12.*exp( -4.4800E+02./T) ); 
Gstr{i,   1}='ISOP';Gstr{i,   2}='NO3';
fISOP(i)=fISOP(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO2(i)=fNO2(i)+  0.350;fNTR2(i)=fNTR2(i)+  0.650;fXO2H(i)=fXO2H(i)+  0.640;fXO2(i)=fXO2(i)+  0.330;fXO2N(i)=fXO2N(i)+  0.030;fRO2(i)=fRO2(i)+  1.000;fHCHO(i)=fHCHO(i)+  0.350;fISPD(i)=fISPD(i)+  0.350;fISOPRXN(i)=fISOPRXN(i)+  1.000;

% 158, <R158>
i=i+1;
Rnames{ 158} = 'ISPD +  HO = 0.02200*XO2N +  0.52100*XO2 +  0.11500*MGLY +  0.11500*MO2 +  0.26900*GLYD +  0.26900*ACO3 +  0.45700*OPO3 +  0.11700*PAR +  0.13700*ACET +  0.13700*CO +  0.13700*HO2 +  0.65800*RO2 ';
k(:,i) = (  5.5800E-12.*exp(  5.1100E+02./T) ); 
Gstr{i,   1}='ISPD';Gstr{i,   2}='HO';
fISPD(i)=fISPD(i)-1.0;fHO(i)=fHO(i)-1.0;
fXO2N(i)=fXO2N(i)+  0.022;fXO2(i)=fXO2(i)+  0.521;fMGLY(i)=fMGLY(i)+  0.115;fMO2(i)=fMO2(i)+  0.115;fGLYD(i)=fGLYD(i)+  0.269;fACO3(i)=fACO3(i)+  0.269;fOPO3(i)=fOPO3(i)+  0.457;fPAR(i)=fPAR(i)+  0.117;fACET(i)=fACET(i)+  0.137;fCO(i)=fCO(i)+  0.137;fHO2(i)=fHO2(i)+  0.137;fRO2(i)=fRO2(i)+  0.658;

% 159, <R159>
i=i+1;
Rnames{ 159} = 'ISPD + O3 = 0.04000*ALD2 +  0.23100*HCHO +  0.53100*MGLY +  0.17000*GLY +  0.17000*ACET +  0.54300*CO +  0.46100* HO +  0.15000*FACD +  0.39800*HO2 +  0.14300*ACO3 ';
k(:,i) = (  3.8800E-15.*exp( -1.7700E+03./T) ); 
Gstr{i,   1}='ISPD';Gstr{i,   2}='O3';
fISPD(i)=fISPD(i)-1.0;fO3(i)=fO3(i)-1.0;
fALD2(i)=fALD2(i)+  0.040;fHCHO(i)=fHCHO(i)+  0.231;fMGLY(i)=fMGLY(i)+  0.531;fGLY(i)=fGLY(i)+  0.170;fACET(i)=fACET(i)+  0.170;fCO(i)=fCO(i)+  0.543;fHO(i)=fHO(i)+  0.461;fFACD(i)=fFACD(i)+  0.150;fHO2(i)=fHO2(i)+  0.398;fACO3(i)=fACO3(i)+  0.143;

% 160, <R160>
i=i+1;
Rnames{ 160} = 'ISPD + NO3 = 0.71700*HNO3 +  0.14200*NTR2 +  0.14200*NO2 +  0.14200*XO2 +  0.14200*XO2H +  0.11300*GLYD +  0.11300*MGLY +  0.71700*PAR +  0.71700*CXO3 +  0.28400*RO2 ';
k(:,i) = (  4.1000E-12.*exp( -1.8600E+03./T) ); 
Gstr{i,   1}='ISPD';Gstr{i,   2}='NO3';
fISPD(i)=fISPD(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fHNO3(i)=fHNO3(i)+  0.717;fNTR2(i)=fNTR2(i)+  0.142;fNO2(i)=fNO2(i)+  0.142;fXO2(i)=fXO2(i)+  0.142;fXO2H(i)=fXO2H(i)+  0.142;fGLYD(i)=fGLYD(i)+  0.113;fMGLY(i)=fMGLY(i)+  0.113;fPAR(i)=fPAR(i)+  0.717;fCXO3(i)=fCXO3(i)+  0.717;fRO2(i)=fRO2(i)+  0.284;

% 161, <R161>
i=i+1;
Rnames{ 161} = 'ISPD = 0.76000*HO2 +  0.34000*XO2H +  0.16000*XO2 +  0.34000*MO2 +  0.20800*ACO3 +  0.26000*HCHO +  0.24000*OLE +  0.24000*PAR +  0.17000*ACET +  0.12800*GLYD +  0.84000*RO2 ';
k(:,i) = (JISPD ); 
Gstr{i,   1}='ISPD';
fISPD(i)=fISPD(i)-1.0;
fHO2(i)=fHO2(i)+  0.760;fXO2H(i)=fXO2H(i)+  0.340;fXO2(i)=fXO2(i)+  0.160;fMO2(i)=fMO2(i)+  0.340;fACO3(i)=fACO3(i)+  0.208;fHCHO(i)=fHCHO(i)+  0.260;fOLE(i)=fOLE(i)+  0.240;fPAR(i)=fPAR(i)+  0.240;fACET(i)=fACET(i)+  0.170;fGLYD(i)=fGLYD(i)+  0.128;fRO2(i)=fRO2(i)+  0.840;

% 162, <R162>
i=i+1;
Rnames{ 162} = 'ISPX +  HO = 0.90400*EPOX +  0.93300* HO +  0.06700*ISO2 +  0.06700*RO2 +  0.02900*IOLE +  0.02900*ALDX ';
k(:,i) = (  2.2300E-11.*exp(  3.7200E+02./T) ); 
Gstr{i,   1}='ISPX';Gstr{i,   2}='HO';
fISPX(i)=fISPX(i)-1.0;fHO(i)=fHO(i)-1.0;
fEPOX(i)=fEPOX(i)+  0.904;fHO(i)=fHO(i)+  0.933;fISO2(i)=fISO2(i)+  0.067;fRO2(i)=fRO2(i)+  0.067;fIOLE(i)=fIOLE(i)+  0.029;fALDX(i)=fALDX(i)+  0.029;

% 163, <R163>
i=i+1;
Rnames{ 163} = 'HPLD =  HO + ISPD ';
k(:,i) = (JHPALD ); 
Gstr{i,   1}='HPLD';
fHPLD(i)=fHPLD(i)-1.0;
fHO(i)=fHO(i)+  1.000;fISPD(i)=fISPD(i)+  1.000;

% 164, <R164>
i=i+1;
Rnames{ 164} = 'HPLD + NO3 = HNO3 + ISPD ';
k(:,i) = (  6.0000E-12.*exp( -1.8600E+03./T) ); 
Gstr{i,   1}='HPLD';Gstr{i,   2}='NO3';
fHPLD(i)=fHPLD(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;fISPD(i)=fISPD(i)+  1.000;

% 165, <R165>
i=i+1;
Rnames{ 165} = 'EPOX +  HO = EPX2 + RO2 ';
k(:,i) = (  5.7800E-11.*exp( -4.0000E+02./T) ); 
Gstr{i,   1}='EPOX';Gstr{i,   2}='HO';
fEPOX(i)=fEPOX(i)-1.0;fHO(i)=fHO(i)-1.0;
fEPX2(i)=fEPX2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 166, <R166>
i=i+1;
Rnames{ 166} = 'EPX2 + HO2 = 0.27500*GLYD +  0.27500*GLY +  0.27500*MGLY +  1.12500* HO +  0.82500*HO2 +  0.37500*HCHO +  0.07400*FACD +  0.25100*CO +  2.17500*PAR ';
k(:,i) = (  7.4300E-13.*exp(  7.0000E+02./T) ); 
Gstr{i,   1}='EPX2';Gstr{i,   2}='HO2';
fEPX2(i)=fEPX2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fGLYD(i)=fGLYD(i)+  0.275;fGLY(i)=fGLY(i)+  0.275;fMGLY(i)=fMGLY(i)+  0.275;fHO(i)=fHO(i)+  1.125;fHO2(i)=fHO2(i)+  0.825;fHCHO(i)=fHCHO(i)+  0.375;fFACD(i)=fFACD(i)+  0.074;fCO(i)=fCO(i)+  0.251;fPAR(i)=fPAR(i)+  2.175;

% 167, <R167>
i=i+1;
Rnames{ 167} = 'EPX2 + NO = 0.27500*GLYD +  0.27500*GLY +  0.27500*MGLY +  0.12500* HO +  0.82500*HO2 +  0.37500*HCHO + NO2 +  0.25100*CO +  2.17500*PAR ';
k(:,i) = (  2.3900E-12.*exp(  3.6500E+02./T) ); 
Gstr{i,   1}='EPX2';Gstr{i,   2}='NO';
fEPX2(i)=fEPX2(i)-1.0;fNO(i)=fNO(i)-1.0;
fGLYD(i)=fGLYD(i)+  0.275;fGLY(i)=fGLY(i)+  0.275;fMGLY(i)=fMGLY(i)+  0.275;fHO(i)=fHO(i)+  0.125;fHO2(i)=fHO2(i)+  0.825;fHCHO(i)=fHCHO(i)+  0.375;fNO2(i)=fNO2(i)+  1.000;fCO(i)=fCO(i)+  0.251;fPAR(i)=fPAR(i)+  2.175;

% 168, <R168>
i=i+1;
Rnames{ 168} = 'EPX2 + ACO3 = 0.22000*GLYD +  0.22000*GLY +  0.22000*MGLY +  0.10000* HO +  0.66000*HO2 +  0.30000*HCHO +  0.20000*CO +  1.74000*PAR +  0.80000*MO2 +  0.20000*AACD +  0.80000*RO2 ';
k(:,i) = (k(:,  58) ); 
Gstr{i,   1}='EPX2';Gstr{i,   2}='ACO3';
fEPX2(i)=fEPX2(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fGLYD(i)=fGLYD(i)+  0.220;fGLY(i)=fGLY(i)+  0.220;fMGLY(i)=fMGLY(i)+  0.220;fHO(i)=fHO(i)+  0.100;fHO2(i)=fHO2(i)+  0.660;fHCHO(i)=fHCHO(i)+  0.300;fCO(i)=fCO(i)+  0.200;fPAR(i)=fPAR(i)+  1.740;fMO2(i)=fMO2(i)+  0.800;fAACD(i)=fAACD(i)+  0.200;fRO2(i)=fRO2(i)+  0.800;

% 169, <R169>
i=i+1;
Rnames{ 169} = 'EPX2 + RO2 = 0.27500*GLYD +  0.27500*GLY +  0.27500*MGLY +  0.12500* HO +  0.82500*HO2 +  0.37500*HCHO +  0.25100*CO +  2.17500*PAR + RO2 ';
k(:,i) = (k(:,  70) ); 
Gstr{i,   1}='EPX2';Gstr{i,   2}='RO2';
fEPX2(i)=fEPX2(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fGLYD(i)=fGLYD(i)+  0.275;fGLY(i)=fGLY(i)+  0.275;fMGLY(i)=fMGLY(i)+  0.275;fHO(i)=fHO(i)+  0.125;fHO2(i)=fHO2(i)+  0.825;fHCHO(i)=fHCHO(i)+  0.375;fCO(i)=fCO(i)+  0.251;fPAR(i)=fPAR(i)+  2.175;fRO2(i)=fRO2(i)+  1.000;

% 170, <R170>
i=i+1;
Rnames{ 170} = 'INTR +  HO = 0.63000*XO2 +  0.37000*XO2H + RO2 +  0.44400*NO2 +  0.18500*NO3 +  0.10400*INTR +  0.59200*HCHO +  0.33100*GLYD +  0.18500*FACD +  2.70000*PAR +  0.09800*OLE +  0.07800*ALDX +  0.26600*NTR2 ';
k(:,i) = (  3.1000E-11 ); 
Gstr{i,   1}='INTR';Gstr{i,   2}='HO';
fINTR(i)=fINTR(i)-1.0;fHO(i)=fHO(i)-1.0;
fXO2(i)=fXO2(i)+  0.630;fXO2H(i)=fXO2H(i)+  0.370;fRO2(i)=fRO2(i)+  1.000;fNO2(i)=fNO2(i)+  0.444;fNO3(i)=fNO3(i)+  0.185;fINTR(i)=fINTR(i)+  0.104;fHCHO(i)=fHCHO(i)+  0.592;fGLYD(i)=fGLYD(i)+  0.331;fFACD(i)=fFACD(i)+  0.185;fPAR(i)=fPAR(i)+  2.700;fOLE(i)=fOLE(i)+  0.098;fALDX(i)=fALDX(i)+  0.078;fNTR2(i)=fNTR2(i)+  0.266;

% 171, <R171>
i=i+1;
Rnames{ 171} = 'TERP + O = 0.15000*ALDX +  5.12000*PAR + TRPRXN ';
k(:,i) = (  3.6000E-11 ); 
Gstr{i,   1}='TERP';Gstr{i,   2}='O';
fTERP(i)=fTERP(i)-1.0;fO(i)=fO(i)-1.0;
fALDX(i)=fALDX(i)+  0.150;fPAR(i)=fPAR(i)+  5.120;fTRPRXN(i)=fTRPRXN(i)+  1.000;

% 172, <R172>
i=i+1;
Rnames{ 172} = 'TERP +  HO = 0.75000*XO2H +  0.50000*XO2 +  0.25000*XO2N +  1.50000*RO2 +  0.28000*HCHO +  1.66000*PAR +  0.47000*ALDX + TRPRXN ';
k(:,i) = (  1.5000E-11.*exp(  4.4900E+02./T) ); 
Gstr{i,   1}='TERP';Gstr{i,   2}='HO';
fTERP(i)=fTERP(i)-1.0;fHO(i)=fHO(i)-1.0;
fXO2H(i)=fXO2H(i)+  0.750;fXO2(i)=fXO2(i)+  0.500;fXO2N(i)=fXO2N(i)+  0.250;fRO2(i)=fRO2(i)+  1.500;fHCHO(i)=fHCHO(i)+  0.280;fPAR(i)=fPAR(i)+  1.660;fALDX(i)=fALDX(i)+  0.470;fTRPRXN(i)=fTRPRXN(i)+  1.000;

% 173, <R173>
i=i+1;
Rnames{ 173} = 'TERP + O3 = 0.57000* HO +  0.07000*XO2H +  0.69000*XO2 +  0.18000*XO2N +  0.94000*RO2 +  0.24000*HCHO +  0.00100*CO +  7.00000*PAR +  0.21000*ALDX +  0.39000*CXO3 + TRPRXN ';
k(:,i) = (  1.2000E-15.*exp( -8.2100E+02./T) ); 
Gstr{i,   1}='TERP';Gstr{i,   2}='O3';
fTERP(i)=fTERP(i)-1.0;fO3(i)=fO3(i)-1.0;
fHO(i)=fHO(i)+  0.570;fXO2H(i)=fXO2H(i)+  0.070;fXO2(i)=fXO2(i)+  0.690;fXO2N(i)=fXO2N(i)+  0.180;fRO2(i)=fRO2(i)+  0.940;fHCHO(i)=fHCHO(i)+  0.240;fCO(i)=fCO(i)+  0.001;fPAR(i)=fPAR(i)+  7.000;fALDX(i)=fALDX(i)+  0.210;fCXO3(i)=fCXO3(i)+  0.390;fTRPRXN(i)=fTRPRXN(i)+  1.000;

% 174, <R174>
i=i+1;
Rnames{ 174} = 'TERP + NO3 = 0.47000*NO2 +  0.28000*XO2H +  0.75000*XO2 +  0.25000*XO2N +  1.28000*RO2 +  0.47000*ALDX +  0.53000*NTR2 + TERPNRO2 ';
k(:,i) = (  3.7000E-12.*exp(  1.7500E+02./T) ); 
Gstr{i,   1}='TERP';Gstr{i,   2}='NO3';
fTERP(i)=fTERP(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO2(i)=fNO2(i)+  0.470;fXO2H(i)=fXO2H(i)+  0.280;fXO2(i)=fXO2(i)+  0.750;fXO2N(i)=fXO2N(i)+  0.250;fRO2(i)=fRO2(i)+  1.280;fALDX(i)=fALDX(i)+  0.470;fNTR2(i)=fNTR2(i)+  0.530;fTERPNRO2(i)=fTERPNRO2(i)+  1.000;

% 175, <R171a>
i=i+1;
Rnames{ 175} = 'APIN + O = 0.15000*ALDX +  5.12000*PAR + TRPRXN ';
k(:,i) = (  3.6000E-11 ); 
Gstr{i,   1}='APIN';Gstr{i,   2}='O';
fAPIN(i)=fAPIN(i)-1.0;fO(i)=fO(i)-1.0;
fALDX(i)=fALDX(i)+  0.150;fPAR(i)=fPAR(i)+  5.120;fTRPRXN(i)=fTRPRXN(i)+  1.000;

% 176, <R172a>
i=i+1;
Rnames{ 176} = 'APIN +  HO = 0.75000*XO2H +  0.50000*XO2 +  0.25000*XO2N +  1.50000*RO2 +  0.28000*HCHO +  1.66000*PAR +  0.47000*ALDX + TRPRXN ';
k(:,i) = (  1.5000E-11.*exp(  4.4900E+02./T) ); 
Gstr{i,   1}='APIN';Gstr{i,   2}='HO';
fAPIN(i)=fAPIN(i)-1.0;fHO(i)=fHO(i)-1.0;
fXO2H(i)=fXO2H(i)+  0.750;fXO2(i)=fXO2(i)+  0.500;fXO2N(i)=fXO2N(i)+  0.250;fRO2(i)=fRO2(i)+  1.500;fHCHO(i)=fHCHO(i)+  0.280;fPAR(i)=fPAR(i)+  1.660;fALDX(i)=fALDX(i)+  0.470;fTRPRXN(i)=fTRPRXN(i)+  1.000;

% 177, <R173a>
i=i+1;
Rnames{ 177} = 'APIN + O3 = 0.57000* HO +  0.07000*XO2H +  0.69000*XO2 +  0.18000*XO2N +  0.94000*RO2 +  0.24000*HCHO +  0.00100*CO +  7.00000*PAR +  0.21000*ALDX +  0.39000*CXO3 + TRPRXN ';
k(:,i) = (  1.2000E-15.*exp( -8.2100E+02./T) ); 
Gstr{i,   1}='APIN';Gstr{i,   2}='O3';
fAPIN(i)=fAPIN(i)-1.0;fO3(i)=fO3(i)-1.0;
fHO(i)=fHO(i)+  0.570;fXO2H(i)=fXO2H(i)+  0.070;fXO2(i)=fXO2(i)+  0.690;fXO2N(i)=fXO2N(i)+  0.180;fRO2(i)=fRO2(i)+  0.940;fHCHO(i)=fHCHO(i)+  0.240;fCO(i)=fCO(i)+  0.001;fPAR(i)=fPAR(i)+  7.000;fALDX(i)=fALDX(i)+  0.210;fCXO3(i)=fCXO3(i)+  0.390;fTRPRXN(i)=fTRPRXN(i)+  1.000;

% 178, <R174a>
i=i+1;
Rnames{ 178} = 'APIN + NO3 = 0.47000*NO2 +  0.28000*XO2H +  0.75000*XO2 +  0.25000*XO2N +  1.28000*RO2 +  0.47000*ALDX +  0.53000*NTR2 ';
k(:,i) = (  3.7000E-12.*exp(  1.7500E+02./T) ); 
Gstr{i,   1}='APIN';Gstr{i,   2}='NO3';
fAPIN(i)=fAPIN(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO2(i)=fNO2(i)+  0.470;fXO2H(i)=fXO2H(i)+  0.280;fXO2(i)=fXO2(i)+  0.750;fXO2N(i)=fXO2N(i)+  0.250;fRO2(i)=fRO2(i)+  1.280;fALDX(i)=fALDX(i)+  0.470;fNTR2(i)=fNTR2(i)+  0.530;

% 179, <R175>
i=i+1;
Rnames{ 179} = 'BENZENE +  HO = 0.53000*CRES +  0.35200*BZO2 +  0.35200*RO2 +  0.11800*OPEN +  0.11800* HO +  0.53000*HO2 + BENZRO2 ';
k(:,i) = (  2.3000E-12.*exp( -1.9000E+02./T) ); 
Gstr{i,   1}='BENZENE';Gstr{i,   2}='HO';
fBENZENE(i)=fBENZENE(i)-1.0;fHO(i)=fHO(i)-1.0;
fCRES(i)=fCRES(i)+  0.530;fBZO2(i)=fBZO2(i)+  0.352;fRO2(i)=fRO2(i)+  0.352;fOPEN(i)=fOPEN(i)+  0.118;fHO(i)=fHO(i)+  0.118;fHO2(i)=fHO2(i)+  0.530;fBENZRO2(i)=fBENZRO2(i)+  1.000;

% 180, <R176>
i=i+1;
Rnames{ 180} = 'BZO2 + NO = 0.91800*NO2 +  0.08200*NTR2 +  0.91800*GLY +  0.91800*OPEN +  0.91800*HO2 ';
k(:,i) = (  2.7000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='BZO2';Gstr{i,   2}='NO';
fBZO2(i)=fBZO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  0.918;fNTR2(i)=fNTR2(i)+  0.082;fGLY(i)=fGLY(i)+  0.918;fOPEN(i)=fOPEN(i)+  0.918;fHO2(i)=fHO2(i)+  0.918;

% 181, <R177>
i=i+1;
Rnames{ 181} = 'BZO2 + ACO3 = GLY + OPEN + HO2 + MO2 + RO2 ';
k(:,i) = (k(:,  58) ); 
Gstr{i,   1}='BZO2';Gstr{i,   2}='ACO3';
fBZO2(i)=fBZO2(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fGLY(i)=fGLY(i)+  1.000;fOPEN(i)=fOPEN(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 182, <R178>
i=i+1;
Rnames{ 182} = 'BZO2 + HO2 =';
k(:,i) = (  1.9000E-13.*exp(  1.3000E+03./T) ); 
Gstr{i,   1}='BZO2';Gstr{i,   2}='HO2';
fBZO2(i)=fBZO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;


% 183, <R179>
i=i+1;
Rnames{ 183} = 'BZO2 + RO2 = GLY + OPEN + HO2 + RO2 ';
k(:,i) = (k(:,  70) ); 
Gstr{i,   1}='BZO2';Gstr{i,   2}='RO2';
fBZO2(i)=fBZO2(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fGLY(i)=fGLY(i)+  1.000;fOPEN(i)=fOPEN(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 184, <R180>
i=i+1;
Rnames{ 184} = 'TOL +  HO = 0.18000*CRES +  0.65000*TO2 +  0.72000*RO2 +  0.10000*OPEN +  0.10000* HO +  0.07000*XO2H +  0.18000*HO2 + TOLRO2 ';
k(:,i) = (  1.8000E-12.*exp(  3.4000E+02./T) ); 
Gstr{i,   1}='TOL';Gstr{i,   2}='HO';
fTOL(i)=fTOL(i)-1.0;fHO(i)=fHO(i)-1.0;
fCRES(i)=fCRES(i)+  0.180;fTO2(i)=fTO2(i)+  0.650;fRO2(i)=fRO2(i)+  0.720;fOPEN(i)=fOPEN(i)+  0.100;fHO(i)=fHO(i)+  0.100;fXO2H(i)=fXO2H(i)+  0.070;fHO2(i)=fHO2(i)+  0.180;fTOLRO2(i)=fTOLRO2(i)+  1.000;

% 185, <R181>
i=i+1;
Rnames{ 185} = 'TO2 + NO = 0.86000*NO2 +  0.14000*NTR2 +  0.41700*GLY +  0.44300*MGLY +  0.66000*OPEN +  0.20000*XOPN +  0.86000*HO2 ';
k(:,i) = (  2.7000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='TO2';Gstr{i,   2}='NO';
fTO2(i)=fTO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  0.860;fNTR2(i)=fNTR2(i)+  0.140;fGLY(i)=fGLY(i)+  0.417;fMGLY(i)=fMGLY(i)+  0.443;fOPEN(i)=fOPEN(i)+  0.660;fXOPN(i)=fXOPN(i)+  0.200;fHO2(i)=fHO2(i)+  0.860;

% 186, <R182>
i=i+1;
Rnames{ 186} = 'TO2 + ACO3 = 0.48000*GLY +  0.52000*MGLY +  0.77000*OPEN +  0.23000*XOPN + HO2 + MO2 + RO2 ';
k(:,i) = (k(:,  58) ); 
Gstr{i,   1}='TO2';Gstr{i,   2}='ACO3';
fTO2(i)=fTO2(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fGLY(i)=fGLY(i)+  0.480;fMGLY(i)=fMGLY(i)+  0.520;fOPEN(i)=fOPEN(i)+  0.770;fXOPN(i)=fXOPN(i)+  0.230;fHO2(i)=fHO2(i)+  1.000;fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 187, <R183>
i=i+1;
Rnames{ 187} = 'TO2 + HO2 =';
k(:,i) = (  1.9000E-13.*exp(  1.3000E+03./T) ); 
Gstr{i,   1}='TO2';Gstr{i,   2}='HO2';
fTO2(i)=fTO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;


% 188, <R184>
i=i+1;
Rnames{ 188} = 'TO2 + RO2 = 0.48000*GLY +  0.52000*MGLY +  0.77000*OPEN +  0.23000*XOPN + HO2 + RO2 ';
k(:,i) = (k(:,  70) ); 
Gstr{i,   1}='TO2';Gstr{i,   2}='RO2';
fTO2(i)=fTO2(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fGLY(i)=fGLY(i)+  0.480;fMGLY(i)=fMGLY(i)+  0.520;fOPEN(i)=fOPEN(i)+  0.770;fXOPN(i)=fXOPN(i)+  0.230;fHO2(i)=fHO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 189, <R185>
i=i+1;
Rnames{ 189} = 'XYLMN +  HO = 0.15500*CRES +  0.54400*XLO2 +  0.60200*RO2 +  0.24400*XOPN +  0.24400* HO +  0.05800*XO2H +  0.15500*HO2 + XYLRO2 ';
k(:,i) = (  1.8500E-11 ); 
Gstr{i,   1}='XYLMN';Gstr{i,   2}='HO';
fXYLMN(i)=fXYLMN(i)-1.0;fHO(i)=fHO(i)-1.0;
fCRES(i)=fCRES(i)+  0.155;fXLO2(i)=fXLO2(i)+  0.544;fRO2(i)=fRO2(i)+  0.602;fXOPN(i)=fXOPN(i)+  0.244;fHO(i)=fHO(i)+  0.244;fXO2H(i)=fXO2H(i)+  0.058;fHO2(i)=fHO2(i)+  0.155;fXYLRO2(i)=fXYLRO2(i)+  1.000;

% 190, <R185a>
i=i+1;
Rnames{ 190} = 'NAPH +  HO = 0.15500*CRES +  0.54400*XLO2 +  0.60200*RO2 +  0.24400*XOPN +  0.24400* HO +  0.05800*XO2H +  0.15500*HO2 + PAHRO2 ';
k(:,i) = (  1.8500E-11 ); 
Gstr{i,   1}='NAPH';Gstr{i,   2}='HO';
fNAPH(i)=fNAPH(i)-1.0;fHO(i)=fHO(i)-1.0;
fCRES(i)=fCRES(i)+  0.155;fXLO2(i)=fXLO2(i)+  0.544;fRO2(i)=fRO2(i)+  0.602;fXOPN(i)=fXOPN(i)+  0.244;fHO(i)=fHO(i)+  0.244;fXO2H(i)=fXO2H(i)+  0.058;fHO2(i)=fHO2(i)+  0.155;fPAHRO2(i)=fPAHRO2(i)+  1.000;

% 191, <R186>
i=i+1;
Rnames{ 191} = 'XLO2 + NO = 0.86000*NO2 +  0.14000*NTR2 +  0.22100*GLY +  0.67500*MGLY +  0.30000*OPEN +  0.56000*XOPN +  0.86000*HO2 ';
k(:,i) = (  2.7000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='XLO2';Gstr{i,   2}='NO';
fXLO2(i)=fXLO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  0.860;fNTR2(i)=fNTR2(i)+  0.140;fGLY(i)=fGLY(i)+  0.221;fMGLY(i)=fMGLY(i)+  0.675;fOPEN(i)=fOPEN(i)+  0.300;fXOPN(i)=fXOPN(i)+  0.560;fHO2(i)=fHO2(i)+  0.860;

% 192, <R187>
i=i+1;
Rnames{ 192} = 'XLO2 + HO2 =';
k(:,i) = (  1.9000E-13.*exp(  1.3000E+03./T) ); 
Gstr{i,   1}='XLO2';Gstr{i,   2}='HO2';
fXLO2(i)=fXLO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;


% 193, <R188>
i=i+1;
Rnames{ 193} = 'XLO2 + ACO3 = 0.26000*GLY +  0.77000*MGLY +  0.35000*OPEN +  0.65000*XOPN + HO2 + MO2 + RO2 ';
k(:,i) = (k(:,  58) ); 
Gstr{i,   1}='XLO2';Gstr{i,   2}='ACO3';
fXLO2(i)=fXLO2(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fGLY(i)=fGLY(i)+  0.260;fMGLY(i)=fMGLY(i)+  0.770;fOPEN(i)=fOPEN(i)+  0.350;fXOPN(i)=fXOPN(i)+  0.650;fHO2(i)=fHO2(i)+  1.000;fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 194, <R189>
i=i+1;
Rnames{ 194} = 'XLO2 + RO2 = 0.26000*GLY +  0.77000*MGLY +  0.35000*OPEN +  0.65000*XOPN + HO2 + RO2 ';
k(:,i) = (k(:,  70) ); 
Gstr{i,   1}='XLO2';Gstr{i,   2}='RO2';
fXLO2(i)=fXLO2(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fGLY(i)=fGLY(i)+  0.260;fMGLY(i)=fMGLY(i)+  0.770;fOPEN(i)=fOPEN(i)+  0.350;fXOPN(i)=fXOPN(i)+  0.650;fHO2(i)=fHO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 195, <R190>
i=i+1;
Rnames{ 195} = 'CRES +  HO = 0.02500*GLY +  0.02500*OPEN + HO2 +  0.20000*CRO +  0.73200*CAT1 +  0.02000*XO2N +  0.02000*RO2 ';
k(:,i) = (  1.7000E-12.*exp(  9.5000E+02./T) ); 
Gstr{i,   1}='CRES';Gstr{i,   2}='HO';
fCRES(i)=fCRES(i)-1.0;fHO(i)=fHO(i)-1.0;
fGLY(i)=fGLY(i)+  0.025;fOPEN(i)=fOPEN(i)+  0.025;fHO2(i)=fHO2(i)+  1.000;fCRO(i)=fCRO(i)+  0.200;fCAT1(i)=fCAT1(i)+  0.732;fXO2N(i)=fXO2N(i)+  0.020;fRO2(i)=fRO2(i)+  0.020;

% 196, <R191>
i=i+1;
Rnames{ 196} = 'CRES + NO3 = 0.30000*CRO + HNO3 +  0.48000*XO2 +  0.12000*XO2H +  0.24000*GLY +  0.24000*MGLY +  0.48000*OPO3 +  0.10000*XO2N +  0.70000*RO2 ';
k(:,i) = (  1.4000E-11 ); 
Gstr{i,   1}='CRES';Gstr{i,   2}='NO3';
fCRES(i)=fCRES(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fCRO(i)=fCRO(i)+  0.300;fHNO3(i)=fHNO3(i)+  1.000;fXO2(i)=fXO2(i)+  0.480;fXO2H(i)=fXO2H(i)+  0.120;fGLY(i)=fGLY(i)+  0.240;fMGLY(i)=fMGLY(i)+  0.240;fOPO3(i)=fOPO3(i)+  0.480;fXO2N(i)=fXO2N(i)+  0.100;fRO2(i)=fRO2(i)+  0.700;

% 197, <R192>
i=i+1;
Rnames{ 197} = 'CRO + NO2 = CRON ';
k(:,i) = (  2.1000E-12 ); 
Gstr{i,   1}='CRO';Gstr{i,   2}='NO2';
fCRO(i)=fCRO(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fCRON(i)=fCRON(i)+  1.000;

% 198, <R193>
i=i+1;
Rnames{ 198} = 'CRO + HO2 = CRES ';
k(:,i) = (  5.5000E-12 ); 
Gstr{i,   1}='CRO';Gstr{i,   2}='HO2';
fCRO(i)=fCRO(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fCRES(i)=fCRES(i)+  1.000;

% 199, <R194>
i=i+1;
Rnames{ 199} = 'CRON +  HO = NTR2 +  0.50000*CRO ';
k(:,i) = (  1.5300E-12 ); 
Gstr{i,   1}='CRON';Gstr{i,   2}='HO';
fCRON(i)=fCRON(i)-1.0;fHO(i)=fHO(i)-1.0;
fNTR2(i)=fNTR2(i)+  1.000;fCRO(i)=fCRO(i)+  0.500;

% 200, <R195>
i=i+1;
Rnames{ 200} = 'CRON + NO3 = NTR2 +  0.50000*CRO + HNO3 ';
k(:,i) = (  3.8000E-12 ); 
Gstr{i,   1}='CRON';Gstr{i,   2}='NO3';
fCRON(i)=fCRON(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNTR2(i)=fNTR2(i)+  1.000;fCRO(i)=fCRO(i)+  0.500;fHNO3(i)=fHNO3(i)+  1.000;

% 201, <R196>
i=i+1;
Rnames{ 201} = 'CRON = HONO + HO2 + HCHO + OPEN ';
k(:,i) = (JNTR_IUPAC10 ); 
Gstr{i,   1}='CRON';
fCRON(i)=fCRON(i)-1.0;
fHONO(i)=fHONO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fHCHO(i)=fHCHO(i)+  1.000;fOPEN(i)=fOPEN(i)+  1.000;

% 202, <R197>
i=i+1;
Rnames{ 202} = 'XOPN = 0.40000*GLY + XO2H +  0.70000*HO2 +  0.70000*CO +  0.30000*ACO3 ';
k(:,i) = (  5.0000E-02.*JNO2_IUPAC10 ); 
Gstr{i,   1}='XOPN';
fXOPN(i)=fXOPN(i)-1.0;
fGLY(i)=fGLY(i)+  0.400;fXO2H(i)=fXO2H(i)+  1.000;fHO2(i)=fHO2(i)+  0.700;fCO(i)=fCO(i)+  0.700;fACO3(i)=fACO3(i)+  0.300;

% 203, <R198>
i=i+1;
Rnames{ 203} = 'XOPN +  HO = MGLY +  0.40000*GLY +  2.00000*XO2H +  2.00000*RO2 ';
k(:,i) = (  9.0000E-11 ); 
Gstr{i,   1}='XOPN';Gstr{i,   2}='HO';
fXOPN(i)=fXOPN(i)-1.0;fHO(i)=fHO(i)-1.0;
fMGLY(i)=fMGLY(i)+  1.000;fGLY(i)=fGLY(i)+  0.400;fXO2H(i)=fXO2H(i)+  2.000;fRO2(i)=fRO2(i)+  2.000;

% 204, <R199>
i=i+1;
Rnames{ 204} = 'XOPN + O3 = 1.20000*MGLY +  0.50000* HO +  0.60000*ACO3 +  0.10000*ALD2 +  0.50000*CO +  0.30000*XO2H +  0.30000*RO2 ';
k(:,i) = (  1.0800E-16.*exp( -5.0000E+02./T) ); 
Gstr{i,   1}='XOPN';Gstr{i,   2}='O3';
fXOPN(i)=fXOPN(i)-1.0;fO3(i)=fO3(i)-1.0;
fMGLY(i)=fMGLY(i)+  1.200;fHO(i)=fHO(i)+  0.500;fACO3(i)=fACO3(i)+  0.600;fALD2(i)=fALD2(i)+  0.100;fCO(i)=fCO(i)+  0.500;fXO2H(i)=fXO2H(i)+  0.300;fRO2(i)=fRO2(i)+  0.300;

% 205, <R200>
i=i+1;
Rnames{ 205} = 'XOPN + NO3 = 0.50000*NO2 +  0.50000*NTR2 +  0.45000*XO2H +  0.45000*XO2 +  0.10000*XO2N + RO2 +  0.25000*OPEN +  0.25000*MGLY ';
k(:,i) = (  3.0000E-12 ); 
Gstr{i,   1}='XOPN';Gstr{i,   2}='NO3';
fXOPN(i)=fXOPN(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO2(i)=fNO2(i)+  0.500;fNTR2(i)=fNTR2(i)+  0.500;fXO2H(i)=fXO2H(i)+  0.450;fXO2(i)=fXO2(i)+  0.450;fXO2N(i)=fXO2N(i)+  0.100;fRO2(i)=fRO2(i)+  1.000;fOPEN(i)=fOPEN(i)+  0.250;fMGLY(i)=fMGLY(i)+  0.250;

% 206, <R201>
i=i+1;
Rnames{ 206} = 'OPEN = OPO3 + HO2 + CO ';
k(:,i) = (  2.8000E-02.*JNO2_IUPAC10 ); 
Gstr{i,   1}='OPEN';
fOPEN(i)=fOPEN(i)-1.0;
fOPO3(i)=fOPO3(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fCO(i)=fCO(i)+  1.000;

% 207, <R202>
i=i+1;
Rnames{ 207} = 'OPEN +  HO = 0.60000*OPO3 +  0.40000*XO2H +  0.40000*RO2 +  0.40000*GLY ';
k(:,i) = (  4.4000E-11 ); 
Gstr{i,   1}='OPEN';Gstr{i,   2}='HO';
fOPEN(i)=fOPEN(i)-1.0;fHO(i)=fHO(i)-1.0;
fOPO3(i)=fOPO3(i)+  0.600;fXO2H(i)=fXO2H(i)+  0.400;fRO2(i)=fRO2(i)+  0.400;fGLY(i)=fGLY(i)+  0.400;

% 208, <R203>
i=i+1;
Rnames{ 208} = 'OPEN + O3 = 1.40000*GLY +  0.24000*MGLY +  0.50000* HO +  0.12000*ACO3 +  0.08000*HCHO +  0.02000*ALD2 +  1.98000*CO +  0.56000*HO2 ';
k(:,i) = (  5.4000E-17.*exp( -5.0000E+02./T) ); 
Gstr{i,   1}='OPEN';Gstr{i,   2}='O3';
fOPEN(i)=fOPEN(i)-1.0;fO3(i)=fO3(i)-1.0;
fGLY(i)=fGLY(i)+  1.400;fMGLY(i)=fMGLY(i)+  0.240;fHO(i)=fHO(i)+  0.500;fACO3(i)=fACO3(i)+  0.120;fHCHO(i)=fHCHO(i)+  0.080;fALD2(i)=fALD2(i)+  0.020;fCO(i)=fCO(i)+  1.980;fHO2(i)=fHO2(i)+  0.560;

% 209, <R204>
i=i+1;
Rnames{ 209} = 'OPEN + NO3 = OPO3 + HNO3 ';
k(:,i) = (  3.8000E-12 ); 
Gstr{i,   1}='OPEN';Gstr{i,   2}='NO3';
fOPEN(i)=fOPEN(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fOPO3(i)=fOPO3(i)+  1.000;fHNO3(i)=fHNO3(i)+  1.000;

% 210, <R205>
i=i+1;
Rnames{ 210} = 'CAT1 +  HO = 0.14000*HCHO +  0.20000*HO2 +  0.50000*CRO ';
k(:,i) = (  5.0000E-11 ); 
Gstr{i,   1}='CAT1';Gstr{i,   2}='HO';
fCAT1(i)=fCAT1(i)-1.0;fHO(i)=fHO(i)-1.0;
fHCHO(i)=fHCHO(i)+  0.140;fHO2(i)=fHO2(i)+  0.200;fCRO(i)=fCRO(i)+  0.500;

% 211, <R206>
i=i+1;
Rnames{ 211} = 'CAT1 + NO3 = CRO + HNO3 ';
k(:,i) = (  1.7000E-10 ); 
Gstr{i,   1}='CAT1';Gstr{i,   2}='NO3';
fCAT1(i)=fCAT1(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fCRO(i)=fCRO(i)+  1.000;fHNO3(i)=fHNO3(i)+  1.000;

% 212, <R207>
i=i+1;
Rnames{ 212} = 'OPO3 + NO = NO2 +  0.50000*GLY +  0.50000*CO +  0.80000*HO2 +  0.20000*CXO3 ';
k(:,i) = (  1.0000E-11 ); 
Gstr{i,   1}='OPO3';Gstr{i,   2}='NO';
fOPO3(i)=fOPO3(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO2(i)=fNO2(i)+  1.000;fGLY(i)=fGLY(i)+  0.500;fCO(i)=fCO(i)+  0.500;fHO2(i)=fHO2(i)+  0.800;fCXO3(i)=fCXO3(i)+  0.200;

% 213, <R208>
i=i+1;
Rnames{ 213} = 'OPO3 + NO2 = OPAN ';
k(:,i) = (k(:,  54) ); 
Gstr{i,   1}='OPO3';Gstr{i,   2}='NO2';
fOPO3(i)=fOPO3(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fOPAN(i)=fOPAN(i)+  1.000;

% 214, <R209>
i=i+1;
Rnames{ 214} = 'OPAN = OPO3 + NO2 ';
k(:,i) = (k(:,  55) ); 
Gstr{i,   1}='OPAN';
fOPAN(i)=fOPAN(i)-1.0;
fOPO3(i)=fOPO3(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

% 215, <R210>
i=i+1;
Rnames{ 215} = 'OPO3 + HO2 = 0.41000*PACD +  0.15000*AACD +  0.15000*O3 +  0.44000*ALDX +  0.44000*XO2H +  0.44000*RO2 +  0.44000* HO ';
k(:,i) = (k(:,  57) ); 
Gstr{i,   1}='OPO3';Gstr{i,   2}='HO2';
fOPO3(i)=fOPO3(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fPACD(i)=fPACD(i)+  0.410;fAACD(i)=fAACD(i)+  0.150;fO3(i)=fO3(i)+  0.150;fALDX(i)=fALDX(i)+  0.440;fXO2H(i)=fXO2H(i)+  0.440;fRO2(i)=fRO2(i)+  0.440;fHO(i)=fHO(i)+  0.440;

% 216, <R211>
i=i+1;
Rnames{ 216} = 'OPO3 + ACO3 = MO2 + XO2 + ALDX +  2.00000*RO2 ';
k(:,i) = (k(:,  59) ); 
Gstr{i,   1}='OPO3';Gstr{i,   2}='ACO3';
fOPO3(i)=fOPO3(i)-1.0;fACO3(i)=fACO3(i)-1.0;
fMO2(i)=fMO2(i)+  1.000;fXO2(i)=fXO2(i)+  1.000;fALDX(i)=fALDX(i)+  1.000;fRO2(i)=fRO2(i)+  2.000;

% 217, <R212>
i=i+1;
Rnames{ 217} = 'OPO3 + RO2 = 0.80000*XO2H +  0.80000*ALDX +  1.80000*RO2 +  0.20000*AACD ';
k(:,i) = (k(:,  58) ); 
Gstr{i,   1}='OPO3';Gstr{i,   2}='RO2';
fOPO3(i)=fOPO3(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fXO2H(i)=fXO2H(i)+  0.800;fALDX(i)=fALDX(i)+  0.800;fRO2(i)=fRO2(i)+  1.800;fAACD(i)=fAACD(i)+  0.200;

% 218, <R213>
i=i+1;
Rnames{ 218} = 'OPAN +  HO = 0.50000*NO2 +  0.50000*GLY + CO +  0.50000*NTR2 ';
k(:,i) = (  3.6000E-11 ); 
Gstr{i,   1}='OPAN';Gstr{i,   2}='HO';
fOPAN(i)=fOPAN(i)-1.0;fHO(i)=fHO(i)-1.0;
fNO2(i)=fNO2(i)+  0.500;fGLY(i)=fGLY(i)+  0.500;fCO(i)=fCO(i)+  1.000;fNTR2(i)=fNTR2(i)+  0.500;

% 219, <R214>
i=i+1;
Rnames{ 219} = 'PANX +  HO = ALD2 + NO2 ';
k(:,i) = (  3.0000E-12 ); 
Gstr{i,   1}='PANX';Gstr{i,   2}='HO';
fPANX(i)=fPANX(i)-1.0;fHO(i)=fHO(i)-1.0;
fALD2(i)=fALD2(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

% 220, <R216>
i=i+1;
Rnames{ 220} = 'ECH4 +  HO = MO2 + RO2 ';
k(:,i) = (  1.8500E-12.*exp( -1.6900E+03./T) ); 
Gstr{i,   1}='ECH4';Gstr{i,   2}='HO';
fECH4(i)=fECH4(i)-1.0;fHO(i)=fHO(i)-1.0;
fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 221, <R217>
i=i+1;
Rnames{ 221} = 'XPRP = XO2N + RO2 ';
xko =   2.3700E-21.*M.*exp(  0.0000E+00./T).*(T./300).^  0.0000E+00;
xkinf =   4.3000E-01.*exp(  0.0000E+00./T).*(T./300).^ -8.0000E+00;
xn =   1.0000E+00;
F =   4.1000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='XPRP';
fXPRP(i)=fXPRP(i)-1.0;
fXO2N(i)=fXO2N(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 222, <R218>
i=i+1;
Rnames{ 222} = 'XPRP = 0.73200*ACET +  0.26800*ALDX +  0.26800*PAR + XO2H + RO2 ';
k(:,i) = (  1.0000E+00 ); 
Gstr{i,   1}='XPRP';
fXPRP(i)=fXPRP(i)-1.0;
fACET(i)=fACET(i)+  0.732;fALDX(i)=fALDX(i)+  0.268;fPAR(i)=fPAR(i)+  0.268;fXO2H(i)=fXO2H(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 223, <R219>
i=i+1;
Rnames{ 223} = 'XPAR = XO2N + RO2 ';
xko =   4.8100E-20.*M.*exp(  0.0000E+00./T).*(T./300).^  0.0000E+00;
xkinf =   4.3000E-01.*exp(  0.0000E+00./T).*(T./300).^ -8.0000E+00;
xn =   1.0000E+00;
F =   4.1000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='XPAR';
fXPAR(i)=fXPAR(i)-1.0;
fXO2N(i)=fXO2N(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 224, <R220>
i=i+1;
Rnames{ 224} = 'XPAR = 0.12600*ALDX +  0.87400*ROR +  0.12600*XO2H +  0.87400*XO2 + RO2 -  0.12600*PAR ';
k(:,i) = (  1.0000E+00 ); 
Gstr{i,   1}='XPAR';
fXPAR(i)=fXPAR(i)-1.0;
fALDX(i)=fALDX(i)+  0.126;fROR(i)=fROR(i)+  0.874;fXO2H(i)=fXO2H(i)+  0.126;fXO2(i)=fXO2(i)+  0.874;fRO2(i)=fRO2(i)+  1.000;fPAR(i)=fPAR(i)-  0.126;

% 225, <CL1>
i=i+1;
Rnames{ 225} = 'CL2 = 2.00000*CL ';
k(:,i) = 1E-10;%(JCL2_IUPAC04 ); 
Gstr{i,   1}='CL2';
fCL2(i)=fCL2(i)-1.0;
fCL(i)=fCL(i)+  2.000;

% 226, <CL2>
i=i+1;
Rnames{ 226} = 'HOCL =  HO + CL ';
k(:,i) = 1E-10;%(JHOCL_IUPAC04 ); 
Gstr{i,   1}='HOCL';
fHOCL(i)=fHOCL(i)-1.0;
fHO(i)=fHO(i)+  1.000;fCL(i)=fCL(i)+  1.000;

% 227, <CL3>
i=i+1;
Rnames{ 227} = 'CL + O3 = CLO ';
k(:,i) = (  2.3000E-11.*exp( -2.0000E+02./T) ); 
Gstr{i,   1}='CL';Gstr{i,   2}='O3';
fCL(i)=fCL(i)-1.0;fO3(i)=fO3(i)-1.0;
fCLO(i)=fCLO(i)+  1.000;

% 228, <CL4>
i=i+1;
Rnames{ 228} = 'CLO + CLO = 0.30000*CL2 +  1.40000*CL ';
k(:,i) = (  1.6300E-14 ); 
Gstr{i,   1}='CLO';Gstr{i,   2}='CLO';
fCLO(i)=fCLO(i)-1.0;fCLO(i)=fCLO(i)-1.0;
fCL2(i)=fCL2(i)+  0.300;fCL(i)=fCL(i)+  1.400;

% 229, <CL5>
i=i+1;
Rnames{ 229} = 'CLO + NO = CL + NO2 ';
k(:,i) = (  6.4000E-12.*exp(  2.9000E+02./T) ); 
Gstr{i,   1}='CLO';Gstr{i,   2}='NO';
fCLO(i)=fCLO(i)-1.0;fNO(i)=fNO(i)-1.0;
fCL(i)=fCL(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

% 230, <CL6>
i=i+1;
Rnames{ 230} = 'CLO + HO2 = HOCL ';
k(:,i) = (  2.2000E-12.*exp(  3.4000E+02./T) ); 
Gstr{i,   1}='CLO';Gstr{i,   2}='HO2';
fCLO(i)=fCLO(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHOCL(i)=fHOCL(i)+  1.000;

% 231, <CL7>
i=i+1;
Rnames{ 231} = 'CLO + MO2 = CL + HCHO + HO2 ';
k(:,i) = (  3.2000E-12.*exp( -1.1000E+02./T) ); 
Gstr{i,   1}='CLO';Gstr{i,   2}='MO2';
fCLO(i)=fCLO(i)-1.0;fMO2(i)=fMO2(i)-1.0;
fCL(i)=fCL(i)+  1.000;fHCHO(i)=fHCHO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

% 232, <CL8>
i=i+1;
Rnames{ 232} = ' HO + FMCL = CL + CO ';
k(:,i) = (  5.0000E-13 ); 
Gstr{i,   1}='HO';Gstr{i,   2}='FMCL';
fHO(i)=fHO(i)-1.0;fFMCL(i)=fFMCL(i)-1.0;
fCL(i)=fCL(i)+  1.000;fCO(i)=fCO(i)+  1.000;

% 233, <CL9>
i=i+1;
Rnames{ 233} = 'FMCL = CL + CO + HO2 ';
k(:,i) = 1E-10;%(JFMCL_IUPAC04 ); 
Gstr{i,   1}='FMCL';
fFMCL(i)=fFMCL(i)-1.0;
fCL(i)=fCL(i)+  1.000;fCO(i)=fCO(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;

% 234, <CL10>
i=i+1;
Rnames{ 234} = 'CL + CH4 = HCL + MO2 + RO2 ';
k(:,i) = (  6.6000E-12.*exp( -1.2400E+03./T) ).*CH4; 
Gstr{i,   1}='CL';
fCL(i)=fCL(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fMO2(i)=fMO2(i)+  1.000;fRO2(i)=fRO2(i)+  1.000;

% 235, <CL11>
i=i+1;
Rnames{ 235} = 'CL + PAR = HCL + XPAR ';
k(:,i) = (  5.0000E-11 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='PAR';
fCL(i)=fCL(i)-1.0;fPAR(i)=fPAR(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fXPAR(i)=fXPAR(i)+  1.000;

% 236, <CL12>
i=i+1;
Rnames{ 236} = 'CL + PRPA = HCL + ACET +  0.97000*XO2H +  0.03000*XO2N + RO2 ';
k(:,i) = (  1.4000E-10 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='PRPA';
fCL(i)=fCL(i)-1.0;fPRPA(i)=fPRPA(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fACET(i)=fACET(i)+  1.000;fXO2H(i)=fXO2H(i)+  0.970;fXO2N(i)=fXO2N(i)+  0.030;fRO2(i)=fRO2(i)+  1.000;

% 237, <CL13>
i=i+1;
Rnames{ 237} = 'CL + ETHA = HCL +  0.99100*ALD2 +  0.99100*XO2H +  0.00900*XO2N + RO2 ';
k(:,i) = (  8.3000E-11.*exp( -1.0000E+02./T) ); 
Gstr{i,   1}='CL';Gstr{i,   2}='ETHA';
fCL(i)=fCL(i)-1.0;fETHA(i)=fETHA(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fALD2(i)=fALD2(i)+  0.991;fXO2H(i)=fXO2H(i)+  0.991;fXO2N(i)=fXO2N(i)+  0.009;fRO2(i)=fRO2(i)+  1.000;

% 238, <CL14>
i=i+1;
Rnames{ 238} = 'CL + ETH = FMCL +  2.00000*XO2 + HO2 + HCHO ';
k(:,i) = (  1.0700E-10 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='ETH';
fCL(i)=fCL(i)-1.0;fETH(i)=fETH(i)-1.0;
fFMCL(i)=fFMCL(i)+  1.000;fXO2(i)=fXO2(i)+  2.000;fHO2(i)=fHO2(i)+  1.000;fHCHO(i)=fHCHO(i)+  1.000;

% 239, <CL15>
i=i+1;
Rnames{ 239} = 'CL + OLE = FMCL +  0.33000*ALD2 +  0.67000*ALDX +  2.00000*XO2 + HO2 - PAR ';
k(:,i) = (  2.5000E-10 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='OLE';
fCL(i)=fCL(i)-1.0;fOLE(i)=fOLE(i)-1.0;
fFMCL(i)=fFMCL(i)+  1.000;fALD2(i)=fALD2(i)+  0.330;fALDX(i)=fALDX(i)+  0.670;fXO2(i)=fXO2(i)+  2.000;fHO2(i)=fHO2(i)+  1.000;fPAR(i)=fPAR(i)-  1.000;

% 240, <CL16>
i=i+1;
Rnames{ 240} = 'CL + IOLE = 0.30000*HCL +  0.70000*FMCL +  0.45000*ALD2 +  0.55000*ALDX +  0.30000*OLE +  0.30000*PAR +  1.70000*XO2 + HO2 ';
k(:,i) = (  3.5000E-10 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='IOLE';
fCL(i)=fCL(i)-1.0;fIOLE(i)=fIOLE(i)-1.0;
fHCL(i)=fHCL(i)+  0.300;fFMCL(i)=fFMCL(i)+  0.700;fALD2(i)=fALD2(i)+  0.450;fALDX(i)=fALDX(i)+  0.550;fOLE(i)=fOLE(i)+  0.300;fPAR(i)=fPAR(i)+  0.300;fXO2(i)=fXO2(i)+  1.700;fHO2(i)=fHO2(i)+  1.000;

% 241, <CL17>
i=i+1;
Rnames{ 241} = 'CL + ISOP = FMCL + ISPD +  0.96000*XO2H +  0.04000*XO2N + RO2 ';
k(:,i) = (  4.3000E-10 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='ISOP';
fCL(i)=fCL(i)-1.0;fISOP(i)=fISOP(i)-1.0;
fFMCL(i)=fFMCL(i)+  1.000;fISPD(i)=fISPD(i)+  1.000;fXO2H(i)=fXO2H(i)+  0.960;fXO2N(i)=fXO2N(i)+  0.040;fRO2(i)=fRO2(i)+  1.000;

% 242, <CL18>
i=i+1;
Rnames{ 242} = 'CL + HCHO = HCL + HO2 + CO ';
k(:,i) = (  8.2000E-11.*exp( -3.4000E+01./T) ); 
Gstr{i,   1}='CL';Gstr{i,   2}='HCHO';
fCL(i)=fCL(i)-1.0;fHCHO(i)=fHCHO(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fCO(i)=fCO(i)+  1.000;

% 243, <CL19>
i=i+1;
Rnames{ 243} = 'CL + ALD2 = HCL + ACO3 ';
k(:,i) = (  7.9000E-11 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='ALD2';
fCL(i)=fCL(i)-1.0;fALD2(i)=fALD2(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fACO3(i)=fACO3(i)+  1.000;

% 244, <CL20>
i=i+1;
Rnames{ 244} = 'CL + ALDX = HCL + CXO3 ';
k(:,i) = (  1.3000E-10 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='ALDX';
fCL(i)=fCL(i)-1.0;fALDX(i)=fALDX(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fCXO3(i)=fCXO3(i)+  1.000;

% 245, <CL21>
i=i+1;
Rnames{ 245} = 'CL + ME HO = HCL + HO2 + HCHO ';
k(:,i) = (  5.5000E-11 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='MEOH';
fCL(i)=fCL(i)-1.0;fMEOH(i)=fMEOH(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fHCHO(i)=fHCHO(i)+  1.000;

% 246, <CL22>
i=i+1;
Rnames{ 246} = 'CL + ET HO = HCL + HO2 + ALD2 ';
k(:,i) = (  8.2000E-11.*exp(  4.5000E+01./T) ); 
Gstr{i,   1}='CL';Gstr{i,   2}='ETOH';
fCL(i)=fCL(i)-1.0;fETOH(i)=fETOH(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fHO2(i)=fHO2(i)+  1.000;fALD2(i)=fALD2(i)+  1.000;

% 247, <CL23>
i=i+1;
Rnames{ 247} = 'HCL +  HO = CL ';
k(:,i) = (  6.5800E-13.*exp(  5.8000E+01./T).*(T./300).^(  1.1600E+00 ) ); 
Gstr{i,   1}='HCL';Gstr{i,   2}='HO';
fHCL(i)=fHCL(i)-1.0;fHO(i)=fHO(i)-1.0;
fCL(i)=fCL(i)+  1.000;

% 248, <CL24>
i=i+1;
Rnames{ 248} = 'CL + TOL = HCL +  0.18000*CRES +  0.65000*TO2 +  0.72000*RO2 +  0.10000*OPEN +  0.10000* HO +  0.07000*XO2H +  0.18000*HO2 + TOLRO2 ';
k(:,i) = (  6.1000E-11 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='TOL';
fCL(i)=fCL(i)-1.0;fTOL(i)=fTOL(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fCRES(i)=fCRES(i)+  0.180;fTO2(i)=fTO2(i)+  0.650;fRO2(i)=fRO2(i)+  0.720;fOPEN(i)=fOPEN(i)+  0.100;fHO(i)=fHO(i)+  0.100;fXO2H(i)=fXO2H(i)+  0.070;fHO2(i)=fHO2(i)+  0.180;fTOLRO2(i)=fTOLRO2(i)+  1.000;

% 249, <CL25>
i=i+1;
Rnames{ 249} = 'CL + XYLMN = HCL +  0.15500*CRES +  0.54400*XLO2 +  0.60200*RO2 +  0.24400*XOPN +  0.24400* HO +  0.05800*XO2H +  0.15500*HO2 + XYLRO2 ';
k(:,i) = (  1.2000E-10 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='XYLMN';
fCL(i)=fCL(i)-1.0;fXYLMN(i)=fXYLMN(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fCRES(i)=fCRES(i)+  0.155;fXLO2(i)=fXLO2(i)+  0.544;fRO2(i)=fRO2(i)+  0.602;fXOPN(i)=fXOPN(i)+  0.244;fHO(i)=fHO(i)+  0.244;fXO2H(i)=fXO2H(i)+  0.058;fHO2(i)=fHO2(i)+  0.155;fXYLRO2(i)=fXYLRO2(i)+  1.000;

% 250, <CL26>
i=i+1;
Rnames{ 250} = 'CL + NAPH = HCL +  0.15500*CRES +  0.54400*XLO2 +  0.60200*RO2 +  0.24400*XOPN +  0.24400* HO +  0.05800*XO2H +  0.15500*HO2 + PAHRO2 ';
k(:,i) = (  1.2000E-10 ); 
Gstr{i,   1}='CL';Gstr{i,   2}='NAPH';
fCL(i)=fCL(i)-1.0;fNAPH(i)=fNAPH(i)-1.0;
fHCL(i)=fHCL(i)+  1.000;fCRES(i)=fCRES(i)+  0.155;fXLO2(i)=fXLO2(i)+  0.544;fRO2(i)=fRO2(i)+  0.602;fXOPN(i)=fXOPN(i)+  0.244;fHO(i)=fHO(i)+  0.244;fXO2H(i)=fXO2H(i)+  0.058;fHO2(i)=fHO2(i)+  0.155;fPAHRO2(i)=fPAHRO2(i)+  1.000;

% 251, <CL27>
i=i+1;
Rnames{ 251} = 'CLNO2 = CL + NO2 ';
k(:,i) = 1E-10;%(JCLNO2_IUPAC13 ); 
Gstr{i,   1}='CLNO2';
fCLNO2(i)=fCLNO2(i)-1.0;
fCL(i)=fCL(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

% 252, <CL28>
i=i+1;
Rnames{ 252} = 'CLO + NO2 = CLNO3 ';
xko =   1.8000E-31.*M.*exp(  0.0000E+00./T).*(T./300).^ -3.4000E+00;
xkinf =   1.5000E-11.*exp(  0.0000E+00./T).*(T./300).^ -1.9000E+00;
xn =   1.0000E+00;
F =   6.0000E-01;
G=1.0./(1.0+(log10(xko./xkinf)./xn).^2);
k(:,i) = (xko./( 1.0+xko./xkinf).*F.^G ); 
Gstr{i,   1}='CLO';Gstr{i,   2}='NO2';
fCLO(i)=fCLO(i)-1.0;fNO2(i)=fNO2(i)-1.0;
fCLNO3(i)=fCLNO3(i)+  1.000;

% 253, <CL30>
i=i+1;
Rnames{ 253} = 'CLNO3 = CLO + NO2 ';
k(:,i) = 1E-10;%(JCLONO2_1 ); 
Gstr{i,   1}='CLNO3';
fCLNO3(i)=fCLNO3(i)-1.0;
fCLO(i)=fCLO(i)+  1.000;fNO2(i)=fNO2(i)+  1.000;

% 254, <CL31>
i=i+1;
Rnames{ 254} = 'CLNO3 = CL + NO3 ';
k(:,i) = 1E-10;%(JCLONO2_2 ); 
Gstr{i,   1}='CLNO3';
fCLNO3(i)=fCLNO3(i)-1.0;
fCL(i)=fCL(i)+  1.000;fNO3(i)=fNO3(i)+  1.000;

% 255, <HET_CLNO3_WAJ>
i=i+1;
Rnames{ 255} = 'CLNO3 = HOCL + HNO3 ';
k(:,i) = (K_HETERO_CLNO3_WAJ ); 
Gstr{i,   1}='CLNO3';
fCLNO3(i)=fCLNO3(i)-1.0;
fHOCL(i)=fHOCL(i)+  1.000;fHNO3(i)=fHNO3(i)+  1.000;

% 256, <SA01>
i=i+1;
Rnames{ 256} = 'TOLRO2 + NO = NO +  0.01600*SVAVB2 +  0.05100*SVAVB3 +  0.04700*SVAVB4 ';
k(:,i) = (  2.7000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='TOLRO2';Gstr{i,   2}='NO';
fTOLRO2(i)=fTOLRO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO(i)=fNO(i)+  1.000;fSVAVB2(i)=fSVAVB2(i)+  0.016;fSVAVB3(i)=fSVAVB3(i)+  0.051;fSVAVB4(i)=fSVAVB4(i)+  0.047;

% 257, <SA02>
i=i+1;
Rnames{ 257} = 'TOLRO2 + HO2 = HO2 +  0.14000*SVAVB1 ';
k(:,i) = (  1.9000E-13.*exp(  1.3000E+03./T) ); 
Gstr{i,   1}='TOLRO2';Gstr{i,   2}='HO2';
fTOLRO2(i)=fTOLRO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fSVAVB1(i)=fSVAVB1(i)+  0.140;

% 258, <SA03>
i=i+1;
Rnames{ 258} = 'XYLRO2 + NO = NO +  0.01500*SVAVB2 +  0.02300*SVAVB3 +  0.06000*SVAVB4 ';
k(:,i) = (  2.7000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='XYLRO2';Gstr{i,   2}='NO';
fXYLRO2(i)=fXYLRO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO(i)=fNO(i)+  1.000;fSVAVB2(i)=fSVAVB2(i)+  0.015;fSVAVB3(i)=fSVAVB3(i)+  0.023;fSVAVB4(i)=fSVAVB4(i)+  0.060;

% 259, <SA04>
i=i+1;
Rnames{ 259} = 'XYLRO2 + HO2 = HO2 +  0.19300*SVAVB1 ';
k(:,i) = (  1.9000E-13.*exp(  1.3000E+03./T) ); 
Gstr{i,   1}='XYLRO2';Gstr{i,   2}='HO2';
fXYLRO2(i)=fXYLRO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fSVAVB1(i)=fSVAVB1(i)+  0.193;

% 260, <SA06>
i=i+1;
Rnames{ 260} = 'BENZRO2 + NO = NO +  0.03400*SVAVB2 +  0.39200*SVAVB4 ';
k(:,i) = (  2.7000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='BENZRO2';Gstr{i,   2}='NO';
fBENZRO2(i)=fBENZRO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO(i)=fNO(i)+  1.000;fSVAVB2(i)=fSVAVB2(i)+  0.034;fSVAVB4(i)=fSVAVB4(i)+  0.392;

% 261, <SA07>
i=i+1;
Rnames{ 261} = 'BENZRO2 + HO2 = HO2 +  0.14600*SVAVB1 ';
k(:,i) = (  1.9000E-13.*exp(  1.3000E+03./T) ); 
Gstr{i,   1}='BENZRO2';Gstr{i,   2}='HO2';
fBENZRO2(i)=fBENZRO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fSVAVB1(i)=fSVAVB1(i)+  0.146;

% 262, <SA08>
i=i+1;
Rnames{ 262} = 'SESQ + O3 = O3 + SESQRXN ';
k(:,i) = (  1.1600E-14 ); 
Gstr{i,   1}='SESQ';Gstr{i,   2}='O3';
fSESQ(i)=fSESQ(i)-1.0;fO3(i)=fO3(i)-1.0;
fO3(i)=fO3(i)+  1.000;fSESQRXN(i)=fSESQRXN(i)+  1.000;

% 263, <SA09>
i=i+1;
Rnames{ 263} = 'SESQ +  HO =  HO + SESQRXN ';
k(:,i) = (  1.9700E-10 ); 
Gstr{i,   1}='SESQ';Gstr{i,   2}='HO';
fSESQ(i)=fSESQ(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fSESQRXN(i)=fSESQRXN(i)+  1.000;

% 264, <SA10>
i=i+1;
Rnames{ 264} = 'SESQ + NO3 = NO3 + SESQRXN ';
k(:,i) = (  1.9000E-11 ); 
Gstr{i,   1}='SESQ';Gstr{i,   2}='NO3';
fSESQ(i)=fSESQ(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;fSESQRXN(i)=fSESQRXN(i)+  1.000;

% 265, <SA11>
i=i+1;
Rnames{ 265} = 'PAHRO2 + NO = NO +  0.02800*SVAVB2 +  0.22500*SVAVB3 +  0.19100*SVAVB4 ';
k(:,i) = (  2.7000E-12.*exp(  3.6000E+02./T) ); 
Gstr{i,   1}='PAHRO2';Gstr{i,   2}='NO';
fPAHRO2(i)=fPAHRO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO(i)=fNO(i)+  1.000;fSVAVB2(i)=fSVAVB2(i)+  0.028;fSVAVB3(i)=fSVAVB3(i)+  0.225;fSVAVB4(i)=fSVAVB4(i)+  0.191;

% 266, <SA12>
i=i+1;
Rnames{ 266} = 'PAHRO2 + HO2 = HO2 +  0.47300*SVAVB1 ';
k(:,i) = (  1.9000E-13.*exp(  1.3000E+03./T) ); 
Gstr{i,   1}='PAHRO2';Gstr{i,   2}='HO2';
fPAHRO2(i)=fPAHRO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fSVAVB1(i)=fSVAVB1(i)+  0.473;

% 267, <SA13>
i=i+1;
Rnames{ 267} = 'SOAALK +  HO =  HO +  0.00600*SVAVB2 +  0.05200*SVAVB3 +  0.08100*SVAVB4 ';
k(:,i) = (  2.7000E-12.*exp(  3.7400E+02./T) ); 
Gstr{i,   1}='SOAALK';Gstr{i,   2}='HO';
fSOAALK(i)=fSOAALK(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fSVAVB2(i)=fSVAVB2(i)+  0.006;fSVAVB3(i)=fSVAVB3(i)+  0.052;fSVAVB4(i)=fSVAVB4(i)+  0.081;

% 268, <HET_NTR2>
i=i+1;
Rnames{ 268} = 'NTR2 = HNO3 ';
k(:,i) = (  1.4000E+00.*K_HETERO_NTR2 ); 
Gstr{i,   1}='NTR2';
fNTR2(i)=fNTR2(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;

% 269, <HET_N2O5IJ>
i=i+1;
Rnames{ 269} = 'N2O5 = HNO3 + H2NO3PIJ ';
k(:,i) = (K_HETERO_N2O5IJ ); 
Gstr{i,   1}='N2O5';
fN2O5(i)=fN2O5(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;fH2NO3PIJ(i)=fH2NO3PIJ(i)+  1.000;

% 270, <HET_N2O5K>
i=i+1;
Rnames{ 270} = 'N2O5 = HNO3 + H2NO3PK ';
k(:,i) = (K_HETERO_N2O5K ); 
Gstr{i,   1}='N2O5';
fN2O5(i)=fN2O5(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;fH2NO3PK(i)=fH2NO3PK(i)+  1.000;

% 271, <HET_H2NO3PIJA>
i=i+1;
Rnames{ 271} = 'H2NO3PIJ = HNO3 ';
k(:,i) = (K_HETERO_H2NO3PAIJ ); 
Gstr{i,   1}='H2NO3PIJ';
fH2NO3PIJ(i)=fH2NO3PIJ(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;

% 272, <HET_H2NO3PKA>
i=i+1;
Rnames{ 272} = 'H2NO3PK = HNO3 ';
k(:,i) = (K_HETERO_H2NO3PAK ); 
Gstr{i,   1}='H2NO3PK';
fH2NO3PK(i)=fH2NO3PK(i)-1.0;
fHNO3(i)=fHNO3(i)+  1.000;

% 273, <HET_H2NO3PIB>
i=i+1;
Rnames{ 273} = 'H2NO3PIJ + ACLI = CLNO2 ';
k(:,i) = (K_HETERO_H2NO3PBIJ ); 
Gstr{i,   1}='H2NO3PIJ';Gstr{i,   2}='ACLI';
fH2NO3PIJ(i)=fH2NO3PIJ(i)-1.0;fACLI(i)=fACLI(i)-1.0;
fCLNO2(i)=fCLNO2(i)+  1.000;

% 274, <HET_H2NO3PJB>
i=i+1;
Rnames{ 274} = 'H2NO3PIJ + ACLJ = CLNO2 ';
k(:,i) = (K_HETERO_H2NO3PBIJ ); 
Gstr{i,   1}='H2NO3PIJ';Gstr{i,   2}='ACLJ';
fH2NO3PIJ(i)=fH2NO3PIJ(i)-1.0;fACLJ(i)=fACLJ(i)-1.0;
fCLNO2(i)=fCLNO2(i)+  1.000;

% 275, <HET_H2NO3PKB>
i=i+1;
Rnames{ 275} = 'H2NO3PK + ACLK = CLNO2 ';
k(:,i) = (K_HETERO_H2NO3PBK ); 
Gstr{i,   1}='H2NO3PK';Gstr{i,   2}='ACLK';
fH2NO3PK(i)=fH2NO3PK(i)-1.0;fACLK(i)=fACLK(i)-1.0;
fCLNO2(i)=fCLNO2(i)+  1.000;

% 276, <HET_N02>
i=i+1;
Rnames{ 276} = 'NO2 = 0.50000*HONO +  0.50000*HNO3 ';
k(:,i) = (K_HETERO_NO2 ); 
Gstr{i,   1}='NO2';
fNO2(i)=fNO2(i)-1.0;
fHONO(i)=fHONO(i)+  0.500;fHNO3(i)=fHNO3(i)+  0.500;

% 277, <HAL_Ozone>
i=i+1;
Rnames{ 277} = 'O3 =';
ILLUMINATED =  ( SZA > 0.0 );
OPEN_OCEAN  = SZA;
OPEN_OCEAN  = true;
Patm = 0.001.*P;
a =  6.701E-11.*exp( 1.074E+01.*Patm) + 3.415E-08.*exp(-6.713E-01.*Patm);
b =  2.000E-06;
a(a>b) = b;
k(:,i) = a.*ILLUMINATED.*OPEN_OCEAN;
Gstr{i,   1}='O3';
fO3(i)=fO3(i)-1.0;


% 278, <HET_IEPOX>
i=i+1;
Rnames{ 278} = 'EPOX = IEPOXP ';
k(:,i) = (K_HETERO_IEPOX ); 
Gstr{i,   1}='EPOX';
fEPOX(i)=fEPOX(i)-1.0;
fIEPOXP(i)=fIEPOXP(i)+  1.000;

% 279, <HET_IEPOXOS>
i=i+1;
Rnames{ 279} = 'IEPOXP + ASO4J = AISO3J ';
k(:,i) = (K_HETERO_IEPOXOS ); 
Gstr{i,   1}='IEPOXP';Gstr{i,   2}='ASO4J';
fIEPOXP(i)=fIEPOXP(i)-1.0;fASO4J(i)=fASO4J(i)-1.0;
fAISO3J(i)=fAISO3J(i)+  1.000;

% 280, <HET_TETROL>
i=i+1;
Rnames{ 280} = 'IEPOXP = AISO3J ';
k(:,i) = (K_HETERO_TETROL ); 
Gstr{i,   1}='IEPOXP';
fIEPOXP(i)=fIEPOXP(i)-1.0;
fAISO3J(i)=fAISO3J(i)+  1.000;

% 281, <HET_GLY>
i=i+1;
Rnames{ 281} = 'GLY = AGLYJ ';
k(:,i) = (K_HETERO_GLY ); 
Gstr{i,   1}='GLY';
fGLY(i)=fGLY(i)-1.0;
fAGLYJ(i)=fAGLYJ(i)+  1.000;

% 282, <HET_MGLY>
i=i+1;
Rnames{ 282} = 'MGLY = AGLYJ ';
k(:,i) = (K_HETERO_MGLY ); 
Gstr{i,   1}='MGLY';
fMGLY(i)=fMGLY(i)-1.0;
fAGLYJ(i)=fAGLYJ(i)+  1.000;

% 283, <BL18a>
i=i+1;
Rnames{ 283} = 'TERPNRO2 + NO = NO +  0.68800*MTNO3 ';
k(:,i) = (  2.6000E-12.*exp(  3.8000E+02./T) ); 
Gstr{i,   1}='TERPNRO2';Gstr{i,   2}='NO';
fTERPNRO2(i)=fTERPNRO2(i)-1.0;fNO(i)=fNO(i)-1.0;
fNO(i)=fNO(i)+  1.000;fMTNO3(i)=fMTNO3(i)+  0.688;

% 284, <BL18b>
i=i+1;
Rnames{ 284} = 'TERPNRO2 + HO2 = HO2 + MTNO3 ';
k(:,i) = (  2.6500E-13.*exp(  1.3000E+03./T) ); 
Gstr{i,   1}='TERPNRO2';Gstr{i,   2}='HO2';
fTERPNRO2(i)=fTERPNRO2(i)-1.0;fHO2(i)=fHO2(i)-1.0;
fHO2(i)=fHO2(i)+  1.000;fMTNO3(i)=fMTNO3(i)+  1.000;

% 285, <BL18c>
i=i+1;
Rnames{ 285} = 'TERPNRO2 + NO3 = NO3 +  0.42200*MTNO3 ';
k(:,i) = (  2.3000E-12 ); 
Gstr{i,   1}='TERPNRO2';Gstr{i,   2}='NO3';
fTERPNRO2(i)=fTERPNRO2(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;fMTNO3(i)=fMTNO3(i)+  0.422;

% 286, <BL18d>
i=i+1;
Rnames{ 286} = 'TERPNRO2 + RO2 = RO2 +  0.71100*MTNO3 ';
k(:,i) = (  3.5000E-14 ); 
Gstr{i,   1}='TERPNRO2';Gstr{i,   2}='RO2';
fTERPNRO2(i)=fTERPNRO2(i)-1.0;fRO2(i)=fRO2(i)-1.0;
fRO2(i)=fRO2(i)+  1.000;fMTNO3(i)=fMTNO3(i)+  0.711;

% 287, <CP07mtp>
i=i+1;
Rnames{ 287} = 'MTNO3 + CL = CL +  0.37000*MTNO3 ';
k(:,i) = (  1.9200E-10 ); 
Gstr{i,   1}='MTNO3';Gstr{i,   2}='CL';
fMTNO3(i)=fMTNO3(i)-1.0;fCL(i)=fCL(i)-1.0;
fCL(i)=fCL(i)+  1.000;fMTNO3(i)=fMTNO3(i)+  0.370;

% 288, <BP70mtp>
i=i+1;
Rnames{ 288} = 'MTNO3 +  HO =  HO +  0.24000*MTNO3 ';
k(:,i) = (  7.2000E-12 ); 
Gstr{i,   1}='MTNO3';Gstr{i,   2}='HO';
fMTNO3(i)=fMTNO3(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fMTNO3(i)=fMTNO3(i)+  0.240;

% 289, <BP71mtp>
i=i+1;
Rnames{ 289} = 'MTNO3 =';
k(:,i) = 1E-10;%(JIC3ONO2 ); 
Gstr{i,   1}='MTNO3';
fMTNO3(i)=fMTNO3(i)-1.0;


% 290, <HYD_MT>
i=i+1;
Rnames{ 290} = 'AMTNO3J = AMTHYDJ ';
k(:,i) = (  9.2590E-05 ); 
Gstr{i,   1}='AMTNO3J';
fAMTNO3J(i)=fAMTNO3J(i)-1.0;
fAMTHYDJ(i)=fAMTHYDJ(i)+  1.000;

% 291, <OLIG_AROMATIC1>
i=i+1;
Rnames{ 291} = 'AAVB2J = 0.90700*AOLGAJ ';
k(:,i) = (  9.4882E-06 ); 
Gstr{i,   1}='AAVB2J';
fAAVB2J(i)=fAAVB2J(i)-1.0;
fAOLGAJ(i)=fAOLGAJ(i)+  0.907;

% 292, <OLIG_AROMATIC2>
i=i+1;
Rnames{ 292} = 'AAVB3J = 0.92500*AOLGAJ ';
k(:,i) = (  9.4882E-06 ); 
Gstr{i,   1}='AAVB3J';
fAAVB3J(i)=fAAVB3J(i)-1.0;
fAOLGAJ(i)=fAOLGAJ(i)+  0.925;

% 293, <OLIG_AROMATIC3>
i=i+1;
Rnames{ 293} = 'AAVB4J = 0.94300*AOLGAJ ';
k(:,i) = (  9.4882E-06 ); 
Gstr{i,   1}='AAVB4J';
fAAVB4J(i)=fAAVB4J(i)-1.0;
fAOLGAJ(i)=fAOLGAJ(i)+  0.943;

% 294, <OLIG_ISOPRENE1>
i=i+1;
Rnames{ 294} = 'AISO1J = 0.50000*AOLGBJ ';
k(:,i) = (  9.4882E-06 ); 
Gstr{i,   1}='AISO1J';
fAISO1J(i)=fAISO1J(i)-1.0;
fAOLGBJ(i)=fAOLGBJ(i)+  0.500;

% 295, <OLIG_ISOPRENE2>
i=i+1;
Rnames{ 295} = 'AISO2J = 0.50000*AOLGBJ ';
k(:,i) = (  9.4882E-06 ); 
Gstr{i,   1}='AISO2J';
fAISO2J(i)=fAISO2J(i)-1.0;
fAOLGBJ(i)=fAOLGBJ(i)+  0.500;

% 296, <OLIG_SESQT1>
i=i+1;
Rnames{ 296} = 'ASQTJ = 1.50000*AOLGBJ ';
k(:,i) = (  9.4882E-06 ); 
Gstr{i,   1}='ASQTJ';
fASQTJ(i)=fASQTJ(i)-1.0;
fAOLGBJ(i)=fAOLGBJ(i)+  1.500;

% 297, <RPOAGEPI>
i=i+1;
Rnames{ 297} = 'APOCI +  HO = 1.25000*APNCOMI + APOCI +  HO ';
k(:,i) = (  2.5000E-12 ); 
Gstr{i,   1}='APOCI';Gstr{i,   2}='HO';
fAPOCI(i)=fAPOCI(i)-1.0;fHO(i)=fHO(i)-1.0;
fAPNCOMI(i)=fAPNCOMI(i)+  1.250;fAPOCI(i)=fAPOCI(i)+  1.000;fHO(i)=fHO(i)+  1.000;

% 298, <RPOAGELI>
i=i+1;
Rnames{ 298} = 'APNCOMI +  HO =  HO ';
k(:,i) = (K_HETERO_PNCOMLI ); 
Gstr{i,   1}='APNCOMI';Gstr{i,   2}='HO';
fAPNCOMI(i)=fAPNCOMI(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;

% 299, <RPOAGEPJ>
i=i+1;
Rnames{ 299} = 'APOCJ +  HO = 1.25000*APNCOMJ + APOCJ +  HO ';
k(:,i) = (  2.5000E-12 ); 
Gstr{i,   1}='APOCJ';Gstr{i,   2}='HO';
fAPOCJ(i)=fAPOCJ(i)-1.0;fHO(i)=fHO(i)-1.0;
fAPNCOMJ(i)=fAPNCOMJ(i)+  1.250;fAPOCJ(i)=fAPOCJ(i)+  1.000;fHO(i)=fHO(i)+  1.000;

% 300, <RPOAGELJ>
i=i+1;
Rnames{ 300} = 'APNCOMJ +  HO =  HO ';
k(:,i) = (K_HETERO_PNCOMLJ ); 
Gstr{i,   1}='APNCOMJ';Gstr{i,   2}='HO';
fAPNCOMJ(i)=fAPNCOMJ(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;

% 301, <PCSOA>
i=i+1;
Rnames{ 301} = 'PCVOC +  HO =  HO + PCSOARXN ';
k(:,i) = (  1.2500E-11 ); 
Gstr{i,   1}='PCVOC';Gstr{i,   2}='HO';
fPCVOC(i)=fPCVOC(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fPCSOARXN(i)=fPCSOARXN(i)+  1.000;

% 302, <POA_AGE1>
i=i+1;
Rnames{ 302} = 'VLVPO1 +  HO =  HO +  0.48570*VLVPO1 +  0.00620*VSVPO1 +  0.00250*VSVPO2 +  0.00260*VSVPO3 +  0.00230*VIVPO1 +  0.29440*VLVOO1 +  0.20210*VLVOO2 +  0.00190*VSVOO2 +  0.00230*VSVOO3 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VLVPO1';Gstr{i,   2}='HO';
fVLVPO1(i)=fVLVPO1(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVPO1(i)=fVLVPO1(i)+  0.486;fVSVPO1(i)=fVSVPO1(i)+  0.006;fVSVPO2(i)=fVSVPO2(i)+  0.003;fVSVPO3(i)=fVSVPO3(i)+  0.003;fVIVPO1(i)=fVIVPO1(i)+  0.002;fVLVOO1(i)=fVLVOO1(i)+  0.294;fVLVOO2(i)=fVLVOO2(i)+  0.202;fVSVOO2(i)=fVSVOO2(i)+  0.002;fVSVOO3(i)=fVSVOO3(i)+  0.002;

% 303, <POA_AGE2>
i=i+1;
Rnames{ 303} = 'VSVPO1 +  HO =  HO +  0.30030*VLVPO1 +  0.28620*VSVPO1 +  0.00410*VSVPO2 +  0.00350*VSVPO3 +  0.22390*VLVOO1 +  0.18200*VLVOO2 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VSVPO1';Gstr{i,   2}='HO';
fVSVPO1(i)=fVSVPO1(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVPO1(i)=fVLVPO1(i)+  0.300;fVSVPO1(i)=fVSVPO1(i)+  0.286;fVSVPO2(i)=fVSVPO2(i)+  0.004;fVSVPO3(i)=fVSVPO3(i)+  0.004;fVLVOO1(i)=fVLVOO1(i)+  0.224;fVLVOO2(i)=fVLVOO2(i)+  0.182;

% 304, <POA_AGE3>
i=i+1;
Rnames{ 304} = 'VSVPO2 +  HO =  HO +  0.38560*VLVPO1 +  0.09500*VSVPO1 +  0.13730*VSVPO2 +  0.00050*VSVPO3 +  0.20510*VLVOO1 +  0.17640*VLVOO2 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VSVPO2';Gstr{i,   2}='HO';
fVSVPO2(i)=fVSVPO2(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVPO1(i)=fVLVPO1(i)+  0.386;fVSVPO1(i)=fVSVPO1(i)+  0.095;fVSVPO2(i)=fVSVPO2(i)+  0.137;fVSVPO3(i)=fVSVPO3(i)+  0.001;fVLVOO1(i)=fVLVOO1(i)+  0.205;fVLVOO2(i)=fVLVOO2(i)+  0.176;

% 305, <POA_AGE4>
i=i+1;
Rnames{ 305} = 'VSVPO3 +  HO =  HO +  0.21810*VLVPO1 +  0.30630*VSVPO1 +  0.01530*VSVPO2 +  0.10430*VSVPO3 +  0.18930*VLVOO1 +  0.16680*VLVOO2 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VSVPO3';Gstr{i,   2}='HO';
fVSVPO3(i)=fVSVPO3(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVPO1(i)=fVLVPO1(i)+  0.218;fVSVPO1(i)=fVSVPO1(i)+  0.306;fVSVPO2(i)=fVSVPO2(i)+  0.015;fVSVPO3(i)=fVSVPO3(i)+  0.104;fVLVOO1(i)=fVLVOO1(i)+  0.189;fVLVOO2(i)=fVLVOO2(i)+  0.167;

% 306, <POA_AGE5>
i=i+1;
Rnames{ 306} = 'VIVPO1 +  HO =  HO +  0.24120*VLVPO1 +  0.20890*VSVPO1 +  0.30000*VSVPO2 +  0.20280*VLVOO1 +  0.04710*VLVOO2 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VIVPO1';Gstr{i,   2}='HO';
fVIVPO1(i)=fVIVPO1(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVPO1(i)=fVLVPO1(i)+  0.241;fVSVPO1(i)=fVSVPO1(i)+  0.209;fVSVPO2(i)=fVSVPO2(i)+  0.300;fVLVOO1(i)=fVLVOO1(i)+  0.203;fVLVOO2(i)=fVLVOO2(i)+  0.047;

% 307, <POA_AGE6>
i=i+1;
Rnames{ 307} = 'VLVOO1 +  HO =  HO +  0.66640*VLVOO1 +  0.01430*VLVOO2 +  0.01230*VSVOO1 +  0.12390*VSVOO2 +  0.18310*VSVOO3 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VLVOO1';Gstr{i,   2}='HO';
fVLVOO1(i)=fVLVOO1(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVOO1(i)=fVLVOO1(i)+  0.666;fVLVOO2(i)=fVLVOO2(i)+  0.014;fVSVOO1(i)=fVSVOO1(i)+  0.012;fVSVOO2(i)=fVSVOO2(i)+  0.124;fVSVOO3(i)=fVSVOO3(i)+  0.183;

% 308, <POA_AGE7>
i=i+1;
Rnames{ 308} = 'VLVOO2 +  HO =  HO +  0.28580*VLVOO1 +  0.39310*VLVOO2 +  0.01390*VSVOO1 +  0.10270*VSVOO2 +  0.20450*VSVOO3 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VLVOO2';Gstr{i,   2}='HO';
fVLVOO2(i)=fVLVOO2(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVOO1(i)=fVLVOO1(i)+  0.286;fVLVOO2(i)=fVLVOO2(i)+  0.393;fVSVOO1(i)=fVSVOO1(i)+  0.014;fVSVOO2(i)=fVSVOO2(i)+  0.103;fVSVOO3(i)=fVSVOO3(i)+  0.204;

% 309, <POA_AGE8>
i=i+1;
Rnames{ 309} = 'VSVOO1 +  HO =  HO +  0.33030*VLVOO1 +  0.22720*VLVOO2 +  0.26070*VSVOO1 +  0.07020*VSVOO2 +  0.11160*VSVOO3 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VSVOO1';Gstr{i,   2}='HO';
fVSVOO1(i)=fVSVOO1(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVOO1(i)=fVLVOO1(i)+  0.330;fVLVOO2(i)=fVLVOO2(i)+  0.227;fVSVOO1(i)=fVSVOO1(i)+  0.261;fVSVOO2(i)=fVSVOO2(i)+  0.070;fVSVOO3(i)=fVSVOO3(i)+  0.112;

% 310, <POA_AGE9>
i=i+1;
Rnames{ 310} = 'VSVOO2 +  HO =  HO +  0.34440*VLVOO1 +  0.27490*VLVOO2 +  0.04910*VSVOO1 +  0.25770*VSVOO2 +  0.07390*VSVOO3 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VSVOO2';Gstr{i,   2}='HO';
fVSVOO2(i)=fVSVOO2(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVOO1(i)=fVLVOO1(i)+  0.344;fVLVOO2(i)=fVLVOO2(i)+  0.275;fVSVOO1(i)=fVSVOO1(i)+  0.049;fVSVOO2(i)=fVSVOO2(i)+  0.258;fVSVOO3(i)=fVSVOO3(i)+  0.074;

% 311, <POA_AGE10>
i=i+1;
Rnames{ 311} = 'VSVOO3 +  HO =  HO +  0.38860*VLVOO1 +  0.24210*VLVOO2 +  0.06400*VSVOO1 +  0.03850*VSVOO2 +  0.26670*VSVOO3 ';
k(:,i) = (  4.0000E-11 ); 
Gstr{i,   1}='VSVOO3';Gstr{i,   2}='HO';
fVSVOO3(i)=fVSVOO3(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fVLVOO1(i)=fVLVOO1(i)+  0.389;fVLVOO2(i)=fVLVOO2(i)+  0.242;fVSVOO1(i)=fVSVOO1(i)+  0.064;fVSVOO2(i)=fVSVOO2(i)+  0.038;fVSVOO3(i)=fVSVOO3(i)+  0.267;

% 312, <T01>
i=i+1;
Rnames{ 312} = 'HCHO_PRIMARY +  HO =  HO ';
k(:,i) = (  5.4000E-12.*exp(  1.3500E+02./T) ); 
Gstr{i,   1}='HCHO_PRIMARY';Gstr{i,   2}='HO';
fHCHO_PRIMARY(i)=fHCHO_PRIMARY(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;

% 313, <T02>
i=i+1;
Rnames{ 313} = 'HCHO_PRIMARY + NO3 = NO3 ';
k(:,i) = (  5.5000E-16 ); 
Gstr{i,   1}='HCHO_PRIMARY';Gstr{i,   2}='NO3';
fHCHO_PRIMARY(i)=fHCHO_PRIMARY(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;

% 314, <T03>
i=i+1;
Rnames{ 314} = 'HCHO_PRIMARY + O = O ';
k(:,i) = (  3.4000E-11.*exp( -1.6000E+03./T) ); 
Gstr{i,   1}='HCHO_PRIMARY';Gstr{i,   2}='O';
fHCHO_PRIMARY(i)=fHCHO_PRIMARY(i)-1.0;fO(i)=fO(i)-1.0;
fO(i)=fO(i)+  1.000;

% 315, <T04>
i=i+1;
Rnames{ 315} = 'HCHO_PRIMARY =';
k(:,i) = (JFORM_R_IUPAC10 ); 
Gstr{i,   1}='HCHO_PRIMARY';
fHCHO_PRIMARY(i)=fHCHO_PRIMARY(i)-1.0;


% 316, <T05>
i=i+1;
Rnames{ 316} = 'HCHO_PRIMARY =';
k(:,i) = (JFORM_M_IUPAC10 ); 
Gstr{i,   1}='HCHO_PRIMARY';
fHCHO_PRIMARY(i)=fHCHO_PRIMARY(i)-1.0;


% 317, <TCL1>
i=i+1;
Rnames{ 317} = 'HCHO_PRIMARY + CL = CL ';
k(:,i) = (  8.2000E-11.*exp( -3.4000E+01./T) ); 
Gstr{i,   1}='HCHO_PRIMARY';Gstr{i,   2}='CL';
fHCHO_PRIMARY(i)=fHCHO_PRIMARY(i)-1.0;fCL(i)=fCL(i)-1.0;
fCL(i)=fCL(i)+  1.000;

% 318, <T06>
i=i+1;
Rnames{ 318} = 'ALD2_PRIMARY +  HO =  HO ';
k(:,i) = (  4.7000E-12.*exp(  3.4500E+02./T) ); 
Gstr{i,   1}='ALD2_PRIMARY';Gstr{i,   2}='HO';
fALD2_PRIMARY(i)=fALD2_PRIMARY(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;

% 319, <T07>
i=i+1;
Rnames{ 319} = 'ALD2_PRIMARY + NO3 = NO3 ';
k(:,i) = (  1.4000E-12.*exp( -1.8600E+03./T) ); 
Gstr{i,   1}='ALD2_PRIMARY';Gstr{i,   2}='NO3';
fALD2_PRIMARY(i)=fALD2_PRIMARY(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;

% 320, <T08>
i=i+1;
Rnames{ 320} = 'ALD2_PRIMARY + O = O ';
k(:,i) = (  1.8000E-11.*exp( -1.1000E+03./T) ); 
Gstr{i,   1}='ALD2_PRIMARY';Gstr{i,   2}='O';
fALD2_PRIMARY(i)=fALD2_PRIMARY(i)-1.0;fO(i)=fO(i)-1.0;
fO(i)=fO(i)+  1.000;

% 321, <T09>
i=i+1;
Rnames{ 321} = 'ALD2_PRIMARY =';
k(:,i) = (JALD2_R_IUPAC10 ); 
Gstr{i,   1}='ALD2_PRIMARY';
fALD2_PRIMARY(i)=fALD2_PRIMARY(i)-1.0;


% 322, <TCL2>
i=i+1;
Rnames{ 322} = 'ALD2_PRIMARY + CL = CL ';
k(:,i) = (  7.9000E-11 ); 
Gstr{i,   1}='ALD2_PRIMARY';Gstr{i,   2}='CL';
fALD2_PRIMARY(i)=fALD2_PRIMARY(i)-1.0;fCL(i)=fCL(i)-1.0;
fCL(i)=fCL(i)+  1.000;

% 323, <T10>
i=i+1;
Rnames{ 323} = 'BUTADIENE13 +  HO =  HO +  0.58000*ACROLEIN ';
k(:,i) = (  1.4800E-11.*exp(  4.4800E+02./T) ); 
Gstr{i,   1}='BUTADIENE13';Gstr{i,   2}='HO';
fBUTADIENE13(i)=fBUTADIENE13(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;fACROLEIN(i)=fACROLEIN(i)+  0.580;

% 324, <T11>
i=i+1;
Rnames{ 324} = 'BUTADIENE13 + O3 = O3 +  0.52000*ACROLEIN ';
k(:,i) = (  1.3400E-14.*exp( -2.2830E+03./T) ); 
Gstr{i,   1}='BUTADIENE13';Gstr{i,   2}='O3';
fBUTADIENE13(i)=fBUTADIENE13(i)-1.0;fO3(i)=fO3(i)-1.0;
fO3(i)=fO3(i)+  1.000;fACROLEIN(i)=fACROLEIN(i)+  0.520;

% 325, <T12>
i=i+1;
Rnames{ 325} = 'BUTADIENE13 + NO3 = NO3 +  0.04500*ACROLEIN ';
k(:,i) = (  1.7900E-13 ); 
Gstr{i,   1}='BUTADIENE13';Gstr{i,   2}='NO3';
fBUTADIENE13(i)=fBUTADIENE13(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;fACROLEIN(i)=fACROLEIN(i)+  0.045;

% 326, <TCL3>
i=i+1;
Rnames{ 326} = 'BUTADIENE13 + CL = CL +  0.58000*ACROLEIN ';
k(:,i) = (  2.5100E-10 ); 
Gstr{i,   1}='BUTADIENE13';Gstr{i,   2}='CL';
fBUTADIENE13(i)=fBUTADIENE13(i)-1.0;fCL(i)=fCL(i)-1.0;
fCL(i)=fCL(i)+  1.000;fACROLEIN(i)=fACROLEIN(i)+  0.580;

% 327, <T13>
i=i+1;
Rnames{ 327} = 'ACRO_PRIMARY +  HO =  HO ';
k(:,i) = (  2.0000E-11 ); 
Gstr{i,   1}='ACRO_PRIMARY';Gstr{i,   2}='HO';
fACRO_PRIMARY(i)=fACRO_PRIMARY(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;

% 328, <T14>
i=i+1;
Rnames{ 328} = 'ACRO_PRIMARY + O3 = O3 ';
k(:,i) = (  2.6100E-19 ); 
Gstr{i,   1}='ACRO_PRIMARY';Gstr{i,   2}='O3';
fACRO_PRIMARY(i)=fACRO_PRIMARY(i)-1.0;fO3(i)=fO3(i)-1.0;
fO3(i)=fO3(i)+  1.000;

% 329, <T15>
i=i+1;
Rnames{ 329} = 'ACRO_PRIMARY + NO3 = NO3 ';
k(:,i) = (  1.1500E-15 ); 
Gstr{i,   1}='ACRO_PRIMARY';Gstr{i,   2}='NO3';
fACRO_PRIMARY(i)=fACRO_PRIMARY(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;

% 330, <T16>
i=i+1;
Rnames{ 330} = 'ACRO_PRIMARY =';
k(:,i) = (JACRO_09 ); 
Gstr{i,   1}='ACRO_PRIMARY';
fACRO_PRIMARY(i)=fACRO_PRIMARY(i)-1.0;


% 331, <TCL4>
i=i+1;
Rnames{ 331} = 'ACRO_PRIMARY + CL = CL ';
k(:,i) = (  2.3700E-10 ); 
Gstr{i,   1}='ACRO_PRIMARY';Gstr{i,   2}='CL';
fACRO_PRIMARY(i)=fACRO_PRIMARY(i)-1.0;fCL(i)=fCL(i)-1.0;
fCL(i)=fCL(i)+  1.000;

% 332, <T17>
i=i+1;
Rnames{ 332} = 'ACROLEIN +  HO =  HO ';
k(:,i) = (  2.0000E-11 ); 
Gstr{i,   1}='ACROLEIN';Gstr{i,   2}='HO';
fACROLEIN(i)=fACROLEIN(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;

% 333, <T18>
i=i+1;
Rnames{ 333} = 'ACROLEIN + O3 = O3 ';
k(:,i) = (  2.6100E-19 ); 
Gstr{i,   1}='ACROLEIN';Gstr{i,   2}='O3';
fACROLEIN(i)=fACROLEIN(i)-1.0;fO3(i)=fO3(i)-1.0;
fO3(i)=fO3(i)+  1.000;

% 334, <T19>
i=i+1;
Rnames{ 334} = 'ACROLEIN + NO3 = NO3 ';
k(:,i) = (  1.1500E-15 ); 
Gstr{i,   1}='ACROLEIN';Gstr{i,   2}='NO3';
fACROLEIN(i)=fACROLEIN(i)-1.0;fNO3(i)=fNO3(i)-1.0;
fNO3(i)=fNO3(i)+  1.000;

% 335, <T20>
i=i+1;
Rnames{ 335} = 'ACROLEIN =';
k(:,i) = (JACRO_09 ); 
Gstr{i,   1}='ACROLEIN';
fACROLEIN(i)=fACROLEIN(i)-1.0;


% 336, <TCL5>
i=i+1;
Rnames{ 336} = 'ACROLEIN + CL = CL ';
k(:,i) = (  2.3700E-10 ); 
Gstr{i,   1}='ACROLEIN';Gstr{i,   2}='CL';
fACROLEIN(i)=fACROLEIN(i)-1.0;fCL(i)=fCL(i)-1.0;
fCL(i)=fCL(i)+  1.000;

% 337, <T21>
i=i+1;
Rnames{ 337} = 'TOLU +  HO =  HO ';
k(:,i) = (  1.8000E-12.*exp(  3.4000E+02./T) ); 
Gstr{i,   1}='TOLU';Gstr{i,   2}='HO';
fTOLU(i)=fTOLU(i)-1.0;fHO(i)=fHO(i)-1.0;
fHO(i)=fHO(i)+  1.000;

% 338, <TCL6>
i=i+1;
Rnames{ 338} = 'TOLU + CL = CL ';
k(:,i) = (  6.1000E-11 ); 
Gstr{i,   1}='TOLU';Gstr{i,   2}='CL';
fTOLU(i)=fTOLU(i)-1.0;fCL(i)=fCL(i)-1.0;
fCL(i)=fCL(i)+  1.000;

% 339, <HG1>
i=i+1;
Rnames{ 339} = 'HG + O3 = 0.50000*HGIIAER +  0.50000*HGIIGAS + O3 ';
k(:,i) = (  2.1100E-18.*exp( -1.2565E+03./T) ); 
Gstr{i,   1}='HG';Gstr{i,   2}='O3';
fHG(i)=fHG(i)-1.0;fO3(i)=fO3(i)-1.0;
fHGIIAER(i)=fHGIIAER(i)+  0.500;fHGIIGAS(i)=fHGIIGAS(i)+  0.500;fO3(i)=fO3(i)+  1.000;

% 340, <HG2>
i=i+1;
Rnames{ 340} = 'HG + CL2 = HGIIGAS + CL2 ';
k(:,i) = (  2.6000E-18 ); 
Gstr{i,   1}='HG';Gstr{i,   2}='CL2';
fHG(i)=fHG(i)-1.0;fCL2(i)=fCL2(i)-1.0;
fHGIIGAS(i)=fHGIIGAS(i)+  1.000;fCL2(i)=fCL2(i)+  1.000;

% 341, <HG3>
i=i+1;
Rnames{ 341} = 'HG + H2O2 = HGIIGAS + H2O2 ';
k(:,i) = (  8.5000E-19 ); 
Gstr{i,   1}='HG';Gstr{i,   2}='H2O2';
fHG(i)=fHG(i)-1.0;fH2O2(i)=fH2O2(i)-1.0;
fHGIIGAS(i)=fHGIIGAS(i)+  1.000;fH2O2(i)=fH2O2(i)+  1.000;

% 342, <HG4>
i=i+1;
Rnames{ 342} = 'HG +  HO = 0.50000*HGIIAER +  0.50000*HGIIGAS +  HO ';
k(:,i) = (  7.7000E-14 ); 
Gstr{i,   1}='HG';Gstr{i,   2}='HO';
fHG(i)=fHG(i)-1.0;fHO(i)=fHO(i)-1.0;
fHGIIAER(i)=fHGIIAER(i)+  0.500;fHGIIGAS(i)=fHGIIGAS(i)+  0.500;fHO(i)=fHO(i)+  1.000;

% 343, <HG5>
i=i+1;
Rnames{ 343} = 'HG + CL + M = 0.50000*HG +  0.50000*HGIIGAS + CL ';
k(:,i) = (  2.2500E-33.*exp(  6.8000E+02./T) ).*M; 
Gstr{i,   1}='HG';Gstr{i,   2}='CL';
fHG(i)=fHG(i)-1.0;fCL(i)=fCL(i)-1.0;
fHG(i)=fHG(i)+  0.500;fHGIIGAS(i)=fHGIIGAS(i)+  0.500;fCL(i)=fCL(i)+  1.000;

