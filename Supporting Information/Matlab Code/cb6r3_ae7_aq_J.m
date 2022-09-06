function J = CB6R3_AE7_AQ_J(Met,Jmethod)
% Calculates photolysis frequencies for the cb6r3_ae7_aq mechanism in the CMAQ model
% Met: structure containing required meteorological constraints. Required vars depend on Jmethod.
%       Met.SZA: solar zenith angle in degrees
%       Met.ALT: altitude, meters
%       Met.O3col: overhead ozone column, DU
%       Met.albedo: surface reflectance, 0-1 (unitless)
%       Met.T: temperature, T
%       Met.P: pressure, mbar
%       Met.LFlux: name of a text file containing an actinic flux spectrum
%
% Jmethod: numeric flag or string specifying how to calculate J-values. Default is 'MCM'.
%       0 or 'MCM':      use MCMv3.3.1 parameterization.
%                         Some reactions are not included in MCM. For these, 'HYBRID' values are used.
%                         Required Met fields: SZA
%       1 or 'BOTTOMUP': bottom-up integration of cross sections/quantum yields.
%                         See J_BottomUp.m for more info.
%                         Required Met fields: LFlux, T, P
%       2 or 'HYBRID':   Interpolation of hybrid J-values from TUV solar spectra.
%                         See J_TUVhybrid.m for more info.
%                         Required Met fields: SZA, ALT, O3col, albedo
%
% OUTPUTS:
% J: structure of J-values.
%
% INPUTS
struct2var(Met)

if nargin<2
    Jmethod = 'MCM';
elseif ischar(Jmethod)
    Jmethod = upper(Jmethod);
end

% J-Values
switch Jmethod
    case {0,'MCM'}
        error(['MCM option not functional for cb6r3_ae7_aq mechanism.'])

    case {1,'BOTTOMUP'}
        Jmcm = J_BottomUp(LFlux,T,P);

    case {2,'HYBRID'}
        Jmcm = J_Hybrid(SZA,ALT,O3col,albedo);

    otherwise
        fprintf('Jmethod = %f\n',Jmethod);
        error(['MCMv331_J: invalid Jmethod option selected'])

end
%rename
J=struct;
J.JNO2_IUPAC10        = Jmcm.J4;
J.JNO3NO_06     = Jmcm.J5;
J.JNO3NO2_06    = Jmcm.J6;
J.JHNO3_IUPAC10       = Jmcm.J8;
J.JO3_O3P_IUPAC10        = Jmcm.J2;
J.JO3_O1D_IUPAC10       = Jmcm.J1;
J.JHONO_IUPAC10       = Jmcm.J7;
J.JH2O2_IUPAC10       = Jmcm.J3;
J.JPNA_IUPAC10 = Jmcm.J3;
J.JMEPX_IUPAC10 = Jmcm.J41;
J.JNTR_IUPAC10       = Jmcm.J54;
J.JALD2_R_IUPAC10 = Jmcm.J13;
J.JALDX_R_IUPAC10 = Jmcm.J14;
J.JFORM_R_IUPAC10 = Jmcm.J11;
J.JFORM_M_IUPAC10     = Jmcm.J12;
J.JMGLY_IUPAC10       = Jmcm.J34;
J.JHPALD       = Jmcm.J20;
J.JGLY_R_IUPAC10        = Jmcm.J31 + Jmcm.J32 + Jmcm.J33;
J.JACET_IUPAC10       = Jmcm.J21;
J.JMEK        = Jmcm.J22;
J.JISPD       = Jmcm.J23 + Jmcm.J24;

% NO DIRECT MCM ANALOGUES
J.JN2O5_IUPAC10       = Jmcm.Jn19 + Jmcm.Jn20;
J.JPAN_IUPAC10        = Jmcm.Jn14 + Jmcm.Jn15;
J.JACRO_09       = Jmcm.Jn11;
J.JGLYD_IUPAC10       = Jmcm.Jn9;
J.JKET_IUPAC10 = Jmcm.Jn3;
