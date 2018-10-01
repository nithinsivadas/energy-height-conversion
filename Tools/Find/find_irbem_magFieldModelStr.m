function magFieldModelStr=find_irbem_magFieldModelStr(magFieldModel)
%% find_irbem_magFieldModelStr.m Generates the name of the HDF5 data-structure
%                                of the magnetic field model specified by
%                                the input integer. 
%
% magFieldModel: Flat specifying which magnetic field to use. Default:11
%                 0   = no external field 
%                 1   = Mead & Fairfield [1975] 
%                     (uses 0?Kp?9 - Valid for rGEO?17. Re) 
%                 2   = Tsyganenko short [1987] 
%                     (uses 0?Kp?9 - Valid for rGEO?30. Re) 
%                 3   = Tsyganenko long [1987] 
%                     (uses 0?Kp?9 - Valid for rGEO?70. Re) 
%                 4   = Tsyganenko [1989c] 
%                     (uses 0?Kp?9 - Valid for rGEO?70. Re) 
%                 5   = Olson & Pfitzer quiet [1977] 
%                     (default - Valid for rGEO?15. Re) 
%                 6   = Olson & Pfitzer dynamic [1988] 
%                     (uses 5.?dens?50., 300.?velo?500., 
%                     -100.?Dst?20. - Valid for rGEO?60. Re) 
%                 7   = Tsyganenko [1996] 
%                     (uses -100.?Dst (nT)?20., 0.5?Pdyn (nPa)?10., 
%                     |ByIMF| (nT)?10., |BzIMF| (nT)?10. - Valid for rGEO?40. Re) 
%                 8   = Ostapenko & Maltsev [1997] 
%                     (uses dst,Pdyn,BzIMF, Kp) 
%                 9   = Tsyganenko [2001] 
%                     (uses -50.?Dst (nT)?20., 0.5?Pdyn (nPa)?5., |ByIMF| (nT)?5., 
%                     |BzIMF| (nT)?5., 0.?G1?10., 0.?G2?10. - Valid for xGSM?-15. Re) 
%                 10 =Tsyganenko [2001] storm  
%                     (uses Dst, Pdyn, ByIMF, BzIMF, G2, G3 
%                     - there is no upper or lower limit 
%                     for those inputs - Valid for xGSM?-15. Re) 
%                 11 =Tsyganenko [2004] storm  
%                     (uses Dst, Pdyn, ByIMF, BzIMF, 
%                     W1, W2, W3, W4, W5, W6 
%                     - there is no upper or lower limit for those inputs 
%                     - Valid for xGSM?-15. Re) 
%                 12 =Alexeev [2000] 
%                     - also known as Paraboloid model - 
%                     Submitted to ISO  (uses Dens, velo, Dst, BzIMF, AL)

if nargin<1
    magFieldModel = 11;
end

    switch magFieldModel 
        case 0
            magFieldModelStr = 'NoExternalField';
        case 1 
            magFieldModelStr = 'MF75';
        case 2
            magFieldModelStr = 'TS87short';
        case 3
            magFieldModelStr = 'TS87long';
        case 4
            magFieldModelStr = 'TS89';
        case 5
            magFieldModelStr = 'OP77quiet';
        case 6
            magFieldModelStr = 'OP88dynamic';
        case 7
            magFieldModelStr = 'TS96';
        case 8
            magFieldModelStr = 'OM97';
        case 9
            magFieldModelStr = 'TS01';
        case 10
            magFieldModelStr = 'TS01storm';
        case 11
            magFieldModelStr = 'TS04storm';
        case 12
            magFieldModelStr = 'Alexeev2000'; 
    end
    
end