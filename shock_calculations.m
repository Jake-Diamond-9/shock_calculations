clear
clc
clf
close all

% Jake Diamond, Johns Hopkins University, jdiamo15@jhu.edu

% This code solves the jump conditions for impact of two known materials
% also calculates the temperature rise due to the shock

% a material library is coded into the program.
% to see what materials are available, run the code 
% and then copy
% "sort(convertCharsToStrings(EOS.Properties.VariableNames.'))" into the
% command window

% create the table with the EOS params
EOS = EOS_table();

%% User input values

% Flyer
flyer_name = 'W';                  % name of flyer material
imp_vels = 1000;                   % impact velocity or velocities
T0_f = 298;                        % reference temperature
h_f = 3e-3;                      % flyer thickness

% Target
target_name = 'SS_304';             % name of target material
T0_t = 298;                     % reference temperature


%% Calculate particle velocity, shock velocity, and pressure

% get EOS params from the table
rho_0f = EOS{"rho0", flyer_name};                          % initial flyer density
C_f = [EOS{"C0", flyer_name}, EOS{"S", flyer_name}];                    % EOS coefficients
gamma0_f = EOS{"gamma0", flyer_name};                              % gruneisen parameter 
Cv_f = EOS{"Cv", flyer_name};                                 % heat capacity at constant volume 

rho_0t = EOS{"rho0", target_name};                          % initial target density
C_t = [EOS{"C0", target_name}, EOS{"S", target_name}];                    % EOS coefficients
gamma0_t = EOS{"gamma0", target_name};                              % gruneisen parameter 
Cv_t = EOS{"Cv", target_name};                                 % heat capacity at constant volume


for j = 1:length(imp_vels)

    U_imp = imp_vels(j);

    % Symbolic variable for particle velocity in the target
    syms U_pt;
    
    % Build the EOS for the target using the user input coefficients.
    % I originally did it this way so I could use higher order hugoniots,
    % like quadratic, cubic, etc. But the temperature calculations are for
    % linear hugoniots only, so I have restricted to EOS material library
    % to only have linear hugioniots. 
    U_st = C_t(1);
    for i = 2:length(C_t)
        U_st = U_st + C_t(i)*U_pt^(i-1);
    end
    
    % Define particle velocity in the flyer as a function of particle velocity
    % in the target
    % See "Dynamic Behavior of Materials" by Meyers, Section 4.3 "Impact"
    U_pf = U_imp-U_pt;
    
    % Build the EOS for the flyer as a function of the particle velocity in the
    % target using the user input coefficients
    U_sf = C_f(1);
    for i = 2:length(C_f)
        U_sf = U_sf + C_f(i)*U_pf^(i-1);
    end
    
    % Symbolic equation to be solved for the particle velocity in the target
    eqn1 = U_pt == (rho_0f*U_pf*U_sf)/(rho_0t*U_st);
    
    % Solve for particle velocity in the target and use that solution to get
    % other useful parameters
    U_pt_sol = double(solve(eqn1,U_pt));
    U_pf_sol = double(subs(U_pf,U_pt,U_pt_sol));
    U_st_sol = double(subs(U_st,U_pt,U_pt_sol));
    U_sf_sol = double(subs(U_sf,U_pt,U_pt_sol));
    P_t_sol  = rho_0t.*U_pt_sol.*U_st_sol;
    P_f_sol  = rho_0f.*U_pf_sol.*U_sf_sol;
    
    % Only take the real solutions in case there are imaginary roots
    U_pt_real = U_pt_sol(imag(U_pt_sol)==0);
    U_pf_real = U_pf_sol(imag(U_pf_sol)==0);
    U_st_real = U_st_sol(imag(U_st_sol)==0);
    U_sf_real = U_sf_sol(imag(U_sf_sol)==0);
    P_t_real = P_t_sol(imag(P_t_sol)==0);
    P_f_real = P_f_sol(imag(P_f_sol)==0);
    
    % Only take the solution in which all parameters have a positive output
    pos_idx = (U_pt_real > 0)&(U_pf_real > 0)&(U_st_real > 0)&(U_sf_real > 0)&(P_t_real > 0)&(P_f_real > 0);
    U_pt = U_pt_real(pos_idx);
    U_pf = U_pf_real(pos_idx);
    U_st = U_st_real(pos_idx);
    U_sf = U_sf_real(pos_idx);
    P_t = P_t_real(pos_idx);
    P_f = P_f_real(pos_idx);
    
    % calculate the shocked density of the flyer and target
    rho_f = (U_sf/(U_sf-U_pf))*rho_0f;    % Davison Eqn. 3.3
    rho_t = (U_st/(U_st-U_pt))*rho_0t;

    % calculate the longitudional green-lagrange strain in the flyer and target
    E11_f = 0.5*((rho_0f/rho_f)^2 - 1);    % Davison Eqns. 2.46 and 2.49
    E11_t = 0.5*((rho_0t/rho_t)^2 - 1);

    % calculate the temperature of the flyer and target
    theta_f = hugoniot_temp(rho_0f, rho_f, C_f(1), C_f(2), gamma0_f, Cv_f, T0_f);
    theta_t = hugoniot_temp(rho_0t, rho_t, C_t(1), C_t(2), gamma0_t, Cv_t, T0_t);

    % Display the solutions in the command window
    disp(' ')
    disp('                          Flyer:   ' + string(flyer_name))
    disp('         Impact velocity  (m/s):   ' + string(U_imp))
    disp('        Shock Drive Time    (s):   ' + string(((2*h_f)/U_sf)))
    disp(' Flyer Particle Velocity  (m/s):   ' + string(U_pf))
    disp('    Flyer Shock Velocity  (m/s):   ' + string(U_sf))
    disp('            Flyer Stress  (GPa):   ' + string(P_f/1e9))
    disp('       Flyer Compression    (-):   ' + string(rho_f/rho_0f))
    disp('            Flyer Strain    (-):   ' + string(E11_f))
    disp('       Flyer Temperature (degC):   ' + string(theta_f - 273))
    disp(' ')
    disp('                         Target:   ' + string(target_name))
    disp('Target Particle Velocity  (m/s):   ' + string(U_pt))
    disp('   Target Shock Velocity  (m/s):   ' + string(U_st))
    disp('           Target Stress  (GPa):   ' + string(P_t/1e9))
    disp('      Target Compression    (-):   ' + string(rho_t/rho_0t))
    disp('           Target Strain    (-):   ' + string(E11_t))
    disp('      Target Temperature (degC):   ' + string(theta_t - 273))
    disp(' ')
    
end

function theta = hugoniot_temp(rho0, rho, C0, S, gamma0, Cv, T0)
    
    % From "Fundamentals of Shock Wave Propagation in Solids" by Davison
    % Free on Springer Link
    % https://doi.org/10.1007/978-3-540-74569-3

    V0 = 1/rho0;
    V = 1/rho;

    % Define needed helper functions
    % Pressure on the Hugoniot with a linear Us-Up relationship
    P_VH = @(V) (rho0.*C0.^2.*((V0 - V)./V0))./(1-S.*((V0 - V)./V0)).^2;    % Davison 3.12
    
    % Slope of the pressure-volume Hugoniot
    dPdVH = @(V) -((C0^2 .* (V0 - S.*V + S.*V0))./(V0 + S.*V - S.*V0).^3);
    
    % Kappa term for the temperature calculation
    kappa_temp = @(V) P_VH(V) + (V0 - V).*dPdVH(V);    % Davison 5.137
    
    % Chi term for the temperature calculation
    chi_temp = @(V) exp((gamma0./V0).*(V0 - V));    % Davison 5.95
    
    % Integrand for the temperature calculation as a single function
    temp_integrand = @(V) kappa_temp(V)./chi_temp(V);
    
    % Function to calculate temperature on the Hugoniot
    temp_H = @(V) chi_temp(V).*(T0 + (1./(2.*Cv)).*integral(temp_integrand, V0, V));    % Davison 5.143
    
    theta = temp_H(V);
end

function EOS = EOS_table()
    
    % Additional values can be found in the following appendix, but they
    % have not been added here. 
    % "Modern Impact and Penetration Mechanics"
    % https://doi.org/10.1017/9781108684026
    % https://doi.org/10.1017/9781108684026.018
     
    rows = ["rho0" "C0" "S" "gamma0" "Cv"].';

    % Kinslow (https://doi.org/10.1016/B978-0-124-08950-1.X5001-0)
    Al_2024 = [2785, 5328, 1.338, 2.0, 838].';
    Cu = [8930, 3940, 1.489, 1.99, 371.8].';
    Be = [1851, 7998, 1.124, 1.16, 1755].';
    Mg = [1740, 4492, 1.263, 1.425, 993.3].';
    Mg_AZ31B = [1780, 4522, 1.242, 1.434, 963.8].';
    Ti = [4528, 5220, 0.767, 1.09, 518.6].';
    Zr = [6505, 3757, 1.018, 1.09, 272.1].';
    Hf = [12885, 2954, 1.121, 0.98, 149.3].';
    V = [6100, 5077, 1.201, 1.29, 479.1].';
    Nb = [8586, 4438, 1.207, 1.470, 265.3].';
    Ta = [16654, 3414, 1.2, 1.6, 137.4].';
    Cr = [7117, 5173, 1.473, 1.19, 444.4].';
    Mo = [10206, 5124, 1.233, 1.52, 245.7].';
    W = [19224, 4029, 1.237, 1.54, 131].';
    Re = [21020, 4184, 1.367, 2.44, 135.6].';
    Rh = [12428, 4807, 1.376, 1.88, 238.8].';
    Pd = [11991, 3948, 1.588, 2.26, 233.4].';
    Ir = [22484, 3916, 1.457, 1.97, 132.8].';
    Pt = [21419, 3598, 1.544, 2.4, 127].';
    SiC = [3120, 8000, 0.95, 1.25, 646.3].';
    WC = [15020, 4920, 1.339, 1.5, 207.7].';
    SS_304 = [7896, 4569, 1.490, 2.17, 441.3].';
    Neoprene = [1439, 2785, 1.419, 1.39, 1649].';
    Polyrubber = [1010, 852, 1.865, 1.5, 1990].';
    LDPE = [915, 2901, 1.481, 1.644, 2300].';
    Polyurethane = [1265, 2486, 1.577, 1.55, 1806].';
    PS = [1044, 2746, 1.319, 1.18, 1345].';
    PMMA = [1185, 2572, 1.536, 0.97, 1469].';
    Paraffin = [900, 2960, 2.531, 1.18, 2890].';
    Epoxy_828 = [1198, 2678, 1.52, 1.13, 1042].';
    PVC = [1681, 1952, 1.66, 0.87, 2090].';
    PTFE = [2147, 1682, 1.819, 0.59, 1048].';

    % Meyers (http://dx.doi.org/10.1002/9780470172278)
    % The heat capacities here are for constant pressure. 
    % Here I assume that the constant pressure heat capacity is close to 
    % the constant volume heat capacity. 
    % Constant volume heat capacity is what is needed for the calculations.
    Ag = [10490, 3230, 1.6, 2.5, 240].';
    Au = [19240, 3060, 1.57, 3.1, 130].';
    Bi = [9840, 1830, 1.47, 1.1, 120].';
    Ca = [1550, 3600, 0.95, 1.1, 660].';
    Cs = [1830, 1050, 1.04, 1.5, 240].';
    Fe = [7850, 3570, 1.92, 1.8, 450].'; % above phase transition
    Hg = [13540, 1490, 2.05, 3.0, 140].';
    K = [860, 1970, 1.18, 1.4, 760].';
    Li = [530, 4650, 1.13, 0.9, 3410].';
    Na = [970, 2580, 1.24, 1.3, 1230].';
    Ni = [8870, 4600, 1.44, 2.0, 440].';
    Pb = [11350, 2050, 1.46, 2.8, 130].';
    Rb = [1530, 1130, 1.27, 1.9, 360].';
    Sn = [7290, 2610, 1.49, 2.3, 220].';
    U = [18950, 2490, 2.20, 2.1, 120].';
    Zn = [7140, 3010, 1.58, 2.1, 390].';
    KCl = [1990, 2150, 1.54, 1.3, 0.68].'; % above phase transition
    LiF = [2640, 5150, 1.35, 2.0, 1500].';
    NaCl = [2160, 3530, 1.34, 1.6, 870].'; % below phase transition
    Al_6061 = [2700, 5350, 1.34, 2.0, 890].';
    Brass = [8450, 3730, 1.43, 2.0, 380].';
    H2O = [1000, 1650, 1.92, 0.1, 4190].';

    % Carter and Marsh and https://doi.org/10.1016/j.ijimpeng.2012.09.001
    Al_1100 = [2714, 5390, 1.34, 1.97, 898.7].';

    % density from manufacturer spec sheet (lexan 8010)
    % hugoniot params from a linear fit of literature values up to 3GPa
    % gruneisen param assume similar to PMMA
    % heat capacity from https://apps.dtic.mil/sti/citations/tr/ADA566369
    PC = [1200, 2232, 1.737, 1.0, 1300].';

    % values for niobium used for Nicolo's paper
    Nb_Nicolo = [8579, 4500, 1.2, 1.47, 265].';

    % Steinberg "Equation of State and Strength Properties of Selected Materials"
    % The heat capacities here are for constant pressure. 
    % Here I assume that the constant pressure heat capacity is close to 
    % the constant volume heat capacity. 
    % Constant volume heat capacity is what is needed for the calculations.
    Cd = [8639, 2480, 1.64, 2.5, 231].';
    Al_1100_O = [2707, 5386, 1.339, 1.97, 884].';
    Al_7075_T6 = [2804, 5200, 1.36, 2.2, 848].';
    Cu_OFHC = [8930, 3940, 1.489, 2.02, 383].';
    Ni_Monel = [8810, 4190, 1.54, 1.95, 411].';
    Steel_Vascomax_250 = [8129, 3980, 1.58, 1.6, 408].';
    Steel_4340_RC_38 = [7810, 4578, 1.33, 1.67, 448].';
    Th = [11700, 2130, 1.278, 1.45, 112].';
    Ti_6Al_4V = [4419, 5130, 1.028, 1.23, 525].';
    LiF_100_singleCrystal = [2638, 5150, 1.35, 1.69, 1560].';

    % Marsh "LASL Shock Hugoniot Data"
    % Does not include gruneisen parameter or heat capacity
    C_diamond = [3191, 7810, 1.43, NaN, NaN].';
    C_foam = [560, 360, 1.22, NaN, NaN].';
    C_graphite = [1000, 790, 1.30, NaN, NaN].';
    Co = [8820, 4770, 1.28, NaN, NaN].';
    Er = [8300, 1760, 1.19, NaN, NaN].';
    In = [7278, 2490, 1.5, NaN, NaN].';
    O2_liquid = [1202, 1880, 1.34, NaN, NaN].';
    Tl = [11840, 1900, 1.48, NaN, NaN].';
    Al_921_T = [2828, 5150, 1.37, NaN, NaN].';
    Steel_347 = [7910, 4620, 1.42, NaN, NaN].';
    NbC = [7500, 5380, 1.46, NaN, NaN].'; % niobium carbide
    Olivine = [3214, 8220, 0.83, NaN, NaN].';
    Periclase_singleCrystal = [3584, 6600, 1.37, NaN, NaN].';
    TaC = [12600, 3320, 1.49, NaN, NaN].'; % tantalum carbide
    Eclogite = [3551, 6260, 1.02, NaN, NaN].';
    Oil_shale = [2192, 3760, 1.15, NaN, NaN].';
    Adiprene = [1094, 2330, 1.54, NaN, NaN].';
    Estane = [1186, 2320, 1.7, NaN, NaN].';
    Kel_F = [2122, 2030, 1.64, NaN, NaN].';
    Silastic_RTV_521 = [1372, 1540, 1.44, NaN, NaN].';
    Birch = [693, 650, 1.44, NaN, NaN].';
    Mahogany = [412, 330, 1.34, NaN, NaN].';
    Pine = [450, 450, 1.33, NaN, NaN].';
    NH3_ammonia = [726, 2000, 1.51, NaN, NaN].';
    Ethylene_glycol = [1112, 2150, 1.55, NaN, NaN].';
    Composition_B = [1715, 3060, 2.01, NaN, NaN].';
    Nitromethane = [1125, 1650, 1.64, NaN, NaN].';
    TNT_liquid = [1472, 2140, 1.57, NaN, NaN].';

    % Davison (https://doi.org/10.1007/978-3-540-74569-3)
    % Does not include gruneisen parameter or heat capacity
    Al2O3 = [3988, 11186, 1.05, NaN, NaN].';


    EOS = table(Al_2024, Cu, Be, Mg, Mg_AZ31B, Ti, Zr, Hf, V, Nb, Ta, Cr, Mo, W, Re, Rh, Pd, Ir, Pt, SiC, WC, SS_304, Neoprene, Polyrubber, LDPE, Polyurethane, PS, PMMA, Paraffin, Epoxy_828, PVC, PTFE,...
        Ag, Au, Bi, Ca, Cs, Fe, Hg, K, Li, Na, Ni, Pb, Rb, Sn, U, Zn, KCl, LiF, NaCl, Al_6061, Brass, H2O, Al_1100, PC, Nb_Nicolo,...
        Cd, Al_1100_O, Al_7075_T6, Cu_OFHC, Ni_Monel, Steel_Vascomax_250, Steel_4340_RC_38, Th, Ti_6Al_4V, LiF_100_singleCrystal,...
        C_diamond, C_foam, C_graphite,  Co, Er, In, O2_liquid, Tl, Al_921_T, Steel_347, NbC, Olivine, Periclase_singleCrystal, TaC, Eclogite, Oil_shale, Adiprene,...
        Estane, Kel_F, Silastic_RTV_521, Birch, Mahogany, Pine, NH3_ammonia, Ethylene_glycol, Composition_B, Nitromethane, TNT_liquid, Al2O3, 'RowNames',rows);

end


