function R = canard_roll_envelope(Pover, mode)
%CANARD_ROLL_ENVELOPE  Effect of canard roll control vs altitude and velocity.
%
%   R = canard_roll_envelope()            run with the YJSP vehicle (default_params)
%   R = canard_roll_envelope(Pover)       override any parameter, e.g.
%       canard_roll_envelope(struct('thrust',4000,'eta',0.2))
%   F = canard_roll_envelope(Pover, 'tools')   return the parameter set and
%       handles to the model functions (aero_roll, fly_trajectory, ...)
%       without running anything; used by roll_gain_design.m
%
%   PART 1  Frozen-coefficient envelope maps over altitude x velocity:
%           dynamic pressure, canard roll authority, uncontrolled spin
%           rate, pitch natural frequency, resonance amplification, canard
%           trim angle, and the closed-loop response of the PI loop
%           (gain-scheduled vs fixed gains) to a fin-misalignment step.
%   PART 2  A 2-D point-mass trajectory flown three ways (canards locked,
%           scheduled PI, fixed-gain PI); its path is drawn on the maps.
%
%   Roll plant (from the 6-DOF derivation):   pdot = -a*p + b*delta + dL
%       a  roll damping            [1/s]      fins + canards, strip theory
%       b  canard effectiveness    [1/s^2]    eta * canard roll moment / Ix
%       dL misalignment forcing    [rad/s^2]  fin cant, all fins same sense
%   Pitch/yaw (non-rolling frame): wn^2 = -Cma*qbar*S*d/Iy, damping from
%   CNa and Cmq. Roll-pitch resonance amplification:
%       |A|/eps = wn^2 / sqrt((wn^2-p^2)^2 + (2*zeta*wn*p)^2)  -> 1/(2 zeta)
%
%   Controller: PI on roll rate, pole placement s^2 + 2 zc wc s + wc^2.
%   Scheduled: gains recomputed from (h, V) with wc = clamp(wc_ratio*a,
%   wc, wc_max). Fixed: gains frozen at the design point (ref_h, ref_V).
%   Phase margin includes servo lag, half-sample hold and sensor delay.
%   The map's misalignment STEP is a worst case: in flight dL builds up
%   smoothly with q, so trajectory results (Part 2) are more realistic.
%
%   default_params() holds the YJSP pressure-fed vehicle (CoDR Candidate 3,
%   fins and triangular canard updated 11/12); each number is tagged with
%   its source file, and ASSUMED marks what is not in the data yet (servo,
%   canard planform orientation). Override any field through Pover.
%
%   Model limits: flat Earth, ISA to 84 km (geopotential altitude), unswept
%   fins, Barrowman subsonic lift slope, supersonic linear theory with tip
%   loss, linear blend across Mach 0.8-1.2, constant eta. Map points hold
%   altitude, velocity and mass fixed (frozen coefficients); the
%   misalignment is applied as a step at t = 0 (conservative vs. a smooth
%   build-up).
%
%   Runs in MATLAB (no toolboxes needed) and GNU Octave.

P = default_params();
if nargin >= 1 && ~isempty(Pover)
    fn = fieldnames(Pover);
    for i = 1:numel(fn), P.(fn{i}) = Pover.(fn{i}); end
end
P = derive_geometry(P);
if nargin >= 2 && strcmpi(mode, 'tools')
    R = struct('P', P, 'mass_props', @mass_props, 'aero_roll', @aero_roll, ...
               'isa_atm', @isa_atm, 'pi_gains', @pi_gains, 'phase_margin', @phase_margin, ...
               'fly_trajectory', @fly_trajectory, 'roll_sim', @roll_sim);
    return
end
rng(1);                                   % repeatable gyro noise

%% ===================== PART 2: trajectory (needed for overlays) ========
T   = fly_trajectory(P);
N   = numel(T.t);
mpT = mass_props(P, T.f);
Etr = aero_roll(P, T.h, T.V, mpT, P.eta);         % true plant
Emd = aero_roll(P, T.h, T.V, mpT, P.eta_model);   % controller's model
on  = T.offrail & (T.qbar >= P.qbar_min);

[KpS, KiS] = pi_gains(P, Emd.a, Emd.b);           % scheduled gains
mpRef = mass_props(P, P.map_f);                   % fixed-gain design point
Eref  = aero_roll(P, P.ref_h, P.ref_V, mpRef, P.eta_model);
[KpF, KiF] = pi_gains(P, Eref.a, Eref.b);

bT = Etr.b(:).';  dLT = Etr.dL(:).';
bT(~T.offrail) = 0;  dLT(~T.offrail) = 0;         % rail buttons stop roll
C.a  = repmat(Etr.a(:).', 3, 1);
C.b  = repmat(bT, 3, 1);
C.dL = repmat(dLT, 3, 1);
C.Kp = [zeros(1,N); KpS(:).'; KpF*ones(1,N)];
C.Ki = [zeros(1,N); KiS(:).'; KiF*ones(1,N)];
C.on = [false(1,N); on(:).'; on(:).'];
Str  = roll_sim(P, C, N, true);                   % rows: locked, sched, fixed

wnT  = Etr.wn(:).';  zT = Etr.zeta(:).';
ampT = zeros(3, N);
for i = 1:3, ampT(i,:) = amp_factor(abs(Str.p(i,:)), wnT, zT); end
ks   = find(on(:).');  ks = ks(1:20:end);          % phase margin along flight
pmTS = phase_margin(P, Etr.a(ks), Etr.b(ks), KpS(ks), KiS(ks));
pmTF = phase_margin(P, Etr.a(ks), Etr.b(ks), KpF, KiF);

%% ===================== PART 1: envelope maps ===========================
hv = linspace(0, P.map_hmax, 61);
Vv = linspace(P.map_Vmin, P.map_Vmax, 67);
[VV, HH] = meshgrid(Vv, hv);
mpM = mass_props(P, P.map_f);
E   = aero_roll(P, HH, VV, mpM, P.eta);
Em  = aero_roll(P, HH, VV, mpM, P.eta_model);

pss    = E.dL ./ E.a;                    % uncontrolled steady spin [rad/s]
dtrim  = E.dL ./ E.b * 180/pi;           % canard angle to cancel it [deg]
auth   = E.b * P.dmax;                   % max roll acceleration [rad/s^2]
ampmax = 1 ./ (2*E.zeta);                % amplification at p = wn (peak ~ this for small zeta)
ctrl   = E.qbar >= P.qbar_min;           % controller enabled

n  = numel(HH);  Nst = round(P.t_step/P.dt);
[KpM, KiM] = pi_gains(P, Em.a, Em.b);
Cm = struct('a', E.a(:), 'b', E.b(:), 'dL', E.dL(:), ...
            'Kp', KpM(:), 'Ki', KiM(:), 'on', ctrl(:));
Ms = roll_sim(P, Cm, Nst, false);        % scheduled gains
Cm.Kp = KpF*ones(n,1);  Cm.Ki = KiF*ones(n,1);
Mf = roll_sim(P, Cm, Nst, false);        % fixed gains

pkS  = reshape(Ms.pk, size(HH));       pkS(~ctrl) = NaN;
pkF  = reshape(Mf.pk, size(HH));       pkF(~ctrl) = NaN;
satS = reshape(Ms.satfrac, size(HH));  satS(~ctrl) = NaN;
% phase margin with servo lag, half-sample hold and sensor delay
pmS  = reshape(phase_margin(P, E.a, E.b, KpM, KiM), size(HH));  pmS(~ctrl) = NaN;
pmF  = reshape(phase_margin(P, E.a, E.b, KpF, KiF), size(HH));  pmF(~ctrl) = NaN;

%% ===================== Figures =========================================
f1 = figure('Name', 'Envelope: physics', 'Color', 'w', 'Position', [40 40 1350 760]);
subplot(2,3,1); map_panel(VV, HH, E.qbar/1e3, E.M, T, 'Dynamic pressure', 'log_{10} q [kPa]', true);
subplot(2,3,2); map_panel(VV, HH, auth, E.M, T, 'Canard roll authority b\delta_{max}', 'log_{10} [rad/s^2]', true);
subplot(2,3,3); map_panel(VV, HH, pss, E.M, T, 'Uncontrolled spin p_{ss} (red: p_{ss} = \omega_n)', ...
                          'log_{10} [rad/s]', true, pss./E.wn, 1);
subplot(2,3,4); map_panel(VV, HH, E.wn, E.M, T, 'Pitch/yaw natural frequency \omega_n', '[rad/s]', false);
subplot(2,3,5); map_panel(VV, HH, ampmax, E.M, T, 'Resonance amplification 1/(2\zeta)', 'log_{10} [-]', true);
subplot(2,3,6); map_panel(VV, HH, dtrim, E.M, T, 'Canard trim angle (red: \delta_{max})', '[deg]', false, ...
                          dtrim, P.dmax*180/pi);

f2 = figure('Name', 'Envelope: controller', 'Color', 'w', 'Position', [60 60 1350 760]);
subplot(2,3,1); map_panel(VV, HH, pkS, E.M, T, 'Scheduled PI: peak |p|, misalignment step', 'log_{10} [rad/s]', true);
subplot(2,3,2); map_panel(VV, HH, pkS./E.wn, E.M, T, 'Scheduled PI: peak|p| / \omega_n (red: 1)', ...
                          'log_{10} [-]', true, pkS./E.wn, 1);
subplot(2,3,3); map_panel(VV, HH, pmS, E.M, T, 'Scheduled PI: phase margin (red: 30 deg)', '[deg]', false, pmS, 30);
subplot(2,3,4); map_panel(VV, HH, pkF, E.M, T, 'Fixed gains: peak |p|, misalignment step', 'log_{10} [rad/s]', true);
subplot(2,3,5); map_panel(VV, HH, pkF./E.wn, E.M, T, 'Fixed gains: peak|p| / \omega_n (red: 1)', ...
                          'log_{10} [-]', true, pkF./E.wn, 1);
subplot(2,3,6); map_panel(VV, HH, pmF, E.M, T, 'Fixed gains: phase margin (red: 30 deg)', '[deg]', false, pmF, 30);

tt = T.t;  lw = 1.3;  cs = {'k', 'b', [0.85 0.4 0]};
nm = {'canards locked', 'scheduled PI', 'fixed-gain PI'};
f3 = figure('Name', 'Trajectory: time histories', 'Color', 'w', 'Position', [80 80 1200 760]);
subplot(2,2,1);
plot(tt, T.h/1e3, 'k', tt, T.M, 'b', tt, T.qbar/1e4, 'r', 'LineWidth', lw); grid on
burn_line(P); xlabel('t [s]'); title('Flight');
legend('altitude [km]', 'Mach', 'q [10 kPa]', 'Location', 'best');
subplot(2,2,2);  rate_plot(tt, Str.p, wnT, cs, nm, lw);  burn_line(P);
xlabel('t [s]');  title('Roll rate vs pitch frequency');
subplot(2,2,3);
hl = plot(tt, Str.d(2,:)*180/pi, 'Color', cs{2}, 'LineWidth', lw);  hold on
hl = [hl, plot(tt, Str.d(3,:)*180/pi, 'Color', cs{3}, 'LineWidth', lw)];
plot(tt([1 end]), P.dmax*180/pi*[1 1], 'r:', tt([1 end]), -P.dmax*180/pi*[1 1], 'r:');
grid on; burn_line(P); xlabel('t [s]'); ylabel('\delta [deg]'); title('Canard angle');
legend(hl, nm(2:3), 'Location', 'best');
subplot(2,2,4);
hl = [];
for i = 1:3, hl = [hl, semilogy(tt, ampT(i,:), 'Color', cs{i}, 'LineWidth', lw)]; hold on; end
grid on; burn_line(P); xlabel('t [s]'); ylabel('|A|/\epsilon');
title('Roll-pitch resonance amplification'); legend(hl, nm, 'Location', 'best');

f4 = figure('Name', 'Trajectory: effect vs altitude and velocity', 'Color', 'w', 'Position', [100 100 1150 450]);
subplot(1,2,1);  rate_plot(T.h/1e3, Str.p, wnT, cs, nm, lw);
xlabel('Altitude [km]');  title('Roll rate vs altitude (boost + coast)');
subplot(1,2,2);  rate_plot(T.V, Str.p, wnT, cs, nm, lw);
xlabel('Velocity [m/s]');  title('Roll rate vs velocity (boost + coast)');

if P.save_figs
    figs = [f1 f2 f3 f4];
    for i = 1:4, print(figs(i), '-dpng', '-r110', sprintf('canard_fig%d.png', i)); end
end

%% ===================== Summary =========================================
dt = P.dt;  [~, ib] = min(abs(T.t - P.t_burn));
fprintf('\n--- Trajectory ---\n');
fprintf('Burnout: V = %.0f m/s, Mach %.2f, h = %.1f km\n', T.V(ib), T.M(ib), T.h(ib)/1e3);
fprintf('Apogee: %.1f km at t = %.1f s, max q = %.0f kPa\n', T.h(end)/1e3, T.t(end), max(T.qbar)/1e3);
fprintf('Fixed gains designed at h = %.0f m, V = %.0f m/s: Kp = %.4f s, Ki = %.3f\n', ...
        P.ref_h, P.ref_V, KpF, KiF);
fprintf('Phase margin while control is on: scheduled %.0f-%.0f deg, fixed %.0f-%.0f deg\n', ...
        min(pmTS), max(pmTS), min(pmTF), max(pmTF));
% statistics over the aerodynamically relevant part of flight (q >= qbar_min);
% near apogee q -> 0 and wn -> 0, so |A|/eps there has no structural meaning
aero = T.qbar(:).' >= P.qbar_min;
fprintf('Statistics for q >= %.0f kPa (%.1f s of flight):\n', P.qbar_min/1e3, sum(aero)*dt);
fprintf('%-16s %10s %12s %14s %12s\n', 'case', 'peak|p|', 'max|A|/eps', 't near wn [s]', 'max|delta|');
for i = 1:3
    pr  = abs(Str.p(i,:));
    nr  = sum(aero & pr > 0.8*wnT & pr < 1.2*wnT) * dt;
    fprintf('%-16s %8.2f r/s %12.2f %14.2f %9.1f deg\n', nm{i}, max(pr(aero)), max(ampT(i,aero)), ...
            nr, max(abs(Str.d(i,:)))*180/pi);
end
dS   = Str.d(2,:);
nsat = sum(aero & abs(dS) >= P.dmax - 1e-4) * dt;
Hfl  = Etr.Hd(:).' .* dS;                               % scheduled PI, in flight
fprintf('Scheduled PI: canard at travel limit for %.1f s (limit %.0f deg)\n', nsat, P.dmax*180/pi);
fprintf('Hinge moment per canard (hinge at %.0f%% MAC): max in flight %.1f N m, at delta_max & max q %.1f N m\n', ...
        100*P.hinge_frac, max(abs(Hfl)), max(abs(Etr.Hd))*P.dmax);

R = struct('P', P, 'T', T, 'traj', Str, 'wn_traj', wnT, 'amp_traj', ampT, ...
           'hinge_traj', Hfl, 'sat_time', nsat, ...
           'map', struct('V', VV, 'h', HH, 'E', E, 'pss', pss, 'dtrim', dtrim, ...
           'auth', auth, 'ampmax', ampmax, 'peak_sched', pkS, 'peak_fixed', pkF, ...
           'sat_sched', satS, 'pm_sched', pmS, 'pm_fixed', pmF));
end

%% =======================================================================
function P = default_params()
% YJSP pressure-fed vehicle (CoDR "Pressure Fed Candidate 3"). Sources:
%   [Y]  CODR_pressure_fed.yaml / fin_can block of 11/12 (fins, boattail)
%   [F]  flight_pressure_fed_regulator.csv (mass, CG, inertia, thrust, drag)
%   [C]  Canard_Calculations_Updated_version_1Fin.csv (triangular canard)
%   [T]  canard_roll_force_1.mat (canard force table)
%   ASSUMED = not in the data yet.
% ---- body and aerodynamic surfaces --------------------------------------
P.d       = 0.31115;            % body diameter [m]                         [Y]
P.canard  = [0.08652 0.300 0];  % exposed span, root, tip chord [m]: triangular [C]
P.canard_ycp = 0.0306;          % spanwise CP from the root [m] (empty -> from planform) [C]
P.Nc      = 2;                  % number of canards (opposite each other)
P.Nc_pitch = 2;                 % canards lifting in the analysed pitch plane
                                %   (both lie in one plane: worst case is that plane)
P.fin     = [0.3485911021958517 0.49238255144153215 0.19119283595291195];
                                % exposed span, root, tip chord [m]           [Y]
                                % (sweep_fraction 0.69 and the boattail taper,
                                %  aft dia 0.2773 m over 0.5183 m, are not modelled)
P.Nf      = 4;                  % number of fins                              [Y]
P.Kc_roll = 1.494;              % canard-body interference factor in roll: the
                                %   body_factor of the team's force table     [T]
                                %   (with it the model's canard slope is within
                                %    ~+-20% of that table over Mach 0.3-5)
P.hinge_frac = 0.175;           % hinge line, fraction of MAC aft of the MAC leading
                                %   edge: hinge 135 mm aft of the root LE [C], MAC =
                                %   2/3 c_r = 200 mm starting 100 mm aft of the root LE
                                %   -- ASSUMES a swept leading edge, straight trailing edge
P.eta       = 0.5;              % canard-tail interaction (net_roll_multiplier) [Y]
P.eta_model = 0.5;              % value assumed by the gain schedule
P.fin_cant  = 0.00174533;       % fin misalignment used by the sizer, 0.1 deg [rad] [Y]
% ---- component CP stations from nose tip [m] ---------------------------
P.x_nose   = 0.466*1.855;       % von Karman nose, fineness 5.963 -> 1.855 m  [Y]
                                %   (0.466 L is the ogive value: approximate)
P.x_canard = 6.372;             % canard CP station: where the sizer applies the
                                %   canard normal force (pitch moment / force) [F]
P.x_fin    = 9.515 - 0.518 + 0.236;  % fin CP: root LE at the front of the boattail +
                                %   Barrowman offset (ESTIMATE; the sizer's body+fin CP
                                %   implies ~9.25 m). Used only for pitch damping:
                                %   the whole-vehicle CP comes from P.aero.xcp below
% ---- mass properties -----------------------------------------------------
P.m_dry   = 217.11;  P.m_prop  = 302.71;   % burnout mass, mass burned [kg]   [F]
P.x_dry   = 4.9141;  P.x_prop  = 5.6864;   % CG stations [m]: 5.364 m full -> 4.914 m empty [F]
P.Ix_dry  = 4.7529;  P.Ix_prop = 7.3267;   % roll inertia 12.08 -> 4.75 kg m^2 [F]
P.Iy_full = 2930.2;  P.Iy_empty = 2116.7;  % pitch inertia [kg m^2]           [F]
% ---- motor and launch ----------------------------------------------------
P.thrust = 26060;  P.t_burn = 26.56;    % total impulse 692.1 kN s / burn time [N], [s]
P.launch_angle = 90*pi/180;  P.L_rail = 15;                                 % [Y]
P.h0 = 1191;                    % launch site altitude [m]                    [Y]
% ---- servo (ASSUMED until the actuator is selected and measured) ---------
P.tau_s     = 0.020;            % small-signal time constant [s]
P.servo_delay = 0.003;          % command-to-motion dead time [s]
                                %   (half PWM frame + pulse width + processing)
P.dmax      = 10*pi/180;        % travel limit [rad]                          [Y]
P.dmax_rate = 250*pi/180;       % rate limit under load [rad/s]
P.backlash  = 0;                % total free play servo -> canard [rad] (0 = none)
P.pos_sensor = 'servo';         % canard angle feedback: 'servo' (before the
                                %   free play, e.g. servo pot) or 'canard'
% ---- controller and sensors ---------------------------------------------
P.f_loop   = 960;               % control loop rate [Hz] (7680 Hz IMU / 8)
P.wc       = 15;  P.zc = 0.7;   % design closed-loop natural frequency [rad/s] / damping
P.wc_ratio = 2;                 % schedule: wc >= wc_ratio * a (beat open-loop pole)
P.wc_max   = 25;                % ... but no faster than servo + delay allow
P.p_ref    = 0;                 % target roll rate [rad/s]
P.qbar_min = 2000;              % control enabled above this q [Pa]
P.Tt_min   = 0.05;              % anti-windup tracking time floor [s]
P.antiwindup = 'measured';      % 'measured' (delta_m - v) or 'command' (sat(v) - v)
P.sensor_delay = 0.002;         % gyro filter + read latency [s]
P.gyro_range   = 2000*pi/180;   % [rad/s]
P.gyro_bias    = 0.5*pi/180;    % [rad/s]
P.gyro_noise   = 0.05*pi/180;   % rms per sample [rad/s]
P.ref_h = 1500;  P.ref_V = 200; % fixed-gain design point
% ---- simulation and maps -------------------------------------------------
P.dt = 1e-3;  P.t_max = 300;  P.t_step = 1.5;
P.map_hmax = 60e3;  P.map_Vmin = 50;  P.map_Vmax = 1800;
P.map_f = 0;                    % propellant fraction on maps (0 = burnout)
P.save_figs = false;
% ---- aero tables from the sizer's flight output (see aero_roll) -----------
% CD : whole-vehicle drag coefficient, power on, times the sizer's 1.12 drag
%      factor (multiplier + canard drag)                                    [F]
% xcp: whole-vehicle centre of pressure, canards included (derived.total_cp_m)
%      [F]. NOTE: that run had the previous fins (span 0.298 m); the 11/12
%      fins (span 0.349 m) move the CP slightly aft -> refresh after a rerun.
P.aero.M   = [0.06 0.43 0.84 1.05 1.28 1.77 2.35 3.01 3.64 4.38 5.18];
P.aero.CD  = [0.425 0.425 0.414 0.468 0.450 0.352 0.253 0.210 0.183 0.159 0.143]*1.12;
P.aero.xcp = [7.639 7.755 7.757 7.896 7.815 7.494 7.060 6.589 6.283 6.172 5.552];
end

function P = derive_geometry(P)
r = P.d/2;
if isempty(P.Nc_pitch), P.Nc_pitch = P.Nc/2; end
[s, cr, ct] = deal(P.canard(1), P.canard(2), P.canard(3));
P.Sc    = 0.5*(cr + ct)*s;                      % canard panel area
P.AR_c  = 2*s^2/P.Sc;
P.rc    = r + s/3*(cr + 2*ct)/(cr + ct);        % spanwise CP radius
if isfield(P, 'canard_ycp') && ~isempty(P.canard_ycp), P.rc = r + P.canard_ycp; end
P.Kc_fb = 1 + r/(s + r);                        % Barrowman interference
P.mac_c = 2/3*(cr + ct - cr*ct/(cr + ct));      % canard mean aerodynamic chord
[s, cr, ct] = deal(P.fin(1), P.fin(2), P.fin(3));
P.Sf    = 0.5*(cr + ct)*s;
P.AR_f  = 2*s^2/P.Sf;
P.rf    = r + s/3*(cr + 2*ct)/(cr + ct);
P.Kf_fb = 1 + r/(s + r);
end

function mp = mass_props(P, f)          % f = propellant fraction remaining
mp.m   = P.m_dry + f*P.m_prop;
mp.xcg = (P.m_dry*P.x_dry + f*P.m_prop*P.x_prop) ./ mp.m;
mp.Ix  = P.Ix_dry + f*P.Ix_prop;
mp.Iy  = P.Iy_empty + f*(P.Iy_full - P.Iy_empty);
end

function E = aero_roll(P, h, V, mp, eta)
[rho, asnd] = isa_atm(h);
Vs     = max(V, 1);
E.M    = V ./ asnd;
E.qbar = 0.5 * rho .* V.^2;
CNc = lift_slope(E.M, P.AR_c);          % per canard panel [1/rad]
CNf = lift_slope(E.M, P.AR_f);          % per fin panel   [1/rad]
% ---- roll channel
E.b  = eta .* P.Nc .* P.Kc_roll .* E.qbar .* P.Sc .* CNc .* P.rc ./ mp.Ix;
E.a  = 0.5 * rho .* Vs .* (P.Nf*P.Sf*P.rf^2*CNf + P.Nc*P.Sc*P.rc^2*CNc) ./ mp.Ix;
E.a  = max(E.a, 1e-6);
E.dL = P.Nf * P.Sf * P.rf * P.fin_cant * E.qbar .* CNf ./ mp.Ix;
% ---- hinge moment per canard per radian of deflection [N m/rad]
% CP moves from 25% MAC (M <= 0.8) to 50% MAC (M >= 1.2); + = CP aft of
% hinge (restoring), - = CP ahead of hinge (overbalanced)
xcp  = 0.25 + 0.25*min(max((E.M - 0.8)/0.4, 0), 1);
E.Hd = E.qbar .* P.Sc .* CNc .* (xcp - P.hinge_frac) * P.mac_c;
% ---- pitch plane, component build-up
S    = pi*P.d^2/4;
CNn  = 2*ones(size(E.M));                               % nose
CNcb = P.Kc_fb*P.Nc_pitch*(P.Sc/S)*CNc;                 % canards lifting in this plane
CNfb = P.Kf_fb*(P.Nf/2)*(P.Sf/S)*CNf;                   % fins in plane
CN   = CNn + CNcb + CNfb;
xcp  = (CNn*P.x_nose + CNcb*P.x_canard + CNfb*P.x_fin) ./ CN;
Cmq  = -2*(CNn.*(P.x_nose - mp.xcg).^2 + CNcb.*(P.x_canard - mp.xcg).^2 ...
         + CNfb.*(P.x_fin - mp.xcg).^2) / P.d^2;
% ---- measured / CFD / RASAero values override the geometry estimates
%   P.aero.M     Mach breakpoints (ascending); every other field is a vector
%                of the same length and is optional. Reference area
%                S = pi d^2/4, reference length d, all angles in radians.
%   .Cld   canard roll-control derivative per rad of canard deflection,
%          canards alone (tail interaction is applied through eta; if the
%          value already includes it, e.g. full-configuration CFD, set eta = 1)
%   .Clp   roll-damping derivative per (p d / 2V)            (negative)
%   .Cl0   roll moment coefficient from fin misalignment at the
%          expected misalignment (replaces fin_cant)
%   .CNa   whole-vehicle normal-force slope [1/rad]
%   .xcp   whole-vehicle CP station from nose tip [m], positive aft
%   .Cmq   pitch-damping derivative per (q d / 2V) about the CG (negative)
%   .CD    drag coefficient (used by fly_trajectory)
%   .xcp_canard  canard CP as fraction of its MAC (hinge moment)
%   Outside the table the end values are held.
if isfield(P, 'aero') && ~isempty(P.aero)
    A  = P.aero;  Mq = min(max(E.M, A.M(1)), A.M(end));
    tb = @(f) reshape(interp1(A.M(:), A.(f)(:), Mq(:), 'linear'), size(E.M));
    if isfield(A, 'Clp'), E.a  = max(-rho.*Vs*S*P.d^2.*tb('Clp')./(4*mp.Ix), 1e-6); end
    if isfield(A, 'Cld'), E.b  = eta .* E.qbar*S*P.d .* tb('Cld') ./ mp.Ix; end
    if isfield(A, 'Cl0'), E.dL = E.qbar*S*P.d .* tb('Cl0') ./ mp.Ix; end
    if isfield(A, 'CNa'), CN   = tb('CNa'); end
    if isfield(A, 'xcp'), xcp  = tb('xcp'); end
    if isfield(A, 'Cmq'), Cmq  = tb('Cmq'); end
    if isfield(A, 'xcp_canard')
        E.Hd = E.qbar .* P.Sc .* CNc .* (tb('xcp_canard') - P.hinge_frac) * P.mac_c;
    end
end
E.SM = (xcp - mp.xcg) / P.d;                            % static margin [cal]
Cma  = -CN .* E.SM;
wn2  = -Cma .* E.qbar * S * P.d ./ mp.Iy;
E.wn = sqrt(max(wn2, 0));
E.wn(wn2 <= 0) = NaN;                                   % statically unstable
E.zeta = (E.qbar*S.*CN./(mp.m.*Vs) - E.qbar*S*P.d^2.*Cmq./(2*Vs.*mp.Iy)) ./ (2*E.wn);
end

function c = lift_slope(M, AR)          % unswept panel, on panel area
c   = zeros(size(M));
sub = M <= 0.8;  sup = M >= 1.2;  tr = ~sub & ~sup;
c(sub) = slope_sub(M(sub), AR);
c(sup) = slope_sup(M(sup), AR);
c08 = slope_sub(0.8, AR);  c12 = slope_sup(1.2, AR);
c(tr) = c08 + (M(tr) - 0.8)/0.4*(c12 - c08);
end
function c = slope_sub(M, AR)           % Barrowman / Diederich
c = 2*pi*AR ./ (2 + sqrt(4 + (sqrt(1 - M.^2)*AR).^2));
end
function c = slope_sup(M, AR)           % linear theory + rectangular tip loss
B = sqrt(M.^2 - 1);
c = 4./B .* max(1 - 1./(2*B*AR), 0.5);
end

function [rho, a, Tk, p] = isa_atm(h)   % ISA, geopotential altitude, 0-84 km
% Layers (base altitude [km], lapse rate [K/km]):
%   0 -6.5 | 11 0 | 20 +1.0 | 32 +2.8 | 47 0 | 51 -2.8 | 71 -2.0 | top 84 km
% Base temperatures and pressures are built once by integrating upward,
% so T, p and rho are continuous across layer boundaries.
persistent hb Lb Tb pb
g0 = 9.80665;  Rg = 287.05;
if isempty(hb)
    hb = [0 11 20 32 47 51 71]*1e3;
    Lb = [-6.5 0 1.0 2.8 0 -2.8 -2.0]*1e-3;
    Tb = zeros(size(hb));  pb = Tb;  Tb(1) = 288.15;  pb(1) = 101325;
    for j = 1:numel(hb) - 1
        dz = hb(j+1) - hb(j);
        Tb(j+1) = Tb(j) + Lb(j)*dz;
        if Lb(j) == 0
            pb(j+1) = pb(j)*exp(-g0*dz/(Rg*Tb(j)));
        else
            pb(j+1) = pb(j)*(Tb(j+1)/Tb(j))^(-g0/(Rg*Lb(j)));
        end
    end
end
h  = min(max(h, 0), 84000);
Tk = zeros(size(h));  p = Tk;
for j = 1:numel(hb)                     % later layers overwrite lower ones
    i = h >= hb(j);
    if ~any(i(:)), break; end
    Tk(i) = Tb(j) + Lb(j)*(h(i) - hb(j));
    if Lb(j) == 0
        p(i) = pb(j)*exp(-g0*(h(i) - hb(j))/(Rg*Tb(j)));
    else
        p(i) = pb(j)*(Tk(i)/Tb(j)).^(-g0/(Rg*Lb(j)));
    end
end
rho = p ./ (Rg*Tk);
a   = sqrt(1.4*Rg*Tk);
end

function [Kp, Ki, wc] = pi_gains(P, a, b)
% Pole placement on pdot = -a p + b delta:  s^2 + 2 zc wc s + wc^2.
% Bandwidth schedule: at least wc_ratio x the open-loop roll pole a,
% clamped to [wc, wc_max].
wc = min(max(P.wc, P.wc_ratio*a), P.wc_max);
b  = max(b, 1e-6);
Kp = max(2*P.zc*wc - a, 0) ./ b;
Ki = wc.^2 ./ b;
end

function pm = phase_margin(P, a, b, Kp, Ki)
% Loop L = b(Kp s + Ki) / (s (s+a)(tau_s s+1)) * exp(-s Td),
% Td = half-sample hold + sensor delay + servo dead time. Returns PM [deg].
a  = a(:);  b = b(:);
Kp = Kp(:).*ones(size(a));  Ki = Ki(:).*ones(size(a));
w  = logspace(-1, 3.5, 500);
Td = 0.5/P.f_loop + P.sensor_delay + P.servo_delay;
mag = b.*sqrt((Kp.*w).^2 + Ki.^2) ./ (w.*sqrt(w.^2 + a.^2).*sqrt(1 + (w*P.tau_s).^2));
ph  = atan2(Kp.*w, Ki) - pi/2 - atan(w./a) - atan(w*P.tau_s) - w*Td;
pm  = nan(size(a));
for i = 1:numel(a)
    k = find(mag(i,:) < 1, 1);
    if isempty(k) || k == 1, continue; end
    s  = log(mag(i,k-1)) / (log(mag(i,k-1)) - log(mag(i,k)));   % interpolate crossover
    pm(i) = 180 + (ph(i,k-1) + s*(ph(i,k) - ph(i,k-1)))*180/pi;
end
end

function A = amp_factor(p, wn, zeta)
A = wn.^2 ./ sqrt((wn.^2 - p.^2).^2 + (2*zeta.*wn.*p).^2);
end

function cd = drag_coef(M)              % crude whole-rocket drag curve
cd = 0.45*ones(size(M));
i = M > 0.8 & M <= 1.1;  cd(i) = 0.45 + (M(i) - 0.8)/0.3*0.30;
i = M > 1.1;             cd(i) = 0.30 + 0.45*1.1./M(i);
end

function T = fly_trajectory(P)          % 2-D point mass, gravity turn
dt = P.dt;  Nmax = round(P.t_max/dt);
S  = pi*P.d^2/4;  g = 9.80665;  th0 = P.launch_angle;
t = zeros(Nmax,1); h = t; V = t; f = t; offr = false(Nmax,1);
vx = 0; vh = 0; hh = P.h0; Vm = 0; srail = 0; kend = Nmax;   % start at the launch site altitude
for k = 1:Nmax
    tk = (k - 1)*dt;
    fk = max(1 - tk/P.t_burn, 0);
    m  = P.m_dry + fk*P.m_prop;
    Th = P.thrust*(tk < P.t_burn);
    [rho, asnd] = isa_atm(hh);
    if isfield(P, 'aero') && isfield(P.aero, 'CD')
        cdk = interp1(P.aero.M, P.aero.CD, min(max(Vm/asnd, P.aero.M(1)), P.aero.M(end)));
    else
        cdk = drag_coef(Vm/asnd);
    end
    D  = 0.5*rho*Vm^2*S*cdk;
    onrail = srail < P.L_rail;
    if onrail
        Vm = Vm + max((Th - D)/m - g*sin(th0), 0)*dt;
        srail = srail + Vm*dt;
        vx = Vm*cos(th0);  vh = Vm*sin(th0);
    else
        F  = (Th - D)/m;
        vx = vx + F*vx/Vm*dt;
        vh = vh + (F*vh/Vm - g)*dt;
        Vm = hypot(vx, vh);
    end
    hh = hh + vh*dt;
    t(k) = tk; h(k) = hh; V(k) = Vm; f(k) = fk; offr(k) = ~onrail;
    if tk > P.t_burn && vh <= 0, kend = k; break; end
end
T.t = t(1:kend); T.h = h(1:kend); T.V = V(1:kend); T.f = f(1:kend);
T.offrail = offr(1:kend);
[rho, asnd] = isa_atm(T.h);
T.M = T.V ./ asnd;  T.qbar = 0.5*rho.*T.V.^2;
end

function S = roll_sim(P, C, N, keep)
% Vectorised over rows of C: plant pdot = -a p + b delta + dL, first-order
% servo with dead time, rate and travel limits, gyro with delay/bias/noise/range, and
% a PI loop at P.f_loop with zero-order hold and back-calculation
% anti-windup (P.antiwindup: 'measured' uses the measured canard angle,
% 'command' uses sat(v) - v). Optional free play P.backlash (total, rad;
% C.backlash per row overrides it) between servo and canard, with no aero
% preload (worst case): the canard moves only when the servo has taken up
% the play. P.pos_sensor selects whether the measured angle is the servo
% ('servo') or the canard ('canard') side of the play.
n    = size(C.a, 1);  dt = P.dt;
if isfield(C, 'backlash'), Bh = C.backlash(:)/2; else, Bh = P.backlash/2*ones(n,1); end
useC = strcmpi(P.pos_sensor, 'canard');
awC  = strcmpi(P.antiwindup, 'command');
nsub = max(1, round(1/(P.f_loop*dt)));  Ts = nsub*dt;
nd   = round(P.sensor_delay/dt);
na   = round(P.servo_delay/dt);
col  = @(X, k) X(:, min(k, size(X, 2)));
p = zeros(n,1); dl = p; cn = p; xI = p; dcmd = p; pk = p; nsat = p; ntick = 0;
buf  = zeros(n, nd + 1);                              % gyro delay line
cbuf = zeros(n, na + 1);                              % servo command delay line
if keep, S.p = zeros(n, N); S.d = S.p; S.ds = S.p; end
for k = 1:N
    a = col(C.a, k);  b = col(C.b, k);  dL = col(C.dL, k);
    buf(:, mod(k - 1, nd + 1) + 1) = p;
    if mod(k - 1, nsub) == 0                          % controller tick
        on = col(C.on, k);  Kp = col(C.Kp, k);  Ki = col(C.Ki, k);
        pm = buf(:, mod(k - 1 - nd, nd + 1) + 1) + P.gyro_bias + P.gyro_noise*randn(n,1);
        pm = min(max(pm, -P.gyro_range), P.gyro_range);
        e  = P.p_ref - pm;
        v  = Kp.*e + xI;
        v(~on) = xI(~on);            % control off: hold last trim (0 before first enable)
        dcmd = min(max(v, -P.dmax), P.dmax);
        Tt = max(Kp./Ki, P.Tt_min);
        if awC, dm = dcmd; elseif useC, dm = cn; else, dm = dl; end
        xIn = xI + Ts*(Ki.*e + (dm - v)./Tt);
        xI(on) = xIn(on);            % integrator frozen while off
        nsat = nsat + (on & abs(v) > P.dmax);  ntick = ntick + 1;
    end
    cbuf(:, mod(k - 1, na + 1) + 1) = dcmd;
    dapp = cbuf(:, mod(k - 1 - na, na + 1) + 1);      % command after servo dead time
    rate = min(max((dapp - dl)/P.tau_s, -P.dmax_rate), P.dmax_rate);
    dl   = min(max(dl + rate*dt, -P.dmax), P.dmax);    % servo output
    cn   = min(max(cn, dl - Bh), dl + Bh);            % canard, through the free play
    phi  = exp(-a*dt);                                % exact over dt
    p    = phi.*p + (1 - phi)./a.*(b.*cn + dL);
    pk   = max(pk, abs(p));
    if keep, S.p(:, k) = p; S.d(:, k) = cn; S.ds(:, k) = dl; end
end
S.pk = pk;  S.satfrac = nsat/max(ntick, 1);
end

function map_panel(VV, HH, Z, MM, T, ttl, cbl, uselog, Zt, lev)
% Filled map of Z (log10 if uselog), Mach iso-lines (dotted), the
% trajectory (white) and an optional red threshold line Zt = lev.
Zp = Z;  if uselog, Zp = log10(Z); end
Zp(~isfinite(Zp)) = NaN;
contourf(VV, HH/1e3, Zp, 20, 'LineColor', 'none');  hold on
cl = caxis;
[c, hc] = contour(VV, HH/1e3, MM, [1 2 3 4 5], 'k:');
clabel(c, hc, 'FontSize', 7);
if nargin >= 10
    Zt(~isfinite(Zt)) = NaN;
    contour(VV, HH/1e3, Zt, [lev lev], 'r-', 'LineWidth', 1.6);
end
caxis(cl);                                   % keep colour scale on Z only
plot(T.V, T.h/1e3, 'w-', 'LineWidth', 1.8);
cb = colorbar;  ylabel(cb, cbl);
xlabel('Velocity [m/s]');  ylabel('Altitude [km]');  title(ttl);
end

function rate_plot(x, p, wn, cs, nm, lw)
hl = [];
for i = 1:size(p, 1)
    hl = [hl, semilogy(x, max(abs(p(i,:)), 1e-3), 'Color', cs{i}, 'LineWidth', lw)];  hold on
end
hl = [hl, semilogy(x, wn, 'r--', 'LineWidth', lw)];
grid on;  ylabel('|p| [rad/s]');  legend(hl, [nm, {'\omega_n'}], 'Location', 'best');
end

function burn_line(P)
hold on;  yl = ylim;  plot(P.t_burn*[1 1], yl, 'k:', 'HandleVisibility', 'off');  ylim(yl);
end