function R = roll_rl_flightdata(csvfile, matfile, spec)
%ROLL_RL_FLIGHTDATA  PI roll-rate gains by root locus from the team's flight data.
%
%   R = roll_rl_flightdata('flight_pressure_fed_regulator.csv', 'canard_roll_force_1.mat')
%   R = roll_rl_flightdata(csv, mat, struct('wz', 6, 'zeta', 0.7))
%
%   Inputs
%     csvfile  trajectory output of the vehicle sizer (one row per time step)
%     matfile  canard force table: force(Mach, dynp, defl)
%
%   STEP 1  Roll plant  pdot = -a p + b delta + dL  at every time step
%     canard moment  Lc = q * G * CN(M, delta)       CN from the .mat table,
%                    G  = max_roll_moment / (q CN(M, 10 deg)) from the CSV
%                    (one constant; the CSV moments are reproduced exactly)
%     b  = q G dCN/ddelta(M, delta_ref) / Ixx        control effectiveness
%     dL = -control_roll_moment / Ixx                fin-cant disturbance (the
%                    sizer's canards cancel it exactly, so it is minus their moment)
%     a  = (fin moment per rad of cant) r_eff / (V Ixx)   roll damping, strip
%                    theory on the same fins: ESTIMATE, not in the data files
%
%   STEP 2  Open loop with the loop gain K = Kp*b pulled out:
%       L(s) = K (s + wz) / ( s (s + a) (tau_s s + 1) ) * Pade2(Td)
%
%   STEP 3  Root locus: fix the PI zero wz, sweep K, take the largest K for
%     which every oscillatory closed-loop pole still has zeta >= spec.zeta.
%
%   STEP 4  Gains:  Kp = K / b,  Ki = wz * Kp   (b changes ~300x in flight;
%     K and wz hardly change because a << the loop bandwidth).
%
%   STEP 5  Checks: phase/gain margin, and a time simulation of the whole
%     flight (needs canard_roll_envelope.m on the path): roll rate with the
%     canards locked vs controlled, against the pitch natural frequency.
%
%   Files written: roll_gain_schedule_flight.csv (per time step),
%   roll_gain_table_mach.csv (Mach -> Kp q/Ixx, Ki q/Ixx for the flight computer)
%   and roll_gain_table.mat (struct roll_gain_table: the Mach table, the design
%   constants and the per-time-step schedule; R.table is the same struct).
%
%   spec fields (defaults)
%     wz 6 [rad/s] PI zero, or 'auto' (largest K*wz that still reaches zeta)
%     zeta 0.7             damping target on the locus
%     tau_s 0.020 [s]      servo time constant            (ASSUMED)
%     Td  [] -> 0.5/f_loop + 0.002 + 0.003 s              (ASSUMED)
%     f_loop 960 [Hz]; qbar_min 2000 [Pa] control enabled above this
%     delta_ref 3 [deg]    deflection at which dCN/ddelta is taken
%     fin_cant 0.00174533 [rad] cant used by the sizer (YAML roll.fin_cant)
%     fin [span root tip] = [0.2981 0.508 0.2032] m, fin_r0 0.1475 m (root radius)
%     dmax 10 [deg] canard travel; plot true; quiet false
%   Verification (see roll_rl_verify.m): the gains stay as designed while the
%   "true" vehicle in the checks and in the flight simulation is perturbed:
%     b_scale 1      true canard effectiveness / model (e.g. 2 if eta is 1.0, not 0.5)
%     dL_scale 1     true roll disturbance / sizer value (2.5 = 0.25 deg fin cant)
%     tau_true []    true servo time constant [s] (default tau_s)
%     delay_extra 0  extra loop delay [s]
%     backlash_deg 0 canard free play, total [deg]
%     gyro_bias_dps [] residual gyro bias [deg/s] (default 0.5)
%   Runs in MATLAB (no toolboxes) and GNU Octave.

if nargin < 3, spec = struct(); end
def = struct('wz', 6, 'zeta', 0.7, 'tau_s', 0.020, 'Td', [], 'f_loop', 960, 'qbar_min', 2000, ...
             'delta_ref', 3, 'fin_cant', 0.00174533, 'fin', [0.2980845 0.508 0.2032], ...
             'fin_r0', 0.1475, 'dmax', 10, 'plot', true, 'quiet', false, ...
             'b_scale', 1, 'dL_scale', 1, 'tau_true', [], 'delay_extra', 0, ...
             'backlash_deg', 0, 'gyro_bias_dps', []);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(spec, fn{i}), spec.(fn{i}) = def.(fn{i}); end, end
if isempty(spec.Td), spec.Td = 0.5/spec.f_loop + 0.002 + 0.003; end
if isempty(spec.tau_true), spec.tau_true = spec.tau_s; end
tau = spec.tau_s;  Td = spec.Td;

%% ---- STEP 1: data -> roll plant ---------------------------------------
want = {'state.kinematics.t', 'state.kinematics.h', 'state.atmosphere.Ma', 'state.atmosphere.q', ...
        'state.wind.airspeed', 'state.mass_properties.Ixx', ...
        'state.plant.aero.canard.deflection_rad', 'state.plant.aero.canard.control_roll_moment_nm', ...
        'state.plant.aero.canard.max_roll_moment_nm', 'state.plant.aero.pitch_hz'};
D = read_csv_columns(csvfile, want);
t = D(:,1);  h = D(:,2);  Ma = D(:,3);  q = D(:,4);  V = D(:,5);  Ixx = D(:,6);
dtrim = D(:,7);  Lc = D(:,8);  Lmax = D(:,9);  wn = 2*pi*D(:,10);
[t, iu] = unique(t);                                    % drop repeated time stamps
h = h(iu); Ma = Ma(iu); q = q(iu); V = V(iu); Ixx = Ixx(iu);
dtrim = dtrim(iu); Lc = Lc(iu); Lmax = Lmax(iu); wn = wn(iu);

S  = load(matfile);  tb = S.table;
Mg = squeeze(tb.mach(:,1,1));  Dg = squeeze(tb.defl(1,1,:))*pi/180;
CN = squeeze(tb.force(:,end,:)) / (tb.dynp(1,end,1)*1e3*tb.canard_area_m2);   % CN(Mach, defl)
cn = @(M, d) interp2(Dg(:).', Mg(:), CN, d, min(max(M, Mg(1)), Mg(end)));

ok   = isfinite(Lmax) & isfinite(Lc) & q > 50;
Gk   = Lmax ./ (q .* cn(Ma, spec.dmax*pi/180 + 0*Ma));
G    = median(Gk(ok));                                  % geometry factor [m^3]
chk  = max(abs(q(ok)*G.*cn(Ma(ok), dtrim(ok)) ./ Lc(ok) - 1));   % table vs CSV moment

dd   = 0.1*pi/180;  dref = spec.delta_ref*pi/180;
CNd  = (cn(Ma, dref + dd + 0*Ma) - cn(Ma, dref - dd + 0*Ma)) / (2*dd);   % dCN/ddelta [1/rad]
b    = q .* G .* CNd ./ Ixx;                            % [1/s^2]
dL   = -Lc ./ Ixx;                                      % [rad/s^2]
% roll damping from the fin-cant moment (same fin lift slope), strip theory
y    = linspace(0, spec.fin(1), 2001);  c = spec.fin(2) + (spec.fin(3) - spec.fin(2))*y/spec.fin(1);
Y    = spec.fin_r0 + y;
reff = trapz(y, c.*Y.^2) / trapz(y, c.*Y);              % <y^2>/<y> [m]
a    = (-Lc/spec.fin_cant) * reff ./ (V .* Ixx);        % [1/s]
pss  = dL ./ a;                                         % canards locked: steady spin [rad/s]

on   = ok & q >= spec.qbar_min;                         % controller enabled
k    = find(on);  k = k(:);
gap  = find(diff(k) > 1, 1);                            % ascent only: first unbroken block
if ~isempty(gap), k = k(1:gap); end
on   = false(size(on));  on(k) = true;
[~, iq] = max(q);                                       % max-q row

%% ---- STEP 2-3: root locus -> loop gain K and PI zero wz ----------------
% the design depends on the flight condition only through a: do it on a
% small grid of a and interpolate
ag = linspace(min(a(k)), max(a(k)), 9);
Kg = zeros(size(ag));  wzg = Kg;
for j = 1:numel(ag), [Kg(j), wzg(j)] = rl_design(ag(j), tau, Td, spec.zeta, spec.wz); end
K  = interp1(ag, Kg, a(k));  wz = interp1(ag, wzg, a(k));

%% ---- STEP 4: gains -----------------------------------------------------
Kp = K ./ b(k);  Ki = wz .* Kp;

%% ---- STEP 5: checks ----------------------------------------------------
[pm, gm, wc] = margins(a(k), K, wz, tau, Td);
% margins of the perturbed ("true") loop with the designed gains
[pmT, gmT] = margins(a(k), K*spec.b_scale, wz, spec.tau_true, Td + spec.delay_extra);
jq = find(k == iq, 1);  if isempty(jq), [~, jq] = max(q(k)); end
kv = logspace(-1, log10(400), 800);                     % locus at max q
rloc = zeros(5, numel(kv));
for j = 1:numel(kv), rloc(:, j) = cl_poles(a(k(jq)), kv(j), wz(jq), tau, Td); end
pdes = cl_poles(a(k(jq)), K(jq), wz(jq), tau, Td);
kcr  = k_unstable(a(k(jq)), wz(jq), tau, Td);

sim = [];
if exist('canard_roll_envelope', 'file') == 2
    sim = fly_roll(t, q, a, b, dL, wn, on, k, Kp, Ki, spec);
end

%% ---- files -------------------------------------------------------------
fid = fopen('roll_gain_schedule_flight.csv', 'w');
fprintf(fid, 't_s,h_m,mach,qbar_Pa,Ixx_kgm2,a_1ps,b_1ps2,dL_radps2,Kp_s,Ki\n');
fprintf(fid, '%.3f,%.1f,%.4f,%.1f,%.4f,%.5f,%.5f,%.5f,%.6g,%.6g\n', ...
        [t(k) h(k) Ma(k) q(k) Ixx(k) a(k) b(k) dL(k) Kp Ki].');
fclose(fid);
% flight-computer table: Kp = f(M) * Ixx / qbar  with  f = K / (G dCN/ddelta)
Mt  = (0.1:0.1:floor(10*max(Ma(k)))/10).';
CNt = (cn(Mt, dref + dd + 0*Mt) - cn(Mt, dref - dd + 0*Mt)) / (2*dd);
Kt  = median(K);  wzt = median(wz);
fid = fopen('roll_gain_table_mach.csv', 'w');
fprintf(fid, 'mach,dCN_ddelta_per_rad,Kp_times_qbar_over_Ixx,Ki_times_qbar_over_Ixx\n');
fprintf(fid, '%.2f,%.5f,%.6g,%.6g\n', [Mt, CNt, Kt./(G*CNt), wzt*Kt./(G*CNt)].');
fclose(fid);
% the same table as a struct in roll_gain_table.mat (variable roll_gain_table)
roll_gain_table = struct( ...
    'mach',             Mt, ...                         % [-]      table breakpoints
    'dCN_ddelta',       CNt, ...                        % [1/rad]  canard force slope at delta_ref
    'Kp_qbar_over_Ixx', Kt./(G*CNt), ...                % Kp = this * Ixx / qbar   [s]
    'Ki_qbar_over_Ixx', wzt*Kt./(G*CNt), ...            % Ki = this * Ixx / qbar   [-]
    'K', Kt, 'wz', wzt, 'zeta', spec.zeta, ...          % loop gain [1/s], PI zero [rad/s], damping target
    'G', G, 'delta_ref_deg', spec.delta_ref, ...        % canard geometry factor [m^3]
    'tau_s', tau, 'Td', Td, 'f_loop', spec.f_loop, ...  % assumed servo time constant, loop delay [s]
    'qbar_min', spec.qbar_min, 'dmax_deg', spec.dmax, ...
    'usage', 'Kp = interp1(mach, Kp_qbar_over_Ixx, M) * Ixx / qbar;  Ki = wz * Kp  (SI units, angles in rad)', ...
    'schedule', struct('t', t(k), 'h', h(k), 'mach', Ma(k), 'qbar', q(k), 'Ixx', Ixx(k), ...
                       'a', a(k), 'b', b(k), 'dL', dL(k), 'Kp', Kp, 'Ki', Ki), ...   % along this flight
    'source_csv', csvfile, 'source_mat', matfile, 'created', datestr(now, 'yyyy-mm-dd HH:MM'));
if exist('OCTAVE_VERSION', 'builtin')
    save('-v7', 'roll_gain_table.mat', 'roll_gain_table');
else
    save('roll_gain_table.mat', 'roll_gain_table', '-v7');
end

%% ---- report ------------------------------------------------------------
if ~spec.quiet
fprintf('\nSTEP 1  plant from the data (G = %.5f m^3, table reproduces CSV canard moments within %.3f%%)\n', G, 100*chk);
fprintf('        roll damping a is an ESTIMATE (strip theory, r_eff = %.3f m)\n', reff);
fprintf('%6s %7s %6s %8s %7s %8s %9s %9s %9s %8s %8s\n', 't[s]', 'h[km]', 'Mach', 'q[kPa]', 'Ixx', ...
        'a[1/s]', 'b[1/s2]', 'dL[r/s2]', 'p_ss[r/s]', 'wn[r/s]', 'trim[deg]');
sel = k(round(linspace(1, numel(k), 9)));
pts = unique([sel(:); iq]);
for i = pts.'
    fprintf('%6.1f %7.2f %6.2f %8.1f %7.2f %8.3f %9.2f %9.3f %9.2f %8.2f %8.2f\n', t(i), h(i)/1e3, Ma(i), ...
            q(i)/1e3, Ixx(i), a(i), b(i), dL(i), pss(i), wn(i), dtrim(i)*180/pi);
end
fprintf('\nSTEP 2-3  root locus (tau_s = %.0f ms, Td = %.2f ms, zeta target %.2f)\n', 1e3*tau, 1e3*Td, spec.zeta);
fprintf('  at max q (t = %.1f s): open-loop poles 0, %.3f, %.0f; PI zero %.2f rad/s\n', t(iq), -a(iq), -1/tau, -wz(jq));
fprintf('  loop gain K = Kp*b = %.2f 1/s  ->  closed-loop poles:', K(jq));
ps = pdes(abs(pdes) < 0.3*max(abs(pdes)));
fprintf(' %.1f%+.1fj', [real(ps(:)) imag(ps(:))].');  fprintf('\n');
fprintf('  locus enters the right half plane at K = %.0f (gain margin %.1f dB)\n', kcr, 20*log10(kcr/K(jq)));
fprintf('  over the flight: K %.2f-%.2f 1/s, wz %.2f-%.2f rad/s\n', min(K), max(K), min(wz), max(wz));
fprintf('\nSTEP 4  gains Kp = K/b, Ki = wz*Kp\n');
fprintf('%6s %6s %8s %9s %9s %9s %7s %7s %7s\n', 't[s]', 'Mach', 'q[kPa]', 'b[1/s2]', 'Kp[s]', 'Ki[-]', 'PM', 'GM[dB]', 'fc[Hz]');
for i = pts.'
    j = find(k == i, 1);  if isempty(j), continue; end
    fprintf('%6.1f %6.2f %8.1f %9.2f %9.4f %9.4f %7.1f %7.1f %7.2f\n', t(i), Ma(i), q(i)/1e3, b(i), ...
            Kp(j), Ki(j), pm(j), gm(j), wc(j)/2/pi);
end
fprintf('  controlled from t = %.1f s to %.1f s; Kp %.3g-%.3g s (factor %.0f)\n', t(k(1)), t(k(end)), ...
        min(Kp), max(Kp), max(Kp)/min(Kp));
fprintf('  margins over the flight: PM %.1f-%.1f deg, GM %.1f-%.1f dB, crossover %.2f-%.2f Hz\n', ...
        min(pm), max(pm), min(gm), max(gm), min(wc)/2/pi, max(wc)/2/pi);
if ~isempty(sim)
    fprintf('\nSTEP 5  flight simulation, fin cant %.2f deg (q >= %.0f kPa)\n', spec.fin_cant*180/pi, spec.qbar_min/1e3);
    fprintf('  canards locked : peak |p| %.2f rad/s, max |p|/wn %.2f\n', sim.pk(1), sim.ratio(1));
    fprintf('  root-locus PI  : peak |p| %.3f rad/s, max |p|/wn %.3f, max |delta| %.2f deg of %.0f\n', ...
            sim.pk(2), sim.ratio(2), sim.dmax_used, spec.dmax);
end
fprintf('  wrote roll_gain_schedule_flight.csv, roll_gain_table_mach.csv and roll_gain_table.mat\n');
end

R = struct('t', t, 'h', h, 'M', Ma, 'qbar', q, 'Ixx', Ixx, 'a', a, 'b', b, 'dL', dL, 'pss', pss, 'wn', wn, ...
           'dtrim', dtrim, 'G', G, 'reff', reff, 'on', on, 'idx', k, 'K', K, 'wz', wz, 'Kp', Kp, 'Ki', Ki, ...
           'pm', pm, 'gm', gm, 'wc', wc, 'pm_true', pmT, 'gm_true', gmT, 'iq', iq, 'locus_k', kv, 'locus', rloc, 'design_poles', pdes, ...
           'k_unstable', kcr, 'sim', sim, 'spec', spec, 'table', roll_gain_table);

if ~spec.plot, return; end
%% ---- figures -----------------------------------------------------------
figure('Name', 'Root locus at max q', 'Color', 'w');
plot(real(rloc(:)), imag(rloc(:)), 'b.', 'MarkerSize', 4);  hold on
plot(real(pdes), imag(pdes), 'rs', 'MarkerSize', 9, 'LineWidth', 1.5);
plot([0 -a(iq) -1/tau], [0 0 0], 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(-wz(jq), 0, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
th = acos(spec.zeta);
plot([0 -200*cos(th)], [0 200*sin(th)], 'k:', [0 -200*cos(th)], [0 -200*sin(th)], 'k:', [0 0], [-200 200], 'k-');
axis equal;  axis([-70 15 -50 50]);  grid on
xlabel('Re(s) [rad/s]');  ylabel('Im(s) [rad/s]');
title(sprintf('Root locus at max q: K = %.1f, w_z = %.1f rad/s', K(jq), wz(jq)));
legend('locus', 'chosen K', 'open-loop poles', 'PI zero', sprintf('\\zeta = %.1f', spec.zeta), 'Location', 'best');

figure('Name', 'Gain schedule', 'Color', 'w');
subplot(2,1,1);  semilogy(t(k), Kp, 'b-', t(k), Ki, 'r-', 'LineWidth', 1.3);  grid on
xlabel('t [s]');  legend('K_p [s]', 'K_i [-]');  title('Scheduled PI gains');
subplot(2,1,2);  semilogy(t(k), b(k), 'k-', 'LineWidth', 1.3);  grid on
xlabel('t [s]');  ylabel('b [1/s^2]');  title('Canard effectiveness');

if ~isempty(sim)
    figure('Name', 'Flight simulation', 'Color', 'w');
    subplot(2,1,1);
    semilogy(sim.t, max(abs(sim.p(1,:)), 1e-4), 'k-', sim.t, max(abs(sim.p(2,:)), 1e-4), 'b-', ...
             sim.t, sim.wn, 'r--', 'LineWidth', 1.3);  grid on
    xlabel('t [s]');  ylabel('[rad/s]');  legend('|p| canards locked', '|p| root-locus PI', '\omega_n pitch');
    title('Roll rate vs pitch natural frequency');
    subplot(2,1,2);
    plot(sim.t, sim.d(2,:)*180/pi, 'b-', sim.t([1 end]), spec.dmax*[1 1], 'r:', ...
         sim.t([1 end]), -spec.dmax*[1 1], 'r:', 'LineWidth', 1.3);  grid on
    xlabel('t [s]');  ylabel('\delta [deg]');  title('Canard deflection');
end
end

%% =======================================================================
function D = read_csv_columns(file, want)
% Numeric columns by header name; non-numeric cells (True/False, empty) -> NaN.
fid = fopen(file, 'r');
if fid < 0, error('cannot open %s', file); end
hdr = regexp(strtrim(fgetl(fid)), ',', 'split');
idx = zeros(1, numel(want));
for j = 1:numel(want)
    i = find(strcmp(hdr, want{j}), 1);
    if isempty(i), fclose(fid); error('column "%s" not found in %s', want{j}, file); end
    idx(j) = i;
end
D = zeros(0, numel(want));
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    c = regexp(ln, ',', 'split');
    if numel(c) < max(idx), continue; end
    D(end+1, :) = str2double(c(idx)); %#ok<AGROW>
end
fclose(fid);
end

function r = cl_poles(a, K, wz, tau, Td)
% roots of  s (s+a)(tau s+1) Dden(s) + K (s+wz) Dnum(s),  2nd-order Pade delay
pd = [Td^2/12, Td/2, 1];  pn = [Td^2/12, -Td/2, 1];
den = conv(conv([1 0], [1 a]), conv([tau 1], pd));
num = K*conv([1 wz], pn);
r = roots(den + [zeros(1, numel(den) - numel(num)), num]);
end

function z = min_zeta(r)
% smallest damping ratio among the physical closed-loop poles (the Pade pair
% is an artefact of the delay approximation); -1 if any pole is unstable
if any(real(r) > 0), z = -1; return; end
r = r(abs(r) < 0.3*max(abs(r)));
c = r(abs(imag(r)) > 1e-9);
if isempty(c), z = 1; else, z = min(-real(c)./abs(c)); end
end

function [K, wz] = rl_design(a, tau, Td, zt, wzopt)
% Largest loop gain K with min damping >= zt, for a fixed PI zero; with
% wzopt = 'auto' the zero is the one that maximises K*wz (integral gain).
kv = logspace(0, log10(300), 500);
if ischar(wzopt), cand = 1:0.25:15; else, cand = wzopt; end
best = -1;  K = NaN;  wz = NaN;
for w = cand
    z = zeros(size(kv));
    for j = 1:numel(kv), z(j) = min_zeta(cl_poles(a, kv(j), w, tau, Td)); end
    j = find(z >= zt, 1, 'last');
    if isempty(j), continue; end
    if j < numel(kv)                                   % refine the crossing
        lo = kv(j);  hi = kv(j+1);
        for it = 1:30
            m = sqrt(lo*hi);
            if min_zeta(cl_poles(a, m, w, tau, Td)) >= zt, lo = m; else, hi = m; end
        end
        Kw = lo;
    else
        Kw = kv(j);
    end
    if Kw*w > best, best = Kw*w;  K = Kw;  wz = w; end
end
if isnan(K), error('zeta = %.2f cannot be reached with this PI zero; lower wz or zeta', zt); end
end

function kc = k_unstable(a, wz, tau, Td)
lo = 1;  hi = 1e4;
for it = 1:60
    m = sqrt(lo*hi);
    if any(real(cl_poles(a, m, wz, tau, Td)) > 0), hi = m; else, lo = m; end
end
kc = sqrt(lo*hi);
end

function [pm, gm, wgc] = margins(a, K, wz, tau, Td)
w = logspace(-1, 4, 4000);
a = a(:);  K = K(:);  wz = wz(:);
mag = K.*sqrt(w.^2 + wz.^2) ./ (w.*sqrt(w.^2 + a.^2).*sqrt(1 + (w*tau).^2));
ph  = atan2(w + 0*wz, wz + 0*w) - pi/2 - atan(w./a) - atan(w*tau) - w*Td;
pm = nan(size(a));  gm = pm;  wgc = pm;
for i = 1:numel(a)
    j = find(mag(i,:) < 1, 1);
    if ~isempty(j) && j > 1
        s = log(mag(i,j-1)) / (log(mag(i,j-1)) - log(mag(i,j)));
        pm(i)  = 180 + (ph(i,j-1) + s*(ph(i,j) - ph(i,j-1)))*180/pi;
        wgc(i) = w(j-1) * (w(j)/w(j-1))^s;
    end
    j = find(ph(i,:) < -pi, 1);
    if ~isempty(j) && j > 1
        s = (ph(i,j-1) + pi) / (ph(i,j-1) - ph(i,j));
        gm(i) = -20*log10(mag(i,j-1) * (mag(i,j)/mag(i,j-1))^s);
    end
end
end

function sim = fly_roll(t, q, a, b, dL, wn, on, k, Kp, Ki, spec)
% Whole-flight roll simulation with the nonlinear loop of
% canard_roll_envelope.m (servo lag, dead time, rate and travel limits, gyro
% delay/bias/noise, 960 Hz PI with anti-windup): canards locked vs PI.
F  = canard_roll_envelope(struct('dmax', spec.dmax*pi/180, 'f_loop', spec.f_loop, ...
                                 'tau_s', spec.tau_true, 'antiwindup', 'command', ...
                                 'backlash', spec.backlash_deg*pi/180), 'tools');
P  = F.P;
P.sensor_delay = P.sensor_delay + spec.delay_extra;
if ~isempty(spec.gyro_bias_dps), P.gyro_bias = spec.gyro_bias_dps*pi/180; end
i0 = find(isfinite(a) & isfinite(b) & q > 50, 1);       % first usable row (off the rail)
i1 = find(isfinite(a) & isfinite(b) & q > 50, 1, 'last');
ts = (t(i0):P.dt:t(i1));  N = numel(ts);
g  = isfinite(a) & isfinite(b);
ip = @(x) interp1(t(g), x(g), ts, 'linear');
A  = max(ip(a), 1e-6);  B = spec.b_scale*ip(b);  DL = spec.dL_scale*ip(dL);  WN = ip(wn);
ON = interp1(t, double(on), ts, 'nearest') > 0.5;
KP = interp1(t(k), Kp, ts, 'linear', 'extrap');  KI = interp1(t(k), Ki, ts, 'linear', 'extrap');
KP(~ON) = 0;  KI(~ON) = 0;
C  = struct('a', [A; A], 'b', [B; B], 'dL', [DL; DL], 'Kp', [0*KP; KP], 'Ki', [0*KI; KI], ...
            'on', [false(1, N); ON]);
rng(1);  S = F.roll_sim(P, C, N, true);
aero = ip(q) >= spec.qbar_min;                          % where resonance matters
sim.t = ts;  sim.p = S.p;  sim.d = S.d;  sim.wn = WN;  sim.on = ON;
sim.pk    = max(abs(S.p(:, aero)), [], 2);
sim.ratio = max(abs(S.p(:, aero)) ./ WN(aero), [], 2);
sim.dmax_used = max(abs(S.d(2, :)))*180/pi;
sim.sat_time  = sum(abs(S.d(2, :)) >= spec.dmax*pi/180 - 1e-5)*P.dt;   % [s] at the travel limit
end