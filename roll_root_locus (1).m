function RL = roll_root_locus(Pover, spec)
%ROLL_ROOT_LOCUS  Root locus of the PI roll-rate loop.
%
%   RL = roll_root_locus(yjsp_params())
%   RL = roll_root_locus(yjsp_params(), struct('zeta', 0.7, 'plot', true))
%
%   Loop at a frozen flight point (canard_linearization.pdf, Sec. 17):
%       L(s) = k * Kp * b (s + wz) / ( s (s + a) (tau_s s + 1) ) * D(s)
%   D(s) = 2nd-order Pade approximation of the loop delay exp(-s Td).
%   The PI zero wz = Ki/Kp is held fixed and the gain multiplier k is swept
%   (k = 1 is the design from roll_gain_design). Closed-loop poles are the
%   roots of  s (s+a)(tau_s s+1) Dden(s) + k Kp b (s + wz) Dnum(s) = 0.
%
%   Outputs / figures
%   1. Root locus at max q: branches for k = 0.01 ... 300, design poles
%      (k = 1), zeta lines, and the k where a branch crosses into the right
%      half plane (= gain margin; compared with the Bode value).
%   2. Root-locus gain selection: for the same wz, the largest k that keeps
%      the dominant pair at zeta >= spec.zeta, and the best zeta this PI
%      zero can give, at several flight points.
%   3. Closed-loop dominant poles along the trajectory: scheduled gains
%      (stay in place) vs fixed gains designed at (P.ref_h, P.ref_V), which
%      slide along the locus as b ~ qbar changes by ~80x.

if nargin < 1, Pover = []; end
if nargin < 2, spec = struct(); end
if ~isfield(spec, 'zeta'), spec.zeta = 0.7;  end
if ~isfield(spec, 'plot'), spec.plot = true; end
spec.quiet = true;

F  = canard_roll_envelope(Pover, 'tools');  P = F.P;
G  = roll_gain_design(Pover, spec);
tau = P.tau_s;  Td = G.Td;
kv = logspace(-2, log10(300), 700);

%% 1. locus at max q
[~, iq] = max(G.qbar);
[r1, kc] = locus(G.a(iq), G.b(iq), G.Kp(iq), G.wz(iq), tau, Td, kv);
p1  = locus(G.a(iq), G.b(iq), G.Kp(iq), G.wz(iq), tau, Td, 1);
gmL = 20*log10(kc);

%% 2. root-locus gain selection at several points
[~, i08] = min(abs(G.M - 0.8));
pts  = [1, i08, iq, find(G.h >= 20e3, 1), numel(G.t)];
nms  = {'control on', 'Mach 0.8', 'max q', '20 km', 'control off'};
ksel = nan(size(pts));  zdes = ksel;  zmx = ksel;  kmx = ksel;
for j = 1:numel(pts)
    i = pts(j);
    [ksel(j), zmx(j), kmx(j)] = k_for_zeta(G.a(i), G.b(i), G.Kp(i), G.wz(i), tau, Td, ...
                                           spec.zeta, 0.3*2*pi*G.spec.fc);
    [~, zdes(j)] = dominant(locus(G.a(i), G.b(i), G.Kp(i), G.wz(i), tau, Td, 1));
end

%% 3. dominant poles along the trajectory: scheduled vs fixed gains
c   = 1:max(1, round(numel(G.t)/400)):numel(G.t);
mpR = F.mass_props(P, P.map_f);
ER  = F.aero_roll(P, P.ref_h, P.ref_V, mpR, P.eta_model);
[KpF, wzF] = shape(ER.a, ER.b, tau, Td, 2*pi*G.spec.fc, G.spec.pm);
dS = zeros(numel(c), 1);  dF = dS;
for j = 1:numel(c)
    i = c(j);
    dS(j) = dominant(locus(G.a(i), G.b(i), G.Kp(i), G.wz(i), tau, Td, 1));
    dF(j) = dominant(locus(G.a(i), G.b(i), KpF, wzF, tau, Td, 1));
end
zS = -real(dS)./abs(dS);  zF = -real(dF)./abs(dF);
unstF = sum(real(dF) > 0)/numel(c);

%% report
fprintf('\nRoot locus of the PI roll loop (tau_s %.0f ms, Td %.2f ms, 2nd-order Pade)\n', 1e3*tau, 1e3*Td);
fprintf('Max q (t = %.1f s, M %.2f, q %.0f kPa): open-loop poles 0, %.2f, %.1f; PI zero %.2f rad/s\n', ...
        G.t(iq), G.M(iq), G.qbar(iq)/1e3, -G.a(iq), -1/tau, -G.wz(iq));
fprintf('  design poles (k = 1):');  fprintf(' %.1f%+.1fi', [real(p1(:)) imag(p1(:))].');  fprintf('\n');
fprintf('  locus crosses into RHP at k = %.1f -> gain margin %.1f dB (Bode: %.1f dB)\n', ...
        kc, gmL, G.gm(find(G.check_idx >= iq, 1)));
fprintf('\nRoot-locus gain selection, same PI zero, target zeta = %.2f:\n', spec.zeta);
fprintf('%-12s %7s %14s %10s %12s %16s\n', 'point', 'Mach', 'zeta at k = 1', 'k(zeta)', 'Kp_RL [s]', 'best zeta (at k)');
for j = 1:numel(pts)
    i = pts(j);
    fprintf('%-12s %7.2f %14.2f %10.2f %12.4f %10.2f (%.2f)\n', nms{j}, G.M(i), zdes(j), ksel(j), ...
            ksel(j)*G.Kp(i), zmx(j), kmx(j));
end
fprintf('  (k(zeta) = NaN: this PI zero cannot reach zeta = %.2f; move the zero, i.e. change Ki/Kp)\n', spec.zeta);
fprintf('\nAlong the trajectory, dominant closed-loop pair:\n');
fprintf('  scheduled gains: |s| %.1f-%.1f rad/s, zeta %.2f-%.2f\n', min(abs(dS)), max(abs(dS)), min(zS), max(zS));
fprintf('  fixed gains (designed at h = %.0f m, V = %.0f m/s): |s| %.1f-%.1f rad/s, zeta %.2f-%.2f, unstable %.0f%% of the time\n', ...
        P.ref_h, P.ref_V, min(abs(dF)), max(abs(dF)), min(zF), max(zF), 100*unstF);

RL = struct('k', kv, 'roots_maxq', r1, 'design_poles_maxq', p1, 'k_cross', kc, 'gm_dB', gmL, ...
            'ol_poles', [0, -G.a(iq), -1/tau], 'ol_zero', -G.wz(iq), 'Td', Td, 'tau', tau, ...
            'pts', pts, 'names', {nms}, 'k_zeta', ksel, 'zeta_design', zdes, 'zeta_best', zmx, 'k_best', kmx, ...
            't', G.t(c), 'qbar', G.qbar(c), 'dom_sched', dS, 'dom_fixed', dF, 'spec', spec);

if ~spec.plot, return; end
%% figures
figure('Name', 'Root locus at max q', 'Color', 'w', 'Position', [60 60 1150 480]);
for sp = 1:2
    subplot(1, 2, sp);
    plot(real(r1(:)), imag(r1(:)), '.', 'Color', [0.2 0.4 0.8], 'MarkerSize', 4);  hold on
    plot(real(p1), imag(p1), 'rs', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot([0 -G.a(iq) -1/tau], [0 0 0], 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
    plot(-G.wz(iq), 0, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
    R = 400;  if sp == 2, R = 40; end
    th = acos(spec.zeta);
    plot([0 -R*cos(th)], [0 R*sin(th)], 'k:', [0 -R*cos(th)], [0 -R*sin(th)], 'k:');
    plot([0 0], [-R R], 'k-', 'LineWidth', 0.5);
    axis equal;  grid on;
    if sp == 1, axis([-120 30 -90 90]);  title('Root locus at max q (k = 0.01 ... 300)');
    else,       axis([-30 5 -20 20]);    title('Zoom: dominant pair'); end
    xlabel('Re(s) [rad/s]');  ylabel('Im(s) [rad/s]');
end
legend({'locus', 'design (k = 1)', 'open-loop poles', 'PI zero', sprintf('\\zeta = %.1f', spec.zeta)}, ...
       'Location', 'best');

figure('Name', 'Dominant poles along the trajectory', 'Color', 'w', 'Position', [80 80 700 520]);
plot(real(dF), imag(dF), '.', 'Color', [0.85 0.4 0], 'MarkerSize', 8);  hold on
plot(real(dS), imag(dS), '.', 'Color', [0.1 0.3 0.8], 'MarkerSize', 8);
th = acos(spec.zeta);
plot([0 -60*cos(th)], [0 60*sin(th)], 'k:', [0 0], [-5 60], 'k-');
grid on;  xlabel('Re(s) [rad/s]');  ylabel('Im(s) [rad/s]');
legend('fixed gains', 'scheduled gains', sprintf('\\zeta = %.1f', spec.zeta), 'Location', 'best');
title('Dominant closed-loop poles over the controlled flight');
end

%% -----------------------------------------------------------------------
function [r, kc] = locus(a, b, Kp, wz, tau, Td, kv)
% Closed-loop roots for each gain multiplier in kv (columns), and the
% smallest multiplier at which any root enters the right half plane.
pd = [Td^2/12, Td/2, 1];  pn = [Td^2/12, -Td/2, 1];     % 2nd-order Pade
den = conv(conv([1 0], [1 a]), conv([tau 1], pd));
num = b*Kp*conv([1 wz], pn);
num = [zeros(1, numel(den) - numel(num)), num];
r = zeros(numel(den) - 1, numel(kv));
for j = 1:numel(kv), r(:, j) = roots(den + kv(j)*num); end
kc = NaN;
j = find(max(real(r), [], 1) > 0, 1);
if ~isempty(j)
    if j == 1, kc = kv(1); return; end
    lo = kv(j-1);  hi = kv(j);                          % bisect the crossing
    for it = 1:40
        m = sqrt(lo*hi);
        if max(real(roots(den + m*num))) > 0, hi = m; else, lo = m; end
    end
    kc = sqrt(lo*hi);
end
end

function [d, z] = dominant(r)
% Slowest oscillatory pole (upper half plane) of the physical loop, or the
% slowest real pole if there is none. The 2nd-order Pade pair (|s| ~ 2 sqrt(3)/Td,
% zeta = 0.87) is an artefact of the delay approximation and is excluded.
r = r(abs(r) < 0.3*max(abs(r)));
c = r(imag(r) > 1e-9);  if isempty(c), c = r; end
[~, j] = min(abs(c));  d = c(j);  z = -real(d)/abs(d);
end

function [k, zmax, kmax] = k_for_zeta(a, b, Kp, wz, tau, Td, zt, wmin)
% Root-locus gain selection with the PI zero fixed. Along the locus the
% dominant pair leaves the real axis (zeta drops), swings round the PI zero
% (zeta rises to a maximum) and finally bends toward the right half plane
% (zeta falls again). Ignoring the slow breakaway region (|s| < wmin), k is
% the LARGEST multiplier that still gives zeta >= zt (fastest loop meeting
% the damping target). If zt is never reached, k = NaN and zmax / kmax give
% the best damping this PI zero allows: then move the zero (change Ki/Kp).
kv = logspace(-2, 2, 400);  z = nan(size(kv));  w = z;
for j = 1:numel(kv)
    [d, z(j)] = dominant(locus(a, b, Kp, wz, tau, Td, kv(j)));  w(j) = abs(d);
end
z(w < wmin) = NaN;
[zmax, jm] = max(z);  kmax = kv(jm);
j = find(z(1:end-1) >= zt & z(2:end) < zt, 1, 'last');
if isempty(j) || zmax < zt, k = NaN; return; end
k = kv(j) * (kv(j+1)/kv(j))^((z(j) - zt)/(z(j) - z(j+1)));
end

function [Kp, wz] = shape(a, b, tau, Td, wgc, pm)
% loop-shaping PI design (same as roll_gain_design)
thz = max(-pi/2 + pm*pi/180 + atan(wgc/a) + atan(wgc*tau) + wgc*Td, 1e-3);
wz  = wgc / tan(thz);
Kp  = wgc*sqrt(wgc^2 + a^2)*sqrt(1 + (wgc*tau)^2) / (b*sqrt(wgc^2 + wz^2));
end
