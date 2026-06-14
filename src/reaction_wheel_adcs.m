function reaction_wheel_adcs()
% REACTION_WHEEL_ADCS  Quaternion-feedback attitude control of a rigid
% spacecraft with three reaction wheels.
%
%   Personal project (Ali Murtaza). A rigid satellite performs a large-
%   angle (60 deg) rest-to-rest eigenaxis slew under quaternion-feedback
%   PD control. Control torque is realised by three body-axis reaction
%   wheels with torque and momentum limits; a small constant external
%   disturbance drives a secular wheel-momentum build-up, illustrating
%   why momentum management (desaturation) is needed.
%
%   Quaternion convention: q = [q1 q2 q3 q4], q4 scalar. No toolboxes
%   required.

clc; close all;
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, '..', 'assets');
if ~exist(outdir,'dir'); mkdir(outdir); end

%% ------------------------------------------------------------- Plant
J  = diag([10, 12, 8]);                 % spacecraft inertia [kg m^2]
Jinv = inv(J);
u_max = 0.5;                            % per-wheel torque limit [N m]
h_max = 5.0;                            % per-wheel momentum limit [N m s]
tau_d = [5e-3; -2e-3; 1e-3];           % constant external disturbance [N m]

%% -------------------------------------------------------- Controller
% Quaternion PID:  u = Kp*qe_vec + Ki*int(qe_vec) - Kd*omega
% The integral term cancels the steady-state pointing bias produced by the
% constant external disturbance (pure PD would leave Kp*qe = -tau_d offset).
Kp = 0.80; Ki = 0.040; Kd = 7.2;       % near-critical damping (zeta ~ 0.9)
ei = [0;0;0]; ei_lim = 0.5;            % integral state + anti-windup clamp
i_gate = deg2rad(5);                    % conditional integration band

%% ---------------------------------------------------------- Scenario
dt = 0.02; T = 120; t = (0:dt:T).'; N = numel(t);
n_eig = [1;1;1]/sqrt(3);                % slew eigenaxis
ang   = deg2rad(60);                    % slew angle
q_cmd = [n_eig*sin(ang/2); cos(ang/2)];% commanded attitude

q = [0;0;0;1];                         % initial attitude (identity)
w = [0;0;0];                           % initial body rate
h = [0;0;0];                           % wheel momenta

Q=zeros(N,4); W=zeros(N,3); H=zeros(N,3); U=zeros(N,3); ERR=zeros(N,1);

for k = 1:N
    % --- Attitude error (current -> commanded)
    qe = quatmul(quatconj(q), q_cmd);
    if qe(4) < 0; qe = -qe; end          % shortest path
    ERR(k) = 2*acos( min(1,abs(qe(4))) );% eigenaxis error angle [rad]

    % --- Control law (PID) with conditional integration (anti-windup):
    %     the integrator only acts once the slew is nearly complete, so it
    %     never winds up during the large-angle transient.
    if ERR(k) < i_gate
        ei = ei + qe(1:3)*dt;
        ei = max(min(ei, ei_lim), -ei_lim);
    end
    u = Kp*qe(1:3) + Ki*ei - Kd*w;       % drives qe_vec -> 0
    u = max(min(u, u_max), -u_max);
    % wheel momentum saturation: if a wheel is saturated, clip its torque
    for i=1:3
        if (h(i)>=h_max && -u(i)>0) || (h(i)<=-h_max && -u(i)<0)
            u(i)=0;
        end
    end

    % --- Log
    Q(k,:)=q.'; W(k,:)=w.'; H(k,:)=h.'; U(k,:)=u.';

    % --- Dynamics: J*wdot = -w x (J*w + h) + u + tau_d ;  hdot = -u
    f = @(qq,ww,hh) deal( 0.5*qmatrix(ww)*qq, ...
                          Jinv*(-cross(ww, J*ww + hh) + u + tau_d), ...
                          -u );
    [qd1,wd1,hd1]=f(q,w,h);
    [qd2,wd2,hd2]=f(q+0.5*dt*qd1, w+0.5*dt*wd1, h+0.5*dt*hd1);
    [qd3,wd3,hd3]=f(q+0.5*dt*qd2, w+0.5*dt*wd2, h+0.5*dt*hd2);
    [qd4,wd4,hd4]=f(q+dt*qd3,     w+dt*wd3,     h+dt*hd3);
    q = q + (dt/6)*(qd1+2*qd2+2*qd3+qd4); q = q/norm(q);
    w = w + (dt/6)*(wd1+2*wd2+2*wd3+wd4);
    h = h + (dt/6)*(hd1+2*hd2+2*hd3+hd4);
end

%% --------------------------------------------------------- Metrics
errd = rad2deg(ERR);
first_below = find(errd<1, 1, 'first');
max_post = max(errd(first_below:end));
fprintf('\n=== Reaction-Wheel Satellite ADCS ===\n');
fprintf('  Commanded slew:            %.0f deg about [1 1 1]/sqrt(3)\n', rad2deg(ang));
fprintf('  Slew settling (first <1 deg):  %.1f s\n', t(first_below));
fprintf('  Max post-slew excursion:       %.2f deg (disturbance transient)\n', max_post);
fprintf('  Final pointing error:          %.3f deg\n', errd(end));
fprintf('  Peak body rate:            %.2f deg/s\n', max(rad2deg(vecnorm(W,2,2))));
fprintf('  Wheel momentum at T:       [%.2f %.2f %.2f] N m s (limit %.1f)\n', H(end,:), h_max);
fprintf('  Secular momentum rate:     ~%.1e N m s / s (from %.0e N m disturbance)\n', ...
        norm(H(end,:)-H(round(N/2),:))/(T/2), norm(tau_d));

%% --------------------------------------------------------------- Plots
co=[.85 .1 .1; .1 .45 .85; .1 .65 .30];

f1=figure('Color','w','Position',[80 80 920 640]);
subplot(2,1,1); hold on; grid on; box on;
plot(t,errd,'Color',[.2 .2 .2],'LineWidth',2);
yline(1,'--k'); xlabel('time [s]'); ylabel('eigenaxis error [deg]');
title('Large-angle slew: attitude error to commanded 60 deg','FontWeight','bold');
legend('pointing error','1 deg threshold','Location','northeast');
subplot(2,1,2); hold on; grid on; box on;
for i=1:3; plot(t,rad2deg(W(:,i)),'Color',co(i,:),'LineWidth',1.6); end
xlabel('time [s]'); ylabel('body rate [deg/s]');
title('Body angular rates','FontWeight','bold');
legend('\omega_x','\omega_y','\omega_z','Location','northeast');
exportgraphics(f1, fullfile(outdir,'01_slew_response.png'),'Resolution',200);

f2=figure('Color','w','Position',[80 80 920 420]); hold on; grid on; box on;
plot(t,Q(:,1),'Color',co(1,:),'LineWidth',1.5);
plot(t,Q(:,2),'Color',co(2,:),'LineWidth',1.5);
plot(t,Q(:,3),'Color',co(3,:),'LineWidth',1.5);
plot(t,Q(:,4),'Color',[.2 .2 .2],'LineWidth',1.5);
yline(q_cmd(1),':','Color',co(1,:)); yline(q_cmd(4),':','Color',[.2 .2 .2]);
xlabel('time [s]'); ylabel('quaternion');
title('Attitude quaternion converging to command','FontWeight','bold');
legend('q_1','q_2','q_3','q_4','Location','east'); xlim([0 T]);
exportgraphics(f2, fullfile(outdir,'02_quaternion.png'),'Resolution',200);

f3=figure('Color','w','Position',[80 80 920 640]);
subplot(2,1,1); hold on; grid on; box on;
for i=1:3; plot(t,H(:,i),'Color',co(i,:),'LineWidth',1.6); end
yline(h_max,'--k'); yline(-h_max,'--k');
ylabel('wheel momentum [N m s]');
title('Reaction-wheel momentum - secular build-up under disturbance','FontWeight','bold');
legend('h_x','h_y','h_z','limit','Location','northwest');
subplot(2,1,2); hold on; grid on; box on;
for i=1:3; plot(t,U(:,i),'Color',co(i,:),'LineWidth',1.4); end
yline(u_max,':k'); yline(-u_max,':k');
xlabel('time [s]'); ylabel('wheel torque [N m]');
title('Commanded wheel torques (with saturation)','FontWeight','bold');
legend('u_x','u_y','u_z','Location','northeast');
exportgraphics(f3, fullfile(outdir,'03_wheel_momentum.png'),'Resolution',200);

save(fullfile(outdir,'adcs_results.mat'),'t','Q','W','H','U','ERR','q_cmd');
fprintf('\nFigures and results written to %s\n', outdir);

% Export quaternion trajectory (decimated) for the web 3D viewer
dec = 1:5:N;
writematrix([t(dec) Q(dec,:)], fullfile(outdir,'adcs_quaternion.csv'));
end

% --------------------------------------------------------------- locals
function r = quatmul(a,b)   % Hamilton product, scalar-last
av=a(1:3); bv=b(1:3); as=a(4); bs=b(4);
r=[as*bv + bs*av + cross(av,bv); as*bs - dot(av,bv)];
end
function c = quatconj(q); c=[-q(1:3); q(4)]; end
function M = qmatrix(w)     % qdot = 0.5*M*q for body rate w (scalar-last)
M=[ 0     w(3) -w(2) w(1);
   -w(3)  0     w(1) w(2);
    w(2) -w(1)  0    w(3);
   -w(1) -w(2) -w(3) 0   ];
end
