clear; close all; clc;
%%
saving = false;

addpath('sub_direct')
addpath('trajectory_generator')
addpath('results')

% obj = 'tracking'; % or regulation
% obj = 'tracking2';
% obj = 'regulation2';
obj = 'regulation';
%% initialization
if strcmp(obj,'tracking')
    t_end = 10;
elseif strcmp(obj,'regulation')
    t_end = 5;
elseif strcmp(obj,'tracking2')
    t_end = 5;
elseif strcmp(obj,'regulation2')
    t_end = 15;
end
t = 0 : 0.001 : t_end;
N = length(t);

S = zeros(12,N);
% S(1:6,1) = [0.1721, -1.0447, 1.6729, -0.6282, 0.1721, pi/3];
S(1:6,1) = [0.2, -0.5, 0.4, 0.6, -0.5, 0.2];

if strcmp(obj,'tracking2')
%    S(1:6,1) = [0.1721, -1.0447, 1.6729, -0.6282, 0.1721, 0];

   g_0 = g_st_fun(S(1,1),S(2,1),S(3,1),S(4,1),S(5,1),S(6,1));
   p_i = g_0(1:3,4);
   p_f = [-0.5, 0.3, 0.5]';
   R_i = g_0(1:3,1:3);

   w = [1,1,1]';
   w_norm = norm(w,2);
   w = w/w_norm;
   b = 0;
   
   t_end_traj = 2;
   
   R_f = rodrigues(w,pi-0.01) * R_i;
   [param_x,param_y,param_z,param_w1,param_w2,param_w3] = trajectory_calculator(p_i,p_f,R_i,R_f,t(1),t_end_traj);
   param_p = [param_x;  param_y;  param_z];
   param_w = [param_w1; param_w2; param_w3];

%    S(1:6,1) = [0.4, -pi/6, 0.4, 0.6, -0.5, 0.2];
elseif strcmp(obj,'regulation2')
    S(1:6,1) = [0.1721, -1.0447, 1.6729, -0.6282, 0.1721, 0];
end

% g_0 = g_st_fun(0,0,0,0,0,0);
g_0 = g_st_fun(S(1,1),S(2,1),S(3,1),S(4,1),S(5,1),S(6,1));

T_arr = zeros(6,N);
T1_arr = zeros(6,N);
g_se_arr = zeros(N,4,4);
pd_arr = zeros(3,N);
g_d_arr = zeros(N,4,4);
F_arr = zeros(6,N);
kinetic_arr = zeros(1,N);
log_err_arr = zeros(6,N);

lyap_arr = zeros(1,N);
potential_arr = zeros(1,N);

sing_arr = zeros(1,N);

err_arr = zeros(1,N);
err_rot_arr = zeros(1,N);

ev_arr = zeros(6,N);

r = zeros(3,N);

%% simulation
for k = 1 : N-1
    %% input calculation
    gain_t = 100;
    gain_o = 100;
    
    %%% bakje
%     kt1 = 100; kt2 = 30; kt3 = 40;
%     ko1 = 100; ko2 = 30; ko3 = 40;
%     kt1 = 200; kt2 = 60; kt3 = 80;
%     ko1 = 10; ko2 = 30; ko3 = 100;

%     kt1 = 100; kt2 = 60; kt3 = 80;

%     kt1 = 100; kt2 = 10; kt3 = 20;
    kt1 = gain_t; kt2 = gain_t; kt3 = gain_t;
%     kt1 = 200; kt2 = 60; kt3 = 80;
    ko1 = gain_o; ko2 = gain_o; ko3 = gain_o;
%     ko1 = 10; ko2 = 100; ko3 = 30;

%     kt1 = 200; kt2 = 60; kt3 = 80;
%     ko1 = 10; ko2 = 30; ko3 = 100;
    
    kv = 50;
    
    Kt = diag([kt1,kt2,kt3]);

    Ko = diag([ko1, ko2, ko3]);

    T_max = [150, 100, 100, 100, 100, 100]';
    
    Kg = blkdiag(Kt,Ko);
    Kxi = kv * eye(6);

    Kxi = 5 * sqrt(Kg);
    
    %%
    q1 = S(1,k); q2 = S(2,k); q3 = S(3,k); 
    q4 = S(4,k); q5 = S(5,k); q6 = S(6,k);
    dq1 = S(7,k); dq2 = S(8,k); dq3 = S(9,k); 
    dq4 = S(10,k); dq5 = S(11,k); dq6 = S(12,k);
    
    q = [q1, q2, q3, q4, q5, q6]';
    dq = [dq1, dq2, dq3, dq4, dq5, dq6]';
    %%
    M = M_fun(q1,q2,q3,q4,q5,q6);
    C = C_fun(q1,q2,q3,q4,q5,q6,dq1,dq2,dq3,dq4,dq5,dq6);
    Jb = Jb_fun(q1,q2,q3,q4,q5,q6);
    Jb_dot = Jb_dot_fun(q1,q2,q3,q4,q5,q6,dq1,dq2,dq3,dq4,dq5,dq6);
    G = G_fun(q1,q2,q3,q4,q5,q6);
    g_se = g_st_fun(q1,q2,q3,q4,q5,q6);
    
    g_se_arr(k,:,:) = g_se;
    
    R = g_se(1:3,1:3); p = g_se(1:3,4);
    xi = Jb * dq;
    v = xi(1:3); w = xi(4:6);
    
    %%
    if strcmp(obj,'tracking')
        [Rd,pd,xi_d,dxi_d] = desired_trajectory(t(k),g_0,'geo');
    elseif strcmp(obj,'regulation')
        [Rd,pd,xi_d,dxi_d] = desired_trajectory_regulation(t(k),g_0);
    elseif strcmp(obj,'tracking2')
        [Rd,pd,xi_d,dxi_d] = desired_trajectory2(t(k),g_0,param_p, param_w,'geo',t_end_traj);
    elseif strcmp(obj,'regulation2')
        [Rd,pd,xi_d,dxi_d] = desired_trajectory3(t(k),g_0,t(end),'geo');
    end
    
    vd = xi_d(1:3); wd = xi_d(4:6);
    g_ed = [R'*Rd, R'*(p - pd);
            zeros(1,3),1];
        
    R_ed = g_ed(1:3,1:3);
    p_ed = g_ed(1:3,4);
    
    pd_arr(:,k) = pd;
    g_d_arr(k,:,:) = [Rd, pd;
                      zeros(1,3),1];
         
    g_d = [Rd, pd;
           zeros(1,3),1];
    g_de = inv(g_d) * g_se;
    log_err = log_map_SE3(g_de);

    xi_d_t = Adj_map(g_ed) * xi_d;
    
    f_g = Kg * log_err;
    e_xi = xi - xi_d_t;
    %%
    R_ed_dot = -hat_map(w)*R'*Rd + R'*Rd*hat_map(wd);
    p_ed_dot = -hat_map(w)*R'*(p - pd) + v - R'*Rd*vd;
    Ad_ged_dot = [R_ed_dot, hat_map(p_ed_dot)*R_ed + hat_map(p_ed)*R_ed_dot;
                  zeros(3,3), R_ed_dot];
    
    P = Ad_ged_dot * xi_d + Adj_map(g_ed) * dxi_d;
    %%
    Jb_inv = inv(M) * Jb' * pinv(Jb * inv(M) * Jb'); 
    M_tilde = Jb_inv' * M * Jb_inv;
%     M_tilde = pinv(Jb)' * M * pinv(Jb); 
    C_tilde = (Jb_inv)' * (C - M * Jb_inv * Jb_dot) * Jb_inv; 
    
    if abs(det(Jb)) < 0.01
        sing_arr(k) = 1;
        F = - f_g - Kxi * e_xi;
    else
        F = C_tilde * xi_d_t -  f_g - Kxi * e_xi + M_tilde * P;
    end
    
    T1 = Jb'* F;
    T = T1 + G';
    T_arr(:,k) = T;
    T1_arr(:,k) = T1;
    F_arr(:,k) = F;
    ev_arr(:,k) = dq;
    log_err_arr(:,k) = log_err;
    
    initVal = S(:,k);
    [Time,S_] = ode15s(@(Time,S_) robot_dynamics(Time,S_,T), [t(k),t(k+1)], initVal ,' ');
    [Nn,~] = size(S_);
    S(:,k+1) = S_(Nn,:);
    
    
    potential = potential_fun(g_se, g_d_arr(k,:,:), Kt,Ko);
    kinetic = 1/2 * e_xi' * M_tilde * e_xi;
    kinetic_arr(k) = 1/2 * dq' * M * dq;

    if (strcmp(obj,'regulation') || strcmp(obj,'regulation2'))
        kinetic = 1/2 * (dq)' * dq;
    end

    err_arr(k) = err_fun(g_se_arr(k,:,:),g_d_arr(k,:,:));
    err_rot_arr(k) = err_fun_rot(g_se_arr(k,:,:),g_d_arr(k,:,:)); 

    lyap_arr(k) = potential + kinetic;
    potential_arr(k) = potential;
    
    r(1,k) = 1 - Rd(:,1)' * R(:,1);
    r(2,k) = 1 - Rd(:,2)' * R(:,2);
    r(3,k) = 1 - Rd(:,3)' * R(:,3);
    
    if mod(k,1000) == 0
        disp(k/1000)
    end
end

%%
t_ = t(1:end-1);

figure
subplot(3,1,1);
plot(t_,g_se_arr(1:end-1,1,4)); hold on; plot(t_,pd_arr(1,1:end-1))
subplot(3,1,2);
plot(t_,g_se_arr(1:end-1,2,4)); hold on; plot(t_,pd_arr(2,1:end-1))
subplot(3,1,3);
plot(t_,g_se_arr(1:end-1,3,4)); hold on; plot(t_,pd_arr(3,1:end-1))

%%
p = zeros(3,N);
R = zeros(3,3,N);
for k = 1 : N
    p(:,k) = g_se_arr(k,1:3,4);
    R(:,:,k) = g_se_arr(k,1:3,1:3);
end

figure
plot3(g_se_arr(1:end-1,1,4),g_se_arr(1:end-1,2,4),g_se_arr(1:end-1,3,4)); hold on; grid on;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

nk = 500;
for k = 1 : N/nk
    idx = nk*(k-1) + 1;
    scale = 0.05;
    plot3([p(1,idx), p(1,idx) + scale * R(1,1,idx)],...
          [p(2,idx), p(2,idx) + scale * R(2,1,idx)],...
          [p(3,idx), p(3,idx) + scale * R(3,1,idx)],'r');
    plot3([p(1,idx), p(1,idx) + scale * R(1,2,idx)],...
          [p(2,idx), p(2,idx) + scale * R(2,2,idx)],...
          [p(3,idx), p(3,idx) + scale * R(3,2,idx)],'b');
    plot3([p(1,idx), p(1,idx) + scale * R(1,3,idx)],...
          [p(2,idx), p(2,idx) + scale * R(2,3,idx)],...
          [p(3,idx), p(3,idx) + scale * R(3,3,idx)],'g');
end
axis equal
%%
if saving == true
result_log.S = S;
result_log.g = g_se_arr;
result_log.g_d = g_d_arr;
result_log.T = T;
result_log.T1 = T1;
result_log.t = t;
result_log.lyap = lyap_arr;
result_log.dq = ev_arr;
result_log.potential = potential_arr;

    if strcmp(obj,"tracking")
        save('results/result_log_tracking.mat','result_log');
    elseif strcmp(obj,'regulation')
        save('results/result_log_regulation.mat','result_log');
    elseif strcmp(obj,'tracking2')
        save('results/result_log_tracking2.mat','result_log');
    elseif strcmp(obj,'regulation2')
        save('results/result_log_regulation2.mat','result_log');
    end
end
%%
T_rms = rms(T,2);
T1_rms = rms(T1,2);
F_rms = rms(F,2);
pos_err_rms = rms(p - pd_arr,2);
lyap_rms = rms(lyap_arr);

err_rms = rms(err_arr);
err_rot_rms = rms(err_rot_arr);

disp('Log T rms')
disp(T_rms')
disp(F_rms');

disp(norm(T,'fro')); disp(norm(T1,'fro'));

fprintf('Log pos error rms %f \n', pos_err_rms);
fprintf('lyap rms %f \n', lyap_rms);
fprintf('err rms %f\n', err_rms);
fprintf('err rot rms %f \n', err_rot_rms);

%%
figure
plot(lyap_arr);

%%
figure
plot(r(1,:)); hold on; plot(r(2,:)); plot(r(3,:),'-');