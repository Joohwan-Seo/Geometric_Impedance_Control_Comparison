clear; close all; clc;
addpath('sub_direct')
%%
R_d = eye(3);
N = 100;
t_e = linspace(-pi,pi,N);

w_log = zeros(1,N);
w_norm = zeros(1,N);
phi_list = zeros(1,N);
w_list = zeros(3,N);
mat_tr = zeros(1,N);

w_frob = zeros(1,N);
w_vec_fro = zeros(3,N);

for k = 1 : N
    R_e = [cos(t_e(k)), -sin(t_e(k)), 0;
           sin(t_e(k)), cos(t_e(k)), 0;
           0, 0, 1];
    R_de = transpose(R_d) * R_e;   
    phi = acos((trace(R_de) - 1)/2);
    
    mat_tr(k) = trace(R_de);
    w_frob(k) = norm(eye(3) - R_de,'fro');

    w = log_map(R_de);
    w_log(k) = w(3);
    w_norm(k) = norm(w,2);
    phi_list(k) = phi;
    w_list(:,k) = w;
    
    err = (R_de - transpose(R_de));
    w_vec_fro(:,k) = [-err(2,3), err(1,3), -err(1,2)]';
end
%%
% V_twist = w_norm.^2;
V_twist = 0.5 * w_log.^2;
V_frob = 0.5 * w_frob.^2;
%%
figure(1)
plot(t_e, V_frob, 'r'); hold on; grid on;
plot(t_e, V_twist, 'k--'); 
xlim([-3.15, 3.15])
ax = gca;
ax.FontSize = 13;
xlabel('$\theta$ (rad)', 'Interpreter', 'latex', 'FontSize', 13); 
ylabel('Error Functions on SO(2)','Fontsize', 13);
legend('$\Psi_1$','$\Psi_2$', 'Interpreter', 'latex', 'FontSize', 12)

figure(2)
plot(t_e, w_vec_fro(3,:),'r'); hold on; grid on;
plot(t_e, w_log, 'k--'); 
xlim([-3.15, 3.15])
ylim([-3.15, 3.15])
ax = gca;
ax.FontSize = 13;
xlabel('$\theta$ (rad)', 'Interpreter', 'latex', 'FontSize', 13); 
ylabel('$f_{_{G,i}}$', 'Interpreter', 'latex', 'FontSize', 13);
legend('$f_{_{G,1}}$', '$f_{_{G,2}}$', 'Interpreter', 'latex', 'FontSize', 12);

figure(3)
subplot(1,2,1)
plot(t_e, V_frob, 'r'); hold on; grid on;
plot(t_e, V_twist, 'k--'); 
subplot(1,2,2);
plot(t_e, w_vec_fro(3,:),'r'); hold on; grid on;
plot(t_e, w_log, 'k--'); 
