clear; close all; clc;

addpath('results')
addpath('sub_direct')
addpath('trajectory_generator')

% obj = 'tracking'; % or regulation
% obj = 'tracking2';
% obj = 'regulation2';
obj = 'regulation';
%%
if strcmp(obj,'regulation')
    load('result_log_regulation.mat');
    load('result_geo2_regulation.mat');
elseif strcmp(obj,'tracking')
    load('result_log_tracking.mat');
    load('result_geo2_tracking.mat');
elseif strcmp(obj,'tracking2')
    load('result_log_tracking2.mat');
    load('result_geo2_tracking2.mat');
elseif strcmp(obj,'regulation2')
    load('result_log_regulation2.mat');
    load('result_geo2_regulation2.mat');
end

%downsample data by factor
down = 10; % default 1

t = result_log.t;

N = length(t);
t_ = t(1:end-1);

gd = result_log.g_d;

if down > 1
    gd = downsample(gd, down);
    N = ceil(N/down);
    t = downsample(t, down);
    t_ = downsample(t_, down);
    result_log.S = downsample(result_log.S', down)';
    result_log.g = downsample(result_log.g, down);
    result_log.g_d = downsample(result_log.g_d, down);
    result_log.t = downsample(result_log.t, down);
    result_geo2.S = downsample(result_geo2.S', down)';
    result_geo2.g = downsample(result_geo2.g, down);
    result_geo2.g_d = downsample(result_geo2.g_d, down);
    result_geo2.t = downsample(result_geo2.t, down);
    
    result_log.lyap = downsample(result_log.lyap, down);
    result_geo2.lyap = downsample(result_geo2.lyap, down);
end

%%
p_log = zeros(3,N-1);
R_log = zeros(3,3,N-1);
p_geo = zeros(3,N-1);
R_geo = zeros(3,3,N-1);
p_des = zeros(3,N-1);

for k = 1 : N-1
    p_log(:,k) = result_log.g(k,1:3,4);
    R_log(:,:,k) = result_log.g(k,1:3,1:3);
    p_geo(:,k) = result_geo2.g(k,1:3,4);
    R_geo(:,:,k) = result_geo2.g(k,1:3,1:3);
    p_des(:,k) = result_log.g_d(k,1:3,4);
end

h = 1:N-1/down;

down2 = 20;

figure(1)
plot3(p_log(1,h),p_log(2,h),p_log(3,h),'k'); hold on; grid on;
plot3(p_geo(1,h),p_geo(2,h),p_geo(3,h),'m');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
nk = 500/down2; 
% legend('log','geo')

for k = 1 : length(h)/nk
    idx = nk*(k-1) + 1;
    scale = 0.05;
    plot3([p_log(1,idx), p_log(1,idx) + scale * R_log(1,1,idx)],...
          [p_log(2,idx), p_log(2,idx) + scale * R_log(2,1,idx)],...
          [p_log(3,idx), p_log(3,idx) + scale * R_log(3,1,idx)],'r');
    plot3([p_log(1,idx), p_log(1,idx) + scale * R_log(1,2,idx)],...
          [p_log(2,idx), p_log(2,idx) + scale * R_log(2,2,idx)],...
          [p_log(3,idx), p_log(3,idx) + scale * R_log(3,2,idx)],'b');
    plot3([p_log(1,idx), p_log(1,idx) + scale * R_log(1,3,idx)],...
          [p_log(2,idx), p_log(2,idx) + scale * R_log(2,3,idx)],...
          [p_log(3,idx), p_log(3,idx) + scale * R_log(3,3,idx)],'g');
      
    plot3([p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,1,idx)], ...
          [p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,1,idx)], ...
          [p_geo(3,idx), p_geo(3,idx) + scale * R_geo(3,1,idx)],'r');
    plot3([p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,2,idx)], ...
          [p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,2,idx)], ...
          [p_geo(3,idx), p_geo(3,idx) + scale * R_geo(3,2,idx)],'b');
    plot3([p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,3,idx)], ...
          [p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,3,idx)], ...
          [p_geo(3,idx), p_geo(3,idx) + scale * R_geo(3,3,idx)],'g');
end
axis equal
%%
figure(2)
subplot(3,1,1);
plot(t_,gd(1:end-1,1,4),'k:'); hold on; grid on; ylabel('x (m)');
plot(t_,result_geo2.g(1:end-1,1,4),'r');
plot(t_,result_log.g(1:end-1,1,4),'b--'); 
legend('desired','geo','log')
subplot(3,1,2);
plot(t_,gd(1:end-1,2,4),'k:'); hold on; grid on; ylabel('y (m)');
plot(t_,result_geo2.g(1:end-1,2,4),'r');
plot(t_,result_log.g(1:end-1,2,4),'b--'); 
subplot(3,1,3);
plot(t_,gd(1:end-1,3,4),'k:'); hold on; grid on; 
plot(t_,result_geo2.g(1:end-1,3,4),'r');
plot(t_,result_log.g(1:end-1,3,4),'b--'); 
ylabel('z (m)'); xlabel('t (s)');

%%
err_log = zeros(1,N);
err_geo = zeros(1,N);

err_log_rot = zeros(1,N);
err_geo_rot = zeros(1,N);

for k = 1 : N-1
    err_log(k) = err_fun(result_log.g(k,:,:),result_log.g_d(k,:,:));
    err_geo(k) = err_fun(result_geo2.g(k,:,:),result_geo2.g_d(k,:,:));
    
    err_log_rot(k) = err_fun_rot(result_log.g(k,:,:),result_log.g_d(k,:,:));
    err_geo_rot(k) = err_fun_rot(result_geo2.g(k,:,:),result_geo2.g_d(k,:,:));
end

%%
figure(3)
plot(t_,err_geo(1:end-1),'r'); hold on; grid on;
plot(t_,err_log(1:end-1),'b--'); 
xlabel('t (s)'); ylabel('Error function \Psi')
legend('geo','log')

figure(4)
plot(t_,result_geo2.lyap(1:end-1),'r'); hold on; grid on;
plot(t_,result_log.lyap(1:end-1),'b--'); 
xlabel('t (s)'); ylabel('Dynamic Cost \Phi');
legend('geo','log')
% ylim([0,1])
%%
RMS_p_geo = rms(p_log - p_des,2);
RMS_p_imp = rms(p_geo - p_des,2);

fprintf('LOG P RMS: %f, %f, %f\n',RMS_p_geo(1),RMS_p_geo(2),RMS_p_geo(3))
fprintf('GIC P RMS: %f, %f, %f\n',RMS_p_imp(1),RMS_p_imp(2),RMS_p_imp(3))

fprintf('LOG ERR RMS: %f, ROT RMS:%f\n', rms(err_log), rms(err_log_rot));
fprintf('GIC ERR RMS: %f, ROT RMS:%f\n', rms(err_geo), rms(err_geo_rot));

fprintf('LOG Dynamic cost RMS: %f\n', rms(result_log.lyap));
fprintf('GIC Dynamic cost RMS: %f\n', rms(result_geo2.lyap));

%% Testing for rightful rotation matrix
det_log = zeros(1,N-1);
det_geo = zeros(1,N-1);

for k = 1 : N-1
    det_log(k) = det(R_geo(:,:,k));
    det_geo(k) = det(R_log(:,:,k));
end

