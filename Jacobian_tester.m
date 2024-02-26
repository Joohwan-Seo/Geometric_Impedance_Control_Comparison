clear; close all; clc;

addpath('sub_direct');

q = randn(6,1);

q1 = q(1); q2 = q(2); q3 = q(3);
q4 = q(4); q5 = q(5); q6 = q(6);

Jb = Jb_fun(q1,q2,q3,q4,q5,q6);
Je = Je_fun(q1,q2,q3,q4,q5,q6);

g = g_st_fun(q1,q2,q3,q4,q5,q6);

Jv = Je(1:3,:); Jw = Je(4:6,:);

R = g(1:3,1:3);

Jb_recal = [R' * Jv; R' * Jw];

Jb_recal - Jb