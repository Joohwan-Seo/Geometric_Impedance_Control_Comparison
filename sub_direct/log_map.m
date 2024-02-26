function w = log_map(R)
phi = acos((trace(R) - 1)/2);

if abs(phi) < 0.001
    log_R = phi / (2*sin(phi) + 0.0001) * (R - transpose(R));
else
    log_R = phi / (2*sin(phi)) * (R - transpose(R));
end
w = [-log_R(2,3), log_R(1,3), -log_R(1,2)]';