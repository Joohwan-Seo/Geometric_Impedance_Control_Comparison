function xi = log_map_SE3(g)
    R = g(1:3,1:3);
    p = g(1:3,4);
    
    w = log_map(R);
    w_hat = hat_map(w);
    w_norm = norm(w,2);
    
    A_inv = eye(3) - 1/2 * w_hat ...
          + (2 * sin(w_norm) - w_norm * (1 + cos(w_norm)))/(2*w_norm^2*sin(w_norm) + 0.001) * w_hat^2;
      
    v = A_inv * p;
    
    xi = [v;w];
end