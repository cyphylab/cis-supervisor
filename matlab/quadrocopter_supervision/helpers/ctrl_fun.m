function u = ctrl_fun(K, err)
    u = zeros(3,1);
   
    u_p = [K(1) * err(1:3)];
    u_v = [K(2) * err(4:6)];
    u_a = [K(3) * err(7:9)];

    for (i = 1:3) 
       u(i) = u_p(i) + u_v(i) + u_a(i);
    end
end