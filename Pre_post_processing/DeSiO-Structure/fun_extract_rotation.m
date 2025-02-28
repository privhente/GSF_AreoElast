% script for calculating rotor speed from director-based formulation
function [phi] = fun_extract_rotation(q,t,node) 

% indizes for node
inz = [12*(node-1)+1:12*(node-1)+12];

d1_0 = q(1,inz(4:6));
d2_0 = q(1,inz(7:9));
d3_0 = q(1,inz(10:12));

phi = zeros(length(t),3);
for j = 2:length(t)
    d1_n = q(j,inz(4:6));
    d2_n = q(j,inz(7:9));
    d3_n = q(j,inz(10:12));

    delta_d1 = d1_n - d1_0;
    delta_d2 = d2_n - d2_0;
    delta_d3 = d3_n - d3_0;

    phi(j,1:3) = phi(j-1,1:3) + 0.5*( cross(d1_n,delta_d1) + cross(d2_n,delta_d2) + cross(d3_n,delta_d3) );

    d1_0 = d1_n;
    d2_0 = d2_n;
    d3_0 = d3_n;
end
return

