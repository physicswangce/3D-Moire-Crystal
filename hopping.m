function t = hopping(r1,r2,Vx,Vy,Vz,K,eta)
% K size [3,1]
t = 0;
t_ = -(4/sqrt(pi))*eta^(3/4)*exp(-2*sqrt(eta));
R1 = r1;
for qx = -3 : 1 : 3
for qy = -3 : 1 : 3
for qz = -3 : 1 : 3
shift =   qx*Vx + qy*Vy + qz*Vz;
R2 = r2 + shift;
if isneighbor(R1,R2) == 1
t = t + t_*exp(-1i*(shift'*K));
 end
end
end
end




end

