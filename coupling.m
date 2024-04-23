function  omega = coupling(r1,r2,Vx,Vy,Vz,K,eta,beta)
omega = 0;

R1 = r1;

for qx = -3 : 1 : 3
for qy = -3 : 1 : 3
for qz = -3 : 1 : 3
shift =   qx*Vx + qy*Vy + qz*Vz;
R2 = r2 + shift;
D = (R1(1)-R2(1))^2 + (R1(2)-R2(2))^2+(R1(3)-R2(3))^2;
omega = omega + beta*exp(-D*sqrt(eta)*pi^2/4)*exp(-1i*(shift'*K));
end
end
end
end








