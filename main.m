% This code generate the tight-binding Hamiltonian for Morie lattice
% As an example, we consider the lattice for m = 3 , n = 1, lx = ly = lz, 
% a \pi/3 rotation along (1,1,1) axis.


%The first step is to determine the unit cell with the mathematica code 
Vx =[-1;1;0];
Vy = [1;0;-1];
Vz = [1;1;1];
%The second step is to determine the position of all the sites in the unit cell
% as Alist_ and Blist_.  (rs is the displacement)
m = 3; n = 1;
lx = 1; ly = 1; lz = 1; ll = [lx;ly;lz];
D = lx^2 + ly^2 +lz^2;
c = (m^2 - D*n^2)/(m^2 + D*n^2); 
s = 2*n*m*sqrt(D)/(m^2 + D*n^2);
ux = lx/sqrt(D); uy = ly/sqrt(D);  uz = lz/sqrt(D); 
S = [ c + ux^2 * ( 1- c), ux*uy*(1-c)-uz*s, ux*uz*(1-c) + uy*s;
      ux*uy*(1-c) + uz *s, c + uy^2 *(1-c), uy*uz*(1-c) - ux*s;
      ux*uz*(1-c) - uy *s, uy*uz*(1-c) + ux*s, c + uz^2*(1-c) ];

Ux = round(S'*Vx);
Uy = round(S'*Vy);
Uz = round(S'*Vz);
rs = (0.1) * [0;-1;1];
Rs = rs * ones(1,3);
Alist = cellposition(Vx,Vy,Vz);
Alist_ = Alist + Rs;
Blist = cellposition(Ux,Uy,Uz);
Blist_ = S*Blist;
%The third step is to determine the high symmetry points in reciprocal
%space, and a path to present the band structure according to those high
%symmetry points.

%primitive reciprocal lattice vectors
absV = abs(det([Vx,Vy,Vz]));
Gx = 2*pi*cross(Vy,Vz)/absV;
Gy = 2*pi*cross(Vz,Vx)/absV;
Gz = 2*pi*cross(Vx,Vy)/absV;

%High symmetry points hexagonal crystal system
G_gamma = [0;0;0];
G_M = 0.5*Gx;
G_K = Gx/3  + Gy/3;
G_A = Gz/2;
G_L = Gx/2 + Gz/2;
G_H = Gx/3 + Gy/3 + Gz/2;

nk = 50;

%Path  gamma - M - K - gamma  - A - L - H - A
Gpath{1} = G_gamma; Gpath{2} = G_M; Gpath{3} = G_K;
Gpath{4} = G_gamma; Gpath{5} = G_A; Gpath{6} = G_L;
Gpath{7} = G_H;     Gpath{8} = G_A;
Glength = zeros(7,1);
for p = 1 : 7
   Glength(p) =  sqrt((Gpath{p} - Gpath{p+1})'*(Gpath{p} - Gpath{p+1}));
end
nklist = round(Glength*100);
position = zeros(8,1);
position(1) = 0;
for p = 1 : 7
    position(1+p) = position(p) + nklist(p);
end

Kpath = zeros(3,sum(nklist));
ns = 0;
for p = 1 : 7
    nkplist  = linspace(0,1,nklist(p));
for k = 1 : nklist(p)
Kpath(:,k+ns) = Gpath{p}*(1-nkplist(k))+Gpath{p+1}*nkplist(k);
end
   ns = ns + nklist(p);
end





%The fourth step is to find the band structure on the path within the tight-binding
% approximation, for given system parameters, eta for 'V / Er', beta for 'Omega / Er' 
% , delta for 'delta /Er ' .
eta = 6; beta = 1;  delta = 0.00;
t_ = -(4/sqrt(pi))*eta^(3/4)*exp(-2*sqrt(eta)); % hopping strength 
Dim = size(Alist_,2);
nb = Dim*2;
Specturm = zeros(size(Kpath,2),nb);
for k = 1 : size(Kpath,2)
K = Kpath(:,k);
H = Hamiltonian(Alist_,Blist_,Vx,Vy,Vz,eta,beta,delta,K,Dim);
[uu,vv]=eig(H);
Specturm(k,:)= sort(diag(vv));
end


%plot the band structure
plot_band



