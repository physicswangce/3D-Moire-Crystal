function H = Hamiltonian(Alist_,Blist_,Vx,Vy,Vz,eta,beta,delta,K,Dim)

H = zeros(2*Dim,2*Dim);
for k = 1 : Dim
  for j = 1 : Dim
    if j == k 
      H(k,j) =  - delta;
    
    else
      H(k,j) = hopping(Alist_(:,k),Alist_(:,j),Vx,Vy,Vz,K,eta);
    end
  end
end
 
for k = 1 : Dim
  for j = 1 : Dim
    if j == k 
      H(k+Dim,j+Dim) =   delta;
    
    else
      H(k+Dim,j+Dim) = hopping(Blist_(:,k),Blist_(:,j),Vx,Vy,Vz,K,eta);
    end
  end
end

for k = 1 : Dim
  for j = 1 : Dim
    
      H(k,j+Dim) = coupling(Alist_(:,k),Blist_(:,j),Vx,Vy,Vz,K,eta,beta);
    
  end
end

for k = 1 : Dim
  for j = 1 : Dim
    
      H(k+Dim,j) = coupling(Blist_(:,k),Alist_(:,j),Vx,Vy,Vz,K,eta,beta);
    
  end
end
H = (H + H')/2;


end

