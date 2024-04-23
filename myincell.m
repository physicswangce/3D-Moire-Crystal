function  m = myincell(R,V1,V2,V3)
%vectors with size[3,1]
m = 0;
m3 = ((V3 - R)'*cross(V1,V2))*(R'*cross(V1,V2));
m2 = ((V2 - R)'*cross(V3,V1))*(R'*cross(V3,V1));
m1 = ((V1 - R)'*cross(V2,V3))*(R'*cross(V2,V3));

a3 = R'*cross(V1,V2);
a1 = R'*cross(V2,V3);
a2 = R'*cross(V3,V1);

if (a1 == 0 || m1 > 0) && (a2 == 0 || m2 > 0) && (a3 == 0 || m3 > 0)
 m = 1;
end

end

