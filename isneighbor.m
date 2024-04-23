function m = isneighbor(r1,r2)
m = 0; 
c = (r1(1)-r2(1))^2 + (r1(2)-r2(2))^2 + (r1(3)-r2(3))^2;

if abs(c-1)<0.001
 m = 1;
end

end

