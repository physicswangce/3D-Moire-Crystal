function  plist = cellposition(Vx,Vy,Vz)
V000 = [0;0;0];
V100 = Vx; V010 = Vy; V001 = Vz;
V110 = Vx+Vy; V011 = Vy+Vz; V101 = Vx+Vz;
V111 = Vx+Vy+Vz;

V = [V000,V100,V010,V001,V110,V101,V011,V111];
Vmax_x = max(V(1,:));
Vmax_y = max(V(2,:));
Vmax_z = max(V(3,:));
Vmin_x = min(V(1,:));
Vmin_y = min(V(2,:));
Vmin_z = min(V(3,:));
%unit cell 
plist = zeros(3,1);
k = 1;
for x = Vmin_x : 1 : Vmax_x
    for y = Vmin_y : 1 : Vmax_y
        for z = Vmin_z : 1 : Vmax_z
            ifin = myincell([x;y;z],Vx,Vy,Vz);
            if ifin == 1 
               plist(:,k)=[x;y;z];
                k = k + 1;
            end
            
        end
    end
end


end

