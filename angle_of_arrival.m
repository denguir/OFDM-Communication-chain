function [AoA] = angle_of_arrival(h_mat, phi, theta, beta, dist)
%ANGLE_OF_ARRIVAL Summary of this function goes here
%   Detailed explanation goes here

[xx, yy, zz, tap_max] = size(h_mat);
AoA = zeros(length(phi),length(theta),tap_max);

for p=1:length(phi)
    for t=1:length(theta)
        for X=0:xx-1
            for Y=yy-1:-1:0
                for Z=0:zz-1
                    BETA = beta*[sin(theta(t))*cos(phi(p)),sin(theta(t))*sin(phi(p)),cos(theta(t))];
                    R = [X Y Z] * dist;
                    Bi = exp(-1i*dot(BETA, R));
                    for n=1:tap_max
                        hn = h_mat(10-X, Y+1, 10-Z, n);
                        AoA(p,t,n) = AoA(p,t,n) + hn*conj(Bi);
                    end
                end
            end
        end
    end
end
AoA = AoA/(xx*yy*zz);

end

