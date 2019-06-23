function [h_simo] = get_simo_ch(AoA, index_simo, phi, theta, h_mat_rect, beta, dist)

[xx, yy, zz, tap_max] = size(h_mat_rect);
h_simo = zeros(xx, yy, zz, tap_max);
for j=1:tap_max
    for X=0:xx-1
        for Y=yy-1:-1:0
            for Z=0:zz-1
                R = [X Y Z] * dist;
                for i=1:size(index_simo{1,j},2)
                    idx1 = index_simo{1,j}(1,i); % phi index
                    idx2 = index_simo{1,j}(2,i); % theta index
                    BETA = beta*[sin(theta(idx2))*cos(phi(idx1)),sin(theta(idx2))*sin(phi(idx1)),cos(theta(idx2))];
                    angle_rand = random('unif',0,2*pi);
                    expo = exp(1i*(angle_rand-dot(BETA, R)));
                    h_simo(X+1,Y+1,Z+1,j) = h_simo(X+1,Y+1,Z+1,j) + AoA(idx1,idx2,j)*expo;
                end
            end
        end
    end
end

end

