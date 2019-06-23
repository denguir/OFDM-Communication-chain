function [h] = get_impulse_rsp2(Data, xx, yy, zz, fbins, window)
%GET_IMPULSE_RSP returns the impulse response of Data cells
% factor is the downsampling factor to be used
h_func = @(x) ifft(ifftshift(x)); % Impulse response fct
h = zeros(xx, yy, zz, fbins); % Impulse response
rect = zeros(1, fbins);
rect(window) = 1;

for x=1:xx
    for y=1:yy
        for z=1:zz
            Hxyz = Data{1,x}{1,y}{1,z} .* rect;
            h(x,y,z,:) = h_func(Hxyz);
        end
    end
end
end

