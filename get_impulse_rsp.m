function [h] = get_impulse_rsp(Data, xx, yy, zz, fbins, window)
%GET_IMPULSE_RSP returns the impulse response of Data cells
% xx, yy, zz are the number of measures along x,y,z axis
% fbins is the number of frequcency bins used for the frequency response
% window is the filter to apply on the frequency response

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

