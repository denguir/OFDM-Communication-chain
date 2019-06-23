function [Rx,Ry,Rz] = spatial_corr(main_a, beta, dX, dY, dZ)
%SPATIAL_CORR Summary of this function goes here
%   Detailed explanation goes here

Rx = zeros(length(dX), 1);
Ry = zeros(length(dY), 1);
Rz = zeros(length(dZ), 1);
for i=1:length(dX)
    phaseX = exp(1j * beta(1,:) * dX(i));
    phaseY = exp(1j * beta(2,:) * dY(i));
    phaseZ = exp(1j * beta(3,:) * dZ(i));
    Rx(i) = dot(abs(main_a).^2, phaseX);
    Ry(i) = dot(abs(main_a).^2, phaseY);
    Rz(i) = dot(abs(main_a).^2, phaseZ);
end

