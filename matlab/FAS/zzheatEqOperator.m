function [ A ] = zzheatEqOperator( H, S, params )
%NONLINEARDIFFUSIONOPERATOR Summary of this function goes here
%   Detailed explanation goes here

alpha = params.dt/params.dx^2;
Valpha = params.V*params.dt/params.dx;


A = diag(H + Valpha*H + 2*alpha*temperature(H, S, params), 0) ...
    - diag(alpha*temperature(H(1:end-1),S(1:end-1), params)   + Valpha*H(1:end-1) ,-1) ...
    - diag(alpha*temperature(H(1:end-1),S(1:end-1), params),1);


% Linear diffusion
onesArr = ones(length(H), 1);
%A = diag((1 + 2*alpha)*onesArr, 0) ...
%    - diag(alpha* onesArr(1:end-1),-1) ...
%    - diag(alpha* onesArr(1:end-1),1);


%Now sort out BCs
A(1,1) = 1; A(1,2) = 0;
A(end,end) = 1; A(end,end-1) =0;


end

