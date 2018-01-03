function G_N = get_G_N(N)
% GET_G_N The nth Kronecker power of G_2, as described in Section
% 5.3.1.2 of 3GPP TS 38.212
%   G_N = GET_G_N(N) obtains the nth Kronecker power of G_2.
%
%   N should be an integer scalar, which should be a power of 2. It 
%   specifies the number of bits in the input and output of the polar 
%   encoder kernal. 
%
%   G_N will be an N by N binary matrix. The polar encoding kernal process
%   may be performed according to d = mod(u*G_N,2).
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

n = log2(N);
if n ~= round(n)
    error('N should be a power of 2');
end

G_N = 1;
for i=1:n
    G_N = kron(G_N,[1 0; 1 1]);
end