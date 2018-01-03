function N = get_3GPP_N(K,E,n_max)
% GET_3GPP_N Obtain the number of bits in the input and output of the polar
% encoder kernal, according to Section 5.3.1 of 3GPP TS 38.212
%   N = GET_3GPP_N(K,E,n_max) obtains the number of bits in the input and
%   output of the polar encoder kernal N.
%
%   K should be an integer scalar. It specifies the number of bits in the
%   information and CRC bit sequence.
%
%   E should be an integer scalar. It specifies the number of bits in the
%   encoded bit sequence.
%
%   n_max should be an integer scalar. It specifies the log2 of the maximum
%   number of bits in the input and output of the polar encoder kernal.
%   n_max = 9 in the PBCH and PDCCH channels, while n_max = 10 in the PUCCH
%   channel.
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.


if E <= (9/8)*2^(ceil(log2(E))-1) && K/E < 9/16
    n_1=ceil(log2(E))-1;
else
    n_1=ceil(log2(E));
end

R_min=1/8;
n_min=5;
n_2=ceil(log2(K/R_min));
n=max(n_min,min([n_1,n_2,n_max]));

N=2^n;
end