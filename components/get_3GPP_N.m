function N = get_3GPP_N(K,E,n_max)
% GET_3GPP_N Obtain the number of bits in the input and output of the polar
% encoder kernal, according to Section 5.3.1 of 3GPP TS 38.212 V1.0.1...
% http://www.3gpp.org/ftp/TSG_RAN/WG1_RL1/TSGR1_AH/NR_AH_1709/Docs/R1-1716928.zip
%   N = GET_3GPP_N(K,E,n_max) obtains the number of bits in the input and
%   output of the polar encoder kernal N.
%
%   K should be an integer scalar. It specifies the number of bits in the
%   information and CRC bit sequence. It should be less than E and no
%   greater than 2^n_max.
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

if K >= E
    error('K should be less than E.');
end
if K > 2^n_max
    error('K should be no greater than 2^n_max.');
end


if E <= (9/8)*2^(ceil(log2(E))-1) && K/E < 9/16
    n_1=ceil(log2(E))-1;
else
    n_1=ceil(log2(E));
end

R_min=1/8;
n_2=ceil(log2(K/R_min));
n=min([n_1,n_2,n_max]);

n = max(n,5); % This is not stated in 3GPP TS 38.212 V1.0.1, but it is necessary for compatibility with the sub-block interleaver during rate-matching

N=2^n;
end