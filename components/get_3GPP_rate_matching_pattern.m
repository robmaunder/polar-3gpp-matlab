function [rate_matching_pattern, mode] = get_3GPP_rate_matching_pattern(K, N, E)
% GET_3GPP_RATE_MATCHING_PATTERN Get the 3GPP rate matching 
% sequence, as specified in Sections 5.4.1.1 and 5.4.1.2 of TS 38.212 
%   [rate_matching_pattern, mode] = GET_3GPP_RATE_MATCHING_PATTERN(K,N,E)
%   obtains the rate matching sequence.
%
%   K should be an integer scalar. It specifies the number of bits in the 
%   information and CRC bit sequence.
%
%   N should be an integer scalar, which should be a power of 2 and no less 
%   than 32. It specifies the number of bits in the input and output of the 
%   polar encoder kernal. 
%
%   E should be an integer scalar. It specifies the number of bits in the 
%   encoded bit sequence. 
%
%   rate_matching_pattern will be a row vector comprising E number of
%   integers, each having a value in the range 1 to N. Each integer
%   identifies which one of the N outputs from the polar encoder kernal
%   provides the corresponding bit in the encoded bit sequence. Rate
%   matching can be achieved according to e = d(rate_matching_pattern).
%
%   mode will have the value 'repetition', 'puncturing' or 'shortening'.
%   This specifies how the rate matching has been achieved. 'repetition'
%   indicates that some outputs of the polar encoder kernal are repeated in 
%   the encoded bit sequence two or more times. 'puncturing' and 
%   'shortening' indicate that some outputs of the polar encoder kernal 
%   have been excluded from the encoded bit sequence. In the case of
%   'puncturing' these excluded bits could have values of 0 or 1. In the
%   case of 'shortening' these excluded bits are guaranteed to have values
%   of 0.
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
if n < 5
    error('polar_3gpp_matlab:UnsupportedBlockLength','N should be no smaller than 32');
end

P = [0 1 2 4 3 5 6 7 8 16 9 17 10 18 11 19 12 20 13 21 14 22 15 23 24 25 26 28 27 29 30 31];

d = 1:N;

J = zeros(1,N);
y = zeros(1,N);
for n=0:N-1
    i=floor(32*n/N);
    J(n+1)=P(i+1)*(N/32)+mod(n,N/32);
    y(n+1) = d(J(n+1)+1);
end

rate_matching_pattern = zeros(1,E);
if E >= N
    for k = 0:E-1
        rate_matching_pattern(k+1) = y(mod(k,N)+1);
    end
    mode = 'repetition';
else
    if K/E <= 7/16
        for k=0:E-1
            rate_matching_pattern(k+1) = y(k+N-E+1);
        end
        mode = 'puncturing';
    else
        for k = 0:E-1
            rate_matching_pattern(k+1) = y(k+1);
        end
        mode = 'shortening';
    end
end
    


