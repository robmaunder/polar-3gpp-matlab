function info_bit_pattern = get_info_bit_pattern(I, Q_N, rate_matching_pattern)
% GET_INFO_BIT_PATTERN Obtain a generic information bit pattern.
%   info_bit_pattern = GET_INFO_BIT_PATTERN(I, Q_N, rate_matching_pattern)
%   obtains the information bit pattern info_bit_pattern.
%
%   I should be an integer scalar. It specifies the number of bits in the 
%   information and CRC bit sequence. It should be no greater than N or E.
%
%   Q_N should be a row vector comprising N number of unique integers in the 
%   range 1 to N. Each successive element of Q_N provides the index of the
%   next most reliable input to the polar encoder kernal, where the first
%   element of Q_N gives the index of the least reliable bit and the last
%   element gives the index of the most reliable bit.
%
%   rate_matching_pattern should be a row vector comprising E number of
%   integers, each having a value in the range 1 to N. Each integer
%   identifies which one of the N outputs from the polar encoder kernal
%   provides the corresponding bit in the encoded bit sequence e.
%
%   info_bit_pattern will be a row vector comprising N number of logical 
%   elements, each having the value true or false. The number of elements 
%   in info_bit_pattern having the value true will be I. These elements 
%   having the value true identify the positions of the information and 
%   CRC bits (if any) within the input to the polar encoder kernal. The
%   information bit arrangement can be achieved according to
%   u(info_bit_pattern) = a.
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.


N = length(Q_N);
n = log2(N);
E = length(rate_matching_pattern);

if n ~= round(n)
    error('N should be a power of 2');
end
if max(rate_matching_pattern) > N
    error('rate_matching_pattern is not compatible with N');
end
if I > N
    error('polar_3gpp_matlab:UnsupportedBlockLength','I should be no greater than N');
end
if I > E
    error('polar_3gpp_matlab:UnsupportedBlockLength','I should be no greater than E');
end

Q_Itmp_N = intersect(Q_N, rate_matching_pattern, 'stable');

Q_I_N=Q_Itmp_N(end-I+1:end);

info_bit_pattern = false(1,N);
info_bit_pattern(Q_I_N) = true;