function e = polar_encoder(a, info_bit_pattern, rate_matching_pattern)
% POLAR_ENCODER Basic polar encoder.
%   e = POLAR_ENCODER(a, info_bit_pattern, rate_matching_pattern) encodes
%   the information bit sequence a, in order to obtain the encoded bit 
%   sequence e.
%
%   a should be a binary row vector comprising A number of bits, each 
%   having the value 0 or 1. 
%
%   info_bit_pattern should be a row vector comprising N number of logical 
%   elements, each having the value true or false. The number of elements 
%   in info_bit_pattern having the value true should be A. These elements 
%   having the value true identify the positions of the information bits 
%   within the input to the polar encoder kernal.
%
%   rate_matching_pattern should be a row vector comprising E number of
%   integers, each having a value in the range 1 to N. Each integer
%   identifies which one of the N outputs from the polar encoder kernal
%   provides the corresponding bit in the encoded bit sequence e.
%
%   e will be a binary row vector comprising E number of bits, each having
%   the value 0 or 1.
%
%   See also POLAR_DECODER
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.


A = length(a);
N = length(info_bit_pattern);

if log2(N) ~= round(log2(N))
    error('N should be a power of 2');  
end
if sum(info_bit_pattern) ~= A
    error('info_bit_pattern should contain A number of ones.');
end
if max(rate_matching_pattern) > N
    error('rate_matching_pattern is not compatible with N');
end

% Position the information bits within the input to the polar encoder 
% kernal.
u = zeros(1,N);
u(info_bit_pattern) = a;

% Perform the polar encoder kernal operation.
G_N = get_G_N(N);
d = mod(u*G_N,2);

% Extract the encoded bits from the output of the polar encoder kernal.
e = d(rate_matching_pattern);

end
