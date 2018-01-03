function e = PCCA_polar_encoder(a, crc_polynomial_pattern, info_bit_pattern, PC_bit_pattern, PC_circular_buffer_length, rate_matching_pattern)
% PCCA_POLAR_ENCODER Parity-Check and CRC-Aided (PC-CA) polar encoder.
%   e = PCCA_POLAR_ENCODER(a, crc_polynomial_pattern, info_bit_pattern, PC_bit_pattern, rate_matching_pattern) 
%   encodes the information bit sequence a, in order to obtain the encoded 
%   bit sequence e.
%
%   a should be a binary row vector comprising A number of bits, each 
%   having the value 0 or 1. 
%
%   crc_polynomial_pattern should be a binary row vector comprising P+1
%   number of bits, each having the value 0 or 1. These bits parameterise a
%   Cyclic Redundancy Check (CRC) comprising P bits. Each bit provides the
%   coefficient of the corresponding element in the CRC generator
%   polynomial. From left to right, the bits provide the coefficients for
%   the elements D^P, D^P-1, D^P-2, ..., D^2, D, 1.
%
%   info_bit_pattern should be a row vector comprising N number of logical 
%   elements, each having the value true or false. The number of elements 
%   in info_bit_pattern having the value true should be K+n_PC, where K = A+P. 
%   These elements having the value true identify the positions of the 
%   information, CRC and PC bits within the input to the polar encoder kernal.
%
%   PC_bit_pattern should be a row vector comprising N number of logical 
%   elements, each having the value true or false. The number of elements 
%   in PC_bit_pattern having the value true should be n_PC. 
%   These elements having the value true identify the positions of the 
%   PC bits within the input to the polar encoder kernal.
%
%   PC_circular_buffer_length should be an integer scalar. It specifies the
%   length of the circular buffer used to generate the PC bits.
%
%   rate_matching_pattern should be a row vector comprising E number of
%   integers, each having a value in the range 1 to N. Each integer
%   identifies which one of the N outputs from the polar encoder kernal
%   provides the corresponding bit in the encoded bit sequence e.
%
%   e will be a binary row vector comprising E number of bits, each having
%   the value 0 or 1.
%
%   See also PCCA_POLAR_DECODER
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
P = length(crc_polynomial_pattern)-1;
n_PC = sum(PC_bit_pattern);
I = A+P+n_PC;
N = length(info_bit_pattern);

if log2(N) ~= round(log2(N))
    error('N should be a power of 2');
end
if sum(info_bit_pattern) ~= I
    error('info_bit_pattern should contain I number of ones.');
end
if max(rate_matching_pattern) > N
    error('rate_matching_pattern is not compatible with N');
end

% Generate the CRC bits and append them to the information bits.
G_P = get_crc_generator_matrix(A,crc_polynomial_pattern);
b = [a, mod(a*G_P,2)];

% Generate the PC bits and position them together with the information and CRC 
% bits within the input to the polar encoder kernal, according to Section 
% 5.3.1.2 of TS 38.212.
u = zeros(1,N);
k = 0;
y = zeros(1,PC_circular_buffer_length);
for n=0:N-1
    y = [y(2:end),y(1)];
    if info_bit_pattern(n+1)
        if PC_bit_pattern(n+1)
            u(n+1) = y(1);
        else
            u(n+1) = b(k+1);
            k=k+1;
            y(1) = xor(y(1),u(n+1));
        end
    end
end

% Perform the polar encoder kernal operation.
G_N = get_G_N(N);
d = mod(u*G_N,2);

% Extract the encoded bits from the output of the polar encoder kernal.
e = d(rate_matching_pattern);

end
