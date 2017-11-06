function G_P = get_crc_generator_matrix(A, crc_polynomial_pattern)
% GET_CRC_GENERATOR_MATRIX Obtain a Cyclic Redudancy Check (CRC) generator 
% matrix.
%   G_P = GET_CRC_GENERATOR_MATRIX(A, crc_polynomial_pattern) obtains the CRC
%   generator matrix G.
%
%   A should be an integer scalar. It specifies the number of bits in the
%   information bit sequence.
%
%   crc_polynomial_pattern should be a binary row vector comprising P+1
%   number of bits, each having the value 0 or 1. These bits parameterise a
%   Cyclic Redundancy Check (CRC) comprising P bits. Each bit provides the
%   coefficient of the corresponding element in the CRC generator
%   polynomial. From left to right, the bits provide the coefficients for
%   the elements D^P, D^P-1, D^P-2, ..., D^2, D, 1.
%
%   G_P will be a K by P binary matrix. The CRC bits can be generated
%   according to mod(a*G_P,2).
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.


P = length(crc_polynomial_pattern)-1;

if P<1
    error('crc_polynomial_pattern is invalid');
end

G_P = zeros(A,P);

if A>0
    G_P(end,:) = crc_polynomial_pattern(2:end);

    for k = A-1:-1:1
        G_P(k,:) = xor([G_P(k+1,2:end),0],G_P(k+1,1)*crc_polynomial_pattern(2:end));
    end
end