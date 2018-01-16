function f = PBCH_encoder(a, E)
% PBCH_ENCODER Polar encoder for the Public Broadcast Channel (PBCH) of 3GPP New Radio, as
% defined in Section 7.1 of TS38.212. Implements the Cyclic Redudancy
% Check (CRC) attachment of Section 7.1.3, the channel coding of Section 7.1.4
% and the rate matching of Section 7.1.5. Note that this code does not
% implement the payload generation of Section 7.1.1 or the scrambling of
% Section 7.1.2.
%   f = PBCH_ENCODER(a, E) encodes the information bit sequence a, in
%   order to obtain the encoded bit sequence e.
%
%   a should be a binary row vector comprising 32 bits, each
%   having the value 0 or 1. The first input bit corresponds to a'_0 from 
%   Section 7.1.3 of TS38.212, while the last input bit corresponds 
%   to a'_A-1.
%
%   E should be 864. It specifies the number of bits in the
%   encoded bit sequence. Since there is only one valid value for this 
%   parameter, it can be omitted.
%
%   f will be a binary row vector comprising 864 bits, each having
%   the value 0 or 1. The first output bit corresponds to f_0 from Section 
%   7.1.5 of TS38.212, while the last output bit corresponds to 
%   f_E-1.
%
%   See also PBCH_DECODER
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

addpath 'components'

A = length(a);

% A is always 32 in PBCH
if A ~= 32
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be 32.');
end
% E is always 864 in PBCH
if nargin<2
    E = 864;
end
if E ~= 864
    error('polar_3gpp_matlab:UnsupportedBlockLength','E should be 864.');
end

% The CRC polynomial used in 3GPP PBCH and PDCCH channel is
% D^24 + D^23 + D^21 + D^20 + D^17 + D^15 + D^13 + D^12 + D^8 + D^4 + D^2 + D + 1
crc_polynomial_pattern = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];
P = length(crc_polynomial_pattern)-1;

% Determine the number of information and CRC bits.
K = A+P; 

% Determine the number of bits used at the input and output of the polar
% encoder kernal.
N = get_3GPP_N(K,E,9); % n_max = 9 is used in PBCH and PDCCH channels

% Get the 3GPP CRC interleaver pattern.
crc_interleaver_pattern = get_3GPP_crc_interleaver_pattern(K);

% Get the 3GPP rate matching pattern.
[rate_matching_pattern, mode] = get_3GPP_rate_matching_pattern(K,N,E);

% Get the 3GPP sequence pattern.
Q_N = get_3GPP_sequence_pattern(N);

% Get the 3GPP information bit pattern.
info_bit_pattern = get_3GPP_info_bit_pattern(K, Q_N, rate_matching_pattern, mode);

% Perform Distributed-CRC-Aided polar encoding.
f = DCA_polar_encoder(a,crc_polynomial_pattern,crc_interleaver_pattern,info_bit_pattern,rate_matching_pattern);
