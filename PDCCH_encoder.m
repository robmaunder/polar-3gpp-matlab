function f = PDCCH_encoder(a, E, RNTI)
% PDCCH_ENCODER Physical Downlink Control Channel (PDCCH) polar encoder from 3GPP New
% Radio, as specified in Section 7.3 of TS 38.212 v1.0.1...
% http://www.3gpp.org/ftp/TSG_RAN/WG1_RL1/TSGR1_AH/NR_AH_1709/Docs/R1-1716928.zip
%   f = PDCCH_ENCODER(a, E) encodes the information bit sequence a, in
%   order to obtain the encoded bit sequence e.
%
%   a should be a binary row vector comprising A number of bits, each
%   having the value 0 or 1. A should be in the range 1 to 200.
%
%   E should be an integer scalar. It specifies the number of bits in the
%   encoded bit sequence, where E should greater than A.
%
%   RNTI should be a binary row vector comprising 16 bits, each having the
%   value 0 or 1. If this parameter is omitted, then ones(1,16) will be
%   used for the RNTI.
%
%   f will be a binary row vector comprising E number of bits, each having
%   the value 0 or 1.
%
%   See also PDCCH_DECODER
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

if nargin == 2
    RNTI = ones(1,16);
end
if length(RNTI) ~= 16
    error('RNTI length should be 16');
end


A = length(a);

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
f = DSCA_polar_encoder(a,crc_polynomial_pattern, RNTI, crc_interleaver_pattern,info_bit_pattern,rate_matching_pattern);
