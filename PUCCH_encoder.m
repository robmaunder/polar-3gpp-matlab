function f = PUCCH_encoder(a, E)
% PUCCH_ENCODER Physical Uplink Control Channel (PUCCH) polar encoder from 3GPP New
% Radio, as specified in Section 6.3.1 of TS 38.212 v1.1.1
%   f = PUCCH_ENCODER(a, E) encodes the information bit sequence a, in 
%   order to obtain the encoded bit sequence f.
%
%   a should be a binary row vector comprising A number of bits, each 
%   having the value 0 or 1. 
%
%   E should be an integer scalar. It specifies the number of bits in the
%   encoded bit sequence, where E should be greater than A.
%
%   f will be a binary row vector comprising E number of bits, each having
%   the value 0 or 1.
%
%   See also PUCCH_DECODER
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

% NEED TO ADD SEGMENTATION

% NEED TO CONFIRM CRC

% The CRC polynomial used in 3GPP PUCCH channel is
% D^11 + D^10 + D^9 + D^5 + 1
crc_polynomial_pattern = [1 1 1 0 0 0 1 0 0 0 0 1];
P = length(crc_polynomial_pattern)-1;


% Determine the number of information and CRC bits.
K = A+P;

if K-3 < 12
    error('polar_3gpp_matlab:UnsupportedBlockLength','K-3 should be no less than 12.');
end


% Determine the number of bits used at the input and output of the polar
% encoder kernal.
N = get_3GPP_N(K,E,10); % n_max = 10 is used in PUCCH channels

% Get a rate matching pattern.
[rate_matching_pattern, mode] = get_3GPP_rate_matching_pattern(K,N,E);

% Get a sequence pattern.
Q_N = get_3GPP_sequence_pattern(N);

if K-3 <= 22

    n_PC = 3;

    % Get an information bit pattern.
    info_bit_pattern = get_3GPP_info_bit_pattern(K+n_PC, Q_N, rate_matching_pattern, mode);

    % Get a PC bit pattern.
    if E-K+3 > 192
        PC_bit_pattern = get_PC_bit_pattern(info_bit_pattern, Q_N, n_PC, 1);
    else        
        PC_bit_pattern = get_PC_bit_pattern(info_bit_pattern, Q_N, n_PC, 0);
    end
    
    % Perform polar encoding.
    e = PCCA_polar_encoder(a, crc_polynomial_pattern, info_bit_pattern, PC_bit_pattern, 5, rate_matching_pattern);
else
    % Get an information bit pattern.
    info_bit_pattern = get_3GPP_info_bit_pattern(K, Q_N, rate_matching_pattern, mode);

    % Perform polar encoding.
    e = CA_polar_encoder(a,crc_polynomial_pattern,info_bit_pattern,rate_matching_pattern);    
end

% Perform channel interleaving.
channel_interleaver_pattern = get_3GPP_channel_interleaver_pattern(E);
f = e(channel_interleaver_pattern);


% NEED TO ADD SEGMENTATION
