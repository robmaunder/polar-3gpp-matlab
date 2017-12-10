function f = PUCCH_encoder(a, G)
% PUCCH_ENCODER Physical Uplink Control Channel (PUCCH) polar encoder from 3GPP New
% Radio, as specified in Section 6.3.1 of TS 38.212 v1.1.1
%   f = PUCCH_ENCODER(a, G) encodes the information bit sequence a, in
%   order to obtain the encoded bit sequence f.
%
%   a should be a binary row vector comprising A number of bits, each
%   having the value 0 or 1.
%
%   E should be an integer scalar. It specifies the number of bits in the
%   encoded bit sequence, where G should be greater than A.
%
%   f will be a binary row vector comprising G number of bits, each having
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

if A < 12
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be no less than 12.');
elseif A > 1706
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be no greater than 1706.');    
elseif A <= 19 % Use PCCA-polar
    % The CRC polynomial used with PCCA-polar in 3GPP PUCCH channel is
    % D^6 + D^5 + 1
    crc_polynomial_pattern = [1 1 0 0 0 0 1];
    
    % Use one segment
    C = 1;
else % Use CA-polar
    % The CRC polynomial used with CA-polar in 3GPP PUCCH channel is
    % D^11 + D^10 + D^9 + D^5 + 1
    crc_polynomial_pattern = [1 1 1 0 0 0 1 0 0 0 0 1];
    
    if A >= 360 && G >= 1088
        if mod(G,2) ~= 0
            error('polar_3gpp_matlab:UnsupportedBlockLength','G should be divisible by 2 when code block segmentation is used.');
        end
        
        % Use two segments
        C = 2;
        
    else
        % Use one segment
        C = 1;
    end
    
end



if G > 8192*C
    error('polar_3gpp_matlab:UnsupportedBlockLength','G is too long.');
end
if isempty(find(unique([unique(16*(1:2)'*(1:12));unique(24*(4:14)'*[1,2,3,4,5,6,8,9,10,12,15,16]);unique(24*(4:14)'*1./[2,4])])==G,1))
    error('polar_3gpp_matlab:UnsupportedBlockLength','G should be selected from the set supported in PDDCH.');
end


% Determine the number of information and CRC bits.
P = length(crc_polynomial_pattern)-1;
K = ceil(A/C)+P;

% Determine the number of bits used at the input and output of the polar
% encoder kernal.
N = get_3GPP_N(K,G/C,10); % n_max = 10 is used in PUCCH channels

% Get a rate matching pattern.
[rate_matching_pattern, mode] = get_3GPP_rate_matching_pattern(K,N,G/C);

% Get a sequence pattern.
Q_N = get_3GPP_sequence_pattern(N);

% Get the channel interleaving pattern
channel_interleaver_pattern = get_3GPP_channel_interleaver_pattern(G/C);

if A <= 19 % Use PCCA-polar
    % We use 3 PC bits
    n_PC = 3;
    
    % Get an information bit pattern.
    info_bit_pattern = get_3GPP_info_bit_pattern(K+n_PC, Q_N, rate_matching_pattern, mode);
    
    % Get a PC bit pattern.
    if G-K+3 > 192
        PC_bit_pattern = get_PC_bit_pattern(info_bit_pattern, Q_N, n_PC, 1);
    else
        PC_bit_pattern = get_PC_bit_pattern(info_bit_pattern, Q_N, n_PC, 0);
    end
    
    % Perform polar encoding.
    e = PCCA_polar_encoder(a, crc_polynomial_pattern, info_bit_pattern, PC_bit_pattern, 5, rate_matching_pattern);
    
    % Perform channel interleaving.
    f = e(channel_interleaver_pattern);
    
else % Use CA-polar
    
    % Get an information bit pattern.
    info_bit_pattern = get_3GPP_info_bit_pattern(K, Q_N, rate_matching_pattern, mode);
    
    if C == 2

        info_bit_pattern2 = info_bit_pattern;
        if mod(A,2) ~= 0 % Prepend a zero to the first segment during encoding
            % Treat the first information bit as a frozen bit
            info_bit_pattern2(find(info_bit_pattern == 1,1,'first')) = 0;
        end
        
        % Perform polar encoding for first segment.
        e = CA_polar_encoder(a(1:floor(A/C)),crc_polynomial_pattern,info_bit_pattern2,rate_matching_pattern);
        
        % Perform channel interleaving for first segment.
        f = e(channel_interleaver_pattern);
        
        % Perform polar encoding for second segment.
        e = CA_polar_encoder(a(floor(A/C)+1:A),crc_polynomial_pattern,info_bit_pattern,rate_matching_pattern);
        
        % Perform channel interleaving for second segment.
        f = [f,e(channel_interleaver_pattern)];
        
    else
        % Perform polar encoding.
        e = CA_polar_encoder(a,crc_polynomial_pattern,info_bit_pattern,rate_matching_pattern);
        
        % Perform channel interleaving.
        f = e(channel_interleaver_pattern);
    end
end

