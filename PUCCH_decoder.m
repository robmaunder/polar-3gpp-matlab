function a_hat = PUCCH_decoder(f_tilde, A, L, min_sum)
% PUCCH_DECODER Polar decoder for the Physical Uplink Control Channel (PUCCH) and the
% Physical Uplink Shared Channel (PUSCH) of 3GPP New Radio, as defined in
% Section 6.3 of TS38.212. Implements the code block segmentation and
% Cyclic Redudancy Check (CRC) attachment of Sections 6.3.1.2.1 and 6.3.2.2.1,
% the channel coding of Sections 6.3.1.3.1 and 6.3.2.3.2, the rate matching of
% Sections 6.3.1.4.1 and 6.3.2.4.1, as well as the code block concatenation of
% Sections 6.3.1.5.1 and 6.3.2.5.1. Note that this code does not implement the
% UCI bit sequence generation of Sections 6.3.1.1 and 6.3.2.1, the
% determination of the encoded block length E_UCI of Sections 6.3.1.4.1 and
% 6.3.2.4.1, or the multiplexing of Sections 6.3.1.6 and 6.3.2.6. Also, this
% code does not implement the small block lengths, which are detailed in
% Sections 6.3.1.2.2, 6.3.1.3.2, 6.3.1.4.2, 6.3.2.2.2, 6.3.2.3.2 and 6.3.2.4.2.
%   a_hat = PUCCH_DECODER(f_tilde, A, L, min_sum) decodes the encoded LLR sequence
%   f_tilde, in order to obtain the recovered information bit sequence
%   a_hat.
%
%   f_tilde should be a real row vector comprising G number of Logarithmic
%   Likelihood Ratios (LLRS), each having a value obtained as LLR =
%   ln(P(bit=0)/P(bit=1)), where G should be no greater than 8192 if A<360
%   and no greater than 16384 if A>=360. The first LLR corresponds to g_0 
%   from Sections 6.3.1.5 and 6.3.2.5 of TS38.212, while the last 
%   LLR corresponds to g_G-1.
%
%   A should be an integer scalar. It specifies the number of bits in the
%   information bit sequence, where A should be in the range 12 to 1706.
%
%   L should be a scalar integer. It specifies the list size to use during
%   Successive Cancellation List (SCL) decoding.
%
%   min_sum shoular be a scalar logical. If it is true, then the SCL
%   decoding process will be completed using the min-sum approximation.
%   Otherwise, the log-sum-product will be used. The log-sum-product gives
%   better error correction capability than the min-sum, but it has higher
%   complexity.
%
%   a_hat will be a binary row vector comprising A number of bits, each
%   having the value 0 or 1. The first output bit corresponds to a_0 from 
%   Sections 6.3.1.2 and 6.3.2.2 of TS38.212, while the last output 
%   bit corresponds to a_A-1.
%
%   See also PUCCH_ENCODER
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

G = length(f_tilde);

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
    
    if (A >= 360 && G >= 1088) || A >= 1013
        % Use two segments
        C = 2;
    else
        % Use one segment
        C = 1;
    end  
end

% The CRC has P bits. P-min(P2,log2(L)) of these are used for error
% detection, where L is the list size. Meanwhile, min(P2,log2(L)) of
% them are used to improve error correction. So the CRC needs to be
% min(P2,log2(L)) number of bits longer than CRCs used in other codes,
% in order to achieve the same error detection capability.
P = length(crc_polynomial_pattern)-1;
P2 = 3;

% Determine the number of information and CRC bits.
K = ceil(A/C)+P;

E_r = floor(G/C);

if E_r > 8192
    error('polar_3gpp_matlab:UnsupportedBlockLength','G is too long.');
end

% Determine the number of bits used at the input and output of the polar
% encoder kernal.
N = get_3GPP_N(K,E_r,10); % n_max = 10 is used in PUCCH channels

% Get a rate matching pattern.
[rate_matching_pattern, mode] = get_3GPP_rate_matching_pattern(K,N,E_r);

% Get a sequence pattern.
Q_N = get_3GPP_sequence_pattern(N);

% Get the channel interleaving pattern
channel_interleaver_pattern = get_3GPP_channel_interleaver_pattern(E_r);



if A <= 19 % Use PCCA-polar
    
    % Perform channel interleaving
    e_tilde = zeros(1,G);
    e_tilde(channel_interleaver_pattern) = f_tilde;
    
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
    
    % Perform polar decoding.
    a_hat = PCCA_polar_decoder(e_tilde,crc_polynomial_pattern,info_bit_pattern,PC_bit_pattern,5,rate_matching_pattern,mode,L,min_sum,P2);
else % Use CA-polar
        
    % Get an information bit pattern.
    info_bit_pattern = get_3GPP_info_bit_pattern(K, Q_N, rate_matching_pattern, mode);
        
    if C == 2
        
        info_bit_pattern2 = info_bit_pattern;
        if mod(A,2) ~= 0 % If a zero was prepended to the first segment during encoding
            % Treat the first information bit as a frozen bit
            info_bit_pattern2(find(info_bit_pattern == 1,1,'first')) = 0;
        end
        
        % Perform channel interleaving for first segment.
        e_tilde = zeros(1,E_r);        
        e_tilde(channel_interleaver_pattern) = f_tilde(1:E_r);
                
        % Perform polar decoding for first segment.
        a_hat = CA_polar_decoder(e_tilde,crc_polynomial_pattern,info_bit_pattern2,rate_matching_pattern,mode,L,min_sum,P2);
        
        if length(a_hat) == floor(A/C)
        
            % Perform channel interleaving for second segment.
            e_tilde(channel_interleaver_pattern) = f_tilde(E_r+1:2*E_r);

            % Perform polar decoding for second segment.
            a_hat = [a_hat, CA_polar_decoder(e_tilde,crc_polynomial_pattern,info_bit_pattern,rate_matching_pattern,mode,L,min_sum,P2)];
        end
        
        if length(a_hat) ~= A
            a_hat = [];
        end
    else
        % Perform channel interleaving
        e_tilde = zeros(1,G);        
        e_tilde(channel_interleaver_pattern) = f_tilde;
        
        % Perform polar decoding.
        a_hat = CA_polar_decoder(e_tilde,crc_polynomial_pattern,info_bit_pattern,rate_matching_pattern,mode,L,min_sum,P2);
        
    end
end
end
