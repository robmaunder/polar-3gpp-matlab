function a_hat = PDCCH_decoder(f_tilde, A, L, min_sum, RNTI)
% PDCCH_DECODER Polar decoder for the Physical Downlink Control Channel (PDCCH) of 3GPP New
% Radio, as defined in Section 7.3 of TS38.212. Implements the zero-
% padding to increase the length of short payloads to 12 bits of Section 7.3.1,
% the Cyclic Redudancy Check (CRC) attachment of Section 7.3.2, the channel
% coding of Section 7.3.3 and the rate matching of Section 7.3.4. Note that
% this code does not implement the DCI bit sequence generation of Section
% 7.3.1.
%   a_hat = PDCCH_DECODER(f_tilde, A, L, min_sum) decodes the encoded LLR sequence 
%   f_tilde, in order to obtain the recovered information bit sequence 
%   a_hat.
%
%   f_tilde should be a real row vector comprising E number of Logarithmic
%   Likelihood Ratios (LLRS), each having a value obtained as LLR =
%   ln(P(bit=0)/P(bit=1)). The first LLR corresponds to f_0 from Section 
%   7.3.4 of TS38.212, while the last LLR corresponds to 
%   f_E-1.
%
%   A should be an integer scalar. It specifies the number of bits in the
%   information bit sequence, where A should be in the range 1 to 140.
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
%   RNTI should be a binary row vector comprising 16 bits, each having the
%   value 0 or 1. If this parameter is omitted, then ones(1,16) will be
%   used for the RNTI. The first bit corresponds to x_rnti,0 from Section
%   7.3.2 of TS38.212, while the last bit corresponds to x_rnti,15.
%
%   a_hat will be a binary row vector comprising A number of bits, each 
%   having the value 0 or 1. The first output bit corresponds to a_0 from 
%   Section 7.3.1 of TS38.212, while the last output bit corresponds 
%   to a_A-1.
%
%   See also PDCCH_ENCODER
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

if A == 0
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be no less than 1.');
end
if A > 140
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be no greater than 140.');
end

if nargin == 4
    RNTI = ones(1,16);
end
if length(RNTI) ~= 16
    error('RNTI length should be 16');
end

E = length(f_tilde);
if E > 8192
    error('polar_3gpp_matlab:UnsupportedBlockLength','E should be no greater than 8192.');
end

% The CRC polynomial used in 3GPP PBCH and PDCCH channel is
% D^24 + D^23 + D^21 + D^20 + D^17 + D^15 + D^13 + D^12 + D^8 + D^4 + D^2 + D + 1
crc_polynomial_pattern = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];

% The CRC has P bits. P-min(P2,log2(L)) of these are used for error
% detection, where L is the list size. Meanwhile, min(P2,log2(L)) of
% them are used to improve error correction. So the CRC needs to be
% min(P2,log2(L)) number of bits longer than CRCs used in other codes,
% in order to achieve the same error detection capability.
P = length(crc_polynomial_pattern)-1;
P2 = 3;

% Determine the number of information and CRC bits.
if A < 12
    % a has been padded with zeros to increase its length to 12
    K = 12+P;
else
    K = A+P;
end

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

if A < 12
    % We know that a has been padded with zeros
    a_tilde = [NaN(1,A),zeros(1,12-A)];
    
    % Perform Distributed-CRC-Aided polar decoding.
    a_hat = DS1CKA_polar_decoder(f_tilde,crc_polynomial_pattern,RNTI,crc_interleaver_pattern,info_bit_pattern,rate_matching_pattern,mode,L,min_sum,P2,a_tilde);
    
    if ~isempty(a_hat)
        % Remove the padding
        a_hat = a_hat(1:A);
    end    
else
    % Perform Distributed-CRC-Aided polar decoding.
    a_hat = DS1CA_polar_decoder(f_tilde,crc_polynomial_pattern,RNTI,crc_interleaver_pattern,info_bit_pattern,rate_matching_pattern,mode,L,min_sum,P2);
end