function a_hat = PBCH_decoder(f_tilde, A, L, min_sum, a_tilde)
% PCBH_DECODER Polar decoder for the Public Broadcast Channel (PBCH) of 3GPP New Radio, as
% defined in Section 7.1 of TS38.212. Implements the Cyclic Redudancy
% Check (CRC) attachment of Section 7.1.3, the channel coding of Section 7.1.4
% and the rate matching of Section 7.1.5. Note that this code does not
% implement the payload generation of Section 7.1.1 or the scrambling of
% Section 7.1.2.
%   a_hat = PBCH_DECODER(f_tilde, A, L, min_sum, a_tilde) decodes the encoded LLR sequence 
%   f_tilde, in order to obtain the recovered information bit sequence 
%   a_hat.
%
%   f_tilde should be a real row vector comprising 864 Logarithmic
%   Likelihood Ratios (LLRS), each having a value obtained as LLR =
%   ln(P(bit=0)/P(bit=1)). The first LLR corresponds to f_0 from Section 
%   7.1.5 of TS38.212, while the last LLR corresponds to 
%   f_E-1.
%
%   A should be 32. It specifies the number of bits in the
%   information bit sequence.
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
%   a_tilde should be a row vector comprising 32 elements. Elements having
%   the value 0 or 1 indicate that the receiver has prior knowledge that
%   the corresponding information bit has the corresponding value. Elements
%   having any other value indicate that the receiver does not have prior
%   knowledge of the value of the corresponding information bit. This
%   parameter can be omitted, in which case it will be assumed that the
%   receiver has no prior knowledge about the value of any of the
%   information bits. Note that this code does not perform the 
%   scrambling of Section 7.1.1 in TS38.212 V1.1.2, since this depends on 
%   several higher-layer parameters. If scrambling is implemented
%   externally to this code, then a_tilde should pertain to the scrambled
%   bit values. The first input bit corresponds to a'_0 from 
%   Section 7.1.3 of TS38.212, while the last input bit corresponds 
%   to a'_A-1.
%
%   a_hat will be a binary row vector comprising 32 bits, each 
%   having the value 0 or 1. The first output bit corresponds to a'_0 from 
%   Section 7.1.3 of TS38.212, while the last output bit corresponds 
%   to a'_A-1.
%
%
%   See also PBCH_ENCODER
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

E = length(f_tilde);

% A is always 32 in PBCH
if A ~= 32
    error('polar_3gpp_matlab:UnsupportedBlockLength','A should be 32.');
end
% E is always 864 in PBCH
if E ~= 864
    error('polar_3gpp_matlab:UnsupportedBlockLength','E should be 864.');
end

if nargin < 5
    a_tilde = nan(1,A);
end
if length(a_tilde) ~= A
    error('a_tilde should contain A number of elements.');
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

% Perform Distributed-CRC-and-Known-bit-Aided polar decoding.
a_hat = DCKA_polar_decoder(f_tilde,crc_polynomial_pattern,crc_interleaver_pattern,info_bit_pattern,rate_matching_pattern,mode,L,min_sum,P2,a_tilde);


