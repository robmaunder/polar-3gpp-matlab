function a_hat = DS1CKA_polar_decoder(e_tilde, crc_polynomial_pattern, crc_scrambling_pattern, crc_interleaver_pattern, info_bit_pattern, rate_matching_pattern, mode, L, min_sum, P2, a_tilde)
% DS1CKA_POLAR_DECODER Distributed-Scrambled-and-1-initialised-CRC-and-Known-bit-Aided (DS1CKA) polar decoder.
%   a_hat = DS1CKA_POLAR_DECODER(e_tilde, crc_polynomial_pattern, crc_interleaver_pattern, info_bit_pattern, rate_matching_pattern, mode, L, min_sum, P2, a_tilde) 
%   decodes the encoded LLR sequence e_tilde, in order to obtain the
%   recovered information bit sequence a_hat.
%
%   e_tilde should be a real row vector comprising E number of Logarithmic
%   Likelihood Ratios (LLRS), each having a value obtained as LLR =
%   ln(P(bit=0)/P(bit=1)).
%
%   crc_polynomial_pattern should be a binary row vector comprising P+1
%   number of bits, each having the value 0 or 1. These bits parameterise a
%   Cyclic Redundancy Check (CRC) comprising P bits. Each bit provides the
%   coefficient of the corresponding element in the CRC generator
%   polynomial. From left to right, the bits provide the coefficients for
%   the elements D^P, D^P-1, D^P-2, ..., D^2, D, 1.
%
%   crc_scrambling_pattern should be a binary row vector, with each element
%   having the value 0 or 1. This vector is right-aligned with the
%   vector of CRC bits (before CRC interleaving in the encoder), then 
%   applied using XOR operations.
%
%   crc_interleaver_pattern should be a row vector comprising K number of
%   integers, each having a unique value in the range 1 to K. Each integer
%   identifies which one of the K information or CRC bits provides the 
%   corresponding bit in the input to the polar encoder kernal.
%
%   info_bit_pattern should be a row vector comprising N number of logical 
%   elements, each having the value true or false. The number of elements 
%   in info_bit_pattern having the value true should be K, where K = A+P. 
%   These elements having the value true identify the positions of the 
%   information and CRC bits within the input to the polar encoder kernal.
%
%   rate_matching_pattern should be a row vector comprising E number of
%   integers, each having a value in the range 1 to N. Each integer
%   identifies which one of the N outputs from the polar encoder kernal
%   provides the corresponding bit in the encoded bit sequence e.
%
%   mode should have the value 'repetition', 'puncturing' or 'shortening'.
%   This specifies how the rate matching has been achieved. 'repetition'
%   indicates that some outputs of the polar encoder kernal are repeated in 
%   the encoded bit sequence two or more times. 'puncturing' and 
%   'shortening' indicate that some outputs of the polar encoder kernal 
%   have been excluded from the encoded bit sequence. In the case of
%   'puncturing' these excluded bits could have values of 0 or 1. In the
%   case of 'shortening' these excluded bits are guaranteed to have values
%   of 0.
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
%   P2 should be a scalar integer. Although the CRC has P bits, only 
%   P-min(P2,log2(L)) of these are used for error detection. The remaining 
%   min(P2,log2(L)) of the CRC bits are used to improve error correction. 
%   So the CRC needs to be min(P2,log2(L)) number of bits longer than CRCs 
%   used in other codes, in order to achieve the same error detection 
%   capability.
%
%   a_tilde should be a row vector comprising A elements. Elements having
%   the value 0 or 1 indicate that the receiver has prior knowledge that
%   the corresponding information bit has the corresponding value. Elements
%   having any other value indicate that the receiver does not have prior
%   knowledge of the value of the corresponding information bit.
%
%   a_hat will normally be a binary row vector comprising A number of bits, 
%   each having the value 0 or 1. However, in cases where the CRC check 
%   fails, a_hat will be an empty vector.
%
%   See also DS1CA_POLAR_ENCODER
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

E = length(e_tilde);
N = length(info_bit_pattern);
K = length(crc_interleaver_pattern);
P = length(crc_polynomial_pattern)-1;
A = K-P;

if log2(N) ~= round(log2(N))
    error('N should be a power of 2');
end
if sum(info_bit_pattern) ~= K
    error('info_bit_pattern should contain K number of ones.');
end
if max(rate_matching_pattern) > N
    error('rate_matching_pattern is not compatible with N');
end
if strcmp(mode,'repetition') 
    if E < N
        error('mode is not compatible with E');
    end
elseif strcmp(mode,'puncturing')
    if E >= N
        error('mode is not compatible with E');
    end
elseif strcmp(mode,'shortening')
     if E >= N
        error('mode is not compatible with E');
    end
else
    error('Unsupported mode');
end
if P < length(crc_scrambling_pattern)
    error('polar_3gpp_matlab:UnsupportedBlockLength','P should be no less than the length of the scrambing pattern');
end
if length(a_tilde) ~= A
    error('a_tilde should contain A number of elements.');
end

% This global variable is used by the minstar and phi functions.
global approx_minstar
approx_minstar=min_sum;

%% Characterise the interleaved CRC

% Get the CRC generator matrix, which has dimensions A+P by P.
G_P = get_crc_generator_matrix(A+P,crc_polynomial_pattern);

% Extend the CRC generator matrix by append an identity matrix to 
% represent the CRC bits, giving dimenstions K by P.
G_P2 = [G_P(P+1:end,:);eye(P)];

% Interleave the rows of the extended CRC generator matrix, according to
% the CRC interleaver.
G_P3 = G_P2(crc_interleaver_pattern,:);

% Determine where the last 1-valued bit appears in each column. When the 
% SCL decoding process reaches the corresponding interleaved CRC bit, we 
% will terminate the decoding process if all list entries have failed the 
% checks assocated with at least one of this or the preceeding CRC bits.
last_one_index = zeros(1,P);
for p = 1:P
    last_one_index(p) = find(G_P3(:,p) == 1, 1, 'last');
end
 
% Extend the scrambling pattern to match the length of the CRC
extended_crc_scrambling_pattern = [zeros(1,P-length(crc_scrambling_pattern)), crc_scrambling_pattern];

%% Rate matching
if strcmp(mode,'repetition')
    % LLRs for repeated bits are added together.
    d_tilde = zeros(1,N);
    for i=1:E
        d_tilde(rate_matching_pattern(i)) = d_tilde(rate_matching_pattern(i)) + e_tilde(i);
    end
else
    if strcmp(mode,'puncturing')
        % Zero valued LLRs are used for punctured bits, because the decoder
        % doesn't know if they have values of 0 or 1.
        d_tilde = zeros(1,N);
    elseif strcmp(mode,'shortening')
        % Infinite valued LLRs are used for shortened bits, because the 
        % decoder knows that they have values of 0.
        d_tilde = inf(1,N);
    else
        error('Unknown rate matching mode');
    end
    
    d_tilde(rate_matching_pattern) = e_tilde;
end

%% Perform the SCL polar decoder kernal operation.
% This is achieved according to the algorithm described in 
% http://ieeexplore.ieee.org/abstract/document/7114328/

global bits % Matrix to store bit throughout the polar code graph
global bits_updated % Matrix to identify which bits in the polar code graph have been calculated so far
global llrs % Matrix to store LLRs throughout the polar code graph
global llrs_updated % Matrix to identify which LLRs in the polar code graph have been calculated so far

bits = zeros(N, log2(N)+1); % Initialse all bits to zero. The left-most column corresponds to the decoded information, CRC and frozen bits
bits_updated = [~info_bit_pattern',false(N,log2(N))]; % The zero values that have initialised the frozen bits are known to be correct
llrs = [zeros(N,log2(N)),d_tilde']; % Initialse the LLRs. The right-most column corresponds to the received LLRS
llrs_updated = [false(N,log2(N)),true(N,1)]; % The received LLRs have been updated.

PM = zeros(1,1,1); % Initialise the path metrics
L_prime = 1; % Initialise the list size to 1. This will grow as the decoding proceeds

% We will calculate the CRC checksums alongside the SCL decoding process.
% Initialise the checksums with P number of 1s.
crc_checksums = mod(sum(G_P(1:P,:)),2)';
% We will keep track of whether any of the checks associated with the CRC
% bits have failed.
crc_okay = true;
% We need a counter to keep track of which information or CRC bit we are
% working on.
i2=1;
% We will return an empty vector if all list entries have failed the checks 
% assocated with at least one of the CRC bits.
a_hat = [];

% Consider each bit in turn
for i = 1:N
    % Make recursive function calls to perform the XOR, g and f functions
    % necessary to obtain the corresponding LLR
    update_llr(i,1);
    
    if info_bit_pattern(i) == 0 % Frozen bit
        PM = phi(PM, llrs(i,1,:), 0); 
    else % Information or CRC bit
        if crc_interleaver_pattern(i2) <= A && a_tilde(crc_interleaver_pattern(i2)) == 0 % Information bit with known value of 0
            PM = phi(PM, llrs(i,1,:), 0);         
            bits_updated(i,1) = true;
       elseif crc_interleaver_pattern(i2) <= A && a_tilde(crc_interleaver_pattern(i2)) == 1 % Information bit with known value of 1
            PM = phi(PM, llrs(i,1,:), 1); 
            bits(i,1,:) = 1;
            bits_updated(i,1) = true;        
            % We use the interleaved CRC generator matrix to update the CRC 
            % check sums whenever an information bit adopts a value of 1.
            crc_checksums = mod(crc_checksums+repmat(G_P3(i2,:)',[1 1 L_prime]),2);           
        else
            % Double the list size, using 0-valued bits for the first half and 1-valued bits for the other half
            PM = cat(3,phi(PM, llrs(i,1,:), 0), phi(PM, llrs(i,1,:), 1));
            llrs = cat(3,llrs,llrs);
            bits = cat(3,bits,bits);
            bits(i,1,1:L_prime) = 0;
            bits(i,1,L_prime+1:2*L_prime) = 1;
            bits_updated(i,1) = true;

            % We use the interleaved CRC generator matrix to update the CRC 
            % check sums whenever an information or CRC bit adopts a value of
            % 1.
            crc_checksums = cat(3,crc_checksums,mod(crc_checksums+repmat(G_P3(i2,:)',[1 1 L_prime]),2));
            % We need to keep track of whether any of the checks associated 
            % with the previous CRC bits have failed.
            crc_okay = cat(3,crc_okay,crc_okay);

            % If the list size has grown above L, then we need to find and keep only the best L entries in the list
            L_prime = size(bits,3);        
            if L_prime > L
                [~,max_indices] = sort(PM,3);
                PM = PM(:,:,max_indices(1:L));
                bits = bits(:,:,max_indices(1:L));
                llrs = llrs(:,:,max_indices(1:L));
                crc_checksums = crc_checksums(:,:,max_indices(1:L));
                crc_okay = crc_okay(:,:,max_indices(1:L));
                L_prime = L;
            end
        end

        % We check the corresponding CRC checksums whenever we reach the
        % last 1-valued bit in a column of the interleaved CRC generator
        % matrix.
        check_crc_bits = find(last_one_index == i2);
        for crc_bit_index = 1:length(check_crc_bits)
            for list_index = 1:L_prime
                % The checksum should equal the value of the corresponding 
                % CRC scrambling pattern bit. If not, then the CRC
                % check has failed. Note that we should not prune these
                % entries from the list, even though we know that they will
                % fail the CRC. We should continue the decoding of this
                % list entries, otherwise we will damage the error
                % detection capability of the CRC.
                if crc_checksums(check_crc_bits(crc_bit_index),1,list_index) ~= extended_crc_scrambling_pattern(check_crc_bits(crc_bit_index))
                    % We keep track of the failing check.
                    crc_okay(1,1,list_index) = false;
                end
            end
        end
        
        % We terminate the decoding process if all list entries have failed 
        % the checks assocated with at least one of this or the preceeding 
        % CRC bits.
        if ~any(crc_okay)
            return;
        end
        
        % Increment the counter of information and CRC bits
        i2 = i2+1;
    end
end

%% Information bit extraction

% We use the list entry with a passing CRC that has the best metric. But we
% only consider the best min(L,2^P2) entries in the list, to avoid
% degrading the CRC's error detection capability.
[~,max_indices] = sort(PM,3);
b_hat = zeros(1,K);
for list_index = 1:min(L,2^P2)
    % Consider the next best list entry.

    % We already checked the CRC during the SCL decoding process.
    if crc_okay(max_indices(list_index))
        u_hat = bits(:,1,max_indices(list_index))';
 
        % Extract the information bits from the output of the polar decoder 
        % kernal.
        c_hat = u_hat(info_bit_pattern);

        % Deinterleave the information and CRC bits.
        b_hat(crc_interleaver_pattern) = c_hat;
        
        % Remove the CRC and output the information bits.
        a_hat = b_hat(1:end-P);
        return;
    end
end
end