function a_hat = PC_polar_decoder(e_tilde, info_bit_pattern, PC_bit_pattern, PC_circular_buffer_length, rate_matching_pattern, mode, L, min_sum)
% PCCA_POLAR_DECODER Parity Check and CRC-Aided (CA) polar decoder.
%   a_hat = PCCA_POLAR_DECODER(e_tilde, crc_polynomial_pattern, info_bit_pattern, PC_bit_pattern, PC_circular_buffer_length, rate_matching_pattern, mode, L, min_sum, P2) 
%   decodes the encoded LLR sequence e_tilde, in order to obtain the
%   recovered information bit sequence a_hat.
%
%   e_tilde should be a real row vector comprising E number of Logarithmic
%   Likelihood Ratios (LLRS), each having a value obtained as LLR =
%   ln(P(bit=0)/P(bit=1)).
%
%   info_bit_pattern should be a row vector comprising N number of logical 
%   elements, each having the value true or false. The number of elements 
%   in info_bit_pattern having the value true should be A+n_PC. 
%   These elements having the value true identify the positions of the 
%   information and PC bits within the input to the polar encoder kernal.
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
%   a_hat will be a binary row vector comprising A number of bits, 
%   each having the value 0 or 1.
%
%   See also PCCA_POLAR_ENCODER
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
A = sum(info_bit_pattern & ~PC_bit_pattern);

if log2(N) ~= round(log2(N))
    error('N should be a power of 2');
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

% This global variable is used by the minstar and phi functions.
global approx_minstar
approx_minstar=min_sum;

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

bits = zeros(N, log2(N)+1); % Initialse all bits to zero. The left-most column corresponds to the decoded information, PC and frozen bits
bits_updated = [~info_bit_pattern',false(N,log2(N))]; % The zero values that have initialised the frozen bits are known to be correct
llrs = [zeros(N,log2(N)),d_tilde']; % Initialse the LLRs. The right-most column corresponds to the received LLRS
llrs_updated = [false(N,log2(N)),true(N,1)]; % The received LLRs have been updated.

PM = zeros(1,1,1); % Initialise the path metrics
L_prime = 1; % Initialise the list size to 1. This will grow as the decoding proceeds

y = zeros(1,PC_circular_buffer_length);



% Consider each bit in turn
for i = 1:N
    % Rotate the PC bit generators
    y = [y(1,2:end,:),y(1,1,:)];
    
    % Make recursive function calls to perform the XOR, g and f functions
    % necessary to obtain the corresponding LLR
    update_llr(i,1);
    
    if info_bit_pattern(i) == 0 % Frozen bit
        PM = phi(PM, llrs(i,1,:), 0);
    elseif PC_bit_pattern(i) == 1 % PC bit
        % The bit values come from the PC bit generator
        bits(i,1,:) = y(1,1,:);
        PM = phi(PM, llrs(i,1,:), bits(i,1,:));
        bits_updated(i,1) = true;
    else % Information bit
        % Double the list size, using 0-valued bits for the first half and 1-valued bits for the other half
        PM = cat(3,phi(PM, llrs(i,1,:), 0), phi(PM, llrs(i,1,:), 1));
        llrs = cat(3,llrs,llrs);
        bits = cat(3,bits,bits);
        bits(i,1,1:L_prime) = 0;
        bits(i,1,L_prime+1:2*L_prime) = 1;
        bits_updated(i,1) = true;
        
        % Update the PC bit generators
        y = cat(3,y,y);
        y(1,1,:) = xor(y(1,1,:),bits(i,1,:));
        
        % If the list size has grown above L, then we need to find and keep only the best L entries in the list
        L_prime = size(bits,3);        
        if L_prime > L
            [~,max_indices] = sort(PM,3);
            PM = PM(:,:,max_indices(1:L));
            bits = bits(:,:,max_indices(1:L));
            llrs = llrs(:,:,max_indices(1:L));
            y = y(:,:,max_indices(1:L));
            L_prime = L;
        end
    end
end

%% Information bit extraction

% We use the list entry that has the best metric
[~,max_indices] = sort(PM,3);
u_hat = bits(:,1,max_indices(1))';

% Extract the information bits from the output of the polar decoder kernal.
a_hat = u_hat(info_bit_pattern & ~PC_bit_pattern);

end