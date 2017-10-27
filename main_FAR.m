function main_FAR(code, A, E, L, min_sum, target_false_alarms, seed)
% PLOT_BLER_VS_SNR Plots Block Error Rate (BLER) versus Signal to Noise
% Ratio (SNR) for polar codes.
%   plot_BLER_vs_SNR(code, A, E, L, min_sum, target_block_errors, target_BLER, EsN0_start, EsN0_delta, seed)
%   generates the plots.
%
%   code should be a string. This identifies which encoder and decoder
%   functions to call. For example, if code is 'custom1', then the
%   functions custom1_encoder and custom1_decoder will be called. The
%   encoder function should have a format f = custom1_encoder(a, E). The
%   decoder function should have a format
%   a_hat = custom1_decoder(f_tilde, A, L, min_sum). Refer to these
%   functions for explanations of their inputs and outputs. Suitable values
%   for code include 'custom1' and 'PBCH'.
%
%   A should be an integer scalar. It specifies the number of bits in each
%   simulated information bit sequence, before CRC and other redundant bits
%   are included.
%
%   E should be an integer row vector. Each element of E specifies one
%   encoded block length to simulate, where E is the number of bits in each
%   encoded bit sequence.
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
%   target_block_errors should be an integer scalar. The simulation of each
%   SNR for each coding rate will continue until this number of block
%   errors have been observed. A value of 100 is sufficient to obtain
%   smooth BLER plots for most values of A. Higher values will give
%   smoother plots, at the cost of requiring longer simulations.
%
%   target_BLER should be a real scalar, in the range (0, 1). The
%   simulation of each coding rate will continue until the BLER plot
%   reaches this value.
%
%   EsN0_start should be a real row vector, having the same length as the
%   vector of coding rates. Each value specifies the Es/N0 SNR to begin at
%   for the simulation of the corresponding coding rate.
%
%   EsN0_delta should be a real scalar, having a value greater than 0.
%   The Es/N0 SNR is incremented by this amount whenever
%   target_block_errors number of block errors has been observed for the
%   previous SNR. This continues until the BLER reaches target_BLER.
%
%   seed should be an integer scalar. This value is used to seed the random
%   number generator, allowing identical results to be reproduced by using
%   the same seed. When running parallel instances of this simulation,
%   different seeds should be used for each instance, in order to collect
%   different results that can be aggregated together.
%
%   See also CUSTOM1_ENCODER and CUSTOM1_DECODER
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you
% can redistribute it and/or modify it under the terms of the GNU General
% Public License as published by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version. This
% program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.

% Default values
if nargin == 0
    code = 'custom1';
    A = 16;
    E = round(A./[0.8333 0.7500 0.6666 0.5000 0.4000 0.3333 0.2500 0.2000 0.1666 0.1250]);
    L = 1;
    min_sum = true;
    target_false_alarms = 10;
    seed = 0;
end

% Seed the random number generator
rng(seed);

% Consider each coding rate in turn
for E_index = 1:length(E)
    
    % Counters to store the number of bits and errors simulated so far
    block_count=0;
    false_alarm_count=0;
    
    % Skip any encoded block lengths that generate errors
    try
        
        % Continue the simulation until enough block errors have been simulated
        while false_alarm_count(end) < target_false_alarms
            
            f_tilde = randn(1,E(E_index));
            
            % Perform polar decoding
            a_hat = feval([code, '_decoder'],f_tilde,A,L,min_sum);

            
            % Accumulate the number of blocks that have been simulated
            % so far
            block_count = block_count + 1;
 
            % Determine if we have a block error
            if length(a_hat) == A
                false_alarm_count = false_alarm_count + 1;
                
                [false_alarm_count, block_count, false_alarm_count/block_count]
            end
            
            
            
            
            
            
        end
    catch ME
        % Issue any errors as warnings and move on to the next encoded block length.
        warning(['E=',num2str(E(E_index)),' - ',ME.message]);
    end
    
end





