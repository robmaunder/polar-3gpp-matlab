function main_FAR(code, A, E, L, min_sum, target_false_alarms, seed)
% main_FAR calculates False Alarm Rate (FAR) for polar codes.
%   main_FAR(code, A, E, L, min_sum, target_false_alarms, seed)
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
%   target_false_alarms should be an integer scalar. The simulation of each
%   encoded block length will continue until this number of false
%   alarms have been observed. A value of 100 is sufficient to obtain
%   smooth BLER plots for most values of A. Higher values will give
%   smoother plots, at the cost of requiring longer simulations.
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
    code = 'PUCCH';
    A = 32;
    E = [54 108 216 432 864 1728 3456 6912 13824];
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





