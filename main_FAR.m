function main_FAR(code, A, E, L, min_sum, target_false_alarms, seed)
% main_FAR calculates False Alarm Rate (FAR) for polar codes.
%   main_FAR(code, A, E, L, min_sum, target_false_alarms, seed) calculates
%   the FAR.
%
%   code should be a string. This identifies which encoder and decoder
%   functions to call. For example, if code is 'custom1', then the
%   functions custom1_encoder and custom1_decoder will be called. The
%   encoder function should have a format f = custom1_encoder(a, E). The
%   decoder function should have a format
%   a_hat = custom1_decoder(f_tilde, A, L, min_sum). Refer to these
%   functions for explanations of their inputs and outputs. Suitable values
%   for code include 'PBCH', 'PDCCH, 'PUCCH' and 'custom1'.
%
%   A should be an integer row vector. Each element specifies the number of
%   bits in each set of simulated information bit sequences, before CRC and
%   other redundant bits are included.
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
%   reliable FAR results for most values of A. Higher values will give
%   more reliable results, at the cost of requiring longer simulations.
%
%   seed should be an integer scalar. This value is used to seed the random
%   number generator, allowing identical results to be reproduced by using
%   the same seed. When running parallel instances of this simulation,
%   different seeds should be used for each instance, in order to collect
%   different results that can be aggregated together.
%
%   See also MAIN_BLER_VS_SNR and MAIN_SNR_VS_A
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
    A = [16 32 64 128 256 512 1024];
    E = [54 108 216 432 864 1728 3456 6912 13824];
    L = 1;
    min_sum = true;
    target_false_alarms = 10;
    seed = 0;
end

% Open a file to save the results into.
filename = ['results/FAR_',code,'_',num2str(L),'_',num2str(min_sum),'_',num2str(target_false_alarms),'_',num2str(seed)];
fid = fopen([filename,'.txt'],'w');
if fid == -1
    error('Could not open %s.txt',filename);
end


% Seed the random number generator
rng(seed);

fprintf("   A          E        FAs     blocks        FAR\n");

% Consider each information block length in turn
for A_index = 1:length(A)
    
    % Consider each encoded block length in turn
    for E_index = 1:length(E)
        
        % Counters to store the number of bits and errors simulated so far
        block_count=0;
        false_alarm_count=0;
        
        chars_to_erase=0;
        
        % Skip any encoded block lengths that generate errors
        try
            
            % Continue the simulation until enough block errors have been simulated
            while false_alarm_count(end) < target_false_alarms
                
                % Use Gaussian distributed random LLRs
                f_tilde = randn(1,E(E_index));
                
                % Perform polar decoding
                a_hat = feval([code, '_decoder'],f_tilde,A(A_index),L,min_sum);
                                
                % Accumulate the number of blocks that have been simulated
                % so far
                block_count = block_count + 1;
                
                % If the CRC is satisfied, then we have a false alarm
                if length(a_hat) == A(A_index)
                    false_alarm_count = false_alarm_count + 1;
                    
                    fprintf(repmat('\b',1,chars_to_erase));                    
                    msg = sprintf('%4d %10d %10d %10d  %.3e\n', A(A_index), E(E_index), false_alarm_count, block_count, false_alarm_count/block_count);                    
                    fprintf(msg);
                    chars_to_erase = numel(msg);
                end
            end
            
            fprintf(fid,'%d\t%d\t%e\n',A(A_index),E(E_index),false_alarm_count/block_count);

        catch ME
            if strcmp(ME.identifier, 'polar_3gpp_matlab:UnsupportedBlockLength')
                warning('polar_3gpp_matlab:UnsupportedBlockLength','%s does not support the combination of block lengths A=%d and E=%d. %s',code,A(A_index),E(E_index), getReport(ME, 'basic', 'hyperlinks', 'on' ));
                continue
            else
                rethrow(ME);
            end
        end        
    end
end

fclose(fid);


