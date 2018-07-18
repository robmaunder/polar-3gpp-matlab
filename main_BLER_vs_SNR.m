function main_BLER_vs_SNR(code, A, E, L, min_sum, target_block_errors, target_BLER, EsN0_start, EsN0_delta, seed)
% MAIN_BLER_VS_SNR Plots Block Error Rate (BLER) versus Signal to Noise
% Ratio (SNR) for polar codes.
%   main_BLER_vs_SNR(code, A, E, L, min_sum, target_block_errors, target_BLER, EsN0_start, EsN0_delta, seed)
%   generates the plots.
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
%   min_sum should be a scalar logical. If it is true, then the SCL
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
%   See also MAIN_SNR_VS_A and MAIN_FAR
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
    target_block_errors = 10;
    target_BLER = 1e-1;
    EsN0_start = -25;
    EsN0_delta = 0.5;
    seed = 0;
end

% Seed the random number generator
rng(seed);

% Consider each information block length in turn
for A_index = 1:length(A)
    
    % Create a figure to plot the results.
    figure
    axes1 = axes('YScale','log');
    title([code, ' polar code, A = ',num2str(A(A_index)),', L = ',num2str(L),', minsum = ',num2str(min_sum),', errors = ',num2str(target_block_errors),', QPSK, AWGN']);
    ylabel('BLER');
    xlabel('E_s/N_0 [dB]');
    ylim([target_BLER,1]);
    hold on
    drawnow
    
    % Consider each encoded block length in turn
    for E_index = 1:length(E)
        
        % Create the plot
        plot1 = plot(nan,'Parent',axes1);
        legend(cellstr(num2str(E(1:E_index)', 'E=%d')),'Location','southwest');
        
        % Counters to store the number of bits and errors simulated so far
        block_counts=[];
        block_error_counts=[];
        EsN0s = [];
        
        % Open a file to save the results into.
        filename = ['results/BLER_vs_SNR_',code,'_',num2str(A(A_index)),'_',num2str(E(E_index)),'_',num2str(L),'_',num2str(min_sum),'_',num2str(target_block_errors),'_',num2str(seed)];
        fid = fopen([filename,'.txt'],'w');
        if fid == -1
            error('Could not open %s.txt',filename);
        end
        
        % Initialise the BLER and SNR
        BLER = 1;
        EsN0 = EsN0_start;
        
        found_start = false;
        
        % Skip any encoded block lengths that generate errors
        try
            % Loop over the SNRs
            while BLER > target_BLER
                
                % Convert from SNR (in dB) to noise power spectral density
                N0 = 1/(10^(EsN0/10));
                
                % Start new counters
                block_counts(end+1) = 0;
                block_error_counts(end+1) = 0;
                EsN0s(end+1) = EsN0;
                
                keep_going = true;
                
                % Continue the simulation until enough block errors have been simulated
                while keep_going && block_error_counts(end) < target_block_errors
                    
                    % Generate a random block of bits
                    a = round(rand(1,A(A_index)));
                    
                    % Perform polar encoding
                    f = feval([code,'_encoder'], a, E(E_index));
                    
                    % QPSK modulation
                    f2 = [f,zeros(1,mod(-length(f),2))];
                    tx = sqrt(1/2)*(2*f2(1:2:end)-1)+1i*sqrt(1/2)*(2*f2(2:2:end)-1);
                    
                    % Simulate transmission
                    rx = tx + sqrt(N0/2)*(randn(size(tx))+1i*randn(size(tx)));
                    
                    % QPSK demodulation
                    f2_tilde = zeros(size(f2));
                    f2_tilde(1:2:end) = -4*sqrt(1/2)*real(rx)/N0;
                    f2_tilde(2:2:end) = -4*sqrt(1/2)*imag(rx)/N0;
                    f_tilde = f2_tilde(1:length(f));
                    
                    % Perform polar decoding
                    a_hat = feval([code, '_decoder'],f_tilde,A(A_index),L,min_sum);
                    
                    
                    if found_start == false && ~isequal(a,a_hat)
                        keep_going = false;
                        BLER = 1;
                    else
                        found_start = true;
                        
                        % Determine if we have a block error
                        if ~isequal(a,a_hat)
                            block_error_counts(end) = block_error_counts(end) + 1;
                        end
                        
                        % Accumulate the number of blocks that have been simulated
                        % so far
                        block_counts(end) = block_counts(end) + 1;
                        
                        % Calculate the BLER and save it in the file
                        BLER = block_error_counts(end)/block_counts(end);
                        
                        % Plot the BLER vs SNR results
                        set(plot1,'XData',EsN0s);
                        set(plot1,'YData',block_error_counts./block_counts);
                        drawnow
                    end
                end
                
                if BLER < 1
                    fprintf(fid,'%f\t%e\n',EsN0,BLER);
                end
                
                % Update the SNR, ready for the next loop
                EsN0 = EsN0 + EsN0_delta;
                
            end
        catch ME
            if strcmp(ME.identifier, 'polar_3gpp_matlab:UnsupportedBlockLength')
                warning('polar_3gpp_matlab:UnsupportedBlockLength','%s does not support the combination of block lengths A=%d and E=%d. %s',code,A(A_index),E(E_index), getReport(ME, 'basic', 'hyperlinks', 'on' ));
                continue
            else
                rethrow(ME);
            end
        end
        
        % Close the file
        fclose(fid);
    end
end




