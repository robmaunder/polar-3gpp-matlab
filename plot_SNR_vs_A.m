function plot_SNR_vs_A(A, R, L, n_max, new_approx_minstar, target_block_errors, target_BLER, EsN0_start, EsN0_delta, seed)
% Copyright © 2017 Robert G. Maunder. This program is free software: you 
% can redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version. This 
% program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
% more details.

if nargin == 0
    A = 16:16:768;
    R = [0.8333 0.7500 0.6666 0.5000 0.4000 0.3333 0.2500 0.2000 0.1666 0.1250];
    L = 1;
    n_max = 9;
    new_approx_minstar = 1;
    target_block_errors = 100;
    target_BLER = 1e-1;
    EsN0_start = [-1, -2, -3, -4, -5, -6, -7, -8, -9, -10];
    EsN0_delta = 0.5;
    seed = 0;
end


global approx_minstar
approx_minstar=new_approx_minstar;

rng(seed);


figure
axes1 = axes;
title(['L = ', num2str(L), ' n_{max} = ',num2str(n_max)]);
ylabel(['E_s/N_0 [dB] at BLER = ',num2str(target_BLER)]);
xlabel('A');
xlim([0,max(A)]);
grid on
hold on
drawnow

crc_polynomial_pattern = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];
P = length(crc_polynomial_pattern)-1;
L2 = 8;

for R_index = 1:length(R)
    plots(R_index) = plot(nan,'Parent',axes1);
    set(plots(R_index),'XData',A);
    legend(cellstr(num2str(R(1:R_index)', 'R=%0.3f')),'Location','eastoutside');
    
    EsN0s = nan(1,length(A));
    
    for A_index = 1:length(A)
        
        K = A(A_index)+P;
        
        E = round(A(A_index)/R(R_index));
        if K <= E
            N = get_3GPP_N(K,E,n_max);
            crc_interleaver_pattern = get_3GPP_crc_interleaver_pattern(K);            
            info_bit_pattern = get_3GPP_info_bit_pattern(K,N,E);
            [rate_matching_pattern, mode] = get_3GPP_rate_matching_pattern(K,N,E);
            
            
            
            BLER=1;
            prev_BLER = nan;
            EsN0 = EsN0_start(R_index);
            prev_EsN0 = nan;
            
            
            while BLER > target_BLER
                
                block_errors = 0;
                blocks = 0;
                while block_errors < target_block_errors
                    
                    % Generate a random frame of bits
                    a = round(rand(1,A(A_index)));
                    
                    % Perform polar encoding
%                    e = CA_polar_encoder(a,crc_polynomial_pattern, info_bit_pattern, rate_matching_pattern);
                    e = DCA_polar_encoder(a,crc_polynomial_pattern,crc_interleaver_pattern,info_bit_pattern,rate_matching_pattern);
                    
                    % QPSK modulation
                    e2 = [e,zeros(1,mod(-length(e),2))];
                    tx = sqrt(1/2)*(2*e2(1:2:end)-1)+1i*sqrt(1/2)*(2*e2(2:2:end)-1);
                    
                    % Simulate transmission
                    N0 = 1/(10^(EsN0/10));
                    rx = tx + sqrt(N0/2)*(randn(size(tx))+1i*randn(size(tx)));
                    
                    % QPSK demodulation
                    e2_tilde = zeros(size(e2));
                    e2_tilde(1:2:end) = -4*sqrt(1/2)*real(rx)/N0;
                    e2_tilde(2:2:end) = -4*sqrt(1/2)*imag(rx)/N0;
                    e_tilde = e2_tilde(1:length(e));
                    
                    % Perform polar decoding
%                    a_hat = CA_polar_decoder(e_tilde,crc_polynomial_pattern, info_bit_pattern, rate_matching_pattern, mode, L, L2);
                    a_hat = DCA_polar_decoder(e_tilde,crc_polynomial_pattern,crc_interleaver_pattern,info_bit_pattern,rate_matching_pattern,mode,L,L2);
                    
                    if ~isequal(a, a_hat)
                        block_errors = block_errors+1;
                    end
                    blocks = blocks+1;
                end
                prev_BLER = BLER;
                BLER = block_errors/blocks;
                prev_EsN0 = EsN0;
                EsN0 = EsN0 + EsN0_delta;
            end
            EsN0s(A_index) = interp1(log10([prev_BLER, BLER]),[prev_EsN0,EsN0],log10(target_BLER));
            set(plots(R_index),'YData',EsN0s);
            drawnow;
        end
    end
end
