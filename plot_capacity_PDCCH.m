addpath '../spectre/bi-awgn'

A = [1:15, 16:2:30, 32:4:60, 64:8:120, 128:16:140];
E = [108 216 432 864 1728];


K = A+21;

target_BLER = 0.001;

% Create a figure to plot the results.
figure
axes1 = axes;
ylabel('Required E_s/N_0 [dB]');
xlabel('A');
xt = 0:11;
set(gca, 'XTick', xt);
set (gca, 'XTickLabel', 2.^xt);
grid on
hold on
drawnow

for E_index = 1:length(E)
    
    % Create the plot
    plots(E_index) = plot(nan,'Parent',axes1);
    legend(cellstr(num2str(E(1:E_index)', 'E=%d, capacity')),'Location','eastoutside');

    % Open a file to save the results into.
    filename = ['results/SNR_vs_A_PDCCH_',num2str(target_BLER),'_',num2str(E(E_index)),'_cap'];
    fid = fopen([filename,'.txt'],'w');
    if fid == -1
        error('Could not open %s.txt',filename);
    end
    
    
    EsN0s = nan(1,length(A));
    
    plot_EsN0s = [];
    plot_As = [];
    for A_index = 1:length(A)
        if K(A_index)+3 <= E(E_index)
            EsN0 = converse_mc(E(E_index), target_BLER, K(A_index)/E(E_index),'On2','error');
            if isempty(plot_EsN0s) || EsN0 > plot_EsN0s(end)
                plot_EsN0s(end+1) = EsN0;
                plot_As(end+1) = A(A_index);
                
                % Plot the SNR vs A results
                set(plots(E_index),'XData',log2(plot_As));
                set(plots(E_index),'YData',plot_EsN0s);
                
                xlim auto;
                xl = xlim;
                xlim([floor(xl(1)), ceil(xl(2))]);
                
                drawnow;
                
                fprintf(fid,'%d\t%f\n',A(A_index),EsN0);
                
            end
            
        end
    end
    % Close the file
    fclose(fid);

end