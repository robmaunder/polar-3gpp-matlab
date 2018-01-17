addpath '../spectre/bi-awgn'

A = 32;
K = A + 21;
E = 864;

    target_BLER = 1e-3;
    points = 100;


% Create a figure to plot the results.
figure
axes1 = axes('YScale','log');
ylabel('BLER');
xlabel('E_s/N_0 [dB]');
ylim([target_BLER,1]);
hold on
drawnow

    % Create the plot
    plot1 = plot(nan,'Parent',axes1);
    legend('capacity','Location','southwest');

    % Open a file to save the results into.
    filename = ['results/BLER_vs_SNR_PBCH_',num2str(A),'_',num2str(E),'_cap'];
    fid = fopen([filename,'.txt'],'w');
    if fid == -1
        error('Could not open %s.txt',filename);
    end

    BLERs = 10.^(log10(0.1):(log10(target_BLER)-log10(0.1))/points:log10(target_BLER));
    EsN0s = nan(size(BLERs));
for BLER_index = 1:length(BLERs)
    EsN0s(BLER_index) = converse_mc(E, BLERs(BLER_index), K/E, 'On2','error');
                    fprintf(fid,'%f\t%e\n',EsN0s(BLER_index),BLERs(BLER_index));
                    % Plot the BLER vs SNR results
                    set(plot1,'XData',EsN0s);
                    set(plot1,'YData',BLERs);                
                    drawnow

end

fclose(fid);
