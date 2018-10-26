% This is the strong scaling - i.e. how the execution time varies with
% number of processors for fixed total problem size.


clear all;

close all;
baseDir = getDataDir('parallelTest/');
%prefix = 'CR6.0RaC35Le200ChiCubedPermeabilitypts128';
%prefix = 'CR6.0RaC5000Le200KozenyPermeabilitypts128';

prefixes = {'CR6.0RaC35Le200ChiCubedPermeabilitypts128', ...
   % 'CR6.0RaC5000Le200KozenyPermeabilitypts128', ...
    'CR6.0RaC5000Le200KozenyPermeabilitypts128', ...
    %'CR20RaC-20000.0Le200KozenyPermeabilitypts512',...
    'CR20RaC-20000.0Le200KozenyPermeabilitypts512'};

leg = {'Darcy 128x256', ...
    'Darcy-Brinkman 128x256',...
    'Darcy-Brinkman 512x512'};

nodes = [1 1 1];

nproc = [1 2 4 8 16];
%nproc = [1];
t = NaN*ones(length(prefixes), length(nproc));
for i = 1:length(nproc)
    for p_i = 1:length(prefixes)
        
    if nodes(p_i) > 1
        node_str = ['NODES', num2str(nodes(p_i))];
    else
        node_str = '';
    end
    
    folder = [baseDir, prefixes{p_i}, 'NPROC', num2str(nproc(i)), node_str, '-0/'];
    
    timingsFile = [folder, 'time.table.0'];
    
    test = exist(timingsFile, 'file');
    
    if exist(timingsFile, 'file') == 2
        %fprintf(['Loading', file, '\n']);
        file_contents = fileread(timingsFile);
    else
        continue;
    end
    
    
    pattern = ['main[^\n]*\n', ...
        '[^\n]*\n', ...
        '(^\n)*\n'];
    
    pattern = ['main[^\n]*[\n]',...
        '([^AMR::run]*)'];
    %   57.0%   19.07869        1 AMR::conclude [1] f=0
    %   37.5%   12.54402        1 AMR::run [4] f=102750144
    %    0.3%    0.10057        1 AMR::setupForRestart [164] f=0
    %    0.2%    0.06587        1 AMR::define [189] f=0
    %    0.0%    0.00708        2 AMR::clearMemory [586] f=0
    %    0.0%    0.00000        1 AMR::blockFactor [2772] f=0
    %    0.0%    0.00000        1 AMR::maxGridSize [2858] f=0
    %    0.0%    0.00000        1 AMR::fillRatio [2861] f=0
    %    0.0%    0.00000        1 AMR::gridBufferSize [2862] f=0
    %   95.0%                  Total
    %---------------------------------------------------------'];
    
    [mat,tok,ext]  = regexp(file_contents, pattern, 'match', ...
        'tokens', 'tokenExtents');
    
    if length(tok) > 0
        
        match = tok{1}{1};
        
        temp = 0;
        parts = strsplit(match, '%');
        parts2 = strsplit(strtrim(parts{2}), ' ');
        seconds = str2num(parts2{1});
        
        t(p_i, i) = seconds;
    end
    
    end % end loop over prefixes
end % loop over nproc


t_perfect = t;
for i = 2:length(nproc)
    t_perfect(:, i) = t_perfect(:, i-1)*nproc(i-1)/nproc(i);
    
end

h = figure();
set(h, 'Position', [300 300 900 600]);
hold on;
plot(log10(nproc), log10(t), '-x');

ax = gca;
ax.ColorOrderIndex = 1;

plot(log10(nproc), log10(t_perfect), '--');
hold off;
xlabel('log10($N_{proc}$)');
ylabel('log10($t$)');
title('Execution time in AMR::run');

box on;


for i=1:length(prefixes)
    leg{end+1} = 'Perfect scaling';
end

legend(leg, 'Location', 'eastoutside');
