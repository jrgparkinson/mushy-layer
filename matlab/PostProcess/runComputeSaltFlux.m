
close all; clear all;

data_dir = getDataDir('middleton/');
%simulation = 'CR1.179RaC500Le200ChiCubedPermeabilitypts256-S3.5-TB-15.0-domWidth1.0-0';
%output_dir= [data_dir, simulation, '/'];
prefix = 'plt';
frames = 100:500:11000;

% sims = {'CR1.179RaC8Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.00013-domWidth2.5-stoppedEarly', ...
%     'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.0013-domWidth2.5-stoppedEarly',...
%     'CR1.179RaC8Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.013-domWidth2.5-stoppedEarlyVerySmooth',...
%     'CR1.179RaC800Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-domWidth2.5-stoppedEarlyBitDodgy', ...
%     'CR1.179RaC800Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R1.3-domWidth2.5-stoppedEarlyBitDodgy', ...
%     'CR1.179RaC800Le200ChiCubedPermeabilitypts512-S3.5-TB-20.0-R0.013-domWidth2.5-stoppedEarlyGood', ...
%     'CR1.179RaC80Le200KozenyPermeabilitypts512-S3.5-TB-20.0-R0.0013-domWidth2.5-stoppedEarlyBitDodgy',...
%     'CR1.179RaC80Le200KozenyPermeabilitypts512-S3.5-TB-20.0-R0.005-domWidth2.5-stoppedEarlyBitDodgy'    };
% maxFrame = [12000 11800 9800 9200 7200 4800 3100 2500];

 sims = {  'CR1.179RaC80Le200KozenyPermeabilitypts256-S3.5-TB-20.0-R0.13-domWidth2.5-stoppedEarlyWavey'  };
% maxFrame = [12000 11800 9800 9200 7200 4800 3100 2500];


for i=1:length(sims)
    
    frames = 500:500:8000;

computeIceOceanSaltFlux([data_dir,sims{i},'/'], prefix, frames)

end