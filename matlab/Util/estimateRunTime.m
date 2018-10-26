function estimateRunTime(Nx, Nz, maxU, lx, hours)
dx = lx/Nx;
cfl= 0.1;
nonDimTimestep = cfl*dx/maxU


timescale = 800;

%
%
%
%
%
wallclockSecondsPerTimestep = Nx*Nz*180/(256*1024); % num seconds
secondsPerTimestep = nonDimTimestep*timescale;
wallcockSecondsPerSimSecond = wallclockSecondsPerTimestep*secondsPerTimestep


%wallClockHoursPerSimHour = wallcockSecondsPerSimSecond*hours

%realTime = 60*60; % 1 hour of ice formation





end
