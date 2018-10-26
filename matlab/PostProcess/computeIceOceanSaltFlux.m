function computeIceOceanSaltFlux(output_dir, prefix, frames)

outfile = [output_dir, 'saltFluxes.mat'];


if exist(outfile, 'file') == 2
    load(outfile)
    
else
    
    time = NaN*frames;
    F = time;
    h = time;
    Si = time;
    processedFrames = time;
    
end

for f=1:length(frames)
    
    frame = frames(f);
    fprintf('Processing frame %d \n', frame);
    
    if find(frame==processedFrames)
        fprintf('Frame already done \n');
        continue
    end
    
    dim =2;
    subcycled = true;
    
    ml = MushyLayerOutput(dim, frame, output_dir, prefix, subcycled);
    
    if length(ml.levelArray) == 0 
        continue
    end
    
    chi =ml.dataForComp(ml.components.Porosity).';
    Sl =ml.dataForComp(ml.components.Liquidconcentration).';
    V =ml.dataForComp(ml.components.yAdvectionvelocity).';
    S =ml.dataForComp(ml.components.Bulkconcentration).';
    
    [mushyi,mushyj] = find(chi<1);
    Si(f) = nanmean(nanmean(S(mushyi,mushyj)));
    
    % [~,Z] = ml.grid();
    domHeight = size(chi, 1);
    
    
    z_i = getIceBase(chi);
    
    if ~isnan(z_i)
        h(f) = (domHeight-z_i)*ml.levelArray(1).dx;
        F(f) = computeFlux(Sl, V, max(z_i-10, 1));
    end
    
    processedFrames(f) = frame;
    
    
    time(f) = ml.t;
    
end

figure();
m = 2; n=2;
subplot(m, n, 1);
velScale = 1.25e-6;
plot(time, -F*velScale*100, '-');

ylabel('$F (cm/s)$');
xlabel('$t$');

subplot(m, n, 2);
plot(time, h, '-');
ylabel('$h (cm)$');
xlabel('$t$');

subplot(m, n, 3);
Se = 233;
deltaS = 233-35;
plot(time, Se+deltaS*Si, '-');
ylabel('$S_i (g/kg)$');
xlabel('$t$');

save(outfile, 'time', 'F', 'Si', 'h', 'processedFrames');


end

function F = computeFlux(Sl, V, z_i)

flux = Sl.*V;
F = sum(squeeze(flux(z_i, :)));
%TODO

end

function base_i = getIceBase(chi)
[mush_i,~] = find(chi<1);
base_i = min(mush_i);

if length(base_i) < 1
    base_i = NaN;
end

%z = NaN;
%TODO
end




