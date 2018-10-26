function plotGrids
close all;

tags = fileread('grids');

p = parseTags(tags);

tagged = zeros(32,128);

for i=1:length(p)
    coords = p(i, :);
    tagged(coords(1)+1, coords(2)+1) = 1;    
end

figure();
pcolor(tagged.');
colorbar();

end


function p = parseTags(tagsStr)

pattern = ['\s+?\d+:\s+(\((\d+),(\d+)\) \((\d+),(\d+)\) \(0,0\)\)\n'];
                [mat,tok,ext]  = regexp(tagsStr, pattern, 'match', ...
                    'tokens', 'tokenExtents');
p = [];

for i=1:length(tok)   
    temp = tok{i};          
    x = str2num(temp{1});
    y = str2num(temp{2});
    
    xe = str2num(temp{3});
    ye = str2num(temp{4});
    
    %p(end+1, :) = [x y];
    
    xi = x;
    
    while xi <= xe
       % p(end+1, :) = [xi y];
       yi = y;
     while yi <= ye
        p(end+1, :) = [xi yi];
        yi = yi + 1;
     end
        xi = xi+1;
    end
    
    
    
    
end


end