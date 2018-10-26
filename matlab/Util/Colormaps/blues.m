function cm = blues(m)
%BLUES blues colormap

if nargin < 1
    m = 100;
end

cm = makeColorMap([8 48 107]/255, [104, 172, 212]/255, [1 1 1], m);
cm = flipud(cm);

end

