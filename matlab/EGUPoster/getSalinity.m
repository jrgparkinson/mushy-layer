function S = getSalinity(dimlessSalinity)
Se = 233;
deltaS = 233-35;

S = dimlessSalinity*deltaS + Se;
end