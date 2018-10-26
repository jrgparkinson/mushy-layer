function res = getRes(folder)
res = NaN;

[tokens,matches] = regexp(folder,'-(\d+)-','tokens','match');

if length(matches) > 0
    res = str2num(tokens{1}{1});
end


end