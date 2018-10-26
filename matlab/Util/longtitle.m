function longtitle(title_str, lineWidth)
%LONGTITLE Display long plot titles
%   Break long plot titles over multiple lines
if nargin < 2
    lineWidth = 15;
end


if length(title_str) <= lineWidth
    title(title_str);
    return;
end

% Try and break title sensibly
% Search backwards from max lineWidth character and find somewhere to break
% title

new_title = {};

remaining_title = title_str;

while length(remaining_title) > lineWidth
   [split, remaining_title] = splitTitle(remaining_title, lineWidth);
   new_title{end+1} = split;
end

new_title{end+1} = remaining_title;
    
title(new_title);

end

function [split, remaining_title] = splitTitle(remaining_title, lineWidth)

for i = lineWidth:-1:1
    splitHere = false;
    
   if remaining_title(i) == '-'
       splitHere = true;
   end
   
   % other possibilities...
   %fprintf('Consider %s%s \n',  remaining_title(i),remaining_title(i+1));
   if isnumber(remaining_title(i)) && ~isnumber(remaining_title(i+1))
       %fprintf('Split here \n'); 
       splitHere = true;
   end
   
   
   if splitHere
       [split, remaining_title] = splitAt(remaining_title, i);
       return;
   end
   
end

% If we make it this far, just split at lineWidth
[split, remaining_title] = splitAt(remaining_title, lineWidth);


end

function [split, remaining_title] = splitAt(remaining_title, i)
split = remaining_title(1:i);
remaining_title = remaining_title(i+1:end);
end


function isnum = isnumber(test_str)
if num2str(str2num(test_str)) == test_str
    isnum = true;
else
    isnum = false;
end

end