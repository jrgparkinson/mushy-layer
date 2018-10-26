function numCells = getNumCells(folderName)

B = regexp(folderName,'-(\d+)-','Match');

numCells = NaN;

if length(B) > 0
    
    if  ~isempty(B{1})
      num = B{1};
      num = strrep(num, '-', '');
      numCells = str2num(num);

    end
end

end