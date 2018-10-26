classdef ChComparePout
    properties
        file
        err
        
    end
    
    
    methods
        function obj = ChComparePout(loc)
            obj.file = loc;
            
            
        end
        
        function err = getErr(obj)
            err = struct();
            
            fileID = fopen(obj.file,'r');
            
            test = exist(obj.file, 'file');
            
            if exist(obj.file, 'file') == 2
                %fprintf(['Loading', file, '\n']);
                file_contents = fileread(obj.file);
            else
                return;
            end
            
            re = '([\w\s]+):([^,]*), ([^,]*), ([^,]*), ([^,\n]*)\n';
            [m, t] = regexp(file_contents,re,'tokens','match');
            
            for i=1:length(m)
                comp = m{i}{1}; 
                
                if strcmp(comp, ' step 0')
                    continue
                end
                
                if str2num(comp) == 0
                    continue
                end
                
                comp = strrep(comp, ' ', '');
                
                L1 = str2num(m{i}{2}); 
                L2 = str2num(m{i}{3}); 
                Max = str2num(m{i}{4}); 
                Sum = str2num(m{i}{5}); 
                
                err.(comp).L1 = L1;
                err.(comp).L2 = L2;
                err.(comp).Max = Max;
                err.(comp).Sum = Sum;
            end
        end
    end
    
end