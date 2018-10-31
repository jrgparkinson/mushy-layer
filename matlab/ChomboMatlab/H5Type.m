classdef H5Type
   properties
      Value
   end
   methods
      function obj = H5Type(v)
         if nargin > 0
            obj.Value = v;
         end
      end
   end
end