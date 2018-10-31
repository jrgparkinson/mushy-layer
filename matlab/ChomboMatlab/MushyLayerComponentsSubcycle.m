classdef MushyLayerComponentsSubcycle
    properties
        enthalpy
        bulkConcentration
        temperature
        porosity
        liquidConcentration
        solidConcentration
        pressure
        permeability
        lambda
        temperatureAnalytic
        
        xDarcyVel
        yDarcyVel
        xUPorosity
        yUPorosity
        xUstar
        yUstar
        xadvVel
        yadvVel
        xViscousSrc
        yViscousSrc
         xDarcyVelErr
        yDarcyVelErr
        
        
    end
    methods
        function obj = MushyLayerComponentsSubcycle(comp_names)
            
            if nargin > 0
                obj = obj.parseCompList(comp_names);
            end
            
        end
        
        function obj = parseCompList(obj, comp_names)
        
            
           
                
            obj.enthalpy = find(strcmp(comp_names,'Enthalpy')); 
            obj.bulkConcentration = find(strcmp(comp_names,'Bulk concentration')); 
            obj.temperature = find(strcmp(comp_names,'Temperature')); 
            obj.porosity = find(strcmp(comp_names,'Porosity')); 
            obj.permeability = find(strcmp(comp_names,'Permeability')); 
            obj.pressure = find(strcmp(comp_names,'Pressure')); 
            obj.liquidConcentration = find(strcmp(comp_names,'Liquid concentration')); 
            obj.solidConcentration = find(strcmp(comp_names,'Solid concentration')); 
            obj.lambda = find(strcmp(comp_names,'lambda')); 
            obj.temperatureAnalytic = find(strcmp(comp_names,'T analytic')); 
            
            obj.xDarcyVel = find(strcmp(comp_names,'xDarcy velocity')); 
            obj.yDarcyVel = find(strcmp(comp_names,'yDarcy velocity')); 
            obj.xUPorosity = find(strcmp(comp_names,'xU divided by porosity')); 
            obj.yUPorosity = find(strcmp(comp_names,'yU divided by porosity')); 
             obj.xadvVel = find(strcmp(comp_names,'xAdvection velocity')); 
              obj.yadvVel = find(strcmp(comp_names,'yAdvection velocity')); 
               obj.xUstar = find(strcmp(comp_names,'xU star')); 
                obj.yUstar = find(strcmp(comp_names,'yU star'));
                 obj.xViscousSrc = find(strcmp(comp_names,'xViscous solve src')); 
                  obj.yViscousSrc = find(strcmp(comp_names,'yViscous solve src'));
                   obj.xDarcyVelErr = find(strcmp(comp_names,'xDarcy Vel Error')); 
            obj.yDarcyVelErr = find(strcmp(comp_names,'yDarcy Vel Error'));
%            obj. = find(strcmp(comp_names,'')); 

                
           
        end
    end
end