close all;

x = 0:0.01:1;
y = 0:0.01:1;

u_theta = zeros(length(x), length(y));

u_x = u_theta;
u_y = u_theta;

R = 0.2;
L = 1.0;
xc = [0.5 0.5];

for i=1:length(x)
    for j=1:length(y)
        
        r = sqrt((x(i)-xc(1))^2 + (y(j)-xc(2))^2);
        theta = atan(  ((x(i)-xc(1))/(y(j)-xc(2)))  );
        
        if r < R
            u_theta(i,j) = L*( (8/3)*r^4/R^5 - 5*r^3/R^4 + (10/3)*r/R^2);
        else
            u_theta(i,j) = L/r;
        end
        
        %x_rel =x(i)-xc(1);
        %sign = -x_rel/abs(x_rel);
        
        u_x(i,j) = -sign*r*u_theta(i,j)*sin(theta);
        u_y(i,j) = sign*r*u_theta(i,j)*cos(theta);

        
    end
end

h = figure();
set(h, 'Position', [200 200 1000 1000]);
quiver(x,y, u_x, u_y);