function [x,y]=myrectangle(t,type)
%rotate
 a=sin(pi/4);
 radius = 0.5;
 if type==1
     x1 = radius*( cos(t).^3 + cos(t) );   %% primal function
     y1 = radius*( sin(t).^3 + sin(t) );
     x = a*x1-a*y1;
     y = a*x1+a*y1+10;
 else if type==2
        x1 = -radius*sin(t).*(  3*cos(t).^2 + 1  );  %% derivative of order one
        y1 = radius*cos(t).*(   3*sin(t).^2 + 1 );
        x = a*x1-a*y1;
        y = a*x1+a*y1;
     else if type==3     %% derivative of order two
             x1 = -radius*( 3*cos(t).^3 + cos(t) -6*cos(t).*sin(t).^2 );
             y1 =  radius*( -(3*sin(t).^3 + sin(t)) + 6*sin(t).*cos(t).^2 ) ;
             x = a*x1-a*y1;
             y = a*x1+a*y1;
         end
     end
 end

 return