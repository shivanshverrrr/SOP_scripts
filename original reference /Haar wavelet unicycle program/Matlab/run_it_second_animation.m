tt=collocation(k) ;

figure(1)
% h = animatedline;
 % plot(xxb,xxc)
for k = 1:length(tt)
    % hold off
      scatter(x1(k),x2(k), 'o', 'MarkerFaceColor', 'r')

      hold on
circle(x_center_circle1,y_center_circle1,r_circle1,1,0,0,0.3,2)

  hold on
circle(x_center_circle2,y_center_circle2,r_circle2,0,1,0,0.3,2)

hold on
circle(x_center_circle3,y_center_circle3,r_circle3,0,0,1,0.3,2)

hold on
circle(x_center_circle4,y_center_circle4,r_circle4,1,1,0,0.3,2)

  axis([ -0.3 1.1 0 1.1])
   pause(0.1)

end




function h = circle(x,y,r,a,b,c,d,e)
th = 0:pi/50:2*pi;
xunit = r.* cos(th) + x;
yunit = r .* sin(th) + y;
h = plot(xunit, yunit,'Color',[a b c d],'LineWidth',e);
end

