% Book: Optimal Control Applied to Biological Models
% Problem on page 44, Example 3.6

clear 
clc

% addpath('/home/amit/Dropbox/Documents/casadi-linux-matlabR2014b-v3.5.5')

addpath('/home/amit/Documents/casadi-3.6.4-linux64-matlab2018b')

import casadi.* 

opti = casadi.Opti();

k=32;
alpha=1; %alpha always equal to 1 for first order differential
b=1;


x1_0=0;
x2_0=0;



x1_1=1;
x2_1=1;



H=haar_matrix(k);
palpha = operational_matrix( k , alpha ,b, H);

x = opti.variable(1,4*k); 
C1T=x(1,1:k);           % used to express state variable x
C2T=x(1,(k+1):2*k);      %  used to express state variable control y
C3T=x(1,(2*k+1):3*k);      %  used to express control u1
C4T=x(1,(3*k+1):4*k); 



helement=@(k,t)haar_column_element(k, t);

x1dot=@(t)C1T*haar_column_vector(k, t);
x2dot=@(t)C2T*haar_column_vector(k, t);


x1state=@(t)C1T*palpha*haar_column_vector(k, t) + x1_0;
x2state=@(t)C2T*palpha*haar_column_vector(k, t) + x2_0;
v_control=@(t)C3T*haar_column_vector(k, t);
theta_control=@(t)C4T*haar_column_vector(k, t);

 % V=10000;

galerkinEq1=@(i,t) (helement(i,t)).*( x1dot(t)- (v_control(t)).*cos(theta_control(t))  );
galerkinEq2=@(i,t) (helement(i,t)).*( x2dot(t)- (v_control(t)).*sin(theta_control(t))  );

 coll=collocation(k);




A=MX.zeros(k,k);
B=MX.zeros(k,k);


E=MX.zeros(k,1);
F=MX.zeros(k,1);

    for i=1:k
        for j=1:k
        A(i,j)=  galerkinEq1(i,coll(j));
        B(i,j)=  galerkinEq2(i,coll(j));
        end
        E(i)=(1/k)*sum(A(i,1:k));
        F(i)=(1/k)*sum(B(i,1:k));
    end
         


F1=E;
F2=F;




xxa=MX.zeros(1,k);
for i=1:k
aaa=C3T*H(:,i)   ;
bbb=C4T*H(:,i)   ;
% xxa(i)=(aaa).^2  + (bbb).^2 ;
xxa(i)=(aaa).^2;
end
xxa;

cost=(1/k)*sum(xxa);

% cost=tf;

 
 x1one=C1T*palpha*H(:,k) + x1_0  -  x1_1  ;   %   equality constraint i.e. boundary condtion for x1 at x=1  
 x2one=C2T*palpha*H(:,k) + x2_0  -  x2_1  ; 


xxx=C1T*palpha*H + x1_0;ss
yyy=C2T*palpha*H + x2_0;
vvvv=C3T*H;ttv
ttheta=C4T*H;

% Four obstacles



x_center_circle1=0.3;
y_center_circle1=0.4;
r_circle1=0.11;

x_center_circle2=0.6;
y_center_circle2=0.5;
r_circle2=0.1;

x_center_circle3=0;
y_center_circle3=0.3;
r_circle3=0.16;

x_center_circle4=0.7;
y_center_circle4=0.7;
r_circle4=0.1;



c1=((xxx-x_center_circle1).^2 + (yyy-y_center_circle1).^2 - r_circle1.^2) ;
 c2=((xxx-x_center_circle2).^2 + (yyy-y_center_circle2).^2 - r_circle2.^2 );
c3=((xxx-x_center_circle3).^2 + (yyy-y_center_circle3).^2 - r_circle3.^2 );
 c4=((xxx-x_center_circle4).^2 + (yyy-y_center_circle4).^2 - r_circle4.^2 );


opti.minimize(cost);
opti.subject_to(F1==0);
opti.subject_to(F2==0);


opti.subject_to(x1one==0);
opti.subject_to(x2one==0);

opti.subject_to(c1>=0);
 opti.subject_to(c2>=0);
opti.subject_to(c3>=0);
 opti.subject_to(c4>=0);


 opti.subject_to(abs(vvvv)<=1.5);

opti.subject_to(ttheta<=pi/2);
opti.subject_to(ttheta>=-pi/2);


 opti.solver('ipopt');
 sol = opti.solve();
 sol.value(x)

c1t=sol.value(C1T);
c2t=sol.value(C2T);
c3t=sol.value(C3T);
c4t=sol.value(C4T);
 
 xxxxa=zeros(1,k);
for i=1:k
aaaaa=c3t*H(:,i)   ;
bbbbb=c4t*H(:,i)   ;
% xxxxa(i)=(aaaaa).^2  + (bbbbb).^2 ;
xxxxa(i)=(aaaaa).^2  ;
end
xxxxa;

cost_value=(1/k)*sum(xxxxa);

save('c1t_c2t_c3t_c4t_k32.mat','c1t','c2t','c3t','c4t','cost_value' )
 
 
 
 tt=collocation(k) ;
   x1=zeros(1,k);
   x2=zeros(1,k);
   u=zeros(1,k);
    for i=1:k
    x1(i)= c1t*palpha*H(:,i)  +  x1_0;
    x2(i)= c2t*palpha*H(:,i)  +  x2_0;
    vv(i) = c3t*H(:,i)  ;
    theta(i) = c4t*H(:,i)  ;
    end
    x1;
    x2;
    vv  ;
    theta;
    
    
tt=collocation(k) ;
   x1=zeros(1,k);
   x2=zeros(1,k);
   u=zeros(1,k);
    for i=1:k
    % this is for mobile robot
    x1(i)= c1t*palpha*H(:,i)  +  x1_0;
    x2(i)= c2t*palpha*H(:,i)  +  x2_0;
    vv(i) = c3t*H(:,i)  ;
    theta(i) = c4t*H(:,i)  ;
    end
    x1;
    x2;
    vv  ;
    theta;
    









































%figures after this

 figure(1)
 scatter(x1,x2,'filled')
    hold on
circle(x_center_circle1,y_center_circle1,r_circle1,1,0,0,0.3,2)

  hold on
circle(x_center_circle2,y_center_circle2,r_circle2,0,1,0,0.3,2)

hold on
circle(x_center_circle3,y_center_circle3,r_circle3,0,0,1,0.3,2)

hold on
circle(x_center_circle4,y_center_circle4,r_circle4,1,1,0,0.3,2)


figure(2)
  scatter(tt,abs(vv),'blue','filled')
  hold on
    plot(tt,abs(vv),'blue')


  hold off

   figure(3)
  scatter(tt,theta*180/pi,'blue','filled')
  hold on
    plot(tt,theta*180/pi,'blue')

    hold off




function h = circle(x,y,r,a,b,c,d,e)
th = 0:pi/50:2*pi;
xunit = r.* cos(th) + x;
yunit = r .* sin(th) + y;
h = plot(xunit, yunit,'Color',[a b c d],'LineWidth',e);
end
