clc
clear all
format shortG

%% Que 1
%%

%given data
x = [0 0.02 0.04 0.06 0.08 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.0];
y = [0 0.134 0.23 0.304 0.365 0.418 0.579 0.665 0.729 0.779 0.825 0.87 0.915 0.958 0.979 1];

%% 1
curve_fun = @(a,x) (a(1).*x)./(1 + a(2).*x +a(3).*(x.^2));
a_ig = [0,0,0];
a = lsqcurvefit(curve_fun,a_ig,x,y);

f_curve = @(x) a(1).*x./(1 + a(2).*x + a(3).*(x.^2));
X = linspace(0,1,100);

figure
plot(X,f_curve(X), LineWidth=1.5);
hold on
plot(x,y,marker = "*",MarkerSize=4)
legend('experimental','simulation', Location='best');
title('Equilibrium Curve',FontWeight='normal');
xlabel('X'); ylabel('Y');
grid minor
axis auto
hold off

%% 3

%overall mass balance
eq1 = [1 1 ; 0.97 0.02];
eq2 = [450; 190];
sol = linsolve(eq1,eq2);
D = sol(1); W = sol(2); %kmol/hr

%% 4,5,6,7,8

%Feed
zf=0.45 ; F=500 ; q=0.8;

%Distilate
xd=0.97; xw=0.2;

%Side stream 
S=50; xs=0.70; qs=1;

%Point of contact
feed_intr = fsolve(@feed_eq, [1 1]); %poc of eqbm curve and feed line
stream_intr=(f_curve(0.7)); %poc of eqbm curve and stream line

%Feed Line
feed_x = linspace(0.45,feed_intr(1),100); %x of feed line
feed_y = feed_x.*(q./(q-1))-zf./(q-1); %y of feed line

%Side-Stream line
stream_x = linspace(0.7,0.7,100); %x of stream line
stream_y = linspace(0.70,stream_intr(1),100); %y of stream line

%Pinch point
Pinch_p=[0.7,stream_intr];

%Reflux ratio
eq3 = @(a) a./(a+1)-(0.97-stream_intr)./(0.27);
rr_min = fsolve(eq3,1);
rr_actual = 2.5*rr_min;

%actual liquid & vapor flow rates
%Section I
x4 = linspace(0.7,0.97,100);
L1 = rr_actual*D;
V1 = (rr_actual+1)*D;
m1 = L1/V1;
y4 = m1.*(x4-xd)+xd; %operating line for section 1
y1_intr = y4(1);
x1_intr = x4(1);

%Section II
L2 = L1-50;
V2 = V1;
m2 = L2/V2;
A1 = [q/(q-1) -1; m2 -1];
B1 = [(zf/(q-1)) ; (m2*x1_intr-y1_intr)];
sol_2 = linsolve(A1,B1);
x2_intr = sol_2(1);
y2_intr = sol_2(2);
x2 = linspace(x2_intr,0.7,100);
y2 = m2.*(x2-x1_intr)+y1_intr; %operating line for section 2

%Section III
x3 = linspace(xw,x2_intr,100);
L3 = L2 + 500*q;
V3 = V2 - 500*(1-q);
m3 = L3/V3;
y3 = m3.*(x3-xw)+xw; %operating line for section 3

figure
plot(X,f_curve(X),'--k',LineWidth=1.5)
hold on

plot(xd, xd, 'o', MarkerFaceColor='red', MarkerSize=6)
plot(xw, xw, 'o', MarkerFaceColor='magenta', MarkerSize=6)
plot(xs, xs, 'o', MarkerFaceColor='cyan', MarkerSize=6)
plot(zf, zf, 'o', MarkerFaceColor='blue', MarkerSize=6)

plot(X,X,'--b',linewidth=1.5)
plot(feed_x, feed_y,LineWidth=1.5) %feed line eqn
plot(stream_x,stream_y,LineWidth=1.5) %part 5 side stream line eqn

plot(Pinch_p(1), Pinch_p(2), 'o', MarkerFaceColor='green', MarkerSize=6)

plot(x4,y4,LineWidth=1.5) %op line for sec 1
plot(x2,y2,LineWidth=1.5) %op line for sec 2
plot(x3,y3,LineWidth=1.5) %op line for sec 3

%% 9

y_old = xd;
fun = @(x) (y_old - a(1)*x/(1+a(2)*x+a(3)*x^2));
x_old = fsolve(fun,xd);
plot([xd x_old],[y_old y_old],'--b',LineWidth=1)
hold on
i=1;

while (y_old >= xw)
    y_new = y_old;
    fun2 = @(x1) (-y_old+ q/(q-1)*x1-zf/(q-1));
    x4 = fsolve(fun2,0.5);
    y4 = rr_actual/(rr_actual+1)*xs + xd/(rr_actual+1);
    
    if (x_old > xs)
       y_old = rr_actual/(rr_actual+1)*x_old + xd/(rr_actual+1);
       plot([x_old x_old],[y_old y_new],'--r',LineWidth=1)
       i=i+1;
    end
    
    if (x_old <= xs && x_old > x4 && xw-L3/V3*(xw-x_old)>q/(q-1)*x_old-zf/(q-1))
        y_old = y4-L2/V2*(xs-x_old);
        plot([x_old x_old],[y_old y_new],'--r',LineWidth=1)
        i=i+1;
    end
    
    if (x_old <= x4 || xw-L3/V3*(xw-x_old) <= q/(q-1)*x_old-zf/(q-1))
        y_old = xw-L3/V3*(xw-x_old);
        plot([x_old x_old],[y_old y_new],'--r',LineWidth=1)
        i=i+1;
    end
    x_new=x_old;
    fun3 = @(x) (y_old - a(1)*x/(1+a(2)*x + a(3)*x^2));
    x_old=fsolve(fun3,x_new);
    plot([x_new x_old],[y_old y_old],'--m',LineWidth=1); 
end

xlabel('X'); ylabel('Y');
title('Equilibrium Curve of Distillation Column',FontWeight='normal')
legend('Eqbm curve','D','W','S','F','y=x line','feed line','side-stream line','Pinch Point','op line 1','op line 2', 'op line 3',location='southeast');
grid minor
hold off

% Number of plates
N = i+(y_new-xw)/(y_new-y_old)
