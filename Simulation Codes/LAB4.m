
x=[0 0.02 0.04 0.06 0.08 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.0];
y=[0 0.134 0.23 0.304 0.365 0.418 0.579 0.665 0.729 0.779 0.825 0.87 0.915 0.958 0.979 1];


F_fit=@(a,x) (a(1).*x)./(1+a(2).*x +a(3).*x.*x);
a0=[1 1 1];
[a,resonorm]=lsqcurvefit(F_fit,a0,x,y);
F_eqll=@(x) 7.9673.*x./(1+9.4109.*x -2.4383.*(x.^2));


%OVERALL MATERIAL AND SPECIES BALANCES
%MATERIAL BALANCE D+W+50-500;
%SPECIES BALANCE 500*0.45 -0.7*50-0.97*D-0.02*W;
A=[1 1;0.97 0.02];
B=[450;500*0.45-0.7*50];
sol=linsolve(A,B);
D=sol(1);W=sol(2);
zf=0.45;q=0.8;xs=0.7;
feed_cut=fsolve(@Feedeqll,[1 1]);%intersection of eqll curve and feed line
stream_eqll=@(m) 7.9673.*m./(1+9.4109.*m -2.4383.*(m.^2));
stream_cut=(stream_eqll(0.7));%intersection of eqll curve and stream line
x2=linspace(0.45,feed_cut(1),100);%for adjusting limits of feed line 
y3=linspace(0.70,stream_cut(1),100);%for adjusting stream
x3=linspace(0.7,0.7,100);
x4=linspace(0.01,1,100);%grid for plotting D and W 
%AS SIR TOLD IN CLASS RMIN WHEN  CUTS LINE JOINING D AND STREAM LINE CUTS
%EQUILLIBRIUM CURVE SO
f3=@(k) k./(k+1) -(0.97-stream_cut)./(0.97-0.7);
R_min=fsolve(f3,1);
R_actual=2.5*R_min;
Pinch_point=[0.7,stream_cut];
xD=0.97;
xW=0.02;
%for section-1
x1_grid=linspace(0.7,0.97,101);
L1=R_actual*D;V1=(R_actual+1)*D;
m_1=L1/V1;
y_1=m_1.*(x1_grid-xD)+xD;%operating line for section 1
y1_intersect=y_1(1);
x1_intersect=x1_grid(1);
%for section-2
L2=L1-50;%STREAM WITHDRAWN
V2=V1;
m_2=L2/V2;
A1=[q/(q-1) -1;m_2 -1];
B1=[(zf/(q-1)) ; (m_2*x1_intersect-y1_intersect)];
sol_2=linsolve(A1,B1);
x2_intersect=sol_2(1);y2_intersect=sol_2(2);
x2_grid=linspace(x2_intersect,0.7,101);
y_2=m_2.*(x2_grid-x1_intersect)+y1_intersect;
%for section-3
x3_grid=linspace(xW,x2_intersect,101);
L3=L2+500*q;V3=V2-500*(1-q);
m_3=L3/V3;
y_3=m_3.*(x3_grid-xW)+xW;

%SOLVING FOR STAGES NUMBER 
figure
plot(x,y,'blue',Marker='o',LineWidth=2);
hold on
plot(x,F_eqll(x),'--r',LineWidth=2);
grid minor
legend('experimental','simulation',Location='best');
title('EQUILLIBRIUM CURVE');

figure
plot(x,F_eqll(x),'--k',LineWidth=2);
hold on 
plot(x2,q.*x2./(q-1)-zf./(q-1),linewidth=2,MarkerSize=10,Marker="o",MarkerIndices=1);
hold on
plot(x3,y3,LineWidth=2,MarkerSize=10,Marker="o",MarkerIndices=1,MarkerFaceColor='white',MarkerEdgeColor='black');
hold on
plot(xD,xD,MarkerSize=20,Marker="o",MarkerFaceColor='red',MarkerEdgeColor='green');
hold on
plot(xW,xW,MarkerSize=20,Marker="o",MarkerFaceColor='white',MarkerEdgeColor='black');
hold on
plot(x,x,'--k',LineWidth=2);
hold on
plot(Pinch_point(1),Pinch_point(2),Marker="o",MarkerSize=10,MarkerEdgeColor='black',MarkerFaceColor='red');
hold on
plot(x1_grid,y_1,'blue',LineWidth=2,MarkerSize=10,Marker="o",MarkerIndices=101,MarkerFaceColor='white',MarkerEdgeColor='black');
hold on
plot(x2_grid,y_2,'k',LineWidth=2,MarkerSize=10,Marker="o",MarkerIndices=101,MarkerFaceColor='white',MarkerEdgeColor='black');
hold on
plot(x3_grid,y_3,'red',LineWidth=2,MarkerSize=10,Marker="o",MarkerIndices=101,MarkerFaceColor='white',MarkerEdgeColor='black');
legend('fitted data','feed line','stream line','D(xD,xD)','W(xW,xW)','y=x line','Pinch Point','Op-line for sec-1','Op-line for sec-2','Op-line for sec-3',Location='best');
grid minor
axis tight
title('Equillibrium Curve');

X=xD;
Y=xD;
n_stage=-1;
%for section-1
while(X>x1_intersect)
     x_old=X;
      y_old=Y;
   
 f1=@(x) 7.9673.*x./(1+9.4109.*x -2.4383.*(x.^2))-y_old;
    k21=fsolve(f1,1);
  
    while(k21<x1_intersect)
        X_1=k21;
        Y_1=F_eqll(X_1);
        break
    end
       X=k21;
    Y=m_1.*(X-xD)+xD;
         n_stage=n_stage+1;
end

 X=X_1;Y=Y_1;
%for section-2
while(X>x2_intersect)
     x_old=X;
      y_old=Y;

 f1=@(x) 7.9673.*x./(1+9.4109.*x -2.4383.*(x.^2))-y_old;
    k21=fsolve(f1,1);
    disp(k21);
    while(k21<x2_intersect)
        X_1=k21;
        Y_1=F_eqll(X_1);

        break
    end

    X=k21;
   Y=m_2.*(X-x1_intersect)+y1_intersect;
   n_stage=n_stage+1;
end
X=X_1;Y=Y_1;
%for section-3
while(X>xW)
     x_old=X;
      y_old=Y;

 f1=@(x) 7.9673.*x./(1+9.4109.*x -2.4383.*(x.^2))-y_old;
    k21=fsolve(f1,1);
    disp(k21);
    while(k21<xW)
        X_1=k21;
        Y_1=F_eqll(X_1);
         break
    end
     X=k21;
  Y=m_3.*(X-xW)+xW;
   n_stage=n_stage+1;
end
disp(n_stage);
n_stage=n_stage+(y_old-xW)/(y_old-Y);
fprintf('Number of Stages Required is %0.2f\n',n_stage);

