omega=pi;
kp=pi;
ks=2*pi;
y=[2 ;1];
x=[5;9];
T1= TractionGreenTensor_2D_e1(omega,kp,ks,y,x);
T2= TractionGreenTensor_2D(omega,kp,ks,y,x);
T3=TractionGreenTensor_2D_fre(omega,kp,ks,y,x);

[T3(1:4,:)-T1].'
[T3(5:8,:)-T2].'