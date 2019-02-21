function [sys,x0,str,ts] =DB(t,x,u,flag)
switch flag,
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 1,
    sys=[];
  case 2,
    sys=[];
  case 3,
    sys=mdlOutputs(t,x,u);
  case 4,
    sys=[];
  case 9,
    sys=[];
   otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
end
% mdlInitializeSizes
function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 15;
sizes.NumInputs      = 18;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [0.00001 0];     % 采样时间0.0001s,100000 Hz
end
% mdlOutputs

function sys=mdlOutputs(t,x,u)


Us_A=u(1);Us_B=u(2); Us_C=u(3);
Is_A=u(4); Is_B=u(5); Is_C=u(6);
Ui_A=u(7); Ui_B=u(8); Ui_C=u(9);%%%%%%kshik
Io_A=u(10); Io_B=u(11); Io_C=u(12);
Io_A_g=u(13); Io_B_g=u(14); Io_C_g=u(15);
Us_1_A=u(16);Us_1_B=u(17);Us_1_C=u(18);


persistent dc;
if isempty(dc)%%为第一次做准备
 dc=[140 140 0.5 0.5 0 0 0.5 0.5];%%Vdc_1 Vdc_2 da db Idc_1 Idc_2 d1 d2
end
persistent h;
if isempty(h)%%为第一次做准备
 h=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%%Vdc_1 Vdc_2 da db Idc_1 Idc_2 d1 d2
end
persistent Ua;
if isempty(Ua)%%为第一次做准备
Ua=[1 1 1 0 0 0];
end
persistent Ub;
if isempty(Ub)%%为第一次做准备
   Ub=[1 1 1 0 0 0];
end 
persistent Ia;
if isempty(Ia)%%为第一次做准备
Ia=[1 0 0 1 0 0];
end
persistent Ib;
if isempty(Ib)%%为第一次做准备
   Ib=[0 0 1 0 0 1];
end 

% persistent aa;
% if isempty(aa)%%为第一次做准备
%    aa=1;
% end 

Ts=0.0001;Fs=60;Us_g=150;Rf=3; Lf=0.0003; Cf=20e-6; Rs=2; Ls=0.03;Io_g=8;
% Ts=0.00002;Fs=50;Us_g=105;Io_g=4.5;Rf=0.5; Lf=0.0059; Cf=10e-6; Rs=10; Ls=0.015;

% % K1=0;K2=0; %%%%%K1 reactive power K2 source current
    Ui_alphar=(2/3)* (Ui_A-(1/2)*Ui_B-(1/2)*Ui_C);      Ui_beta=(2/3)*(sqrt(3)/2)*(Ui_B-Ui_C);
    Us_alphar=(2/3)* (Us_A-(1/2)*Us_B-(1/2)*Us_C);      Us_beta=(2/3)*(sqrt(3)/2)*(Us_B-Us_C);
    Is_alphar=(2/3)* (Is_A-(1/2)*Is_B-(1/2)*Is_C);      Is_beta=(2/3)*(sqrt(3)/2)*(Is_B-Is_C);
	Io_g_alphar=(2/3)* (Io_A_g-(1/2)*Io_B_g-(1/2)*Io_C_g);      Io_g_beta=(2/3)*(sqrt(3)/2)*(Io_B_g-Io_C_g);
    Io_alphar=(2/3)* (Io_A-(1/2)*Io_B-(1/2)*Io_C);      Io_beta=(2/3)*(sqrt(3)/2)*(Io_B-Io_C);
    Us_1_alphar=(2/3)* (Us_1_A-(1/2)*Us_1_B-(1/2)*Us_1_C);      Us_1_beta=(2/3)*(sqrt(3)/2)*(Us_1_B-Us_1_C);
%  sita=atan(Io_beta_g/Io_alphar_g);
%     if  Io_beta_g>0     
%         A0=1;
%     else
%         A0=0;
%     end
%     if  0.866*Io_alphar_g-0.5*Io_beta_g>0    
%         A1=1;        % sin(sita+2*pi/3)
%     else
%         A1=0;
%     end    
%     if  -0.866*Io_alphar_g-0.5*Io_beta_g>0     
%         A2=1;       % sin(sita-2*pi/3)
%     else
%         A2=0;
%     end
%     P1=4*A2+2*A1+A0;
% 
%     if  P1==1            
%         sectorn=2;
%     elseif  P1==2        
%         sectorn=6;
%     elseif  P1==3        
%         sectorn=1;
%     elseif  P1==4       
%         sectorn=4;
%     elseif  P1==5        
%         sectorn=3;
%     elseif  P1==6       
%         sectorn=5;
%     else    
%     end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Uo_A=(1/3)*(2*Ua(1)-1*Ua(2)-1*Ua(3))*dc(1)*dc(3)+(1/3)*(2*Ub(1)-1*Ub(2)-1*Ub(3))*dc(2)*dc(4);%%对Uo（k）进行矢量合成，两种开关状态对应两种Uo(k)
    Uo_B=(1/3)*(-1*Ua(1)+2*Ua(2)-1*Ua(3))*dc(1)*dc(3)+(1/3)*(-1*Ub(1)+2*Ub(2)-1*Ub(3))*dc(2)*dc(4);
    Uo_C=(1/3)*(-1*Ua(1)-1*Ua(2)+2*Ua(3))*dc(1)*dc(3)+(1/3)*(-1*Ub(1)-1*Ub(2)+2*Ub(3))*dc(2)*dc(4);
    Uo_alphar=(2/3)* (Uo_A-(1/2)*Uo_B-(1/2)*Uo_C);      Uo_beta=(2/3)*(sqrt(3)/2)*(Uo_B-Uo_C);
    Uo_1_alphar=(Ls/Ts)*(Io_g_alphar-(1-(Rs*Ts)/Ls)^2*Io_alphar-(1-(Rs*Ts)/Ls)*(Ts/Ls)*Uo_alphar);
    Uo_1_beta=(Ls/Ts)*(Io_g_beta-(1-(Rs*Ts)/Ls)^2*Io_beta-(1-(Rs*Ts)/Ls)*(Ts/Ls)*Uo_beta);
    
    
 if  Uo_1_beta>0     
        A0=1;
    else
        A0=0;
 end
    if  0.866*Uo_1_alphar-0.5*Uo_1_beta>0    
        A1=1;        % sin(sita+2*pi/3)
    else
        A1=0;
    end    
    if  -0.866*Uo_1_alphar-0.5*Uo_1_beta>0     
        A2=1;       % sin(sita-2*pi/3)
    else
        A2=0;
    end
    P1=4*A2+2*A1+A0;

    if  P1==1            
        nn=2;
    elseif  P1==2        
        nn=6;
    elseif  P1==3        
        nn=1;
    elseif  P1==4       
        nn=4;
    elseif  P1==5        
        nn=3;
    elseif  P1==6       
        nn=5;
    else    
    end
 if (nn==1)
    
   sita=atan(Uo_1_beta/Uo_1_alphar);
   d1=sin(sita);
   d2=sin(pi/3-sita);
   Uc=[0 0 0 1 1 1];
   Ua=[1 0 0 0 1 1];
   Ub=[1 1 0 0 0 1];
   Ud=[1 1 1 0 0 0];
elseif (nn==2&&Uo_1_alphar>0) 
    sita=atan(Uo_1_beta/Uo_1_alphar);
    d1=sin(-pi/3+sita);
    d2=sin(2*pi/3-sita);
   Uc=[1 1 1 0 0 0];
   Ua=[1 1 0 0 0 1];
   Ub=[0 1 0 1 0 1];
   Ud=[0 0 0 1 1 1];
 elseif (nn==2&&Uo_1_alphar<=0) 
    sita=atan(Uo_1_beta/Uo_1_alphar);
    sita=sita+pi;
    d1=sin(-pi/3+sita);
    d2=sin(2*pi/3-sita);
   Uc=[1 1 1 0 0 0];
   Ua=[1 1 0 0 0 1];
   Ub=[0 1 0 1 0 1];
   Ud=[0 0 0 1 1 1];
elseif (nn==3)
    sita=atan(Uo_1_beta/Uo_1_alphar);
    sita=sita+pi;
    d1=sin(sita-pi*2/3);
    d2=sin(pi-sita);
   Uc=[0 0 0 1 1 1]; 
   Ua=[0 1 0 1 0 1];
   Ub=[0 1 1 1 0 0];
   Ud=[1 1 1 0 0 0];
elseif (nn==4)
    sita=atan(Uo_1_beta/Uo_1_alphar);
    sita=sita+pi;
    d1=sin(sita-pi);
    d2=sin(4*pi/3-sita);
   Uc=[1 1 1 0 0 0]; 
   Ua=[0 1 1 1 0 0];
   Ub=[0 0 1 1 1 0];
   Ud=[0 0 0 1 1 1];
elseif (nn==5&&Uo_1_alphar<=0)
    sita=atan(Uo_1_beta/Uo_1_alphar);
    sita=sita+pi;
    d1=sin(sita-4*pi/3);
    d2=sin(5*pi/3-sita);
   Uc=[0 0 0 1 1 1];
   Ua=[0 0 1 1 1 0];
   Ub=[1 0 1 0 1 0];
   Ud=[1 1 1 0 0 0];
elseif (nn==5&&Uo_1_alphar>0)
    sita=atan(Uo_1_beta/Uo_1_alphar);
    d1=sin(sita-4*pi/3);
    d2=sin(5*pi/3-sita);
   Uc=[0 0 0 1 1 1]; 
   Ua=[0 0 1 1 1 0];
   Ub=[1 0 1 0 1 0];
   Ud=[1 1 1 0 0 0];
 elseif (nn==6)
    sita=atan(Uo_1_beta/Uo_1_alphar);
    d1=sin(sita-5*pi/3);
    d2=sin(2*pi-sita);
   Uc=[1 1 1 0 0 0];
   Ua=[1 0 1 0 1 0];
   Ub=[1 0 0 0 1 1];
   Ud=[0 0 0 1 1 1];
 end
    
    Idc_1=(Ua(1)*Io_A+Ua(2)*Io_B+Ua(3)*Io_C);
    Idc_2=(Ub(1)*Io_A+Ub(2)*Io_B+Ub(3)*Io_C);
    d3=1-d1-d2;
%%%%%%%%%参数计算%%%%%    
AA=[0 1/Cf;-1/Lf -Rf/Lf];
BB=[0 1/-Cf;1/Lf 0];
CC=expm(AA*Ts);
% DD=(CC-eye(2))/AA*BB;
% DD=(CC-eye(2))*BB/AA;
DD=inv(AA)*(CC-eye(2))*BB;

    a=CC(1,1);b=CC(1,2);c=CC(2,1);d=CC(2,2);%%%%%Ts=0.0001
    A=DD(1,1);B=DD(1,2);C=DD(2,1);D=DD(2,2);
     
%         a=0.994151383821788;b=0.463560667578935;c=-0.0243369350478941;d=0.945477513726000;%%%%Ts=0.00001
%     A=0.00584861617821161;B=-0.475257899935359;C=0.0243369350478941;D=0.00584861617821164;


%             a=-0.0480;b=-0.3425;c=0.0180;d=-0.0121;%%%%Ts=0.001
%     A=1.048;B=-1.7536;C=-0.0180;D=1.048;
    
   %%%计算Is*%%%%%%%%% 
    landa=0.9;
    kesi=1-8*pi^2*Fs^2*Cf*Lf;
    k1=-kesi*Us_g;
    k2=sqrt((kesi*Us_g)^2-4*Rs*Rf*Io_g/landa);
    k3=-2*Rf*kesi;
    Is_g=(k1-k2)/k3;%%输入电流的幅值应小于输出电流
%     Is_g=8;
    Is_g_A=Is_g*sin(2*pi*Fs*t);Is_g_C=Is_g*sin(2*pi*Fs*t+2/3*pi); Is_g_B=Is_g*sin(2*pi*Fs*t+4/3*pi);
    Is_g_alphar=(2/3)* (Is_g_A-(1/2)*Is_g_B-(1/2)*Is_g_C);      Is_g_beta=(2/3)*((sqrt(3)/2)*Is_g_B-(sqrt(3)/2)*Is_g_C);
 %%%%%%%%%%%%%%%%%%%%%   
    Ii_A=dc(7)*(Ia(1)-Ia(4))*Idc_1+dc(8)*(Ib(1)-Ib(4))*Idc_2;%%Ii矢量合成
    Ii_B=dc(7)*(Ia(2)-Ia(5))*Idc_1+dc(8)*(Ib(2)-Ib(5))*Idc_2;
    Ii_C=dc(7)*(Ia(3)-Ia(6))*Idc_1+dc(8)*(Ib(3)-Ib(6))*Idc_2;
    Ii_alphar=(2/3)* (Ii_A-(1/2)*Ii_B-(1/2)*Ii_C);      Ii_beta=(2/3)*(sqrt(3)/2)*(Ii_B-Ii_C);
    Is_1_A=c*Ui_A+d*Is_A+C*Us_A+D*Ii_A;
    Iinput_A=1/D*(Is_g_A-c*Ui_A-d*Is_A-C*Us_A);
    Ui_1_A=a*Ui_A+b*Is_A+A*Us_A+B*Ii_A;
    Ii_1_alphar=1/D*(Is_g_alphar-((c*a+d*c)*Ui_alphar+(c*b+d*d)*Is_alphar+(c*A+d*C)*Us_alphar+(c*B+d*D)*Ii_alphar+C*Us_1_alphar));
    Ii_1_beta=1/D*(Is_g_beta-((c*a+d*c)*Ui_beta+(c*b+d*d)*Is_beta+(c*A+d*C)*Us_beta+(c*B+d*D)*Ii_beta+C*Us_1_beta));
    Ii_1_A=Ii_1_alphar;Ii_1_B=(-1/2)*Ii_1_alphar+sqrt(3)/2*Ii_1_beta;Ii_1_C=(-1/2)*Ii_1_alphar-sqrt(3)/2*Ii_1_beta;
       

if Ii_1_A>0&&Ii_1_B<0&&Ii_1_C<=0  
   zz=1; 
   da=-Ii_1_B/Ii_1_A;
   db=-Ii_1_C/Ii_1_A;
   Ia=[1 0 0 0 0 1];
   Ib=[1 0 0 0 1 0];
elseif Ii_1_A>0&&Ii_1_B>=0&&Ii_1_C<0
    zz=2; 
    da=-Ii_1_A/Ii_1_C;
    db=-Ii_1_B/Ii_1_C;
    Ia=[0 1 0 0 0 1];
    Ib=[1 0 0 0 0 1];
elseif Ii_1_B>0&&Ii_1_A<=0&&Ii_1_C<0
    zz=3;
    da=-Ii_1_C/Ii_1_B;
    db=-Ii_1_A/Ii_1_B;
    Ia=[0 1 0 1 0 0];
    Ib=[0 1 0 0 0 1];
elseif Ii_1_A<0&&Ii_1_B>0&&Ii_1_C>=0
    zz=4;
    da=-Ii_1_B/Ii_1_A;
    db=-Ii_1_C/Ii_1_A;
    Ia=[0 0 1 1 0 0];
    Ib=[0 1 0 1 0 0];
elseif Ii_1_C>0&&Ii_1_B<=0&&Ii_1_A<0
    zz=5;
    da=-Ii_1_A/Ii_1_C;
    db=-Ii_1_B/Ii_1_C;
    Ia=[0 0 1 0 1 0];
    Ib=[0 0 1 1 0 0];   
 elseif Ii_1_C>0&&Ii_1_B<0&&Ii_1_A>=0
    zz=6;
    da=-Ii_1_C/Ii_1_B;
    db=-Ii_1_A/Ii_1_B;
    Ia=[1 0 0 0 1 0];
    Ib=[0 0 1 0 1 0];   
end
    Vdc_1=((Ia(1)-Ia(4))*Ui_A+(Ia(2)-Ia(5))*Ui_B+(Ia(3)-Ia(6))*Ui_C);
    Vdc_2=((Ib(1)-Ib(4))*Ui_A+(Ib(2)-Ib(5))*Ui_B+(Ib(3)-Ib(6))*Ui_C);

    dc=[Vdc_1 Vdc_2 da db Idc_1 Idc_2 d1 d2];
    
    
    
%     y1=d1*da*Ts;
%     y2=d2*da*Ts;
%     y3=d2*db*Ts;
%     y4=d1*db*Ts;
% if rem(t,Ts)>=0&rem(t,Ts)<(y1)
%      h=[Ia Ua 1 Is_g_B Is_g_C];
% elseif rem(t,Ts)>=(y1)&rem(t,Ts)<(y1+y2)
%      h=[Ia Ub 2 Is_g_B Is_g_C];             
% elseif rem(t,Ts)>=(y1+y2)&rem(t,Ts)<(y1+y2+y3)
%      h=[Ib Ub 3 Is_g_B Is_g_C]; 
% elseif rem(t,Ts)>=(y1+y2+y3)&rem(t,Ts)<(y1+y2+y3+y4)
%      h=[Ib Ua 4 Is_g_B Is_g_C]; 
% end  

    y1=1/4*d3*da*Ts;
    y2=d1*da*Ts;
    y3=d2*da*Ts;
    y4=1/4*d3*da*Ts;
    y5=1/4*d3*db*Ts;
    y6=d2*db*Ts;
    y7=d1*db*Ts;
    y8=1/4*d3*db*Ts;
if rem(t,Ts)>=0&rem(t,Ts)<(y1)
     h=[Ia Uc Is_1_A Ui_1_A Is_g_A];
elseif rem(t,Ts)>=(y1)&rem(t,Ts)<(y1+y2)
     h=[Ia Ua Is_1_A Ui_1_A Is_g_A];             
elseif rem(t,Ts)>=(y1+y2)&rem(t,Ts)<(y1+y2+y3)
     h=[Ia Ub Is_1_A Ui_1_A Is_g_A]; 
elseif rem(t,Ts)>=(y1+y2+y3)&rem(t,Ts)<(y1+y2+y3+y4)
     h=[Ia Ud Is_1_A Ui_1_A Is_g_A]; 
elseif rem(t,Ts)>=(y1+y2+y3+y4)&rem(t,Ts)<(y1+y2+y3+y4+y5)
     h=[Ib Ud Is_1_A Ui_1_A Is_g_A];
elseif rem(t,Ts)>=(y1+y2+y3+y4+y5)&rem(t,Ts)<(y1+y2+y3+y4+y5+y6)
     h=[Ib Ub Is_1_A Ui_1_A Is_g_A];             
elseif rem(t,Ts)>=(y1+y2+y3+y4+y5+y6)&rem(t,Ts)<(y1+y2+y3+y4+y5+y6+y7)
     h=[Ib Ua Is_1_A Ui_1_A Is_g_A]; 
elseif rem(t,Ts)>=(y1+y2+y3+y4+y5+y6+y7)&rem(t,Ts)<(y1+y2+y3+y4+y5+y6+y7+y8)
     h=[Ib Uc Is_1_A Ui_1_A Is_g_A]; 
end 

    
    sys=h;
end



