function [sys,x0,str,ts] = MMPC(t,x,u,flag)%%%整流级有零矢量1.零矢量为100100以及2.000000
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
sizes.NumInputs      = 15;
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
Uc_A_1=u(7); Uc_B_1=u(8); Uc_C_1=u(9);%%%%%%kshik
Io_A_1=u(10); Io_B_1=u(11); Io_C_1=u(12);
Io_A_g=u(13); Io_B_g=u(14); Io_C_g=u(15);
h=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
persistent Vdc;

if isempty(Vdc)%%为第一次Vdc做准备
Vdc=180;
end

Ts=0.00005;Rf=2; Lf=0.0004; Cf=21e-6; Rs=2; Ls=0.02;%%若修改周期求取Is（K+1)的状态方程系数也需相应的修改
% % K1=0;K2=0; %%%%%K1 reactive power K2 source current
% 
% 
% 
% 
% %%%%%%%%%%%%%%更新Uc Is Io
% %%%%%%时延补偿%%%%%%%%%%%%%%%%%%%%%用k时刻的开关状态估计在k+1时刻的Vs,Is,Io,Uo,Ic,Uc
% % if isempty(buffer)%%判断上一时刻的开关状态是否为空（为第一次采样而设计）
% % a=1; b=0;c=0;d=0;e=1;f=0;
% % g=1;h=1;i=1;j=0;k=0;l=0;
% % end
% % Vdc_prev=(a-d)*Uc_A+(b-e)*Uc_B+(c-f)*Uc_C;
% % Uo_A=Vdc_prev*(g-j); Uo_B=Vdc_prev*(h-k); Uo_C=Vdc_prev*(1/3)*(i-l);%%使用k-1时刻的开关状态求取k时刻输出电压
% %  Uo_A_c=Uo_A-(Uo_A+Uo_B+Uo_C)/3;Uo_B_c=Uo_B-(Uo_A+Uo_B+Uo_C)/3;Uo_C_c=Uo_C-(Uo_A+Uo_B+Uo_C)/3;%%%%%%%%%%%%%消除共模电压
% %  Idc_prev=g*Io_A+h*Io_B+i*Io_C;%%使用K-1时刻的开关状态求取k时刻的输入变换器电流
% %  Ic_A=(a-d)*Idc_prev;Ic_B=(b-e)*Idc_prev;Ic_C=(c-f)*Idc_prev; 
% % 
% % 
% % Is_A_1=-0.0729*Uc_A+0.9188*Is_A+0.0729*Us_A+0.0729*Ic_A;%%%%Is(k+1)=Is(k)+Uc(k)+Us(k)+Ic(k);
% % Is_B_1=-0.0729*Uc_B+0.9188*Is_B+0.0729*Us_B+0.0739*Ic_B;
% % Is_C_1=-0.0729*Uc_C+0.9188*Is_C+0.0729*Us_C+0.0739*Ic_C;
% %  Io_A_1=1/(Ls)*((Ls-Ts*Rs)*Io_A+Ts*Uo_A_c);
% %  Io_B_1=1/(Ls)*((Ls-Ts*Rs)*Io_B+Ts*Uo_B_c);%%使用k-1时刻的开关状态求取k时刻输出电流
% %  Io_C_1=1/(Ls)*((Ls-Ts*Rs)*Io_C+Ts*Uo_C_c);
% % Uc_A_1=0.9261*Uc_A+1.9431*Is_A+0.0739*Us_A-1.9505*Ic_A;
% % Uc_B_1=0.9261*Uc_B+1.9431*Is_B+0.0739*Us_B-1.9505*Ic_B;%%Io(k+1)=Io(k)+Uo(k)
% % Uc_C_1=0.9261*Uc_C+1.9431*Is_C+0.0739*Us_C-1.9505*Ic_C;
% %%%%%时延补偿
Uc_alphar_1=(2/3)* (Uc_A_1-(1/2)*Uc_B_1-(1/2)*Uc_C_1);      Uc_beta_1=(2/3)*(sqrt(3)/2)*(Uc_B_1-Uc_C_1);
Is_alphar=(2/3)* (Is_A-(1/2)*Is_B-(1/2)*Is_C);      Is_beta=(2/3)*(sqrt(3)/2)*(Is_B-Is_C);
sita=atan(Uc_beta_1/Uc_alphar_1);

% % Is_alphar_1=(2/3)* (Is_A_1-(1/2)*Is_B_1-(1/2)*Is_C_1);      Is_beta_1=(2/3)*((sqrt(3)/2)*Is_B_1-(sqrt(3)/2)*Is_C_1);
% % Io_g_alphar=(2/3)* (Io_A_g-(1/2)*Io_B_g-(1/2)*Io_C_g);      Io_g_beta=(2/3)*((sqrt(3)/2)*Io_B_g-(sqrt(3)/2)*Io_C_g);
% % 

% %%%%%%%%%%%%%%%%%%%%%计算给定的Is%%%%%%%%%%%%
% 
% % Pi=(1-8*pi^2*50^2*Cf*Lf)*(Is_A_1*(Us_A-Rf*Is_A_1)+Is_B_1*(Us_B-Rf*Is_A_1)+Is_C_1*(Us_C-Rf*Is_C_1));
% % Po=3*Rs*5^2;
% % landa=0.9;
% % kesi=1-8*pi^2*60^2*Cf*Lf;
% % Is_g=((-kesi*220)+sqrt((kesi*220)^2-4*kesi*Rf*Rs*(30)^2/landa))/(-2*kesi*Rf);
% % %%电流电压同相位
% % %  Is_g=2.4;
% %  Is_g_A=Is_g*sin(2*pi*60*t);
% %  Is_g_B=Is_g*sin(2*pi*60*t+2/3*pi);
% %  Is_g_C=Is_g*sin(2*pi*60*t-2/3*pi);
% %  Is_g_alphar=(2/3)* (Is_g_A-(1/2)*Is_g_B-(1/2)*Is_g_C);      Is_g_beta=(2/3)*((sqrt(3)/2)*Is_g_B-(sqrt(3)/2)*Is_g_C);
% % 
% % %  Pi_A=(1-8*pi*50^2*Cf*Lf)
% % Po_g=1.5*Rs*(5^2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %    Is_beta_2=-0.0248*Uc_beta_1+0.9817*Is_beta_1+0.0248*Us_beta+0.0059*Ic_beta_1;
% %    Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);%%%%使用k时刻的电流电压求取K+1时刻的电流
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%逆变级的遍历
     %%%%%%%S1 S5 S6 以及S1 S2 S6 开启%%%%
    Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(2-1*0-1*0)*Vdc;%%%%逆变级S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia1=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%逆变级S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib1=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%逆变级零适量
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);
      Ic1=[0 0 0 1 1 1]; Id1=[1 1 1 0 0 0];%%常规逆变级零矢量
% Ic1=[1 0 1 0 1 0];Id1=[0 1 0 1 0 1];%%逆变端由两个反矢量代替零矢量结果
      Da1=0.66*(gb1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Db1=0.66*(ga1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Dc1=1-Da1-Db1;
    Q1=(ga1*gb1*gc1)/(ga1*gb1+gb1*gc1+gc1*ga1);%%%常规逆变级带有零矢量方式的作用时间


% Da1=gb1/(ga1+gb1);Db1=ga1/(ga1+gb1);Q1=(ga1*gb1)/(ga1+gb1);k1=1;
    
         %%%%%%%S1 S2 S6 以及S4 S2 S6 开启
    Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%逆变级S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%逆变级S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%逆变级零适量
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);
      Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];%%常规逆变级零矢量
%       Ic2=[1 0 0 0 1 1];Id2=[0 1 1 1 0 0];%%逆变端由两个反矢量代替零矢量结果
      Da2=0.66*(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=0.66*(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
% Da2=gb2/(ga2+gb2);Db2=ga2/(ga2+gb2);Q2=(ga2*gb2)/(ga2+gb2);k2=2;
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;
    end
   
         %%%%%%%S4 S2 S6 以及S4 S2 S3 开启
    Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%逆变级S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%逆变级S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%逆变级零适量
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);
     Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];%%常规逆变级零矢量
%       Ic2=[1 1 0 0 0 1];Id2=[0 0 1 1 1 0];%%逆变端由两个反矢量代替零矢量结果
      Da2=0.66*(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=0.66*(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
% Da2=gb2/(ga2+gb2);Db2=ga2/(ga2+gb2);Q2=(ga2*gb2)/(ga2+gb2);k2=3;
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;
    end    
 
          
         %%%%%%%S4 S2 S3 以及S4 S5 S3 开启
    Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%逆变级S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%逆变级S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%逆变级零适量
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);
      Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];%%常规逆变级零矢量
%       Ic2=[0 1 0 1 0 1];Id2=[1 0 1 0 1 0];%%逆变端由两个反矢量代替零矢量结果
      
      
      Da2=0.66*(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=0.66*(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
% Da2=gb2/(ga2+gb2);Db2=ga2/(ga2+gb2);Q2=(ga2*gb2)/(ga2+gb2);k2=4;
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;
    end  
    
         %%%%%%%S4 S5 S3 以及S1 S5 S3 开启
    Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%逆变级S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%逆变级S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%逆变级零适量
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);
      Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];%%常规逆变级零矢量
%       Ic2=[0 1 1 1 0 0];Id2=[1 0 0 0 1 1];%%逆变端由两个反矢量代替零矢量结果
      Da2=0.66*(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=0.66*(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
%     Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
% Da2=gb2/(ga2+gb2);Db2=ga2/(ga2+gb2);Q2=(ga2*gb2)/(ga2+gb2);k2=5;
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;
    end 
    
         %%%%%%%S1 S5 S3 以及S1 S5 S6 开启
    Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%逆变级S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%逆变级S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%逆变级零适量
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);
      Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];%%常规逆变级零矢量
% Ic2=[0 0 1 1 1 0];Id2=[1 1 0 0 0 1];%%逆变端由两个反矢量代替零矢量结果

      Da2=0.66*(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=0.66*(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
% Da2=gb2/(ga2+gb2);Db2=ga2/(ga2+gb2);Q2=(ga2*gb2)/(ga2+gb2);k2=6;
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;
    end  

Idc=Da1*(Ia1(1)*Io_A_1+Ia1(2)*Io_B_1+Ia1(3)*Io_C_1)+Db1*(Ib1(1)*Io_A_1+Ib1(2)*Io_B_1+Ib1(3)*Io_C_1);


%%整流级也运用mpc选择的思想，将输入电压的1扇区变化为（0，3/Π）%%

%judge sector%
Us_alphar=(2/3)* (Us_A-(1/2)*Us_B-(1/2)*Us_C);      Us_beta=(2/3)*(sqrt(3)/2)*(Us_B-Us_C);
    if  Uc_beta_1>0     
        A0=1;
    else
        A0=0;
    end
    if  0.866*Uc_alphar_1-0.5*Uc_beta_1>0    
        A1=1;        % sin(sita+2*pi/3)
    else
        A1=0;
    end    
    if  -0.866*Uc_alphar_1-0.5*Uc_beta_1>0     
        A2=1;       % sin(sita-2*pi/3)
    else
        A2=0;
    end
    P1=4*A2+2*A1+A0;

    if  P1==1            
        sector=2;
    elseif  P1==2        
        sector=6;
    elseif  P1==3        
        sector=1;
    elseif  P1==4       
        sector=4;
    elseif  P1==5        
        sector=3;
    elseif  P1==6       
        sector=5;
    else    
    end
    
    
    %%Ic=[整流]*Idc;Is(k+1)=A*Is(测量）+B*Us(测量）+C*Uc(测量）+D*Ic(计算）;g=abs(Qs*-Qs)^2;da=ga/(ga+gb);！！！！！
   
    Is_alphar_c_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar;%%整流端寻找绝对值最小的零矢量，把此零矢量的gc也投入到调至过程中
    Is_beta_c_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta;%%无论零矢量选择Uaa，Ubb，Ucc，Ic=0*Idc;
    gc=(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1)^2;
if(sector==1)%%为了保证直流环节为postive，可选择的整流级开关状态为（S2,S6),(S1,S6)以及(S1,S5),(S2,S6)；（S1S5) (S1S6)
    Va1=[0 1 0 0 0 1];Vb1=[1 0 0 0 0 1];%%整流级选择26 16
    Vdc_a=Uc_A_1+(-1)*Uc_B_1;
    Vdc_b=Uc_A_1+(-1)*Uc_C_1;
    Ic_A_a=1*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=0*Idc;%%%在整流级为010001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=0;Ic_C_b=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
     Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
     Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);    
    
    if(abs(Uc_A_1)<=abs(Uc_B_1))
        Vc1=[1 0 0 1 0 0];%%零矢量选择14，这样开关顺序应该为26 16 14 16 26
        da=gb*gc/(ga*gc+gb*gc+ga*gb);db=ga*gc/(ga*gc+gb*gc+ga*gb);
        Vdc=da*Vdc_a+db*Vdc_b;
    elseif(abs(Uc_A_1)>abs(Uc_B_1))
        Vc1=[0 1 0 0 1 0];%%零矢量选择25，开关顺序为16 26 25 26 16
    Va1=[1 0 0 0 0 1];Vb1=[0 1 0 0 0 1];
    db=gb*gc/(ga*gc+gb*gc+ga*gb);da=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc=db*Vdc_a+da*Vdc_b;
    end
    
    
    
    Va2=[0 1 0 0 0 1];Vb2=[1 0 0 0 1 0];%%整流级选择 26 15
    Vdc_b=Uc_A_1+(-1)*Uc_B_1;
    Vdc_a=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=0*Idc;Ic_B_a=(1)*Idc;Ic_C_a=(-1)*Idc;%%%在整流级为010001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=0*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);
    
    if(abs(Uc_A_1)<=abs(Uc_B_1)&&abs(Uc_A_1)<=abs(Uc_C_1))
        Vc2=[1 0 0 1 0 0];%%A相最小时，零矢量选择14，开关顺序为26 15 14 15 26
    da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;    
    elseif(abs(Uc_B_1)<=abs(Uc_A_1)&&abs(Uc_B_1)<=abs(Uc_A_1))
        Vc2=[0 1 0 0 1 0];%%C像最小是，零矢量选择25，开关顺序为26 15 25 15 26
        da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b; 
    elseif(abs(Uc_C_1)<abs(Uc_B_1)&&abs(Uc_C_1)<abs(Uc_A_1))
        Vc2=[0  0 1 0 0 1];%%C相最小，零矢量选择36，开关顺序为15 26 36 26 15
        Va2=[1 0 0 0 1 0];Vb2=[0 1 0 0 0 1];
    db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db2*Vdc_a+da2*Vdc_b; 
    end
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;Vc1=Vc2;
    end
    
    
        Va2=[1 0 0 0 0 1];Vb2=[1 0 0 0 1 0];%%整流级选择 16 15
    Vdc_a=Uc_A_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(0)*Idc;Ic_C_a=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为100001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);
    
    
    
    if(abs(Uc_B_1)<=abs(Uc_C_1))
        Vc2=[0 1 0 0 1 0];%%零矢量选择25，开关顺序为16 15 25 15 16
        da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
        Vdc2=da2*Vdc_a+db2*Vdc_b; 
    elseif(abs(Uc_B_1)>abs(Uc_C_1))
        Vc2=[0 0 1 0 0 1];%%零矢量选择36，开关顺序为15 16 36 16 15
        Va2=[1 0 0 0 1 0];Vb2=[1 0 0 0 0 1];
        db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
        Vdc2=da2*Vdc_b+db2*Vdc_a;  %%Vdc对应原先的矢量，dadb相反
    end
      
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;Vc1=Vc2;
    end
    
    
    
    
    
    %%sita存在已无价值，可以将2扇区合并为一个扇区
    
elseif(sector==2)%%为了保证直流环节为postive，可选择的整流级开关状态为（S2,S4),(S2,S6)以及(S2,S4),(S1,S6)AND（S2S6)(S1S6)
    Va1=[0 1 0 1 0 0];Vb1=[0 1 0 0 0 1];%%整流级选择24 26
    Vdc_a=Uc_B_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=1*Idc;Ic_B_a=(1)*Idc;Ic_C_a=0*Idc;%%%在整流级为010100时的Ic
    Ic_A_b=0*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为010001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    
    if(abs(Uc_A_1)>=abs(Uc_C_1))
        Vc1=[0 0 1 0 0 1];%%零矢量选择36，开关顺序为24 26 36 26 24
        da=gb*gc/(ga*gc+gb*gc+ga*gb);db=ga*gc/(ga*gc+gb*gc+ga*gb);
        Vdc=da*Vdc_a+db*Vdc_b;
    elseif(abs(Uc_C_1)>abs(Uc_A_1))
        Vc1=[1 0 0 1 0 0];%%零矢量选择14,开关顺序为26 24 14 24 26
        Va1=[0 1 0 0 0 1];Vb1=[0 1 0 1 0 0];
        db=gb*gc/(ga*gc+gb*gc+ga*gb);da=ga*gc/(ga*gc+gb*gc+ga*gb);
        Vdc=db*Vdc_a+da*Vdc_b;
    end
    
    Va2=[0 1 0 1 0 0];Vb2=[1 0 0 0 0 1];%%整流级选择24 16
    Vdc_a=Uc_B_1+(-1)*Uc_A_1;
    Vdc_b=Uc_A_1+(-1)*Uc_C_1;

    Ic_A_a=-1*Idc;Ic_B_a=(1)*Idc;Ic_C_a=0*Idc;%%%在整流级为010100时的Ic
    Ic_A_b=1*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    if(abs(Uc_A_1)<=abs(Uc_B_1)&&abs(Uc_A_1)<=abs(Uc_C_1))
        Vc2=[1 0 0 1 0 0];%%A相最小时，零矢量选择14，开关顺序为24 16 14 16 24
    da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;    
    elseif(abs(Uc_C_1)<=abs(Uc_B_1)&&abs(Uc_C_1)<=abs(Uc_A_1))
        Vc2=[0 0 1 0 0 1];%%C像最小是，零矢量选择36，开关顺序为24 16 36 16 24
        da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b; 
    elseif(abs(Uc_B_1)<abs(Uc_C_1)&&abs(Uc_B_1)<abs(Uc_A_1))
        Vc2=[0 1 0 0 1 0];%%B相最小，零矢量选择25，开关顺序为16 24 25 24 16
        Va2=[1 0 0 0 0 1];Vb2=[0 1 0 1 0 0];
    db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db2*Vdc_a+da2*Vdc_b; 
    end
        
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;Vc1=Vc2;
    end
        
    
        Va2=[0 1 0 0 0 1];Vb2=[1 0 0 0 0 1];%%整流级选择26 16
    Vdc_a=Uc_B_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_C_1;
    Ic_A_a=0*Idc;Ic_B_a=(1)*Idc;Ic_C_a=-1*Idc;%%%在整流级为010001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    
    if(abs(Uc_A_1)<=abs(Uc_B_1))
        Vc2=[1 0 0 1 0 0];%%零矢量选择14，开关顺序为26 16 14 16 25
        da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
        Vdc2=da2*Vdc_a+db2*Vdc_b;
    elseif(abs(Uc_A_1)>abs(Uc_B_1))
        Vc2=[0 1 0 0 1 0];
        Va2=[1 0 0 0 0 1];Vb2=[0 1 0 0 0 1];%%零矢量选择25，开关顺序为16 26 25 26 16
        db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
        Vdc2=db2*Vdc_a+da2*Vdc_b;
    end
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;Vc1=Vc2;
    end

    
elseif(sector==3)%%为了保证直流环节为postive，可选择的整流级开关状态为（S3,S4),(S2,S4)以及(S3,S4),(S2,S6)AND (S2S4)(S2S6)

    Va1=[0 0 1 1 0 0];Vb1=[0 1 0 1 0 0];%%整流级选择34 24
    Vdc_a=Uc_C_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_A_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(0)*Idc;Ic_C_a=1*Idc;%%%在整流级为001100时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为010100时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    if(abs(Uc_B_1)<=abs(Uc_C_1))
        Vc1=[0 1 0 0 1 0];%%零矢量选择25，开关顺序选择34 24 25 24 34
     da=gb*gc/(ga*gc+gb*gc+ga*gb);db=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc=da*Vdc_a+db*Vdc_b;
    elseif(abs(Uc_B_1)>abs(Uc_C_1))
        Vc1=[0 0 1 0 0 1];%%零矢量选择36，开关顺序为 24 34 36 34 24
        Va1=[0 1 0 1 0 0];Vb1=[0 0 1 1 0 0];
    db=gb*gc/(ga*gc+gb*gc+ga*gb);da=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc=db*Vdc_a+da*Vdc_b;
    end
    
    
    Va2=[0 0 1 1 0 0];Vb2=[0 1 0 0 0 1];%%整流级选择34 26
    Vdc_a=Uc_C_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(0)*Idc;Ic_C_a=1*Idc;%%%在整流级为001100时的Ic
    Ic_A_b=0*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为010001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);
    if(abs(Uc_B_1)<=abs(Uc_A_1)&&abs(Uc_B_1)<=abs(Uc_C_1))
        Vc2=[0 1 0 0 1 0];%%B相最小时，零矢量选择25，开关顺序为34 26 25 26 34
    da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;    
    elseif(abs(Uc_C_1)<=abs(Uc_B_1)&&abs(Uc_C_1)<=abs(Uc_A_1))
        Vc2=[0 0 1 0 0 1];%%C像最小是，零矢量选择36，开关顺序为34 26 36 26 34
        da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b; 
    elseif(abs(Uc_A_1)<abs(Uc_C_1)&&abs(Uc_A_1)<abs(Uc_B_1))
        Vc2=[1 0 0 1 0 0];%%A相最小，零矢量选择14，开关顺序为26 34 14 34 26
        Va2=[0 1 0 0 0 1];Vb2=[0 0 1 1 0 0];
    db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db2*Vdc_a+da2*Vdc_b; 
    end
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;Vc1=Vc2;
    end
    
        Va2=[0 1 0 1 0 0];Vb2=[0 1 0 0 0 1];%%整流级选择24 26
    Vdc_a=Uc_B_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(1)*Idc;Ic_C_a=0*Idc;%%%在整流级为010100时的Ic
    Ic_A_b=0*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为010001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    if(abs(Uc_A_1)>=abs(Uc_C_1))
        Vc2=[0 0 1 0 0 1];%%零矢量选择36 ，开关顺序为 24 26 36 26 24
        da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;
    elseif(abs(Uc_A_1)<abs(Uc_C_1))
        Vc2=[1 0 0 1 0 0];%%零矢量选择14，开关顺序为26 24 14 24 26
        Vb2=[0 1 0 1 0 0];Va2=[0 1 0 0 0 1];
        db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db2*Vdc_a+da2*Vdc_b;
    end
        
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;Vc1=Vc2;
    end
    
elseif(sector==4)%%为了保证直流环节为postive，可选择的整流级开关状态为（S3,S5),(S3,S4)以及(S3,S5),(S2,S4)AND(S3S4)(S2S4)
    Va1=[0 0 1 0 1 0];Vb1=[0 0 1 1 0 0];%%整流级选择35 34
    Vdc_a=Uc_C_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_A_1;
    Ic_A_a=(0)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=1*Idc;%%%在整流级为001010时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=0;Ic_C_b=(1)*Idc;%%%在整流级为001100时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    if(abs(Uc_A_1)<=abs(Uc_B_1))
        Vc1=[1 0 0 1 0 0];%%零矢量选择14，开关顺序为 35 34 14 34 35
    da=gb*gc/(ga*gc+gb*gc+ga*gb);db=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc=da*Vdc_a+db*Vdc_b;
    elseif(abs(Uc_A_1)>abs(Uc_B_1))
        Vc1=[0 1 0 0 1 0];%%零矢量选择25，开关顺序为 34 35 25 35 34
        Va1=[0 0 1 1 0 0];Vb1=[0 0 1 0 1 0];
     db=gb*gc/(ga*gc+gb*gc+ga*gb);da=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc=db*Vdc_a+da*Vdc_b;
    end
    
    
    
    Va2=[0 0 1 0 1 0];Vb2=[0 1 0 1 0 0];%%整流级选择35 24
    Vdc_a=Uc_C_1+(-1)*Uc_B_1;
    Vdc_b=Uc_B_1+(-1)*Uc_A_1;
    Ic_A_a=(0)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=1*Idc;%%%在整流级为0010100时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为01010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);
    
    if(abs(Uc_A_1)<=abs(Uc_B_1)&&abs(Uc_A_1)<=abs(Uc_C_1))
        Vc2=[1 0 0 1 0 0];%%A相最小时，零矢量选择14，开关顺序为35 24 14 24 35
    da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;    
    elseif(abs(Uc_B_1)<=abs(Uc_A_1)&&abs(Uc_B_1)<=abs(Uc_C_1))
        Vc2=[0 1 0 0 1 0];%%B像最小是，零矢量选择25，开关顺序为35 24 25 24 35
        da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b; 
    elseif(abs(Uc_C_1)<abs(Uc_B_1)&&abs(Uc_C_1)<abs(Uc_A_1))
        Vc2=[0 0 1 0 0 1];%%c相最小，零矢量选择25，开关顺序为24 35 36 35 24
        Va2=[0 1 0 1 0 0];Vb2=[0 0 1 0 1 0];
    db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db2*Vdc_a+da2*Vdc_b; 
    end
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;Vc1=Vc2;
    end
    
    Va2=[0 0 1 1 0 0];Vb2=[0 1 0 1 0 0];%%整流级选择34 24
    Vdc_a=Uc_C_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_A_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(0)*Idc;Ic_C_a=1*Idc;%%%在整流级为001100时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为01010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    if(abs(Uc_C_1)>abs(Uc_B_1))
        Vc2=[0 1 0 0 1 0];%%零矢量选择25，开关顺序为34 24 25 24 34
    da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;   
    elseif(abs(Uc_C_1)<=abs(Uc_B_1))
        Vc2=[0 0 1 0 0 1];%%零矢量选择36，开关顺序为24 34 36 34 24
        Vb2=[0 0 1 1 0 0];Va2=[0 1 0 1 0 0];
       db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db2*Vdc_a+da2*Vdc_b;    
    end
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;Vc1=Vc2;
    end
    
    
    %%sita作用不再，可以将5扇区的两个小部分整合
elseif(sector==5)%%为了保证直流环节为postive，可选择的整流级开关状态为（S1,S5),(S3,S5)以及(S1,S5),(S3,S4)AND(S3S5)((S3S4)

    Va1=[1 0 0 0 1 0];Vb1=[0 0 1 0 1 0];%%整流级选择15 35
    Vdc_a=Uc_A_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_B_1;
    Ic_A_a=(1)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=0;%%%在整流级为100010时的Ic
    Ic_A_b=(0)*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(1)*Idc;%%%在整流级为001010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    if(abs(Uc_A_1)>=abs(Uc_C_1))
        Vc1=[0 0 1 0 0 1];%%零矢量选择36，开关顺序为15 35 36 35 15
      da=gb*gc/(ga*gc+gb*gc+ga*gb);db=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc=da*Vdc_a+db*Vdc_b;
    elseif(abs(Uc_A_1)<abs(Uc_C_1))
        Vc1=[1 0 0 1 0 0];%%零矢量选择14，开关顺序为35 15 14 15 35
    Va1=[0 0 1 0 1 0];Vb1=[1 0 0 0 1 0];
    db=gb*gc/(ga*gc+gb*gc+ga*gb);da=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc=db*Vdc_a+da*Vdc_b;
    end
    
    
    
    Va2=[1 0 0 0 1 0];Vb2=[0 0 1 1 0 0];%%整流级选择15 34
    Vdc_a=Uc_A_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_A_1;
    Ic_A_a=(1)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=0;%%%在整流级为010100时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(1)*Idc;%%%在整流级为001010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);
    if(abs(Uc_A_1)<=abs(Uc_B_1)&&abs(Uc_A_1)<=abs(Uc_C_1))
        Vc2=[1 0 0 1 0 0];%%A相最小时，零矢量选择14，开关顺序为15 34 14 34 15
    da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;    
    elseif(abs(Uc_C_1)<=abs(Uc_B_1)&&abs(Uc_C_1)<=abs(Uc_A_1))
        Vc2=[0 0 1 0 0 1];%%C像最小是，零矢量选择36，开关顺序为15 34 36 34 15
        da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b; 
    elseif(abs(Uc_B_1)<abs(Uc_C_1)&&abs(Uc_B_1)<abs(Uc_A_1))
        Vc2=[0 1 0 0 1 0];%%B相最小，零矢量选择25，开关顺序为34 15 25 15 34
        Va2=[0 0 1 1 0 0];Vb2=[1 0 0 0 1 0];
    db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db2*Vdc_a+da2*Vdc_b; 
    end
    
        if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;Vc1=Vc2;
        end
    
Va2=[0 0 1 0 1 0];Vb2=[0 0 1 1 0 0];%%整流级选择35 34
    Vdc_a=Uc_C_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_A_1;
    Ic_A_a=(0)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=1*Idc;%%%在整流级为001010时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(1)*Idc;%%%在整流级为001010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    
     if(abs(Uc_A_1)<=abs(Uc_B_1))
        Vc1=[1 0 0 1 0 0];%%零矢量选择14，开关顺序为35 34 14 34 35
    da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;
    elseif(abs(Uc_A_1)>abs(Uc_B_1))
        Vc1=[0 1 0 0 1 0];%%零矢量选择25，开关顺序为34 35 25 35 34
    Va1=[0 0 1 1 0 0];Vb1=[0 0 1 0 1 0];
    db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db*Vdc_a+da*Vdc_b;
    end
    
        if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;Vc1=Vc2;
        end
        
elseif(sector==6)%%为了保证直流环节为postive，可选择的整流级开关状态为（S1,S6),(S1,S5)以及(S1,S6),(S3,S5)and(s1s5)(s3s5)
    Va1=[1 0 0 0 0 1];Vb1=[1 0 0 0 1 0];%%整流级选择16 15
    Vdc_a=Uc_A_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(0)*Idc;Ic_C_a=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为100010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    if(abs(Uc_B_1)<=abs(Uc_C_1))
        Vc1=[0 1 0 0 1 0];%%整流级零矢量选择25，开关顺序为16 15 25 15 16
    da=gb*gc/(ga*gc+gb*gc+ga*gb);db=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc=da*Vdc_a+db*Vdc_b;
    elseif(abs(Uc_C_1)<abs(Uc_B_1))
        Vc1=[0 0 1 0 0 1];%%零矢量选择36，开关顺序为15 16 36 16 15
        Va1=[1 0 0 0 1 0];Vb1=[1 0 0 0 0 1];
        db=gb*gc/(ga*gc+gb*gc+ga*gb);da=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc=db*Vdc_a+da*Vdc_b;
    end
    
    Va2=[1 0 0 0 0 1];Vb2=[0 0 1 0 1 0];%%整流级选择16 35
    Vdc_a=Uc_A_1+(-1)*Uc_C_1;
    Vdc_b=Uc_C_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(0)*Idc;Ic_C_a=-1*Idc;%%%在整流级为100001时的Ic
    Ic_A_b=(0)*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(1)*Idc;%%%在整流级为001010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);
    if(abs(Uc_B_1)<=abs(Uc_A_1)&&abs(Uc_B_1)<=abs(Uc_C_1))
        Vc2=[0 1 0 0 1 0];%%B相最小时，零矢量选择25，开关顺序为16 35 25 35 16
    da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;    
    elseif(abs(Uc_C_1)<=abs(Uc_B_1)&&abs(Uc_C_1)<=abs(Uc_A_1))
        Vc2=[0 0 1 0 0 1];%%C像最小是，零矢量选择36，开关顺序为16 35 36 35 16
        da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b; 
    elseif(abs(Uc_A_1)<abs(Uc_C_1)&&abs(Uc_A_1)<abs(Uc_B_1))
        Vc2=[1 0 0 1 0 0];%%A相最小，零矢量选择14，开关顺序为35 16 14 16 35
        Va2=[0 0 1 0 1 0];Vb2=[1 0 0 0 0 1];
    db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db2*Vdc_a+da2*Vdc_b; 
    end
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;Vc1=Vc2;
    end
    
        Va2=[1 0 0 0 1 0];Vb2=[0 0 1 0 1 0];%%整流级选择15 35
    Vdc_a=Uc_A_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=0*Idc;%%%在整流级为100010时的Ic
    Ic_A_b=(0)*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(1)*Idc;%%%在整流级为001010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_a;
     Is_beta_a_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1052*Uc_alphar_1+0.6558*Is_alphar+0.1052*Us_alphar+0.1338*Ic_alphar_b;
    Is_beta_b_1=-0.1052*Uc_beta_1+0.6558*Is_beta+0.1052*Us_beta+0.1338*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb*gc)/(ga*gc+gb*gc+ga*gb);

    
    if(abs(Uc_C_1)<abs(Uc_A_1))
        Vc2=[0 0 1 0 0 1];%%零矢量选择36，开关顺序为15 35 36 35 15
    da2=gb*gc/(ga*gc+gb*gc+ga*gb);db2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;
    elseif(abs(Uc_C_1)>=abs(Uc_A_1))
    Vc2=[1 0 0 1 0 0];%%零矢量选择14，开关顺序为35 15 14 15 35
    Va2=[0 0 1 0 1 0];Vb2=[1 0 0 0 1 0];
    db2=gb*gc/(ga*gc+gb*gc+ga*gb);da2=ga*gc/(ga*gc+gb*gc+ga*gb);
    Vdc2=db2*Vdc_a+da2*Vdc_b;
    end
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;Vc1=Vc2;
    end
    
end
dc=1-da-db;

if(rem(t,Ts)>=0&&rem(t,Ts)<Ts/2*(da*Dc1/8))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%16steps
    h=[Va1 Ic1 1 Vdc sector];
elseif(rem(t,Ts)>=Ts/2*(da*Dc1/8)&&rem(t,Ts)<Ts/2*(da*Dc1/8+da*Da1))
    h=[Va1 Ia1 2 Vdc sector];
elseif(rem(t,Ts)>=Ts/2*(da*Dc1/8+da*Da1)&&rem(t,Ts)<Ts/2*(da*Dc1/8+da*Da1+da*Db1))
    h=[Va1 Ib1 3 Vdc sector];
elseif(rem(t,Ts)>=Ts/2*(da*Dc1/8+da*Da1+da*Db1)&&rem(t,Ts)<Ts/2*(da*Dc1/8+da*Da1+da*Db1+da*Dc1/8))
    h=[Va1 Id1 4 Vdc sector];
elseif(rem(t,Ts)>=Ts/2*(da*Dc1/8+da*Da1+da*Db1+da*Dc1/8)&&rem(t,Ts)<Ts/2*(da*Dc1/8+da*Da1+da*Db1+Dc1/8))
    h=[Vb1 Id1 5 Vdc sector];
elseif(rem(t,Ts)>=Ts/2*(da*Dc1/8+da*Da1+da*Db1+Dc1/8)&&rem(t,Ts)<Ts/2*(da*Dc1/8+da*Da1+da*Db1+Dc1/8+db*Da1))
    h=[Vb1 Ib1 6 Vdc sector];   
elseif(rem(t,Ts)>=Ts/2*(da*Dc1/8+da*Da1+da*Db1+Dc1/8+db*Da1)&&rem(t,Ts)<Ts/2*(da*Dc1/8+da*Da1+da*Db1+Dc1/8+db*Da1+db*Db1))
    h=[Vb1 Ia1 7 Vdc sector]; 
elseif(rem(t,Ts)>=Ts/2*(da*Dc1/8+da*Da1+da*Db1+Dc1/8+db*Da1+db*Db1)&&rem(t,Ts)<(Ts/2))
    h=[Vb1 Ic1 8 Vdc sector];        
elseif(rem(t,Ts)>=(Ts/2)&&rem(t,Ts)<(Ts/2+Ts/2*db*Dc1/8))
    h=[Vb1 Ic1 9 Vdc sector];   
elseif(rem(t,Ts)>=(Ts/2+Ts/2*db*Dc1/8)&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/8+db*Db1)))
    h=[Vb1 Ia1 10 Vdc sector];       
elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/8+db*Db1))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1)))
    h=[Vb1 Ib1 11 Vdc sector];  
elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1+Dc1/8)))
    h=[Vb1 Id1 12 Vdc sector]; 
elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1+db*Dc1/8))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1+Dc1/8)))
    h=[Va1 Id1 13 Vdc sector];
elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1+Dc1/8))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1+Dc1/8+da*Db1)))
    h=[Va1 Ib1 14 Vdc sector];
elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1+Dc1/8+da*Db1))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1+Dc1/8+da*Db1+da*Da1)))
    h=[Va1 Ia1 15 Vdc sector];
elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/8+db*Db1+db*Da1+Dc1/8+da*Db1+da*Da1))&&rem(t,Ts)<(Ts))
    h=[Va1 Ic1 16 Vdc sector];
end
% 

% if(rem(t,Ts)>=0&&rem(t,Ts)<Ts*(Dc1/4))%%%%%%%%%%%%%%%%%%%%%%%%%%eight steps 可供逆变级零矢量，逆变级两个方向相反大作用时间相等的矢量合成
%     h=[Va1 Ic1 1 Dc1 sector];
% elseif(rem(t,Ts)>=Ts*(Dc1/4)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1))
%     h=[Va1 Ia1 2 Dc1 sector];
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1))
%     h=[Va1 Ib1 3 Dc1 sector];
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1+Dc1/4))
%     h=[Va1 Id1 4 Dc1 sector];
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1+Dc1/4)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2))
%     h=[Vb1 Id1 5 Dc1 sector];
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2+db*Da1))
%     h=[Vb1 Ib1 6 Dc1 sector];   
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2+db*Da1)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2+db*Da1+db*Db1))
%     h=[Vb1 Ia1 7 Dc1 sector]; 
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2+db*Da1+db*Db1)&&rem(t,Ts)<(Ts))
%     h=[Vb1 Ic1 8 Dc1 sector]; 
% end






%%为了满足如一扇区整流级零矢量作用时，逆变级为100011，作用时间要严格遵循T逆a/T逆b=Da/Db
%%作用时间重新分配
% a=Ts*dc/4;b=(1/2*da*Da1)*Ts+a;c=b+da*Db1*Ts;d=c+(1/2*da*Da1)*Ts;e=d+1/2*dc*Ts;
% f=e+Ts*(1/2*db*Da1);g=f+db*Db1*Ts;h1=g+(1/2*db*Da1)*Ts;i=h1+dc/4*Ts;
% if(rem(t,Ts)>=0&&rem(t,Ts)<a)
%     h=[Vc1 Ia1 1 b-a d-c];
% elseif(rem(t,Ts)>=a&&rem(t,Ts)<b)
%     h=[Va1 Ia1 2 b-a d-c];
% elseif(rem(t,Ts)>=b&&rem(t,Ts)<c)    
%     h=[Va1 Ib1 3 b-a d-c];
% elseif(rem(t,Ts)>=c&&rem(t,Ts)<d)  
%     h=[Va1 Ia1 4 b-a d-c];
% elseif(rem(t,Ts)>=d&&rem(t,Ts)<e)  
%     h=[Vd1 Ia1 5 b-a d-c];
% elseif(rem(t,Ts)>=e&&rem(t,Ts)<f)  
%     h=[Vb1 Ia1 6 b-a d-c]; 
% elseif(rem(t,Ts)>=f&&rem(t,Ts)<g)  
%     h=[Vb1 Ib1 7 b-a d-c];
% elseif(rem(t,Ts)>=g&&rem(t,Ts)<h1)  
%     h=[Vb1 Ia1 8 b-a d-c];  
%  elseif(rem(t,Ts)>=h1&&rem(t,Ts)<i)  
%     h=[Vc1 Ia1 9 b-a d-c];
% end
sys=h;
end
