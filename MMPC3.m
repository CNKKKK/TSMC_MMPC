function [sys,x0,str,ts] = MMPC3(t,x,u,flag)%%逆变级双矢量，整流级零矢量选择[000000]
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
 h=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
persistent Vdc;

if isempty(Vdc)%%为第一次Vdc做准备
Vdc=180;
end

Ts=0.0001;Rf=2; Lf=0.0004; Cf=21e-6; Rs=2; Ls=0.02;%%若修改周期求取Is（K+1)的状态方程系数也需相应的修改
% % K1=0;K2=0; %%%%%K1 reactive power K2 source current
landa=0.8;k1=0.01;
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


 Io_alphar_g=(2/3)* (Io_A_g-(1/2)*Io_B_g-(1/2)*Io_C_g);      Io_beta_g=(2/3)*(sqrt(3)/2)*(Io_B_g-Io_C_g);
 Io_alphar_1=(2/3)* (Io_A_1-(1/2)*Io_B_1-(1/2)*Io_C_1);      Io_beta_1=(2/3)*(sqrt(3)/2)*(Io_B_1-Io_C_1);
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
%%%%%逆变级的遍历%%%%%%
 %%%%%%%S1 S5 S6 以及S1 S2 S6 开启以及S1 S5 S3%%%%
%     Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(2-1*0-1*0)*Vdc;%%%%逆变级S1 S5 S6
     Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(0-1);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);    
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia1=[1 0 0 0 1 1];
%      Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%逆变级S1 S2 S6
     Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1);    
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib1=[1 1 0 0 0 1];
k=1;
      Q1=(ga1*gb1)/(ga1+gb1);
%       Db1=1-landa*cos(sita)+landa/sqrt(3)*sin(sita);
%       Da1=-1+2*landa*cos(sita);
%       Dc1=1-landa*cos(sita)-landa/sqrt(3)*sin(sita);
      Da1=gb1/(ga1+gb1);
      Db1=ga1/(gb1+ga1);
      
               %%%%%%%S1 S2 S6 以及S4 S2 S6 以及S1 S5 S6开启
%     Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%逆变级S1 S2 S6
     Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[1 1 0 0 0 1];
%      Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%逆变级S4 S2 S6
     Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib2=[0 1 0 1 0 1];
      Q2=(ga2*gb2)/(ga2+gb2);
%         sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
      Da2=gb2/(gb2+ga2);
      Db2=ga2/(gb2+ga2);
    if(Q1<=Q2)
        k=1;
    else
        Ia1=Ia2;Ib1=Ib2;Da1=Da2;Db1=Db2;Q1=Q2;k=2;
    end
    
      
        %%%%%%%S4 S2 S6 以及S4 S2 S3 以及S1 S2 S6开启
%     Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%逆变级S4 S2 S6
      Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1);
      Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[0 1 0 1 0 1];
%      Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%逆变级S4 S2 S3
     Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib2=[0 1 1 1 0 0];
      Q2=(ga2*gb2)/(ga2+gb2);
%      sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
      
      Da2=gb2/(gb2+ga2);
      Db2=ga2/(gb2+ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Da1=Da2;Db1=Db2;Q1=Q2;k=3;
    end 
      

         %%%%%%%S4 S2 S3 以及S4 S5 S3 以及S4 S2 S6开启
%     Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%逆变级S4 S2 S3
     Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[0 1 1 1 0 0];
%      Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%逆变级S4 S5 S3
     Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib2=[0 0 1 1 1 0];
      Q2=(ga2*gb2)/(ga2+gb2);
%      sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
          
     Da2=gb2/(gb2+ga2);
     Db2=ga2/(gb2+ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Da1=Da2;Db1=Db2;Q1=Q2;k=4;
    end  
       
      
         %%%%%%%S4 S5 S3 以及S1 S5 S3 以及S4 S2 S3开启
%     Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%逆变级S4 S5 S3
     Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[0 0 1 1 1 0];
%      Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%逆变级S1 S5 S3
     Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib2=[1 0 1 0 1 0];
      Q2=(ga2*gb2)/(ga2+gb2);
%      sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
          
     Da2=gb2/(gb2+ga2);
     Db2=ga2/(gb2+ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Da1=Da2;Db1=Db2;Q1=Q2;k=5;
    end  
      

         %%%%%%%S1 S5 S3 以及S1 S5 S6 以及S4 S5 S3开启
%     Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%逆变级S1 S5 S3
     Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[1 0 1 0 1 0];
%      Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%逆变级S1 S5 S6
     Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(0-1);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib2=[1 0 0 0 1 1];
      Q2=(ga2*gb2)/(ga2+gb2);
%      sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
          
     Da2=gb2/(gb2+ga2);
     Db2=ga2/(gb2+ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Da1=Da2;Db1=Db2;Q1=Q2;k=6;
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
    
    
    %%Ic=[整流]*Idc;Is(k+1)=A*Is(测量）+B*Us(测量）+C*Uc(测量）+D*Ic(计算）;g=abs(Qs*-Qs)^2;da=(gb*gc)/(ga*gb+ga*gc+gb*gc);！！！！！
   

if(sector==1)%%为了保证直流环节为postive，可选择的整流级开关状态为（S2,S6),(S1,S6)以及(S1,S5),(S2,S6)；（S1S5) (S1S6)
    Va1=[0 1 0 0 0 1];Vb1=[1 0 0 0 0 1];%%整流级选择26 16
    Vdc_a=Uc_B_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_C_1;
    Ic_A_a=0*Idc;Ic_B_a=(1)*Idc;Ic_C_a=(-1)*Idc;%%%在整流级为010001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=0;Ic_C_b=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_alphar_a=(2/3)*(Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)*(Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
     Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
     Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
     Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
     Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1);    
    g1=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);    

    db=(ga*gc)/(ga*gb+ga*gc+gb*gc);da=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc=1-da-db;
    Vdc=da*Vdc_a+db*Vdc_b;
    j1=1;
    
    
    
    Va2=[0 1 0 0 0 1];Vb2=[1 0 0 0 1 0];%%整流级选择 26 15
    Vdc_b=Uc_A_1+(-1)*Uc_B_1;
    Vdc_a=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=0*Idc;Ic_B_a=(1)*Idc;Ic_C_a=(-1)*Idc;%%%在整流级为010001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为100001时的Ic
    Ic_alphar_a=(2/3)*(Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)*(Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc); 
    

    db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b;
j2=2;
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;g1=g2;j1=j2;
    end
    
    
        Va2=[1 0 0 0 0 1];Vb2=[1 0 0 0 1 0];%%整流级选择 16 15
    Vdc_a=Uc_A_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(0)*Idc;Ic_C_a=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为100010时的Ic
    Ic_alphar_a=(2/3)*(Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)*(Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
    
    
    j2=3;

        db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
        Vdc2=da2*Vdc_a+db2*Vdc_b;  
      
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;j1=j2;
    end
    
    
    
    
    
    %%sita存在已无价值，可以将2扇区合并为一个扇区
    
elseif(sector==2)%%为了保证直流环节为postive，可选择的整流级开关状态为（S2,S4),(S2,S6)以及(S2,S4),(S1,S6)AND（S2S6)(S1S6)
    Va1=[0 1 0 1 0 0];Vb1=[0 1 0 0 0 1];%%整流级选择24 26
    Vdc_a=Uc_B_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(1)*Idc;Ic_C_a=0*Idc;%%%在整流级为010100时的Ic
    Ic_A_b=0*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为010001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g1=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
j1=1;

        db=(ga*gc)/(ga*gb+ga*gc+gb*gc);da=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc=1-da-db;
        Vdc=da*Vdc_a+db*Vdc_b;

    
    Va2=[0 1 0 1 0 0];Vb2=[1 0 0 0 0 1];%%整流级选择24 16
    Vdc_a=Uc_B_1+(-1)*Uc_A_1;
    Vdc_b=Uc_A_1+(-1)*Uc_C_1;

    Ic_A_a=-1*Idc;Ic_B_a=(1)*Idc;Ic_C_a=0*Idc;%%%在整流级为010100时的Ic
    Ic_A_b=1*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
j2=2;

    db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b;

        
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;g1=g2;j1=j2;
    end
        
    
        Va2=[0 1 0 0 0 1];Vb2=[1 0 0 0 0 1];%%整流级选择26 16
    Vdc_a=Uc_B_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_C_1;
    Ic_A_a=0*Idc;Ic_B_a=(1)*Idc;Ic_C_a=-1*Idc;%%%在整流级为010001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);

    j2=3;

        db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
        Vdc2=da2*Vdc_a+db2*Vdc_b;

    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;j1=j2;
    end

    
elseif(sector==3)%%为了保证直流环节为postive，可选择的整流级开关状态为（S3,S4),(S2,S4)以及(S3,S4),(S2,S6)AND (S2S4)(S2S6)

    Va1=[0 0 1 1 0 0];Vb1=[0 1 0 1 0 0];%%整流级选择34 24
    Vdc_a=Uc_C_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_A_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(0)*Idc;Ic_C_a=1*Idc;%%%在整流级为001100时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为010100时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g1=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
j1=1;

    da=(gb*gc)/(ga*gb+ga*gc+gb*gc);db=(ga*gc)/(ga*gb+ga*gc+gb*gc);dc=1-da-db;
    Vdc=da*Vdc_a+db*Vdc_b;

    
    
    Va2=[0 0 1 1 0 0];Vb2=[0 1 0 0 0 1];%%整流级选择34 26
    Vdc_a=Uc_C_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(0)*Idc;Ic_C_a=1*Idc;%%%在整流级为001100时的Ic
    Ic_A_b=0*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为010001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);

    db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    j2=2;
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;g1=g2;j1=j2;
    end
    
        Va2=[0 1 0 1 0 0];Vb2=[0 1 0 0 0 1];%%整流级选择24 26
    Vdc_a=Uc_B_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(1)*Idc;Ic_C_a=0*Idc;%%%在整流级为010100时的Ic
    Ic_A_b=0*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(-1)*Idc;%%%在整流级为010001时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);


        db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b;

        j2=3;
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;j1=j2;
    end
    
elseif(sector==4)%%为了保证直流环节为postive，可选择的整流级开关状态为（S3,S5),(S3,S4)以及(S3,S5),(S2,S4)AND(S3S4)(S2S4)
    Va1=[0 0 1 0 1 0];Vb1=[0 0 1 1 0 0];%%整流级选择35 34
    Vdc_a=Uc_C_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_A_1;
    Ic_A_a=(0)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=1*Idc;%%%在整流级为001010时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=0;Ic_C_b=(1)*Idc;%%%在整流级为001100时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g1=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
j1=1;

     db=(ga*gc)/(ga*gb+ga*gc+gb*gc);da=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc=1-da-db;
    Vdc=da*Vdc_a+db*Vdc_b;

    
    
    
    Va2=[0 0 1 0 1 0];Vb2=[0 1 0 1 0 0];%%整流级选择35 24
    Vdc_a=Uc_C_1+(-1)*Uc_B_1;
    Vdc_b=Uc_B_1+(-1)*Uc_A_1;
    Ic_A_a=(0)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=1*Idc;%%%在整流级为0010100时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为01010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
    j2=2;

    db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;g1=g2;j1=j2;
    end
    
    Va2=[0 0 1 1 0 0];Vb2=[0 1 0 1 0 0];%%整流级选择34 24
    Vdc_a=Uc_C_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_A_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(0)*Idc;Ic_C_a=1*Idc;%%%在整流级为001100时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为01010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
j2=3;
db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b;   

    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;j1=j2;
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
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g1=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
j1=1;

    da=(gb*gc)/(ga*gb+ga*gc+gb*gc);db=(ga*gc)/(ga*gb+ga*gc+gb*gc);dc=1-da-db;
    Vdc=da*Vdc_a+db*Vdc_b;

    
    
    
    Va2=[1 0 0 0 1 0];Vb2=[0 0 1 1 0 0];%%整流级选择15 34
    Vdc_a=Uc_A_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_A_1;
    Ic_A_a=(1)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=0;%%%在整流级为010100时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(1)*Idc;%%%在整流级为001010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);

    db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b; 

    j2=2;
        if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;g1=g2;j1=j2;
        end
    
Va2=[0 0 1 0 1 0];Vb2=[0 0 1 1 0 0];%%整流级选择35 34
    Vdc_a=Uc_C_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_A_1;
    Ic_A_a=(0)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=1*Idc;%%%在整流级为001010时的Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(1)*Idc;%%%在整流级为001010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);

    j2=3;

    db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    
        if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;j1=j2;
        end
        
elseif(sector==6)%%为了保证直流环节为postive，可选择的整流级开关状态为（S1,S6),(S1,S5)以及(S1,S6),(S3,S5)and(s1s5)(s3s5)
    Va1=[1 0 0 0 0 1];Vb1=[1 0 0 0 1 0];%%整流级选择16 15
    Vdc_a=Uc_A_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(0)*Idc;Ic_C_a=(-1)*Idc;%%%在整流级为100001时的Ic
    Ic_A_b=1*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(0)*Idc;%%%在整流级为100010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g1=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
j1=1;

        db=(ga*gc)/(ga*gb+ga*gc+gb*gc);da=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc=1-da-db;
    Vdc=da*Vdc_a+db*Vdc_b;

    
    Va2=[1 0 0 0 0 1];Vb2=[0 0 1 0 1 0];%%整流级选择16 35
    Vdc_a=Uc_A_1+(-1)*Uc_C_1;
    Vdc_b=Uc_C_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(0)*Idc;Ic_C_a=-1*Idc;%%%在整流级为100001时的Ic
    Ic_A_b=(0)*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(1)*Idc;%%%在整流级为001010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);

    db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    j2=2;
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;g1=g2;j1=j2;
    end
    
        Va2=[1 0 0 0 1 0];Vb2=[0 0 1 0 1 0];%%整流级选择15 35
    Vdc_a=Uc_A_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=0*Idc;%%%在整流级为100010时的Ic
    Ic_A_b=(0)*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(1)*Idc;%%%在整流级为001010时的Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=abs(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1);
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=abs(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1);
    Is_alphar_c_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*0;
    Is_beta_c_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*0;
    gc=abs(Us_alphar*Is_beta_c_1-Us_beta*Is_alphar_c_1); 
    g2=(ga*gb*gc)/(ga*gb+ga*gc+gb*gc);
j2=3;
    

    db2=(ga*gc)/(ga*gb+ga*gc+gb*gc);da2=(gb*gc)/(ga*gb+ga*gc+gb*gc);dc2=1-da2-db2;
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;dc=dc2;Vdc=Vdc2;j1=j2;
    end
    
end
Vc1=[0 0 0 0 0 0];
%%%矢量切换点的计算%%%%%%
a=dc/4*Ts;
b=a+1/2*Ts*da*Da1;
c=b+Ts*da*Db1;
d=c+1/2*Ts*da*Da1;
e=d+1/2*dc*Ts;
f=e+1/2*Ts*db*Da1;
g=f+Ts*db*Db1;
h=g+1/2*Ts*db*Da1;
% if(rem(t,Ts)>=0&&rem(t,Ts)<Ts/2*(da*Dc1/2))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%16steps
%     h=[Va1 Ic1 Da1 e-f];
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2))
%     h=[Va1 Ia1 Da1 e-f];
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2+da*Da1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2))
%     h=[Va1 Ib1 Da1 e-f];
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2))
%     h=[Vb1 Ib1 Da1 e-f];
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2+db*Da1/2))
%     h=[Vb1 Ia1 Da1 e-f];    
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2+db*Da1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2+db*Da1/2+db*Dc1/2))
%     h=[Vb1 Ic1 Da1 e-f]; 
% elseif(rem(t,Ts)>=Ts/2&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2)))
%     h=[Vb1 Ic1 Da1 e-f];  
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2)))
%     h=[Vb1 Ia1 Da1 e-f];  
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2)))
%     h=[Vb1 Ib1 Da1 e-f]; 
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2+da*Db1/2)))
%     h=[Va1 Ib1 Da1 e-f];
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2+da*Db1/2))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2+da*Db1/2+da*Da1/2)))
%     h=[Va1 Ia1 Da1 e-f];
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2+da*Db1/2+da*Da1/2))&&rem(t,Ts)<(Ts/2+Ts/2))
%     h=[Va1 Ic1 Da1 e-f];
% end


if(rem(t,Ts)>=0&&rem(t,Ts)<a)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%6steps
    h=[Vc1 Ia1 Vng sector j1];
elseif(rem(t,Ts)>=a&&rem(t,Ts)<b)
    h=[Va1 Ia1 Vng sector j1];
elseif(rem(t,Ts)>=b&&rem(t,Ts)<c)
    h=[Va1 Ib1 Vng sector j1];
elseif(rem(t,Ts)>=c&&rem(t,Ts)<d)
    h=[Va1 Ia1 Vng sector j1];
elseif(rem(t,Ts)>=d&&rem(t,Ts)<e)
    h=[Vc1 Ia1 Vng sector j1];    
elseif(rem(t,Ts)>=e&&rem(t,Ts)<f)
    h=[Vb1 Ia1 Vng sector j1];
elseif(rem(t,Ts)>=f&&rem(t,Ts)<g)
    h=[Vb1 Ib1 Vng sector j1];
elseif(rem(t,Ts)>=g&&rem(t,Ts)<h)
    h=[Vb1 Ia1 Vng sector j1];
else
    h=[Vc1 Ia1 Vng sector j1];
end







% if(rem(t,Ts)>=0&&rem(t,Ts)<Ts*(da*Dc1/2))%%%%%%%%%%%%%%%%%%%%%%%%%%eight steps 可供逆变级零矢量，逆变级两个方向相反大作用时间相等的矢量合成
%     h=[Va1 Ic1 sector k1 j1];
% elseif(rem(t,Ts)>=Ts*(da*Dc1/2)&&rem(t,Ts)<Ts*(da*Dc1/2+da*Da1))
%     h=[Va1 Ia1 sector k1 j1];
% elseif(rem(t,Ts)>=Ts*(da*Dc1/2+da*Da1)&&rem(t,Ts)<Ts*(da*Dc1/2+da*Da1+da*Db1))
%     h=[Va1 Ib1 sector k1 j1];
% elseif(rem(t,Ts)>=Ts*(da*Dc1/2+da*Da1+da*Db1)&&rem(t,Ts)<Ts*(da*Dc1/2+da*Da1+da*Db1+da*Dc1/2))
%     h=[Va1 Id1 sector k1 j1];
% elseif(rem(t,Ts)>=Ts*(da*Dc1/2+da*Da1+da*Db1+da*Dc1/2)&&rem(t,Ts)<Ts*(da*Dc1/2+da*Da1+da*Db1+Dc1/2))
%     h=[Vb1 Id1 sector k1 j1];
% elseif(rem(t,Ts)>=Ts*(da*Dc1/2+da*Da1+da*Db1+Dc1/2)&&rem(t,Ts)<Ts*(da*Dc1/2+da*Da1+da*Db1+Dc1/2+db*Da1))
%     h=[Vb1 Ib1 sector k1 j1];   
% elseif(rem(t,Ts)>=Ts*(da*Dc1/2+da*Da1+da*Db1+Dc1/2+db*Da1)&&rem(t,Ts)<Ts*(da*Dc1/2+da*Da1+da*Db1+Dc1/2+db*Da1+db*Db1))
%     h=[Vb1 Ia1 sector k1 j1]; 
% elseif(rem(t,Ts)>=Ts*(da*Dc1/2+da*Da1+da*Db1+Dc1/2+db*Da1+db*Db1)&&rem(t,Ts)<(Ts))
%     h=[Vb1 Ic1 sector k1 j1]; 
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



