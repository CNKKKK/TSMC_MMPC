function [sys,x0,str,ts] = MMPC2(t,x,u,flag)%%��伶������ʸ���ϳ�
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
ts  = [0.00001 0];     % ����ʱ��0.0001s,100000 Hz
end
% mdlOutputs

function sys=mdlOutputs(t,x,u)


Us_A=u(1);Us_B=u(2); Us_C=u(3);
Is_A=u(4); Is_B=u(5); Is_C=u(6);
Uc_A_1=u(7); Uc_B_1=u(8); Uc_C_1=u(9);%%%%%%kshik
Io_A_1=u(10); Io_B_1=u(11); Io_C_1=u(12);
Io_A_g=u(13); Io_B_g=u(14); Io_C_g=u(15);
 h=[0 0 0 0 0 0 1 1 1 1 1 1 0 0 0];
persistent Vdc;

if isempty(Vdc)%%Ϊ��һ��Vdc��׼��
Vdc=120;
end

Ts=0.0001;Rf=2; Lf=0.0004; Cf=21e-6; Rs=2; Ls=0.02;%%���޸�������ȡIs��K+1)��״̬����ϵ��Ҳ����Ӧ���޸�
% % K1=0;K2=0; %%%%%K1 reactive power K2 source current
k1=0.01;
% 
% 
% 
% 
% %%%%%%%%%%%%%%����Uc Is Io
% %%%%%%ʱ�Ӳ���%%%%%%%%%%%%%%%%%%%%%��kʱ�̵Ŀ���״̬������k+1ʱ�̵�Vs,Is,Io,Uo,Ic,Uc
% % if isempty(buffer)%%�ж���һʱ�̵Ŀ���״̬�Ƿ�Ϊ�գ�Ϊ��һ�β�������ƣ�
% % a=1; b=0;c=0;d=0;e=1;f=0;
% % g=1;h=1;i=1;j=0;k=0;l=0;
% % end
% % Vdc_prev=(a-d)*Uc_A+(b-e)*Uc_B+(c-f)*Uc_C;
% % Uo_A=Vdc_prev*(g-j); Uo_B=Vdc_prev*(h-k); Uo_C=Vdc_prev*(1/3)*(i-l);%%ʹ��k-1ʱ�̵Ŀ���״̬��ȡkʱ�������ѹ
% %  Uo_A_c=Uo_A-(Uo_A+Uo_B+Uo_C)/3;Uo_B_c=Uo_B-(Uo_A+Uo_B+Uo_C)/3;Uo_C_c=Uo_C-(Uo_A+Uo_B+Uo_C)/3;%%%%%%%%%%%%%������ģ��ѹ
% %  Idc_prev=g*Io_A+h*Io_B+i*Io_C;%%ʹ��K-1ʱ�̵Ŀ���״̬��ȡkʱ�̵�����任������
% %  Ic_A=(a-d)*Idc_prev;Ic_B=(b-e)*Idc_prev;Ic_C=(c-f)*Idc_prev; 
% % 
% % 
% % Is_A_1=-0.0729*Uc_A+0.9188*Is_A+0.0729*Us_A+0.0729*Ic_A;%%%%Is(k+1)=Is(k)+Uc(k)+Us(k)+Ic(k);
% % Is_B_1=-0.0729*Uc_B+0.9188*Is_B+0.0729*Us_B+0.0739*Ic_B;
% % Is_C_1=-0.0729*Uc_C+0.9188*Is_C+0.0729*Us_C+0.0739*Ic_C;
% %  Io_A_1=1/(Ls)*((Ls-Ts*Rs)*Io_A+Ts*Uo_A_c);
% %  Io_B_1=1/(Ls)*((Ls-Ts*Rs)*Io_B+Ts*Uo_B_c);%%ʹ��k-1ʱ�̵Ŀ���״̬��ȡkʱ���������
% %  Io_C_1=1/(Ls)*((Ls-Ts*Rs)*Io_C+Ts*Uo_C_c);
% % Uc_A_1=0.9261*Uc_A+1.9431*Is_A+0.0739*Us_A-1.9505*Ic_A;
% % Uc_B_1=0.9261*Uc_B+1.9431*Is_B+0.0739*Us_B-1.9505*Ic_B;%%Io(k+1)=Io(k)+Uo(k)
% % Uc_C_1=0.9261*Uc_C+1.9431*Is_C+0.0739*Us_C-1.9505*Ic_C;
% %%%%%ʱ�Ӳ���
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
%%%%%��伶�ı���%%%%%%
 %%%%%%%S1 S5 S6 �Լ�S1 S2 S6 �����Լ�S1 S5 S3%%%%
%     Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(2-1*0-1*0)*Vdc;%%%%��伶S1 S5 S6
    Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(0-1);
    Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);    
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia1=[1 0 0 0 1 1];
%      Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
    Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1); 
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib1=[1 1 0 0 0 1];
%      Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶������
    Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ic1=[1 0 1 0 1 0];
      Q1=(ga1*gb1*gc1)/(ga1*gb1+ga1*gc1+gb1*gc1);
%       Db1=1-landa*cos(sita)+landa/sqrt(3)*sin(sita);
%       Da1=-1+2*landa*cos(sita);
%       Dc1=1-landa*cos(sita)-landa/sqrt(3)*sin(sita);
      Dc1=gb1*ga1/(ga1*gb1+ga1*gc1+gc1*gb1);
      Da1=gb1*gc1/(ga1*gb1+ga1*gc1+gc1*gb1);
      Db1=gc1*ga1/(ga1*gb1+ga1*gc1+gc1*gb1);
      k=1;
      
               %%%%%%%S1 S2 S6 �Լ�S4 S2 S6 �Լ�S1 S5 S6����
%     Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
    Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1);
    Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[1 1 0 0 0 1];
%      Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
    Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ic2=[0 1 0 1 0 1];
%      Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%��伶������
    Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(0-1);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib2=[1 0 0 0 1 1];
      Q2=(ga2*gb2*gc2)/(ga2*gb2+ga2*gc2+gb2*gc2);
%         sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
      Dc2=gb2*ga2/(ga2*gb2+ga2*gc2+gc2*gb2);
      Da2=gb2*gc2/(ga2*gb1+ga2*gc2+gc2*gb2);
      Db2=gc2*ga2/(ga2*gb2+ga2*gc2+gc2*gb2);
    if(Q1<=Q2)
        k=1;
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Da1=Da2;Db1=Db2;Dc1=Dc2;Q1=Q2;k=2;
    end
    
 
      
        %%%%%%%S4 S2 S6 �Լ�S4 S2 S3 �Լ�S1 S2 S6����
%     Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
    Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1);
    Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[0 1 0 1 0 1];
%      Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
    Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ic2=[0 1 1 1 0 0];
%      Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶������
    Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib2=[1 1 0 0 0 1];
      Q2=(ga2*gb2*gc2)/(ga2*gb2+ga2*gc2+gb2*gc2);
%      sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
      Dc2=gb2*ga2/(ga2*gb2+ga2*gc2+gc2*gb2);
      Da2=gb2*gc2/(ga2*gb1+ga2*gc2+gc2*gb2);
      Db2=gc2*ga2/(ga2*gb2+ga2*gc2+gc2*gb2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Da1=Da2;Db1=Db2;Dc1=Dc2;Q1=Q2;k=3;
    end 
    
 

         %%%%%%%S4 S2 S3 �Լ�S4 S5 S3 �Լ�S4 S2 S6����
%     Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
    Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(1-0);
    Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[0 1 1 1 0 0];
%      Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
    Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ic2=[0 0 1 1 1 0];
%      Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶������
    Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(0-1);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib2=[0 1 0 1 0 1];
      Q2=(ga2*gb2*gc2)/(ga2*gb2+ga2*gc2+gb2*gc2);
%      sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
          Dc2=ga2*gb2/(ga2*gb2+gb2*gc2+ga2*gc2);
     Da2=gc2*gb2/(ga2*gb2+gb2*gc2+ga2*gc2);
     Db2=gc2*ga2/(ga2*gb2+gb2*gc2+ga2*gc2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Da1=Da2;Db1=Db2;Dc1=Dc2;Q1=Q2;k=4;
    end  
    
            
      
         %%%%%%%S4 S5 S3 �Լ�S1 S5 S3 �Լ�S4 S2 S3����
%     Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
    Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
    Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[0 0 1 1 1 0];
%      Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
    Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ic2=[1 0 1 0 1 0];
%      Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶������
    Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(1-0);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ib2=[0 1 1 1 0 0];
      Q2=(ga2*gb2*gc2)/(ga2*gb2+ga2*gc2+gb2*gc2);
%      sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
          Dc2=ga2*gb2/(ga2*gb2+gb2*gc2+ga2*gc2);
     Da2=gc2*gb2/(ga2*gb2+gb2*gc2+ga2*gc2);
     Db2=gc2*ga2/(ga2*gb2+gb2*gc2+ga2*gc2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Da1=Da2;Db1=Db2;Dc1=Dc2;Q1=Q2;k=5;
    end  
    
      
         %%%%%%%S1 S5 S3 �Լ�S1 S5 S6 �Լ�S4 S5 S3����
%     Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
    Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
    Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ia2=[1 0 1 0 1 0];
%      Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%��伶S1 S5 S6
    Uo_A_1=Vdc*(1-0);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(0-1);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;Ic2=[1 0 0 0 1 1];
%      Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶������
    Uo_A_1=Vdc*(0-1);Uo_B_1=Vdc*(0-1);Uo_C_1=Vdc*(1-0);
     Vng=abs(Uo_A_1+Uo_B_1+Uo_C_1)/3;
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);Ib2=[0 0 1 1 1 0];
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2)+k1*Vng;
      Q2=(ga2*gb2*gc2)/(ga2*gb2+ga2*gc2+gb2*gc2);
%      sita1=pi/3-sita;
%       Db2=1-landa*cos(sita1)+landa/sqrt(3)*sin(sita1);
%       Da2=-1+2*landa*cos(sita1);
%       Dc2=1-landa*cos(sita1)-landa/sqrt(3)*sin(sita1);
          Dc2=ga2*gb2/(ga2*gb2+gb2*gc2+ga2*gc2);
     Da2=gc2*gb2/(ga2*gb2+gb2*gc2+ga2*gc2);
     Db2=gc2*ga2/(ga2*gb2+gb2*gc2+ga2*gc2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Da1=Da2;Db1=Db2;Dc1=Dc2;Q1=Q2;k=6;
    end  
    
    



Idc=Da1*(Ia1(1)*Io_A_1+Ia1(2)*Io_B_1+Ia1(3)*Io_C_1)+Db1*(Ib1(1)*Io_A_1+Ib1(2)*Io_B_1+Ib1(3)*Io_C_1)+Dc1*(Ic1(1)*Io_A_1+Ic1(2)*Io_B_1+Ic1(3)*Io_C_1);


%%������Ҳ����mpcѡ���˼�룬�������ѹ��1�����仯Ϊ��0��3/����%%

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
    
    
    %%Ic=[����]*Idc;Is(k+1)=A*Is(������+B*Us(������+C*Uc(������+D*Ic(���㣩;g=abs(Qs*-Qs)^2;da=ga/(ga+gb);����������
   

if(sector==1)%%Ϊ�˱�ֱ֤������Ϊpostive����ѡ�������������״̬Ϊ��S2,S6),(S1,S6)�Լ�(S1,S5),(S2,S6)����S1S5) (S1S6)
    Va1=[0 1 0 0 0 1];Vb1=[1 0 0 0 0 1];%%������ѡ��26 16
    Vdc_a=Uc_B_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_C_1;
    Ic_A_a=0*Idc;Ic_B_a=(1)*Idc;Ic_C_a=(-1)*Idc;%%%��������Ϊ010001ʱ��Ic
    Ic_A_b=1*Idc;Ic_B_b=0;Ic_C_b=(-1)*Idc;%%%��������Ϊ100001ʱ��Ic
    Ic_alphar_a=(2/3)*(Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)*(Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
     Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
     Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=(ga*gb)/(ga+gb);    

    db=ga/(ga+gb);da=gb/(ga+gb);
    Vdc=da*Vdc_a+db*Vdc_b;
    j1=1;
    
    
    
    Va2=[0 1 0 0 0 1];Vb2=[1 0 0 0 1 0];%%������ѡ�� 26 15
    Vdc_b=Uc_A_1+(-1)*Uc_B_1;
    Vdc_a=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=0*Idc;Ic_B_a=(1)*Idc;Ic_C_a=(-1)*Idc;%%%��������Ϊ010001ʱ��Ic
    Ic_A_b=1*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(0)*Idc;%%%��������Ϊ100001ʱ��Ic
    Ic_alphar_a=(2/3)*(Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)*(Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=(ga*gb)/(ga+gb); 
    

    db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;
j2=2;
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;j1=j2;
    end
    
    
        Va2=[1 0 0 0 0 1];Vb2=[1 0 0 0 1 0];%%������ѡ�� 16 15
    Vdc_a=Uc_A_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(0)*Idc;Ic_C_a=(-1)*Idc;%%%��������Ϊ100001ʱ��Ic
    Ic_A_b=1*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(0)*Idc;%%%��������Ϊ100010ʱ��Ic
    Ic_alphar_a=(2/3)*(Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)*(Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);
    
    
    j2=3;

        db2=ga/(ga+gb);da2=gb/(ga+gb);
        Vdc2=da2*Vdc_a+db2*Vdc_b;  
      
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;j1=j2;
    end
    
    
    
    
    
    %%sita�������޼�ֵ�����Խ�2�����ϲ�Ϊһ������
    
elseif(sector==2)%%Ϊ�˱�ֱ֤������Ϊpostive����ѡ�������������״̬Ϊ��S2,S4),(S2,S6)�Լ�(S2,S4),(S1,S6)AND��S2S6)(S1S6)
    Va1=[0 1 0 1 0 0];Vb1=[0 1 0 0 0 1];%%������ѡ��24 26
    Vdc_a=Uc_B_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(1)*Idc;Ic_C_a=0*Idc;%%%��������Ϊ010100ʱ��Ic
    Ic_A_b=0*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(-1)*Idc;%%%��������Ϊ010001ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=ga*gb/(ga+gb);
j1=1;

        db=ga/(ga+gb);da=gb/(ga+gb);
        Vdc=da*Vdc_a+db*Vdc_b;

    
    Va2=[0 1 0 1 0 0];Vb2=[1 0 0 0 0 1];%%������ѡ��24 16
    Vdc_a=Uc_B_1+(-1)*Uc_A_1;
    Vdc_b=Uc_A_1+(-1)*Uc_C_1;

    Ic_A_a=-1*Idc;Ic_B_a=(1)*Idc;Ic_C_a=0*Idc;%%%��������Ϊ010100ʱ��Ic
    Ic_A_b=1*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(-1)*Idc;%%%��������Ϊ100001ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);
j2=2;

    db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;

        
    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;j1=j2;
    end
        
    
        Va2=[0 1 0 0 0 1];Vb2=[1 0 0 0 0 1];%%������ѡ��26 16
    Vdc_a=Uc_B_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_C_1;
    Ic_A_a=0*Idc;Ic_B_a=(1)*Idc;Ic_C_a=-1*Idc;%%%��������Ϊ010001ʱ��Ic
    Ic_A_b=1*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(-1)*Idc;%%%��������Ϊ100001ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);

    j2=3;

        db2=ga/(ga+gb);da2=gb/(ga+gb);
        Vdc2=da2*Vdc_a+db2*Vdc_b;

    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;j1=j2;
    end

    
elseif(sector==3)%%Ϊ�˱�ֱ֤������Ϊpostive����ѡ�������������״̬Ϊ��S3,S4),(S2,S4)�Լ�(S3,S4),(S2,S6)AND (S2S4)(S2S6)

    Va1=[0 0 1 1 0 0];Vb1=[0 1 0 1 0 0];%%������ѡ��34 24
    Vdc_a=Uc_C_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_A_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(0)*Idc;Ic_C_a=1*Idc;%%%��������Ϊ001100ʱ��Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(0)*Idc;%%%��������Ϊ010100ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=ga*gb/(ga+gb);
j1=1;

    da=gb/(ga+gb);db=ga/(ga+gb);
    Vdc=da*Vdc_a+db*Vdc_b;

    
    
    Va2=[0 0 1 1 0 0];Vb2=[0 1 0 0 0 1];%%������ѡ��34 26
    Vdc_a=Uc_C_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(0)*Idc;Ic_C_a=1*Idc;%%%��������Ϊ001100ʱ��Ic
    Ic_A_b=0*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(-1)*Idc;%%%��������Ϊ010001ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);

    db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    j2=2;
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;j1=j2;
    end
    
        Va2=[0 1 0 1 0 0];Vb2=[0 1 0 0 0 1];%%������ѡ��24 26
    Vdc_a=Uc_B_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_C_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(1)*Idc;Ic_C_a=0*Idc;%%%��������Ϊ010100ʱ��Ic
    Ic_A_b=0*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(-1)*Idc;%%%��������Ϊ010001ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);


        db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;

        j2=3;
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;j1=j2;
    end
    
elseif(sector==4)%%Ϊ�˱�ֱ֤������Ϊpostive����ѡ�������������״̬Ϊ��S3,S5),(S3,S4)�Լ�(S3,S5),(S2,S4)AND(S3S4)(S2S4)
    Va1=[0 0 1 0 1 0];Vb1=[0 0 1 1 0 0];%%������ѡ��35 34
    Vdc_a=Uc_C_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_A_1;
    Ic_A_a=(0)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=1*Idc;%%%��������Ϊ001010ʱ��Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=0;Ic_C_b=(1)*Idc;%%%��������Ϊ001100ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=ga*gb/(ga+gb);
j1=1;

     db=ga/(ga+gb);da=gb/(ga+gb);
    Vdc=da*Vdc_a+db*Vdc_b;

    
    
    
    Va2=[0 0 1 0 1 0];Vb2=[0 1 0 1 0 0];%%������ѡ��35 24
    Vdc_a=Uc_C_1+(-1)*Uc_B_1;
    Vdc_b=Uc_B_1+(-1)*Uc_A_1;
    Ic_A_a=(0)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=1*Idc;%%%��������Ϊ0010100ʱ��Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(0)*Idc;%%%��������Ϊ01010ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);
    j2=2;

    db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;j1=j2;
    end
    
    Va2=[0 0 1 1 0 0];Vb2=[0 1 0 1 0 0];%%������ѡ��34 24
    Vdc_a=Uc_C_1+(-1)*Uc_A_1;
    Vdc_b=Uc_B_1+(-1)*Uc_A_1;
    Ic_A_a=(-1)*Idc;Ic_B_a=(0)*Idc;Ic_C_a=1*Idc;%%%��������Ϊ001100ʱ��Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(1)*Idc;Ic_C_b=(0)*Idc;%%%��������Ϊ01010ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);
j2=3;
db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;   

    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;j1=j2;
    end
    
    
    %%sita���ò��٣����Խ�5����������С��������
elseif(sector==5)%%Ϊ�˱�ֱ֤������Ϊpostive����ѡ�������������״̬Ϊ��S1,S5),(S3,S5)�Լ�(S1,S5),(S3,S4)AND(S3S5)((S3S4)

    Va1=[1 0 0 0 1 0];Vb1=[0 0 1 0 1 0];%%������ѡ��15 35
    Vdc_a=Uc_A_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_B_1;
    Ic_A_a=(1)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=0;%%%��������Ϊ100010ʱ��Ic
    Ic_A_b=(0)*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(1)*Idc;%%%��������Ϊ001010ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=ga*gb/(ga+gb);
j1=1;

    da=gb/(ga+gb);db=ga/(ga+gb);
    Vdc=da*Vdc_a+db*Vdc_b;

    
    
    
    Va2=[1 0 0 0 1 0];Vb2=[0 0 1 1 0 0];%%������ѡ��15 34
    Vdc_a=Uc_A_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_A_1;
    Ic_A_a=(1)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=0;%%%��������Ϊ010100ʱ��Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(1)*Idc;%%%��������Ϊ001010ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);

    db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b; 

    j2=2;
        if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;j1=j2;
        end
    
Va2=[0 0 1 0 1 0];Vb2=[0 0 1 1 0 0];%%������ѡ��35 34
    Vdc_a=Uc_C_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_A_1;
    Ic_A_a=(0)*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=1*Idc;%%%��������Ϊ001010ʱ��Ic
    Ic_A_b=(-1)*Idc;Ic_B_b=(0)*Idc;Ic_C_b=(1)*Idc;%%%��������Ϊ001010ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);

    j2=3;

    db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    
        if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;j1=j2;
        end
        
elseif(sector==6)%%Ϊ�˱�ֱ֤������Ϊpostive����ѡ�������������״̬Ϊ��S1,S6),(S1,S5)�Լ�(S1,S6),(S3,S5)and(s1s5)(s3s5)
    Va1=[1 0 0 0 0 1];Vb1=[1 0 0 0 1 0];%%������ѡ��16 15
    Vdc_a=Uc_A_1+(-1)*Uc_C_1;
    Vdc_b=Uc_A_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(0)*Idc;Ic_C_a=(-1)*Idc;%%%��������Ϊ100001ʱ��Ic
    Ic_A_b=1*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(0)*Idc;%%%��������Ϊ100010ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g1=ga*gb/(ga+gb);
j1=1;

        db=ga/(ga+gb);da=gb/(ga+gb);
    Vdc=da*Vdc_a+db*Vdc_b;

    
    Va2=[1 0 0 0 0 1];Vb2=[0 0 1 0 1 0];%%������ѡ��16 35
    Vdc_a=Uc_A_1+(-1)*Uc_C_1;
    Vdc_b=Uc_C_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(0)*Idc;Ic_C_a=-1*Idc;%%%��������Ϊ100001ʱ��Ic
    Ic_A_b=(0)*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(1)*Idc;%%%��������Ϊ001010ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);

    db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    j2=2;
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;g1=g2;j1=j2;
    end
    
        Va2=[1 0 0 0 1 0];Vb2=[0 0 1 0 1 0];%%������ѡ��15 35
    Vdc_a=Uc_A_1+(-1)*Uc_B_1;
    Vdc_b=Uc_C_1+(-1)*Uc_B_1;
    Ic_A_a=1*Idc;Ic_B_a=(-1)*Idc;Ic_C_a=0*Idc;%%%��������Ϊ100010ʱ��Ic
    Ic_A_b=(0)*Idc;Ic_B_b=(-1)*Idc;Ic_C_b=(1)*Idc;%%%��������Ϊ001010ʱ��Ic
    Ic_alphar_a=(2/3)* (Ic_A_a-(1/2)*Ic_B_a-(1/2)*Ic_C_a);      Ic_beta_a=(2/3)*((sqrt(3)/2)*Ic_B_a-(sqrt(3)/2)*Ic_C_a);
    Ic_alphar_b=(2/3)* (Ic_A_b-(1/2)*Ic_B_b-(1/2)*Ic_C_b);      Ic_beta_b=(2/3)*((sqrt(3)/2)*Ic_B_b-(sqrt(3)/2)*Ic_C_b);
     Is_alphar_a_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_a;
     Is_beta_a_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_a;
     ga=(Us_alphar*Is_beta_a_1-Us_beta*Is_alphar_a_1)^2;
    Is_alphar_b_1=-0.1601*Uc_alphar_1+0.2192*Is_alphar+0.1601*Us_alphar+0.4606*Ic_alphar_b;
    Is_beta_b_1=-0.1601*Uc_beta_1+0.2192*Is_beta+0.1601*Us_beta+0.4606*Ic_beta_b;
    gb=(Us_alphar*Is_beta_b_1-Us_beta*Is_alphar_b_1)^2;
    g2=ga*gb/(ga+gb);
j2=3;
    

    db2=ga/(ga+gb);da2=gb/(ga+gb);
    Vdc2=da2*Vdc_a+db2*Vdc_b;

    
    if(g1>=g2)
        Va1=Va2;Vb1=Vb2;da=da2;db=db2;Vdc=Vdc2;j1=j2;
    end
    
end
% if(rem(t,Ts)>=0&&rem(t,Ts)<Ts/2*(da*Dc1/2))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%16steps
%     h=[Va1 Ic1 Da1 Db1 Dc1];
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2))
%     h=[Va1 Ia1 Da1 Db1 Dc1];
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2+da*Da1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2))
%     h=[Va1 Ib1 Da1 Db1 Dc1];
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2))
%     h=[Vb1 Ib1 Da1 Db1 Dc1];
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2+db*Da1/2))
%     h=[Vb1 Ia1 Da1 Db1 Dc1];    
% elseif(rem(t,Ts)>=Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2+db*Da1/2)&&rem(t,Ts)<Ts/2*(da*Dc1/2+da*Da1/2+da*Db1/2+db*Db1/2+db*Da1/2+db*Dc1/2))
%     h=[Vb1 Ic1 Da1 Db1 Dc1]; 
% elseif(rem(t,Ts)>=Ts/2&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2)))
%     h=[Vb1 Ic1 Da1 Db1 Dc1];  
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2)))
%     h=[Vb1 Ia1 Da1 Db1 Dc1];  
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2)))
%     h=[Vb1 Ib1 Da1 Db1 Dc1]; 
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2+da*Db1/2)))
%     h=[Va1 Ib1 Da1 Db1 Dc1];
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2+da*Db1/2))&&rem(t,Ts)<(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2+da*Db1/2+da*Da1/2)))
%     h=[Va1 Ia1 Da1 Db1 Dc1];
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(db*Dc1/2+db*Da1/2+db*Db1/2+da*Db1/2+da*Da1/2))&&rem(t,Ts)<(Ts/2+Ts/2))
%     h=[Va1 Ic1 Da1 Db1 Dc1];
% end


if(rem(t,Ts)>=0&&rem(t,Ts)<Ts*(da*Dc1))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%6steps
    h=[Va1 Ic1 k Vdc Vng];
elseif(rem(t,Ts)>=Ts*(da*Dc1)&&rem(t,Ts)<Ts*(da*Dc1+da*Da1))
    h=[Va1 Ia1 k Vdc Vng];
elseif(rem(t,Ts)>=Ts*(da*Dc1+da*Da1)&&rem(t,Ts)<Ts*(da*Dc1+da*Da1+da*Db1))
    h=[Va1 Ib1 k Vdc Vng];
elseif(rem(t,Ts)>=Ts*(da*Dc1+da*Da1+da*Db1)&&rem(t,Ts)<Ts*(da*Dc1+da*Da1+da*Db1+db*Db1))
    h=[Vb1 Ib1 k Vdc Vng];
elseif(rem(t,Ts)>=Ts*(da*Dc1+da*Da1+da*Db1+db*Db1)&&rem(t,Ts)<Ts*(da*Dc1+da*Da1+da*Db1+db*Db1+db*Da1))
    h=[Vb1 Ia1 k Vdc Vng];    
elseif(rem(t,Ts)>=Ts*(da*Dc1+da*Da1+da*Db1+db*Db1+db*Da1)&&rem(t,Ts)<Ts*(da*Dc1+da*Da1+da*Db1+db*Db1+db*Da1+db*Dc1))
    h=[Vb1 Ic1 k Vdc Vng];
end







% if(rem(t,Ts)>=0&&rem(t,Ts)<Ts*(da*Dc1/2))%%%%%%%%%%%%%%%%%%%%%%%%%%eight steps �ɹ���伶��ʸ������伶���������෴������ʱ����ȵ�ʸ���ϳ�
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






%%Ϊ��������һ������������ʸ������ʱ����伶Ϊ100011������ʱ��Ҫ�ϸ���ѭT��a/T��b=Da/Db
%%����ʱ�����·���
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


