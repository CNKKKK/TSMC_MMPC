function [sys,x0,str,ts] = MMPC(t,x,u,flag)
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
Is_A_1=u(4); Is_B_1=u(5); Is_C_1=u(6);
Uc_A_1=u(7); Uc_B_1=u(8); Uc_C_1=u(9);%%%%%%kshik
Io_A_1=u(10); Io_B_1=u(11); Io_C_1=u(12);
Io_A_g=u(13); Io_B_g=u(14); Io_C_g=u(15);
h=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];

Ts=0.00005;Rf=0.5; Lf=0.0004; Cf=21e-6; Rs=2; Ls=0.02;%%���޸�������ȡIs��K+1)��״̬����ϵ��Ҳ����Ӧ���޸�
% % K1=0;K2=0; %%%%%K1 reactive power K2 source current
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

sita=atan(Uc_beta_1/Uc_alphar_1);

% % Is_alphar_1=(2/3)* (Is_A_1-(1/2)*Is_B_1-(1/2)*Is_C_1);      Is_beta_1=(2/3)*((sqrt(3)/2)*Is_B_1-(sqrt(3)/2)*Is_C_1);
% % Io_g_alphar=(2/3)* (Io_A_g-(1/2)*Io_B_g-(1/2)*Io_C_g);      Io_g_beta=(2/3)*((sqrt(3)/2)*Io_B_g-(sqrt(3)/2)*Io_C_g);
% % 
% % %judge sector%
% % Us_alphar=(2/3)* (Us_A-(1/2)*Us_B-(1/2)*Us_C);      Us_beta=(2/3)*(sqrt(3)/2)*(Us_B-Us_C);
% %     if  Us_beta>0     
% %         A0=1;
% %     else
% %         A0=0;
% %     end
% %     if  0.866*Us_alphar-0.5*Us_beta>0    
% %         A1=1;        % sin(sita+2*pi/3)
% %     else
% %         A1=0;
% %     end    
% %     if  -0.866*Us_alphar-0.5*Us_beta>0     
% %         A2=1;       % sin(sita-2*pi/3)
% %     else
% %         A2=0;
% %     end
% %     P1=4*A2+2*A1+A0;
% % 
% %     if  P1==1            
% %         sector=2;
% %     elseif  P1==2        
% %         sector=6;
% %     elseif  P1==3        
% %         sector=1;
% %     elseif  P1==4       
% %         sector=4;
% %     elseif  P1==5        
% %         sector=3;
% %     elseif  P1==6       
% %         sector=5;
% %     else    
% %     end
% %%%%%%%%%%%%%%%%%********%%%%%%%%%%%end
% %%%%%%%%%%%%%%%%%%%%%���������Is%%%%%%%%%%%%
% 
% % Pi=(1-8*pi^2*50^2*Cf*Lf)*(Is_A_1*(Us_A-Rf*Is_A_1)+Is_B_1*(Us_B-Rf*Is_A_1)+Is_C_1*(Us_C-Rf*Is_C_1));
% % Po=3*Rs*5^2;
% % landa=0.9;
% % kesi=1-8*pi^2*60^2*Cf*Lf;
% % Is_g=((-kesi*220)+sqrt((kesi*220)^2-4*kesi*Rf*Rs*(30)^2/landa))/(-2*kesi*Rf);
% % %%������ѹͬ��λ
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
% %    Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1_c);%%%%ʹ��kʱ�̵ĵ�����ѹ��ȡK+1ʱ�̵ĵ���
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Uc_A_1>0&&Uc_B_1<0&&Uc_C_1<=0 )

    %������ѡ�� S1 S5 �Լ� S1 S6
    Vdc_a=1*Uc_A_1+(-1)*Uc_B_1;%%S1 S5��״̬�µ�Vdc
    Vdc_b=1*Uc_A_1+(-1)*Uc_C_1;%%S1 S6��״̬�µ�Vdc
    da=(sin(sita-(-pi/6)));db=(sin(pi/3-(sita-(-pi/6)))); dc=1-da-db;
%     da=(-1)*(Uc_B_1/Uc_A_1);db=(-1)*(Uc_C_1/Uc_A_1);
     Vdc=da*Vdc_a+db*Vdc_b;
    Va=[1 0 0 0 1 0];Vb=[1 0 0 0 0 1];
    if(abs(Uc_B_1)<=abs(Uc_C_1)) Vc=[0 1 0 0 1 0];
    else 
        Vc=[0 0 1 0 0 1];
    end
        k=1;
     %%%%%%%S1 S5 S6 �Լ�S1 S2 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(2-1*0-1*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia1=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib1=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic1=[0 0 0 1 1 1]; Id1=[1 1 1 0 0 0];
      Da1=(gb1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Db1=(ga1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Dc1=1-Da1-Db1;
    Q1=(ga1*gb1*gc1)/(ga1*gb1+gb1*gc1+gc1*ga1);
    
         %%%%%%%S1 S2 S6 �Լ�S4 S2 S6 ����
    Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end
   
         %%%%%%%S4 S2 S6 �Լ�S4 S2 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end    
 
          
         %%%%%%%S4 S2 S3 �Լ�S4 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end  
    
         %%%%%%%S4 S5 S3 �Լ�S1 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end 
    
         %%%%%%%S1 S5 S3 �Լ�S1 S5 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end  
elseif Uc_A_1>0&&Uc_B_1>=0&&Uc_C_1<0

    %������ѡ�� S2 S6 �Լ� S1 S6
    Vdc_a=1*Uc_A_1+(-1)*Uc_C_1;%%S1 S6��״̬�µ�Vdc
    Vdc_b=1*Uc_B_1+(-1)*Uc_C_1;%%S2 S6��״̬�µ�Vdc
    da=(sin(sita-pi/6));db=(sin(pi/2-(sita)));dc=1-da-db;
%     da=(-1)*(Uc_A_1/Uc_C_1);db=(-1)*(Uc_B_1/Uc_C_1);
     Vdc=da*Vdc_a+db*Vdc_b;
    Va=[1 0 0 0 0 1];Vb=[0 1 0 0 0 1];
    if(abs(Uc_A_1)<=abs(Uc_B_1))
        Vc=[1 0 0 1 0 0];
    else Vc=[0 1 0 0 1 0];
    end
    k=2;
     %%%%%%%S1 S5 S6 �Լ�S1 S2 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(2-1*0-1*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia1=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib1=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic1=[0 0 0 1 1 1]; Id1=[1 1 1 0 0 0];
      Da1=(gb1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Db1=(ga1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Dc1=1-Da1-Db1;
    Q1=(ga1*gb1*gc1)/(ga1*gb1+gb1*gc1+gc1*ga1);
    
         %%%%%%%S1 S2 S6 �Լ�S4 S2 S6 ����
    Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end
   
         %%%%%%%S4 S2 S6 �Լ�S4 S2 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end    
 
          
         %%%%%%%S4 S2 S3 �Լ�S4 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end  
    
         %%%%%%%S4 S5 S3 �Լ�S1 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end 
    
         %%%%%%%S1 S5 S3 �Լ�S1 S5 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end
           
elseif Uc_B_1>0&&Uc_A_1<=0&&Uc_C_1<0  
  
    %������ѡ�� S2 S6 �Լ� S2 S4
    Vdc_a=1*Uc_B_1+(-1)*Uc_C_1;%%S2 S6��״̬�µ�Vdc
    Vdc_b=1*Uc_B_1+(-1)*Uc_A_1;%%S2 S4��״̬�µ�Vdc
%     da=(-1)*(Uc_C_1/Uc_B_1);db=(-1)*(Uc_A_1/Uc_B_1); 
    sita=pi+sita;
    da=(sin(sita-pi/2));db=(sin(5*pi/6-(sita)));dc=1-da-db;
     Vdc=da*Vdc_a+db*Vdc_b;
    Va=[0 1 0 0 0 1];Vb=[0 1 0 1 0 0];
    if(abs(Uc_C_1)<=abs(Uc_A_1))
        Vc=[0 0 1 0 0 1];
    else
        Vc=[1 0 0 1 0 0];
    end
    k=3;
     %%%%%%%S1 S5 S6 �Լ�S1 S2 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(2-1*0-1*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia1=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib1=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic1=[0 0 0 1 1 1]; Id1=[1 1 1 0 0 0];
      Da1=(gb1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Db1=(ga1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Dc1=1-Da1-Db1;
    Q1=(ga1*gb1*gc1)/(ga1*gb1+gb1*gc1+gc1*ga1);
    
         %%%%%%%S1 S2 S6 �Լ�S4 S2 S6 ����
    Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end
   
         %%%%%%%S4 S2 S6 �Լ�S4 S2 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end    
 
          
         %%%%%%%S4 S2 S3 �Լ�S4 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end  
    
         %%%%%%%S4 S5 S3 �Լ�S1 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end 
    
         %%%%%%%S1 S5 S3 �Լ�S1 S5 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end

elseif Uc_A_1<0&&Uc_B_1>0&&Uc_C_1>=0
    sita=pi+sita;
    %������ѡ�� S2 S4 �Լ� S3 S4
    Vdc_a=1*Uc_B_1+(-1)*Uc_A_1;%%S2 S4��״̬�µ�Vdc
    Vdc_b=1*Uc_C_1+(-1)*Uc_A_1;%%S3 S4��״̬�µ�Vdc
%     da=(-1)*(Uc_B_1/Uc_A_1);db=(-1)*(Uc_C_1/Uc_A_1);
da=(sin(sita-(5*pi/6)));db=(sin(7*pi/6-(sita))); dc=1-da-db;
     Vdc=da*Vdc_a+db*Vdc_b;
    Va=[0 1 0 1 0 0];Vb=[0 0 1 1 0 0];
    if(abs(Uc_B_1)<=abs(Uc_C_1))
        Vc=[0 1 0 0 1 0 ];
    else
        Vc=[0 0 1 0 0 1];
    end
        k=4;
     %%%%%%%S1 S5 S6 �Լ�S1 S2 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(2-1*0-1*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia1=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib1=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic1=[0 0 0 1 1 1]; Id1=[1 1 1 0 0 0];
      Da1=(gb1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Db1=(ga1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Dc1=1-Da1-Db1;
    Q1=(ga1*gb1*gc1)/(ga1*gb1+gb1*gc1+gc1*ga1);
    
         %%%%%%%S1 S2 S6 �Լ�S4 S2 S6 ����
    Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end
   
         %%%%%%%S4 S2 S6 �Լ�S4 S2 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end    
 
          
         %%%%%%%S4 S2 S3 �Լ�S4 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end  
    
         %%%%%%%S4 S5 S3 �Լ�S1 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end 
    
         %%%%%%%S1 S5 S3 �Լ�S1 S5 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end

elseif Uc_C_1>0&&Uc_B_1<=0&&Uc_A_1<0
    sita=pi+sita;
    %������ѡ�� S3 S4 �Լ� S3 S5
    Vdc_a=1*Uc_C_1+(-1)*Uc_A_1;%%S3 S4��״̬�µ�Vdc
    Vdc_b=1*Uc_C_1+(-1)*Uc_B_1;%%S3 S5��״̬�µ�Vdc
%     da=(-1)*(Uc_A_1/Uc_C_1);db=(-1)*(Uc_B_1/Uc_C_1);
   da=(sin(sita-7*pi/6)); db=sin(3*pi/2-(sita)); dc=1-da-db;
     Vdc=da*Vdc_a+db*Vdc_b;
    Va=[0 0 1 1 0 0];Vb=[0 0 1 0 1 0];
    if(abs(Uc_A_1)<=abs(Uc_B_1))
        Vc=[1 0 0 1 0 0];
    else
        Vc=[0 1 0 0 1 0];
    end
    k=5;
     %%%%%%%S1 S5 S6 �Լ�S1 S2 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(2-1*0-1*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia1=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib1=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic1=[0 0 0 1 1 1]; Id1=[1 1 1 0 0 0];
      Da1=(gb1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Db1=(ga1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Dc1=1-Da1-Db1;
    Q1=(ga1*gb1*gc1)/(ga1*gb1+gb1*gc1+gc1*ga1);
    
         %%%%%%%S1 S2 S6 �Լ�S4 S2 S6 ����
    Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end
   
         %%%%%%%S4 S2 S6 �Լ�S4 S2 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end    
 
          
         %%%%%%%S4 S2 S3 �Լ�S4 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end  
    
         %%%%%%%S4 S5 S3 �Լ�S1 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end 
    
         %%%%%%%S1 S5 S3 �Լ�S1 S5 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end  

  
else

    %������ѡ�� S3 S5 �Լ� S1 S5
    Vdc_a=1*Uc_C_1+(-1)*Uc_B_1;%%S3 S5��״̬�µ�Vdc
    Vdc_b=1*Uc_A_1+(-1)*Uc_B_1;%%S1 S5��״̬�µ�Vdc
%     da=(-1)*(Uc_C_1/Uc_B_1);db=(-1)*(Uc_A_1/Uc_B_1);
da=(sin(pi/2+(sita))); db=(sin(pi/3-(pi/2+(sita)))); dc=1-da-db;
     Vdc=da*Vdc_a+db*Vdc_b;
    Va=[0 0 1 0 1 0];Vb=[1 0 0 0 1 0];
    if(abs(Uc_C_1)<=abs(Uc_A_1))
        Vc=[0 0 1 0 0 1];
    else
        Vc=[1 0 0 1 0 0];
    end
    k=6;
     %%%%%%%S1 S5 S6 �Լ�S1 S2 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(2-1*0-1*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia1=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib1=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc1=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic1=[0 0 0 1 1 1]; Id1=[1 1 1 0 0 0];
      Da1=(gb1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Db1=(ga1*gc1)/(gb1*gc1+ga1*gc1+gb1*ga1);Dc1=1-Da1-Db1;
    Q1=(ga1*gb1*gc1)/(ga1*gb1+gb1*gc1+gc1*ga1);
    
         %%%%%%%S1 S2 S6 �Լ�S4 S2 S6 ����
    Uo_A_1=(1/3)*(2-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1+2*0)*Vdc;%%%%��伶S1 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 1 0 0 0 1];
     Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end
   
         %%%%%%%S4 S2 S6 �Լ�S4 S2 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1*0)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1*0)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2*0)*Vdc;%%%%��伶S4 S2 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 0 1 0 1];
     Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end    
 
          
         %%%%%%%S4 S2 S3 �Լ�S4 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1+2)*Vdc;%%%%��伶S4 S2 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 1 1 1 0 0];
     Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end  
    
         %%%%%%%S4 S5 S3 �Լ�S1 S5 S3 ����
    Uo_A_1=(1/3)*(2*0-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1*0+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1*0-1*0+2)*Vdc;%%%%��伶S4 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[0 0 1 1 1 0];
     Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ic2=[0 0 0 1 1 1]; Id2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end 
    
         %%%%%%%S1 S5 S3 �Լ�S1 S5 S6 ����
    Uo_A_1=(1/3)*(2-1*0-1)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2)*Vdc;%%%%��伶S1 S5 S3
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      ga2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ia2=[1 0 1 0 1 0];
     Uo_A_1=(1/3)*(2-1*0-1*0)*Vdc;Uo_B_1=(1/3)*(-1+2*0-1*0)*Vdc;Uo_C_1=(1/3)*(-1-1*0+2*0)*Vdc;%%%%��伶S1 S5 S6
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gb2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Ib2=[1 0 0 0 1 1];
     Uo_A_1=(1/3)*(2-1-1)*Vdc;Uo_B_1=(1/3)*(-1+2-1)*Vdc;Uo_C_1=(1/3)*(-1-1+2)*Vdc;%%%%��伶������
     Io_A_2=1/(Ls)*((Ls-Ts*Rs)*Io_A_1+Ts*Uo_A_1);
     Io_B_2=1/(Ls)*((Ls-Ts*Rs)*Io_B_1+Ts*Uo_B_1);
     Io_C_2=1/(Ls)*((Ls-Ts*Rs)*Io_C_1+Ts*Uo_C_1);
      gc2=abs(Io_A_g-Io_A_2)+abs(Io_B_g-Io_B_2)+abs(Io_C_g-Io_C_2);Id2=[0 0 0 1 1 1]; Ic2=[1 1 1 0 0 0];
      Da2=(gb2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Db2=(ga2*gc2)/(gb2*gc2+ga2*gc2+gb2*ga2);Dc2=1-Da2-Db2;
    Q2=(ga2*gb2*gc2)/(ga2*gb2+gb2*gc2+gc2*ga2);
    if(Q1<=Q2)
    else
        Ia1=Ia2;Ib1=Ib2;Ic1=Ic2;Id1=Id2;Da1=Da2;Db1=Db2;Dc1=Dc2;
    end     
   

end

% if(rem(t,Ts)>=0&&rem(t,Ts)<Ts/2*(Dc1/8))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%16steps
%     h=[Va Ic1 1 Vdc k];
% elseif(rem(t,Ts)>=Ts/2*(Dc1/8)&&rem(t,Ts)<Ts/2*(Dc1/8+da*Da1))
%     h=[Va Ia1 2 Vdc k];
% elseif(rem(t,Ts)>=Ts/2*(Dc1/8+da*Da1)&&rem(t,Ts)<Ts/2*(Dc1/8+da*Da1+da*Db1))
%     h=[Va Ib1 3 Vdc k];
% elseif(rem(t,Ts)>=Ts/2*(Dc1/8+da*Da1+da*Db1)&&rem(t,Ts)<Ts/2*(Dc1/8+da*Da1+da*Db1+Dc1/8))
%     h=[Va Id1 4 Vdc k];
% elseif(rem(t,Ts)>=Ts/2*(Dc1/8+da*Da1+da*Db1+Dc1/8)&&rem(t,Ts)<Ts/2*(Dc1/8+da*Da1+da*Db1+Dc1/4))
%     h=[Vb Id1 5 Vdc k];
% elseif(rem(t,Ts)>=Ts/2*(Dc1/8+da*Da1+da*Db1+Dc1/4)&&rem(t,Ts)<Ts/2*(Dc1/8+da*Da1+da*Db1+Dc1/4+db*Da1))
%     h=[Vb Ib1 6 Vdc k];   
% elseif(rem(t,Ts)>=Ts/2*(Dc1/8+da*Da1+da*Db1+Dc1/4+db*Da1)&&rem(t,Ts)<Ts/2*(Dc1/8+da*Da1+da*Db1+Dc1/4+db*Da1+db*Db1))
%     h=[Vb Ia1 7 Vdc k]; 
% elseif(rem(t,Ts)>=Ts/2*(Dc1/8+da*Da1+da*Db1+Dc1/4+db*Da1+db*Db1)&&rem(t,Ts)<(Ts/2))
%     h=[Vb Ic1 8 Vdc k];        
% elseif(rem(t,Ts)>=(Ts/2)&&rem(t,Ts)<(Ts/2+Ts/2*Dc1/8))
%     h=[Vb Ic1 9 Vdc k];   
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*Dc1/8)&&rem(t,Ts)<(Ts/2+Ts/2*(Dc1/8+db*Db1)))
%     h=[Vb Ia1 10 Vdc k];       
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(Dc1/8+db*Db1))&&rem(t,Ts)<(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1)))
%     h=[Vb Ib1 11 Vdc k];  
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1))&&rem(t,Ts)<(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1+Dc1/8)))
%     h=[Vb Id1 12 Vdc k]; 
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1+Dc1/8))&&rem(t,Ts)<(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1+Dc1/4)))
%     h=[Va Id1 13 Vdc k];
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1+Dc1/4))&&rem(t,Ts)<(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1+Dc1/8+da*Db1)))
%     h=[Va Ib1 14 Vdc k];
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1+Dc1/8+da*Db1))&&rem(t,Ts)<(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1+Dc1/8+da*Db1+da*Da1)))
%     h=[Va Ia1 15 Vdc k];
% elseif(rem(t,Ts)>=(Ts/2+Ts/2*(Dc1/8+db*Db1+db*Da1+Dc1/8+da*Db1+da*Da1))&&rem(t,Ts)<(Ts))
%     h=[Va Ic1 16 Vdc k];
% end
% 

% if(rem(t,Ts)>=0&&rem(t,Ts)<Ts*(Dc1/4))%%%%%%%%%%%%%%%%%%%%%%%%%%eight steps 
%     h=[Va Ic1 da sita k];
% elseif(rem(t,Ts)>=Ts*(Dc1/4)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1))
%     h=[Va Ia1 da sita k];
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1))
%     h=[Va Ib1 da sita k];
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1+Dc1/4))
%     h=[Va Id1 da sita k];
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1+Dc1/4)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2))
%     h=[Vb Id1 da sita k];
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2+db*Da1))
%     h=[Vb Ib1 da sita k];   
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2+db*Da1)&&rem(t,Ts)<Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2+db*Da1+db*Db1))
%     h=[Vb Ia1 da sita k]; 
% elseif(rem(t,Ts)>=Ts*(Dc1/4+da*Da1+da*Db1+Dc1/2+db*Da1+db*Db1)&&rem(t,Ts)<(Ts))
%     h=[Vb Ic1 da sita k]; 
% end

if(rem(t,2*Ts)>=0&&rem(t,2*Ts)<Ts*(da*Da1))
    h=[Va Ia1 da 1 k];
elseif(rem(t,2*Ts)>=(Ts*(da*Da1))&&rem(t,2*Ts)<Ts*(da*Da1+da*Db1))
    h=[Va Ib1 da 2 k];
elseif(rem(t,2*Ts)>=(Ts*(da*Da1+da*Db1))&&rem(t,2*Ts)<Ts*(da*Da1+da*Db1+db*Db1))
    h=[Vb Ib1 da 3 k];
elseif(rem(t,2*Ts)>=(Ts*(da*Da1+da*Db1+db*Db1))&&rem(t,2*Ts)<Ts*(da*Da1+da*Db1+db*Db1+db*Da1))
    h=[Vb Ia1 da 4 k];
elseif(rem(t,2*Ts)>=Ts*(da*Da1+da*Db1+db*Db1+db*Da1)&&rem(t,2*Ts)<Ts)
    h=[Vc Ic1 da 5 k];
elseif(rem(t,2*Ts)>=Ts&&rem(t,2*Ts)<Ts+Ts*dc)
    h=[Vc Id1 da 6 k];
elseif(rem(t,2*Ts)>=(Ts+Ts*dc)&&rem(t,2*Ts)<(Ts+Ts*(dc+db*Da1)))
    h=[Vb Ib1 da 7 k];
elseif(rem(t,2*Ts)>=(Ts+Ts*(dc+db*Da1))&&rem(t,2*Ts)<(Ts+Ts*(dc+db*Da1+db*Db1)))
    h=[Vb Ia1 da 8 k];
elseif(rem(t,2*Ts)>=(Ts+Ts*(dc+db*Da1+db*Db1))&&rem(t,2*Ts)<(Ts+Ts*(dc+db*Da1+db*Db1+da*Db1)))
    h=[Va Ia1 da 9 k];
elseif(rem(t,2*Ts)>=(Ts+Ts*(dc+db*Da1+db*Db1+da*Db1))&&rem(t,2*Ts)<(Ts+Ts))
    h=[Va Ib1 da 10 k];
end
sys=h;
end