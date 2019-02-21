function [sys,x0,str,ts] =Us_1(t,x,u,flag)%%%预测Us(K+1)
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
sizes.NumOutputs     = 3;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [0.000005 0];     % 采样时间0.0001s,100000 Hz
end
% mdlOutputs
function sys=mdlOutputs(t,x,u)
persistent Us_p;
if isempty(Us_p)
    Us_p=[0 0 0];
end

persistent Us_pp;
if isempty(Us_pp)
    Us_pp=[0 0 0];
end

Us=[u(1) u(2) u(3)];

sys=[3*Us(1)-3*Us_p(1)+Us_pp(1) 3*Us(2)-3*Us_p(2)+Us_pp(2) 3*Us(3)-3*Us_p(3)+Us_pp(3)];
% sys=[Us(1) Us_p(1) Us_pp(1)];

Us_pp=[Us_p(1) Us_p(2) Us_p(3)];
Us_p=[Us(1) Us(2) Us(3)];
end