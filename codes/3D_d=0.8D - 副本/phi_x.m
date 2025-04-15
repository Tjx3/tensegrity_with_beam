function output=phi_x(xyz,L_i)
global I f E phi_x
I=pi*(xyz^4-(0.8*xyz)^4)/64;
mju_x=1; % mju计算长度系数，两端铰接取1.0
l_x=mju_x*L_i;  % 构件对截面主轴x和y的计算长度
lamda_x=l_x/sqrt(I/xyz); % 构件的有效长细比
lamda_nx=(lamda_x/pi)*sqrt(f/E);        
if lamda_nx<=0.215
    phi_x=1-0.41*lamda_nx^2; % 轴心受压构件的稳定系数
else
    phi_x=(1/(2*lamda_nx^2))*((0.986+0.152*lamda_nx+lamda_nx^2)-sqrt(0.986+0.152*lamda_nx+lamda_nx^2));
end
output= phi_x;
end