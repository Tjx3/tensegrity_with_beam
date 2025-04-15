function output=N_Ex(xyz,L_i)
global I E
I=pi*(xyz^4-(0.8*xyz)^4)/64;
mju_x=1; % mju计算长度系数，两端铰接取1.0
l_x=mju_x*L_i;  % 构件对截面主轴x和y的计算长度
lamda_x=l_x/sqrt(I/xyz); % 构件的有效长细比
N_Ex=pi^2*E*xyz/(1.1*lamda_x^2);
output= N_Ex;
end
