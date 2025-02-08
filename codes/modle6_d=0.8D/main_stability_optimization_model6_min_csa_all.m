clc,clear,close all
global T_diag V2 A ne w L z d_min3 d1 phi_x   
[~,N,C_s_in,C_b_in,nodes,groups,num_qa] = input_var;
%% 手动输入参数N（全部坐标节点的坐标向量x,y,z） 
%% 手动输入参数C_s_in和C_b_in得到C（索1，2，3，4和杆1，2的连接索引）
C_b = tenseg_ind2C(C_b_in,N);C_s = tenseg_ind2C(C_s_in,N);
C=[C_s;C_b];    % C：关系矩阵
[ne,nn]=size(C);    % ne：杆/梁单元数量;nn：构件节点数量
q=eye(nn);
% num=numel(N);   % num：节点数量
n=N(:); % 全部节点的非广义坐标向量（含x,y,z）
B=N*C'; % B:两根杆件方向向量，保证B中向量顺序与C_b_in相符合
L=diag(B'*B).^0.5; % L：杆件长度向量，与L=diag(diag(B'*B)).^0.5*ones(pole,1);等效
%% T_diag：各杆件单元坐标转换矩阵组成的对角矩阵的转置
for i=1:ne  % ne：同cell_nums=length(L); % cell_nums：杆/梁单元数量    
    l_tem=L(i);    
    The_rod_element_direction_vector_tem=kron(C(i,:),eye(3))*n;
    dx=The_rod_element_direction_vector_tem(1);
    dy=The_rod_element_direction_vector_tem(2);
    T_tem=1/l_tem*[dx,-dy,0,0,0,0;
               dy,dx,0,0,0,0;
               0,0,l_tem,0,0,0;
               0,0,0,dx,-dy,0;
               0,0,0,dy,dx,0;
               0,0,0,0,0,l_tem];    % T_tem：单元坐标转换矩阵    
    T_cell_all{i,:}=T_tem;
end
T_diag=blkdiag(T_cell_all{:})';
%% 铰节点处理
num=numel(nodes);
E=eye(length(nodes));I=zeros(length(nodes),length(6));
for i=1:length(groups)
    group=cell2mat(groups(i));
    for j=1:length(cell2mat(groups(i)))
        I(:,j)=E(:,group(j));
    end
    FF{:,i}=I;
end
C_i=cell2mat(FF);
%% Eqa
E_tem=eye(num);     % E_tem=eye(9);
% num_qa=[7 8 9 11 12 13]';   % num_qa：自由节点，节点总数（重码后）=9是这次错误的来源（*****）
% num_qb=setdiff([1:num]',num_qa);    % num_qb：约束节点
Eqa=E_tem(:,num_qa);
%% A1：整体坐标系下的节点平衡方程
A1=Eqa'*C_i;
%% A2：整体坐标系下的杆件单元平衡方程(内力平衡)
for i=1:ne
    % if i<=3
        A2_tem=[1 0 0 1 0 0
                0 1 0 0 1 0
                0 0 1 0 1 1];
        A2{i,:}=A2_tem;
end
A2_diag=blkdiag(A2{:});
A2=A2_diag*T_diag;
%% A：整体坐标系下的平衡方程
A=[A1;A2];
%% 奇异值分解
[U,S,V] = svd(A);
r=rank(A); 
U1=U(:,1:r);
U2=U(:,r+1:end);        % U1 is C(A_1g); U2 is N(A_1g') mechanism mode
% S1=S(1:r,1:r);                      % S1 is singular value of A_1g
V1=V(:,1:r);
%% V2：整体坐标系下各节点的力
V2=V(:,r+1:end);      % V1 is C(A_1g'); V2 is N(A_1g) self stress mode
columns=size(V2,2);   % columns：V2的列数，即z的列数
% d=fsolve(@root,double(0));
%% 最小质量优化
[z,volumes]=fmincon(@stability_test_func,zeros(columns,1),[],[],[],[]);volumes=volumes*1e-9; % 体积单位转换为立方米
Phi_x=phi_x;D_out=d1;A_c=d_min3';F=w; % D_out：环形圆杆外径，A_c：杆件截面积最小值(mm2)，F：力(N)


%% Plot结构轴力图
N=N';F_f=F*5e-5;%%
mm=1;nn=4;
mmm=1;
for i=1:11
    % 结构图
    for j=1:2
        var1(:,j)=C_b_in(i,j);
    end
    var2=var1(1);var3=var1(2);
    var4=N(var2,:);var5=N(var3,:);
    var6=var4(1);var7=var4(2);var8=var5(1);var9=var5(2);
    x=[var6,var8];y=[var7,var9];
    plot(x,y,'k','LineWidth',2)
    hold on
   % 轴力图
   % 坐标转换
    l_tem=L(i);    
    The_rod_element_direction_vector_tem=kron(C(i,:),eye(3))*n;
    dx=The_rod_element_direction_vector_tem(1);
    dy=The_rod_element_direction_vector_tem(2);
    cos=1/l_tem*dx;sin=1/l_tem*dy;
   % 横线
    xn1=F_f(mm)*sin+x(1);xn2=-F_f(nn)*sin+x(2);
    fn1=-F_f(mm)*cos+y(1);fn2=F_f(nn)*cos+y(2);
    mm=mm+6;nn=nn+6;
    xn=[xn1,xn2];
    fn=[fn1,fn2];
    plot(xn,fn,'b','LineWidth',1)
    hold on;
   % 边缘竖线
    xy_x=[x(1),xn(1)];xy_y=[y(1),fn(1)];
    plot(xy_x,xy_y,'b','LineWidth',1)
    hold on;
    yx_x=[x(2),xn(2)];yx_y=[y(2),fn(2)];
    plot(yx_x,yx_y,'b','LineWidth',1)
   % 已知两点求直线方程
   x=[var6,var8];y=[var7,var9];
   if x(1) ~= x(2)
       X=x(1):[x(2)-x(1)]/10:x(2);
       slope=(y(2)-y(1))/(x(2)-X(1));
       intercept=y(1)-slope*x(1);
       Y=slope*X+intercept;
       % plot(X,Y,'r','LineWidth',1)
       % hold on
   else
       X=[x(1),x(1),x(1),x(1),x(1),x(1),x(1),x(1),x(1),x(1),x(1)];
       Y=y(1):[y(2)-y(1)]/10:y(2);
       % plot(X,Y,'r','LineWidth',1)
       % hold on
   end
   % 中间竖线
   for i=2:length(X)-1
       xx=F_f(mmm)*sin+X(i);
       yy=-F_f(mmm)*cos+Y(i);

       xx_x=[X(i),xx];yy_y=[Y(i),yy];
       plot(xx_x,yy_y,'b','LineWidth',1)
       hold on;
   end
   % 标明轴力大小
   X_x=F_f(mmm)*sin+X(6);Y_y=-F_f(mmm)*cos+Y(6);
   str=mat2str(round(-F(mmm),2))+"N";
   text(X_x,Y_y,str,'Color','r','FontSize',8,'FontWeight','bold');
   mmm=mmm+6;
end
% Plot标识
xlabel({'L/m','轴力图'});ylabel('Fn/N');
title("张弦梁");
grid
1
% x=1;
% for i=1:10
%     % str1='杆件';
%     % str2=i;
%     % str3=str1+' '+str2;
%     str="杆件"+num2str(i);
%     figure('Name',str,'NumberTitle','off')
%     fn=F(x);l=L(i);
%     x=x+6;
%     m=0:l:l;
%     n=-fn+m-m;% fn左负轴拉力
%     plot(m,n,'-go', 'LineWidth', 1);
%     xlabel({'L/m','轴力图'});ylabel('Fn/N');
%     title("杆"+num2str(i));
%     grid
% end

