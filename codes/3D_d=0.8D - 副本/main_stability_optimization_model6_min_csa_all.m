clc,clear,close all
global T_diag V2 A ne w L z d_min3 d1 phi_x   
[~,N,C_s_in,C_b_in,nodes,groups,num_qa] = input_var;
%% 手动输入参数N（全部坐标节点的坐标向量x,y,z） 
%% 手动输入参数C_s_in和C_b_in得到C（索1，2，3，4和杆1，2的连接索引）
C_b = tenseg_ind2C(C_b_in,N);C_s = tenseg_ind2C(C_s_in,N);
C=[C_s;C_b];    % C：关系矩阵
[ne,bb]=size(C);    % ne：杆/梁单元数量;nn：构件节点数量
q=eye(bb);
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
    dz=The_rod_element_direction_vector_tem(3);
    if dx&&dy==0
        T_tem=1/l_tem*[dx,-dy*l_tem/sqrt(dx^2+dy^2),-dx*dz/sqrt(dx^2+dy^2),0,0,0,0,0,0,0,0,0;
                   dy,dx*l_tem/sqrt(dx^2+dy^2),-dy*dz/sqrt(dx^2+dy^2),0,0,0,0,0,0,0,0,0;
                   dz,0,sqrt(dx^2+dy^2),0,0,0,0,0,0,0,0,0;
                   0,0,0,dx,-dy*l_tem/sqrt(dx^2+dy^2),-dx*dz/sqrt(dx^2+dy^2),0,0,0,0,0,0;
                   0,0,0,dy,dx*l_tem/sqrt(dx^2+dy^2),-dy*dz/sqrt(dx^2+dy^2),0,0,0,0,0,0;
                   0,0,0,dz,0,sqrt(dx^2+dy^2),0,0,0,0,0,0;
                   0,0,0,0,0,0,dx,-dy*l_tem/sqrt(dx^2+dy^2),-dx*dz/sqrt(dx^2+dy^2),0,0,0;
                   0,0,0,0,0,0,dy,dx*l_tem/sqrt(dx^2+dy^2),-dy*dz/sqrt(dx^2+dy^2),0,0,0;
                   0,0,0,0,0,0,dz,0,sqrt(dx^2+dy^2),0,0,0;
                   0,0,0,0,0,0,0,0,0,dx,-dy*l_tem/sqrt(dx^2+dy^2),-dx*dz/sqrt(dx^2+dy^2);
                   0,0,0,0,0,0,0,0,0,dy,dx*l_tem/sqrt(dx^2+dy^2),-dy*dz/sqrt(dx^2+dy^2);
                   0,0,0,0,0,0,0,0,0,dz,0,sqrt(dx^2+dy^2)];    % T_tem：单元坐标转换矩阵
    else
        R=[0 0 1
           1 0 0
           0 -1 0];
        T_tem=blkdiag(R,R,R,R);
    end
    T_cell_all{i,:}=T_tem;
end
T_diag=blkdiag(T_cell_all{:})';
%% 铰节点处理
num=numel(nodes);
E=eye(length(nodes));I=zeros(length(nodes),1);
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
        A2_tem=[1 0 0 0 0 0 1 0 0 0 0 0
                0 1 0 0 0 0 0 1 0 0 0 0
                0 0 1 0 0 0 0 0 1 0 0 0
                0 0 0 1 0 0 0 0 0 1 0 0
                0 0 0 0 1 0 0 0 -L(i) 0 1 0
                0 0 0 0 0 1 0 L(i) 0 0 0 1
                ];
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

sa1=A1*V2   % sa1:等于0
ss=T_diag*V2;  % ss:局部坐标系下的杆端力
sa2=A2_diag*ss  % sa2:等于0

columns=size(V2,2);   % columns：V2的列数，即z的列数
% d=fsolve(@root,double(0));
%% 最小质量优化
% [z,volumes]=fmincon(@stability_test_func,zeros(columns,1),[],[],[],[]);volumes=volumes*1e-9; % 体积单位转换为立方米
[z,volumes]=fmincon(@stability_test_func,[1;1;1;1],[],[],[],[]);volumes=volumes*1e-9;
Phi_x=phi_x;D_out=d1;A_c=d_min3';F=w; % D_out：环形圆杆外径，A_c：杆件截面积最小值(mm2)，F：力(N)

1
