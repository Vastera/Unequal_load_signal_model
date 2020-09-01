《行星轮不均载的行星齿轮箱振动信号模型》论文信号模型代码 作者：马浩群，冯志鹏
CopyRight:vastera@163.com

%%%% 主程序 %%%%%%
Simulation_Script.m

%%%% 参数设置 %%%%
T_s：输入负载大小 单位Nm 默认T_s=3e2;
epsilon_i：行星轮轴的周向位置偏差，向量长度代表行星轮数，每一个位置的值分别对应一个行星轮周向偏差大小 单位是弧度
例如：epsilon_i=[0,0,0,0/180*pi,0,0];%6个行星轮的正常情况
或者：epsilon_i=[0,0.5/180*pi,0,0];%4个行星轮中第二个有误差

%%%% 函数介绍 %%%%
FigureSetting.m 调整图形显示
III.m 生成冲击函数
Load_sharing_coef.m 计算载荷分布比
Resonant_frequency.m 生成固有振动
Transfer_length.m 时变传递路径长度


Code for '"A vibration signal model of planetary gearboxes under
% uneven load sharing among planets" author: Haoqun Ma,Zhipeng Feng. email:vastera@163.com
%%%% Main program %%%%%
Simulation_Script.m

%%%% Parameter setting %%%%
T_s: Input load size unit Nm default T_s=3e2;
epsilon_i：Circumferential position error,The vector length stands for the number of planetary rounds, with each position corresponding to a circumferential error of one planet in radians
For example: epsilon_i = [0,0,0,0/180 * PI, 0, 0];%Normal condition of 6 planetary wheels
Or: epsilon_i = [0,0.5/180 * PI, 0, 0];% The second of the four planet has an angular error

%%%% Function introduction %%%%
FigureSetting.m Adjust the graph display
III.m generates Impact function
Load_sharing_coef.m Calculates the load distribution ratio
Amplitude.m Generated natural vibrations
Transfer_leng.m Time-varying transfer path length