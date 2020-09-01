�������ֲ����ص����ǳ��������ź�ģ�͡������ź�ģ�ʹ��� ���ߣ����Ⱥ����־��
CopyRight:vastera@163.com

%%%% ������ %%%%%%
Simulation_Script.m

%%%% �������� %%%%
T_s�����븺�ش�С ��λNm Ĭ��T_s=3e2;
epsilon_i���������������λ��ƫ��������ȴ�������������ÿһ��λ�õ�ֵ�ֱ��Ӧһ������������ƫ���С ��λ�ǻ���
���磺epsilon_i=[0,0,0,0/180*pi,0,0];%6�������ֵ��������
���ߣ�epsilon_i=[0,0.5/180*pi,0,0];%4���������еڶ��������

%%%% �������� %%%%
FigureSetting.m ����ͼ����ʾ
III.m ���ɳ������
Load_sharing_coef.m �����غɷֲ���
Resonant_frequency.m ���ɹ�����
Transfer_length.m ʱ�䴫��·������


Code for '"A vibration signal model of planetary gearboxes under
% uneven load sharing among planets" author: Haoqun Ma,Zhipeng Feng. email:vastera@163.com
%%%% Main program %%%%%
Simulation_Script.m

%%%% Parameter setting %%%%
T_s: Input load size unit Nm default T_s=3e2;
epsilon_i��Circumferential position error,The vector length stands for the number of planetary rounds, with each position corresponding to a circumferential error of one planet in radians
For example: epsilon_i = [0,0,0,0/180 * PI, 0, 0];%Normal condition of 6 planetary wheels
Or: epsilon_i = [0,0.5/180 * PI, 0, 0];% The second of the four planet has an angular error

%%%% Function introduction %%%%
FigureSetting.m Adjust the graph display
III.m generates Impact function
Load_sharing_coef.m Calculates the load distribution ratio
Amplitude.m Generated natural vibrations
Transfer_leng.m Time-varying transfer path length