function mpc = Extented_case()
% 14-bus : 400kV slack + 63kV meshed network + 20kV distribution loads and generators
% Without shunt susceptances for the initial network
% Distribution loads and generators placed directly on the secondary bus of each distribution OLTC
% No explicit 20kV feeder lines: aggregated loads and generators share the OLTC secondary bus

mpc.version = '2';
mpc.baseMVA = 100;

%% BUS DATA
% bus  type  Pd    Qd    Gs  Bs  area  Vm      Va   baseKV  zone  Vmax  Vmin
mpc.bus = [
% ----------------- Transmission network -----------------

%--- 400 kV network
 1    3     0     0     0   0   1     1.0125   0   400      1     1.05  0.95;
 2    1     0     0     0   0   1     1.000    0   400      1     1.05  0.95;
 3    1     0     0     0   0   1     1.000    0   400      1     1.05  0.95; 

%--- 63 kV meshed network
4    1     0     0     0   0   1     1.000   0    63      1     1.05  0.95;
5    1     0     0     0   0   1     1.000   0    63      1     1.05  0.95;
6    1     0     0     0   0   1     1.000   0    63      1     1.05  0.95;
7    1     0     0     0   0   1     1.000   0    63      1     1.05  0.95;
8    1     0     0     0   0   1     1.000   0    63      1     1.05  0.95;
9    1     0     0     0   0   1     1.000   0    63      1     1.05  0.95;
  
 %--- Generator connected via transformer from bus 8 (63kV)
10    2     0     0     0   0   1     1.010   0    20      1     1.05  0.95;

% ----------------- Distribution network -----------------

%--- 20 kV network (OLTC secondary buses — loads and generators placed directly here)

%--- Distribution connected to bus 4 (20kV)
11    1    50     14    0   0  1     1.000   0    20      1     1.05  0.95; % Aggregated load (n=5)

%--- Distribution connected to bus 6 (20kV)
12    2    30    10    0   0  1     1.010   0    20      1     1.05  0.95; % Aggregated load (n=5) + generator (n=2)

%--- Distribution connected to bus 7 (20kV)
13    2    40    12    0   0  1     1.010   0    20      1     1.05  0.95; % Aggregated load (n=5) + generator (n=2)

%--- Distribution connected to bus 9 (20kV)
14    1    60     20    0   0  1     1.000   0    20      1     1.05  0.95; % Aggregated load (n=5)
];

%% GENERATOR DATA
% bus   Pg    Qg   Qmax  Qmin   Vg      mBase  status  Pmax   Pmin
mpc.gen = [
1       30    0   600  -600   1.0125  100     1      1000     0;   % 400kV slack
10      60    0    50   -30   1.010   100     1       100    10;   % Subtransmission generator (single unit)
12      40    13   20   -20   1.010   100     1        60     0;   % Aggregated distribution generator (n=2)
13      50    17   20   -20   1.010   100     1        60     0;   % Aggregated distribution generator (n=2)
];

%% BRANCH DATA
% from_bus to_bus   r        x        b       rateA  rateB  rateC  ratio(Vf/Vt)  angle  status
mpc.branch = [

% ----------------- Transmission network -----------------

%--- 400 kV double-circuit lines (80km)
1   2    0.001  0.015  0.384  1000   1000   1000     0      0      1;
1   2    0.001  0.015  0.384  1000   1000   1000     0      0      1;

1   3    0.001  0.015  0.384  1000   1000   1000     0      0      1;
1   3    0.001  0.015  0.384  1000   1000   1000     0      0      1;

%--- OLTC (400 -> 63 kV)
2    4    0.003    0.130    0.001        600    600    600    1.00   0      1;
3    6    0.002    0.150    0.001        600    600    600    1.00   0      1;
 
 %--- 63 kV double-circuit lines (20km)
4    5    0.06   0.21  0.0022  80  80   80   0      0      1;
4    5    0.06   0.21  0.0022  80  80   80   0      0      1;

4    7    0.06   0.21  0.0022  80  80   80   0      0      1;
4    7    0.06   0.21  0.0022  80  80   80   0      0      1;

5    6    0.06   0.21  0.0022  80  80   80   0      0      1;
5    6    0.06   0.21  0.0022  80  80   80   0      0      1;

6    9    0.06   0.21  0.0022  80  80   80   0      0      1;
6    9    0.06   0.21  0.0022  80  80   80   0      0      1;

7    8    0.06   0.21  0.0022  80  80   80   0      0      1;
7    8    0.06   0.21  0.0022  80  80   80   0      0      1;

8    9    0.06   0.21  0.0022  80  80   80   0      0      1;
8    9    0.06   0.21  0.0022  80  80   80   0      0      1;

%--- Transformer (63 -> 20 kV)
8    10    0.004    0.09    0.004  120    120    120    1.00   0      1;

% ----------------- Distribution network -----------------

%--- OLTC (63 -> 20 kV)
4    11    0.004    0.09    0.004  120    120    120    1.00   0      1;
6    12    0.004    0.09    0.004  120    120    120    1.00   0      1;
7    13    0.004    0.09    0.004  120    120    120    1.00   0      1;
9    14    0.004    0.09    0.004  120    120    120    1.00   0      1;
];


end
