clc;
clear;
close all;
exp_LTP_raw=[0	7.86695E-6
1	2.24891E-6
2	8.84676E-6
3	9.96004E-6
4	1.06904E-5
5	1.1348E-5
6	1.19418E-5
7	1.24855E-5
8	1.29958E-5
9	1.34718E-5
10	1.39187E-5
11	1.47366E-5
12	1.47575E-5
13	1.51459E-5
14	1.5533E-5
15	1.58968E-5
16	1.62486E-5
17	1.65903E-5
18	1.69211E-5
19	1.72441E-5
20	1.75601E-5
21	1.82655E-5
22	1.81647E-5
23	1.84538E-5
24	1.87396E-5
25	1.90188E-5
26	1.92893E-5
27	1.95577E-5
28	1.98181E-5
29	2.00733E-5
30	2.03237E-5
31	2.09803E-5
32	2.08138E-5
33	2.10462E-5
34	2.12819E-5
35	2.15112E-5
36	2.17392E-5
37	2.19587E-5
38	2.21761E-5
39	2.23893E-5
40	2.25965E-5
41	2.32201E-5
42	2.30113E-5
43	2.32147E-5
44	2.34092E-5
45	2.36061E-5
46	2.38026E-5
47	2.39904E-5
48	2.4179E-5
49	2.43629E-5
50	2.45477E-5
51	2.51389E-5
52	2.49062E-5
53	2.50819E-5
54	2.52532E-5
55	2.54266E-5
56	2.55928E-5
57	2.57581E-5
58	2.59232E-5
59	2.60873E-5
60	2.62477E-5
61	2.68261E-5
62	2.65641E-5
63	2.67229E-5];

exp_LTD_raw=[  63	2.68749E-5
62	2.82636E-5
61	2.49027E-5
60	2.24763E-5
59	2.06688E-5
58	1.92506E-5
57	1.81048E-5
56	1.71485E-5
55	1.63384E-5
54	1.56547E-5
53	1.50517E-5
52	1.45233E-5
51	1.40539E-5
50	1.36337E-5
49	1.32635E-5
48	1.29236E-5
47	1.26223E-5
46	1.23462E-5
45	1.20853E-5
44	1.18515E-5
43	1.1635E-5
42	1.14354E-5
41	1.12495E-5
40	1.10751E-5
39	1.09134E-5
38	1.07593E-5
37	1.06166E-5
36	1.04847E-5
35	1.03537E-5
34	1.02334E-5
33	1.0114E-5
32	1.00027E-5
31	9.89675E-6
30	9.8013E-6
29	9.70227E-6
28	9.61106E-6
27	9.52352E-6
26	9.43806E-6
25	9.35211E-6
24	9.26972E-6
23	9.19897E-6
22	9.12143E-6
21	9.04898E-6
20	8.97929E-6
19	8.91202E-6
18	8.84679E-6
17	8.78381E-6
16	8.72225E-6
15	8.65993E-6
14	8.60153E-6
13	8.54144E-6
12	8.48448E-6
11	8.43188E-6
10	8.37765E-6
9	8.33127E-6
8	8.27874E-6
7	8.22878E-6
6	8.18359E-6
5	8.12955E-6
4	8.08398E-6
3	8.038E-6
2	7.99381E-6
1	7.9541E-6
0	7.91249E-6];
xf_ltp = floor(max(exp_LTP_raw(:,1)));
xf_ltd = floor(max(exp_LTD_raw(:,1)));
x_step_ltp = 1/xf_ltp;
x_step_ltd = 1/xf_ltd;
yf_ltp = max(exp_LTP_raw(:,2));
yf_ltd = max(exp_LTD_raw(:,2));
yi_ltp = min(exp_LTP_raw(:,2));
yi_ltd = min(exp_LTD_raw(:,2));

exp_LTP(:,1) = exp_LTP_raw(:,1)/xf_ltp;
exp_LTD(:,1) = exp_LTD_raw(:,1)/xf_ltd;
exp_LTP(:,2) = (exp_LTP_raw(:,2)-yi_ltp)/(yf_ltp-yi_ltp);
exp_LTD(:,2) = (exp_LTD_raw(:,2)-yi_ltd)/(yf_ltd-yi_ltd);

plot(exp_LTP(:,1), exp_LTP(:,2), 'bo', 'LineWidth', 1);
hold on;
plot(exp_LTD(:,1), exp_LTD(:,2), 'ro', 'LineWidth', 1);

xf = 1;
A_LTP = 0.4;
B_LTP = 1./(1-exp(-1./A_LTP));
A_LTD = -0.2;
B_LTD = 1./(1-exp(-1./A_LTD));

% LTP fitting
var_amp = 0.001;    % LTP cycle-to-cycle variation
rng(103);
x_ltp(1) = 0;
y_ltp(1) = 0;
for n=1:1/x_step_ltp+1
    x_ltp(n+1) = x_ltp(n)+x_step_ltp;
    y_ltp(n+1) = B_LTP(1)*(1-exp(-x_ltp(n+1)/A_LTP(1)));
    delta_y = (y_ltp(n+1)-y_ltp(n)) + randn*var_amp;
    y_ltp(n+1) = y_ltp(n) + delta_y;   
    if y_ltp(n+1)>=1
        y_ltp(n+1)=1;
    elseif y_ltp(n+1)<=0
        y_ltp(n+1)=0;
    end
    x_ltp(n+1) = -A_LTP(1)*log(1-(y_ltp(n+1))/B_LTP(1));
end
plot((0:n-1)/(n-1), y_ltp(1:n), 'b', 'linewidth', 2);

% LTD fitting
var_amp = 0.001;    % LTD cycle-to-cycle variation
rng(898);
x_ltd(1) = 1;
y_ltd(1) = 1;
for n=1:1/x_step_ltd+1
    x_ltd(n+1) = x_ltd(n)-x_step_ltd;
    y_ltd(n+1) = B_LTD(1)*(1-exp(-x_ltd(n+1)/A_LTD(1)));
    delta_y = (y_ltd(n+1)-y_ltd(n)) + randn*var_amp;
    y_ltd(n+1) = y_ltd(n) + delta_y;
    if y_ltd(n+1)>=1
        y_ltd(n+1)=1;
    elseif y_ltd(n+1)<=0
        y_ltd(n+1)=0;
    end
    x_ltd(n+1) = -A_LTD(1)*log(1-(y_ltd(n+1))/B_LTD(1));
end
x_start = numel(x_ltd(:));
x_end = numel(x_ltd(:)) - n;
plot((n-1:-1:0)/(n-1), y_ltd(1:n), 'r', 'linewidth', 2);

xlabel('Normalized Pulse #');
ylabel('Normalized Conductance');
legend('Exp data (LTP)','Exp data (LTD)', 'Fit (LTP)', 'Fit (LTD)', 'location', 'southeast');