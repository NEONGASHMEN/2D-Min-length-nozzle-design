clear
clc
close all

gamma = 1.4;                        %%Ratio of sp. heats
Po = 2.27e6;                        %%Combustion chamber pressure in Pa
To = 1200;                          %%Combustion chamber temperature in K
Pa = 39368;                         %%Ambient pressure
throat_rad = 0.2;                   %%Radius of throat
rsltn = 20;                         %%No of wall divisions

disp('Please wait, the program is running......');

pr = Pa/Po;                         %%Pressure ratio
Me = mach_frm_pr(pr,gamma);         %%Exit mach no
theta_max = 0.5*pm(Me,gamma);       %%maximum wall angle (init wall angle in ML nozzle
unit_theta = theta_max/(rsltn);     %%Each angle of expansion wave

crdnts = zeros(length(0:unit_theta:theta_max),2);
chara = cell(3*rsltn,11);
%%%1 > pt_ref
%%%2 > x_cordinate
%%%3 > y_coordinate
%%%4 > theta
%%%5 > nu
%%%6 > mach
%%%7 > mu
%%%8 > C+
%%%9 > C-
%%%10 > K+
%%%11 > K-

%%Getting characteristic lines
i = 1;
for theta=unit_theta:unit_theta:theta_max       
    chara{i,1} = strcat('0,',num2str(i));
    chara{i,4} = theta*180/pi;
    
    nu = theta;
    chara{i,5} = nu*180/pi;
    
    syms x;
    M = double(abs(vpasolve(pm(x,gamma) == nu,x)));
    chara{i,6} = M;
    
    mu = mach_angle(M);
    chara{i,7} = mu*180/pi;
    
    C_plus = tan(theta + mu);
    C_minus = tan(theta - mu);
    chara{i,8} = C_plus;
    chara{i,9} = C_minus;
    
    K_minus = theta + nu;
    K_plus = theta - nu;
    chara{i,10} = K_plus*180/pi;
    chara{i,11} = K_minus*180/pi;
    
    chara{i,2} = 0;
    chara{i,3} = throat_rad;
    
    i = i + 1;
end

%%Getting axis points
j = 1;
for i = length(unit_theta:unit_theta:theta_max)+1:length(unit_theta:unit_theta:theta_max)*2
    chara{i,1} = j;
    
    chara{i,3} = 0;
    theta = 0;
    chara{i,4} = theta;
    K_minus = chara{j,11};
    chara{i,11} = K_minus;
    
    nu = K_minus - theta;
    chara{i,5} = nu;
    
    syms x;
    M = double(abs(vpasolve(pm(x,gamma) == nu*pi/180,x)));
    chara{i,6} = M;
    
    mu = mach_angle(M);
    chara{i,7} = mu*180/pi;
    
    C_plus = tan(theta + mu);
    C_minus = tan(theta - mu);
    chara{i,8} = C_plus;
    chara{i,9} = C_minus;
    
    K_plus = theta - nu;
    chara{i,10} = K_plus;
    
    syms x_c;
    x_crdnt = double(vpasolve(((0 - throat_rad)/(x_c - 0)) == C_minus,x_c));
    chara{i,2} = x_crdnt; 
    
    j = j+1;
end

k = 2*rsltn;
int_Kp = zeros(rsltn,rsltn);        %%wave index by wave no.    
int_Km = zeros(rsltn,rsltn);
int_theta = zeros(rsltn,rsltn);
int_nu = zeros(rsltn,rsltn);
int_M = zeros(rsltn,rsltn);
int_mu = zeros(rsltn,rsltn);
int_x = zeros(rsltn,rsltn);
int_y = zeros(rsltn,rsltn);

int_Cm = zeros(rsltn,rsltn);
int_Cp = zeros(rsltn,rsltn);

%%Getting interior points
for i = (rsltn):-1:1                %%i > no of pts in each wave
    
    l = rsltn - i + 1;              %%wave index
    for j = 1:i                     %%the pt in each wave
        
        if j ~= 1
            
            if (j == 2)
                int_Kp(l,j) = chara{l+rsltn,10};
                int_Cp(l,j) = tand(0.5*(chara{l+rsltn,4}+int_theta(l,j)) + 0.5*(chara{l+rsltn,7}+int_mu(l,j)));
                if l == 1
                    int_Km(l,j) = chara{j,11};
                    int_Cm(l,j) = tand(0.5*(chara{j,4}+int_theta(l,j)) - 0.5*(chara{j,7}+int_mu(l,j)));
                    syms x y
                    eqns = [(y - chara{j,3})/(x - chara{j,2})==int_Cm(l,j),(y - chara{l+rsltn,3})/(x - chara{l+rsltn,2})==int_Cp(l,j)];
                    sln = solve(eqns,x,y);
                    int_x(l,j) = sln.x;
                    int_y(l,j) = sln.y;
                else
                    int_Km(l,j) = int_Km(l-1,j+1);
                    int_Cm(l,j) = tand(0.5*(int_theta(l-1,j+1)+int_theta(l,j)) - 0.5*(int_mu(l-1,j+1)+int_mu(l,j)));
                    syms x y
                    eqns = [(y - int_y(l-1,j+1))/(x - int_x(l-1,j+1))==int_Cm(l,j),(y - chara{l+rsltn,3})/(x - chara{l+rsltn,2})==int_Cp(l,j)];
                    sln = solve(eqns,x,y);
                    int_x(l,j) = sln.x;
                    int_y(l,j) = sln.y;
                end
            else
                int_Kp(l,j) = int_Kp(l,j-1);
                int_Cp(l,j) = tand(0.5*(int_theta(l,j-1)+int_theta(l,j)) + 0.5*(int_mu(l,j-1)+int_mu(l,j)));
                if l == 1
                    int_Km(l,j) = chara{j,11};
                    int_Cm(l,j) = tand(0.5*(chara{j,4}+int_theta(l,j)) - 0.5*(chara{j,7}+int_mu(l,j)));
                    syms x y
                    eqns = [(y - chara{j,3})/(x - chara{j,2})==int_Cm(l,j),(y - int_y(l,j-1))/(x - int_x(l,j-1))==int_Cp(l,j)];
                    sln = solve(eqns,x,y);
                    int_x(l,j) = sln.x;
                    int_y(l,j) = sln.y;
                else
                    int_Km(l,j) = int_Km(l-1,j+1);
                    int_Cm(l,j) = tand(0.5*(int_theta(l-1,j+1) + int_theta(l,j)) - 0.5*(int_mu(l-1,j+1) + int_mu(l,j)));
                    syms x y
                    eqns = [(y - int_y(l-1,j+1))/(x - int_x(l-1,j+1))==int_Cm(l,j),(y - int_y(l,j-1))/(x - int_x(l,j-1))==int_Cp(l,j)];
                    sln = solve(eqns,x,y);
                    int_x(l,j) = sln.x;
                    int_y(l,j) = sln.y;
                end
            end
            
            int_theta(l,j) = 0.5*(int_Kp(l,j) + int_Km(l,j));
            int_nu(l,j) = 0.5*(int_Km(l,j) - int_Kp(l,j));
            syms var;
            int_M(l,j) = double(abs(vpasolve(pm(var,gamma) == int_nu(l,j)*pi/180,var)));
            int_mu(l,j) = mach_angle(int_M(l,j))*180/pi;
            
        end
        
    end 
end

%%Getting wall points
j = 1;
k = rsltn;
for i = length(unit_theta:unit_theta:theta_max)*2+1:length(unit_theta:unit_theta:theta_max)*3
    
    chara{i,1} = strcat('W',num2str(j));
    
    if j == rsltn
        chara{i,4} = chara{2*rsltn,4};
        chara{i,5} = chara{2*rsltn,5};
        chara{i,6} = chara{2*rsltn,6};
        chara{i,7} = chara{2*rsltn,7};
        chara{i,10} = chara{2*rsltn,10};
        chara{i,11} = chara{2*rsltn,11};
    else
        chara{i,4} = int_theta(j,k);
        chara{i,5} = int_nu(j,k);           %%Since theta_wall = theta_awaypt, K+ chara being equal => nu is same
        chara{i,6} = int_M(j,k);
        chara{i,7} = int_mu(j,k);
        chara{i,10} = int_Kp(j,k);
        chara{i,11} = int_Km(j,k);
    end
    chara{i,8} = tand(chara{i,4} + chara{i,7});
    
    j = j + 1;
    k = k - 1;
end

j = rsltn;
for i = 1:rsltn
    
    if i == 1
        syms x y
        eqns = [(y - int_y(i,j))/(x - int_x(i,j)) == chara{i+2*rsltn,8},(y - throat_rad)/x == tand(0.5*(theta_max*180/pi + chara{i+2*rsltn,4}))];
        sln = solve(eqns,x,y);
        chara{i+2*rsltn,2} = double(sln.x);
        chara{i+2*rsltn,3} = double(sln.y);
    else
        syms x y
        eqns = [(y - int_y(i,j))/(x - int_x(i,j)) == chara{i+2*rsltn,8},(y - chara{i+2*rsltn-1,3})/(x - chara{i+2*rsltn-1,2}) == tand(0.5*(chara{i+2*rsltn-1,4} + chara{i+2*rsltn,4}))];
        sln = solve(eqns,x,y);
        chara{i+2*rsltn,2} = double(sln.x);
        chara{i+2*rsltn,3} = double(sln.y);
    end
    
    
    j = j - 1;
end
    
%%Print stuff
disp('  ');
disp('Done !');

%%Plot stuff
fig1 = figure();
for i = length(unit_theta:unit_theta:theta_max)+1:length(unit_theta:unit_theta:theta_max)*2
    line([0 chara{i,2}],[throat_rad 0],'Color','blue')
end
hold on;
for i = 1:rsltn
    line([chara{i+rsltn,2} chara{i+2*rsltn,2}],[chara{i+rsltn,3} chara{i+2*rsltn,3}],'Color','red');
end
hold on;
line([0 chara{2*rsltn+1,2}],[throat_rad chara{2*rsltn+1,3}],'Color','k','LineWidth',2);
hold on;
for i = 1:(rsltn-1)
    line([chara{2*rsltn+i,2} chara{2*rsltn+i+1,2}],[chara{2*rsltn+i,3} chara{2*rsltn+i+1,3}],'Color','k','LineWidth',2);
end

%%Write stuff
crdnts = zeros(rsltn+1,3);
crdnts(:,3) = 0;
crdnts(1,1) = 0;
crdnts(1,2) = throat_rad;
for i = 1:rsltn
    crdnts(i+1,1) = chara{2*rsltn+i,2};
    crdnts(i+1,2) = chara{2*rsltn+i,3};
end
save 'raw_crdnts.dat' crdnts -ascii;
