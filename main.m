clear
clc
close all

gamma = 1.4;                        %%Ratio of sp. heats
Po = 2.27e6;                        %%Combustion chamber pressure in Pa
To = 1200;                          %%Combustion chamber temperature in K
Pa = 39365;                         %%Ambient pressure
Ta = 243;
throat_rad = 0.2;                   %%Radius of throat
rsltn = 30;                         %%No of wall divisions

pr = Pa/Po;                         %%Pressure ratio
Me = mach_frm_pr(pr,gamma);         %%Exit mach no 
theta_max = 0.5*pm(Me,gamma);       %%maximum wall angle (init wall angle in ML nozzle
unit_theta = theta_max/(rsltn);     %%Each angle of expansion wave

chara = cell(rsltn,11);
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
bar = waitbar(0,'Calculating flow properties...');
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
    waitbar(i/rsltn,bar);
end
close(bar);

int_Kp = nan(rsltn,rsltn);        %%wave index by wave no.    
int_Km = nan(rsltn,rsltn);
int_theta = nan(rsltn,rsltn);
int_nu = nan(rsltn,rsltn);
int_M = nan(rsltn,rsltn);
int_mu = nan(rsltn,rsltn);
int_x = nan(rsltn,rsltn);
int_y = nan(rsltn,rsltn);

int_Cm = nan(rsltn,rsltn);
int_Cp = nan(rsltn,rsltn);

k = rsltn;
bar = waitbar(0,'Computing Characteristic mesh...');
for i = 1:rsltn                             %%C+ wave number - visualisation only
    
    for j = 1:k                             %%j th point on i th wave
        
        if i == 1
            
            if j == 1                       %%theta is zero at axis pnts
                int_theta(i,j) = 0;
                int_Km(i,j) = chara{j,11};
                int_Kp(i,j) = (-1)*int_Km(i,j);
                int_nu(i,j) = int_Km(i,j);

                syms M v;
                M = abs(double(vpasolve(pm(v,gamma) == int_nu(i,j)*pi/180,v)));
                int_M(i,j) = M;
                
                int_mu(i,j) = mach_angle(int_M(i,j))*180/pi;
                
                int_Cp(i,j) = nan;
                avg_theta_m = 0.5*(int_theta(i,j) + chara{j,4});
                avg_mu_m = 0.5*(int_mu(i,j) + chara{j,7});
                int_Cm(i,j) = tand(avg_theta_m - avg_mu_m);
                
                int_y(i,j) = 0;
                syms x v;
                x = double(solve((int_y(i,j)-chara{j,3})/(v-chara{j,2}) == int_Cm(i,j),v));
                int_x(i,j) = x;
            else                              %%Equate K+ & K- chara with previous pts on chara lines
                int_Kp(i,j) = int_Kp(i,j-1);
                int_Km(i,j) = chara{j,11};
                int_theta(i,j) = 0.5*(int_Km(i,j) + int_Kp(i,j));
                int_nu(i,j) = 0.5*(int_Km(i,j) - int_Kp(i,j));
                
                syms M v;
                M = abs(double(vpasolve(pm(v,gamma) == int_nu(i,j)*pi/180,v)));
                int_M(i,j) = M;
                
                int_mu(i,j) = mach_angle(int_M(i,j))*180/pi;
                
                avg_theta_m = 0.5*(int_theta(i,j) + chara{j,4});
                avg_mu_m = 0.5*(int_mu(i,j) + chara{j,7});
                int_Cm(i,j) = tand(avg_theta_m - avg_mu_m);
                avg_theta_p = 0.5*(int_theta(i,j) + int_theta(i,j-1));
                avg_mu_p = 0.5*(int_mu(i,j) + int_mu(i,j-1));
                int_Cp(i,j) = tand(avg_theta_p + avg_mu_p);
                
                syms x y;
                eqn_m = (y - chara{j,3})/(x - chara{j,2}) == int_Cm(i,j);
                eqn_p = (y - int_y(i,j-1))/(x - int_x(i,j-1)) == int_Cp(i,j);
                sltn = solve([eqn_m,eqn_p],x,y);
                int_x(i,j) = sltn.x;
                int_y(i,j) = sltn.y;
            end
            
        else
            
            if j == 1
                int_theta(i,j) = 0;
                int_Km(i,j) = int_Km(i-1,j+1);
                int_Kp(i,j) = (-1)*int_Km(i,j);
                int_nu(i,j) = int_Km(i,j);
                
                syms M v;
                M = abs(double(vpasolve(pm(v,gamma) == int_nu(i,j)*pi/180,v)));
                int_M(i,j) = M;
                
                int_mu(i,j) = mach_angle(int_M(i,j))*180/pi;
                
                avg_theta_m = 0.5*(int_theta(i,j) + int_theta(i-1,j+1));
                avg_mu_m = 0.5*(int_mu(i,j) + int_mu(i-1,j+1));
                int_Cm(i,j) = tand(avg_theta_m - avg_mu_m);
                int_Cp(i,j) = nan;
                
                int_y(i,j) = 0;
                syms x v;
                x = double(solve((int_y(i,j)-int_y(i-1,j+1))/(v-int_x(i-1,j+1)) == int_Cm(i,j),v));
                int_x(i,j) = x;
            else
                int_Kp(i,j) = int_Kp(i,j-1);
                int_Km(i,j) = int_Km(i-1,j+1);
                int_theta(i,j) = 0.5*(int_Km(i,j) + int_Kp(i,j));
                int_nu(i,j) = 0.5*(int_Km(i,j) - int_Kp(i,j));
                
                syms M v;
                M = abs(double(vpasolve(pm(v,gamma) == int_nu(i,j)*pi/180,v)));
                int_M(i,j) = M;
                
                int_mu(i,j) = mach_angle(int_M(i,j))*180/pi;
                
                avg_theta_m = 0.5*(int_theta(i,j) + int_theta(i-1,j+1));
                avg_mu_m = 0.5*(int_mu(i,j) + int_mu(i-1,j+1));
                int_Cm(i,j) = tand(avg_theta_m - avg_mu_m);
                avg_theta_p = 0.5*(int_theta(i,j) + int_theta(i,j-1));
                avg_mu_p = 0.5*(int_mu(i,j) + int_mu(i,j-1));
                int_Cp(i,j) = tand(avg_theta_p + avg_mu_p);
                
                syms x y;
                eqn_m = (y-int_y(i-1,j+1))/(x-int_x(i-1,j+1)) == int_Cm(i,j);
                eqn_p = (y-int_y(i,j-1))/(x-int_x(i,j-1)) == int_Cp(i,j);
                sltn = solve([eqn_m,eqn_p],x,y);
                int_x(i,j) = sltn.x;
                int_y(i,j) = sltn.y;
            end
            
        end
        
    end
    k = k - 1;
    waitbar(i/rsltn,bar);
    
end
close(bar);

j = rsltn;
bar = waitbar(0,'Obtaining contour points...');
for i = 1:rsltn                                     %%Slope of a contour line is tan of avrg theta at two end pts of the line
                                                    %%Theta & nu same for the wall pt as well as the most extrr intrr pt
    wall.theta(i) = int_theta(i,j);
    wall.Kp(i) = int_Kp(i,j);
    wall.Km(i) = int_Km(i,j);
    wall.nu(i) = int_nu(i,j);
    wall.M(i) = int_M(i,j);
    wall.mu(i) = int_mu(i,j);
    
    wall.Cp(i) = tand(wall.theta(i) + wall.mu(i));
    wall.Cm(i) = nan;
    
    if i == 1
        syms x y;
        eqn_lin = (y-chara{i,3})/(x-chara{i,2}) == tand(0.5*(theta_max*(180/pi)+wall.theta(i)));
        eqn_p = (y - int_y(i,j))/(x - int_x(i,j)) == wall.Cp(i);
        sltn = solve([eqn_lin,eqn_p],x,y);
        wall.x(i) = double(sltn.x);
        wall.y(i) = double(sltn.y);
    else
        syms x y;
        eqn_lin = (y-wall.y(i-1))/(x-wall.x(i-1)) == tand(0.5*(wall.theta(i-1)+wall.theta(i)));
        eqn_p = (y - int_y(i,j))/(x - int_x(i,j)) == wall.Cp(i);
        sltn = solve([eqn_lin,eqn_p],x,y);
        wall.x(i) = double(sltn.x);
        wall.y(i) = double(sltn.y);
    end
   
    j = j - 1;
    waitbar(i/rsltn,bar);
end
close(bar);

crdnts = zeros(rsltn,3);
crdnts(1,1) = 0;
crdnts(1,2) = throat_rad;
for i = 1:rsltn
    crdnts(i+1,1) = wall.x(i);
    crdnts(i+1,2) = wall.y(i);
end

%%Plot stuff
fig1 = figure('Name','Minimum Length Nozzle','NumberTitle','off');
sonic_line = line([0 0],[0 throat_rad],'Color','green','LineWidth',1);
hold on;
line([0 wall.x(1)],[throat_rad wall.y(1)],'Color','k','LineWidth',2);
hold on;
for i = 1:rsltn-1
    line([wall.x(i) wall.x(i+1)],[wall.y(i) wall.y(i+1)],'Color','k','LineWidth',2);
end
k = rsltn;
for i = 1:rsltn
    
    for j = 1:k
        
        if i == 1
            
            if j == k
                Cm_line = line([0 int_x(i,j)],[throat_rad int_y(i,j)],'Color','b');
                Cp_line = line([int_x(i,j) wall.x(i)],[int_y(i,j) wall.y(i)],'Color','r');
            else
                line([0 int_x(i,j)],[throat_rad int_y(i,j)],'Color','b');
                line([int_x(i,j) int_x(i,j+1)],[int_y(i,j) int_y(i,j+1)],'Color','r');
            end
            
        else
            
            if j == k
                line([int_x(i,j) int_x(i-1,j+1)],[int_y(i,j) int_y(i-1,j+1)],'Color','b');
                line([int_x(i,j) wall.x(i)],[int_y(i,j) wall.y(i)],'Color','r');
            else
                line([int_x(i,j) int_x(i-1,j+1)],[int_y(i,j) int_y(i-1,j+1)],'Color','b');
                line([int_x(i,j) int_x(i,j+1)],[int_y(i,j) int_y(i,j+1)],'Color','r');
            end
            
        end
    end
    k = k - 1;
    
end
lgnd = legend([sonic_line Cm_line Cp_line],'Sonic Line','C- Characteristics','C+ Characteristics');
lgnd.Location = 'northwest';

%%Save Data
save 'raw_crdnts.dat' crdnts -ascii;
