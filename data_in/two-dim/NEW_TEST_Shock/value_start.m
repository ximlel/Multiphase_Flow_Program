line=10;
column=200;
shock=95;

%%
gamma=1.4;
M=6;

f=1/(2/(gamma+1)/M/M+(gamma-1)/(gamma+1));
g=2*gamma/(gamma+1)*M*M-(gamma-1)/(gamma+1);


rho_L=1
qt_L=1
p_L=1/gamma/M/M
rho_R=f
qt_R=1/f
p_R=g/gamma/M/M
%%


u_R=qt_R*sqrt(2.0)/2.0
v_R=qt_R*sqrt(2.0)/2.0

u_L=qt_L*sqrt(2.0)/2.0
v_L=qt_L*sqrt(2.0)/2.0

rho_M=(rho_R+rho_L)/2.0
u_M=u_R*rho_R/(rho_R+rho_L)+u_L*rho_L/(rho_R+rho_L)
v_M=v_R*rho_R/(rho_R+rho_L)+v_L*rho_L/(rho_R+rho_L)
p_M=p_R

rho=zeros(column,1);
fid = fopen('RHO.txt','wt');
for j=1:line
for i=1:(line+shock-j)
    rho(i)=rho_L;
end
rho(line+shock-j+1)=rho_M;
for i=(line+shock-j+2):column
    rho(i)=rho_R;
end
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
fid = fopen('U.txt','wt');
for j=1:line
for i=1:(line+shock-j)
    u(i)=u_L;
end
u(line+shock-j+1)=u_M;
for i=(line+shock-j+2):column
    u(i)=u_R;
end
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
end
fclose(fid);

v=zeros(column,1);
fid = fopen('V.txt','wt');
for j=1:line
for i=1:(line+shock-j)
    v(i)=v_L;
end
v(line+shock-j+1)=v_M;
for i=(line+shock-j+2):column
    v(i)=v_R;
end
fprintf(fid,'%12.10f\t',v);
fprintf(fid,'\n');
end
fclose(fid);

p=zeros(column,1);
fid = fopen('P.txt','wt');
for j=1:line
for i=1:(line+shock-j)
    p(i)=p_L;
end
p(line+shock-j+1)=p_M;
for i=(line+shock-j+2):column
    p(i)=p_R;
end
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
end
fclose(fid);
