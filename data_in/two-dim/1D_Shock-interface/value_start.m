line=1;
column=1600;
interface=800;
shock=160;

rho_1=1
u_1=0
p_1=1
phi_1=0
rho_2=1.9
u_2=0
p_2=1
phi_2=1
rho_3=2.7647
u_3=1.4833
p_3=4.4468
phi_3=0

rho=zeros(column,1);
for i=1:shock
    rho(i)=rho_3;
end
for i=(shock+1):interface
    rho(i)=rho_1;
end
for i=(interface+1):column
    rho(i)=rho_2;
end
fid = fopen('RHO.txt','wt');
for j=1:line
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
for i=1:shock
    u(i)=u_3;
end
for i=(shock+1):interface
    u(i)=u_1;
end
for i=(interface+1):column
    u(i)=u_2;
end
fid = fopen('U.txt','wt');
for j=1:line
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
end
fclose(fid);

v=zeros(column,1);
fid = fopen('V.txt','wt');
for j=1:line
fprintf(fid,'%12.10f\t',v);
fprintf(fid,'\n');
end
fclose(fid);

p=zeros(column,1);
for i=1:shock
    p(i)=p_3;
end
for i=(shock+1):interface
    p(i)=p_1;
end
for i=(interface+1):column
    p(i)=p_2;
end
fid = fopen('P.txt','wt');
for j=1:line
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
end
fclose(fid);

phi=zeros(column,1);
for i=1:shock
    phi(i)=phi_3;
end
for i=(shock+1):interface
    phi(i)=phi_1;
end
for i=(interface+1):column
    phi(i)=phi_2;
end
fid = fopen('PHI.txt','wt');
for j=1:line
fprintf(fid,'%12.10f\t',phi);
fprintf(fid,'\n');
end
fclose(fid);
