clear
close all
N_s = 1000;
N = 10000;
generate_snaps(N_s,N);
eigfig = figure;
modefig = figure;


% 1 Shot answer
[Lam_min,eigenvalues] = mPOD(N_s,'snapshot',1,1e-8,N_s,N_s)
eigenvalues = real(eigenvalues);
figure(eigfig);
semilogy(eigenvalues)
hold on

figure(modefig);
for i=1:5
    data = myload(['mode1.' num2str(i)]);
    plot(data,'k')
    hold on
end
delete mode*

% Case that doesn't hit maxsize
[Lam_min,eigenvalues] = mPOD(1000,'snapshot',1,1e-8,N_s/10,N_s)
figure(eigfig);
eigenvalues = real(eigenvalues);
hold on
semilogy(eigenvalues)

figure(modefig);
for i=1:5
    data = myload(['mode1.' num2str(i)]);
    plot(data,'r')
    hold on
end
delete mode*


% Case that hits maxsize limit
[Lam_min,eigenvalues] = mPOD(1000,'snapshot',1,1e-8,5,10)
figure(eigfig);
eigenvalues = real(eigenvalues);
semilogy(eigenvalues)

figure(modefig);
for i=1:5
    data = myload(['mode1.' num2str(i)]);
    plot(data,'g')
    hold on
end
delete mode*


