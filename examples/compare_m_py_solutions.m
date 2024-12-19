function compare_m_py_solutions(s,path_py_nc)
% Quick function to compare results between 
% s: the MATLAB solution (structure outputted from "run_fjord.m")
% path_py_nc: the Python solution (full path to .nc file outputted in Python implementation)
if nargin < 2
    path_py_nc = '/Users/mmeb1/fjordrpm_coupling/test_inputs/example2_intermediary/outputs_example2_intermediary.nc';
end
sp.t = ncread(path_py_nc,'time');
sp.t = sp.t-sp.t(1);
sp.T = ncread(path_py_nc,'T');
sp.S = ncread(path_py_nc,'S');
sp.Ts = ncread(path_py_nc,'Ts');
sp.Ss = ncread(path_py_nc,'Ss');
sp.QVs = ncread(path_py_nc,'QVs');
sp.QVp = ncread(path_py_nc,'QVp');
sp.QVk = ncread(path_py_nc,'QVk');
sp.QVi = ncread(path_py_nc,'QVi');
sp.QVv = ncread(path_py_nc,'QVv');
sp.z = ncread(path_py_nc,'depth');

%% add mixing and iceberg melting as third row
figure("Position",[40 40 1500 700]); 
subplot(2,5,1); 
hold on; box on;
plot(s.t,s.T(1,:),'linewidth',2); 
plot(sp.t,sp.T(:,1),'linewidth',2);
xlabel('Time')
ylabel('Surf. temp')
legend('Matlab','Python','location','best')
set(gca,'fontsize',14)
subplot(2,5,2); 
hold on; box on;
plot(s.T(:,1),s.z,'linestyle',':','linewidth',2); 
plot(sp.T(1,:),-sp.z,'linestyle',':','linewidth',2);
plot(s.T(:,end),s.z,'linewidth',2); 
plot(sp.T(end,:),-sp.z,'linewidth',2);
ylabel('Depth')
xlabel('Temperature')
legend('','','Matlab','Python','location','best')
set(gca,'fontsize',14)
subplot(2,5,3); 
hold on; box on;
plot(s.S(:,1),s.z,'linestyle',':','linewidth',2); 
plot(sp.S(1,:),-sp.z,'linestyle',':','linewidth',2);
plot(s.S(:,end),s.z,'linewidth',2); 
plot(sp.S(end,:),-sp.z,'linewidth',2);
ylabel('Depth')
xlabel('Salinity')
legend('','','Matlab','Python','location','best')
set(gca,'fontsize',14)
subplot(2,5,4); 
hold on; box on;
plot(s.t,s.Ts,'linewidth',2); 
plot(sp.t,sp.Ts,'linestyle',':','linewidth',2);
xlabel('Time')
ylabel('Matlab Shelf Temp')
set(gca,'fontsize',14)
subplot(2,5,5); 
hold on; box on;
% plot(sp.t,s.Ts-sp.Ts','linewidth',2);
xlabel('Time')
ylabel('Shelf Temp diff.')
set(gca,'fontsize',14)
subplot(2,5,6); 
hold on;  box on;
plot(s.QVs(:,1),s.z,'linestyle',':','linewidth',2); 
plot(sp.QVs(1,:),-sp.z,'linestyle',':','linewidth',2);
plot(s.QVs(:,end),s.z,'linewidth',2); 
plot(sp.QVs(end,:),-sp.z,'linewidth',2);
ylabel('Depth')
xlabel('Shelf vol. fluxes')
legend('','','Matlab','Python','location','best')
set(gca,'fontsize',14)
subplot(2,5,7); 
hold on; box on;
plot(squeeze(s.QVp(:,1)),s.z,'linestyle',':','linewidth',2); 
plot(squeeze(sp.QVp(1,:)),-sp.z,'linestyle',':','linewidth',2);
plot(squeeze(s.QVp(:,end)),s.z,'linewidth',2); 
plot(squeeze(sp.QVp(end,:)),-sp.z,'linewidth',2);
ylabel('Depth')
xlabel('Plume vol. fluxes')
legend('','','Matlab','Python','location','best')
set(gca,'fontsize',14)
subplot(2,5,8); 
hold on; box on;
plot(squeeze(s.QVk(:,1)),s.z,'linestyle',':','linewidth',2); 
plot(squeeze(sp.QVk(1,:)),-sp.z,'linestyle',':','linewidth',2);
plot(squeeze(s.QVk(:,end)),s.z,'linewidth',2); 
plot(squeeze(sp.QVk(end,:)),-sp.z,'linewidth',2);
ylabel('Depth')
xlabel('Mixing vol. fluxes')
legend('','','Matlab','Python','location','best')
set(gca,'fontsize',14)
subplot(2,5,9); 
hold on; box on;
plot(squeeze(s.QVi(:,1)),s.z,'linestyle',':','linewidth',2); 
plot(squeeze(sp.QVi(1,:)),-sp.z,'linestyle',':','linewidth',2);
plot(squeeze(s.QVi(:,end)),s.z,'linewidth',2); 
plot(squeeze(sp.QVi(end,:)),-sp.z,'linewidth',2);
ylabel('Depth')
xlabel('Iceberg vol. fluxes')
legend('','','Matlab','Python','location','best')
set(gca,'fontsize',14)
subplot(2,5,10); 
hold on; box on;
plot(squeeze(s.QVv(:,1)),s.z,'linestyle',':','linewidth',2); 
plot(squeeze(sp.QVv(1,:)),-sp.z,'linestyle',':','linewidth',2);
plot(squeeze(s.QVv(:,end)),s.z,'linewidth',2); 
plot(squeeze(sp.QVv(end,:)),-sp.z,'linewidth',2);
ylabel('Depth')
xlabel('Vertical vol. fluxes')
legend('','','Matlab','Python','location','best')
set(gca,'fontsize',14)