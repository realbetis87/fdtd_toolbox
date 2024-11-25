function result = FDTD_3D_function_TFSF(domain, materials, surfaces, source, sensor)

e0=8.85e-12; % vacuum's permittivity
m0=4*pi*1e-7; % vacuum's permeability
c0=1/sqrt(e0*m0); % speed of light
f=(domain.fmin+domain.fmax)/2; % reference frequency
BW=domain.fmax-domain.fmin; % bandwidth
ta=3./(pi*BW);
lambda=c0/f; % reference wavelength

% spatial parameters
dx = lambda/domain.wl_ratio;
dy = lambda/domain.wl_ratio;
dz = lambda/domain.wl_ratio;
N=round(domain.length(1)/dx);
M=round(domain.length(2)/dy);
K=round(domain.length(3)/dz);


[X,Y,Z] = meshgrid(linspace(0,1,size(domain.voxel_model,2)),linspace(0,1,size(domain.voxel_model,1)),linspace(0,1,size(domain.voxel_model,3)));
[Xq,Yq,Zq] = meshgrid(linspace(0,1,M),linspace(0,1,N),linspace(0,1,K));

voxel_box = interp3(X,Y,Z,domain.voxel_model,Xq,Yq,Zq,'nearest'); clear X Y Z Xq Yq Zq

% Source
ExSource = [];
EySource = [];
EzSource = [];


% round(source.EzLocation(1,1)*N),round(source.EzLocation(1,2)*M),round(source.EzLocation(1,3)*K)

for i=1:size(source.ExLocation,1)
    ExSource(i) = sub2ind([N M K],round(source.ExLocation(i,1)*N),round(source.ExLocation(i,2)*M),round(source.ExLocation(i,3)*K));
end
for i=1:size(source.EyLocation,1)
    EySource(i) = sub2ind([N M K],round(source.EyLocation(i,1)*N),round(source.EyLocation(i,2)*M),round(source.EyLocation(i,3)*K));
end
for i=1:size(source.EzLocation,1)
    EzSource(i) = sub2ind([N M K],round(source.EzLocation(i,1)*N),round(source.EzLocation(i,2)*M),round(source.EzLocation(i,3)*K));
end


%time parameters
Ttotal=domain.pulse_no*6*ta;
dt=0.99/(c0*sqrt(1/(dx^2)+1/(dy^2)+1/(dz^2)));
T=round(Ttotal/dt);

% signal
n=0:T;
Js = exp(-((dt * n - 3 * ta) / ta) .* ((dt * n - 3 * ta) / ta)) .* sin(2 * pi * f * (dt * n - 3 * ta));

% PML parameters
pmlc=8;
r=1e-3;
nn=3.;
smax=-log(r)*(nn+1)*e0*c0/2/dx/pmlc;


for i=0:pmlc-1
    e_s1(i+1)=(smax*((i+.0)^nn)/(pmlc-1)^nn); % electric field PML coefficient for small indices
    h_s1(i+1)=(smax*((i+.5)^nn)/(pmlc-1)^nn*m0/e0); % magnetic field PML coefficient for small indices
    
    e_s2(i+1)=(smax*((i+.5)^nn)/(pmlc-1)^nn); % electric field PML coefficient for large indices
    h_s2(i+1)=(smax*((i+.0)^nn)/(pmlc-1)^nn*m0/e0); % magnetic field PML coefficient for large indices
end

pHyl=zeros(pmlc,M,K);
pHzl=zeros(pmlc,M,K);
pEyl=zeros(pmlc,M,K);
pEzl=zeros(pmlc,M,K);

pHyr=zeros(pmlc,M,K);
pHzr=zeros(pmlc,M,K);
pEyr=zeros(pmlc,M,K);
pEzr=zeros(pmlc,M,K);

pHxf=zeros(N,pmlc,K);
pHzf=zeros(N,pmlc,K);
pExf=zeros(N,pmlc,K);
pEzf=zeros(N,pmlc,K);

pHxc=zeros(N,pmlc,K);
pHzc=zeros(N,pmlc,K);
pExc=zeros(N,pmlc,K);
pEzc=zeros(N,pmlc,K);

pHxt=zeros(N,M,pmlc);
pHyt=zeros(N,M,pmlc);
pExt=zeros(N,M,pmlc);
pEyt=zeros(N,M,pmlc);

pHxb=zeros(N,M,pmlc);
pHyb=zeros(N,M,pmlc);
pExb=zeros(N,M,pmlc);
pEyb=zeros(N,M,pmlc);



% material initialization
epsilon = zeros(N,M,K);
mu = zeros(N,M,K);
sigma = zeros(N,M,K);
Ca = zeros(N,M,K);
Cb = zeros(N,M,K);
Da = zeros(N,M,K);
Db = zeros(N,M,K);

nn=1;

for i=1:N
    for j=1:M
        for k=1:K
            if strcmp(materials.material(voxel_box(i,j,k)).type,'Standard')
                epsilon(i,j,k)=e0*materials.material(voxel_box(i,j,k)).er;
                mu(i,j,k)=m0*materials.material(voxel_box(i,j,k)).mr;
                sigma(i,j,k)=materials.material(voxel_box(i,j,k)).sigma;
                Ca(i,j,k)=(1.-sigma(i,j,k)*dt/2./epsilon(i,j,k))./(1.+sigma(i,j,k)*dt/2./epsilon(i,j,k));
                Cb(i,j,k)=(dt./epsilon(i,j,k))./(1.+sigma(i,j,k)*dt/2./epsilon(i,j,k));
                Da(i,j,k)=1.;
                Db(i,j,k)=(dt./mu(i,j,k));
            end
            if strcmp(materials.material(voxel_box(i,j,k)).type,'Lorentz')
                for ii=1:materials.material(voxel_box(i,j,k)).poles
                    % Lorentz material
                    e_s = materials.material(voxel_box(i,j,k)).epsilon_s(ii);
                    e_inf = materials.material(voxel_box(i,j,k)).epsilon_inf(ii);
                    omega_0 = materials.material(voxel_box(i,j,k)).omega_0(ii);
                    delta = materials.material(voxel_box(i,j,k)).delta(ii);
                    Gamma_0(nn,ii) = (e_s-e_inf)*omega_0^2/(sqrt(omega_0^2-delta^2));
                    gamma_0(nn,ii) = delta-1i*sqrt(omega_0^2-delta^2);
                end
                Ca(i,j,k)=(2*e0*e_inf-e0*real(sum(Gamma_0(nn,:)))*dt)/(2*e0*e_inf+e0*real(sum(Gamma_0(nn,:)))*dt);
                Cb(i,j,k)=2*dt/(2*e0*e_inf+e0*real(sum(Gamma_0(nn,:)))*dt);
                Da(i,j,k)=1.;
                Db(i,j,k)=(dt./m0);
                id_mat(nn,1)=i;id_mat(nn,2)=j;id_mat(nn,3)=k; nn=nn+1;
            end
            if strcmp(materials.material(voxel_box(i,j,k)).type,'Debye')
                for ii=1:materials.material(voxel_box(i,j,k)).poles
                    % Debye material
                    e_s = materials.material(voxel_box(i,j,k)).epsilon_s(ii);
                    e_inf = materials.material(voxel_box(i,j,k)).epsilon_inf(ii);
                    tau_0 = materials.material(voxel_box(i,j,k)).tau_0(ii);
                    Gamma_0(nn,ii) = (e_s-e_inf)/tau_0;
                    gamma_0(nn,ii) = 1/tau_0;
                end
                Ca(i,j,k)=(2*e0*e_inf-e0*real(sum(Gamma_0(nn,:)))*dt)/(2*e0*e_inf+e0*real(sum(Gamma_0(nn,:)))*dt);
                Cb(i,j,k)=2*dt/(2*e0*e_inf+e0*real(sum(Gamma_0(nn,:)))*dt);
                Da(i,j,k)=1.;
                Db(i,j,k)=(dt./m0);
                id_mat(nn,1)=i;id_mat(nn,2)=j;id_mat(nn,3)=k; nn=nn+1;
            end
        end
    end
end


psi_x = zeros(nn-1,10);
psi_y = zeros(nn-1,10);
psi_z = zeros(nn-1,10);

% Planar material initialization

nx=1;ny=1;nz=1;
for ii=1:surfaces.total_surfaces
    
    if surfaces.surface(ii).orientation == 'x'
        for j=round(surfaces.surface(ii).ydim(1)*M+1):round(surfaces.surface(ii).ydim(2)*M)
            for k=round(surfaces.surface(ii).zdim(1)*K+1):round(surfaces.surface(ii).zdim(2)*K)
                for jj=1:surfaces.surface(ii).poles
                    if isreal(surfaces.surface(ii).amplitude(jj))
                        amc_x(nx,jj) = surfaces.surface(ii).amplitude(jj)/2;
                        gamma_x(nx,jj) = surfaces.surface(ii).gamma(jj)/2;
                    else
                        amc_x(nx,jj) = surfaces.surface(ii).amplitude(jj);
                        gamma_x(nx,jj) = surfaces.surface(ii).gamma(jj);
                    end
                    if surfaces.surface(ii).time_varying
                        surface_tv_x(nx) = true;
                        mod_strength_x(nx) = surfaces.surface(ii).modulation_strength;
                        mod_freq_x(nx) = surfaces.surface(ii).modulation_frequency;
                        amcon_x(nx) = surfaces.surface(ii).amplitude(jj)/2;
                    else
                        surface_tv_x(nx) = false;
                    end
                    
                end
                id_surx(nx,1)=round(surfaces.surface(ii).xdim*N);id_surx(nx,2)=j;id_surx(nx,3)=k; nx=nx+1;
            end
        end
    end
    
    if surfaces.surface(ii).orientation == 'y'
        for i=round(surfaces.surface(ii).xdim(1)*N+1):round(surfaces.surface(ii).xdim(2)*N)
            for k=round(surfaces.surface(ii).zdim(1)*K+1):round(surfaces.surface(ii).zdim(2)*K)
                for jj=1:surfaces.surface(ii).poles
                    if isreal(surfaces.surface(ii).amplitude(jj))
                        amc_y(ny,jj) = surfaces.surface(ii).amplitude(jj)/2;
                        gamma_y(ny,jj) = surfaces.surface(ii).gamma(jj)/2;
                    else
                        amc_y(ny,jj) = surfaces.surface(ii).amplitude(jj);
                        gamma_y(ny,jj) = surfaces.surface(ii).gamma(jj);
                    end
                    if surfaces.surface(ii).time_varying
                        surface_tv_y(ny) = true;
                        mod_strength_y(ny) = surfaces.surface(ii).modulation_strength;
                        mod_freq_y(ny) = surfaces.surface(ii).modulation_frequency;
                        amcon_y(ny) = surfaces.surface(ii).amplitude(jj)/2;
                    else
                        surface_tv_y(ny) = false;
                    end
                end
                id_sury(ny,1)=i;id_sury(ny,2)=round(surfaces.surface(ii).ydim*N);id_sury(ny,3)=k; ny=ny+1;
            end
        end
    end
    
    if surfaces.surface(ii).orientation == 'z'
        for i=round(surfaces.surface(ii).xdim(1)*N+1):round(surfaces.surface(ii).xdim(2)*N)
            for j=round(surfaces.surface(ii).ydim(1)*M+1):round(surfaces.surface(ii).ydim(2)*M)
                for jj=1:surfaces.surface(ii).poles
                    if isreal(surfaces.surface(ii).amplitude(jj))
                        amc_z(nz,jj) = surfaces.surface(ii).amplitude(jj)/2;
                        gamma_z(nz,jj) = surfaces.surface(ii).gamma(jj)/2;
                    else
                        amc_z(nz,jj) = surfaces.surface(ii).amplitude(jj);
                        gamma_z(nz,jj) = surfaces.surface(ii).gamma(jj);
                    end
                    if surfaces.surface(ii).time_varying
                        surface_tv_z(nz) = true;
                        mod_strength_z(nz) = surfaces.surface(ii).modulation_strength;
                        mod_freq_z(nz) = surfaces.surface(ii).modulation_frequency;
                        amcon_z(nz) = surfaces.surface(ii).amplitude(jj)/2;
                    else
                        surface_tv_z(nz) = false;
                    end
                end
                id_surz(nz,1)=i;id_surz(nz,2)=j;id_surz(nz,3)=round(surfaces.surface(ii).zdim*K); nz=nz+1;
            end
        end
    end
    
end

Jy_x = zeros(nx-1,2);
Jz_x = zeros(nx-1,2);

Jx_y = zeros(ny-1,2);
Jz_y = zeros(ny-1,2);

Jx_z = zeros(nz-1,2);
Jy_z = zeros(nz-1,2);


% Electromagnetic field parameter initialization

Ex=zeros(N,M,K);
Ey=zeros(N,M,K);
Ez=zeros(N,M,K);
Hx=zeros(N,M,K);
Hy=zeros(N,M,K);
Hz=zeros(N,M,K);


% Sensor
timeSensorLocation = [];

for i=1:size(sensor.timeLocation,1)
    timeSensorLocation(i) = sub2ind([N M K],round(sensor.timeLocation(i,1)*N),round(sensor.timeLocation(i,2)*M),round(sensor.timeLocation(i,3)*K));
end

timeSensor = zeros(T,6*size(sensor.timeLocation,1)+2);

fEx=zeros(N,M,K,length(sensor.frequencyElectric));
fEy=zeros(N,M,K,length(sensor.frequencyElectric));
fEz=zeros(N,M,K,length(sensor.frequencyElectric));
fHx=zeros(N,M,K,length(sensor.frequencyMagnetic));
fHy=zeros(N,M,K,length(sensor.frequencyMagnetic));
fHz=zeros(N,M,K,length(sensor.frequencyMagnetic));

freq_normalizationE = zeros(length(sensor.frequencyElectric),1);
freq_normalizationH = zeros(length(sensor.frequencyMagnetic),1);


% TF/SF parameters
if source.TFSF
    tot_sc=3;
    m_0=2;
    psi=0*pi/180;
    phi=source.phi*pi/180;
    theta=source.theta*pi/180;
    kx=sin(theta)*cos(phi);
    ky=sin(theta)*sin(phi);
    kz=cos(theta);
    N1=N-2-2*pmlc-2*tot_sc;
    M1=M-2-2*pmlc-2*tot_sc;
    K1=K-2-2*pmlc-2*tot_sc;
    NMK=2*(N1+M1+K1);

    % phase velocity correction
    tsA=dx*cos(phi)*sin(theta)/2;
    tsB=dx*sin(phi)*sin(theta)/2;
    tsD=dx*cos(theta)/2;
    tsC=(dx/c0/dt)^2*(sin(2*pi*f*dt/2))^2;
    ck=2*pi*f/c0;
    for i=0:10
        ck=ck-((sin(tsA*ck)^2)+(sin(tsB*ck)^2)+(sin(tsD*ck)^2)-tsC)/(tsA*sin(2*tsA*ck)+tsB*sin(2*tsB*ck)+tsD*sin(2*tsD*ck));
    end
    tsA=dx*cos(0)*sin(0)/2;
    tsB=dx*sin(0)*sin(0)/2;
    tsD=dx*cos(0)/2;
    tsC=(dx/c0/dt)^2*(sin(2*pi*f*dt/2))^2;
    ck0=2*pi*f/c0;
    for i=0:10
        ck0=ck0-((sin(tsA*ck0)^2)+(sin(tsB*ck0)^2)+(sin(tsD*ck0)^2)-tsC)/(tsA*sin(2*tsA*ck0)+tsB*sin(2*tsB*ck0)+tsD*sin(2*tsD*ck0));
    end
    vp_corr=ck/ck0;

    Einc=zeros(NMK,1);
    Hinc=zeros(NMK,1);

    Exinc=zeros(NMK,1);
    Eyinc=zeros(NMK,1);
    Ezinc=zeros(NMK,1);

    Hxinc=zeros(NMK,1);
    Hyinc=zeros(NMK,1);
    Hzinc=zeros(NMK,1);


    Hzt=zeros(N1,K1);
    Hzb=zeros(N1,K1);
    Hxt=zeros(N1,K1);
    Hxb=zeros(N1,K1);
    Hyf=zeros(N1,M1);
    Hyc=zeros(N1,M1);
    Hxf=zeros(N1,M1);
    Hxc=zeros(N1,M1);
    Hzl=zeros(M1,K1);
    Hzr=zeros(M1,K1);
    Hyl=zeros(M1,K1);
    Hyr=zeros(M1,K1);

    Ezt=zeros(N1,K1);
    Ezb=zeros(N1,K1);
    Ext=zeros(N1,K1);
    Exb=zeros(N1,K1);
    Eyf=zeros(N1,M1);
    Eyc=zeros(N1,M1);
    Exf=zeros(N1,M1);
    Exc=zeros(N1,M1);
    Ezl=zeros(M1,K1);
    Ezr=zeros(M1,K1);
    Eyl=zeros(M1,K1);
    Eyr=zeros(M1,K1);
end

tic

disp('Simulation started')

for n=0:T
    if ~mod(n, round(T / 10))
        fprintf('Progress: %.0f%%\n', (n / T) * 100);
    end
    
    %Js = exp(-((dt * n - 3 * ta) / ta) * ((dt * n - 3 * ta) / ta)) * sin(2 * pi * f * (dt * n - 3 * ta));
    if source.TFSF
        Einc(0+1)=Js(n+1);
        einc=Einc(NMK-2+1);
        for i=1:NMK-2
            Einc(i+1)=Einc(i+1)+dt/(e0*dx*vp_corr)*(Hinc(i+1)-Hinc(i+2));
            Exinc(i+1)=Einc(i+1)*(cos(psi)*sin(phi)-sin(psi)*cos(theta)*cos(phi));
            Eyinc(i+1)=Einc(i+1)*(-cos(psi)*cos(phi)-sin(psi)*cos(theta)*sin(phi));
            Ezinc(i+1)=Einc(i+1)*(sin(psi)*sin(theta));
        end
        Einc(NMK-1+1)=einc+(c0*dt-dx)/(c0*dt+dx)*(Einc(NMK-2+1)-Einc(NMK-1+1));
        
        for i=1:NMK-1
            Hinc(i+1)=Hinc(i+1)+dt/(m0*dx*vp_corr)*(Einc(i-1+1)-Einc(i+1));
            Hxinc(i+1)=Hinc(i+1)*(sin(psi)*sin(phi)+cos(psi)*cos(theta)*cos(phi));
            Hyinc(i+1)=Hinc(i+1)*(-sin(psi)*cos(phi)+cos(psi)*cos(theta)*sin(phi));
            Hzinc(i+1)=Hinc(i+1)*(-cos(psi)*sin(theta));
        end
        
        
        for i=0:N1-1
            for j=0:M1-1
                rx=i*dx+0.5*dx;
                ry=j*dy;
                rz=0*dz; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Exf(i+1,j+1)=(1-dd)*Exinc(m_0+floor(dis)+1)+dd*Exinc(m_0+floor(dis)+2);
                rz=(K1-1)*dz; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Exc(i+1,j+1)=(1-dd)*Exinc(m_0+floor(dis)+1)+dd*Exinc(m_0+floor(dis)+2);
                
                rx=i*dx;
                ry=j*dy+0.5*dy;
                rz=0*dz; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Eyf(i+1,j+1)=(1-dd)*Eyinc(m_0+floor(dis)+1)+dd*Eyinc(m_0+floor(dis)+2);
                rz=(K1-1)*dz; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Eyc(i+1,j+1)=(1-dd)*Eyinc(m_0+floor(dis)+1)+dd*Eyinc(m_0+floor(dis)+2);
                
                rx=i*dx;
                ry=j*dy+0.5*dy;
                rz=(-0.5)*dz; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hxf(i+1,j+1)=(1-dd)*Hxinc(m_0+floor(dis)+1)+dd*Hxinc(m_0+floor(dis)+2);
                rz=(K1-0.5)*dz; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hxc(i+1,j+1)=(1-dd)*Hxinc(m_0+floor(dis)+1)+dd*Hxinc(m_0+floor(dis)+2);
                
                rx=i*dx+0.5*dx;
                ry=j*dy;
                rz=(-0.5)*dz; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hyf(i+1,j+1)=(1-dd)*Hyinc(m_0+floor(dis)+1)+dd*Hyinc(m_0+floor(dis)+2);
                rz=(K1-0.5)*dz; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hyc(i+1,j+1)=(1-dd)*Hyinc(m_0+floor(dis)+1)+dd*Hyinc(m_0+floor(dis)+2);
            end
        end
        
        for i=0:N1-1
            for j=0:K1-1
                rx=i*dx+0.5*dx;
                rz=j*dz;
                ry=0*dy; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Exb(i+1,j+1)=(1-dd)*Exinc(m_0+floor(dis)+1)+dd*Exinc(m_0+floor(dis)+2);
                ry=(M1-1)*dy; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Ext(i+1,j+1)=(1-dd)*Exinc(m_0+floor(dis)+1)+dd*Exinc(m_0+floor(dis)+2);
                
                rx=i*dx;
                rz=j*dz+0.5*dz;
                ry=0*dy; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Ezb(i+1,j+1)=(1-dd)*Ezinc(m_0+floor(dis)+1)+dd*Ezinc(m_0+floor(dis)+2);
                ry=(M1-1)*dy; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Ezt(i+1,j+1)=(1-dd)*Ezinc(m_0+floor(dis)+1)+dd*Ezinc(m_0+floor(dis)+2);
                
                rx=i*dx;
                rz=j*dz+0.5*dz;
                ry=(-0.5)*dy; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hxb(i+1,j+1)=(1-dd)*Hxinc(m_0+floor(dis)+1)+dd*Hxinc(m_0+floor(dis)+2);
                ry=(M1-0.5)*dy; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hxt(i+1,j+1)=(1-dd)*Hxinc(m_0+floor(dis)+1)+dd*Hxinc(m_0+floor(dis)+2);
                
                rx=i*dx+0.5*dx;
                rz=j*dz;
                ry=(-0.5)*dy; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hzb(i+1,j+1)=(1-dd)*Hzinc(m_0+floor(dis)+1)+dd*Hzinc(m_0+floor(dis)+2);
                ry=(M1-0.5)*dy; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hzt(i+1,j+1)=(1-dd)*Hzinc(m_0+floor(dis)+1)+dd*Hzinc(m_0+floor(dis)+2);
            end
        end
        
        
        for i=0:M1-1
            for j=0:K1-1
                ry=i*dy+0.5*dy;
                rz=j*dz;
                rx=0*dx; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Eyl(i+1,j+1)=(1-dd)*Eyinc(m_0+floor(dis)+1)+dd*Eyinc(m_0+floor(dis)+2);
                rx=(N1-1)*dx; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Eyr(i+1,j+1)=(1-dd)*Eyinc(m_0+floor(dis)+1)+dd*Eyinc(m_0+floor(dis)+2);
                
                ry=i*dy;
                rz=j*dz+0.5*dz;
                rx=0*dx; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Ezl(i+1,j+1)=(1-dd)*Ezinc(m_0+floor(dis)+1)+dd*Ezinc(m_0+floor(dis)+2);
                rx=(N1-1)*dx; dis=(rx*kx+ry*ky+rz*kz)/dx; dd=dis-floor(dis);
                Ezr(i+1,j+1)=(1-dd)*Ezinc(m_0+floor(dis)+1)+dd*Ezinc(m_0+floor(dis)+2);
                
                ry=i*dy;
                rz=j*dz+0.5*dz;
                rx=(-0.5)*dx; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hyl(i+1,j+1)=(1-dd)*Hyinc(m_0+floor(dis)+1)+dd*Hyinc(m_0+floor(dis)+2);
                rx=(N1-0.5)*dx; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hyr(i+1,j+1)=(1-dd)*Hyinc(m_0+floor(dis)+1)+dd*Hyinc(m_0+floor(dis)+2);
                
                ry=i*dy+0.5*dy;
                rz=j*dz;
                rx=(-0.5)*dx; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hzl(i+1,j+1)=(1-dd)*Hzinc(m_0+floor(dis)+1)+dd*Hzinc(m_0+floor(dis)+2);
                rx=(N1-0.5)*dx; dis=(rx*kx+ry*ky+rz*kz)/dx+0.5; dd=dis-floor(dis);
                Hzr(i+1,j+1)=(1-dd)*Hzinc(m_0+floor(dis)+1)+dd*Hzinc(m_0+floor(dis)+2);
            end
        end
    end
    
    
    for i=2:N
        for j=2:M
            for k=2:K
                Hx(i,j,k)=Da(i,j,k)*Hx(i,j,k)+Db(i,j,k)*((Ey(i,j,k)-Ey(i,j,k-1))/dz-(Ez(i,j,k)-Ez(i,j-1,k))/dy);
                Hy(i,j,k)=Da(i,j,k)*Hy(i,j,k)+Db(i,j,k)*((Ez(i,j,k)-Ez(i-1,j,k))/dx-(Ex(i,j,k)-Ex(i,j,k-1))/dz);
                Hz(i,j,k)=Da(i,j,k)*Hz(i,j,k)+Db(i,j,k)*((Ex(i,j,k)-Ex(i,j-1,k))/dy-(Ey(i,j,k)-Ey(i-1,j,k))/dx);
            end
        end
    end
    
    
    % TF/SF magnetic field
    if source.TFSF
        for i=1+pmlc+tot_sc:N-1-pmlc-tot_sc-1
            for j=2+pmlc+tot_sc:M-1-pmlc-tot_sc-1
                Hx(i+1,j+1,1+pmlc+tot_sc+1)=Hx(i+1,j+1,1+pmlc+tot_sc+1)-Db(i+1,j+1,1+pmlc+tot_sc+1)*Eyf(i-1-pmlc-tot_sc+1,j-2-pmlc-tot_sc+1);
                Hx(i+1,j+1,K-1-pmlc-tot_sc+1)=Hx(i+1,j+1,K-1-pmlc-tot_sc+1)+Db(i+1,j+1,K-1-pmlc-tot_sc+1)*Eyc(i-1-pmlc-tot_sc+1,j-2-pmlc-tot_sc+1);
            end
        end
        
        for i=2+pmlc+tot_sc:N-1-pmlc-tot_sc-1
            for j=1+pmlc+tot_sc:M-1-pmlc-tot_sc-1
                Hy(i+1,j+1,1+pmlc+tot_sc+1)=Hy(i+1,j+1,1+pmlc+tot_sc+1)+Db(i+1,j+1,1+pmlc+tot_sc+1)*Exf(i-2-pmlc-tot_sc+1,j-1-pmlc-tot_sc+1);
                Hy(i+1,j+1,K-1-pmlc-tot_sc+1)=Hy(i+1,j+1,K-1-pmlc-tot_sc+1)-Db(i+1,j+1,K-1-pmlc-tot_sc+1)*Exc(i-2-pmlc-tot_sc+1,j-1-pmlc-tot_sc+1);
            end
        end
        
        for j=2+pmlc+tot_sc:M-1-pmlc-tot_sc-1
            for k=1+pmlc+tot_sc:K-1-pmlc-tot_sc-1
                Hz(1+pmlc+tot_sc+1,j+1,k+1)=Hz(1+pmlc+tot_sc+1,j+1,k+1)+Db(1+pmlc+tot_sc+1,j+1,k+1)*Eyl(j-2-pmlc-tot_sc+1,k-1-pmlc-tot_sc+1);
                Hz(N-1-pmlc-tot_sc+1,j+1,k+1)=Hz(N-1-pmlc-tot_sc+1,j+1,k+1)-Db(N-1-pmlc-tot_sc+1,j+1,k+1)*Eyr(j-2-pmlc-tot_sc+1,k-1-pmlc-tot_sc+1);
            end
        end
        
        for j=1+pmlc+tot_sc:M-1-pmlc-tot_sc-1
            for k=2+pmlc+tot_sc:K-1-pmlc-tot_sc-1
                Hy(1+pmlc+tot_sc+1,j+1,k+1)=Hy(1+pmlc+tot_sc+1,j+1,k+1)-Db(1+pmlc+tot_sc+1,j+1,k+1)*Ezl(j-1-pmlc-tot_sc+1,k-2-pmlc-tot_sc+1);
                Hy(N-1-pmlc-tot_sc+1,j+1,k+1)=Hy(N-1-pmlc-tot_sc+1,j+1,k+1)+Db(N-1-pmlc-tot_sc+1,j+1,k+1)*Ezr(j-1-pmlc-tot_sc+1,k-2-pmlc-tot_sc+1);
            end
        end
        
        for i=1+pmlc+tot_sc:N-1-pmlc-tot_sc-1
            for k=2+pmlc+tot_sc:K-1-pmlc-tot_sc-1
                Hx(i+1,1+pmlc+tot_sc+1,k+1)=Hx(i+1,1+pmlc+tot_sc+1,k+1)+Db(i+1,1+pmlc+tot_sc+1,k+1)*Ezb(i-1-pmlc-tot_sc+1,k-2-pmlc-tot_sc+1);
                Hx(i+1,M-1-pmlc-tot_sc+1,k+1)=Hx(i+1,M-1-pmlc-tot_sc+1,k+1)-Db(i+1,M-1-pmlc-tot_sc+1,k+1)*Ezt(i-1-pmlc-tot_sc+1,k-2-pmlc-tot_sc+1);
            end
        end
        for i=2+pmlc+tot_sc:N-1-pmlc-tot_sc-1
            for k=1+pmlc+tot_sc:K-1-pmlc-tot_sc-1
                Hz(i+1,1+pmlc+tot_sc+1,k+1)=Hz(i+1,1+pmlc+tot_sc+1,k+1)-Db(i+1,1+pmlc+tot_sc+1,k+1)*Exb(i-2-pmlc-tot_sc+1,k-1-pmlc-tot_sc+1);
                Hz(i+1,M-1-pmlc-tot_sc+1,k+1)=Hz(i+1,M-1-pmlc-tot_sc+1,k+1)+Db(i+1,M-1-pmlc-tot_sc+1,k+1)*Ext(i-2-pmlc-tot_sc+1,k-1-pmlc-tot_sc+1);
            end
        end
    end
    
    for i=2:pmlc
        for j=2:M
            for k=2:K
                pHyl(i,j,k)=(Da(i,j,k)-dt/m0*h_s1(pmlc-i+1))*pHyl(i,j,k)+dt/m0*h_s1(pmlc-i+1)*Db(i,j,k)*(Ez(i,j,k)-Ez(i-1,j,k))/dx;
                pHzl(i,j,k)=(Da(i,j,k)-dt/m0*h_s1(pmlc-i+1))*pHzl(i,j,k)-dt/m0*h_s1(pmlc-i+1)*Db(i,j,k)*(Ey(i,j,k)-Ey(i-1,j,k))/dx;
                Hy(i,j,k)=Hy(i,j,k)-pHyl(i,j,k);
                Hz(i,j,k)=Hz(i,j,k)-pHzl(i,j,k);
            end
        end
    end
    
    for i=N-pmlc+1:N
        for j=2:M
            for k=2:K
                pHyr(i-N+pmlc,j,k)=(Da(i,j,k)-dt/m0*h_s2(i-N+pmlc))*pHyr(i-N+pmlc,j,k)+dt/m0*h_s2(i-N+pmlc)*Db(i,j,k)*(Ez(i,j,k)-Ez(i-1,j,k))/dx;
                pHzr(i-N+pmlc,j,k)=(Da(i,j,k)-dt/m0*h_s2(i-N+pmlc))*pHzr(i-N+pmlc,j,k)-dt/m0*h_s2(i-N+pmlc)*Db(i,j,k)*(Ey(i,j,k)-Ey(i-1,j,k))/dx;
                Hy(i,j,k)=Hy(i,j,k)-pHyr(i-N+pmlc,j,k);
                Hz(i,j,k)=Hz(i,j,k)-pHzr(i-N+pmlc,j,k);
            end
        end
    end
    
    for i=2:N
        for j=2:pmlc
            for k=2:K
                pHxc(i,j,k)=(Da(i,j,k)-dt/m0*h_s1(pmlc-j+1))*pHxc(i,j,k)-dt/m0*h_s1(pmlc-j+1)*Db(i,j,k)*(Ez(i,j,k)-Ez(i,j-1,k))/dy;
                pHzc(i,j,k)=(Da(i,j,k)-dt/m0*h_s1(pmlc-j+1))*pHzc(i,j,k)+dt/m0*h_s1(pmlc-j+1)*Db(i,j,k)*(Ex(i,j,k)-Ex(i,j-1,k))/dy;
                Hx(i,j,k)=Hx(i,j,k)-pHxc(i,j,k);
                Hz(i,j,k)=Hz(i,j,k)-pHzc(i,j,k);
            end
        end
    end
    
    for i=2:N
        for j=M-pmlc+1:M
            for k=2:K
                pHxf(i,j-M+pmlc,k)=(Da(i,j,k)-dt/m0*h_s2(j-M+pmlc))*pHxf(i,j-M+pmlc,k)-dt/m0*h_s2(j-M+pmlc)*Db(i,j,k)*(Ez(i,j,k)-Ez(i,j-1,k))/dy;
                pHzf(i,j-M+pmlc,k)=(Da(i,j,k)-dt/m0*h_s2(j-M+pmlc))*pHzf(i,j-M+pmlc,k)+dt/m0*h_s2(j-M+pmlc)*Db(i,j,k)*(Ex(i,j,k)-Ex(i,j-1,k))/dy;
                Hx(i,j,k)=Hx(i,j,k)-pHxf(i,j-M+pmlc,k);
                Hz(i,j,k)=Hz(i,j,k)-pHzf(i,j-M+pmlc,k);
            end
        end
    end
    
    for i=2:N
        for j=2:M
            for k=2:pmlc
                pHxb(i,j,k)=(Da(i,j,k)-dt/m0*h_s1(pmlc-k+1))*pHxb(i,j,k)+dt/m0*h_s1(pmlc-k+1)*Db(i,j,k)*(Ey(i,j,k)-Ey(i,j,k-1))/dz;
                pHyb(i,j,k)=(Da(i,j,k)-dt/m0*h_s1(pmlc-k+1))*pHyb(i,j,k)-dt/m0*h_s1(pmlc-k+1)*Db(i,j,k)*(Ex(i,j,k)-Ex(i,j,k-1))/dz;
                Hx(i,j,k)=Hx(i,j,k)-pHxb(i,j,k);
                Hy(i,j,k)=Hy(i,j,k)-pHyb(i,j,k);
            end
        end
    end
    
     for i=2:N
        for j=2:M
            for k=K-pmlc+1:K
                pHxt(i,j,k-K+pmlc)=(Da(i,j,k)-dt/m0*h_s2(k-K+pmlc))*pHxt(i,j,k-K+pmlc)+dt/m0*h_s2(k-K+pmlc)*Db(i,j,k)*(Ey(i,j,k)-Ey(i,j,k-1))/dz;
                pHyt(i,j,k-K+pmlc)=(Da(i,j,k)-dt/m0*h_s2(k-K+pmlc))*pHyt(i,j,k-K+pmlc)-dt/m0*h_s2(k-K+pmlc)*Db(i,j,k)*(Ex(i,j,k)-Ex(i,j,k-1))/dz;
                Hx(i,j,k)=Hx(i,j,k)-pHxt(i,j,k-K+pmlc);
                Hy(i,j,k)=Hy(i,j,k)-pHyt(i,j,k-K+pmlc);
            end
        end
     end
    
     
     for i=1:length(psi_x)
         for j=1:length(gamma_0(i,:))
             psi_x(i,j) = psi_x(i,j)*exp(-gamma_0(i,j)*dt)+Gamma_0(i,j)*(1-exp(-gamma_0(i,j)*dt))*Ex(id_mat(i,1),id_mat(i,2),id_mat(i,3));
             psi_y(i,j) = psi_y(i,j)*exp(-gamma_0(i,j)*dt)+Gamma_0(i,j)*(1-exp(-gamma_0(i,j)*dt))*Ey(id_mat(i,1),id_mat(i,2),id_mat(i,3));
             psi_z(i,j) = psi_z(i,j)*exp(-gamma_0(i,j)*dt)+Gamma_0(i,j)*(1-exp(-gamma_0(i,j)*dt))*Ez(id_mat(i,1),id_mat(i,2),id_mat(i,3));
         end
     end
     
     for i=1:length(Jy_x)
         if surface_tv_x(i)
             amc_x(i,1) = amcon_x(i)*(1+mod_strength_x(i)*cos(2*pi*mod_freq_x(i)*n*dt));
         end
         for j=1:length(gamma_x(i,:))
             Jy_x(i,j)=Jy_x(i,j)*exp(-gamma_x(i,j)*dt)+Ey(id_surx(i,1),id_surx(i,2),id_surx(i,3))*amc_x(i,j)*dt;
             Jz_x(i,j)=Jz_x(i,j)*exp(-gamma_x(i,j)*dt)+Ez(id_surx(i,1),id_surx(i,2),id_surx(i,3))*amc_x(i,j)*dt;
         end
     end
     
     for i=1:length(Jx_y)
         if surface_tv_y(i)
             amc_y(i,1) = amcon_y(i)*(1+mod_strength_y(i)*cos(2*pi*mod_freq_y(i)*n*dt));
         end
         for j=1:length(gamma_y(i,:))
             Jx_y(i,j)=Jx_y(i,j)*exp(-gamma_y(i,j)*dt)+Ex(id_sury(i,1),id_sury(i,2),id_sury(i,3))*amc_y(i,j)*dt;
             Jz_y(i,j)=Jz_y(i,j)*exp(-gamma_y(i,j)*dt)+Ez(id_sury(i,1),id_sury(i,2),id_sury(i,3))*amc_y(i,j)*dt;
         end
     end
     
     for i=1:length(Jx_z)
         if surface_tv_z(i)
             amc_z(i,1) = amcon_z(i)*(1+mod_strength_z(i)*cos(2*pi*mod_freq_z(i)*n*dt));
         end
         for j=1:length(gamma_z(i,:))
             Jx_z(i,j)=Jx_z(i,j)*exp(-gamma_z(i,j)*dt)+Ex(id_surz(i,1),id_surz(i,2),id_surz(i,3))*amc_z(i,j)*dt;
             Jy_z(i,j)=Jy_z(i,j)*exp(-gamma_z(i,j)*dt)+Ey(id_surz(i,1),id_surz(i,2),id_surz(i,3))*amc_z(i,j)*dt;
         end
     end
     
    
    for i=1:N-1
        for j=1:M-1
            for k=1:K-1
                Ex(i,j,k)=Ca(i,j,k)*Ex(i,j,k)+Cb(i,j,k)*((Hz(i,j+1,k)-Hz(i,j,k))/dy-(Hy(i,j,k+1)-Hy(i,j,k))/dz);
                Ey(i,j,k)=Ca(i,j,k)*Ey(i,j,k)+Cb(i,j,k)*((Hx(i,j,k+1)-Hx(i,j,k))/dz-(Hz(i+1,j,k)-Hz(i,j,k))/dx);
                Ez(i,j,k)=Ca(i,j,k)*Ez(i,j,k)+Cb(i,j,k)*((Hy(i+1,j,k)-Hy(i,j,k))/dx-(Hx(i,j+1,k)-Hx(i,j,k))/dy);
            end
        end
    end

    % TF/SF electric field
    if source.TFSF
        for i=2+pmlc+tot_sc:N-1-pmlc-tot_sc-1
            for j=1+pmlc+tot_sc:M-1-pmlc-tot_sc-1
                Ex(i+1,j+1,1+pmlc+tot_sc+1)=Ex(i+1,j+1,1+pmlc+tot_sc+1)+Cb(i+1,j+1,1+pmlc+tot_sc+1)*Hyf(i-2-pmlc-tot_sc+1,j-1-pmlc-tot_sc+1);
                Ex(i+1,j+1,K-2-pmlc-tot_sc+1)=Ex(i+1,j+1,K-2-pmlc-tot_sc+1)-Cb(i+1,j+1,K-2-pmlc-tot_sc+1)*Hyc(i-2-pmlc-tot_sc+1,j-1-pmlc-tot_sc+1);
            end
        end
        for i=1+pmlc+tot_sc:N-1-pmlc-tot_sc-1
            for j=2+pmlc+tot_sc:M-1-pmlc-tot_sc-1
                Ey(i+1,j+1,1+pmlc+tot_sc+1)=Ey(i+1,j+1,1+pmlc+tot_sc+1)-Cb(i+1,j+1,1+pmlc+tot_sc+1)*Hxf(i-1-pmlc-tot_sc+1,j-2-pmlc-tot_sc+1);
                Ey(i+1,j+1,K-2-pmlc-tot_sc+1)=Ey(i+1,j+1,K-2-pmlc-tot_sc+1)+Cb(i+1,j+1,K-2-pmlc-tot_sc+1)*Hxc(i-1-pmlc-tot_sc+1,j-2-pmlc-tot_sc+1);
            end
        end
        
        for j=2+pmlc+tot_sc:M-1-pmlc-tot_sc-1
            for k=1+pmlc+tot_sc:K-1-pmlc-tot_sc-1
                Ey(1+pmlc+tot_sc+1,j+1,k+1)=Ey(1+pmlc+tot_sc+1,j+1,k+1)+Cb(1+pmlc+tot_sc+1,j+1,k+1)*Hzl(j-2-pmlc-tot_sc+1,k-1-pmlc-tot_sc+1);
                Ey(N-2-pmlc-tot_sc+1,j+1,k+1)=Ey(N-2-pmlc-tot_sc+1,j+1,k+1)-Cb(N-2-pmlc-tot_sc+1,j+1,k+1)*Hzr(j-2-pmlc-tot_sc+1,k-1-pmlc-tot_sc+1);
            end
        end
        for j=1+pmlc+tot_sc:M-1-pmlc-tot_sc-1
            for k=2+pmlc+tot_sc:K-1-pmlc-tot_sc-1
                Ez(1+pmlc+tot_sc+1,j+1,k+1)=Ez(1+pmlc+tot_sc+1,j+1,k+1)-Cb(1+pmlc+tot_sc+1,j+1,k+1)*Hyl(j-1-pmlc-tot_sc+1,k-2-pmlc-tot_sc+1);
                Ez(N-2-pmlc-tot_sc+1,j+1,k+1)=Ez(N-2-pmlc-tot_sc+1,j+1,k+1)+Cb(N-2-pmlc-tot_sc+1,j+1,k+1)*Hyr(j-1-pmlc-tot_sc+1,k-2-pmlc-tot_sc+1);
            end
        end
        
        for i=2+pmlc+tot_sc:N-1-pmlc-tot_sc-1
            for k=1+pmlc+tot_sc:K-1-pmlc-tot_sc-1
                Ex(i+1,1+pmlc+tot_sc+1,k+1)=Ex(i+1,1+pmlc+tot_sc+1,k+1)-Cb(i+1,1+pmlc+tot_sc+1,k+1)*Hzb(i-2-pmlc-tot_sc+1,k-1-pmlc-tot_sc+1);
                Ex(i+1,M-2-pmlc-tot_sc+1,k+1)=Ex(i+1,M-2-pmlc-tot_sc+1,k+1)+Cb(i+1,M-2-pmlc-tot_sc+1,k+1)*Hzt(i-2-pmlc-tot_sc+1,k-1-pmlc-tot_sc+1);
            end
        end
        for i=1+pmlc+tot_sc:N-1-pmlc-tot_sc-1
            for k=2+pmlc+tot_sc:K-1-pmlc-tot_sc-1
                Ez(i+1,1+pmlc+tot_sc+1,k+1)=Ez(i+1,1+pmlc+tot_sc+1,k+1)+Cb(i+1,1+pmlc+tot_sc+1,k+1)*Hxb(i-1-pmlc-tot_sc+1,k-2-pmlc-tot_sc+1);
                Ez(i+1,M-2-pmlc-tot_sc+1,k+1)=Ez(i+1,M-2-pmlc-tot_sc+1,k+1)-Cb(i+1,M-2-pmlc-tot_sc+1,k+1)*Hxt(i-1-pmlc-tot_sc+1,k-2-pmlc-tot_sc+1);
            end
        end
    end
    
    for i=1:length(psi_x)
        Ex(id_mat(i,1),id_mat(i,2),id_mat(i,3)) = Ex(id_mat(i,1),id_mat(i,2),id_mat(i,3)) +Cb(id_mat(i,1),id_mat(i,2),id_mat(i,3))*e0*sum(imag(psi_x(i,:)));
        Ey(id_mat(i,1),id_mat(i,2),id_mat(i,3)) = Ey(id_mat(i,1),id_mat(i,2),id_mat(i,3)) +Cb(id_mat(i,1),id_mat(i,2),id_mat(i,3))*e0*sum(imag(psi_y(i,:)));
        Ez(id_mat(i,1),id_mat(i,2),id_mat(i,3)) = Ez(id_mat(i,1),id_mat(i,2),id_mat(i,3)) +Cb(id_mat(i,1),id_mat(i,2),id_mat(i,3))*e0*sum(imag(psi_z(i,:)));
    end
    
    for i=1:pmlc
        for j=1:M-1
            for k=1:K-1
                pEyl(i,j,k)=(Ca(i,j,k)-dt/e0*e_s1(pmlc-i+1))*pEyl(i,j,k)-dt/e0*e_s1(pmlc-i+1)*Cb(i,j,k)*(Hz(i+1,j,k)-Hz(i,j,k))/dx;
                pEzl(i,j,k)=(Ca(i,j,k)-dt/e0*e_s1(pmlc-i+1))*pEzl(i,j,k)+dt/e0*e_s1(pmlc-i+1)*Cb(i,j,k)*(Hy(i+1,j,k)-Hy(i,j,k))/dx;
                Ey(i,j,k)=Ey(i,j,k)-pEyl(i,j,k);
                Ez(i,j,k)=Ez(i,j,k)-pEzl(i,j,k);
            end
        end
    end
    
    for i=N-pmlc+1:N-1
        for j=1:M-1
            for k=1:K-1
                pEyr(i-N+pmlc,j,k)=(Ca(i,j,k)-dt/e0*e_s2(i-N+pmlc))*pEyr(i-N+pmlc,j,k)-dt/e0*e_s2(i-N+pmlc)*Cb(i,j,k)*(Hz(i+1,j,k)-Hz(i,j,k))/dx;
                pEzr(i-N+pmlc,j,k)=(Ca(i,j,k)-dt/e0*e_s2(i-N+pmlc))*pEzr(i-N+pmlc,j,k)+dt/e0*e_s2(i-N+pmlc)*Cb(i,j,k)*(Hy(i+1,j,k)-Hy(i,j,k))/dx;
                Ey(i,j,k)=Ey(i,j,k)-pEyr(i-N+pmlc,j,k);
                Ez(i,j,k)=Ez(i,j,k)-pEzr(i-N+pmlc,j,k);
            end
        end
    end
    
    for i=1:N-1
        for j=1:pmlc
            for k=1:K-1
                pExc(i,j,k)=(Ca(i,j,k)-dt/e0*e_s1(pmlc-j+1))*pExc(i,j,k)+dt/e0*e_s1(pmlc-j+1)*Cb(i,j,k)*(Hz(i,j+1,k)-Hz(i,j,k))/dy;
                pEzc(i,j,k)=(Ca(i,j,k)-dt/e0*e_s1(pmlc-j+1))*pEzc(i,j,k)-dt/e0*e_s1(pmlc-j+1)*Cb(i,j,k)*(Hx(i,j+1,k)-Hx(i,j,k))/dy;
                Ex(i,j,k)=Ex(i,j,k)-pExc(i,j,k);
                Ez(i,j,k)=Ez(i,j,k)-pEzc(i,j,k);
            end
        end
    end
    
    for i=1:N-1
        for j=M-pmlc+1:M-1
            for k=1:K-1
                pExf(i,j-M+pmlc,k)=(Ca(i,j,k)-dt/e0*e_s2(j-M+pmlc))*pExf(i,j-M+pmlc,k)+dt/e0*e_s2(j-M+pmlc)*Cb(i,j,k)*(Hz(i,j+1,k)-Hz(i,j,k))/dy;
                pEzf(i,j-M+pmlc,k)=(Ca(i,j,k)-dt/e0*e_s2(j-M+pmlc))*pEzf(i,j-M+pmlc,k)-dt/e0*e_s2(j-M+pmlc)*Cb(i,j,k)*(Hx(i,j+1,k)-Hx(i,j,k))/dy;
                Ex(i,j,k)=Ex(i,j,k)-pExf(i,j-M+pmlc,k);
                Ez(i,j,k)=Ez(i,j,k)-pEzf(i,j-M+pmlc,k);
            end
        end
    end
    
    
    for i=1:N-1
        for j=1:M-1
            for k=1:pmlc
                pExb(i,j,k)=(Ca(i,j,k)-dt/e0*e_s1(pmlc-k+1))*pExb(i,j,k)-dt/e0*e_s1(pmlc-k+1)*Cb(i,j,k)*(Hy(i,j,k+1)-Hy(i,j,k))/dz;
                pEyb(i,j,k)=(Ca(i,j,k)-dt/e0*e_s1(pmlc-k+1))*pEyb(i,j,k)+dt/e0*e_s1(pmlc-k+1)*Cb(i,j,k)*(Hx(i,j,k+1)-Hx(i,j,k))/dz;
                Ex(i,j,k)=Ex(i,j,k)-pExb(i,j,k);
                Ey(i,j,k)=Ey(i,j,k)-pEyb(i,j,k);
            end
        end
    end
    
    for i=1:N-1
        for j=1:M-1
            for k=K-pmlc+1:K-1
                pExt(i,j,k-K+pmlc)=(Ca(i,j,k)-dt/e0*e_s2(k-K+pmlc))*pExt(i,j,k-K+pmlc)-dt/e0*e_s2(k-K+pmlc)*Cb(i,j,k)*(Hy(i,j,k+1)-Hy(i,j,k))/dz;
                pEyt(i,j,k-K+pmlc)=(Ca(i,j,k)-dt/e0*e_s2(k-K+pmlc))*pEyt(i,j,k-K+pmlc)+dt/e0*e_s2(k-K+pmlc)*Cb(i,j,k)*(Hx(i,j,k+1)-Hx(i,j,k))/dz;
                Ex(i,j,k)=Ex(i,j,k)-pExt(i,j,k-K+pmlc);
                Ey(i,j,k)=Ey(i,j,k)-pEyt(i,j,k-K+pmlc);
            end
        end
    end    
    
    for i=1:length(Jy_x)
        for j=1:length(gamma_x(i,:))
            Ey(id_surx(i,1),id_surx(i,2),id_surx(i,3))=Ey(id_surx(i,1),id_surx(i,2),id_surx(i,3))-dt/e0/dx*2*(real(Jy_x(i,j))+imag(Jy_x(i,j)));
            Ez(id_surx(i,1),id_surx(i,2),id_surx(i,3))=Ez(id_surx(i,1),id_surx(i,2),id_surx(i,3))-dt/e0/dx*2*(real(Jz_x(i,j))+imag(Jz_x(i,j)));
        end
    end
    
    for i=1:length(Jx_y)
        for j=1:length(gamma_y(i,:))
            Ex(id_sury(i,1),id_sury(i,2),id_sury(i,3))=Ex(id_sury(i,1),id_sury(i,2),id_sury(i,3))-dt/e0/dx*2*(real(Jx_y(i,j))+imag(Jx_y(i,j)));
            Ez(id_sury(i,1),id_sury(i,2),id_sury(i,3))=Ez(id_sury(i,1),id_sury(i,2),id_sury(i,3))-dt/e0/dx*2*(real(Jz_y(i,j))+imag(Jz_y(i,j)));
        end
    end
    
    for i=1:length(Jx_z)
        for j=1:length(gamma_z(i,:))
            Ex(id_surz(i,1),id_surz(i,2),id_surz(i,3))=Ex(id_surz(i,1),id_surz(i,2),id_surz(i,3))-dt/e0/dx*2*(real(Jx_z(i,j))+imag(Jx_z(i,j)));
            Ey(id_surz(i,1),id_surz(i,2),id_surz(i,3))=Ey(id_surz(i,1),id_surz(i,2),id_surz(i,3))-dt/e0/dx*2*(real(Jy_z(i,j))+imag(Jy_z(i,j)));
        end
    end
    
    for i=1:length(EzSource)
        Ez(EzSource(i))=Ez(EzSource(i))-Cb(EzSource(i))*Js(n+1);
    end
    for i=1:length(EySource)
        Ey(EySource(i))=Ey(EySource(i))-Cb(EySource(i))*Js(n+1);
    end
    for i=1:length(ExSource)
        Ex(ExSource(i))=Ex(ExSource(i))-Cb(ExSource(i))*Js(n+1);
    end

    
    for i=1:length(timeSensorLocation)
        timeSensor(n+1,1) = n*dt;
        timeSensor(n+1,2) = Js(n+1);
        timeSensor(n+1,6*(i-1)+3) = Ex(timeSensorLocation(i));
        timeSensor(n+1,6*(i-1)+4) = Ey(timeSensorLocation(i));
        timeSensor(n+1,6*(i-1)+5) = Ez(timeSensorLocation(i));
        timeSensor(n+1,6*(i-1)+6) = Hx(timeSensorLocation(i));
        timeSensor(n+1,6*(i-1)+7) = Hy(timeSensorLocation(i));
        timeSensor(n+1,6*(i-1)+8) = Hz(timeSensorLocation(i));
    end
    
    for i=1:length(sensor.frequencyElectric)
        cos1 = cos(2*pi*sensor.frequencyElectric(i)*n*dt);
        sin1 = sin(2*pi*sensor.frequencyElectric(i)*n*dt);
        fEx(:,:,:,i)=fEx(:,:,:,i)+Ex*(cos1-1i*sin1);
        fEy(:,:,:,i)=fEy(:,:,:,i)+Ey*(cos1-1i*sin1);
        fEz(:,:,:,i)=fEz(:,:,:,i)+Ez*(cos1-1i*sin1);
        freq_normalizationE(i) = freq_normalizationE(i) + Cb(floor(N/2),floor(M/2),floor(K/2))*Js(n+1)*(cos1-1i*sin1);
    end
    
    for i=1:length(sensor.frequencyMagnetic)
        cos1 = cos(2*pi*sensor.frequencyMagnetic(i)*n*dt);
        sin1 = sin(2*pi*sensor.frequencyMagnetic(i)*n*dt);
        fHx(:,:,:,i)=fHx(:,:,:,i)+Hx*(cos1-1i*sin1);
        fHy(:,:,:,i)=fHy(:,:,:,i)+Hy*(cos1-1i*sin1);
        fHz(:,:,:,i)=fHz(:,:,:,i)+Hz*(cos1-1i*sin1);
        freq_normalizationH(i) = freq_normalizationH(i) + Cb(floor(N/2),floor(M/2),floor(K/2))*Js(n+1)*(cos1-1i*sin1);
    end
    
    % pcolor(squeeze(Ez(50,:,:)))
    % axis image
    % shading interp
    % caxis([-1e-5 1e-5]);
    % title(['n=',num2str(n)])
    % getframe;
    
end

for i=1:length(sensor.frequencyElectric)
    fEx(:,:,:,i)=fEx(:,:,:,i)/abs(freq_normalizationE(i));
    fEy(:,:,:,i)=fEy(:,:,:,i)/abs(freq_normalizationE(i));
    fEz(:,:,:,i)=fEz(:,:,:,i)/abs(freq_normalizationE(i));
end

for i=1:length(sensor.frequencyMagnetic)
    fHx(:,:,:,i)=fHx(:,:,:,i)/abs(freq_normalizationH(i));
    fHy(:,:,:,i)=fHy(:,:,:,i)/abs(freq_normalizationH(i));
    fHz(:,:,:,i)=fHz(:,:,:,i)/abs(freq_normalizationH(i));
end

toc


if ~isempty(timeSensorLocation)
    result.timeSensor = timeSensor;
end

if ~isempty(sensor.frequencyElectric)
    result.fEx = fEx;
    result.fEy = fEy;
    result.fEz = fEz;
end

if ~isempty(sensor.frequencyMagnetic)
    result.fHx = fHx;
    result.fHy = fHy;
    result.fHz = fHz;
end





