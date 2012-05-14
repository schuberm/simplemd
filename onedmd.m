function onedmd();

 format long;
 global sys;
 
 sys.nseeds=1;
 sys.natoms=4; %Even number
 sys.eq=1:1:sys.natoms;
 sys.eq=sys.eq';
 sys.dt=0.01;
 sys.nsteps=10000;

  %Create vector of m
  %sys.mass=ones(sys.natoms,1);
  sys.mass=repmat([1 1],1,sys.natoms/2)';
  MM=diag(sys.mass);

 %harmonic component of the potential
 %Create vector of k: ( 1 2 1 2... or 1 2 4 8 ... or 1 1 2 3 5 ... etc) 
  %k=1:1:sys.natoms;
  %k=ones(1,sys.natoms);
  k=repmat([1 1],1,sys.natoms/2)';
  sys.K=(diag(k,0)+diag(-k(1:sys.natoms-1),-1)+diag(-k(sys.natoms:sys.natoms),sys.natoms-1));
  DM=diag(k,0)+diag(circshift(k,1),0)+diag(-k(1:sys.natoms-1),-1)+diag(-k(1:sys.natoms-1),1)+diag(-k(sys.natoms:sys.natoms),sys.natoms-1)+diag(-k(sys.natoms:sys.natoms),1-sys.natoms);
 
  %anharmonic component of the potential	 
  sys.a=0.1;
  sys.A=(diag(sys.a*ones(sys.natoms,1))+diag(-sys.a*ones(sys.natoms-1,1),-1)+diag(-sys.a*ones(1,1),sys.natoms-1));
  
  %Eigenvectors for normal mode transformation	
  Kt=MM^-.5*DM*MM^-.5;
  [E,omega]=eig(Kt);

  %number of independent simulations to run
 for n=1:1:sys.nseeds
  %initial positions
  sys.xinit=sys.eq;
  sys.string='';
  for l=1:1:sys.natoms
	sys.string=strcat(sys.string,'% 20.8f');
        sys.xinit(l)=sys.xinit(l)+0.1*(2*rand-1);
  end	  
  sys.x=sys.xinit;

  %initial velocities for no net momentum
  sys.velinit=repmat([1 -1],1,sys.natoms/2);
  sys.vel=sys.velinit';

  sys.F=zeros(sys.natoms,1);
  sys.delta=zeros(sys.natoms,1);
  sys.epot=zeros(sys.natoms,1);
  sys.ekin=zeros(sys.natoms,1);
  sys.etot=zeros(sys.natoms,1);
  sys.nmode=zeros(sys.natoms,1);
  sys.nmodemini=zeros(sys.natoms/2,1);
   

  %Initialize
  force();

  sys.velfile = strcat('vel',int2str(n),'.dat');
  sys.energyfile=strcat('energy',int2str(n),'.dat');
  if sys.a==0
   sys.nmfile=strcat('nmode',int2str(n),'.dat');
   sys.nmminifile=strcat('nmodemini',int2str(n),'.dat');
  else
   sys.nmfile=strcat('nmode_a',int2str(n),'.dat');
   sys.nmminifile=strcat('nmodemini_a',int2str(n),'.dat');
  end

  velf = fopen(sys.velfile, 'w');
  energyf = fopen(sys.energyfile, 'w');
  nmf=fopen(sys.nmfile, 'w');
  nmfmini=fopen(sys.nmminifile, 'w'); 
 
  %main md loop
  for i= 1 : 1 : sys.nsteps
	%velocity verlet
	sys.vel=sys.vel+0.5*sys.dt.*sys.F;
	sys.x=sys.x+sys.dt.*sys.vel;
	force();
	sys.vel=sys.vel+0.5*sys.dt.*sys.F;

	fprintf(velf,strcat('',sys.string,'\n'),sys.vel(:)');
        sys.nmode=inv(E)*sys.vel;
        fprintf(nmf,strcat('',sys.string,'\n'),sys.nmode(:)');
	sys.ekin=0.5*sys.mass*sys.vel.^2;
        sys.etot=sys.epot+sys.ekin;
        fprintf(energyf,strcat('',sys.string,'\n'),sys.etot(:)');
  end
  fclose(velf);
  fclose(energyf);
  fclose(nmf);
  fclose(nmfmini);
 end

%vectorized force calculation
 function force()
	global sys;
	sys.diff=sys.eq-sys.x;
	sys.delta=sys.diff-circshift(sys.diff,-1);
	sys.F=sys.K*sys.delta+sys.A*sys.delta.^3;
        sys.epot=0.5*sys.k*sys.delta.^2+0.25*sys.a*sys.delta.^4;



