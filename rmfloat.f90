program atom_type

real :: a1,a2,a3,b1,b2,b3, c1,c2,c3
real :: a,b,c
real :: alpha,beta,gamma
real :: volume
character (len=5) :: name(50000), expanded_name(500000)
real :: coordinates(50000,3), charge(50000), C_coordinates(50000,3), list(50000)=1.0
real :: expanded_coordinates(500000,3)
real :: X, Y, Z, Dis
integer :: iii, jjj, zzz, Toc, Toc_expanded
integer :: num(8), typeid(50000,8)
integer :: ii, jj, kk, i, j, m
integer :: index
integer :: n, l
character (len=5) :: rank(8), type(50000,8)
real :: distance(8)
real :: R_cut
real :: pi = 3.14159265359
real :: rad 
real :: radii(10000)
character (len=20) :: dummy(10000)



! here the program read a cif format 
OPEN (UNIT = 1, FILE = "./FILES.dat",ACTION='READ',STATUS='OLD');
OPEN (UNIT = 2, FILE = "./FILES.dat.txt",ACTION='WRITE',STATUS='NEW');

!read a/b/c and alpha/beta/gamma values
read(1,*)a
read(1,*)b
read(1,*)c
read(1,*)alpha
read(1,*)beta
read(1,*)gamma
read(1,*)volume
read(1,*)Toc

rad= pi/180

! This is lattice vector
! Lattice vector can be computed from a, b, c and alpha, beta, and gamma --> need your implemnetation

a1=    a    ; b1=    b*COS(gamma*rad) ;  c1=   c*COS(beta*rad)
a2=     0.000     ; b2=   b*SIN(gamma*rad) ;  c2=    c*((COS(alpha*rad)-COS(beta*rad)*COS(gamma*rad))/SIN(gamma*rad))
a3=     0.000     ; b3=    0.000 ;  c3=   volume/(a*b*SIN(gamma*rad))

! expand the structure to a 3-3-3 super cells --> no need to implement periodic anymore
Toc_expanded=Toc*27

!real file first
do iii=1,Toc

   read(1,*) name(iii), coordinates(iii,1), coordinates(iii,2), coordinates(iii,3)
enddo

!expand the structure, according to the lattice vectors
index=1

do ii=-1,1
do jj=-1,1
do kk=-1,1

    
   do iii=1,Toc

   expanded_name(index)=name(iii)
   expanded_coordinates(index,1)= (coordinates(iii,1)+real(ii))*a1+(coordinates(iii,2)+real(jj))*b1+(coordinates(iii,3)+real(kk))*c1
   expanded_coordinates(index,2)= (coordinates(iii,1)+real(ii))*a2+(coordinates(iii,2)+real(jj))*b2+(coordinates(iii,3)+real(kk))*c2
   expanded_coordinates(index,3)= (coordinates(iii,1)+real(ii))*a3+(coordinates(iii,2)+real(jj))*b3+(coordinates(iii,3)+real(kk))*c3
   index=index+1
   enddo
enddo
enddo
enddo

index=1

do iii=1,Toc 
  
C_coordinates(index,1)=(coordinates(iii,1))*a1+(coordinates(iii,2))*b1+(coordinates(iii,3))*c1
C_coordinates(index,2)=(coordinates(iii,1))*a2+(coordinates(iii,2))*b2+(coordinates(iii,3))*c2
C_coordinates(index,3)=(coordinates(iii,1))*a3+(coordinates(iii,2))*b3+(coordinates(iii,3))*c3
index=index+1   
enddo

!start to calculate the distance
do jjj=1,Toc
!initialized the maxima point
distance(1)=100.0; distance(2)=100.0; distance(3)=100.0; distance(4)=100.0; distance(5)=100.0;
distance(6)=100.0; distance(7)=100.0;
rank(1)='null'; rank(2)='null'; rank(3)='null'; rank(4)='null'; rank(5)='null'; rank(6)='null'; rank(7)='null'; rank(8)='null'
num(1)=0; num(2)=0; num(3)=0; num(4)=0; num(5)=0; num(6)=0; num(7)=0; num(8)=0
   do zzz=1,Toc_expanded
  
         X=abs(C_coordinates(jjj,1)-expanded_coordinates(zzz,1))
         Y=abs(C_coordinates(jjj,2)-expanded_coordinates(zzz,2)) 
         Z=abs(C_coordinates(jjj,3)-expanded_coordinates(zzz,3))
   

         Dis=sqrt(X**2+Y**2+Z**2);  
              
			   !find radii
			   !of element
				if (expanded_name(jjj) == "Zn") then
					radii(1) = 1.193;
				elseif (expanded_name(jjj) == "C") then
					radii(1) = 0.757;
				elseif (expanded_name(jjj) == "Cl") then
					radii(1) = 1.044;
				elseif (expanded_name(jjj) == "Br") then
					radii(1) = 1.192;
				elseif (expanded_name(jjj) == "F") then
					radii(1) = 0.668;
				elseif (expanded_name(jjj) == "H") then
					radii(1) = 0.354;
				elseif (expanded_name(jjj) == "I") then
					radii(1) = 1.382;
				elseif (expanded_name(jjj) == "N") then
					radii(1) = 0.7;
				elseif (expanded_name(jjj) == "O") then
					radii(1) = 0.634;
				elseif (expanded_name(jjj) == "P") then
					radii(1) = 1.101;
				elseif (expanded_name(jjj) == "S") then
					radii(1) = 1.064;
				elseif (expanded_name(jjj) == "Ag") then
					radii(1) = 1.386;
				elseif (expanded_name(jjj) == "Al") then
					radii(1) = 1.244;
				elseif (expanded_name(jjj) == "As") then
					radii(1) = 1.211;
				elseif (expanded_name(jjj) == "Au") then
					radii(1) = 1.262;
				elseif (expanded_name(jjj) == "B") then
					radii(1) = 0.838;
				elseif (expanded_name(jjj) == "Ba") then
					radii(1) = 2.277;
				elseif (expanded_name(jjj) == "Be") then
					radii(1) = 1.074;
				elseif (expanded_name(jjj) == "Bi") then
					radii(1) = 1.512;
				elseif (expanded_name(jjj) == "Ca") then
					radii(1) = 1.761;
				elseif (expanded_name(jjj) == "Cd") then
					radii(1) = 1.403;
				elseif (expanded_name(jjj) == "Ce") then
					radii(1) = 1.841;
				elseif (expanded_name(jjj) == "Cl") then
					radii(1) = 1.044;
				elseif (expanded_name(jjj) == "Co") then
					radii(1) = 1.241;
				elseif (expanded_name(jjj) == "Cr") then
					radii(1) = 1.345;
				elseif (expanded_name(jjj) == "Cu") then
					radii(1) = 1.302;
				elseif (expanded_name(jjj) == "Dy") then
					radii(1) = 1.71;	
				elseif (expanded_name(jjj) == "Er") then
					radii(1) = 1.673;
				elseif (expanded_name(jjj) == "Fe") then
					radii(1) = 1.335;
				elseif (expanded_name(jjj) == "Ga") then
					radii(1) = 1.26;
				elseif (expanded_name(jjj) == "Gd") then
					radii(1) = 1.735;
				elseif (expanded_name(jjj) == "Ge") then
					radii(1) = 1.197;
				elseif (expanded_name(jjj) == "Hf") then
					radii(1) = 1.611;
				elseif (expanded_name(jjj) == "Hg") then
					radii(1) = 1.34;
				elseif (expanded_name(jjj) == "Ho") then
					radii(1) = 1.696;
				elseif (expanded_name(jjj) == "In") then
					radii(1) = 1.459;
				elseif (expanded_name(jjj) == "Ir") then
					radii(1) = 1.371;
				elseif (expanded_name(jjj) == "K") then
					radii(1) = 1.953;
				elseif (expanded_name(jjj) == "La") then
					radii(1) = 1.943;
				elseif (expanded_name(jjj) == "Li") then
					radii(1) = 1.336;
				elseif (expanded_name(jjj) == "Lu") then
					radii(1) = 1.671;
				elseif (expanded_name(jjj) == "Mg") then
					radii(1) = 1.421;
				elseif (expanded_name(jjj) == "Mn") then
					radii(1) = 1.382;
				elseif (expanded_name(jjj) == "Mo") then
					radii(1) = 1.47;
				elseif (expanded_name(jjj) == "Na") then
					radii(1) = 1.539;
				elseif (expanded_name(jjj) == "Nb") then
					radii(1) = 1.473;
				elseif (expanded_name(jjj) == "Nd") then
					radii(1) = 1.816;
				elseif (expanded_name(jjj) == "Ni") then
					radii(1) = 1.164;
				elseif (expanded_name(jjj) == "Np") then
					radii(1) = 1.666;
				elseif (expanded_name(jjj) == "Pb") then
					radii(1) = 1.459;
				elseif (expanded_name(jjj) == "Pd") then
					radii(1) = 1.338;
				elseif (expanded_name(jjj) == "Pr") then
					radii(1) = 1.823;
				elseif (expanded_name(jjj) == "Pt") then
					radii(1) = 1.364;
				elseif (expanded_name(jjj) == "Pu") then
					radii(1) = 1.657;
				elseif (expanded_name(jjj) == "Rb") then
					radii(1) = 2.26;
				elseif (expanded_name(jjj) == "Re") then
					radii(1) = 1.343;
				elseif (expanded_name(jjj) == "Rh") then
					radii(1) = 1.332;
				elseif (expanded_name(jjj) == "Ru") then
					radii(1) = 1.478;
				elseif (expanded_name(jjj) == "Sb") then
					radii(1) = 1.407;
				elseif (expanded_name(jjj) == "Sc") then
					radii(1) = 1.513;
				elseif (expanded_name(jjj) == "Se") then
					radii(1) = 1.19;
				elseif (expanded_name(jjj) == "Si") then
					radii(1) = 1.117;
				elseif (expanded_name(jjj) == "Sm") then
					radii(1) = 1.78;
				elseif (expanded_name(jjj) == "Sn") then
					radii(1) = 1.398;
				elseif (expanded_name(jjj) == "Sr") then
					radii(1) = 2.052;
				elseif (expanded_name(jjj) == "Te") then
					radii(1) = 1.386;
				elseif (expanded_name(jjj) == "Th") then
					radii(1) = 1.721;
				elseif (expanded_name(jjj) == "Ti") then
					radii(1) = 1.412;
				elseif (expanded_name(jjj) == "Tm") then
					radii(1) = 1.66;
				elseif (expanded_name(jjj) == "U") then
					radii(1) = 1.684;
				elseif (expanded_name(jjj) == "V") then
					radii(1) = 1.402;
				elseif (expanded_name(jjj) == "W") then
					radii(1) = 1.392;
				elseif (expanded_name(jjj) == "Y") then
					radii(1) = 1.698;
				elseif (expanded_name(jjj) == "Yb") then
					radii(1) = 1.637;
				elseif (expanded_name(jjj) == "He") then
					radii(1) = 0.849;
				elseif (expanded_name(jjj) == "Ne") then
					radii(1) = 0.920;
				elseif (expanded_name(jjj) == "Zr") then
					radii(1) = 1.564;
				elseif (expanded_name(jjj) == "Tc") then
					radii(1) = 1.322;
				elseif (expanded_name(jjj) == "Xe") then
					radii(1) = 1.267;
				elseif (expanded_name(jjj) == "Cs") then
					radii(1) = 2.57;
				elseif (expanded_name(jjj) == "Pm") then
					radii(1) = 1.801;
				elseif (expanded_name(jjj) == "Eu") then
					radii(1) = 1.771;
				elseif (expanded_name(jjj) == "Tb") then
					radii(1) = 1.732;
				elseif (expanded_name(jjj) == "Ta") then
					radii(1) = 1.511;
				elseif (expanded_name(jjj) == "Os") then
					radii(1) = 1.372;
				elseif (expanded_name(jjj) == "Po") then
					radii(1) = 1.50;
				elseif (expanded_name(jjj) == "At") then
					radii(1) = 1.545;
				elseif (expanded_name(jjj) == "Rn") then
					radii(1) = 1.420;
				elseif (expanded_name(jjj) == "Fr") then
					radii(1) = 2.880;
				elseif (expanded_name(jjj) == "Ra") then
					radii(1) = 2.512;
				elseif (expanded_name(jjj) == "Ac") then
					radii(1) = 1.983;
				elseif (expanded_name(jjj) == "Pa") then
					radii(1) = 01.711;
				elseif (expanded_name(jjj) == "Am") then
					radii(1) = 1.660;
				elseif (expanded_name(jjj) == "Cm") then
					radii(1) = 1.801;
				elseif (expanded_name(jjj) == "Bk") then
					radii(1) = 1.761;
				elseif (expanded_name(jjj) == "Cf") then
					radii(1) = 1.750;
				elseif (expanded_name(jjj) == "Es") then
					radii(1) = 1.724;
				elseif (expanded_name(jjj) == "Fm") then
					radii(1) = 1.712;
				elseif (expanded_name(jjj) == "Md") then
					radii(1) = 1.689;
				elseif (expanded_name(jjj) == "No") then
					radii(1) = 1.679;
				elseif (expanded_name(jjj) == "Lw") then
					radii(1) = 1.698;
				else
				    radii(1) = 1.500; 
				endif
				!of connected atom
				if (expanded_name(zzz) == "Zn") then
					radii(2) = 1.193;
				elseif (expanded_name(zzz) == "C") then
					radii(2) = 0.757;
				elseif (expanded_name(zzz) == "Cl") then
					radii(2) = 1.044;
				elseif (expanded_name(zzz) == "Br") then
					radii(2) = 1.192;
				elseif (expanded_name(zzz) == "F") then
					radii(2) = 0.668;
				elseif (expanded_name(zzz) == "H") then
					radii(2) = 0.354;
				elseif (expanded_name(zzz) == "I") then
					radii(2) = 1.382;
				elseif (expanded_name(zzz) == "N") then
					radii(2) = 0.7;
				elseif (expanded_name(zzz) == "O") then
					radii(2) = 0.634;
				elseif (expanded_name(zzz) == "P") then
					radii(2) = 1.101;
				elseif (expanded_name(zzz) == "S") then
					radii(2) = 1.064;
				elseif (expanded_name(zzz) == "Ag") then
					radii(2) = 1.386;
				elseif (expanded_name(zzz) == "Al") then
					radii(2) = 1.244;
				elseif (expanded_name(zzz) == "As") then
					radii(2) = 1.211;
				elseif (expanded_name(zzz) == "Au") then
					radii(2) = 1.262;
				elseif (expanded_name(zzz) == "B") then
					radii(2) = 0.838;
				elseif (expanded_name(zzz) == "Ba") then
					radii(2) = 2.277;
				elseif (expanded_name(zzz) == "Be") then
					radii(2) = 1.074;
				elseif (expanded_name(zzz) == "Bi") then
					radii(2) = 1.512;
				elseif (expanded_name(zzz) == "Ca") then
					radii(2) = 1.761;
				elseif (expanded_name(zzz) == "Cd") then
					radii(2) = 1.403;
				elseif (expanded_name(zzz) == "Ce") then
					radii(2) = 1.841;
				elseif (expanded_name(zzz) == "Cl") then
					radii(2) = 1.044;
				elseif (expanded_name(zzz) == "Co") then
					radii(2) = 1.241;
				elseif (expanded_name(zzz) == "Cr") then
					radii(2) = 1.345;
				elseif (expanded_name(zzz) == "Cu") then
					radii(2) = 1.302;
				elseif (expanded_name(zzz) == "Dy") then
					radii(2) = 1.71;	
				elseif (expanded_name(zzz) == "Er") then
					radii(2) = 1.673;
				elseif (expanded_name(zzz) == "Fe") then
					radii(2) = 1.335;
				elseif (expanded_name(zzz) == "Ga") then
					radii(2) = 1.26;
				elseif (expanded_name(zzz) == "Gd") then
					radii(2) = 1.735;
				elseif (expanded_name(zzz) == "Ge") then
					radii(2) = 1.197;
				elseif (expanded_name(zzz) == "Hf") then
					radii(2) = 1.611;
				elseif (expanded_name(zzz) == "Hg") then
					radii(2) = 1.34;
				elseif (expanded_name(zzz) == "Ho") then
					radii(2) = 1.696;
				elseif (expanded_name(zzz) == "In") then
					radii(2) = 1.459;
				elseif (expanded_name(zzz) == "Ir") then
					radii(2) = 1.371;
				elseif (expanded_name(zzz) == "K") then
					radii(2) = 1.953;
				elseif (expanded_name(zzz) == "La") then
					radii(2) = 1.943;
				elseif (expanded_name(zzz) == "Li") then
					radii(2) = 1.336;
				elseif (expanded_name(zzz) == "Lu") then
					radii(2) = 1.671;
				elseif (expanded_name(zzz) == "Mg") then
					radii(2) = 1.421;
				elseif (expanded_name(zzz) == "Mn") then
					radii(2) = 1.382;
				elseif (expanded_name(zzz) == "Mo") then
					radii(2) = 1.47;
				elseif (expanded_name(zzz) == "Na") then
					radii(2) = 1.539;
				elseif (expanded_name(zzz) == "Nb") then
					radii(2) = 1.473;
				elseif (expanded_name(zzz) == "Nd") then
					radii(2) = 1.816;
				elseif (expanded_name(zzz) == "Ni") then
					radii(2) = 1.164;
				elseif (expanded_name(zzz) == "Np") then
					radii(2) = 1.666;
				elseif (expanded_name(zzz) == "Pb") then
					radii(2) = 1.459;
				elseif (expanded_name(zzz) == "Pd") then
					radii(2) = 1.338;
				elseif (expanded_name(zzz) == "Pr") then
					radii(2) = 1.823;
				elseif (expanded_name(zzz) == "Pt") then
					radii(2) = 1.364;
				elseif (expanded_name(zzz) == "Pu") then
					radii(2) = 1.657;
				elseif (expanded_name(zzz) == "Rb") then
					radii(2) = 2.26;
				elseif (expanded_name(zzz) == "Re") then
					radii(2) = 1.343;
				elseif (expanded_name(zzz) == "Rh") then
					radii(2) = 1.332;
				elseif (expanded_name(zzz) == "Ru") then
					radii(2) = 1.478;
				elseif (expanded_name(zzz) == "Sb") then
					radii(2) = 1.407;
				elseif (expanded_name(zzz) == "Sc") then
					radii(2) = 1.513;
				elseif (expanded_name(zzz) == "Se") then
					radii(2) = 1.19;
				elseif (expanded_name(zzz) == "Si") then
					radii(2) = 1.117;
				elseif (expanded_name(zzz) == "Sm") then
					radii(2) = 1.78;
				elseif (expanded_name(zzz) == "Sn") then
					radii(2) = 1.398;
				elseif (expanded_name(zzz) == "Sr") then
					radii(2) = 2.052;
				elseif (expanded_name(zzz) == "Te") then
					radii(2) = 1.386;
				elseif (expanded_name(zzz) == "Th") then
					radii(2) = 1.721;
				elseif (expanded_name(zzz) == "Ti") then
					radii(2) = 1.412;
				elseif (expanded_name(zzz) == "Tm") then
					radii(2) = 1.66;
				elseif (expanded_name(zzz) == "U") then
					radii(2) = 1.684;
				elseif (expanded_name(zzz) == "V") then
					radii(2) = 1.402;
				elseif (expanded_name(zzz) == "W") then
					radii(2) = 1.392;
				elseif (expanded_name(zzz) == "Y") then
					radii(2) = 1.698;
				elseif (expanded_name(zzz) == "Yb") then
					radii(2) = 1.637;
				elseif (expanded_name(zzz) == "He") then
					radii(2) = 0.849;
				elseif (expanded_name(zzz) == "Ne") then
					radii(2) = 0.920;
				elseif (expanded_name(zzz) == "Zr") then
					radii(2) = 1.564;
				elseif (expanded_name(zzz) == "Tc") then
					radii(2) = 1.322;
				elseif (expanded_name(zzz) == "Xe") then
					radii(2) = 1.267;
				elseif (expanded_name(zzz) == "Cs") then
					radii(2) = 2.57;
				elseif (expanded_name(zzz) == "Pm") then
					radii(2) = 1.801;
				elseif (expanded_name(zzz) == "Eu") then
					radii(2) = 1.771;
				elseif (expanded_name(zzz) == "Tb") then
					radii(2) = 1.732;
				elseif (expanded_name(zzz) == "Ta") then
					radii(2) = 1.511;
				elseif (expanded_name(zzz) == "Os") then
					radii(2) = 1.372;
				elseif (expanded_name(zzz) == "Po") then
					radii(2) = 1.50;
				elseif (expanded_name(zzz) == "At") then
					radii(2) = 1.545;
				elseif (expanded_name(zzz) == "Rn") then
					radii(2) = 1.420;
				elseif (expanded_name(zzz) == "Fr") then
					radii(2) = 2.880;
				elseif (expanded_name(zzz) == "Ra") then
					radii(2) = 2.512;
				elseif (expanded_name(zzz) == "Ac") then
					radii(2) = 1.983;
				elseif (expanded_name(zzz) == "Pa") then
					radii(2) = 1.711;
				elseif (expanded_name(zzz) == "Am") then
					radii(2) = 1.660;
				elseif (expanded_name(zzz) == "Cm") then
					radii(2) = 1.801;
				elseif (expanded_name(zzz) == "Bk") then
					radii(2) = 1.761;
				elseif (expanded_name(zzz) == "Cf") then
					radii(2) = 1.750;
				elseif (expanded_name(zzz) == "Es") then
					radii(2) = 1.724;
				elseif (expanded_name(zzz) == "Fm") then
					radii(2) = 1.712;
				elseif (expanded_name(zzz) == "Md") then
					radii(2) = 1.689;
				elseif (expanded_name(zzz) == "No") then
					radii(2) = 1.679;
				elseif (expanded_name(zzz) == "Lw") then
					radii(2) = 1.698;
				else
				    radii(2) = 1.500; 
				endif
				R_cut = (radii(1) + radii(2))*1.25	
				if (Dis <distance(1) .and. Dis <R_cut) then
                 distance(8)=distance(7); rank(8)=rank(7); num(8)=num(7);
				 distance(7)=distance(6); rank(7)=rank(6); num(7)=num(6);
				 distance(6)=distance(5); rank(6)=rank(5); num(6)=num(5);
				 distance(5)=distance(4); rank(5)=rank(4); num(5)=num(4);
				 distance(4)=distance(3); rank(4)=rank(3); num(4)=num(3);
				 distance(3)=distance(2); rank(3)=rank(2); num(3)=num(2);
				 distance(2)=distance(1); rank(2)=rank(1); num(2)=num(1);
				 distance(1)=Dis; rank(1)=expanded_name(zzz); flo = zzz/Toc; num(1)= zzz - floor(flo)*Toc;
					if (floor(flo)*Toc == zzz) then
						num(1) = Toc
					else
					endif
				 
				else if (Dis <distance(2) .and. Dis <R_cut) then
                 distance(8)=distance(7); rank(8)=rank(7); num(8)=num(7);    
				 distance(7)=distance(6); rank(7)=rank(6); num(7)=num(6);
				 distance(6)=distance(5); rank(6)=rank(5); num(6)=num(5);
				 distance(5)=distance(4); rank(5)=rank(4); num(5)=num(4);
				 distance(4)=distance(3); rank(4)=rank(3); num(4)=num(3);
				 distance(3)=distance(2); rank(3)=rank(2); num(3)=num(2);
				 distance(2)=Dis; rank(2)=expanded_name(zzz); flo = zzz/Toc; num(2)= zzz - floor(flo)*Toc;
				 if (floor(flo)*Toc == zzz) then
						num(2) = Toc
				else
				endif
				else if (Dis <distance(3) .and. Dis <R_cut) then
                 distance(8)=distance(7); rank(8)=rank(7); num(8)=num(7);
				 distance(7)=distance(6); rank(7)=rank(6); num(7)=num(6);
				 distance(6)=distance(5); rank(6)=rank(5); num(6)=num(5);
				 distance(5)=distance(4); rank(5)=rank(4); num(5)=num(4);
				 distance(4)=distance(3); rank(4)=rank(3); num(4)=num(3);
				 distance(3)=Dis; rank(3)=expanded_name(zzz); flo = zzz/Toc; num(3)= zzz - floor(flo)*Toc;
				 if (floor(flo)*Toc == zzz) then
						num(3) = Toc
				else
				endif
				else if (Dis <distance(4) .and. Dis <R_cut) then
                 distance(8)=distance(7); rank(8)=rank(7); num(8)=num(7);
				 distance(7)=distance(6); rank(7)=rank(6); num(7)=num(6);
				 distance(6)=distance(5); rank(6)=rank(5); num(6)=num(5);
				 distance(5)=distance(4); rank(5)=rank(4); num(5)=num(4);
				distance(4)=Dis; rank(4)=expanded_name(zzz); flo = zzz/Toc; num(4)= zzz - floor(flo)*Toc;
				if (floor(flo)*Toc == zzz) then
						num(4) = Toc
				else
				endif
				else if (Dis < distance(5) .and. Dis <R_cut) then
                 distance(8)=distance(7); rank(8)=rank(7); num(8)=num(7);
				 distance(7)=distance(6); rank(7)=rank(6); num(7)=num(6);
				 distance(6)=distance(5); rank(6)=rank(5); num(6)=num(5);
				 distance(5)=Dis; rank(5)=expanded_name(zzz); flo = zzz/Toc; num(5)= zzz - floor(flo)*Toc;
				 if (floor(flo)*Toc == zzz) then
						num(5) = Toc
				else
				endif
				else if (Dis < distance(6) .and. Dis <R_cut) then
                 distance(8)=distance(7); rank(8)=rank(7); num(8)=num(7);
				 distance(7)=distance(6); rank(7)=rank(6); num(7)=num(6);
				 distance(6)=Dis; rank(6)=expanded_name(zzz); flo = zzz/Toc; num(6)= zzz - floor(flo)*Toc;
				 if (floor(flo)*Toc == zzz) then
						num(6) = Toc
				else
				endif
				else if (Dis < distance(7) .and. Dis <R_cut) then
                 distance(8)=distance(7); rank(8)=rank(7); num(8)=num(7);
				 distance(7)=Dis; rank(7)=expanded_name(zzz); flo = zzz/Toc; num(7)= zzz - floor(flo)*Toc;
				 if (floor(flo)*Toc == zzz) then
						num(7) = Toc
				else
				endif
				else if (Dis < distance(8) .and. Dis <R_cut) then
				 distance(8)=Dis; rank(8)=expanded_name(zzz); flo = zzz/Toc; num(8)= zzz - floor(flo)*Toc;
				if (floor(flo)*Toc == zzz) then
						num(8) = Toc
				else
				endif
                                
				endif
				
    enddo
                    
         write (2,*) rank, num
	     !write (2,*) distance
enddo  
close (2)       
OPEN (UNIT = 3, FILE = "./FILES.dat.txt",ACTION='READ',STATUS='OLD');
OPEN (UNIT = 4, FILE = "./FILES.checkfloat.txt",ACTION='WRITE',STATUS='NEW');
OPEN (UNIT = 5, FILE = "./FILES.freefloat.txt",ACTION='WRITE',STATUS='NEW');
do iii=1,Toc

   read(3,*) (type(iii,j), j=1,8), (typeid(iii,j), j = 1,8)

enddo

!check float atoms
do iii=1,Toc
  m = count(typeid(iii,:)/=0)
  l = 0
  do i=2,m
    n = count(typeid(typeid(iii,i),:)/=0)
    l = l+n-1
  enddo
  if (m >= l+1) then
	write(4,*) "float"
	do i=1,m
	index = typeid(iii,i)
	list(index) = 0.0
	enddo
  else
  endif
 enddo
do iii=1,Toc
   if (list(iii) > 0.0) then
   write(5,*) name(iii),name(iii), coordinates(iii,1), coordinates(iii,2), coordinates(iii,3), 0 
   else
   endif
enddo
endprogram
