#!/bin/sh
#!/bin/python
#echo What is the filename?
#read FILES
FILES=`ls *.cif`
#dos2unix *.cif 
rm float.list
 for f in $FILES
 do
 rm ${f}*.dat ${f}*.txt
 echo $f
 awk '$1 == "_cell_length_a" {print $2}' $f > temp
 a=`awk '$1 == "_cell_length_a" {print $2}' $f`
 echo $a
 awk '$1 == "_cell_length_b" {print $2}' $f >> temp
 b=`awk '$1 == "_cell_length_b" {print $2}' $f`
 echo $b
 awk '$1 == "_cell_length_c" {print $2}' $f >> temp
 c=`awk '$1 == "_cell_length_c" {print $2}' $f`
 echo $c
 awk '$1 == "_cell_angle_alpha" {print $2}' $f >> temp
 alpha=`awk '$1 == "_cell_angle_alpha" {print $2}' $f `
 echo $alpha
 awk '$1 == "_cell_angle_beta" {print $2}' $f >> temp
 beta=`awk '$1 == "_cell_angle_beta" {print $2}' $f `
 echo $beta
 awk '$1 == "_cell_angle_gamma" {print $2}' $f >> temp
 gamma=`awk '$1 == "_cell_angle_gamma" {print $2}' $f`
 echo $gamma
 rlabel=`awk '$1 == "_atom_site_label" {print NR}' $f`
 echo $rlabel
 nlabel=1
 nsymbol=`awk '$1 == "_atom_site_type_symbol" {print NR-'$rlabel'+1}' $f`
 echo $nsymbol
 nfrax=`awk '$1 == "_atom_site_fract_x" {print NR-'$rlabel'+1}' $f`
 nfray=`awk '$1 == "_atom_site_fract_y" {print NR-'$rlabel'+1}' $f`
 nfraz=`awk '$1 == "_atom_site_fract_z" {print NR-'$rlabel'+1}' $f`
 echo $nfrax $nfray $nfraz
 awk 'BEGIN{print '$a'*'$b'*'$c'*(1-(cos(atan2(0, -1)*'$alpha'/180))**2-(cos(atan2(0, -1)*'$beta'/180))**2-(cos(atan2(0, -1)*'$gamma'/180))**2+2*cos(atan2(0, -1)*'$alpha'/180)*cos(atan2(0, -1)*'$beta'/180)*cos(atan2(0, -1)*'$gamma'/180))**0.5}' >> temp
 awk 'NR >= 25 {print $'$nsymbol'" "$'$nfrax'" "$'$nfray'" "$'$nfraz'" " "0"}' $f > $f.dat
 awk 'END {print NR}' $f.dat >> temp
 cat $f.dat >> temp
 cp temp $f.dat

 cp secondlayer.f90 temp.f90

 sed -i 's/FILES/'${f}'/g' temp.f90
 gfortran temp.f90 -o temp
 ./temp
 rm temp.f90 temp


 file1=${f}.checkfloat.txt
 awk -vf="$f" 'NR==1 && $1 {print f}' $file1 >> float.list
 
 rm ${f}*.dat ${f}*.txt
 done


 