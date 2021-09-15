mkdir ./Float_Structures/
list=`awk '{print}' float.list`
for f in $list
do
 rm head_${f}.txt
 rm ${f}.dat ${f}.dat.* 
 rm Final_Charge_${f}*
 rm *_Charge_${f}*
 rm tail_*
 rm FINAL_${f}*.cif
 rm ${f}.checkfloat.txt
 rm ${f}.freefloat.txt
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
 nsymbol=$((nsymbol+0))
 if [ $nsymbol -eq 0 ]; then
 nsymbol=1
 fi
 echo $nsymbol
 nfrax=`awk '$1 == "_atom_site_fract_x" {print NR-'$rlabel'+1}' $f`
 nfray=`awk '$1 == "_atom_site_fract_y" {print NR-'$rlabel'+1}' $f`
 nfraz=`awk '$1 == "_atom_site_fract_z" {print NR-'$rlabel'+1}' $f`
 echo $nfrax $nfray $nfraz
 crow=` awk -F"_" '$2=="atom" {a=NR} END{print a+1}' $f`
 awk 'BEGIN{print '$a'*'$b'*'$c'*(1-(cos(atan2(0, -1)*'$alpha'/180))**2-(cos(atan2(0, -1)*'$beta'/180))**2-(cos(atan2(0, -1)*'$gamma'/180))**2+2*cos(atan2(0, -1)*'$alpha'/180)*cos(atan2(0, -1)*'$beta'/180)*cos(atan2(0, -1)*'$gamma'/180))**0.5}' >> temp
 awk 'NR >= '$crow' {print $'$nsymbol'" "$'$nfrax'" "$'$nfray'" "$'$nfraz'" " "0" " "}' $f > $f.dat
 awk 'END {print NR}' $f.dat >> temp
 cat $f.dat >> temp
 mv temp $f.dat

 cp rmfloat.f90 temp.f90

 sed -i 's/FILES/'${f}'/g' temp.f90
 gfortran temp.f90 -o temp
 ./temp
 rm temp.f90 temp
 num=`awk -F"_" '$2=="atom"{print NR-1;exit}' $f`
file2=${f}.freefloat.txt
#num=$(awk '$1 == "_atom_site_fract_z" {print NR}' $f)
head -${num} ${f} > head_${f}.txt

echo " _atom_site_label" >> head_${f}.txt
echo " _atom_site_type_symbol" >> head_${f}.txt
echo " _atom_site_fract_x" >> head_${f}.txt
echo " _atom_site_fract_y" >> head_${f}.txt
echo " _atom_site_fract_z" >> head_${f}.txt
echo " _atom_site_charge" >> head_${f}.txt

#awk '{print $1,$0}' ${file2} > new_tail.txt
#paste -d ' ' tail_1.txt $file2 > new_tail.txt
cat head_${f}.txt ${file2} > Freefloat_${f}
mv ${f} Float_Structures/
done