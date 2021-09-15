####################################################################################################################
#This is the beta version of 2nd-CBAC method for quick charge assignment of MOFs                                   #
#Copyright Li-Chiang Lin, Computational Material Discovery Group, College of Enginnering, The Ohio State University#
####################################################################################################################

#####Make sure cif files in correct format################
dos2unix *.cif 
#####Get the input file list##############################
FILES=`ls *.cif`

#####Read box,angle,atoms,coordinates from cif files######
 for f in $FILES
 do
 #####clean of previously dumped files#####
 rm head_${f}.txt
 rm ${f}.dat ${f}.dat.* 
 rm *_Charge_*${f}*
 rm tail_1.txt
 rm ${f}.checkfloat.txt
 rm new_tail.txt
 rm FINAL_${f}*.cif
 rm FINAL_${f}*.cif_ratio
 
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

 cp secondlayer.f90 temp.f90

 sed -i 's/FILES/'${f}'/g' temp.f90
 ####Build the connectivity file##########
 gfortran temp.f90 -o temp
 ./temp
 rm temp.f90 temp
########Combine connectivity information#######

file1=${f}.dat.0nd.txt
file2=${f}.dat.1st.txt
file3=${f}.dat.2nd.txt
file4=${f}.dat.total.txt
sed -i "s/null //g" ./$file3

for InputName in $file2 $file3
do
cp sort_colum.py sort_temp.py
sed -i "s/InputName/$InputName/g" sort_temp.py
#python sort_temp.py
python3 sort_temp.py
done
sed -i 's/ //g' $file1
sed -i 's/ //g' $file2
sed -i 's/ //g' $file3

paste $file1 $file2 $file3  > $file4


######define the database#########
data1=Database_2nd.txt
data2=Database_1st.txt
data3=Database_0th.txt

######search the 2nd layer CBAC charge##########
awk -F' ' 'NR == FNR { h[$1,$2,$3]=$4; next }; h[$1,$2,$3] {print h[$1,$2,$3];next} {print $1,$2} ' $data1 $file4 > Sec_Charge_${f}.txt
#######serach the 1st layer CABC charge##########
awk -F' ' 'NR == FNR { h[$1,$2]=$3; next }; h[$1,$2] {print h[$1,$2];next} {print $1} ' $data2 Sec_Charge_${f}.txt > First_Charge_${f}.txt
#######search the 0th layer CBAC charge##########
awk -F' ' 'NR == FNR { h[$1]=$2; next }; h[$1] {print h[$1];next} {print $1} ' $data3 First_Charge_${f}.txt > Zero_Charge_${f}.txt       
#######record the ratio of 2nd CBAC assignment ratio######
Num_atom=`awk 'END {print NR}' Sec_Charge_${f}.txt`
Num_assign=`awk '$1==$1+0 {print}' Sec_Charge_${f}.txt | wc -l` 
Mat_Ratio=`awk 'BEGIN {print '$Num_assign'/'$Num_atom'}'`
echo $Mat_Ratio $Num_Mat > ${f}_ratio.txt

#####get the head txt of cif file######
num=`awk -F"_" '$2=="atom"{print NR-1;exit}' $f`
head -${num} ${f} > head_${f}.txt
length=$(awk 'END {print NR}' ${f})
######get the coordinates section of cif file #######
awk 'NR >= '$crow' {print $'$nlabel'" "$'$nsymbol'" "$'$nfrax'" "$'$nfray'" "$'$nfraz'" "}' ${f} > tail_1.txt
########
#######assign charges to the elements not shown in current database#####
Sum_Charge=$(awk -F, '$1+0 == $1 {sum += $1} END {print sum}' Zero_Charge_${f}.txt)
Num_unknown=$(awk -F, '$1+0 != $1 {print}' Zero_Charge_${f}.txt | wc -l)
if [ $Num_unknown -gt 0 ]; then
Unknown_Charge=$(awk 'BEGIN {print (0 - '$Sum_Charge')/'$Num_unknown'}')
echo $f >> unkonwn_element.list
awk -F, '$1+0 != $1{print '$Unknown_Charge';next}1' Zero_Charge_${f}.txt > Final_Charge_${f}.txt  
else
mv Zero_Charge_${f}.txt Final_Charge_${f}.txt 
fi 
##########neturalize the total charge by weighted##############
Abs_sum=`awk 'function abs(v) {return v < 0 ? -v : v} {SUM += abs($1)} END {print SUM}' Final_Charge_${f}.txt`
Real_sum=`awk '{SUM += $1} END {print SUM}' Final_Charge_${f}.txt`
awk -v Abs="$Abs_sum" -v Real="$Real_sum" 'function abs(v) {return v < 0 ? -v : v} {print ($1-(Real*abs($1)/Abs))}' Final_Charge_${f}.txt > Final_Charge_Zeroed_${f}.txt
awk -v Abs="$Abs_sum" -v Real="$Real_sum" 'function abs(v) {return v < 0 ? -v : v} {print (Real*abs($1)/Abs)}' Final_Charge_${f}.txt > Shift_Charge_${f}.txt  
paste -d ' ' tail_1.txt Final_Charge_Zeroed_${f}.txt > new_tail.txt
###########record the total charge after neturalization######
echo $f >> zero_check.txt
awk '{ SUM += $1} END {print SUM}' Final_Charge_Zeroed_${f}.txt >> zero_check.txt
###########make new cif file######################
echo " _atom_site_label" >> head_${f}.txt
echo " _atom_site_type_symbol" >> head_${f}.txt
echo " _atom_site_fract_x" >> head_${f}.txt
echo " _atom_site_fract_y" >> head_${f}.txt
echo " _atom_site_fract_z" >> head_${f}.txt
echo " _atom_site_charge" >> head_${f}.txt
cat head_${f}.txt new_tail.txt > FINAL_${f}
##########optional check if there is floating ions or solvent#########
file5=${f}.checkfloat.txt
awk -vf="$f" 'NR==1 && $1 {print f}' $file5 >> float.list
#########clean the temp files, uncomment if need charge assignment details###########
rm head_${f}.txt
rm ${f}.dat ${f}.dat.* 
rm *_Charge_*${f}*
rm tail_1.txt
rm ${f}.checkfloat.txt
rm new_tail.txt
rm FINAL_${f}*.cif_ratio
rm sort_temp.py
done


 