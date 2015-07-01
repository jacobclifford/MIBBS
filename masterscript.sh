#!/bin/bash

#this script first runs the clustering simulation that resides in the folder 'cluster', once the clustering simulation is complete and the sequece files (fasta's of each cluster of binding sites) have been created the sequence files are copied to six subprograms  (which are numbered below as second subprogram to seventh subprogram, where the 'first subprogram' is the clustering experiment).

# $cluster folder runs the clustering alogorithm to decompose cbe.fa into dc and du and cb fasta files for various cluster parameters (spacer window cutoffs, Twist motif..),

# branch folders: Logodds of predicted sites, and Ranksum ($ranksumproj), Mutual Information I(I,O) ($rocproj),  Mutual Information I(C,P) ($miclasses), ROC ($rocproj, logos ($data_dir, $data_dir2)
 

logosproj=~/Desktop/MIBBS/
cluster=~/Desktop/MIBBS/cluster

logodds=~/Desktop/MIBBS/logodds
#ranksum is also computed in this folder
ranksumproj=~/Desktop/MIBBS/logodds

rocenhancer=~/Desktop/MIBBS/rocenhancer
miclasses=~/Desktop/MIBBS/miclasses
rocproj=~/Desktop/MIBBS/rocproj

data_dir=~/Desktop/MIBBS/data_dir
data_dir2=~/Desktop/MIBBS/data_dir2
#first subprogram (clustering Dorsal loci based on Twist motif and spacer window)
cd $cluster
rm format.tex
 ./seq2exp -bs fas.txt -bi topbot2Int.txt -s efastanotw12.txt  -p synmyout -e expre12.tab -m factordts.wtmx -f factorexpdts2.tab -na 0 -i factorinfodts.txt -o BINS -c coopdt.txt -fo out.txt -oo corr -dc  cbe.fa -du cbe.fa -tw Ebox.txt -sa NEEandREPseq.txt   -pva .0001 -cbo cb_11b.txt   -l 9 -nr 1000 -ip 1  -nu 5 -w 4 -m2 twmotifs842014.wtmx


j=0
q=1  # this is the number of motifs of twist in cluster folder (still need to add this index to files.. in roc ,ranksum etc..)
echo " hi "
i=1
k=2 #10  # this is the nupdates of the spacer in cluster folder
while [ $j -lt $q ]; do
	while [ $i -lt $k ]; do
#file names
cbpartition=$cluster/coreOPT"$j"cb$i.txt
dcpartition=$cluster/coreOPT"$j"dc$i.txt
dupartition=$cluster/coreOPT"$j""du"$i.txt

cd $cluster

dccoreshif=$cluster/coredc"$i"OPT1.txt
ducoreshif=$cluster/coredu"$i"OPT1.txt
cbcoreshif=$cluster/corecb"$i"OPT1.txt  #  i is the update of sapcer, 1 is the shift, dc du cb indicate the type of site..
cp $cbcoreshif $data_dir/
cp $dccoreshif $data_dir/
cp $ducoreshif $data_dir/

cp $cbcoreshif $data_dir2/
cp $dccoreshif $data_dir2/
cp $ducoreshif $data_dir2/

cp $dcpartition $data_dir2/
cp $dupartition $data_dir2/
cp $cbpartition $data_dir2/

cp $dcpartition $data_dir/
cp $dupartition $data_dir/
cp $cbpartition $data_dir/

cp $dcpartition $ranksumproj/
cp $dupartition $ranksumproj/
cp $cbpartition $ranksumproj/
cp $dcpartition $rocproj/
cp $dupartition $rocproj/
cp $cbpartition $rocproj/



cp $dcpartition $logodds/
cp $dupartition $logodds/
cp $cbpartition $logodds/

cp $dcpartition $rocenhancer/
cp $dupartition $rocenhancer/
cp $cbpartition $rocenhancer/

cp $dcpartition $miclasses/
cp $dupartition $miclasses/
cp $cbpartition $miclasses/


cp cbe.fa $miclasses/
cp cbe.fa $rocenhancer/

cp cbe.fa $rocproj/
cp cbe.fa $ranksumproj/
cp cbe.fa $logodds/
cp cbe.fa $rocenhancer/
#: '

#second subprogram (logodds of predicted sites adjaceny to Twist)
cd $logodds
rm format.tex
rm beautydcdu.txt
 ./seq2exp -bs fas.txt -bi topbot2Int.txt -s efastanotw12.txt  -p synmyout -e expre12.tab -m factordts.wtmx -f factorexpdts2.tab -na 0 -i factorinfodts.txt -o BINS -c coopdt.txt -fo out.txt -oo corr -du coreOPT"$j""du"$i.txt  -dc coreOPT"$j"dc$i.txt -tw Ebox.txt -sa  NEEandREPseq.txt   -pva .0005 -cbo coreOPT"$j"cb$i.txt   -l 9 -nr 8 -ip 1  -nu 2 -w 1 -m2 twmotifs.wtmx -fi 20
matlab -nodesktop -r "logoddsPredpval; exit"

cp logoddsPredpvaAlltw.eps $logosproj

#third subprogram (ROC with CRMs as negatives (where Dorsal TFBS loci are masked))
cd $rocenhancer
rm format.tex
rm beautydcdu.txt
rm TPDC.txt
rm TPDU.txt
rm TPCB.txt
rm FPDC.txt
rm FPDU.txt
rm FPCB.txt
 ./seq2exp -bs fas.txt -bi topbot2Int.txt -s efastanotw12.txt  -p synmyout -e expre12.tab -m factordts.wtmx -f factorexpdts2.tab -na 0 -i factorinfodts.txt -o BINS -c coopdt.txt -fo out.txt -oo corr -du coreOPT"$j""du"$i.txt -dc coreOPT"$j"dc$i.txt -tw Eboxn.txt -sa  cbe.fa  -pva .0005 -cbo coreOPT"$j"cb$i.txt    -l 9 -nr 8 -ip 1  -nu 2 -w 1 -m2 twmotifsn.wtmx -fi 20 -saFP NEEandREPseq.txt

matlab -nodesktop -r "rocdiff; exit"
# roc of OR gate vs roc of CB detector
cp rocenha.eps $logosproj/rocenha$i.eps

# mi I(I;O) of Input class (Dorsal or not Dorsal) and the Output predictions of OR gate for a given energy cutoff of class type; and similarly for CB detector
matlab -nodesktop -r "channelinputoutput; exit"
cp mienhancers.eps $logosproj/mienhancers$i.eps

#fourth subprogram (compute I(C,P), for class types (DC, DU) and predictions of class type for each detector DC and DU for a given energy cutoff.  In a sense, this program answers the question: can DC predict DC sites as positives, while also predicting DU sites as negatives? 

cd $miclasses
rm format.tex
rm beautydcdu.txt
rm TPDC.txt
rm TPDU.txt
rm TPCB.txt
rm FPDC.txt
rm FPDU.txt
rm FPCB.txt
 ./seq2exp -bs fas.txt -bi topbot2Int.txt -s efastanotw12.txt  -p synmyout -e expre12.tab -m factordts.wtmx -f factorexpdts2.tab -na 0 -i factorinfodts.txt -o BINS -c coopdt.txt -fo out.txt -oo corr -du coreOPT"$j""du"$i.txt  -dc coreOPT"$j"dc$i.txt -tw Eboxn.txt -sa  cbe.fa  -pva .0005 -cbo coreOPT"$j"cb$i.txt    -l 9 -nr 8 -ip 1  -nu 2 -w 1 -m2 twmotifsn.wtmx -fi 20 -saFP NEEandREPseq.txt

matlab -nodesktop -r "cdfwchannelfptw; exit"
cp miclasses.eps $logosproj/miclasses$i.eps

# the folder 'rocproj' generates a figure (rocp$i.eps) useful for analyis of I(C,P) behavior.  By generating an ROC where DC sites are positives for DC PWM (detector), and DU sites are negatives for DC PWM, similarly DU sites are postives for DU PWM, and DC sites are negatives for DU PWM.
cd $rocproj
rm format.tex
rm beautydcdu.txt
rm TPDC.txt
rm TPDU.txt
rm TPCB.txt
rm FPDC.txt
rm FPDU.txt
rm FPCB.txt
echo " echo rocproj"
#exp1: using dc as pos, du as neg
 ./seq2exp -bs fas.txt -bi topbot2Int.txt -s efastanotw12.txt  -p synmyout -e expre12.tab -m factordts.wtmx -f factorexpdts2.tab -na 0 -i factorinfodts.txt -o BINS -c coopdt.txt -fo out.txt -oo corr -du coreOPT"$j""du"$i.txt  -dc coreOPT"$j"dc$i.txt -tw Eboxn.txt -sa  coreOPT"$j"dc$i.txt  -pva .0005 -cbo coreOPT"$j"cb$i.txt   -l 9 -nr 8 -ip 1  -nu 2 -w 1 -m2 twmotifsn.wtmx -fi 20 -saFP coreOPT"$j""du"$i.txt 
cp TPDC.txt TPDCt.txt
cp FPDC.txt FPDCt.txt
rm format.tex
rm beautydcdu.txt
rm TPDC.txt
rm TPDU.txt
rm TPCB.txt
rm FPDC.txt
rm FPDU.txt
rm FPCB.txt
#exp2: using dc as neg, du as pos
echo "echo rocproj2"
./seq2exp -bs fas.txt -bi topbot2Int.txt -s efastanotw12.txt  -p synmyout -e expre12.tab -m factordts.wtmx -f factorexpdts2.tab -na 0 -i factorinfodts.txt -o BINS -c coopdt.txt -fo out.txt -oo corr -du coreOPT"$j""du"$i.txt  -dc coreOPT"$j"dc$i.txt -tw Eboxn.txt -sa  coreOPT"$j""du"$i.txt    -pva .0005 -cbo coreOPT"$j"cb$i.txt  -l 9 -nr 8 -ip 1  -nu 2 -w 1 -m2 twmotifsn.wtmx -fi 20 -saFP coreOPT"$j"dc$i.txt
cp TPDCt.txt TPDC.txt  #overwrite the DC rates, so as to use exp1 run
cp FPDCt.txt  FPDC.txt
matlab -nodesktop -r "roc; exit" 
cp rocp.eps $logosproj/rocp$i.eps

# fifth subprogram (ranksum test, with pvalue sampling distribution)
cd $ranksumproj
echo "echo ranksumproj"

./seq2exp -bs fas.txt -bi topbot2Int.txt -s efastanotw12.txt  -p synmyout -e expre12.tab -m factordts.wtmx -f factorexpdts2.tab -na 0 -i factorinfodts.txt -o BINS -c coopdt.txt -fo out.txt -oo corr -du coreOPT"$j""du"$i.txt  -dc coreOPT"$j"dc$i.txt  -tw Ebox.txt -sa  NEEandREPseq.txt   -pva .0005 -cbo coreOPT"$j"cb$i.txt   -l 9 -nr 8 -ip 1  -nu 2 -w 1 -m2 twmotifs.wtmx -fi 20

matlab -nodesktop -r "usw; exit"
cp rs.eps $logosproj/rs$i.eps


# sixth subprogram
# this foder flips all sites in an alignment to A-rich seq, and then makes a pwm table from the flipped strand

cd $data_dir 

./seq2exp  -s coreOPT"$j""du"$i.txt -dc coreOPT"$j""du"$i.txt -et 26 -l 9 -nr 1 -ip 0 -nu 1 -w 1 -du du.txt -motifclass 0 -psed .1

s=5
	x=$(awk -v tit=$s '{print $tit}' info.txt )  #dcinfo.txt)  
	r=$(cat info.txt)
	#train_seq_file=DU$i.txt        #be
#echo $i > t.txt
matlab -nodesktop -r "pwm;  exit;"
cp pwm.eps $logosproj/pwmdu$i.eps
cp rcRead.txt $data_dir2/

# seventh subprogram
cd $data_dir2   # this folder makes a weblogo

./weblogo --format EPS --weight .1 --color-scheme classic --errorbars NO --title B --fineprint "" --first-index 0 --annotate "0,1,2,3,4,5,6,7,8" < rcRead.txt > du11m.eps
cp du11m.eps $logosproj/du11m$i.eps

# sixth subprogram is called again for 'extended binding site sequences' where one flanking base has been added to each binding site locus in the alignment.  Notice the input to the program executable 'seq2exp' under the -s parameter is coredu"$i"OPT1.txt, whereas the call to the sixth subprogram above had input for the -s parameter as coreOPT"$j""du"$i.txt (the quotes around "du" is so bash doesn't call its program 'du').  Inside the ExprPredictor.cpp library, in the "alingNEW2" function one will find the file coredu"$i"OPT1.txt is being generated inside a for loop at the end of the alignment routine for the DU motif, as this is the only place where a record of each loci's binding site start site and strand is kept, hence the file must be generated there in order to access the needed information.  If one wants two flanking bases added to each binding site, then use the file coredu"$i"OPT2.txt, and for three flanking bases use file coredu"$i"OPT3.txt, hence the number of flanking bases is in general coredu"$i"OPT$j+1$.txt... (as can be seen,the j index of the 'do while' shell script loop is not in sync with the file names, hence I have manually just placed the appropriate index (which for one flanking base is '1') in this file for now..)

cd $data_dir
./seq2exp  -s coredu"$i"OPT1.txt -dc coredu"$i"OPT1.txt -et 26 -l 9 -nr 1 -ip 0 -nu 1 -w 1 -du du.txt -motifclass 0 -psed .1

s=5
	x=$(awk -v tit=$s '{print $tit}' info.txt )  #dcinfo.txt)  
	r=$(cat info.txt)

matlab -nodesktop -r "pwm;  exit;"
cp pwm.eps $logosproj/pwm1du$i.eps
cp rcRead.txt $data_dir2/


# seventh subprogram is called again for the extended binding site sequences.; after this program call the sixth and seventh programs are called again for DC and CB motifs.
cd $data_dir2   # this folder makes a weblogo
./weblogo --format EPS --weight .1 --color-scheme classic --errorbars NO --title E --fineprint "" --first-index -1 --annotate "-1,0,1,2,3,4,5,6,7,8,9" < rcRead.txt > du11m1.eps
cp du11m1.eps $logosproj/du11m1$i.eps







cd $data_dir
./seq2exp  -s coreOPT"$j"dc$i.txt -dc coreOPT"$j"dc$i.txt -et 26 -l 9 -nr 1 -ip 0 -nu 1 -w 1 -du du.txt -motifclass 0 -psed .1

s=5
	x=$(awk -v tit=$s '{print $tit}' info.txt )  #dcinfo.txt)  
	r=$(cat info.txt)
	#train_seq_file=DU$i.txt        #be
#echo $i > t.txt
matlab -nodesktop -r "pwm;  exit;"
cp pwm.eps $logosproj/pwmdc$i.eps
cp rcRead.txt $data_dir2

cd $data_dir2

./weblogo --format EPS --weight 0 --color-scheme classic --errorbars NO --title A --fineprint "" --first-index 0 --annotate "0,1,2,3,4,5,6,7,8" < rcRead.txt > dc11m.eps
cp dc11m.eps $logosproj/dc11m$i.eps

cd $data_dir
./seq2exp  -s coredc"$i"OPT1.txt -dc coredc"$i"OPT1.txt -et 26 -l 9 -nr 1 -ip 0 -nu 1 -w 1 -du du.txt -motifclass 0 -psed .1

s=5
	x=$(awk -v tit=$s '{print $tit}' info.txt )  #dcinfo.txt)  
	r=$(cat info.txt)

matlab -nodesktop -r "pwm;  exit;"
cp pwm.eps $logosproj/pwm1dc$i.eps
cp rcRead.txt $data_dir2/
cd $data_dir2   # this folder makes a weblogo

./weblogo --format EPS --weight 0 --color-scheme classic --errorbars NO --title D --first-index -1  --fineprint "" --annotate "-1,0,1,2,3,4,5,6,7,8,9" < rcRead.txt > dc11m1.eps
cp dc11m1.eps $logosproj/dc11m1$i.eps








cd $data_dir
echo "echo align2 cb"
./seq2exp  -s coreOPT"$j"cb$i.txt -dc coreOPT"$j"cb$i.txt -et 26 -l 9 -nr 1 -ip 0 -nu 1 -w 1 -du du.txt -motifclass 0 -psed .1
train_seq_file=$data_dir/coreOPT.txt
s=5
	x=$(awk -v tit=$s '{print $tit}' info.txt )  #dcinfo.txt)  
	r=$(cat info.txt)
	#train_seq_file=DU$i.txt        #be
#echo $i > t.txt
matlab -nodesktop -r "pwm;  exit;"
cp pwm.eps $logosproj/pwmcb$i.eps
cp rcRead.txt $data_dir2/
cd $data_dir2

./weblogo --format EPS --weight 0 --color-scheme classic --errorbars NO --title C  --fineprint "" --first-index 0 --annotate "0,1,2,3,4,5,6,7,8" < rcRead.txt > cb11m.eps
cp $train_seq_file Data/CBm.txt
cp cb11m.eps $logosproj/cb11m$i.eps




cd $data_dir
./seq2exp  -s corecb"$i"OPT1.txt -dc corecb"$i"OPT1.txt -et 26 -l 9 -nr 1 -ip 0 -nu 1 -w 1 -du du.txt -motifclass 0 -psed .1

s=5
	x=$(awk -v tit=$s '{print $tit}' info.txt )  #dcinfo.txt)  
	r=$(cat info.txt)
	#train_seq_file=DU$i.txt        #be
#echo $i > t.txt
matlab -nodesktop -r "pwm;  exit;"
cp pwm.eps $logosproj/pwm1cb$i.eps
cp rcRead.txt $data_dir2/
cd $data_dir2   # this folder makes a weblogo

./weblogo --format EPS --weight 0 --color-scheme classic --errorbars NO --title F --fineprint "" --first-index -1 --annotate "-1,0,1,2,3,4,5,6,7,8,9" < rcRead.txt > cb11m1.eps
cp cb11m1.eps $logosproj/cb11m1$i.eps

	let i++
	done
	#fi
        echo "j "
        echo $j
        i=1
let j++
done

