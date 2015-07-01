rm format.tex
rm beautydcdu.txt
rm TPDC.txt
rm TPDU.txt
rm TPCB.txt
rm FPDC.txt
rm FPDU.txt
rm FPCB.txt
 ./seq2exp -bs fas.txt -bi topbot2Int.txt -s efastanotw12.txt  -p synmyout -e expre12.tab -m factordts.wtmx -f factorexpdts2.tab -na 0 -i factorinfodts.txt -o BINS -c coopdt.txt -fo out.txt -oo corr -du "coreOPT0du1.txt"  -dc "coreOPT0dc1.txt" -tw Eboxn.txt -sa  cbe.fa  -pva .0005 -cbo "coreOPT0cb1.txt"   -l 9 -nr 8 -ip 1  -nu 2 -w 1 -m2 twmotifsn.wtmx -fi 20 -saFP NEEandREPseq.txt

