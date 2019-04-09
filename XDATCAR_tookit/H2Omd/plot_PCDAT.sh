awk <PCDAT >PCDAT.dat '
NR==8 { pcskal=$1}
NR==9 { pcfein=$1}
NR>=13 {
 line=line+1
 if (line==257)  {
    print " "
    line=0
 }
 else
    print (line-0.5)*pcfein/pcskal,$1
}
'
cat >plotfile<<!
# set term postscript enhanced colour lw 2 "Helvetica" 20
# set output "pair_correlation.eps"
set title "pair-correlation of H2O at 2000 K"
set xlabel "r [Angstrom]"
set ylabel "g(r)"
plot [0:15] "PCDAT.dat"  w lines
!
gnuplot -persist plotfile
