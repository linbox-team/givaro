set xlabel "Integer bit size"
set ylabel "Speed (M. arith. op./s)" 
set title "MEANING on an MODEL" # tc lt  2
#set key below
set logscale y 10
#set logscale x 2
#set ytics (10,1000)
set xtics ("64" 6,"128" 7,"256" 8,"512" 9,"1024" 10,"2048" 11,"4096" 12, "8192" 13)
#set xtics tc lt 2
#set ytics tc lt 2
#set grid noxtics ytics lt 2
set grid noxtics ytics 
#set border 4095 lt 7
#set style line 5 pt 2



set terminal pdf enhanced color solid lw 2 size 6,4
set output "rint_FUNCTION.pdf"
plot [6:13] "output.rint.FUNCTION" using 7:($5) title "GMP-6" with linespoint lt 2 lc 1

set terminal pdf enhanced color solid lw 2 size 6,4
set output "rint_FUNCTION.pdf"
replot "output.rint.FUNCTION" using 7:($4) title "RecInt" with fsteps lt 3 
