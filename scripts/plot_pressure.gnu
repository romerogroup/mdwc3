set nokey
plot "pressure_volume.mdout" u ($4):($5) with lines
set title "pressure vs time"
set term png
set output "pressure.png"
replot
