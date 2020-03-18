set nokey
plot "pressure_volume.mdout" u ($4):($6) with lines
set title "volume vs time"
set term png
set output "volume.png"
replot
