#!/bin/csh
#

set tsurf=0
set hwind=0
set omega=0
set EXP=230

while ( $EXP <= 241 )
 while ( $tsurf <= 1 )
  while ( $hwind <= 1 )
   while ( $omega <= 1 )
    echo $tsurf $hwind $omega $EXP
    ./run.greb.hydro.csh $tsurf $hwind $omega $EXP
    @ omega++
   end
   set omega=0
   @ hwind++
  end
  set hwind=0
  @ tsurf++
 end
 set tsurf=0
 @ EXP++
 if ( $EXP == 231 ) set EXP=240
end

echo ' '
echo 'Convert output files to netcdf?'
cd output
sh ctl2nc.sh

exit
