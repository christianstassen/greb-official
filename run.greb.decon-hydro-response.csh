#                                                   Tsurf  Eva  Omega  Omegastd  Hwind
./run.greb.decon2xco2.hydro.csh 1 1 1 1 1 1 1 1 1 1 2      2    2      2         2     # Everything is on and responding (forced)           #1
./run.greb.decon2xco2.hydro.csh 1 1 1 1 1 1 1 1 1 1 2      2    3      3         3     # As #1 but only evaporation is responding (forced)  #2
./run.greb.decon2xco2.hydro.csh 1 1 1 1 1 1 1 1 1 1 2      2    3      3         2     # As #2 plus Hwind is responding                     #3
./run.greb.decon2xco2.hydro.csh 1 1 1 1 1 1 1 1 1 1 2      2    2      3         2     # As #3 plus omega is responding
./run.greb.decon2xco2.hydro.csh 1 1 1 1 1 1 1 1 1 1 2      1    2      2         2     # As #1 but evaporation is dynamic                   #5
./run.greb.decon2xco2.hydro.csh 1 1 1 1 1 1 1 1 1 1 3      1    2      2         2     # As #5 but Tsurf is PI                              #6
