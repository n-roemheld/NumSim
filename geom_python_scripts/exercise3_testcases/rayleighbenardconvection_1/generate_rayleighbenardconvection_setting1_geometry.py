#!/usr/bin/env python3 

temperatureTop=291.20
temperatureBottom=294.78

nx = 85
ny = 18

lx = 8.5
ly = 1.

f = open( "rayleighbenardconvection_setting1_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )
#Top
f.write( "NSW;TN:0," )
for i in range(1,nx+1):
  f.write( "NSW;TD:{},".format( temperatureTop ) )
f.write( "NSW;TN:0\n" )

# Center
for j in range(1,ny+1):
  f.write( "NSW;TN:0," )
  for i in range(1,nx+1):
    f.write( "F," )
  f.write( "NSW;TN:0\n" )
  
#Bottom
f.write( "NSW;TN:0," )
for i in range(1,nx+1):
  f.write( "NSW;TD:{},".format( temperatureBottom ) )
f.write( "NSW;TN:0\n" )


f.close()

