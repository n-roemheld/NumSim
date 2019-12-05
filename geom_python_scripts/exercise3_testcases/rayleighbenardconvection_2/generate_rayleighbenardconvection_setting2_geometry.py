#!/usr/bin/env python3 

temperatureTop=1
temperatureBottom=2

nx = 20
ny = 20

lx = 1.
ly = 1.

f = open( "rayleighbenardconvection_setting2_{}x{}.geom".format(nx, ny), "w" )

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

