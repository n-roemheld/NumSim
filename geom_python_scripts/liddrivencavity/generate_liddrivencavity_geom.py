#!/usr/bin/env python3 

velXTop=1.0
temperatureTop=0.0
temperatureBottom=0.0

nx = 20
ny = 20 

lx = 2.
ly = 2.

f = open( "liddriven_cavity_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )
f.write( "NSW;TN:0," )
for i in range(1,nx+1):
  f.write( "IN:1:0;TD:{},".format(temperatureTop) )
f.write( "NSW;TN:0\n" )

for j in range(1,ny+1):
  f.write( "NSW;TN:0," )
  for i in range(1,nx+1):
    f.write( "F,".format(temperatureTop) )
  f.write( "NSW;TN:0\n" )
f.write( "NSW;TN:0," )
for i in range(1,nx+1):
  f.write( "NSW;TD:{},".format(temperatureBottom) )
f.write( "NSW;TN:0\n" )


f.close()

