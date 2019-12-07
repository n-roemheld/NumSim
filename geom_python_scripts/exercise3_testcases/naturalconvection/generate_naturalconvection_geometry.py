#!/usr/bin/env python3 

temperatureLeft=1.0
temperatureRight=0.0

nx = 50
ny = 50

lx = 1.
ly = 1.

f = open( "naturalconvection_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )
#Top
f.write( "NSW;TD:{},".format(temperatureLeft) )
for i in range(1,nx+1):
  f.write( "NSW;TN:0," )
f.write( "NSW;TD:{}\n".format(temperatureRight) )

# Center
for j in range(1,ny+1):
  f.write( "NSW;TD:{},".format(temperatureLeft) )
  for i in range(1,nx+1):
    f.write( "F," )
  f.write( "NSW;TD:{}\n".format(temperatureRight) )
  
#Bottom
f.write( "NSW;TD:{},".format(temperatureLeft) )
for i in range(1,nx+1):
  f.write( "NSW;TN:0," )
f.write( "NSW;TD:{}\n".format(temperatureRight) )


f.close()

