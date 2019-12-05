#!/usr/bin/env python3 

pressureBottom = 5.0
pressureTop = -5.0

nx = 10
ny = 100

lx = 1.
ly = 10.

f = open( "channelflow_vertical_halfheight_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )

#Left
f.write( "NSW;TN:0," )
# Top
for i in range(1,nx+1):
  f.write( "PR:{};TD:0,".format(pressureTop) )
# Right
f.write( "SLW;TN:0\n" )

for j in range(1,ny+1):
  #Left
  f.write( "NSW;TN:0," )
  #Fluid domain
  for i in range(1,nx+1):
    f.write( "F," )
  #Right
  f.write( "SLW;TN:0\n" )
  
#Left
f.write( "NSW;TN:0," )
#Bottom
for i in range(1,nx+1):
  f.write( "PR:{};TD:0,".format(pressureBottom) )
#Right
f.write( "SLW;TN:0\n" )


f.close()

