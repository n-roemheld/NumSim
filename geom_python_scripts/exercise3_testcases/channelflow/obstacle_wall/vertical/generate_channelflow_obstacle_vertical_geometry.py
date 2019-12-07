#!/usr/bin/env python3 

pressureBottom = 5.0
pressureTop = -5.0

nx = 20
ny = 100

lx = 2.
ly = 10.

assert(nx>=6)

f = open( "channelflow_obstacle_vertical_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx + 4 * lx / nx ) )
f.write( "physicalSizeY = {}\n".format( ly  ) )

f.write( "nCellsX = {}\n".format( nx+4 ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )

#Left
f.write( "NSW;TD:0," )
# Top
for i in range(-1,nx+3):
  f.write( "PR:{};TN:0,".format( pressureTop ) )
# Right
f.write( "NSW;TN:0\n" )

for j in range(1,ny+1):
  #Left
  f.write( "NSW;TD:0," )
  #Fluid domain
  for i in range(-1,nx+3):
    if ( i < 1 or i > nx ):
      f.write( "S," )
    else:
      f.write( "F," )
  #Right
  f.write( "NSW;TN:0\n" )
  
#Left
f.write( "NSW;TD:0," )
#Bottom
for i in range(-1,nx+3):
  f.write( "PR:{};TN:0,".format( pressureBottom ) )
#Right
f.write( "NSW;TN:0\n" )


f.close()

