#!/usr/bin/env python3 

pressureLeft = 5.0
pressureRight = -5.0

nx = 25 # 100
ny = 5 # 20

lx = 10.
ly = 2.

f = open( "channelflow_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )

#Left
f.write( "PR:{};TN:0,".format( pressureLeft ) )
# Top
for i in range(1,nx+1):
  f.write( "NSW;TD:0," )
# Right
f.write( "PR:{};TN:0\n".format( pressureRight ) )

for j in range(1,ny+1):
  #Left
  f.write( "PR:{};TN:0,".format( pressureLeft ) )
  #Fluid domain
  for i in range(1,nx+1):
    f.write( "F," )
  #Right
  f.write( "PR:{};TN:0\n".format( pressureRight ) )
  
#Left
f.write( "PR:{};TN:0,".format( pressureLeft ) )
#Bottom
for i in range(1,nx+1):
  f.write( "NSW;TD:0," )
#Right
f.write( "PR:{};TN:0\n".format( pressureRight ) )


f.close()

