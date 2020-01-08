#!/usr/bin/env python3 

velXIn=1.0
velYIn=0.0

pIn = 1.0
pOut = 0.995


xOrigin=-0.5
yOrigin=0.25

nx = 70
ny = 10

lx = 3.5
ly = 0.5

dx = lx / nx
dy = ly / ny

f = open( "forcedconvection_fluid_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "xOrigin = {}\n".format( xOrigin ) )
f.write( "yOrigin = {}\n".format( yOrigin ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )
f.write( "SLW;TN:0," )
for i in range(1,nx+1):
  f.write( "SLW;TN:0,")
f.write( "SLW;TN:0\n" )

for j in range(1,ny+1):
#  f.write( "IN:{}:{};TD:0,".format( velXIn, velYIn ) )
  f.write( "PR:{};TD:0,".format( pIn ) )
  for i in range(1,nx+1):
    f.write( "F," )
  f.write( "PR:{};TN:0\n".format( pOut ) )
#  f.write( "OUT;TN:0\n" )
  
f.write( "SLW;TN:0," )
for i in range(1,nx+1):

  x = dx * 0.5 + (i-1) * dx + xOrigin
 
  if ( x < 0.0 ):
    f.write( "SLW;TN:0," )
  elif ( x >= 0.0 and x <= 1.0 ):
    f.write( "NSW;TPN," )
    print( "Interface at x = {}".format(x) )
  else:
    f.write( "NSW;TN:0," ) 

f.write( "NSW;TN:0\n" )


f.close()

