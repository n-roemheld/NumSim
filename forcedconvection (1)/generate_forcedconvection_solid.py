#!/usr/bin/env python3 

temperatureBottom=1.0

xOrigin=0.0
yOrigin=0.0

nx = 20
ny = 5

lx = 1
ly = 0.25

dx = lx / nx
dy = ly / ny

f = open( "forcedconvection_solid_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "xOrigin = {}\n".format( xOrigin ) )
f.write( "yOrigin = {}\n".format( yOrigin ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )
f.write( "NSW;TN:0," )
for i in range(1,nx+1):
  x = xOrigin + dx * 0.5 + (i-1) * dx
  print( "Interface at x = {}".format(x) )
  f.write( "NSW;TPD,")
f.write( "NSW;TN:0\n" )

for j in range(1,ny+1):
  f.write( "NSW;TN:0," )
  for i in range(1,nx+1):
    f.write( "F," )
  f.write( "NSW;TN:0\n" )
  
f.write( "NSW;TN:0," )
for i in range(1,nx+1):
  f.write( "NSW;TD:{},".format(temperatureBottom) )


f.write( "NSW;TN:0\n" )


f.close()

