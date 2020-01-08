#!/usr/bin/env python3 

xOrigin=0.15
yOrigin=0.15

nx = 25
ny = 25

lx = 0.75
ly = 0.75

dx = lx / nx
dy = ly / ny

f = open( "naturalconvection_fluid_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "xOrigin = {}\n".format( xOrigin ) )
f.write( "yOrigin = {}\n".format( yOrigin ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )
f.write( "NSW;TN:0," )
for i in range(1,nx+1):
  f.write( "NSW;TPN,")
f.write( "NSW;TN:0\n" )

for j in range(1,ny+1):
  f.write( "NSW;TPN," )
  for i in range(1,nx+1):
    f.write( "F," )
  f.write( "NSW;TPN\n" )
  
f.write( "NSW;TN:0," )
for i in range(1,nx+1):
  f.write( "NSW;TPN," ) 

f.write( "NSW;TN:0\n" )


f.close()

