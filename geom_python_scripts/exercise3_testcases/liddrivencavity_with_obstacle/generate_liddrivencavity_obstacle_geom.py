#!/usr/bin/env python3 

velXTop=1.0
temperatureTop=1.0
temperatureBottom=0.0

nx = 20
ny = 20

lx = 1.
ly = 1.

obstacle_thickness = 3

assert( nx % 2 == 0 )
assert( ny % 2 == 0 )

gap_x = (nx - 2*obstacle_thickness)/2
gap_y = (ny - 2*obstacle_thickness)/2

assert( gap_x > 1 )
assert( gap_y > 1 )

print( gap_x, gap_y ) 

f = open( "liddriven_cavity_obstacle_{}_{}x{}.geom".format(obstacle_thickness, nx, ny), "w" )

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
    if ( (i > gap_x and i < (nx - gap_x + 1e-13) ) and ( j > gap_x and j < (ny - gap_y + 1e-13) ) ):
#      print( "i: {}, gap_x: {}, nx-gap_x: {}".format(i, gap_x, nx-gap_x) )
#      print( "j: {}, gap_y: {}, ny-gap_y: {}".format(j, gap_y, ny-gap_y) )
      f.write( "S," )
    else:
      f.write( "F," )
    
  f.write( "NSW;TN:0\n" )
f.write( "NSW;TN:0," )
for i in range(1,nx+1):
  f.write( "NSW;TD:{},".format(temperatureBottom) )
f.write( "NSW;TN:0\n" )


f.close()

