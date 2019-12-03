#!/usr/bin/env python3 

uIn = 1.0
vIn = 0.0
#pressureRight = -5.0

nx = 20 # 10 
ny = 5 # 10 

lx = 2.
ly = 2.

#print( "Ny: ", ny, ny % 5 )
assert (ny % 5) == 0, "Channel height must be a multiple of 5"

obstacle_height = ny / 5

f = open( "karmanvortex_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )

#Left
f.write( "IN:{}:{};TN:0,".format( uIn, vIn ) )
# Top
for i in range(1,nx+1):
  f.write( "NSW;TD:0," )
# Right
f.write( "OUT;TN:0\n" )

for j in range(1,ny+1):
  #Left
  f.write( "IN:{}:{};TN:0,".format( uIn, vIn ) )
  #Fluid domain
  for i in range(1,nx+1):
    if ( j > 2 * obstacle_height and j < 3 * obstacle_height+1 ):
      if ( j == 2 * obstacle_height+1 and ( i > 3 * obstacle_height - 2 ) and (i < 3 * obstacle_height + 1) ):
        f.write( "S," )
      elif ( j == 3 * obstacle_height and ( i > 2 * obstacle_height ) and (i < 2 * obstacle_height + 3 ) ):
        f.write( "S," )
      elif ( j > 2 * obstacle_height+1 and j < 3 * obstacle_height ):
        if ( (j-i) > 0 and (j-i) < 2 ):
          f.write( "S," )
        else:
          f.write( "F," )
      else: 
#      print("i, j: {}, {}".format(i,j) )
        f.write( "F," )
    else:
      f.write( "F," )
  #Right
  f.write( "OUT;TN:0\n" )
  
#Left
f.write( "IN:{}:{};TN:0,".format( uIn, vIn ) )
#Bottom
for i in range(1,nx+1):
  f.write( "NSW;TD:0," )
#Right
f.write( "OUT;TN:0\n" )


f.close()

