#!/usr/bin/env python3 

pressureLeft = 10.
pressureRight = 9.96

nx = 20 # 100
ny = 4 # 20

lx = 10.
ly = 2.

#assert ( abs( ( lx / nx ) - (ly / ny ) ) < 1e-14 )
assert nx % 2 == 0 and ny % 2 == 0, "Number of grid cells has to be even in both dimensions!"

f = open( "backwardfacingstep_{}x{}.geom".format(nx, ny), "w" )

f.write( "physicalSizeX = {}\n".format( lx ) )
f.write( "physicalSizeY = {}\n".format( ly ) )

f.write( "nCellsX = {}\n".format( nx ) )
f.write( "nCellsY = {}\n".format( ny ) )


f.write( "\nMesh =\n" )

#Left
f.write( "NSW;TD:0,".format( pressureLeft ) )
# Top
for i in range(1,nx+1):
  f.write( "NSW;TD:0," )
# Right
f.write( "PR:{};TD:0\n".format( pressureRight ) )

for j in range(1,ny+1):
  #Left
  if ( j > ny / 2):
    f.write( "NSW;TD:0," )
  else:
    f.write( "PR:{};TD:0,".format( pressureLeft ) )
    
  #Fluid domain
  for i in range(1,nx+1):
    if ( i <= ny / 2 and j > ny / 2):
      f.write( "S," )    
    else:
      f.write( "F," )
  #Right
  f.write( "PR:{};TD:0\n".format( pressureRight ) )
  
#Left
f.write( "NSW;TD:0,".format( pressureLeft ) )
#Bottom
for i in range(1,nx+1):
  f.write( "NSW;TD:0," )
#Right
f.write( "PR:{};TD:0\n".format( pressureRight ) )


f.close()

