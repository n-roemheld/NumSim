#!/usr/bin/env python3 

from enum import Enum

class BoundaryLocation(Enum):
  IsOnBoundary = 1
  Corner = 2
  NotOnBoundary = 3  
  
temperatureBottom=1.0

xOrigin=0.0
yOrigin=0.0

nx = 35
ny = 35

lx = 1.05
ly = 1.05

dx = lx / nx
dy = ly / ny

solidThickness = 0.15
fluidThickness = 0.75

#print( "dx: {}, dy:{}".format(dx,dy) )


#Computation of coordinates is shifted by -1 due to the cells indicating
#boundary conditions/ghost cells
def getX( i ):
  return xOrigin + 0.5 * dx + (i-1) * dx

def getY( j ):
  return yOrigin + 0.5 * dy + (j-1) * dy
  
def isOnOuterBoundary( i, j ):
  if ( i > 0 and i < nx+1 and j > 0 and j < ny+1 ):
    return False
  else:
    return True

def isOnLeftOuterDomain( i, j):
  return (i == 0)

def isOnRightOuterDomain( i, j):
  return (i == nx+1)

def isOnBottomOuterDomain( i, j):
  return (j == nx+1 and i != 0 and i != nx+1)
  
def isOnTopOuterDomain( i, j):
  return (j == 0 and i != 0 and i != nx+1)

def isInFluidDomain( i, j ):
  
  x = getX( i )
  y = getY( j )

  if ( x > solidThickness and x < lx-solidThickness and y > solidThickness and y < ly-solidThickness):
    #print( "(x, y): {:04e}, {:04e}".format(x, y) )
    return True
  else:
    return False
    
def isOnInnerBoundary( i, j ):
 
  x = getX( i )
  y = getY( j )

  if (x > solidThickness and x < (solidThickness+dx) ):
    # Left boundary
    return True
  elif (x > fluidThickness+solidThickness-dx and x <fluidThickness+solidThickness):
    #Right boundary
    return True
  elif (y > solidThickness and y < (solidThickness+dy) ):
    #Bottom boundary
    return True
  elif (y > fluidThickness+solidThickness-dy and y < (fluidThickness+solidThickness) ):
    #Top boundary
    return True
  else:
    return False


 
def getBoundaryTypeAtPosition( i, j ):
  x = getX( i )
  y = getY( j )

  if ( isOnOuterBoundary(i, j) ):
      loc = getBoundaryLocation( i, j )


def generateMesh():

  filename = "naturalconvection_solid_{}x{}.geom".format(nx, ny) 

  f = open( filename, "w" )

  f.write( "physicalSizeX = {}\n".format( lx ) )
  f.write( "physicalSizeY = {}\n".format( ly ) )

  f.write( "xOrigin = {}\n".format( xOrigin ) )
  f.write( "yOrigin = {}\n".format( yOrigin ) )

  f.write( "nCellsX = {}\n".format( nx ) )
  f.write( "nCellsY = {}\n".format( ny ) )

  f.write( "\nMesh =\n" )

  for j in range(0,ny+2):
    for i in range(0,nx+2):
      #print( "({:02},{:02}) is OnOuterBoundary? {:<}, OnInnerBoundary? {:1}".format( i, j, isOnOuterBoundary(i, j), isOnInnerBoundary(i, j) ) )
      x = getX( i )
      y = getY( j )
  #    print( "({:02}, {:02}): {:04f}, {:04f}".format(i, j, x, y) )
      if ( isOnOuterBoundary(i, j) ):
        if ( isOnTopOuterDomain(i,j) or isOnBottomOuterDomain(i,j) ):
          f.write( "NSW;TN:0" )
        elif( isOnLeftOuterDomain(i,j) ):
          f.write( "NSW;TD:323" )
        else:
          f.write( "NSW;TD:283" )
      elif( isInFluidDomain(i,j) ):
        if ( isOnInnerBoundary(i,j) ):
          f.write( "NSW;TPD" )
        else:
          f.write( "S" )
      else:
        f.write( "F" )
        
      if ( i < nx+1 ):
        f.write( "," )
      
    f.write( "\n" )
      
  f.close()  
  
  return filename

filename = generateMesh()

print( "Mesh written to {}.".format( filename ) )
