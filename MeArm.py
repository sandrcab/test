#!/usr/bin/python

from servos import ServoSB
from links import MultiLink
import numpy as np
import sympy as sp
import frames
from IPython.display import display

import matplotlib.pyplot as plt

def _s(s):
  return frames.SYMB(s)

class MeArm(object):
  def __init__(self):
    # TODO: connect the servos to their pin and select the range of the angles
    DUMMY_VAL = 1,0,180 # Remove this!
    self.servos = {"base":ServoSB( 7,-90,90 ),
                  "shoulder":ServoSB( 11,0,90),
                  "elbow":ServoSB(13,0,45),
                  "grip":ServoSB( 15,0,90 )}
    
    # Create the body using multilink code and DH parameters
    _l = [.20, #Origin to base 
          .20, #Base to Shoulder
         .80, #Shoulder to elbow length
         .80, #Elbow to wrist length
         .68] #Wrist to hand length
    self._joints = {_s("j1"):0.0, 
                    _s("j2"):0, 
                    _s("j3"):45.0, 
                    _s("j4"):0.0}
    self._body = MultiLink()
    
    _dh = {"d":[.2, .2, .0, .0],
         "theta":["s:j1","s:j2","s:j3","s:j4"],
         "r":[.0, .0, .8, .8],
         "alpha":[0.0, -90.0, .0, .0]
        }
    self._body.fromDH(_dh)
    self._body.compose(*self._joints.keys())
    self._J = self._body.J.copy()
    print self._body.tt

  def base(self, angle):
    self.servos["base"].set_angle(angle)
    print "Moving base to:", angle
      
  def shoulder(self, angle):
    self.servos["shoulder"].set_angle(angle)
    print "Moving shoulder to:", angle        
  
  def elbow(self, angle):
    self.servos["elbow"].set_angle(angle)
    print "Moving elbow to:", angle
  
  def gripper(self, angle):
    self.servos["gripper"].set_angle(angle)
    print "Moving gripper to:", angle   
  
  def openGripper(self):
    self.gripper(90)
    print "Grip open"
  
  def closeGripper(self):
    self.gripper(0)
    print "Grip close"
  
  def clap(self):
    self.openGripper(); self.closeGripper()
    self.openGripper(); self.closeGripper()
    self.openGripper()        
  
  def readJoints(self):
    # TODO: update joints values from motors
    self._body.bindSymbols(self._joints)
    self._joints[_s("j4")] = 360 - (self._joints[_s("j2")]+self._joints[_s("j3")])
    Q = sp.Matrix( [self._joints[_s("j1")],
                    self._joints[_s("j2")], 
                    self._joints[_s("j3")],
                    self._joints[_s("j4")]] )
    return Q
  
  def getJacobian(self):
    self.readJoints()
    return self._body.getJacobian()
  
  def moveJoints(self, Q):
    self._joints[_s("j1")] = Q[0]
    self._joints[_s("j2")] = Q[1] 
    self._joints[_s("j3")] = Q[2]
    self._joints[_s("j4")] = Q[3]
    
    self.base(Q[0])
    self.shoulder(Q[1])
    self.elbow(Q[2])
    #self.clap(Q[3])
      
  def gotoPoint(self, tx, ty, tz, steps=1):
    """ Simple control loop that uses the Jacobian"""
    print "Moving to:",(tx,ty,tz)
    for i in range (steps):
      print "****** Iteration #%d:"%i
      # Read current joint status
      Q = self.readJoints()
      display(Q)
      
      # Draw current robot position
      
      # Get position in 3D space (DK)
      pose = self._body.jointsToPosition(Q)
      display(pose)
      # Get distance from target
      E = sp.Matrix([[tx - pose[0]],
                     [ty - pose[1]],
                     [tz - pose[2]],
                     [1]])
      display(E)
      print "Es:",E.shape
      
      #Compute jacobian
      J = self.getJacobian()
      display(J)
      print "Js:",J.shape
      
      # Calculate joint change
      dQ =  J.T * E*0.005
      display(dQ)
      # and use it to update current values
      Q +=  dQ
      self.moveJoints(Q)
    
if __name__ == '__main__':
  arm = MeArm()
  #arm.clap()
  #Go up and left to grab something
  arm.gotoPoint(-1,1.9,0.7); 
  #arm.closeGripper();
  #Go down, forward and right to drop it
  #arm.gotoPoint(.70,2.00,.10);
  #arm.openGripper();
  #Back to start position
  #arm.gotoPoint(0,1.00,.50);
