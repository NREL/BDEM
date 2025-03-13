import sys
import os
import math
def fprintf(stream, format_spec, *args):
    stream.write(format_spec % args)

f = open("particle_input.dat", "x")

N = 10

minx = -0.001
maxx = 0.0015

ypos = 0.00025
dp = ( maxx - minx )/N
dp_ext = dp*1.00
zpos = 0.002 + dp_ext
pid = 1
fprintf(f,"    %i\n",N)
# for xid in range(N):
#     for yid in range(N):
#         xpos = minx + dp/2.0 + dp*xid
#         ypos = minx + dp/2.0 + dp*yid
#         fprintf(f,"%i %f %f %f %f 2300 0 0 0 \n",pid, xpos, ypos, zpos, dp_ext)
#         pid = pid +1
for xid in range(N):

    xpos = minx + dp/2.0 + dp*xid
    if xpos < 0: 
        vx =  0.05
    else:
        vx = -0.05
    fprintf(f,"%i %f %f %f %f 2300 %f 0 0 1\n",pid, xpos, ypos, zpos, dp_ext/2.0, vx)
    pid = pid +1

f.close()
