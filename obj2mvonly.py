#Yang Liu 2009-9-30
import os
import fileinput
import math
import pdb
import re
import getopt,sys
import copy


def readObjfilelist(filename):
    vcnt = 0

    mfile = fileinput.input(filename)
    vlist = []
    flist = []
    vtlist=[]
    vnlist=[]

    mtllist=[]
    curmtl=-1

    vtflist = []
    for line in mfile:
        ar = line.split()
        if len(ar)>0:

          if ar[0]=='v':
#          pdb.set_trace()
            vid = vcnt+1
            x = float(ar[1])
            y = float(ar[2])
            z = float(ar[3])

            vinst=[x, y, z]     #vcnt starts from 0, serve as canonical index for vertices
            vlist.append(vinst)

          elif ar[0]=='f':
            if len(ar)>4:
              print('may have face with more than 3 vertices')

            if ar[1].rfind('/')>0 and ar[2].rfind('/')>0 and ar[3].rfind('/')>0 :
                v1id = int(ar[1].split('/')[0])
                v2id = int(ar[2].split('/')[0])
                v3id = int(ar[3].split('/')[0])
            else:
                v1id = int(ar[1])
                v2id = int(ar[2])
                v3id = int(ar[3])

#            pdb.set_trace()
            finst=[v1id, v2id, v3id, curmtl]
            flist.append(finst)

#           deal with texture id
            # vertex index/ texture index/ normal index
            if ar[1].rfind('/')>0 and ar[2].rfind('/')>0 and ar[3].rfind('/')>0 :
                if len(ar[1].split('/')[1]) > 0 :
                    v1vtindex = int(ar[1].split('/')[1])
                else :
                    v1vtindex = -1
                if len(ar[2].split('/')[1]) > 0 :
                    v2vtindex = int(ar[2].split('/')[1])
                else :
                    v2vtindex = -1
                if len(ar[3].split('/')[1]) > 0 :
                    v3vtindex = int(ar[3].split('/')[1])
                else :
                    v3vtindex = -1

                vtinst = [v1vtindex, v2vtindex, v3vtindex]
                vtflist.append(vtinst)

          elif ar[0]=='vt':
          # texture data
            if len(ar)>2:
              # 2d texture
                tx=float(ar[1])
                ty=float(ar[2])
            else :
              # 1d texture
                tx=float(ar[1])
                ty=0
            vtinst=[tx, ty]
            vtlist.append(vtinst)

          elif ar[0]=='vn':
          # normal dictionary
            nx = float(ar[1])
            ny = float(ar[2])
            nz = float(ar[3])
            vninst=[nx, ny, nz]
            vnlist.append(vninst)

          elif ar[0]=='usemtl':
          #material
            mtl = ar[1]
            if not mtl in mtllist:
              mtllist.append(mtl)
              print(mtl)

            curmtl=mtllist.index(mtl)

#    pdb.set_trace()
    print("v: "+str(len(vlist)))
    print("f: "+str(len(flist)))
    return [vlist, flist, vtlist, vnlist, vtflist, mtllist]






if len(sys.argv)<3:
  print("Usage: obj2mvonly.py objname mfilename")
  sys.exit()

objname = sys.argv[1]
mname = sys.argv[2]


readout = readObjfilelist(objname)

vlist = readout[0]
flist = readout[1]


of = open(mname,'w')

nv = len(vlist)
nf = len(flist)

for i in range(0, nv):
    line = ''
    line += 'Vertex '
    line += str(i+1)
    line += ' '

    v = vlist[i]
    line += str(v[0])
    line += ' '
    line += str(v[1])
    line += ' '
    line += str(v[2])
    line += '\n'

    of.write(line)

for i in range(0, nf):
    line = ''
    line += 'Face '
    line += str(i+1)
    line += ' '

    f = flist[i]
    line += str(f[0])
    line += ' '
    line += str(f[1])
    line += ' '
    line += str(f[2])
    line += '\n'

    of.write(line)



of.close()

