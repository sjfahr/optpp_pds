import numpy
import os

resultfileList = [
'./workdir/Study0030/495/',
'./workdir/Study0030/497/',
'./workdir/Study0030/491/',
'./workdir/Study0030/496/',
'./workdir/Study0030/490/',
'./workdir/Study0017/378/',
'./workdir/Study0025/438/',
'./workdir/Study0025/435/',
'./workdir/Study0025/440/',
'./workdir/Study0025/436/',
'./workdir/Study0028/466/',
'./workdir/Study0028/468/',
'./workdir/Study0028/471/',
'./workdir/Study0026/447/',
'./workdir/Study0026/457/',
'./workdir/Study0026/455/',
'./workdir/Study0026/453/',
'./workdir/Study0026/450/',
'./workdir/Study0026/451/',
'./workdir/Study0022/418/',
'./workdir/Study0022/417/',
'./workdir/Study0021/409/',
'./workdir/Study0021/414/',
'./workdir/Study0021/415/',
]

with file('datasummary.txt' , 'w') as fileHandle: 
  # write header
  fileHandle.write("iddata,mu_eff,obj\n")
  # loop over files and extract optimal value
  opttype = 'heating'
  for filenamebase in resultfileList:
    filename = '%s/opt/dakota_q_newton_%s.in.tabular.dat' % (filenamebase,opttype)
    dataid = int(filename.split('/')[3])
    mu_eff = numpy.loadtxt(filename,skiprows=1,usecols=(1,))
    obj_fn = numpy.loadtxt(filename,skiprows=1,usecols=(2,))
    if( len(obj_fn.shape) > 0 ) :
      if( obj_fn.shape[0]  > 1 ) :
        idmin = numpy.argmin(obj_fn)
        print idmin, mu_eff[idmin], obj_fn[idmin]
        mu_effopt =  mu_eff[idmin]
        minobjval =  obj_fn[idmin]
    else:
      print "single entry ??", mu_eff, obj_fn
      mu_effopt =  mu_eff
      minobjval =  obj_fn 
    #dataarray = numpy.loadtxt(filename,skiprows=1,usecols=(0,1,2,3,4,6)
    fileHandle.write("%05d,%12.5e,%12.5e\n" %(dataid,mu_effopt,minobjval))
    # FIXME
    runcmd = "vglrun python ./brainsearch.py --param_file  %s/opt/optpp_pds.%s.in.%d %s/opt/optpp_pds.%s.out.%d --vis_out" % (filenamebase,opttype,idmin,filenamebase,opttype,idmin)
    print runcmd
    #os.system( runcmd )
