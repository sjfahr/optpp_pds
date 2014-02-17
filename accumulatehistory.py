import numpy

resultfileList = [
'./workdir/Patient0002/000/',
'./workdir/Patient0002/001/',
'./workdir/Patient0003/002/',
'./workdir/Patient0004/003/',
'./workdir/Patient0005/004/',
'./workdir/Patient0006/005/',
'./workdir/Patient0004/006/',
'./workdir/Patient0006/007/',
'./workdir/Patient0005/008/',
'./workdir/Patient0006/009/',
'./workdir/Patient0002/010/',
'./workdir/Patient0007/011/',
'./workdir/Patient0005/012/',
'./workdir/Patient0003/013/',
'./workdir/Patient0008/014/',
'./workdir/Patient0007/015/',
'./workdir/Patient0008/016/',
'./workdir/Patient0002/017/',
'./workdir/Patient0003/018/',
'./workdir/Patient0006/019/',
'./workdir/Patient0002/020/',
'./workdir/Patient0002/021/',
'./workdir/Patient0002/022/',
'./workdir/Patient0006/023/',
'./workdir/Patient0006/024/',
'./workdir/Patient0003/025/',
'./workdir/Patient0007/026/',
'./workdir/Patient0008/027/']

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
    print "python ./brainsearch.py --param_file  %s/opt/optpp_pds.%s.in.%d %s/opt/optpp_pds.%s.out.%d --vis_out" % (filenamebase,opttype,idmin,filenamebase,opttype,idmin)
