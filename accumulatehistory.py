import numpy

resultfileList = [
'./workdir/Patient0002/000/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0002/001/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0003/002/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0004/003/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0005/004/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0006/005/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0004/006/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0006/007/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0005/008/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0006/009/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0002/010/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0007/011/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0005/012/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0003/013/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0008/014/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0007/015/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0008/016/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0002/017/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0003/018/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0006/019/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0002/020/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0002/021/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0002/022/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0006/023/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0006/024/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0003/025/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0007/026/opt/dakota_q_newton_heating.in.tabular.dat',
'./workdir/Patient0008/027/opt/dakota_q_newton_heating.in.tabular.dat']

with file('datasummary.txt' , 'w') as fileHandle: 
  # write header
  fileHandle.write("iddata,mu_eff,obj\n")
  # loop over files and extract optimal value
  for filename in resultfileList:
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
