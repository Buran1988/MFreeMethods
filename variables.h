dimension pEBC(2,100),npEBC(3,100),npNBC(3,100),pNBC(2,100)
dimension Dmat(3,3)
dimension x(nx,numd),noCell1(ng,ncn),ds(nx,numd)
dimension xc(nx,numdq),noCell(ng,numc)
dimension gauss(nx,nqc),gs(ng,numg)
dimension gpos(nx),nv(numd),ph(10,numd)
dimension ak(2*numd,2*numd),GSPk(4*numd*numd)
dimension ne(2*numd),force(2*numd)
dimension u2(nx,numd),disp(2*numd)
dimension Stressnode(3,numd)
common/para/xlength,ylength,p,young,anu,aimo
common/rpim/ALFC,DC,Q,nRBF
common /basis/mbasis
