# $Id: countop.py,v 1.1.2.1 2005/03/27 12:57:26 mstorti Exp $
ndof=5
nel=8
ndim=3

vals = {'mu':ndof,'alpha':ndof,'rho':ndof,'m':ndim,
        'sigma':ndof,'j':nel,'m':ndim,'l':ndim,'nu':ndof,
        'k':nel}

matrices = [['taua','mu','alpha'],
            ['cpi','alpha','rho'],
            ['A1','m','rho','sigma'],
            ['A2','l','sigma','nu'],
            ['dNdx1','j','m'],
            ['dNdx2','k','l']]

total = 0

indx = 'alpha'
c=0
for m in matrices:
   if indx in m:
       c=c+1
       print m
if c==2:
    new_matrices = []
    c=0
    for m in matrices:
        if indx in m:
            c=c+1
            if c==1:
                m1=m
            elif c==2:
                m2=m
        else:
            new_matrices.append(m)
    new_mat = '('+m1[0]+'*'+m2[0]+')'
    mm=m1[1:]+m2[1:]
    print m1,m2,mm
    

#while matrices.length()>1:
#    # Find lower operation product
    
