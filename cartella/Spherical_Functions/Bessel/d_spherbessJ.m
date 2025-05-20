function J=d_spherbessJ(nu,z)
J=nu*spherbessJ(nu,z)-z.*spherbessJ(nu+1,z);