function H=d_spherbessH(nu,K,z)
H=nu*spherbessH(nu,K,z)-z*spherbessH(nu+1,K,z);