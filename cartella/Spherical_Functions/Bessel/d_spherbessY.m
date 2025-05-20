function Y=d_spherbessY(nu,z)
Y=nu*spherbessY(nu,z)-z*spherbessY(nu+1,z);
