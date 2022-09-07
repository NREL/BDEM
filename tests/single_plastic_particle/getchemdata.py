import numpy as np
import sys

class AMReXParticleHeader(object):
    '''

    This class is designed to parse and store the information 
    contained in an AMReX particle header file. 

    Usage:

        header = AMReXParticleHeader("plt00000/particle0/Header")
        print(header.num_particles)
        print(header.version_string)

    etc...

    '''

    def __init__(self, header_filename):

        self.real_component_names = []
        self.int_component_names = []
        with open(header_filename, "r") as f:
            self.version_string = f.readline().strip()

            particle_real_type = self.version_string.split('_')[-1]
            particle_real_type = self.version_string.split('_')[-1]
            if particle_real_type == 'double':
                self.real_type = np.float64
            elif particle_real_type == 'single':
                self.real_type = np.float32
            else:
                raise RuntimeError("Did not recognize particle real type.")
            self.int_type = np.int32

            self.dim = int(f.readline().strip())
            self.num_int_base = 2
            self.num_real_base = self.dim
            self.num_real_extra = int(f.readline().strip())
            for i in range(self.num_real_extra):
                self.real_component_names.append(f.readline().strip())
            self.num_int_extra = int(f.readline().strip())
            for i in range(self.num_int_extra):
                self.int_component_names.append(f.readline().strip())
            self.num_int = self.num_int_base + self.num_int_extra
            self.num_real = self.num_real_base + self.num_real_extra
            self.is_checkpoint = bool(int(f.readline().strip()))
            self.num_particles = int(f.readline().strip())
            self.max_next_id = int(f.readline().strip())
            self.finest_level = int(f.readline().strip())
            self.num_levels = self.finest_level + 1

            if not self.is_checkpoint:
                self.num_int_base = 0
                self.num_int_extra = 0
                self.num_int = 0

            self.grids_per_level = np.zeros(self.num_levels, dtype='int64')
            for level_num in range(self.num_levels):
                self.grids_per_level[level_num] = int(f.readline().strip())

            self.grids = [[] for _ in range(self.num_levels)]
            for level_num in range(self.num_levels):
                for grid_num in range(self.grids_per_level[level_num]):
                    entry = [int(val) for val in f.readline().strip().split()]
                    self.grids[level_num].append(tuple(entry))

def read_amrex_binary_particle_file(fn, ptype="particles"):

    timestamp=float(fn.split("plt")[1])
    base_fn = fn + "/" + ptype
    header = AMReXParticleHeader(base_fn + "/Header")

    
    idtype = "(%d,)i4" % header.num_int    
    if header.real_type == np.float64:
        fdtype = "(%d,)f8" % header.num_real
    elif header.real_type == np.float32:
        fdtype = "(%d,)f4" % header.num_real
    
    idata = np.empty((header.num_particles, header.num_int ))
    rdata = np.empty((header.num_particles, header.num_real))
    
    ip = 0
    for lvl, level_grids in enumerate(header.grids):
        for (which, count, where) in level_grids:
            if count == 0: continue
            fn = base_fn + "/Level_%d/DATA_%05d" % (lvl, which)

            with open(fn, 'rb') as f:
                f.seek(where)
                if header.is_checkpoint:
                    ints   = np.fromfile(f, dtype = idtype, count=count)
                    idata[ip:ip+count] = ints

                floats = np.fromfile(f, dtype = fdtype, count=count)
                rdata[ip:ip+count] = floats            
            ip += count

    return idata, rdata, timestamp


def getrate(C,Temp):

    nspec=len(C);
    A=2.3e18
    Ea=285.7e3 #J/mol
    Ru=8.314 #J/K/mol
    rate=A*np.exp(-Ea/Ru/Temp);

    dcdt=np.zeros(nspec);

    dcdt[0]=-rate*C[0]
    dcdt[1]=0.6*rate*C[0]
    dcdt[2]=0.4*rate*C[0]

    return(dcdt)

if __name__ == "__main__":
    import glob

    fn_pattern = sys.argv[1]

    fn_list = glob.glob(fn_pattern)
    fn_list.sort()

    outfile=open("avg_quantities.dat","w")

    for i, fn in enumerate(fn_list):
        print("reading %s"%(fn))
        idata, rdata, timestamp = read_amrex_binary_particle_file(fn)
        ppos = rdata[:,0:3]   # assumes 3D
        rad  = rdata[:,3]
        pl = rdata[:,17]
        tr = rdata[:,18]
        gs = rdata[:,19]
        outfile.write("%e\t%e\t%e\t%e\t%e\n"%(timestamp,rad,pl,tr,gs))

    outfile.close()
    
    nspec=3
    rho=1000.0;
    Y=np.zeros(nspec);

    Y[0]=1.0;

    C=rho*Y;
    eps=1.0

    tfinal=100.0;
    tout=0.1;
    dt=0.01

    tnow=0.0;
    toutput=0.0;
    Temp=773.0;

    outfile=open("chemsoln.dat","w")
    outfile.write("%e\t%e\t%e\t%e\t%e\n"%(tnow,Y[0],Y[1],Y[2],eps));

    while(tnow<tfinal):

        C=rho*Y
        dcdt=getrate(C,Temp)

        for i in range(nspec):
            Y[i]+=dcdt[i]*dt/rho

        eps+=dcdt[0]/rho*dt

        tnow+=dt
        toutput+=dt

        if(toutput>tout):
            outfile.write("%e\t%e\t%e\t%e\t%e\n"%(tnow,Y[0],Y[1],Y[2],eps));
            toutput=0.0

    outfile.close()
