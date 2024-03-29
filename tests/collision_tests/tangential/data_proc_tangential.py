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

def get_output_time(fname):

    output_time=1.0
    infile=open(fname,'r')
    lines=infile.readlines()
    for line in lines:
        spltline=line.split("=")
        if(len(spltline)>0):
            if(spltline[0].strip()=="bdem.write_output_time"):
                output_time=float(spltline[1])
                break

    infile.close()
    return(output_time)


if __name__ == "__main__":
    import glob

    fn_pattern = sys.argv[1]
    inpfile=sys.argv[2]
    
    out_time=get_output_time(inpfile)

    fn_list = glob.glob(fn_pattern)
    fn_list.sort()

    outfile=open("particle_data.dat","w")

    for i, fn in enumerate(fn_list):
        print("reading %s"%(fn))
        idata, rdata, timestamp = read_amrex_binary_particle_file(fn)
        overlap_n = rdata[:,17] 
        kn = rdata[:,18]
        kt = rdata[:,19]
        eta_n = rdata[:,20]
        eta_t = rdata[:,21]
        last = rdata[:,22]
        outfile.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(timestamp*out_time,overlap_n[0],
            kn[0],kt[0],eta_n[0],eta_t[0],last[0]))




