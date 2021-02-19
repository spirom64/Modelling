import re
from Scientific.Geometry import Vector
from math import pi, sin, cos, acos, atan2
import argparse as ag

def get_args():
    """Parses command line arguments"""
    parser = ag.ArgumentParser(description='Places cations between oxygens '
                                'according to doi://10.1261/rna.2390311')
    parser.add_argument('-f', '--file', type=ag.FileType('r'),
                        help='PDB file', required=True)

    parser.add_argument('-o', '--output', type=ag.FileType('w+'),
                        help='Output PDB file', required=True)

    parser.add_argument('-c', '--cation-type', type=str,
                        help='Cation type: NA, MG or CA (default NA)',
                        choices = ['NA', 'MG', 'CA'],
                        default='NA')
   
    parser.add_argument('-a', '--acid-type', type=str,
                        help='Amino Acid type: RNA or DNA (default RNA)',
                        choices = ['RNA', 'DNA'],
                        default='RNA')

    parser.add_argument('-r', '--rotation-angle', type=int,
                        help='Rotation angle in degrees (default 5)',
                        default = 5)
   
    parser.add_argument('-n', '--number', type=int,
                        help='Number of added cations, 0 - infinite (default 0)',
                        default = 0)
   
    parser.add_argument('--na-run', action='store_true', 
                        help='Fill with Na cations')
    
    parser.add_argument('--na-limit', type=int,
                        help='Number of Na cations to be added during Na fill'
                        'run. 0 - infinite (default 0)', default = 0)
    
    get_args = parser.parse_args()

    args_dict = {
                'file' : get_args.file,
                'output' : get_args.output,
                'cation' : get_args.cation_type,
                'acid' : get_args.acid_type,
                'angle' : get_args.rotation_angle,
                'number' : get_args.number,
                'na_run' : get_args.na_run,
                'na_limit': get_args.na_limit}

    return args_dict

def na_run(cat_run, na_limit):
    cat_run.cation_type = 'NA'
    if na_limit == 0.:
        cat_run.cation_limit = 0
    else:
        cat_run.cation_limit = cat_run.cation_count + na_limit
    cat_run.test_params()
    cat_run.run()
    return cat_run

class PlaceCation(object):
    def __init__(self, args):
        self.Atoms, self.cations, self.pairs = {}, {}, []
        self.cation_type = args['cation']
        self.cation_limit = args['number']
        self.acid_type = args['acid']
        self.rot_res = args['angle']
        self.output = args['output']
        self.cation_count = 0
       
        coord_file = args['file']
        pdb_i = self.read_pdb(coord_file)
        coord_file.close()
        self.Atoms = pdb_i[0]
        self.nmb = pdb_i[1]
       
        self.R1, self.R2 = '',''
        self.test_params()
       
        self.q_ox = re.compile("(OP1|O1P|OP2|O2P)")
        self.ox_atoms = self.get_ox(self.Atoms, self.q_ox)
        
        self.q_charge = re.compile("(M[Gg]|N[ZH]|N[ED]1)")
        self.pdb_format = ("%-6s%5d %-4s%1s%3s %1s%4d%1s   "
                           "%8.3f%8.3f%8.3f%6.2f%6.2f           %2s%2s%1s")

    def read_pdb(self, pdb):
        """Reads .pdb file and gets atom coordinates"""
        pdb_a = {}
        for line in pdb:
            at = re.compile("(ATOM|HETATM)")
            if at.match(line):
                nm = re.sub(r'\s', '', line[6:12])
                aname = re.sub(r'\s', '', line[12:17])
                ri_c = re.sub(r'\s', '', line[20:27])
                x = re.sub(r'\s', '', line[30:38])
                y = re.sub(r'\s', '', line[38:46])
                z = re.sub(r'\s', '', line[46:55])
                if ri_c and aname and x and y and z:
                    pdb_a[int(nm)] = [aname, Vector(float(x), float(y), float(z)), ri_c]
        return [pdb_a, nm]
   
    def test_params(self):
        """Chooses radii depending on cation and acid types"""
        param = ''.join([self.cation_type, self.acid_type])
        params = {'NARNA':[2.3524, 2.5374],
                  'NADNA':[2.3749, 3.4039],
                  'MGRNA':[2.0040, 2.0046],
                  'MGDNA':[2.0020, 2.0500],
                  'CARNA':[2.2734, 2.3373],
                  'CADNA':[2.2714, 2.3321]}
   
        self.R1, self.R2 = params[param][0], params[param][1]
        return True

    def get_ox(self, all_atoms, q):
        """Extracts oxygens"""
        oxygens = {}
        for at in all_atoms:
            if q.match(all_atoms[at][0]):
                oxygens[at] = all_atoms[at]
        return oxygens
   
    def run(self):
        """Main method""" 
        if self.cation_limit == 0.:
            for ox in self.ox_atoms:
                self.process_atoms(ox)
        else:
            for ox in self.ox_atoms:
                if self.cation_count < self.cation_limit:
                    self.process_atoms(ox)
    
        return True 
    
    def process_atoms(self, ox):
        """Processes atom list"""
        NB = self.find_nb([ox, self.ox_atoms[ox]], self.Atoms,
                self.R1, self.R2)
        NB = self.find_nb_cat(self.ox_atoms[ox], self.cations, NB,
                self.R1, self.R2)
        nb = NB[0]
        check = NB[1]
        
        if nb and not re.search(self.q_charge, check):
            for nb_ox in nb:
                if (self.q_ox.match(nb[nb_ox][0])
                    and self.check_pairs(self.pairs, [ox, nb_ox])):
                    cat_coord = self.place_cat(self.Atoms[ox], nb[nb_ox], nb, self.R1, self.R2, self.rot_res)
                    if not cat_coord:
                        cat_coord = self.place_cat(self.Atoms[ox], nb[nb_ox], nb, self.R2, self.R1, self.rot_res)
                    if cat_coord:
                        self.cation_count = self.cation_count + 1
                        self.cations[self.cation_count] = [self.cation_type, cat_coord]
                        self.pairs.extend([ox, nb_ox])
                        break

        return True

    def find_nb(self, ox1, atoms, r1, r2):
        """finds neighbours in r1 + r2 distance"""
        nb_check = [{}, ""]
        for k in atoms:
            dox = Vector.length(ox1[1][1] - atoms[k][1])
            if (k != ox1[0] and ox1[1][2] != atoms[k][2]  and
                dox <= (r1 + r2)):
                nb_check[0][k] = atoms[k]
                if dox <= r2:
                    nb_check[1] = ''.join([nb_check[1], atoms[k][0]])
        return  nb_check
   
    def find_nb_cat(self, ox1, cations, nb_check, r1, r2):
        """finds cation's neighbours"""
        k = len(nb_check[0])
        if cations:
            for it in cations:
                dcat = Vector.length(ox1[1] - cations[it][1])
                if (dcat <= (r1 + r2)):
                    nb_check[0][k + it] = cations[it]
                    """affects total charge only if in self.R2 distance"""
                    if dcat <= r2:
                        nb_check[1] = ''.join([nb_check[1], cations[it][0]])
        return  nb_check
    
    def check_pairs(self, all_pr, curr):
        """Checks if oxygens already have cation"""
        flag = True
        for pair_ox in all_pr:
            if (curr[0] == pair_ox or curr[1] == pair_ox):
                flag = False
        return flag
   
    def place_cat(self, ox1, ox2, nb_atoms, r1, r2, rs):
        """Places cation between two oxygens"""
        """rotation angle"""
        rot = pi * rs / 180
        """first oxygen"""
        v1 = ox1[1]
        """second oxygen"""
        v2 = ox2[1]
        d = v2 - v1
        k = (r1 - r2 + Vector.length(d)) / (Vector.length(d) * 2)
        """center of ox1-ox2"""
        hd = self.cartesian_to_spherical(d * k)
        corr = acos( Vector.length(d * k) / r1)
        """cross point of two spheres"""
        cp = self.spherical_to_cartesian(r1, hd[1], (hd[2] - corr))
        ev = Vector.normal(d)
        ka = Vector.cross(ev, cp)
        kb = ev*(ev*cp)
        for rt in range(360 / rs):
            ang = rt*rot
            flag = True
            """Rodrigues' formula for vector rotation"""
            cat_vector = cp * cos(ang) + ka*sin(ang) + kb*(1-cos(ang)) + v1
            for test_nb in nb_atoms:
                """distance between cation and other neighbours has to be more
                than 2A"""
                if Vector.length(cat_vector - nb_atoms[test_nb][1]) <= 2.:
                    flag = False
                    break
            if flag:
                return cat_vector
   
    def cartesian_to_spherical(self, v):
        """Converts cartesian coordinates to spherical"""
        x = Vector.x(v)
        y = Vector.y(v)
        z = Vector.z(v)
        r = Vector.length(v)
        phi = atan2(y, x)
        theta = acos(z / r)
       
        return [r, phi, theta]
   
    def spherical_to_cartesian(self, r, phi, theta):
        """Converts spherical coordinates to cartesian"""
        x = r*cos(phi)*sin(theta)
        y = r*sin(phi)*sin(theta)
        z = r*cos(theta)
       
        return Vector(float(x), float(y), float(z))
   
    def print_cat(self):
        """Prints cations in pdb format"""
        for ion in self.cations:
            data_str = ("HETATM", int(self.nmb) + ion, self.cations[ion][0], "",
                    self.cations[ion][0],"", ion, "",
                        Vector.x(self.cations[ion][1]), Vector.y(self.cations[ion][1]),
                        Vector.z(self.cations[ion][1]), 1.00, 0.00,
                        self.cations[ion][0], "","\n")
            out_str = self.pdb_format % data_str
            self.output.write(out_str)
        self.output.close()
        return 1
   
def main():
    """Main function"""
    args = get_args()
    pc = PlaceCation(args)
    pc.run()
    
    if args['na_run'] == True:
        pc = na_run(pc, args['na_limit'])
    
    pc.print_cat()
    return True

if __name__ ==  "__main__":
    main()
