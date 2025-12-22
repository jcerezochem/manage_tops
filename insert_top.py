#!/usr/bin/env python3
'''
App to insert a topoly (e.g derived wit Jopyce) into 
a large molecule itp (e.g. a protein or DNA fragment)
Requires the following files:
    -pt: a template itp with the [ atoms ] of the new topology (includes charges). May
         inlcude a region with backbone to be left from topology and new residue
    -pr: the bonded parameters derived elsewere. The fragment can be larger than the 
         residue to assign paramters, but the atom naming must be unique
    -p : the full itp where the new fragment is added
The new topology is written to stdout. 
The resieud to be updated is designed with -nres
'''

import sys
import re
import numpy as np
import datetime

# Version tag (here comes the version tag)
#
# (here comes the COMMIT info)
# (here comes the DATE info)
try:
    version_tag
except:
    class version_tag:
        COMMIT="Untracked"
        DATE="No date"

# Define classes for atom, bond...
class Atom:
    def __init__(self,*args):
#    def __init__(self,iat,attype,ires,resname,atname,chgr,q,mass):
        if len(args) == 8:
            self.iat     = int(args[0])
            self.attype  = args[1]
            self.ires    = int(args[2])
            self.resname = args[3]
            self.atname  = args[4]
            self.chrg    = int(args[5])
            self.q       = float(args[6])
            self.mass    = float(args[7])
            self.V       = 0.0  # sigma(cr=2,3) or C6(cr=1)
            self.W       = 0.0  # epsilon(cr=2,3) or C12(cr=1)
        elif len(args) == 7:
            # E.g. IB+ ion
            self.iat     = int(args[0])
            self.attype  = args[1]
            self.ires    = int(args[2])
            self.resname = args[3]
            self.atname  = args[4]
            self.chrg    = int(args[5])
            self.q       = float(args[6])
            self.mass    = 0.0  # To be set..
            self.V       = 0.0  # sigma(cr=2,3) or C6(cr=1)
            self.W       = 0.0  # epsilon(cr=2,3) or C12(cr=1)
        else:
            raise TypeError('Wrong number of elements to instaciate atom. Expected 7 or 8, got '+str(len(args)))
    def setLJ(self,V,W):
        self.V       = float(V)
        self.W       = float(W)
    def entryline(self):
        return "%5i %5s %5i %5s %5s %5i  %7.4f  %7.4f"%(self.iat,self.attype,self.ires,self.resname,self.atname,self.chrg,self.q,self.mass)

        
class Bond:
    def __init__(self,*args):
        if len(args) == 0:
            self.i1   = 0
            self.i2   = 0
            self.ft   = 0
            self.r0   = 0.0
            self.kb   = 0.0
            self.prms = ''        
        else:
            self.i1   = int(args[0])
            self.i2   = int(args[1])
            self.ft   = int(args[2])
            self.r0   = None
            self.kb   = None
            self.prms = ' '.join(args[3:])
    def setbond(self,i1,i2,ft,params):
        self.i1   = int(i1)
        self.i2   = int(i2)
        self.ft   = int(ft)
        self.prms = params
    def entryline(self):
        return "%5i %5i %2i %s"%(self.i1,self.i2,self.ft,self.prms)


class Pair:
    def __init__(self,*args):
        if len(args) == 0:
            self.i1   = 0
            self.i2   = 0
            self.ft   = 0
            self.prms = ''        
        else:
            self.i1   = int(args[0])
            self.i2   = int(args[1])
            self.ft   = int(args[2])
            self.prms = ' '.join(args[3:])
    def setpair(self,i1,i2,ft,params):
        self.i1   = int(i1)
        self.i2   = int(i2)
        self.ft   = int(ft)
        self.prms = params
    def entryline(self):
        return "%5i %5i %2i %s"%(self.i1,self.i2,self.ft,self.prms)
        
class Angle:
    def __init__(self,*args):
        if len(args) == 0:
            self.i1   = 0
            self.i2   = 0
            self.i3   = 0
            self.ft   = 0
            self.prms = ''        
        else:
            self.i1   = int(args[0])
            self.i2   = int(args[1])
            self.i3   = int(args[2])
            self.ft   = int(args[3])
            self.prms = ' '.join(args[4:])
    def setangle(self,i1,i2,i3,ft,params):
        self.i1   = int(i1)
        self.i2   = int(i2)
        self.i3   = int(i3)
        self.ft   = int(ft)
        self.prms = params
    def entryline(self):
        return "%5i %5i %5i %2i %s"%(self.i1,self.i2,self.i3,self.ft,self.prms)
        
class Dihed:
    def __init__(self,*args):
        if len(args) == 0:
            self.i1   = 0
            self.i2   = 0
            self.i3   = 0
            self.i4   = 0
            self.ft   = 0
            self.prms = ''        
        else:
            self.i1   = int(args[0])
            self.i2   = int(args[1])
            self.i3   = int(args[2])
            self.i4   = int(args[3])
            self.ft   = int(args[4])
            self.prms = ' '.join(args[5:])
    def setdihed(self,i1,i2,i3,i4,ft,params):
        self.i1   = int(i1)
        self.i2   = int(i2)
        self.i3   = int(i3)
        self.i4   = int(i4)
        self.ft   = int(ft)
        self.prms = params
    def entryline(self):
        return "%5i %5i %5i %5i %2i %s"%(self.i1,self.i2,self.i3,self.i4,self.ft,self.prms)

# INPUT PARSER
def get_args():
    
    # Options and their defaults 
    final_arguments = dict()
    final_arguments["-pt"]="template.itp"
    final_arguments["-pr"]="params.itp"
    final_arguments["-p"]="topol.top"
    final_arguments["-s"]="none"
    final_arguments["-res"]=1
    final_arguments["-h"]=False
    # Description of the options
    arg_description = dict()
    arg_description["-pt"] ="Name of the template topology"
    arg_description["-pr"] ="Name of the topology with parameters"
    arg_description["-p"] ="Name of the topologyto be updated"
    arg_description["-s"] ="Swapp file (if not equal none)"
    arg_description["-res"] ="Number of residue to update"
    arg_description["-h"] ="Show this help"
    # Type for arguments
    arg_type = dict()
    arg_type["-pt"] ="char"
    arg_type["-pr"] ="char"
    arg_type["-p"] ="char"
    arg_type["-s"] ="char"
    arg_type["-res"] ="int"
    arg_type["-h"]    ="-"
    
    # Get list of input args
    input_args_list = []
    iarg = -1
    for s in sys.argv[1:]:
        # get -flag [val] arguments 
        if s[0]=="-":
            iarg=iarg+1
            input_args_list.append([s])
        else:
            input_args_list[iarg].append(s)
            
    # Transform into dict. Associtaing lonely flats to boolean   
    input_args_dict=dict()
    for input_arg in input_args_list:
        if len(input_arg) == 1:
            # Boolean option. Can be -Bool or -noBool
            input_arg.append(True)
            if input_arg[0][1:3] == "no":
                input_arg[0] = "-" + input_arg[0][3:]
                input_arg[1] = not input_arg[1]
        elif len(input_arg) != 2:
            raise BaseException("Sintax error. Too many arguments")

        input_args_dict[input_arg[0]] = input_arg[1]
    
    for key,value in input_args_dict.items():
        # Check it is allowed
        isValid = final_arguments.get(key,None)
        if isValid is None:
            raise BaseException("Sintax error. Unknown label: " + key)
        # If valid, update final argument
        final_arguments[key]=value
        
    if final_arguments.get("-h"):
        
        print(f"""
 ----------------------------------------
        {sys.argv[0]}

   Convert topology based on [ Xtypes ]
   databases to explicitly defined 
   potential terms in the topology
   
   Version(GIT HASH): %s
   Date             : %s
 ----------------------------------------
        """%(version_tag.COMMIT,version_tag.DATE))
        print("    Options:")
        print("    --------")
        print('      {0:<10}  {1:^4}  {2:<41}  {3:<7}'.format("Flag","Type","Description","Value"))
        print('      {0:-<10}  {1:-^4}  {2:-<41}  {3:-<7}'.format("","","",""))
        for key,value in final_arguments.items():
            descr = arg_description[key]
            atype = arg_type[key]
            #atype=str(type(value)).replace("<type '","").replace("'>","")
            print('      {0:<10}  {1:^4}  {2:<41}  {3:<7}'.format(key, atype, descr, str(value)))
        print("")
        
        sys.exit()
        
    return final_arguments


def read_top(topfile, swapfile='none'):

    # Read input file
    with open(topfile) as f:
        topentry=[]
        for line in f:
            line=line.replace('\n','')
            line=line.replace('\t','  ')
            line=line.split(';')[0]
            if len(line.lstrip()) != 0:
                topentry.append(line)
 
    # Read swap info (TBD)
    if swapfile != 'none':
        with open(swapfile) as f:
            swap=dict()
            for line in f:
                line=line.replace('\n','')
                line=line.replace('\t','  ')
                line=line.split(';')[0]
                if len(line.lstrip()) != 0:
                    line=line.split('--')[0]
                    i,j = line.split()
                    swap[int(i)] = int(j)
    
    # Prepare arrays for molecule items
    molecules=[]
    atoms=[]
    bonds=[]
    pairs=[]
    angles=[]
    diheds=[]
    # Prepare dictionaries for data types database
    atom_prms=dict()
    bond_prms=dict()
    pair_prms=dict()
    angle_prms=dict()
    dihed_prms=dict()
    
    
    # LOOP0: get molecules
    # Initialize
    molecules=[]
    section=''
    for line in topentry:
        # Use re to locate sections
        # We use a named group to individuate the section. This is done with
        # "?P<sect>" (sections are the blocks in parenthesis)
        pattern = r'(^( )*\[( )*(?P<section>[a-zA-Z_]+)( )*\]( )*$)'
        match = re.match(pattern,line)
        if match:
            section=match.group('section')
        
        elif section == 'molecules':
            line_strip = line.split(';')[0]
            if (line_strip.replace(' ','')) == 0:
                continue
            else:
                data=line_strip.split()
                molecules.append(data[0])

    
    # LOOP1: get atoms
    # Initialize
    section=''
    for line in topentry:
        # Use re to locate sections
        # We use a named group to individuate the section. This is done with
        # "?P<sect>" (sections are the blocks in parenthesis)
        pattern = r'(^( )*\[( )*(?P<section>[a-zA-Z_]+)( )*\]( )*$)'
        match = re.match(pattern,line)
        if match:
            section=match.group('section')
            continue
        
        elif section == 'atoms':
            # Get all items for the atom in data array
            data = line.split()
            # Form the atoms array with a new atom (class)
            # *data expands the array into the elements
            atoms.append(Atom(*data))
            # Get LJ parameters from database...

        elif section == 'bonds':
            data = line.split()
            bonds.append(Bond(*data))
        elif section == 'pairs':
            data = line.split()
            pairs.append(Pair(*data))
        elif section == 'angles':
            data = line.split()
            angles.append(Angle(*data))
        elif section == 'dihedrals':
            data = line.split()
            diheds.append(Dihed(*data))

        if section == 'moleculetype':
            molname = line.split()[0]
            molecules.append(molname)

        if len(molecules)>1:
            print('ERROR: only one moleculetype supported')
            sys.exit(1)


    return atoms, [bonds, pairs, angles, diheds]


#############################################################
#
#  MAIN
#
#############################################################
if __name__ == '__main__':

    # Get arguments
    args=get_args()

    targetres = int(args.get('-res'))

    # Template topoly
    topfile = args.get('-pt')
    atoms_template,params = read_top(topfile)

    templatetop_dict = {}
    templatetop_core = []
    templatetop_env  = []
    for atom in atoms_template:
        if atom.iat < 0:
            templatetop_core.append(atom.atname)
        else:
            templatetop_env.append(atom.atname)

    # Parameter topoly
    topfile = args.get('-pr')
    atoms_par, params = read_top(topfile)

    paramstop_dict = {}
    for atom in atoms_par:
        paramstop_dict[atom.iat] = atom
    bonds_par = []
    for bond in params[0]:
        check = [ paramstop_dict[i].atname in templatetop_core for i in [ bond.i1, bond.i2 ] ]
        if any(check):
            bonds_par.append(bond)
    pairs_par = []
    for pair in params[1]:
        check = [ paramstop_dict[i].atname in templatetop_core for i in [ pair.i1, pair.i2 ] ]
        if any(check):
            pairs_par.append(pair)
    angles_par = []
    for angle in params[2]:
        check = [ paramstop_dict[i].atname in templatetop_core for i in [ angle.i1, angle.i2, angle.i3 ] ]
        if any(check):
            angles_par.append(angle)
    diheds_par = []
    for dihed in params[3]:
        check = [ paramstop_dict[i].atname in templatetop_core for i in [ dihed.i1, dihed.i2, dihed.i3, dihed.i4 ] ]
        if any(check):
            diheds_par.append(dihed)

    # Full topology to update
    topfile = args.get('-p')
    atoms, params = read_top(topfile)

    originaltop_dict = {}
    originaltop_core = []
    originaltop_env  = []
    for atom in atoms:
        if atom.ires == targetres:
            if atom.atname in templatetop_env:
                originaltop_env.append(atom.atname)
            else:
                originaltop_core.append(atom.atname)

    if len(templatetop_env) != len(originaltop_env):
        print('ERROR: number of environment atoms different in original and template')
        sys.exit(2)
    check = [ a==b for a,b in zip(templatetop_env,originaltop_env) ]
    if not all(check):
        print('ERROR: environment atoms different in original and template')
        sys.exit(2)

    nres_old = len(originaltop_core) + len(originaltop_env)
    nres_new = len(templatetop_core) + len(templatetop_env)

    atoms_upd = []
    atom_mapping = {}
    for iat, atom in enumerate(atoms,1):
        if atom.ires < targetres:
            atom_mapping[atom.iat] = iat
            atom.iat = iat
            atoms_upd.append(atom)
        else:
            break 
    # Manage mapping of target res
    nat = len(atoms_upd)
    iat_upd = []
    iat_old = []
    iat_rm  = []
    for iat, atom in enumerate(atoms_template,nat+1):
        if atom.atname in originaltop_env:
            iat_upd.append(iat)
        
    for iat, atom in enumerate(atoms[nat:nat+nres_old+1],nat+1):
        if atom.atname in originaltop_env:
            iat_old.append(iat)
        else:
            iat_rm.append(iat)
    for iat, iat0 in zip(iat_upd,iat_old):
        atom_mapping[iat0] = iat
    # Paste new res (updating ires)
    iatpar_to_iatupd = {}
    for iat, atom in enumerate(atoms_template,nat+1):
        atom.iat = iat
        atom.ires = targetres
        atoms_upd.append(atom)
        # search over atoms_par
        for atompar in atoms_par:
            if atompar.atname == atom.atname:
                iatpar_to_iatupd[atompar.iat] = iat
                break
    # Continue with remaingin residues, updating iat
    for iat, atom in enumerate(atoms[nat+nres_old:],nat+nres_new+1):
        atom_mapping[atom.iat] = iat
        atom.iat = iat
        atoms_upd.append(atom)

    # Print
    print('[ atoms ]')
    ires = 0
    for atom in atoms_upd:
        if ires != atom.ires:
            ires = atom.ires
            print(f'; residue {ires}')
        print(atom.entryline())
    print('')

    print('[ bonds ]')
    for bond in params[0]:
        check = [ iat in iat_rm for iat in [ bond.i1, bond.i2 ] ]
        if any(check):
            print(';', bond.entryline())
        else:
            bond.i1 = atom_mapping[bond.i1]
            bond.i2 = atom_mapping[bond.i2]
            print(bond.entryline())
    print('; new bonds with updated residue')
    for bond in bonds_par:
        bond.i1 = iatpar_to_iatupd[bond.i1]
        bond.i2 = iatpar_to_iatupd[bond.i2]
        print(bond.entryline())
    print('')

    print('[ pairs ]')
    for pair in params[1]:
        check = [ iat in iat_rm for iat in [ pair.i1, pair.i2 ] ]
        if any(check):
            print(';', pair.entryline())
        else:
            pair.i1 = atom_mapping[pair.i1]
            pair.i2 = atom_mapping[pair.i2]
            print(pair.entryline())
    print('; new pairs with updated residue')
    for pair in pairs_par:
        pair.i1 = iatpar_to_iatupd[pair.i1]
        pair.i2 = iatpar_to_iatupd[pair.i2]
        print(pair.entryline())
    print('')

    print('[ angles ]')
    for angle in params[2]:
        check = [ iat in iat_rm for iat in [ angle.i1, angle.i2, angle.i3 ] ]
        if any(check):
            print(';', angle.entryline())
        else:
            angle.i1 = atom_mapping[angle.i1]
            angle.i2 = atom_mapping[angle.i2]
            angle.i3 = atom_mapping[angle.i3]
            print(angle.entryline())
    print('; new angles with updated residue')
    for angle in angles_par:
        angle.i1 = iatpar_to_iatupd[angle.i1]
        angle.i2 = iatpar_to_iatupd[angle.i2]
        angle.i3 = iatpar_to_iatupd[angle.i3]
        print(angle.entryline())
    print('')

    print('[ dihedrals ]')
    for dihed in params[3]:
        check = [ iat in iat_rm for iat in [ dihed.i1, dihed.i2, dihed.i3, dihed.i4 ] ]
        if any(check):
            print(';', dihed.entryline())
        else:
            dihed.i1 = atom_mapping[dihed.i1]
            dihed.i2 = atom_mapping[dihed.i2]
            dihed.i3 = atom_mapping[dihed.i3]
            dihed.i4 = atom_mapping[dihed.i4]
            print(dihed.entryline())
    print('; new diheds with updated residue')
    for dihed in diheds_par:
        dihed.i1 = iatpar_to_iatupd[dihed.i1]
        dihed.i2 = iatpar_to_iatupd[dihed.i2]
        dihed.i3 = iatpar_to_iatupd[dihed.i3]
        dihed.i4 = iatpar_to_iatupd[dihed.i4]
        print(dihed.entryline())
    print('')


