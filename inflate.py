import MDAnalysis as mda
import sys

# Arguments
strfile=sys.argv[1]
ires=int(sys.argv[2])-1
scale=float(sys.argv[3])
outfile=sys.argv[4]

# Open input file
u = mda.Universe(strfile)
# Select residue
res = u.residues[ires]

# Scale residue only (not backbone)
for atom in res.atoms:
    if 'P' in atom.name or "'" in atom.name:
        #print('Discard: ', atom.name)
        continue
    else:
        #print('Use:', atom.name)
        if atom.name == 'N':
            r0 = atom.position
        atom.position -= r0
        atom.position *= scale
        atom.position += r0

# Write output
f = mda.Writer(outfile, n_atoms=u.atoms.n_atoms)
f.write(u.atoms)
