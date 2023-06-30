#!/usr/bin/env python
<<<<<<< HEAD:src/pmx/scripts/mutate.py
=======

>>>>>>> develop:pmx/extensions/pmx/scripts/mutate.py
# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

<<<<<<< HEAD:src/pmx/scripts/mutate.py
__doc__="""
Program to insert mutated residues in structure files for
free energy simulations (so far unfinished new version).
"""

import sys,os
from pmx import *
from pmx.parser import *
from pmx import library
from pmx.mutdb import *
from pmx.geometry import *

class UnknownResidueError(Exception):
    def __init__(self, s):
        self.s = s
    def __str__(self):
        return repr(self.s)
class RangeCheckError(Exception):
    def __init__(self, s):
        self.s = s
    def __str__(self):
        return repr(self.s)
class mtpError(Exception):
    def __init__(self,s):
        self.s = s
    def __str__(self):
        return repr(self.s)



ext_one_letter = {
    'ALA':'A',
    'ARG':'R',
    'ASN':'N',
    'ASP':'D',
    'ASPH':'B',
    'ASPP':'B',
    'ASH':'B',
    'CYS':'C',
    'CYS2':'C',
    'CYN':'C',
    'CYX':'CX',
    'CYM':'CM',
    'CYSH':'C',
    'GLU':'E',
    'GLUH':'J',
    'GLUP':'J',
    'GLH':'J',
    'GLN':'Q',
    'GLY':'G',
    'HIS':'H',
    'HIE':'X',
    'HISE':'X',
    'HSE':'X',
    'HIP':'Z',
    'HSP':'Z',
    'HISH':'Z',
    'HID':'H',
    'HSD':'H',
    'ILE':'I',
    'LEU':'L',
    'LYS':'K',
    'LYSH':'K',
    'LYP':'K',
    'LYN':'O',
    'LSN':'O',
    'MET':'M',
    'PHE':'F',
    'PRO':'P',
    'SER':'S',
    'SP1':'SP1', # phosphoserine in charmm36
    'SP2':'SP2', # phosphoserine in charmm36
    'THR':'T',
    'TRP':'W',
    'TYR':'Y',
    'VAL':'V',
}

noncanonical_aa = {
    'S2SP1':'SSP1', # serine to pSer1
    'S2SP2':'SSP2', # serine to pSer2
    'SP12S':'SP1S', # pSer1 to serine
    'SP22S':'SP2S', # pSer2 to setine
}

dna_names = {
    'DA5_DT5':'D5K',
    'DA5_DC5':'D5L',
    'DA5_DG5':'D5M',
    'DT5_DA5':'D5N',
    'DT5_DC5':'D5O',
    'DT5_DG5':'D5P',
    'DC5_DA5':'D5R',
    'DC5_DT5':'D5S',
    'DC5_DG5':'D5T',
    'DG5_DA5':'D5X',
    'DG5_DT5':'D5Y',
    'DG5_DC5':'D5Z',
    'DA3_DT3':'D3K',
    'DA3_DC3':'D3L',
    'DA3_DG3':'D3M',
    'DT3_DA3':'D3N',
    'DT3_DC3':'D3O',
    'DT3_DG3':'D3P',
    'DC3_DA3':'D3R',
    'DC3_DT3':'D3S',
    'DC3_DG3':'D3T',
    'DG3_DA3':'D3X',
    'DG3_DT3':'D3Y',
    'DG3_DC3':'D3Z',
# False names to avoid an error
    'DG3_DG3':'FOO',
    'DC3_DC3':'FOO',
    'DA3_DA3':'FOO',
    'DT3_DT3':'FOO',
    'DG5_DG5':'FOO',
    'DC5_DC5':'FOO',
    'DA5_DA5':'FOO',
    'DT5_DT5':'FOO',
    }

rna_names = {
    'RA5_RU5':'R5K',
    'RA5_RC5':'R5L',
    'RA5_RG5':'R5M',
    'RU5_RA5':'R5N',
    'RU5_RC5':'R5O',
    'RU5_RG5':'R5P',
    'RC5_RA5':'R5R',
    'RC5_RU5':'R5S',
    'RC5_RG5':'R5T',
    'RG5_RA5':'R5X',
    'RG5_RU5':'R5Y',
    'RG5_RC5':'R5Z',
    'RA3_RU3':'R3K',
    'RA3_RC3':'R3L',
    'RA3_RG3':'R3M',
    'RU3_RA3':'R3N',
    'RU3_RC3':'R3O',
    'RU3_RG3':'R3P',
    'RC3_RA3':'R3R',
    'RC3_RU3':'R3S',
    'RC3_RG3':'R3T',
    'RG3_RA3':'R3X',
    'RG3_RU3':'R3Y',
    'RG3_RC3':'R3Z',
# False names to avoid an error
    'RG3_RG3':'FOO',
    'RC3_RC3':'FOO',
    'RA3_RA3':'FOO',
    'RU3_RU3':'FOO',
    'RG5_RG5':'FOO',
    'RC5_RC5':'FOO',
    'RA5_RA5':'FOO',
    'RU5_RU5':'FOO',
    }

def check_residue_name( res ):
    if res.resname == 'LYS':
        if res.has_atom( 'HZ3'):
=======
"""Program to insert mutated residues in structure files for
free energy simulations.
"""

import sys
import argparse
from pmx import library
from pmx.model import Model
from pmx.parser import read_and_format
from pmx.utils import get_ff_path, ff_selection
from pmx.alchemy import mutate
from pmx.scripts.cli import check_unknown_cmd

# resinfo
dna_one_letter = {'A': 'adenosine',
                  'C': 'cytosine',
                  'G': 'guanine',
                  'T': 'thymine'}

rna_one_letter = {'A': 'adenosine',
                  'C': 'cytosine',
                  'G': 'guanine',
                  'U': 'uracil'}

# ================
# Helper functions
# ================
def _print_sorted_dict(d):
    for key in sorted(d.iterkeys()):
        print("{0:>5}{1:>15}".format(key, d[key]))


def _int_input():
    inp = input()
    try:
        inp = int(inp)
        return inp
    except:
        return( inp )
#        print('You entered "%s" -> Try again' % inp)
#        return None


def _check_residue_name(res):
    if res.resname == 'LYS':
        if res.has_atom('HZ3'):
>>>>>>> develop:pmx/extensions/pmx/scripts/mutate.py
            res.set_resname('LYP')
    elif res.resname == 'HIS':
        if res.has_atom('HD1') and \
           res.has_atom('HE2'):
            res.set_resname('HIP')
<<<<<<< HEAD:src/pmx/scripts/mutate.py
        elif res.has_atom('HD1') and not \
                 res.has_atom('HE2'):
            res.set_resname( 'HID' )
        elif not res.has_atom('HD1') and  \
                 res.has_atom('HE2'):
            res.set_resname( 'HIE' )
=======
        elif res.has_atom('HD1') and not res.has_atom('HE2'):
            res.set_resname('HID')
        elif not res.has_atom('HD1') and res.has_atom('HE2'):
            res.set_resname('HIE')
>>>>>>> develop:pmx/extensions/pmx/scripts/mutate.py
    elif res.resname == 'ASP':
        if res.has_atom('HD2'):
            res.set_resname('ASH')
    elif res.resname == 'GLU':
        if res.has_atom('HE2'):
            res.set_resname('GLH')
    elif res.resname == 'CYS':
        if not res.has_atom('HG'):
<<<<<<< HEAD:src/pmx/scripts/mutate.py
            print >>sys.stderr,' Cannot mutate SS-bonded Cys %d' % res.id

def check_OPLS_LYS( res ):
    if res.has_atom( 'HZ3'):
        return('K')
    else:
	return('O')

#def get_restype(r):
#    if r.resname in ['DA','DT','DC','DG']:
#        return 'DNA'
#    elif r.resname in ['RA','RU','RC','RG']:
#        return 'RNA'
#    else: return 'PEPTIDE'

def read_script(fn):
    return read_and_format(fn,"is")

def int_input():
    inp = raw_input()
    try:
        inp = int(inp)
        return inp
    except:
        print 'You entered "%s" -> Try again' % inp
        return None

def check_residue_range(m, idx):
    valid_ids = range(1, len(m.residues)+1)
    if idx not in valid_ids: return False
    return True

def select_residue(m):
    valid_ids = range(1, len(m.residues)+1)
    print '\nSelect residue to mutate:'
    for i,r in enumerate(m.residues):
        if r.resname not in library._ions+library._water:
            sys.stdout.write('%6d-%s-%s' % (r.id,r.resname,r.chain_id))
            if r.id % 6 == 0: print
    print
    selected_residue_id = None
    while not selected_residue_id:
        sys.stdout.write('Enter residue number: ')
        selected_residue_id = int_input()
        if selected_residue_id is not None and selected_residue_id not in valid_ids:
            print 'Residue id %d not in range %d-%d -> Try again' % (selected_residue_id,1,len(residues))
            selected_residue_id = None
    return selected_residue_id

def select_mutation(m, selected_residue_id, ffpath):

    residue = m.residues[selected_residue_id - 1]
    if get_restype(residue) == 'PEPTIDE':
        return select_aa_mutation(residue,ffpath)
    elif get_restype(residue) in ['DNA','RNA']:
        return select_nuc_mutation(residue)

def select_nuc_mutation(residue):
    aa = None
    print '\nSelect new base for %s-%s: ' % (residue.id,residue.resname)
    sys.stdout.write('One-letter code: ')
    while aa is None:
        aa = raw_input().upper()
        if get_restype(residue) == 'DNA' and aa not in ['A','C','G','T']:
            sys.stdout.write('Unknown DNA residue "%s"!\nOne-letter code: ' % aa)
            aa = None
        elif get_restype(residue) == 'RNA' and aa not in ['A','C','G','U']:
            sys.stdout.write('Unknown RNA residue "%s"!\nOne-letter code: ' % aa)
            aa = None
        if aa:
            print 'Will apply mutation %s->%s on residue %s-%d' % (residue.resname[1],aa,residue.resname,residue.id)
        return aa

def select_aa_mutation(residue,ffpath):
    check_residue_name( residue )
    print '\nSelect new amino acid for %s-%s: ' % (residue.id,residue.resname)
    sys.stdout.write('Three- or one-letter code (or four-letter for ff specific residues): ')
    if residue.resname in ['HIE','HISE','HSE']: rol = 'X'
    elif residue.resname in ['HIP','HISH','HSP']: rol = 'Z'
    elif residue.resname in ['GLH','GLUH','GLUP']: rol = 'J'
    elif residue.resname in ['ASH','ASPH','ASPP']: rol = 'B'
    elif residue.resname in ['LYN','LYS','LSN']: rol = 'O'
    else:
        rol = library._one_letter[residue.resname]
    aa = None
    ol = library._aacids_dic.keys()
    tl = library._aacids_dic.values()
    ffpathlower = ffpath.lower()
    if('amber' in ffpathlower):
            ol = library._aacids_ext_amber.keys()
            tl = library._aacids_ext_amber.values()
    if('opls' in ffpathlower):
            ol = library._aacids_ext_oplsaa.keys()
            tl = library._aacids_ext_oplsaa.values()+['ASPP','GLUP','LSN']
    if('charmm' in ffpathlower):
            ol = library._aacids_ext_charmm.keys()
            tl = library._aacids_ext_charmm.values()

    while aa is None:
        aa = raw_input().upper()
        if len(aa) != 1 and len(aa)!=3 and len(aa)!=4:
            sys.stdout.write('Nope!\nThree- or one-letter code (or four-letter for ff specific residues): ')
            aa = None
        elif (len(aa) == 1 and aa not in ol+['B','J','O','X','Z']) or (len(aa)==3 and aa not in tl) or (len(aa)==4 and aa not in tl):
            sys.stdout.write('Unknown aa "%s"!\nThree- or one-letter code (or four-letter for ff specific residues): ' % aa)
            aa = None
        if aa and (len(aa)==3 or len(aa)==4): aa = ext_one_letter[aa]
    print 'Will apply mutation %s->%s on residue %s-%d' % (rol,aa,residue.resname,residue.id)
    return aa


def interactive_selection(m,ffpath):
    residue_id = select_residue(m)
    mutation = select_mutation(m, residue_id, ffpath )
    return residue_id, mutation

def ask_next():
    sys.stdout.write('\nApply another mutation [y/n]? ')
    res = raw_input().lower()
    if res == 'y': return True
    elif res == 'n': return False
    else: return ask_next()

def convert_aa_name( aa ):
    if len(aa) == 1: return aa.upper()
    elif len(aa)==2 and aa=='CM': return(aa) # special case for two letter code
    elif len(aa) == 3: return ext_one_letter[aa.upper()]
    elif len(aa) == 4: return ext_one_letter[aa.upper()]
    else: raise UnkownResidueError(aa)

def rename_to_match_library(res):
    name_hash = {}
    atoms = res.atoms
    for atom in atoms:
	foo = atom.name
	## for serine
	if (atom.resname == 'SER') and (atom.name == 'HG1'):
	    atom.name = 'HG'
        if ('S2' in atom.resname) and (atom.name == 'HG1'):
            atom.name = 'HG'
        if ('SP1' in atom.resname) and (atom.name == 'HG1'): # phosphoserine in charmm36
            atom.name = 'HG'
        if ('SP2' in atom.resname) and (atom.name == 'HG1'): # phosphoserine in charmm36
            atom.name = 'HG'
	## for cysteine
        if (atom.resname == 'CYS') and (atom.name == 'HG1'):
            atom.name = 'HG'
        if ('C2' in atom.resname) and (atom.name == 'HG1'):
            atom.name = 'HG'
#	print atom.resname,atom.name
	name_hash[atom.name] = foo
    return name_hash

def rename_back( res, name_hash ):
    for atom in res.atoms:
        atom.name = name_hash[atom.name]

def set_conformation(old_res, new_res, rotdic):
    old_res.get_real_resname()
    dihedrals = library._aa_dihedrals[old_res.real_resname]
    for key, lst in rotdic.items():
        new = new_res.fetchm(lst)
        rotdic[key] = new
    chis = []
    for key in rotdic.keys():
        at1,at2 = key.split('-')
        for d in dihedrals:
            if d[1] == at1 and d[2] == at2 \
                   and d[-1] != -1:
                chis.append(d)
    for d in chis:
        atoms = old_res.fetchm(d[:4])
        phi = atoms[0].dihedral(atoms[1], atoms[2], atoms[3])
        atoms2 = new_res.fetchm(d[:4])
        phi2 = atoms2[0].dihedral(atoms2[1], atoms2[2], atoms2[3])
        diff = phi-phi2
        a1,a2 = new_res.fetchm(d[1:3])
        key= a1.name+'-'+a2.name
        atoms = rotdic[key]
        rot = Rotation(a1.x,a2.x)
        for atom in atoms:
            atom.x = rot.apply(atom.x,diff)
#    sys.exit(0)
    for atom in new_res.atoms:
        if atom.name[0] != 'D':
            atom.x = old_res[atom.name].x

def get_nuc_hybrid_resname(residue,new_nuc_name,bRNA=False):
    firstLetter = 'D'
    if bRNA:
	firstLetter = 'R'

    # identify if the nucleotide is terminal
    for a in residue.atoms:
	if a.name=='H3T':
	    r1 = firstLetter+residue.resname[1]+'3'
	    r2 = firstLetter+new_nuc_name+'3'
	    dict_key = r1+'_'+r2
	    if bRNA:
	        hybrid_residue_name = rna_names[dict_key]
	    else:
	        hybrid_residue_name = dna_names[dict_key]
	    return(hybrid_residue_name,residue.resname[1],new_nuc_name)
	elif a.name=='H5T':
	    r1 = firstLetter+residue.resname[1]+'5'
	    r2 = firstLetter+new_nuc_name+'5'
	    dict_key = r1+'_'+r2
	    if bRNA:
	        hybrid_residue_name = rna_names[dict_key]
	    else:
	        hybrid_residue_name = dna_names[dict_key]
	    return(hybrid_residue_name,residue.resname[1],new_nuc_name)
    hybrid_residue_name = residue.resname+new_nuc_name
    return(hybrid_residue_name,residue.resname[1],new_nuc_name)

def apply_nuc_mutation(m, residue, new_nuc_name, mtp_file, bRNA=False):

#    hybrid_residue_name = residue.resname+new_nuc_name
    hybrid_residue_name,resname1,resname2 = get_nuc_hybrid_resname(residue,new_nuc_name,bRNA)
    print 'log_> Residue to mutate: %d | %s | %s ' % ( residue.id, residue.resname, residue.chain_id)
    print 'log_> Mutation to apply: %s->%s' % (residue.resname[1], new_nuc_name)
    print 'log_> Hybrid residue name: %s' % hybrid_residue_name
    hybrid_res, bonds, imps, diheds, rotdic = get_hybrid_residue(hybrid_residue_name, mtp_file)
#    hybrid_res.nm2a()

    nuc_super( residue, hybrid_res, resname1, resname2 )
    for atom in hybrid_res.atoms:
        if atom.name[0] != 'D':
            atom.x = residue[atom.name].x
    m.replace_residue( residue, hybrid_res)
    print 'log_> Inserted hybrid residue %s at position %d (chain %s)' %\
          (hybrid_res.resname, hybrid_res.id, hybrid_res.chain_id)


def apply_aa_mutation(m, residue, new_aa_name, mtp_file, bStrB, infileB):

    if residue.resname == 'ILE': rename_ile( residue )
    olkey = convert_aa_name( residue.resname )

    # olkey should contain the correct one letter name of the WT residue
    # however, due to the different namings of the residues in the FFs
    # Lys needs to be checked once again: in OPLS Lys is non-protonated, while in the other FFs it is protonated
    if ('opls' in mtp_file) and ('LYS' in residue.resname):
        olkey = check_OPLS_LYS( residue )

    hybrid_residue_name = olkey+'2'+new_aa_name
    if hybrid_residue_name in noncanonical_aa.keys():
        hybrid_residue_name = noncanonical_aa[hybrid_residue_name]
#    if hybrid_residue_name in longname_aa.keys():
#        hybrid_residue_name = longname_aa[hybrid_residue_name]
    print 'log_> Residue to mutate: %d | %s | %s ' % ( residue.id, residue.resname, residue.chain_id)
    print 'log_> Mutation to apply: %s->%s' % (olkey, new_aa_name)
    print 'log_> Hybrid residue name: %s' % hybrid_residue_name
    hybrid_res, bonds, imps, diheds, rotdic = get_hybrid_residue(hybrid_residue_name, mtp_file)
    #hybrid_res.nm2a()
    bb_super(residue, hybrid_res )

    ## VG rename residue atoms
    hash1 = rename_to_match_library(residue)
    hash2 = rename_to_match_library(hybrid_res)
    set_conformation(residue, hybrid_res, rotdic)
    if bStrB:
	print "log_> Set Bstate geometry according to the provided structure"
   	mB = Model(infileB,bPDBTER=True)
   	rename_atoms_to_gromacs( mB )
	mB.nm2a()
	residueB = mB.residues[residue.id-1]
    	bb_super(residue, residueB )
	for atom in hybrid_res.atoms:
            if atom.name[0] == 'D':
	        for atomB in residueB.atoms:
		    if atomB.name == hybrid_res.morphes[atom.name]['n1']:
	 	        atom.x = atomB.x
    rename_back(residue,hash1)
    rename_back(hybrid_res,hash2)
    ## VG rename residue atoms back

    m.replace_residue( residue, hybrid_res)
    print 'log_> Inserted hybrid residue %s at position %d (chain %s)' %\
          (hybrid_res.resname, hybrid_res.id, hybrid_res.chain_id)


def apply_mutation(m, mut, mtp_file, bStrB, infileB, bRNA):
    residue_id = mut[0]
    if not check_residue_range(m, residue_id):
        raise RangeCheckError(residue_id)
    residue = m.residues[residue_id - 1]
    if get_restype(residue) == 'PEPTIDE':
        new_aa_name = convert_aa_name( mut[1] )
        apply_aa_mutation(m, residue, new_aa_name, mtp_file, bStrB, infileB)
    elif get_restype(residue) in ['DNA','RNA']:
        new_nuc_name = mut[1].upper()
        apply_nuc_mutation(m, residue, new_nuc_name, mtp_file, bRNA)


def get_hybrid_residue(residue_name, mtp_file = 'ffamber99sb.mtp'):
    print 'log_> Scanning database for %s ' % residue_name
    resi, bonds, imps, diheds, rotdic = read_mtp_entry(residue_name, filename = mtp_file, version = 'new')
    if len(resi.atoms) == 0:
        raise mtpError("Hybrid residue %s not found in %s" % (residue_name, mtp_file) )
    return resi, bonds, imps, diheds, rotdic



def rename_ile(residue):
    dic = {'CD':'CD1',
           'HD1':'HD11',
           'HD2':'HD12',
           'HD3':'HD13'
           }
    for key, value in dic.items():
        try:
            atom = residue[key]
            atom.name = value
        except:
            pass

def rename_atoms_to_gromacs( m ):
    for atom in m.atoms:
        if atom.name[0].isdigit():
            atom.name =  atom.name[1:]+atom.name[0]




def get_restype(r):
    if r.resname in ['DA','DT','DC','DG','DA3','DT3','DC3','DG3','DA5','DT5','DC5','DG5']:
        return 'DNA'
    elif r.resname in ['RA','RU','RC','RG','RA3','RU3','RC3','RG3','RA5','RU5','RC5','RG5']:
        return 'RNA'
    else: return 'PEPTIDE'


def get_ff_path( ff ):
    ff_path = None
    if not os.path.isdir(ff):
        gmxlib = os.environ.get('GMXLIB')
        p = os.path.join(gmxlib,ff)
        pff = p+'.ff'
        if os.path.isdir(p):
            ff_path = p
        elif os.path.isdir(pff):
            ff_path = pff
        else:
            print >>sys.stderr,' Error: forcefield path "%s" not found' % ff
            sys.exit(0)
    else:
        ff_path = ff
    print 'Opening forcefield: %s' % ff_path
    return ff_path


def main(argv):

   options = [
        Option( "-resinfo", "bool", False, "print a 3-letter -> 1-letter residue list"),
        Option( "-dna", "bool", False, "generate hybrid residue for the DNA nucleotides"),
        Option( "-rna", "bool", False, "generate hybrid residue for the RNA nucleotides"),
##         Option( "-r", "rvec", [1,2,3], "some string"),
##         Option( "-b", "bool", True, "bool"),
##         Option( "-r2", "rvec", [1,2,3], "some vector that does wonderful things and returns always segfaults")
        ]

   files = [
       FileOption("-f", "r",["pdb","gro"],"protein.pdb", "input structure file"),
       FileOption("-fB", "r",["pdb","gro"],"proteinB.pdb", "input structure file of the Bstate (optional)"),
       FileOption("-o", "w",["pdb","gro"],"out.pdb", "output structure file"),
       FileOption("-ff", "dir",["ff"],"amber99sbmut", "path to mutation forcefield"),
       FileOption("-script", "r",["txt"],"mutations.txt", "text file with mutations to insert"),
       ]

   help_text = ('This script applies mutations of residues in a structure file ',
                'for subsequent free energy calculations like FEP, TI, etc.',
                'The mutation information and dummy placements are taken from',
                'the hybrid residue database "mutres.mtp". The best way to use',
                'this script is to take a pdb/gro file that has been written with pdb2gmx',
                'with all hydrogen atoms present.'
                'The program can either be executed interactively or via script.',
                'The script file simply has to consist of "resi_number target_residue." pairs.',
                'The script uses an extended one-letter code for amino acids to account for',
                'different protonation states. Use the -resinfo flag to print the dictionary.',
                'Currently available force fields:',
                '    - amber99sbmut (Hornak et al, 2006)',
                '    - amber99sb-star-ildn-mut (Best & Hummer, 2009; Lindorff-Larsen et al, 2010)',
                '    - charmm22starmut.ff (Piana et al, 2011)',
                '    - charmm36mut (Best et al, 2012)',
                '    - oplsaamut (Jorgensen et al, 1996; Kaminski et al, 2001)',
                '',
                '',
                'Please cite:',
		'Vytautas Gapsys, Servaas Michielssens, Daniel Seeliger and Bert L. de Groot.',
		'Automated Protein Structure and Topology Generation for Alchemical Perturbations.',
		'J. Comput. Chem. 2015, 36, 348-354. DOI: 10.1002/jcc.23804',
		'',
		'Old pmx (pymacs) version:',
                'Daniel Seeliger and Bert L. de Groot. Protein Thermostability Calculations Using',
                'Alchemical Free Energy Simulations, Biophysical Journal, 98(10):2309-2316 (2010)',
                '',
                '',
                '',
                )


   cmdl = Commandline( argv, options = options,
                       fileoptions = files,
                       program_desc = help_text,
                       check_for_existing_files = False )

   bDNA = cmdl['-dna']
   bRNA = cmdl['-rna']

   if cmdl['-resinfo']:
       print 'Residue dictionary:'
       lst = ext_one_letter.items()
       lst.sort(lambda a,b: cmp(a,b))
       for key, val in lst:
           print "%5s %4s" % (key, val)
       sys.exit(0)

   bStrB = False
   infileB = ''
   if cmdl.opt['-fB'].is_set:
	bStrB = True
	infileB = cmdl['-fB']

   ffpath = get_ff_path(cmdl['-ff'])
   if bDNA:
       mtp_file = os.path.join( ffpath,'mutres_dna.mtp')
   elif bRNA:
       mtp_file = os.path.join( ffpath,'mutres_rna.mtp')
   else:
       mtp_file = os.path.join( ffpath,'mutres.mtp')
   infile = cmdl['-f']

   m = Model(infile,bPDBTER=True)

   rename_atoms_to_gromacs( m )
#   m.write('ll.pdb')
   m.nm2a()
#   m.rename_atoms()
   mutation_list = []
   if cmdl.opt['-script'].is_set:
       mutations_to_make = read_script( cmdl['-script'] )
       for mut in mutations_to_make:
	   check_residue_name( m.residues[ mut[0]-1 ] )
           apply_mutation( m, mut, mtp_file, bStrB, infileB, bRNA )
   else:
       do_more = True
       while do_more:
           mutation = interactive_selection(m,ffpath)
           apply_mutation( m, mutation, mtp_file, bStrB, infileB, bRNA )
           if not ask_next(): do_more = False


   m.write(cmdl['-o'],bPDBTER=True)
   print
   print 'mutations done...........'
   print


def entry_point():
    main(sys.argv[1:])


if __name__ == '__main__':
    main(sys.argv)
=======
            print(' Cannot mutate SS-bonded Cys %d' % res.id, file=sys.stderr)


def _ask_next():
    sys.stdout.write('\nApply another mutation [y/n]? ')
    res = input().lower()
    if res == 'y':
        return True
    elif res == 'n':
        return False
    else:
        return _ask_next()


def _match_mutation(m, ref_m, ref_chain, ref_resid):
    """Matches chain/indices of two Models. Given the chain and resid of a
    reference Model (ref_m), return the resid of the Model (m).

    Parameters
    ----------
    m: Model
        model you want to mutate
    ref_m : Model
        reference model
    ref_chain: str
        chain of residue in reference model
    ref_resid: int
        resid of residue in reference model

    Returns
    -------
    resid: int
        residue ID of model that corresponds to the chain/resid in the
        reference.
    """

    # select non-solvent residues
    res = [r for r in m.residues if r.moltype not in ['water', 'ion']]
    ref_res = [r for r in ref_m.residues if r.moltype not in ['water', 'ion']]
    # check they have same len
    assert len(res) == len(ref_res)

    # iterate through all residue pairs
    resmap = {}
    for r, rref in zip(res, ref_res):
        # first, check that the sequence is the same
        if r.resname != rref.resname:
            raise ValueError('residue %s in the input file does not match '
                             'residue %s in the input reference file'
                             % (r.resname, rref.resname))
        # then, create a dict to map (chain_id, res_id) in the reference
        # to the residues ID in the input file
        resmap[(rref.chain_id, rref.id)] = r.id

    resid = resmap[(ref_chain, ref_resid)]
    print('log_> Residue {ref_id} (chain {ref_ch}) in file {ref} mapped to residue '
          '{m_id} in file {m} after renumbering'.format(ref_ch=ref_chain,
                                                        ref_id=ref_resid,
                                                        ref=ref_m.filename,
                                                        m_id=resid,
                                                        m=m.filename))

    return resid


# ===============================
# Class for interactive selection
# ===============================
class InteractiveSelection:
    """Class containing fucntions related to the interactive selection of
    residues to be mutated.

    Parameters
    ----------
    m : Model object
        instance of pmx.model.Model
    ffpath : str
        path to forcefield files

    Attributes
    ----------
    mut_resid : int
        index of residue to be mutated
    mut_resname : str
        one-letter code of target residue

    """

    def __init__(self, m, ff, renumbered=True):
        self.m = m
        self.ffpath = get_ff_path(ff)

        # get selection
        if renumbered is True:
            self.mut_chain = None
        elif renumbered is False:
            self.mut_chain = self.select_chain()

        self.mut_resid = self.select_residue()
        self.mut_resname = self.select_mutation()

    def select_chain(self):
        """Ask for the chain id to mutate.
        """
        # show selection
        valid_ids = [c.id for c in self.m.chains]
        print('\nSelect a chain:')
        for c in self.m.chains:
            print('{0:>6}'.format(c.id))

        # select id
        selected_chain_id = None
        while selected_chain_id is None:
            sys.stdout.write('Enter chain ID: ')
            selected_chain_id = input()
            if selected_chain_id is not None and selected_chain_id not in valid_ids:
                print('Chain id %s not among valid IDs -> Try again' % selected_chain_id)
                selected_chain_id = None
        return selected_chain_id

    def select_residue(self):
        """Ask for the residue id to mutate.
        """
        # show selection if we do not need chain ID
        if self.mut_chain is None:
            valid_ids = [r.id for r in self.m.residues]
            print('\nSelect residue to mutate:')
            for i, r in enumerate(self.m.residues):
                if r.moltype not in ['water', 'ion']:
                    sys.stdout.write('%6d-%s-%s' % (r.id, r.resname, r.chain_id))
                    if (i+1) % 6 == 0:
                        print("")
        elif self.mut_chain is not None:
            valid_ids = [r.id for r in self.m.chdic[self.mut_chain].residues]
            print('\nSelect residue to mutate:')
            for i, r in enumerate(self.m.chdic[self.mut_chain].residues):
                if r.moltype not in ['water', 'ion']:
                    sys.stdout.write('%6s-%s-%s' % (str(r.id), r.resname, r.chain_id))
                    if (i+1) % 6 == 0:
                        print("")
        print("")

        # select id
        selected_residue_id = None
        while not selected_residue_id:
            sys.stdout.write('Enter residue number: ')
            selected_residue_id = _int_input()
            if (selected_residue_id is not None) and (selected_residue_id not in valid_ids):
                print('Residue id %s not among valid IDs -> Try again' % str(selected_residue_id))
                selected_residue_id = None
        return selected_residue_id

    def select_mutation(self):
        """Ask which residue to mutate to.
        """

        residue = self.m.fetch_residue(idx=self.mut_resid, chain=self.mut_chain)
        if residue.moltype == 'protein':
            aa = self.select_aa_mutation(residue)
        elif residue.moltype in ['dna', 'rna']:
            aa = self.select_nuc_mutation(residue)
        return aa

    def select_aa_mutation(self, residue):
        """Selection for protein residues.
        """

        _check_residue_name(residue)
        print('\nSelect new amino acid for %s-%s: ' % (residue.id, residue.resname))
        sys.stdout.write('Three- or one-letter code (or four-letter for ff specific residues): ')
        if residue.resname in ['HIE', 'HISE', 'HSE']:
            rol = 'X'
        elif residue.resname in ['HIP', 'HISH', 'HSP']:
            rol = 'Z'
        elif residue.resname in ['GLH', 'GLUH', 'GLUP']:
            rol = 'J'
        elif residue.resname in ['ASH', 'ASPH', 'ASPP']:
            rol = 'B'
        elif residue.resname in ['LYN', 'LYS', 'LSN']:
            rol = 'O'
        else:
            rol = library._one_letter[residue.resname]
        aa = None
        ol = list(library._aacids_dic.keys())
        tl = list(library._aacids_dic.values())
        ffpathlower = self.ffpath.lower()
        if 'amber' in ffpathlower:
                ol = list(library._aacids_ext_amber.keys())
                tl = list(library._aacids_ext_amber.values())
        if 'opls' in ffpathlower:
                ol = list(library._aacids_ext_oplsaa.keys())
                tl = list(library._aacids_ext_oplsaa.values()) + ['ASPP', 'GLUP', 'LSN']
        if 'charmm' in ffpathlower:
                ol = list(library._aacids_ext_charmm.keys())
                tl = list(library._aacids_ext_charmm.values())

        while aa is None:
            aa = input().upper()
            # some special residues:
            #   CM - deprotonated cysteine
            #   YM - deprotonated tyrosine
            if aa == 'CM':
                sys.stdout.write('Special case for deprotonated residue')
            elif len(aa) != 1 and len(aa) != 3 and len(aa) != 4:
                sys.stdout.write('Nope!\nThree- or one-letter code (or four-letter for ff specific residues): ')
                aa = None
            elif (len(aa) == 1 and aa not in ol+['B', 'J', 'O', 'X', 'Z']) or \
                 (len(aa) == 3 and aa not in tl) or \
                 (len(aa) == 4 and aa not in tl):
                sys.stdout.write('Unknown aa "%s"!\nThree- or one-letter code (or four-letter for ff specific residues): ' % aa)
                aa = None
            if aa and (len(aa) == 3 or len(aa) == 4):
                aa = library._ext_one_letter[aa]
        print('Will apply mutation %s->%s on residue %s-%s'
              % (rol, aa, residue.resname, str(residue.id)))
        return aa

    def select_nuc_mutation(self, residue):
        """Selection for nucleic acids.
        """
        aa = None
        print('\nSelect new base for %s-%s: ' % (residue.id, residue.resname))
        sys.stdout.write('One-letter code: ')
        while aa is None:
            aa = input().upper()
            if residue.moltype == 'dna' and aa not in ['A', 'C', 'G', 'T']:
                sys.stdout.write('Unknown DNA residue "%s"!\nOne-letter code: ' % aa)
                aa = None
            elif residue.moltype == 'rna' and aa not in ['A', 'C', 'G', 'U']:
                sys.stdout.write('Unknown RNA residue "%s"!\nOne-letter code: ' % aa)
                aa = None
            if aa:
                print('Will apply mutation %s->%s on residue %s-%d'
                      % (residue.resname[1], aa, residue.resname, residue.id))
            return aa


# =============
# Input Options
# =============
def parse_options():
    parser = argparse.ArgumentParser(description='''
This script applies mutations of residues in a structure file for subsequent
free energy calculations. It supports mutations to protein, DNA, and RNA
molecules.

The mutation information and dummy placements are taken from the hybrid residue
database "mutres.mtp". The best way to use this script is to take a pdb/gro file
that has been written with pdb2gmx with all hydrogen atoms present.

By default, all residues are renumbered starting from 1, so to have unique
residue IDs. If you want to keep the original residue IDs, you can use the flag
--keep_resid. In this case, you will also need to provide chain information
in order to be able to mutate the desired residue. Alternatively, if you would
like to use the original residue IDs but these have been changed, e.g. by gromacs,
you can provide a reference PDB file (with chain information too) using the --ref
flag. The input structure will be mutated according to the IDs chosen for the
reference structure after having mapped the two residue indices.

The program can either be executed interactively or via script. The script file
simply has to consist of "residue_id target_residue_name" pairs (just with some
space between the id and the name), or "chain_id residue_id target_residue_name"
if you are keeping the original residue IDs or providing a reference structure.

The script uses an extended one-letter code for amino acids to account for
different protonation states. Use the --resinfo flag to print the dictionary.

''', formatter_class=argparse.RawTextHelpFormatter)

    exclus = parser.add_mutually_exclusive_group()

    parser.add_argument('-f',
                        metavar='infile',
                        dest='infile',
                        type=str,
                        help='Input structure file in PDB or GRO format. '
                        'Default is "protein.pdb"',
                        default='protein.pdb')
    parser.add_argument('-fB',
                        metavar='infileB',
                        dest='infileB',
                        type=str,
                        help='Input structure file of the B state in PDB '
                        'or GRO format (optional).',
                        default=None)
    parser.add_argument('-o',
                        metavar='outfile',
                        dest='outfile',
                        type=str,
                        help='Output structure file in PDB or GRO format. '
                        'Default is "mutant.pdb"',
                        default='mutant.pdb')
    parser.add_argument('-ff',
                        metavar='ff',
                        dest='ff',
                        type=str.lower,
                        help='Force field to use. If none is provided, \n'
                        'a list of available ff will be shown.',
                        default=None)
    parser.add_argument('--script',
                        metavar='script',
                        dest='script',
                        type=str,
                        help='Text file with list of mutations (optional).',
                        default=None)
    exclus.add_argument('--keep_resid',
                        dest='renumber',
                        help='Whether to renumber all residues or to keep the\n'
                        'original residue IDs. By default, all residues are\n'
                        'renumbered so to have unique IDs. With this flags set,\n'
                        'the original IDs are kept. Because the IDs might not\n'
                        'be unique anymore, you will also be asked to choose\n'
                        'the chain ID where the residue you want to mutate is.',
                        default=True,
                        action='store_false')
    exclus.add_argument('--ref',
                        metavar='',
                        dest='ref_infile',
                        help='Provide a reference PDB structure from which to map\n'
                        'the chain and residue IDs onto the file to be mutated (-f).\n'
                        'This can be useful when wanting to mutate a file that\n'
                        'has had its residues renumbered or the chain information\n'
                        'removed (e.g. after gmx grompp). As in the --keep_resid\n'
                        'option, if --ref is chosen, you will need to provide chain\n'
                        'information either interactively or via the --script flag.',
                        default=None)
    parser.add_argument('--resinfo',
                        dest='resinfo',
                        help='Show the list of 3-letter -> 1-letter residues',
                        default=False,
                        action='store_true')

    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)

    # ------------------
    # residue dictionary
    # ------------------
    if args.resinfo is True:
        moltype = Model(args.infile).moltype
        if moltype == 'protein':
            print('\n ---------------------------')
            print(' Protein residues dictionary')
            print(' ---------------------------')
            _print_sorted_dict(library._ext_one_letter)
            print(' ---------------------------\n')
        elif moltype == 'dna':
            print('\n -----------------------')
            print(' DNA residues dictionary')
            print(' -----------------------')
            _print_sorted_dict(dna_one_letter)
            print(' -----------------------\n')
        elif moltype == 'rna':
            print('\n -----------------------')
            print(' RNA residues dictionary')
            print(' -----------------------')
            _print_sorted_dict(rna_one_letter)
            print(' -----------------------\n')
        else:
            print('\n ---------------------------')
            print(' Protein residues dictionary')
            print(' ---------------------------')
            _print_sorted_dict(library._ext_one_letter)
            print(' ---------------------------')
            print(' DNA residues dictionary')
            print(' ---------------------------')
            _print_sorted_dict(dna_one_letter)
            print(' ---------------------------')
            print(' RNA residues dictionary')
            print(' ---------------------------')
            _print_sorted_dict(rna_one_letter)
            print(' ---------------------------\n')
        exit()

    # ------------
    # ff selection
    # ------------
    if args.ff is None:
        args.ff = ff_selection()

    return args


def main(args):

    # input variables
    infile = args.infile
    infileB = args.infileB
    outfile = args.outfile
    ff = args.ff
    script = args.script
    renumber = args.renumber
    ref_infile = args.ref_infile

    # initialise Model
    m = Model(infile, renumber_residues=renumber, bPDBTER=True,
              rename_atoms=True, scale_coords='A')

    # if reference structure provided, initialise that Model too
    if ref_infile is not None:
        ref_m = Model(ref_infile, renumber_residues=False,
                      bPDBTER=True, rename_atoms=True, scale_coords='A')

    # if script is provided, do the mutations in that file
    # ----------------------------------------------------
    if script is not None:
        # 1) defualt: renumbering and not providing a ref struct
        if renumber is True and ref_infile is None:
            mutations_to_make = read_and_format(script, "is")
            # modify mut lists in mutations_to_make so that they have same
            # len of 3, both in the case where renumber is True or False
            mutations_to_make = [[None]+x for x in mutations_to_make]
        # 2) not renumbering and not providing a ref struct
        elif renumber is False and ref_infile is None:
            mutations_to_make = read_and_format(script, "sis")
        # 3) renumbering and providing a ref struct
        elif renumber is True and ref_infile is not None:
            # we read mutations according to numbering in the reference
            ref_mutations_to_make = read_and_format(script, "sis")
            mutations_to_make = []
            # we map these onto the file to be mutated
            for ref_mut in ref_mutations_to_make:
                m_resid = _match_mutation(m=m, ref_m=ref_m,
                                          ref_chain=ref_mut[0],
                                          ref_resid=ref_mut[1])
                # mut = [No Chain, resid in m, target residue]
                mut = [None, m_resid, ref_mut[2]]
                mutations_to_make.append(mut)

        # 4) NOT ALLOWED: providing a ref struct while not renumbering. This is
        # not allowed because we want to make sure there are unique indices for
        # each residue

        for mut in mutations_to_make:
            _check_residue_name(m.fetch_residue(idx=mut[1], chain=mut[0]))
            mutate(m=m,
                   mut_chain=mut[0],
                   mut_resid=mut[1],
                   mut_resname=mut[2],
                   ff=ff,
                   refB=infileB,
                   inplace=True,
                   verbose=True)

    # if script not provided, interactive selection
    # ---------------------------------------------
    else:
        do_more = True
        while do_more:
            # if no reference provided, choose from infile (Model m)
            if ref_infile is None:
                sele = InteractiveSelection(m=m, ff=ff, renumbered=renumber)
            # if reference IS provided, then choose from reference, then map
            # mutation onto infile (Model m)
            elif ref_infile is not None:
                sele = InteractiveSelection(m=ref_m, ff=ff, renumbered=False)
                sele.mut_resid = _match_mutation(m=m, ref_m=ref_m,
                                                 ref_chain=sele.mut_chain,
                                                 ref_resid=sele.mut_resid)
                sele.mut_chain = None

            # we have the needed info ==> carry out the mutation
            mutate(m=m,
                   mut_chain=sele.mut_chain,
                   mut_resid=sele.mut_resid,
                   mut_resname=sele.mut_resname,
                   ff=ff,
                   refB=infileB,
                   inplace=True,
                   verbose=True)

            # ask whether to do more mutations or stop
            if not _ask_next():
                do_more = False

    m.write(outfile)
    print('')
    print('mutations done...........')
    print('')


def entry_point():
    args = parse_options()
    main(args)


if __name__ == '__main__':
    entry_point()
>>>>>>> develop:pmx/extensions/pmx/scripts/mutate.py
