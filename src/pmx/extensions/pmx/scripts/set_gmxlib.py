#!/usr/bin/env python
import os
from pmx import model
<<<<<<< HEAD:src/pmx/scripts/set_gmxlib.py
=======
from subprocess import call
>>>>>>> develop:pmx/extensions/pmx/scripts/set_gmxlib.py


def main():
    path = os.path.abspath(model.__file__)
    dir_path = os.path.dirname(path)
<<<<<<< HEAD:src/pmx/scripts/set_gmxlib.py
    gmxlib_proteins = os.path.join(dir_path, 'data/mutff45')
    gmxlib_nuc_acids = os.path.join(dir_path, 'data/mutff45dna')
=======
    gmxlib = os.path.join(dir_path, 'data/mutff')
>>>>>>> develop:pmx/extensions/pmx/scripts/set_gmxlib.py

    print('\n  In order to be able to use the hybrid/alchemical force fields \n'
          '  available in pmx, the environment variable GMXLIB needs to be set.\n')

<<<<<<< HEAD:src/pmx/scripts/set_gmxlib.py
    print('  The path to your pmx force field library for proteins is:')
    print('  %s\n' % gmxlib_proteins)
    print('  The path to your pmx force field library for nucleic acids is:')
    print('  %s\n' % gmxlib_nuc_acids)

    print('  Set the relevant GMXLIB path in your shell session as follows:')
    print('  $ export GMXLIB=%s\n' % gmxlib_proteins)
    print('  Or you can add this directly in your bashrc file.\n')

=======
    print('  The path to your pmx force field library is the following:')
    print('  %s\n' % gmxlib)

    print('  You can either set GMXLIB in your shell session as follows:')
    print('  $ export GMXLIB=%s\n' % gmxlib)
    print('  Or you can add this directly in your bashrc file.\n')

    set_gmxlib = input('  Do you wish pmx to set the GMXLIB variable in '
                       'your ~/.bashrc? [yes|no]\n  >>> ')

    while set_gmxlib not in ['yes', 'no']:
        set_gmxlib = input('  Please enter only either "yes" or "no"\n  >>> ')

    if set_gmxlib == 'yes':
        call('echo "# GMXLIB for using pmx forcefields" >> ~/.bashrc', shell=True)
        call('echo "export GMXLIB=%s" >> ~/.bashrc' % gmxlib, shell=True)
        print('  \n  All done. GMXLIB is set in your ~/.bashrc file.\n')
    elif set_gmxlib == 'no':
        print('  \n  All done. Do not forget to set GMXLIB any time you want to \n'
              '  use pmx or Gromacs with one of the hybrid forcefields in pmx.\n')

>>>>>>> develop:pmx/extensions/pmx/scripts/set_gmxlib.py

def entry_point():
    main()


if __name__ == '__main__':
    entry_point()
