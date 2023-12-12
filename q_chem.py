import psi4
import ase.io
import time


psi4.set_memory(20 * 1000 * 1000 * 1000)

def get_molecule(file,index):

    molecule = ase.io.read(file,index=str(index))

    formatted  = ""
    for atom in molecule:
        formatted  += atom.symbol + " " + " ".join(str(x) for x in atom.position) + "\n"

    formatted += """\n
                 symmetry c1
                 no_reorient
                 no_com
                 """

    return formatted


def get_trimer_dimer(file,index):

    molecule = ase.io.read(file,index=str(index))

    m1  = ""
    for atom in [0,1,6]:
        m1  += molecule[atom].symbol + " " + " ".join(str(x) for x in molecule[atom].position) + "\n"
    m2  = ""
    for atom in [2,3,7]:
        m2  += molecule[atom].symbol + " " + " ".join(str(x) for x in molecule[atom].position) + "\n"
    m3 = ""
    for atom in [4,5,8]:
        m3  += molecule[atom].symbol + " " + " ".join(str(x) for x in molecule[atom].position) + "\n"


    trimer = m1 + "--\n" + m2 + "--\n" + m3
    dimers = [ m1 + "--\n" + m2,
               m1 + "--\n" + m3,
               m2 + "--\n" + m3  ]

    trimer += """\n
                 symmetry c1
                 no_reorient
                 no_com
              """
    for i in range(3):
       dimers[i] = dimers[i] + """\n
                                 symmetry c1
                                 no_reorient
                                 no_com
                                """

    return trimer, dimers

###################################


def run_dft( index, dft_functional, molecule, output_file, be_quiet ):
    psi4.opt_coordinates='cartesian'
    psi4.geometry(molecule)
    psi4.set_options({
        'basis':    'aug-cc-pvtz',
        'e_convergence': 1e-10,
        'd_convergence': 1e-10,
    })

    
    psi4.set_output_file(output_file)
    if be_quiet:
        psi4.core.be_quiet()
    tic = time.perf_counter()
    out = psi4.energy('scf', dft_functional=dft_functional)
    toc = time.perf_counter()

    psi4.core.clean()
    print( dft_functional +  f"{index}, {out:0.14f}, {toc - tic:0.4f}")



def run_dft_3b( index, dft_functional, trimer, dimers, output_file, be_quiet  ):
    psi4.opt_coordinates='cartesian'
    psi4.geometry(trimer)
    psi4.set_options({
        'basis':    'aug-cc-pvtz',
        'e_convergence': 1e-10,
        'd_convergence': 1e-10,
    })

    psi4.set_output_file(output_file)
    if be_quiet:
        psi4.core.be_quiet()
    tic     = time.perf_counter()
    full    = psi4.energy('scf', dft_functional=dft_functional, bsse_type='cp')
    toc     = time.perf_counter()
    elapsed = toc - tic

    interaction  = full
    intermediate = []
    for d in dimers:
        psi4.opt_coordinates='cartesian'
        psi4.geometry(d)
        psi4.set_options({
            'basis':    'aug-cc-pvtz',
            'e_convergence': 1e-10,
            'd_convergence': 1e-10,
        })
        
        if be_quiet:
            psi4.core.be_quiet()
        tic          = time.perf_counter()
        out          = psi4.energy('scf', dft_functional=dft_functional, bsse_type='cp')
        toc          = time.perf_counter()
        elapsed     += toc - tic
        interaction -= out
        intermediate.append(out)

    print(dft_functional +  f"{index}, {full:0.14f}, {intermediate[0]:0.14f}, {intermediate[1]:0.14f}, {intermediate[2]:0.14f}, {interaction:0.14f}, {elapsed:0.4f}, {toc - tic:0.4f}")





def run_ccsd( index, molecule, output_file, be_quiet ):
    psi4.geometry(molecule)
    psi4.set_options({
        'basis':    'aug-cc-pvtz',
        'e_convergence': 1e-10,
        'd_convergence': 1e-10,
    })

    psi4.set_output_file(output_file)
    if be_quiet:
        psi4.core.be_quiet()
    tic = time.perf_counter()
    out = psi4.energy('ccsd(t)')
    toc = time.perf_counter()

    print(f"{index}, {out:0.14f}, {toc - tic:0.4f}")



def run_ccsd_3b( index, trimer, dimers, output_file, be_quiet ):
    psi4.geometry(trimer)
    psi4.set_options({
        'e_convergence': 1e-10,
        'd_convergence': 1e-10,
    })

    psi4.set_output_file(output_file)
    if be_quiet:
        psi4.core.be_quiet()
    tic     = time.perf_counter()
    full    = psi4.energy('ccsd(t)/aug-cc-pvtz', bsse_type='cp')
    toc     = time.perf_counter()
    elapsed = toc - tic

    interaction  = full
    intermediate = []
    for d in dimers:
        psi4.geometry(d)
        psi4.set_options({
            'e_convergence': 1e-10,
            'd_convergence': 1e-10,
        })

        if be_quiet:
            psi4.core.be_quiet()
        tic          = time.perf_counter()
        out          = psi4.energy('ccsd(t)/aug-cc-pvtz', bsse_type='cp')
        toc          = time.perf_counter()
        elapsed     += toc - tic
        interaction -= out
        intermediate.append(out)

    print(f"{index}, {full:0.14f}, {intermediate[0]:0.14f}, {intermediate[1]:0.14f}, {intermediate[2]:0.14f}, {interaction:0.14f}, {elapsed:0.4f}, {toc - tic:0.4f}" )


