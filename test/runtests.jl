using SolventShellInteractions
using Test
using PDBTools

dir="./files"

@testset "Reading files" begin

    pdb = readPDB("$dir/system.pdb")
    solute = select(pdb, "protein")
    solvent = select(pdb, "resname EMI")

    # Indexes of the atoms of the solute and solvent
    solute_indexes = [ at.index_pdb for at in solute ]
    solvent_indexes = [ at.index_pdb for at in solvent ]
   
    # Number of atoms of the solvent molecules
    natoms_per_molecule = length(PDBTools.select(solvent, "resnum $(solvent[1].resnum)"))
    nmolecules = length(solvent) ÷ natoms_per_molecule
   
    # Read force field parameters from topology file
    ff = read_forcefield("$dir/processed.top")

    @test ff[begin] == FFType("N3A", "N3A", "NC", -0.581, 1.0669, 3.25) 
    @test ff[end] == FFType("opls_403", "I", "I", -1.0, 0.29288, 5.4) 
   
    # Assign parameters to atoms
    solute_params = assign_forcefield(solute, ff; warn_aliasing = false)
    solvent_params = assign_forcefield(solvent, ff; warn_aliasing = false)

    @test solute_params[begin] == (q = -0.3, eps = 0.71128, sig = 3.25)
    @test solute_params[end] == (q = -0.8, eps = 0.87864, sig = 2.96) 
    @test solvent_params[begin] == (q = 0.04, eps = 0.29288, sig = 3.55) 
    @test solvent_params[end] == (q = -0.088, eps = 0.0, sig = 0.0)
    @test solvent_params[1:natoms_per_molecule] == solvent_params[natoms_per_molecule+1:2*natoms_per_molecule]

end

@testset "Computing" begin

    # read pdb file
    pdb = readPDB("$dir/system.pdb")

    # solute atoms
    solute = select(pdb, "protein")
    solvent = select(pdb, "resname EMI")

    # trajectory file (Gromacs xtc only)
    trajectory = "$dir/simulation_short.xtc"

    # topology files
    top_files = "$dir/processed.top"

    # distance of the first dip in the distribution
    cutoff = 10.

    for standard_cutoff in [true, false], shift in [true, false]
        # compute standard electrostatic potential without shifting
        q, lj = nonbonded(
            solute,
            solvent,
            cutoff,
            trajectory, 
            top_files;
            standard_cutoff = standard_cutoff,
            shift = shift,
            show_progress = false,
        )
        q_naive, lj_naive = SolventShellInteractions.naive_nonbonded(
               solute,
               solvent,
               cutoff,
               trajectory, 
               top_files;
               standard_cutoff = standard_cutoff,
               shift = shift,
               show_progress = false,
        )
        @test q ≈ q_naive
        @test lj ≈ lj_naive
    end

end
