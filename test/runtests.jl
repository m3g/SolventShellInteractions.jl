using SolventShellInteractions
using Test
using PDBTools

dir="./files"

@testset "Reading files" begin

    pdb = readPDB("$dir/simulacao_EMIMDCA.pdb")
    solute = select(pdb, "protein")
    solvent = select(pdb, "resname EMI")

    # Indexes of the atoms of the solute and solvent
    solute_indexes = [ at.index_pdb for at in solute ]
    solvent_indexes = [ at.index_pdb for at in solvent ]
   
    # Number of atoms of the solvent molecules
    natoms_per_molecule = length(PDBTools.select(solvent, "resnum $(solvent[1].resnum)"))
    nmolecules = length(solvent) ÷ natoms_per_molecule
   
    # Read charges from topology files
    topology_files = [ "$dir/topol.top", "$dir/tip3p.itp" ]
    ffcharges = FFType[]
    for top_file in topology_files
        read_ffcharges!(ffcharges, top_file)
    end

    @test ffcharges[begin] == FFType("N3A", "N3A", "NC", -0.581)
    @test ffcharges[end] == FFType("opls_112", "HW2", "SOL", 0.417)
   
    # Assign charges
    solute_charges = assign_charges(solute, ffcharges; warn_aliasing = false)
    solvent_charges = assign_charges(solvent, ffcharges; warn_aliasing = false)

    @test solute_charges[begin] == -0.3
    @test solute_charges[end] == -0.8
    @test solvent_charges[begin] == 0.04
    @test solvent_charges[end] == -0.088
    @test solvent_charges[1:natoms_per_molecule] == solvent_charges[natoms_per_molecule+1:2*natoms_per_molecule]

end

@testset "Computing" begin

    # read pdb file
    pdb = readPDB("$dir/simulacao_EMIMDCA.pdb")

    # solute atoms
    solute = select(pdb, "protein")
    solvent = select(pdb, "resname EMI")

    # trajectory file (Gromacs xtc only)
    trajectory = "$dir/simulacao_EMIMDCA_curta.xtc"

    # topology files
    top_files = [ "$dir/topol.top", "$dir/tip3p.itp" ]

    # distance of the first dip in the distribution
    cutoff = 10.

    for standard_cutoff in [true, false], switch in [true, false]
        # compute standard electrostatic potential without switching
        u = electrostatic_potential(
            solute,
            solvent,
            cutoff,
            trajectory, 
            top_files;
            standard_cutoff = standard_cutoff,
            switch = switch,
            show_progress = false,
        )
        u_naive = SolventShellInteractions.naive_electrostatic_potential(
               solute,
               solvent,
               cutoff,
               trajectory, 
               top_files;
               standard_cutoff = standard_cutoff,
               switch = switch,
               show_progress = false,
        )
        @test u ≈ u_naive
    end

end
