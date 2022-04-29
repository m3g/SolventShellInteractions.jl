using SolventShellInteractions
using Test
using PDBTools

dir="./files"

@testset "SolventShellInteractions.jl" begin

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
            switch = switch
        )

        u_naive = SolventShellInteractions.naive_electrostatic_potential(
               solute,
               solvent,
               cutoff,
               trajectory, 
               top_files;
               standard_cutoff = standard_cutoff,
               switch = switch
        )
        @test u ≈ u_naive
    end

end
