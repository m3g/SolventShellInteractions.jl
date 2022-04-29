using SolvationInteractionPotential
using Test
using PDBTools

dir="./files"

@testset "SolvationInteractionPotential.jl" begin

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

    # compute standard electrostatic potential without switching
    u = electrostatic_potential(
        solute,
        solvent,
        cutoff,
        trajectory, 
        top_files;
        standard_cutoff = true,
        switch = false
    )

    u_naive = SolvationInteractionPotential.naive_electrostatic_potential(
           solute,
           solvent,
           cutoff,
           trajectory, 
           top_files;
           switch = false
    )
    @test u ≈ u_naive

    # compute standard electrostatic potential with switching
    u = electrostatic_potential(
        solute,
        solvent,
        cutoff,
        trajectory, 
        top_files;
        standard_cutoff = true,
        switch = true
    )

    u_naive = SolvationInteractionPotential.naive_electrostatic_potential(
           solute,
           solvent,
           cutoff,
           trajectory, 
           top_files;
           switch = true
    )
    @test u ≈ u_naive

end
