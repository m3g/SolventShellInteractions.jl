# SolventShellInteractions

Computes the electrostatic potential of the molecules which have at least
one atom within a distance from a solute. 

Currently only works with Gromacs topology and trajectory formats. 

The calculations run in parallel if `julia` is initialized with multi-threading support, i. e. with `julia -t auto`, for example.

## Installation

```julia
julia> import Pkg

julia> Pkg.add(url="http://github.com/m3g/SolventShellInteractions.jl")
```

## Example

```julia
using Plots
using PDBTools
using SolventShellInteractions

dir="./test/files"

# read pdb file
pdb = readPDB("$dir/system.pdb")

# solute atoms
solute = select(pdb, "protein")
solvent = select(pdb, "resname EMI")

# trajectory file (Gromacs xtc only)
trajectory = "$dir/simulation_short.xtc"

# topology file
topology_file = "$dir/processed.top"

# distance of the first dip in the distribution
cutoff = 10.

# compute electrostatic potential
q, lj = nonbonded(
    solute,
    solvent,
    cutoff,
    trajectory, 
    topology_file,
)

plot(
    [ q lj ],
    labels=[ "electrostatic" "Lennard-Jones" ],
    xlabel="step",
    ylabel="potential energy / kJ / mol",
    linewidth=2, 
    framestyle=:box,
)
```

Will produce (for a longer trajectory):

![example.png](./docs/example.png)


### Options

Three keyword arguments can be passed to the `electrostatic_potential` function:

| keyword |  values | default |  meaning  | 
|:-------------:|:---------------:|:--:|:-------------|
| `show_progress` | `true/false` | `true` | Display or not the progress bar.  | 
| `combination_rule` | `:geometric/:arithmetic` | `:geometric` | Type of combination rule for the sigma parameters of the LJ potential.  | 
| `standard_cutoff` | `true/false` | `false` |Use a standard cutoff, and thus no interaction above the cutoff distance will be considered. | 
| `shift` | `true/false` | `false` | Use a shifting function, as described [here](https://www.ks.uiuc.edu/Research/namd/2.10/ug/node23.html).  | 

Note that with `standard_cutoff` and/or `shift` turned on, the calculation is not particularly interesting, as it does not represent anymore the interaction of the complete molecules that have at least one atom within the desired distance range. These options are mostly used for testing purposes.

## Related work:

[ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl):  Investigating the structure of solutions of complex-shaped molecules from a solvent-shell perspective J. Mol. Liq. 117945, 2021.

[MolecularMinimumDistances.jl](https://github.com/m3g/MolecularMinimumDistances.jl): Computes the set of minimum distances between to sets of particles, which can be grouped (like in molecules).  

## References:

This package uses the following libraries, which should be cited independently:

[CellListMap.jl](https://github.com/m3g/CellListMap.jl): Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff. https://doi.org/10.48550/arXiv.2202.06427.

[Chemfiles.jl](https://github.com/chemfiles/Chemfiles.jl): A high-quality library for reading and writing trajectory files created by computational chemistry simulations programs. [10.5281/zenodo.1202207](https://doi.org/10.5281/zenodo.1202207).












