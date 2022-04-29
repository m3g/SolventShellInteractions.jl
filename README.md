# SolventShellInteractions

Computes the electrostatic potential of the molecules which have at least
one atom within a distance from a solute. 

Currently only works with Gromacs topology and trajectory formats. 

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

# compute electrostatic potential
u = electrostatic_potential(
    solute,
    solvent,
    cutoff,
    trajectory, 
    top_files,
)

plot(
    u,
    xlabel="step",
    ylabel="electrostaic potential / kJ / mol",
    linewidth=2, framestyle=:box, label=nothing
)
```

Will produce (for a longer trajectory):

![example.png](./docs/example.png)


### Options

Three keyword arguments can be passed to the `electrostatic_potential` function:

| keyword |  values | default |  meaning  | 
|:-------------:|:---------------:|:--:|:-------------|
| `standard_cutoff` | `true/false` | `false` |Use a standard cutoff, and thus no interaction above the cutoff distance will be considered. | 
| `shift` | `true/false` | `false` | Use a shifting function, as described [here](https://www.ks.uiuc.edu/Research/namd/2.10/ug/node23.html).  | 
| `show_progress` | `true/false` | `true` | Display or not the progress bar.  | 

## Related work:

[ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl):  Investigating the structure of solutions of complex-shaped molecules from a solvent-shell perspective J. Mol. Liq. 117945, 2021.

[MolecularMinimumDistances.jl](https://github.com/m3g/MolecularMinimumDistances.jl): Computes the set of minimum distances between to sets of particles, which can be grouped (like in molecules).  

## References:

This package uses the following libraries, which should be cited independently:

[CellListMap.jl](https://github.com/m3g/CellListMap.jl): Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff. https://doi.org/10.48550/arXiv.2202.06427.

[Chemfiles.jl](https://github.com/chemfiles/Chemfiles.jl): A high-quality library for reading and writing trajectory files created by computational chemistry simulations programs. [10.5281/zenodo.1202207](https://doi.org/10.5281/zenodo.1202207).












