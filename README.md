
<p align=center>
<img height=200px src="./docs/assets/logo.svg">
</p>

# SolventShellInteractions

Computes the non-bonded potential of the molecules which have at least
one atom within a distance from a solute. Used to evaluate the interactions of a solvent
with a solute, for the solvent molecules in coordination shells. The coordination
number is also computed.  

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
cn, q, lj = nonbonded(
    solute,
    solvent,
    cutoff,
    trajectory, 
    topology_file,
)

plot(
    [ cn q lj ],
    xlabel=[ "" "" "step" ],
    ylabel=[ "coordination\n number" "electrostatic\n energy / kJ / mol" "LJ energy / kJ / mol" ],
    linewidth=2, 
    framestyle=:box,
    labels=:none,
    layout=(3,1),
    size=(600,700),
)
```

Will produce (for a longer trajectory):

![example.png](./docs/example.png)

### Options

Some keyword arguments can be passed to the `nonbonded` function:

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

[CellListMap.jl](https://github.com/m3g/CellListMap.jl): Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff. Computer Physics Communications, 279, 108452, 2022. (https://doi.org/10.1016/j.cpc.2022.108452)

[Chemfiles.jl](https://github.com/chemfiles/Chemfiles.jl): A high-quality library for reading and writing trajectory files created by computational chemistry simulations programs. [10.5281/zenodo.1202207](https://doi.org/10.5281/zenodo.1202207).












