module SolvationInteractionPotential

using LinearAlgebra: norm
using MolecularMinimumDistances
using StaticArrays

import CellListMap: wrap_relative_to
import Chemfiles
import PDBTools

export FFType
export read_ffcharges, read_ffcharges!
export assign_charges
export electrostatic_potential

Base.@kwdef struct FFType
    type::String
    name::String
    resname::String
    charge::Float64
end

# Read charges from topology file, for the first time
function read_ffcharges(topology::String)
    ffcharges = FFType[]
    read_ffcharges!(ffcharges, topology)
    return ffcharges
end

# Append new data to the ffcharges array, given a new topology file
function read_ffcharges!(ffcharges::Vector{FFType}, topology::String)
    reading = false
    open(topology) do file
        for line in eachline(file)
            # start reading atoms section
            if occursin("[ atoms ]", line)
                reading = true
                continue
            end
            if !reading 
                continue
            end
            # end of atoms section 
            if length(line) == 0
                reading = false
                continue
            end
            # comment line on atom section
            type_data = split(line)
            if type_data[1] == ";"
                continue
            end
            # read atom type data
            fftype = FFType(
                type = type_data[2], 
                name = type_data[5], 
                resname = type_data[4], 
                charge = parse(Float64, type_data[7])
            )
            # if the atom type is not repeated, add to the list
            if isnothing(findfirst(isequal(fftype), ffcharges))
                push!(ffcharges, fftype)
            end
        end
    end
    return ffcharges
end

function assign_charges(atoms, ffcharges; warn_aliasing = true)
    charges = zeros(length(atoms))
    for (iat, at) in pairs(atoms)
        itype = findfirst( 
            type -> ( type.resname == at.resname && type.name == at.name ),
            ffcharges
        )
        if isnothing(itype)
            # checking for name aliasing for H atoms
            if length(at.name) == 4
                aliased_name = at.name[2:4]*at.name[1]
                itype = findfirst( 
                    type -> ( type.resname == at.resname && type.name == aliased_name ),
                    ffcharges
                )
                if !isnothing(itype)
                    warn_aliasing && println(" Warning: aliasing $(at.resname): $(at.name) to $(aliased_name).")
                end
            end
            if isnothing(itype)
                error("Could not find charges for atom: $(at.name) $(at.resname)")
            end
        end
        charges[iat] = ffcharges[itype].charge
    end
    return charges
end

function electrostatic_potential(
    solute::AbstractVector{PDBTools.Atom}, 
    solvent::AbstractVector{PDBTools.Atom}, 
    cutoff::Real,
    trajectory::String,
    topology_files::Vector{String}
)

    # Indexes of the atoms of the solute and solvent
    solute_indexes = [ at.index_pdb for at in solute ]
    solvent_indexes = [ at.index_pdb for at in solvent ]

    # Number of atoms of the solvent molecules
    natoms_per_molecule = length(PDBTools.select(solvent, "resnum $(solvent[1].resnum)"))

    # Read charges from topology files
    ffcharges = FFType[]
    for top_file in topology_files
        read_ffcharges!(ffcharges, top_file)
    end

    # Assign charges
    solute_charges = assign_charges(solute, ffcharges; warn_aliasing = false)
    solvent_charges = assign_charges(solvent, ffcharges; warn_aliasing = false)

    # Open trajectory with Chemfiles
    traj = Chemfiles.Trajectory(trajectory)

    u = 0.
    for iframe in 1:1
        frame = read(traj)
        coor = reinterpret(reshape, SVector{3,Float64}, Chemfiles.positions(frame))
        unit_cell = Chemfiles.lengths(Chemfiles.UnitCell(frame))
    
        # Compute electrostatic potential
        u += electrostatic_potential(
            solute_indexes,
            solvent_indexes,
            natoms_per_molecule,
            solute_charges,
            solvent_charges,
            coor,
            unit_cell,
            cutoff
        )
    end

    return u
end


function electrostatic_potential(
    solute::AbstractVector{Int}, 
    solvent::AbstractVector{Int}, 
    natoms_per_molecule::Int, 
    solute_charges::AbstractVector{<:Real}, 
    solvent_charges::AbstractVector{<:Real}, 
    coor::AbstractVector{<:SVector}, 
    unit_cell, 
    cutoff::Real
)

    # Define simulation box in this frame
    box = Box(unit_cell, cutoff)

    # Solute coordinates
    xsolute = @view(coor[solute])

    # Solvent coordinates
    xsolvent = @view(coor[solvent])

    # Compute the list of minimum-distances within cutoff
    list = minimum_distances(xsolute, xsolvent, natoms_per_molecule, box)

    # Compute the electrostatic potential between the solute and the solvent 
    # molecules that have some atom within the cutoff 
    electrostatic_potential = 0.
    for (isolvent,md) in pairs(list)
        # skipt if the solvent molecule is not within the cuotff
        if !md.within_cutoff 
            continue
        end
        ifirst = (isolvent-1)*natoms_per_molecule + 1
        ilast = ifirst + natoms_per_molecule - 1
        for isolvent = ifirst:ilast
            x = xsolvent[isolvent]
            qsolvent = solvent_charges[isolvent]
            for (isolute,y) in pairs(xsolute)
                qsolute = solute_charges[isolute]
                y = wrap_relative_to(y,x,box) 
                d = norm(x-y)
                electrostatic_potential += qsolvent*qsolute / d
            end
        end
    end
    return 332.05382e0 * electrostatic_potential
end

end # module
