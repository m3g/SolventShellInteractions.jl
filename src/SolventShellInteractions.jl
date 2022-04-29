module SolventShellInteractions

using LinearAlgebra: norm
using MolecularMinimumDistances
using StaticArrays
using ProgressMeter
using FastPow

import CellListMap: wrap_relative_to
import Chemfiles
import PDBTools

export FFType
export read_forcefield, read_forcefield!
export assign_forcefield
export nonbonded

Base.@kwdef struct FFType
    type::String
    name::String
    resname::String
    charge::Float64
    eps::Float64
    sig::Float64
end

# Read charges from topology file, for the first time
function read_forcefield(topology::String)
    ff = FFType[]
    read_forcefield!(ff, topology)
    return ff
end

# Append new data to the ffcharges array, given a new topology file
function read_forcefield!(ff::Vector{FFType}, topology::String)
    # Read charges and atom types
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
                charge = parse(Float64, type_data[7]),
                eps = -1.0,
                sig = -1.0,
            )
            # if the atom type is not repeated, add to the list
            if isnothing(findfirst(isequal(fftype), ff))
                push!(ff, fftype)
            end
        end
    end

    # Read LJ parameters, given classes
    reading = false
    open(topology) do file
        for line in eachline(file)
            # start reading atoms section
            if occursin("[ atomtypes ]", line)
                reading = true
                continue
            end
            if !reading
                continue
            end
            # end of atoms section 
            type_data = split(strip(line))
            if length(type_data) == 0
                reading = false
                continue
            end
            # comment line on atom section
            if type_data[1] == ";"
                continue
            end
            ilast = findfirst(isequal(";"), type_data)
            if isnothing(ilast)
                ilast = length(type_data) + 1
            end
            for (itype, atom_type) in pairs(ff)
                if atom_type.sig == -1
                    if atom_type.type == type_data[1]
                        ff[itype] = FFType(
                            type = ff[itype].type,
                            name = ff[itype].name,
                            resname = ff[itype].resname,
                            charge = ff[itype].charge,
                            eps = parse(Float64, type_data[ilast-1]),
                            sig = 10 * parse(Float64, type_data[ilast-2]),
                        )
                    end
                end
            end
        end
    end
    return ff
end

function assign_forcefield(atoms, ff; warn_aliasing = true)
    params = Vector{typeof((q = 0.0, eps = 0.0, sig = 0.0))}(undef, length(atoms))
    for (iat, at) in pairs(atoms)
        itype = findfirst(type -> (type.resname == at.resname && type.name == at.name), ff)
        if isnothing(itype)
            # checking for name aliasing for H atoms
            if length(at.name) == 4
                aliased_name = at.name[2:4] * at.name[1]
                itype = findfirst(
                    type -> (type.resname == at.resname && type.name == aliased_name),
                    ff,
                )
                if !isnothing(itype)
                    warn_aliasing && println(
                        " Warning: aliasing $(at.resname): $(at.name) to $(aliased_name).",
                    )
                end
            end
            if isnothing(itype)
                error("Could not find parameters for atom: $(at.name) $(at.resname)")
            else
                if ff[itype].charge == -1
                    error("Could not find charges for atom: $(at.name) $(at.resname)")
                end
                if ff[itype].eps == -1 || ff[itype].sig == -1
                    error(
                        "Could not find Lennard-Jones parameters for atom: $(at.name) $(at.resname)",
                    )
                end
            end
        end
        params[iat] = (q = ff[itype].charge, eps = ff[itype].eps, sig = ff[itype].sig)
    end
    return params
end

function nonbonded(
    solute::AbstractVector{PDBTools.Atom},
    solvent::AbstractVector{PDBTools.Atom},
    cutoff::Real,
    trajectory::String,
    topology_file::String;
    standard_cutoff::Bool = false,
    shift::Bool = false,
    show_progress::Bool = true,
)

    # Indexes of the atoms of the solute and solvent
    solute_indexes = [at.index_pdb for at in solute]
    solvent_indexes = [at.index_pdb for at in solvent]

    # Number of atoms of the solvent molecules
    natoms_per_molecule = length(PDBTools.select(solvent, "resnum $(solvent[1].resnum)"))

    # Read charges from topology files
    ff = read_forcefield(topology_file)

    # Assign charges
    solute_params = assign_forcefield(solute, ff, warn_aliasing = false)
    solvent_params = assign_forcefield(solvent, ff, warn_aliasing = false)

    # Open trajectory with Chemfiles
    traj = Chemfiles.Trajectory(trajectory)

    coordination_number = zeros(Int, length(traj))
    electrostatic_potential = zeros(length(traj))
    lennard_jones = zeros(length(traj))
    show_progress && (p = Progress(length(traj)))
    for iframe = 1:length(traj)
        frame = read(traj)
        coor = reinterpret(reshape, SVector{3,Float64}, Chemfiles.positions(frame))
        unit_cell = Chemfiles.lengths(Chemfiles.UnitCell(frame))

        # Compute electrostatic potential
        coordination_number[iframe],
        electrostatic_potential[iframe],
        lennard_jones[iframe] = nonbonded(
            solute_indexes,
            solvent_indexes,
            natoms_per_molecule,
            solute_params,
            solvent_params,
            coor,
            unit_cell,
            cutoff;
            standard_cutoff = standard_cutoff,
            shift = shift,
        )

        show_progress && next!(p)
    end

    return coordination_number, electrostatic_potential, lennard_jones
end

function lj(eps_x, eps_y, sig_x, sig_y, d; combination_rule = :geometric)
    if combination_rule == :geometric
        sig = sqrt(sig_x * sig_y)
    else
        sig = 0.5 * (sig_x + sig_y)
    end
    lj = @fastpow 4 * sqrt(eps_x * eps_y) * ((sig / d)^12 - (sig / d)^6)
    return lj
end

function nonbonded(
    solute::AbstractVector{Int},
    solvent::AbstractVector{Int},
    natoms_per_molecule::Int,
    solute_params::AbstractVector{<:NamedTuple},
    solvent_params::AbstractVector{<:NamedTuple},
    coor::AbstractVector{<:SVector},
    unit_cell,
    cutoff::Real;
    standard_cutoff::Bool = false,
    shift::Bool = false,
    combination_rule = :geometric,
)

    # Define simulation box in this frame
    box = Box(unit_cell, cutoff)

    # Solute coordinates
    xsolute = @view(coor[solute])

    # Solvent coordinates
    xsolvent = @view(coor[solvent])

    # Compute the list of minimum-distances within cutoff
    list = minimum_distances(xsolvent, xsolute, natoms_per_molecule, box)

    # Compute the electrostatic potential between the solute and the solvent 
    # molecules that have some atom within the cutoff 
    nbatches = Threads.nthreads()
    coordination_number = zeros(Int, nbatches)
    electrostatic_potential = zeros(nbatches)
    lennard_jones = zeros(nbatches)
    Threads.@threads for ibatch = 1:nbatches
        for isolvent = ibatch:nbatches:length(list)
            # skipt if the solvent molecule is not within the cuotff
            if !list[isolvent].within_cutoff
                continue
            end
            coordination_number[ibatch] += 1
            ifirst = (isolvent - 1) * natoms_per_molecule + 1
            ilast = ifirst + natoms_per_molecule - 1
            for iatom_solvent = ifirst:ilast
                x = xsolvent[iatom_solvent]
                qx, eps_x, sig_x = solvent_params[iatom_solvent]
                for (iatom_solute, y) in pairs(xsolute)
                    qy, eps_y, sig_y = solute_params[iatom_solute]
                    y = wrap_relative_to(y, x, box)
                    d = norm(y - x)
                    if d < box.cutoff || !standard_cutoff
                        qpair = qx * qy / d
                        ljpair = lj(
                            eps_x,
                            eps_y,
                            sig_x,
                            sig_y,
                            d;
                            combination_rule = combination_rule,
                        )
                        if shift
                            qpair -= qx * qy / box.cutoff
                            ljpair -= lj(
                                eps_x,
                                eps_y,
                                sig_x,
                                sig_y,
                                box.cutoff;
                                combination_rule = combination_rule,
                            )
                        end
                        electrostatic_potential[ibatch] += qpair
                        lennard_jones[ibatch] += ljpair
                    end
                end
            end
        end
    end
    return (
        sum(coordination_number),
        (332.05382e0 * 4.184) * sum(electrostatic_potential),
        sum(lennard_jones),
    )

end

include("./testing.jl")

end # module
