
function naive_electrostatic_potential(
    solute::AbstractVector{PDBTools.Atom}, 
    solvent::AbstractVector{PDBTools.Atom}, 
    cutoff::Real,
    trajectory::String,
    topology_files::Vector{String};
    standard_cutoff::Bool = false,
    switch::Bool = false,
    show_progress::Bool = true
)

    # Indexes of the atoms of the solute and solvent
    solute_indexes = [ at.index_pdb for at in solute ]
    solvent_indexes = [ at.index_pdb for at in solvent ]

    # Number of atoms of the solvent molecules
    natoms_per_molecule = length(PDBTools.select(solvent, "resnum $(solvent[1].resnum)"))
    nmolecules = length(solvent) ÷ natoms_per_molecule

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

    u = zeros(length(traj))
    show_progress && (p = Progress(length(traj)))
    for iframe in 1:length(traj)
        frame = read(traj)
        coor = reinterpret(reshape, SVector{3,Float64}, Chemfiles.positions(frame))
        unit_cell = Chemfiles.lengths(Chemfiles.UnitCell(frame))
        box = Box(unit_cell, cutoff)

        xsolute = @view(coor[solute_indexes])
        xsolvent = @view(coor[solvent_indexes])
    
        j = 0
        for _ = 1:nmolecules
            qpair = 0.
            dmin = +Inf
            for _ = 1:natoms_per_molecule
                j += 1
                y = xsolvent[j]
                qsolvent = solvent_charges[j]
                for i in eachindex(xsolute)
                    x = xsolute[i]
                    qsolute = solute_charges[i]
                    y = wrap_relative_to(y, x, box)
                    d = norm(y - x)
                    dmin = min(d,dmin)
                    if d < box.cutoff || !standard_cutoff
                        qpair += qsolute * qsolvent / d
                        if switch
                            qpair -= qsolute * qsolvent / box.cutoff
                        end
                    end
                end
            end
            if dmin < box.cutoff
                u[iframe] += qpair
            end
        end

        show_progress && next!(p)
    end

    return (332.05382e0 * 4.184)*u
end
