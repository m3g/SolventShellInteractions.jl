
function naive_nonbonded(
    solute::AbstractVector{PDBTools.Atom},
    solvent::AbstractVector{PDBTools.Atom},
    cutoff::Real,
    trajectory::String,
    topology_file::String;
    standard_cutoff::Bool = false,
    shift::Bool = false,
    show_progress::Bool = true,
    combination_rule = :geometric,
)

    # Indexes of the atoms of the solute and solvent
    solute_indexes = [at.index_pdb for at in solute]
    solvent_indexes = [at.index_pdb for at in solvent]

    # Number of atoms of the solvent molecules
    natoms_per_molecule = length(PDBTools.select(solvent, "resnum $(solvent[1].resnum)"))
    nmolecules = length(solvent) ÷ natoms_per_molecule

    # Read charges from topology files
    ff = read_forcefield(topology_file)

    # Assign charges
    solute_params = assign_forcefield(solute, ff; warn_aliasing = false)
    solvent_params = assign_forcefield(solvent, ff; warn_aliasing = false)

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
        box = Box(unit_cell, cutoff)

        xsolute = @view(coor[solute_indexes])
        xsolvent = @view(coor[solvent_indexes])

        j = 0
        for _ = 1:nmolecules
            qpair = 0.0
            ljpair = 0.0
            dmin = +Inf
            for _ = 1:natoms_per_molecule
                j += 1
                y = xsolvent[j]
                qy, eps_y, sig_y = solvent_params[j]
                for i in eachindex(xsolute)
                    x = xsolute[i]
                    qx, eps_x, sig_x = solute_params[i]
                    y = wrap_relative_to(y, x, box)
                    d = norm(y - x)
                    dmin = min(d, dmin)
                    if d < box.cutoff || !standard_cutoff
                        qpair += qx * qy / d
                        ljpair += lj(
                            eps_x,
                            eps_y,
                            sig_x,
                            sig_y,
                            d,
                            combination_rule = combination_rule,
                        )
                        if shift
                            qpair -= qx * qy / box.cutoff
                            ljpair -= lj(
                                eps_x,
                                eps_y,
                                sig_x,
                                sig_y,
                                box.cutoff,
                                combination_rule = combination_rule,
                            )
                        end
                    end
                end
            end
            if dmin < box.cutoff
                coordination_number[iframe] += 1
                electrostatic_potential[iframe] += qpair
                lennard_jones[iframe] += ljpair
            end
        end

        show_progress && next!(p)
    end

    return coordination_number,
    (332.05382e0 * 4.184) * electrostatic_potential,
    lennard_jones

end
