def calc_torsions(pose):
    phi_angles = []
    psi_angles = []

    for residx in range(1, pose.size()+1):
        phi = pose.phi(resi)
        psi = pose.psi(resi)
        phi_angles.append(phi)
        psi_angles.append(psi)

    if len(phi_angles) != len(psi_angles):
        raise Exception("Phi/Psi dimension mismatch")

    return {
        "phi": phi_angles,
        "psi": psi_angles
    }
