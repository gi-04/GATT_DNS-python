
import numpy as np
import scipy.sparse as sp

def make_matrices(mesh, domain, boundary, num_methods):
    """
    This function creates the matrices that will be used to compute derivatives and filters
    One pair of LHS and RHS will be computed for each type of region in the domain
    """

    # Find uniform regions
    find_derivative_regions()

    # Get finite differences coefficients
    centered_stencil_lhs, centered_stencil_rhs, decentered_stencil_lhs, decentered_stencil_rhs = finite_difference_coefficients(num_methods.spatial_derivs)
    centered_stencil_lhs_b, centered_stencil_rhs_b, decentered_stencil_lhs_b, decentered_stencil_rhs_b = finite_difference_coefficients(num_methods.spatial_derivs_buffer)
    filter_stencil_lhs, filter_stencil_rhs, filter_decentered_stencil_lhs, filter_decentered_stencil_rhs = spatial_filter_coefficients(num_methods.spatial_filter_strength, num_methods.filter_borders)

    # If the Z direction is not unitary but too small to fit the filter stencil, change the stencil just for it
    if mesh.nz > 1 and mesh.nz < 2 * len(filter_stencil_rhs) - 1:
        filter_stencil_lhs_z = np.ones(mesh.nz)
        filter_stencil_rhs_z = np.ones(mesh.nz)
        filter_decentered_stencil_lhs_z = np.ones(mesh.nz)
        filter_decentered_stencil_rhs_z = np.ones(mesh.nz)
    else:
        filter_stencil_lhs_z = filter_stencil_lhs
        filter_stencil_rhs_z = filter_stencil_rhs
        filter_decentered_stencil_lhs_z = filter_decentered_stencil_lhs
        filter_decentered_stencil_rhs_z = filter_decentered_stencil_rhs

    # If the mesh in Z has exactly 4 nodes and has no boundaries, switch it to spectral mode
    if mesh.nz == 4 and len(deriv_starts_z) == 1 and len(deriv_starts_z[0]) == 0 and len(deriv_ends_z[0]) == 0:
        centered_stencil_lhs_z = np.ones(mesh.nz)
        centered_stencil_rhs_z = np.array([0, np.pi / 4])
        centered_stencil_lhs_bz = centered_stencil_lhs_z
        centered_stencil_rhs_bz = centered_stencil_rhs_z
    else:
        centered_stencil_lhs_z = centered_stencil_lhs
        centered_stencil_rhs_z = centered_stencil_rhs
        centered_stencil_lhs_bz = centered_stencil_lhs_b
        centered_stencil_rhs_bz = centered_stencil_rhs_b

    # Make matrices
    lhs_x, rhs_x = make_matrices_each_direction(centered_stencil_lhs, centered_stencil_rhs, decentered_stencil_lhs, decentered_stencil_rhs, deriv_starts_x, deriv_ends_x, mesh.nx, [])
    lhs_y, rhs_y = make_matrices_each_direction(centered_stencil_lhs, centered_stencil_rhs, decentered_stencil_lhs, decentered_stencil_rhs, deriv_starts_y, deriv_ends_y, mesh.ny, [])
    lhs_z, rhs_z = make_matrices_each_direction(centered_stencil_lhs_z, centered_stencil_rhs_z, decentered_stencil_lhs, decentered_stencil_rhs, deriv_starts_z, deriv_ends_z, mesh.nz, [])

    lhs_xb, rhs_xb = make_matrices_each_direction(centered_stencil_lhs_b, centered_stencil_rhs_b, decentered_stencil_lhs_b, decentered_stencil_rhs_b, deriv_starts_x, deriv_ends_x, mesh.nx, mesh.x.buffer)
    lhs_yb, rhs_yb = make_matrices_each_direction(centered_stencil_lhs_b, centered_stencil_rhs_b, decentered_stencil_lhs_b, decentered_stencil_rhs_b, deriv_starts_y, deriv_ends_y, mesh.ny, mesh.y.buffer)
    lhs_zb, rhs_zb = make_matrices_each_direction(centered_stencil_lhs_bz, centered_stencil_rhs_bz, decentered_stencil_lhs_b, decentered_stencil_rhs_b, deriv_starts_z, deriv_ends_z, mesh.nz, mesh.z.buffer)

    f_lhs_x, f_rhs_x = make_matrices_each_direction(filter_stencil_lhs, filter_stencil_rhs, filter_decentered_stencil_lhs, filter_decentered_stencil_rhs, deriv_starts_x, deriv_ends_x, mesh.nx, [])
    f_lhs_y, f_rhs_y = make_matrices_each_direction(filter_stencil_lhs, filter_stencil_rhs, filter_decentered_stencil_lhs, filter_decentered_stencil_rhs, deriv_starts_y, deriv_ends_y, mesh.ny, [])
    f_lhs_z, f_rhs_z = make_matrices_each_direction(filter_stencil_lhs_z, filter_stencil_rhs_z, filter_decentered_stencil_lhs_z, filter_decentered_stencil_rhs_z, deriv_starts_z, deriv_ends_z, mesh.nz, [])

    # Add buffer zones into the derivative matrices
    if not hasattr(num_methods, 'change_order_x') or num_methods.change_order_x:
        lhs_x, rhs_x = add_buffer_to_matrix(lhs_x, lhs_xb, rhs_x, rhs_xb, mesh.nx, mesh.x.buffer)
    if not hasattr(num_methods, 'change_order_y') or num_methods.change_order_y:
        lhs_y, rhs_y = add_buffer_to_matrix(lhs_y, lhs_yb, rhs_y, rhs_yb, mesh.ny, mesh.y.buffer)
    if not hasattr(num_methods, 'change_order_z') or num_methods.change_order_z:
        lhs_z, rhs_z = add_buffer_to_matrix(lhs_z, lhs_zb, rhs_z, rhs_zb, mesh.nz, mesh.z.buffer)

    # Transform due to mesh stretching

    # Chose the spatial derivative method for the metric
    if not hasattr(num_methods, 'metric_method') or num_methods.metric_method is None:  # If unspecified, use SL4
        num_methods.metric_method = 'SL4'
    elif num_methods.metric_method == '':  # If left blank, use same as spatial derivatives
        num_methods.metric_method = num_methods.spatial_derivs

    if mesh.nz == 4:  # Use SL4 for z if only 4 nodes are being used, due to its smaller stencil
        num_methods.metric_method_z = 'SL4'
    else:
        num_methods.metric_method_z = num_methods.metric_method

    lhs_x = apply_metric(lhs_x, mesh.X, mesh.x, domain.xf, num_methods.metric_method)
    lhs_y = apply_metric(lhs_y, mesh.Y, mesh.y, domain.yf, num_methods.metric_method)
    lhs_z = apply_metric(lhs_z, mesh.Z, mesh.z, domain.zf, num_methods.metric_method_z)

    # Remove ends of filters if needed
    if num_methods.filter_borders:
        if hasattr(num_methods, 'filter_borders_start_x') and not num_methods.filter_borders_start_x:
            for i in range(len(f_lhs_x)):
                f_lhs_x[i][:5, :] = 0
                f_rhs_x[i][:5, :] = 0
                f_lhs_x[i][:5, :5] = np.eye(5)
                f_rhs_x[i][:5, :5] = np.eye(5)
        if hasattr(num_methods, 'filter_borders_end_x') and not num_methods.filter_borders_end_x:
            for i in range(len(f_lhs_x)):
                f_lhs_x[i][-4:, :] = 0
                f_rhs_x[i][-4:, :] = 0
                f_lhs_x[i][-4:, -4:] = np.eye(5)
                f_rhs_x[i][-4:, -4:] = np.eye(5)
        if hasattr(num_methods, 'filter_borders_start_y') and not num_methods.filter_borders_start_y:
            for i in range(len(f_lhs_y)):
                f_lhs_y[i][:5, :] = 0
                f_rhs_y[i][:5, :] = 0
                f_lhs_y[i][:5, :5] = np.eye(5)
                f_rhs_y[i][:5, :5] = np.eye(5)
        if hasattr(num_methods, 'filter_borders_end_y') and not num_methods.filter_borders_end_y:
            for i in range(len(f_lhs_y)):
                f_lhs_y[i][-4:, :] = 0
                f_rhs_y[i][-4:, :] = 0
                f_lhs_y[i][-4:, -4:] = np.eye(5)
                f_rhs_y[i][-4:, -4:] = np.eye(5)
        if hasattr(num_methods, 'filter_borders_start_z') and not num_methods.filter_borders_start_z:
            for i in range(len(f_lhs_z)):
                f_lhs_z[i][:5, :] = 0
                f_rhs_z[i][:5, :] = 0
                f_lhs_z[i][:5, :5] = np.eye(5)
                f_rhs_z[i][:5, :5] = np.eye(5)
        if hasattr(num_methods, 'filter_borders_end_z') and not num_methods.filter_borders_end_z:
            for i in range(len(f_lhs_z)):
                f_lhs_z[i][-4:, :] = 0
                f_rhs_z[i][-4:, :] = 0
                f_lhs_z[i][-4:, -4:] = np.eye(5)
                f_rhs_z[i][-4:, -4:] = np.eye(5)

    # Save to output structure

    matrices.x.types = type_map_x
    matrices.y.types = type_map_y
    matrices.z.types = type_map_z

    matrices.x.lhs = lhs_x
    matrices.x.rhs = rhs_x
    matrices.y.lhs = lhs_y
    matrices.y.rhs = rhs_y
    matrices.z.lhs = lhs_z
    matrices.z.rhs = rhs_z

    matrices.x.flhs = f_lhs_x
    matrices.x.frhs = f_rhs_x
    matrices.y.flhs = f_lhs_y
    matrices.y.frhs = f_rhs_y
    matrices.z.flhs = f_lhs_z
    matrices.z.frhs = f_rhs_z

def make_matrices_each_direction(centered_stencil_lhs, centered_stencil_rhs, decentered_stencil_lhs, decentered_stencil_rhs, deriv_starts, deriv_ends, n, buffer_info=None):
    """
    This function creates the matrices for each direction
    """

    # Check if the mesh is large enough for the stencil
    if n != 1 and n < 2 * len(centered_stencil_rhs) - 1:
        raise ValueError(f"Mesh is not large enough for one of the stencils. It has {n} nodes but the stencil needs at least {2 * len(centered_stencil_rhs) - 1}.")

    n_types = len(deriv_starts)

    lhs = [None] * n_types
    rhs = [None] * n_types

    if n == 1:  # If single point in this direction, set derivative to zero and filter to one
        if centered_stencil_rhs[0] == 0:
            for i in range(n_types):
                lhs[i] = np.ones((n, n))
                rhs[i] = np.zeros((n, n))
        else:
            for i in range(n_types):
                lhs[i] = np.ones((n, n))
                rhs[i] = np.ones((n, n))
    else:
        # Create the base stencils
        lhs_base = sp.diags(centered_stencil_lhs[0] * np.ones(n), 0)
        rhs_base = sp.diags(centered_stencil_rhs[0] * np.ones(n), 0)

        invert_stencil = -1 if centered_stencil_rhs[0] == 0 else 1

        for i in range(1, len(centered_stencil_lhs)):
            inds = np.arange(n)[:, np.newaxis] + np.mod(np.arange(n) - 1 + i, n)
            lhs_base[inds] = centered_stencil_lhs[i]
            inds = np.arange(n)[:, np.newaxis] + np.mod(np.arange(n) - 1 - i, n)
            lhs_base[inds] = centered_stencil_lhs[i]

        for i in range(1, len(centered_stencil_rhs)):
            inds = np.arange(n)[:, np.newaxis] + np.mod(np.arange(n) - 1 + i, n)
            rhs_base[inds] = centered_stencil_rhs[i]
            inds = np.arange(n)[:, np.newaxis] + np.mod(np.arange(n) - 1 - i, n)
            rhs_base[inds] = invert_stencil * centered_stencil_rhs[i]

        # Check if this is a buffer zone that needs to be upwind
        # and add that to the list of starts and ends so that the decentered
        # stencil is used
        if buffer_info is not None:
            if hasattr(buffer_info.i, 'upwind') and buffer_info.i.upwind:
                for i in range(n_types):
                    deriv_starts[i] = np.unique(np.concatenate([deriv_starts[i], np.arange(1, buffer_info.i.n + 1)]), axis=0)[::-1]

            if hasattr(buffer_info.f, 'upwind') and buffer_info.f.upwind:
                for i in range(n_types):
                    deriv_ends[i] = np.unique(np.concatenate([deriv_ends[i], np.arange(n - buffer_info.f.n + 1, n + 1)]), axis=0)

        # Add startings and endings

        m_lhs, n_lhs = decentered_stencil_lhs.shape
        m_rhs, n_rhs = decentered_stencil_rhs.shape

        for i in range(n_types):
            lhs_temp = lhs_base.copy()
            rhs_temp = rhs_base.copy()

            for j in range(len(deriv_starts[i])):
                ind_start = deriv_starts[i][j]

                lhs_temp[ind_start:ind_start + m_lhs - 1, :] = 0
                rhs_temp[ind_start:ind_start + m_rhs - 1, :] = 0

                lhs_temp[ind_start:ind_start + m_lhs - 1, ind_start:ind_start + n_lhs - 1] = decentered_stencil_lhs
                rhs_temp[ind_start:ind_start + m_rhs - 1, ind_start:ind_start + n_rhs - 1] = decentered_stencil_rhs

            for j in range(len(deriv_ends[i])):
                ind_end = deriv_ends[i][j]

                lhs_temp[ind_end - m_lhs + 1:ind_end, :] = 0
                rhs_temp[ind_end - m_rhs + 1:ind_end, :] = 0

                lhs_temp[ind_end - m_lhs + 1:ind_end, ind_end - n_lhs + 1:ind_end] = np.flip(decentered_stencil_lhs, axis=0)
                rhs_temp[ind_end - m_rhs + 1:ind_end, ind_end - n_rhs + 1:ind_end] = invert_stencil * np.flip(decentered_stencil_rhs, axis=0)

            lhs[i] = lhs_temp
            rhs[i] = rhs_temp

    return lhs, rhs

def add_buffer_to_matrix(base_matrix_l, buffer_matrix_l, base_matrix_r, buffer_matrix_r, n, buffer_info):
    """
    This function adds the buffer zone to the matrices
    """

    ni = buffer_info.i.n
    nf = buffer_info.f.n

    if hasattr(buffer_info.i, 'transition'):
        nti = int(buffer_info.i.transition * ni)
    else:
        nti = ni
    if hasattr(buffer_info.f, 'transition'):
        ntf = int(buffer_info.f.transition * nf)
    else:
        ntf = nf

    ni1 = ni - nti + 1
    ni2 = ni
    nf1 = n - nf + 1
    nf2 = nf1 + ntf - 1

    # The buffer zone can be computed by a different type of derivatives. The transition is done smoothly.

    eta = np.ones((n, 1))

    if ni > 0:
        eta[:ni1] = 0
        eta[ni1:ni2] = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, ni2 - ni1 + 1))

    if nf > 0:
        eta[nf2:] = 0
        eta[nf1:nf2] = 0.5 + 0.5 * np.cos(np.linspace(0, np.pi, nf2 - nf1 + 1))

    n_types = len(base_matrix_l)
    new_matrix_l = [None] * n_types
    new_matrix_r = [None] * n_types

    eta = np.sqrt(eta)

    for i in range(n_types):
        new_matrix_l[i] = eta * base_matrix_l[i] + (1 - eta) * buffer_matrix_l[i]
        new_matrix_r[i] = eta * base_matrix_r[i] + (1 - eta) * buffer_matrix_r[i]

    return new_matrix_l, new_matrix_r

def apply_metric(lhs, X, mesh_info, xf, method):
    """
    This function applies the metric to the LHS matrix
    """

    # Scale the LHS matrix due to the mesh stretching. SL4 method is used here

    if mesh_info.periodic and mesh_info.fix_periodic_domain_size and len(X) > 1:  # Add temp node for a periodic mesh
        X = np.concatenate([X, [xf]])
        added_temp_node = True
    else:
        added_temp_node = False

    n = len(X)
    if n == 1:
        return

    centered_stencil_lhs, centered_stencil_rhs, decentered_stencil_lhs, decentered_stencil_rhs = finite_difference_coefficients(method)
    lhs_temp, rhs_temp = make_matrices_each_direction(centered_stencil_lhs, centered_stencil_rhs, decentered_stencil_lhs, decentered_stencil_rhs, [1], [n], n, [])

    dXdEta = lhs_temp @ X[:, np.newaxis]

    if added_temp_node:
        dXdEta[-1] = 0

    n_types = len(lhs)
    for i in range(n_types):
        lhs[i] = lhs[i] @ dXdEta


