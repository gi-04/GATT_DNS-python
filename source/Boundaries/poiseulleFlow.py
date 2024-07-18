import numpy as np

# Define flow region
# Find which nodes will actually contain a flow and which ones will be in or at a wall
flow_region = np.ones((mesh.nx, mesh.ny, mesh.nz), dtype=bool)

# Add cavities to the flow region
if hasattr(flow_type, 'cav'):
    for cavity in flow_type.cav:
        x, y, z = cavity.x, cavity.y, cavity.z
        flow_region[(mesh.X > x[0]) & (mesh.X < x[1]) & 
                    (mesh.Y > y[0]) & (mesh.Y < y[1]) & 
                    (mesh.Z > z[0]) & (mesh.Z < z[1])] = True

# Remove roughnesses from the flow
if hasattr(flow_type, 'rug'):
    for roughness in flow_type.rug:
        x, y, z = roughness.x, roughness.y, roughness.z
        flow_region[(mesh.X >= x[0]) & (mesh.X <= x[1]) & 
                    (mesh.Y >= y[0]) & (mesh.Y <= y[1]) & 
                    (mesh.Z >= z[0]) & (mesh.Z <= z[1])] = False

# Add outer walls
flow_region[:, [0, -1], :] = False

# Get walls
find_walls_for_boundaries()

# Add walls to boundary conditions
wall_directions = ['xi', 'xf', 'yi', 'yf', 'zi', 'zf']
wall_positions = [wall_front_limits, wall_back_limits, wall_up_limits, 
                  wall_down_limits, wall_right_limits, wall_left_limits]

for wall_dir, wall_position in zip(wall_directions, wall_positions):
    for wall in wall_position:
        var.append('p')
        type.append('neu')
        dir.append(wall_dir)
        val.append(0)
        xi.append(wall[0])
        xf.append(wall[1])
        yi.append(wall[2])
        yf.append(wall[3])
        zi.append(wall[4])
        zf.append(wall[5])

# Add regions that are inside walls
for wall in inside_walls:
    for variable in ['p', 'u', 'v', 'w']:
        var.append(variable)
        type.append('dir')
        dir.append('yi')
        val.append(P0 if variable == 'p' else 0)
        xi.append(wall[0])
        xf.append(wall[1])
        yi.append(wall[2])
        yf.append(wall[3])
        zi.append(wall[4])
        zf.append(wall[5])
    
    if hasattr(flow_type, 'tWallRelative'):
        e_wall = E0 * flow_type.tWallRelative
    else:
        e_wall = E0
    
    var.append('e')
    type.append('dir')
    dir.append('yi')
    val.append(e_wall)
    xi.append(wall[0])
    xf.append(wall[1])
    yi.append(wall[2])
    yf.append(wall[3])
    zi.append(wall[4])
    zf.append(wall[5])

# Add moving walls
for y_pos, y_val in zip([1, mesh.ny], [mesh.Y[0], mesh.Y[-1]]):
    var.append('u')
    type.append('dir')
    dir.append('yi' if y_pos == 1 else 'yf')
    val.append(y_val / (mesh.Y[-1] - mesh.Y[0]))
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(y_pos)
    yf.append(y_pos)
    zi.append(1)
    zf.append(mesh.nz)


