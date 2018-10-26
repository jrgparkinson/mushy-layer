# This script generates symmetric fixed grids for Chombo, consisting of squares
# nested inside each other
# Set min_box_size equal to main.min_box_size in order for code to work in parallel

import math

max_level = 1
coarse_resolutions = [16, 32, 64, 128]
min_box_size = 8

for coarse_resolution in coarse_resolutions:
    grid_file = 'grids/grids-minBox' + str(min_box_size) + '-res'

    for lev in range(0, max_level+1):
        grid_file = grid_file + '-' + str(coarse_resolution * (pow(2, lev)))

    f = open(grid_file, 'w')
    f.write(str(max_level+1) + '\n')

    # First fine level of grids
    for lev in range(1, max_level+1):
        res = coarse_resolution * pow(2, lev)

        covered_area_size = int(res*pow(0.5,lev)) # because each level covers half the width of the previous level
        actual_box_size = min(min_box_size, covered_area_size)

        num_rows = max(res/(2*actual_box_size), 1)
        num_grids = pow(num_rows, 2)
        f.write(str(num_grids) + '\n')

        offset = int(0.5*(1-1/float(2**lev))*res)

        for grid in range(0, num_grids):
            row = math.floor(grid/num_rows)
            column = grid % num_rows
            bottom_x = offset + column * actual_box_size
            bottom_y = offset + row * actual_box_size
            top_x = offset + ((column+1) * actual_box_size) - 1
            top_y = offset + ((row+1) * actual_box_size) - 1
            line = '((%d, %d) (%d, %d) (0,0)) \n' % (bottom_x, bottom_y, top_x, top_y)
            f.write(line)

    f.close()
