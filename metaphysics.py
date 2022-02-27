#################################
########## METAPHYSICS ##########
####### Liam Hale McCarty #######
#################################

#####
# _Metaphysics_ is a pair of sculptures. One is a domino made of dice, and the other is a die made of dominoes. The latter is the focus of the code here, since constructing it properly requires a series of fairly complex calculations.
#####

#####
### Goals
#
# This code accomplishes two key goals:
# 1. List the sequence of dominoes that cover the die.
# 2. Determine the minimum number of domino sets required to make the die.
# Moreover, it does so in essentially a "scale free" way: the same goals can be accomplished for larger (or smaller) versions of the sculptures by simply changing a few parameters.

# The latter goal (determining `min_num_sets`) necessitates a long sequence of intermediate calcuations. To help you follow the code below, here's a brief summary of the major steps:
# 1. Number the squares from 0 to `num_squares` - 1.
# 2. Calculate `domino` and `number` for a square from the square number `square`.
# 3. Define the dots with coordinates.
# 4. Define the collections of dots for each side.
# 5. Calculate the white areas from the collections of dots.
# 6. Convert the dots and white areas from local into global coordinates.
# 7. Find the full and half dominoes that constitute the dots and white areas.
# 8. Count how many full and half dominoes constitute them.
# 9. Determine `min_num_sets` from these counts.
#####

#####
### Background
#
# The die of dominoes is a cube tiled by a circular domino train in a Hilbert curve pattern. (To be technically precise, it's a polygonal approximation to a Hilbert curve, since the Hilbert curve itself is the infinite limit of such approximations.) The standard Hilbert curve is open and fills the unit square, such that many copies of it fill the real number plane ($\mathbb{R}^2$). But since a cube is topologically closed, the Hilbert curve in this case will also be closed in that it will "loop back on itself". 
#
# There are many ways one could describe such a Hilbert curve, but it's useful to (arbitrarily) pick out "start" and "end" points for clarity.
#
# 'start' of Hilbert curve                           Orientation of dots:
# v___________ ___________                            ___________ ___________
# |           |           |                          |           |           |
# |           |           |                          |           |        +  |
# |  Side     |  Side     |                          |     +     |           |
# |  1        |  2        |                          |           |  +        |
# |___________|___________|___________               |___________|___________|___________
#             |           |           |                          |           |           |
#             |           |           |                          |  +     +  |  +  +  +  |
#             |  Side     |  Side     |                          |           |           |
#             |  4        |  6        |                          |  +     +  |  +  +  +  |
#             |___________|___________|___________               |___________|___________|___________
#                         |           |           |                          |           |           |
#                         |           |           |                          |  +     +  |        +  |
#                         |  Side     |  Side     |                          |     +     |     +     |
#                         |  5        |  3        |                          |  +     +  |  +        |
#                         |___________|___________|                          |___________|___________|
#                                                 ^
#                            'end' of Hilbert curve
#
# Also shown above is the orientation of dots on the die. Note that this is quite a particular orientation. Not only is the order of numbers specific, but so is the way the dots are laid out on each side. For example, the 2-dot side could have dots on the other diagonal instead, but it doesn't. All of this is based on physical dice I used as references. A fun fact I learned in the process of inspecting and researching them: dice are designed such that two opposite sides' dots always sum to 7!
#
# The full Hilbert curve is composed of (traditional, open) Hilbert curves on each side of the cube. Each of these corresponds to either a "type 1" or "type 2" tiling (hence the 1 and 2 indices in the side names). 
#
# For a type 1 tiling, the Hilbert curve starts in the upper left and ends in the upper right corner. For a type 2 tiling, the curve starts in the upper left and ends in the lower left corner.
# 
# start     end                 start
# v___________v                 v___________
# |           |                 |           |
# |           |                 |           |
# |  Type 1   |         x       |  Type 2   |         y
# |           |    |---- >      |           |    |---- >
# |___________|    |            |___________|    |
#                y v            ^              x v
#                               end
#
# Note how this changes the local coordinate axes, as the diagram above indicates.
#
# These are both rotations of the "standard" Hilbert curve, which is the default for the `hilbertcurve` package I leverage below. That standard curve starts in the lower left and ends in the lower right corner:
#
#  ___________        
# |           |      
# |           |      
# |  Standard |  y ^ 
# |           |    | 
# |___________|    |---- >
# ^           ^        x
# start     end
#
# To keep track of things, it's helpful to label the sides and dots. I label the sides, intuitively enough, based on the number of dots they have. For consistency, I list them in the order the appear on the diagram shown above, following standard left to right and top to bottom ordering: Side 1, Side 2, Side 4, etc.
#
# I label each dot with one index corresponding to the side it's on and a second index corresponding to its ordering on the side, left to right and top to bottom from the perspective of the diagram above. And again, I list them in order: Dot 1A, Dot 2A, Dot 2B, Dot 4A, Dot 4B, Dot 4C, Dot 4D, etc. For clarity:
#
# Labeling of dots:
#  ___________ ___________
# |           |        2A |
# |     1A    |        +  |
# |     +     |  2B       |
# |           |  +        |
# |___________|___________|___________
#             |  4A    4B |  6A 6B 6C |
#             |  +     +  |  +  +  +  |
#             |  4C    4D |           |
#             |  +     +  |  +  +  +  |
#             |___________|___________|___________
#                         |  5A    5B |        3A |
#                         |  +  5C +  |     3B +  |
#                         |  5D +  5E |  3C +     |
#                         |  +     +  |  +        |
#                         |___________|___________|
#
# It's also helpful to have coordinates for each side ("local" coordinates) and for the whole cube ("global" coordinates). Local coordinates always start at [0,0], but global coordinates start at different values for different sides so they're always unique. (See elsewhere below for more details.)
#
# Similarly, these local and global coordinates correspond to squares the local and global Hilbert curves pass through. Locally (on each side), the squares are indexed starting at 0. Globally, the squares are indexed starting at 0 on Side 1 and with higher indices across the cube.
#####

import math
from cmath import sqrt
from turtle import st
# for simple data tables
from tabulate import tabulate
# for colors in tables
from colorama import init, Back, Fore
# for Hilbert curve calculations
from hilbertcurve.hilbertcurve import HilbertCurve
# for Hilbert curve diagrams
import matplotlib.pyplot as plt


#####
### Preliminaries
#
# - 'square' indexes the squares the Hilbert curve runs through, starting at 0.
# - 'domino' indexes the dominoes tiling the cube, starting at 1. Note that here the tiling is defined to begin with a half domino. It could begin with a full one — both are valid tilings in line with the Hilbert curve, so it's a matter of choice. (I explain elsewhere why choosing the half domino approach was important for this project.)
# - 'term' indexes the term in my Wallis-like domino train product. Each such term includes two fractions multiplied together.
# - 'number' is the number of a particular half domino on a square.
##### 

#####
#### `get_domino()`
#
# > Given `square`, find which domino tiles it.
#####
# The 2 here is not a variable because only real, standard dominoes (which cover two squares) are considered.
def get_domino(square):
    return math.trunc(math.floor((square + 1) / 2)) + 1

#####
#### `get_term()`
#
# > Given `square`, find which term it corresponds to.
#####
# The 4 here is not a variable because my Wallis-like domino train product always has 2 fractions with 4 numerator/denominator values.
def get_term(square):
    return math.trunc(math.floor((square + 1) / 4)) + 1

#####
#### `get_number()`
#
# > Given `square`, find which domino number covers it.
#####
# The 7s here are not variables because only real, standard half dominoes (which have 7 possible values, from 0 to 6) are considered. The 4 here is not a variable because my Wallis-like domino train product always has 2 fractions with 4 numerator/denominator values.
def get_number(square):
    term = get_term(square)
    # numerators and denominators of first and second fractions
    num_1 = (2 * term - 1) % 7
    den_1 = num_2 = (2 * term + math.trunc(math.floor((term - 1) / 7))) % 7
    den_2 = (2 * term + 1) % 7
    # a condition to pick out which numerator or denominator to set the number of a square to
    # the number is just one of the two sections of a domino
    condition = (square + 1) % 4
    if condition == 0: return num_1
    elif condition == 1: return den_1
    elif condition == 2: return num_2
    else: return den_2



#####
### Hilbert Curve Parameters
#
# **Important:** Note that these are paramaters for the local Hilbert curves on one side of the cube, not the global Hilbert curve covering the whole cube.
#
# This uses the [`hilbertcurve`](https://pypi.org/project/hilbertcurve/) package.
#
# - `iterations` is the number of iterations of (the polygonal approximation to) the Hilbert curve. For _Metaphysics_, this will be 4 for the smallest scale version but greater for the larger scale versions.
# - `dimensions` is the number of spatial dimensions. For _Metaphysics_, this will always be 2, since each local Hilbert curve corresponds to a tiling of one side of a cube (which has 2 dimensions).
# - `num_squares` is the number of squares in a Hilbert curve with so many iterations and of so many dimensions. (Since for _Metaphysics_ I'm always using 2 dimensions, I use the more specific term "squares" rather than the fully general "hypercubes".) In general, a Hilbert curve fills a hypercube with $2^{i \cdot d}$ unit hypercubes contained within in it, where $i$ is `iterations` and $d$ is `dimensions`. For 4 iterations and 2 dimensions, that's a square with $2^8 = 256$ unit squares contained within it.
#
# For _Metaphysics_, I'm using 6 connecting Hilbert curves (in 2 different orientations) to cover the surface of a cube.

# Note that this package produces a Hilbert curve that begins at the lower left and ends at the lower right corner. As a result, no matter the orientation of the Hilbert curve considered here, I pick coordinates such that [0,0] is at the beginning and [`sqrt(num_squares)`,0] is at the end.
#####
iterations = 4
dimensions = 2
hilbert_curve = HilbertCurve(iterations, dimensions)
num_coordinates_per_side = 2 ** iterations
num_squares = num_coordinates_per_side ** dimensions


#####
### Dots
# 
# "Dots" are the dots (sometimes called "pips") on a die.
#
# For _Metaphysics_, I'm using two Hilbert curve orientations, which means there are two types of tiling. These have different coordinate orientations, as described above.
#
# The dots have two indices. The first (1, 2, 3, ...) indicates the side of the die the dot is on. The second (A, B, C, ...) indicates the order of the dot on its side of the die. They're ordered from left to right and top to bottom from the perspective of the diagram above. 
#
# These are currently defined only for a Hilbert curve tiling with 256 squares. Ideally, they'd be defined independently of the number of squares, with variables rather than numbers, but doing this is complicated because the shape of each dot and the spacing between dots and the cube sides should change with the number of squares. So, I'm skipping this for now.
#
# **Important:** The dots are first defined in "local" coordinates, where each side's coordinates goes from [0,0] to [`sqrt(num_squares)` - 1, `sqrt(num_squares)` - 1], i.e. [15,15]. (Confusingly enough, these are "global variables" in the programming sense!) The `get_global_coordinates()` function further below will later transform these local coordinates into global ones, where each side's coordinates start at a multiple of `sqrt(num_squares)` times the index of the side in the ordering shown in the diagram above (starting at 0). In global coordinates, the first side starts at [0,0], the second at [`sqrt(num_squares)`,`sqrt(num_squares)`] i.e. [16,16], the third at [`2 * sqrt(num_squares)`, `2 * sqrt(num_squares)`] i.e. [32,32], and so on.
#####
dot_1A = [[6,7], [6,8], [7,6], [7,7], [7,8], [7,9], [8,6], [8,7], [8,8], [8,9], [9,7], [9,8]]
dot_2A = [[1,12], [1,13], [2,11], [2,12], [2,13], [2,14], [3,11], [3,12], [3,13], [3,14], [4,12], [4,13]]
dot_2B = [[11,2], [11,3], [12,1], [12,2], [12,3], [12,4], [13,1], [13,2], [13,3], [13,4], [14,2], [14,3]]
dot_4A = [[1,2], [1,3], [2,1], [2,2], [2,3], [2,4], [3,1], [3,2], [3,3], [3,4], [4,2], [4,3]]
dot_4B = [[11,2], [11,3], [12,1], [12,2], [12,3], [12,4], [13,1], [13,2], [13,3], [13,4], [14,2], [14,3]]
dot_4C = [[1,12], [1,13], [2,11], [2,12], [2,13], [2,14], [3,11], [3,12], [3,13], [3,14], [4,12], [4,13]]
dot_4D = [[11,12], [11,13], [12,11], [12,12], [12,13], [12,14], [13,11], [13,12], [13,13], [13,14], [14,12], [14,13]]
dot_6A = [[1,2], [1,3], [2,1], [2,2], [2,3], [2,4], [3,1], [3,2], [3,3], [3,4], [4,2], [4,3]]
dot_6B = [[1,7], [1,8], [2,6], [2,7], [2,8], [2,9], [3,6], [3,7], [3,8], [3,9], [4,7], [4,8]]
dot_6C = [[1,12], [1,13], [2,11], [2,12], [2,13], [2,14], [3,11], [3,12], [3,13], [3,14], [4,12], [4,13]]
dot_6D = [[11,2], [11,3], [12,1], [12,2], [12,3], [12,4], [13,1], [13,2], [13,3], [13,4], [14,2], [14,3]]
dot_6E = [[11,7], [11,8], [12,6], [12,7], [12,8], [12,9], [13,6], [13,7], [13,8], [13,9], [14,7], [14,8]]
dot_6F = [[11,12], [11,13], [12,11], [12,12], [12,13], [12,14], [13,11], [13,12], [13,13], [13,14], [14,12], [14,13]]
dot_5A = [[1,2], [1,3], [2,1], [2,2], [2,3], [2,4], [3,1], [3,2], [3,3], [3,4], [4,2], [4,3]]
dot_5B = [[11,2], [11,3], [12,1], [12,2], [12,3], [12,4], [13,1], [13,2], [13,3], [13,4], [14,2], [14,3]]
dot_5C = [[6,7], [6,8], [7,6], [7,7], [7,8], [7,9], [8,6], [8,7], [8,8], [8,9], [9,7], [9,8]]
dot_5D = [[1,12], [1,13], [2,11], [2,12], [2,13], [2,14], [3,11], [3,12], [3,13], [3,14], [4,12], [4,13]]
dot_5E = [[11,12], [11,13], [12,11], [12,12], [12,13], [12,14], [13,11], [13,12], [13,13], [13,14], [14,12], [14,13]]
dot_3A = [[1,12], [1,13], [2,11], [2,12], [2,13], [2,14], [3,11], [3,12], [3,13], [3,14], [4,12], [4,13]]
dot_3B = [[6,7], [6,8], [7,6], [7,7], [7,8], [7,9], [8,6], [8,7], [8,8], [8,9], [9,7], [9,8]]
dot_3C = [[11,2], [11,3], [12,1], [12,2], [12,3], [12,4], [13,1], [13,2], [13,3], [13,4], [14,2], [14,3]]

dots = [dot_1A, dot_2A, dot_2B, dot_4A, dot_4B, dot_4C, dot_4D, dot_6A, dot_6B, dot_6C, dot_6D, dot_6E, dot_6F, dot_5A, dot_5B, dot_5C, dot_5D, dot_5E, dot_3A, dot_3B, dot_3C]
dot_names = ['Dot 1A', 'Dot 2A', 'Dot 2B', 'Dot 4A', 'Dot 4B', 'Dot 4C', 'Dot 4D', 'Dot 6A', 'Dot 6B', 'Dot 6C', 'Dot 6D', 'Dot 6E', 'Dot 6F', 'Dot 5A', 'Dot 5B', 'Dot 5C', 'Dot 5D', 'Dot 5E', 'Dot 3A', 'Dot 3B', 'Dot 3C']


#####
### Sides
#
# "Sides" are the sides of a dice. The variables below list the dots on each side. This is somewhat redundant, since the dot names encode the side they're on (with their first index), but it's useful to have this information consolidated.
#
# **Important:** `generic_side` here is a generic side in local coordinates. This makes it possible to calculate white areas in local coordinates using dots in local coordinates together with this generic side.
#####

side_1_dots = [dot_1A]
side_2_dots = [dot_2A, dot_2B]
side_4_dots = [dot_4A, dot_4B, dot_4C, dot_4D]
side_6_dots = [dot_6A, dot_6B, dot_6C, dot_6D, dot_6E, dot_6F]
side_5_dots = [dot_5A, dot_5B, dot_5C, dot_5D, dot_5E]
side_3_dots = [dot_3A, dot_3B, dot_3C]

sides_dots = [side_1_dots, side_2_dots, side_4_dots, side_6_dots, side_5_dots, side_3_dots]
num_sides_dots = len(sides_dots)
side_names = ['Side 1', 'Side 2', 'Side 4', 'Side 6', 'Side 5', 'Side 3']

# This is defined so that it can be passed along with dots (defined above), since the order and length of the two lists match.
dots_side_list = [side_1_dots, side_2_dots, side_2_dots, side_4_dots, side_4_dots, side_4_dots, side_4_dots, side_6_dots, side_6_dots, side_6_dots, side_6_dots, side_6_dots, side_6_dots, side_5_dots, side_5_dots, side_5_dots, side_5_dots, side_5_dots, side_3_dots, side_3_dots, side_3_dots]

#####
#### `get_generic_side()`
#
# > Calculate a generic side, i.e. one in local coordinates. This includes [0,0], [0,1], [0,2], ..., [1,0], [1,1], [1,2], ..., [sqrt(num_squares),sqrt(num_squares)].
#####
def get_generic_side():
    generic_side = []
    for i in range(0, int(math.sqrt(num_squares))):
        for j in range(0, int(math.sqrt(num_squares))):
            generic_side.append([i,j])
    return generic_side

generic_side = get_generic_side()


#####
### White Areas
#
# A 'white area' is the part of a side that isn't the dots.
#
# Important: The white areas are first defined in local coordinates. The get_global_coordinates() function below will later transform them into global coordinates.
#####

#####
#### `get_white_area()`
#
# > Given a particular side (which is a list of dots), find the list of coordinates for its white area. This is done by removing the coordinates for the dots on the given side.
#####
def get_white_area(side_dots):
    # I find this quite nonintuitive, but white_area = generic_side doesn't work here because that syntax just creates a reference to the original list rather than creating a copy of that list. So, it's necessary to explictly copy the list so we can make changes to the new list values without changing the corresponding original list values. There are many ways to do this: see https://stackoverflow.com/questions/2612802/list-changes-unexpectedly-after-assignment-why-is-this-and-how-can-i-prevent-it. 
    # For some reason that I have been unable to figure out, `white_area = generic_side[:]`, `white_area = generic_side.copy()`, and the like do NOT work. But `white_area = get_generic_side()` does!
    # white_area = generic_side[:]
    white_area = get_generic_side()
    # remove the coordinates corresponding to dots
    for k in range(0, len(side_dots)):
        for l in range(0, len(side_dots[k])):
            # side[k] is a dot, and side[k][l] is a coordinate in that dot
            # list order matters, so this will remove e,g, [2,5] but not [5,2]
            white_area.remove(side_dots[k][l])
    return white_area

white_area_1 = get_white_area(side_1_dots)
white_area_2 = get_white_area(side_2_dots)
white_area_4 = get_white_area(side_4_dots)
white_area_6 = get_white_area(side_6_dots)
white_area_5 = get_white_area(side_5_dots)
white_area_3 = get_white_area(side_3_dots)

white_areas = [white_area_1, white_area_2, white_area_4, white_area_6, white_area_5, white_area_3]
white_area_names = ['White Area for Side 1', 'White Area for Side 2', 'White Area for Side 4', 'White Area for Side 6', 'White Area for Side 5', 'White Area for Side 3']

# This is defined so that it can be passed along with white areas, since the order and length of the two lists match.
white_areas_side_list = [side_1_dots, side_2_dots, side_4_dots, side_6_dots, side_5_dots, side_3_dots]


#####
### Diagrams
#####

#####
#### `create_hilbert_curve_diagram()`
#
# > Create a Hilbert curve diagram.
# >
# > This adapts code from the GitHub repo of the [`hilbertcurve`](https://pypi.org/project/hilbertcurve/) package. The side index is that of the ordering of sides defined above. This function creates a diagram for one side at a time.
# >
# > Note that, currently, this does not adjust the orientation of the Hilbert curve to be type 1 or 2 for a given side (as defined above). All Hilbert curves it produces are in "standard" orientation.
#####
def create_hilbert_curve_diagram(side_index):
    # this has to be at the beginning, not with the other 'plt' statements below
    plt.figure(figsize = (10,10))
    min_coordinate = 0
    max_coordinate = num_coordinates_per_side - 1
    cmin = min_coordinate - 0.5
    cmax = max_coordinate + 0.5
    colors = ['red', 'blue', 'black', 'green', 'purple', 'cyan', 'gray']
    line_widths = [32, 16, 8, 4, 2, 1, 0.5]
    offset = 0
    dx = 0.5
    for i in range(iterations, iterations - 1, -1):
        curve = HilbertCurve(i, dimensions)
        num_coordinates_per_side_i = 2 ** i
        num_points = 2 ** (i * dimensions)
        points = []
        for j in range(num_points):
            points.append(curve.point_from_distance(j))
        points = [
            [(point[0] * num_coordinates_per_side / num_coordinates_per_side_i) + offset,
            (point[1] * num_coordinates_per_side / num_coordinates_per_side_i) + offset]
            for point in points]
        connectors = range(3, num_points, 4)
        color = colors[i - 1]
        # '+ len(line_widths) - iterations' so it starts at a smaller line width (later in the list) when iterations is smaller than the number of line width values
        # Note that to increase iterations beyond this number, more line width values (and colors) should be added
        line_width = line_widths[i - 1 + len(line_widths) - iterations]
        for k in range(num_points - 1):
            if k in connectors:
                line_style = '--'
                alpha = 0.5
            else:
                line_style = '-'
                alpha = 1.0
            plt.plot((points[k][0], points[k + 1][0]), (points[k][1], points[k + 1][1]),
                    color = color, linewidth = line_width, linestyle = line_style, alpha = alpha)
        for l in range(num_points):
            plt.scatter(points[l][0], points[l][1], 60, color = color)
            plt.text(points[l][0] + 0.1, points[l][1] + 0.1, str(l + side_index * num_points), color = color)
        offset += dx
        dx *= 2
    plt.title('Hilbert Curve Pattern for ' + str(side_names[side_index]))
    plt.grid(alpha = 0.3)
    plt.xlim(cmin, cmax)
    plt.ylim(cmin, cmax)
    plt.xlabel('x', fontsize = 16)
    plt.ylabel('y', fontsize = 16)
    plt.tight_layout()
    plt.savefig(str(side_names[side_index]) + ' - ' + str(iterations) + ' iterations, ' + str(dimensions) + ' dimensions.png')

create_hilbert_curve_diagram(0)


# Colors for tables
backs = [Back.LIGHTBLUE_EX, Back.WHITE, Back.GREEN, Back.YELLOW, Back.LIGHTMAGENTA_EX, Back.CYAN, Back.LIGHTRED_EX]
num_colors = len(backs)


#####
### Values and Number Counts
#####

#####
#### `print_values()`
#
# > Print a table of values for square, domino, term, and number.
#####
def print_values():
    column_headers = ['Square', 'Domino', 'Term', 'Number']
    # array of rows
    data = []
    # num_squares is per side, but we want to tile all sides of the cube
    for square in range(0, num_squares * num_sides_dots):
        domino = get_domino(square)
        term = get_term(square)
        number = get_number(square)
        # setting different colors for different numbers (mainly because it's very easy to confuse 0 and 6 when reading the table)
        # see https://compucademy.net/python-tables-for-multiplication-and-addition/
        color = backs[number % num_colors]
        data.append([square, domino, term, f'{color}{number}{Back.RESET}'])
        # add dash to simulate 'fraction' line on domino
        if square % 2 == 1: data.append(['', '', '', '-'])
        # add spacing to make it easier to see pairs of domino numbers
        else: data.append(['........', '........', '......', '........'])
    print(tabulate(data, column_headers, tablefmt = "pretty"))

print_values()


#####
#### `print_number_counts()`
#
# > Print counts of how many times each number (0 through 6) appears.
#####
def print_number_counts():
    column_headers = ['0', '1', '2', '3', '4', '5', '6']
    # array of rows
    data = [[0, 0, 0, 0, 0, 0, 0]]
    for square in range(0, num_squares * num_sides_dots):
        number = get_number(square)
        data[0][number] += 1
    print("Number counts:")
    print(tabulate(data, column_headers))

print_number_counts()


#####
### Coordinates and Squares
#
# In local coordinates, each side's coordinates goes from [0,0] to [`sqrt(num_squares) - 1`, `sqrt(num_squares) - 1`]. 
#
# In global coordinates, each side's coordinates start at a multiple of `sqrt(num_squares)` times the index of the side in the ordering shown in the diagram above (starting at 0). So, the first side starts at [0,0], the second at [`sqrt(num_squares)`,`sqrt(num_squares)`], the third at [`2 * sqrt(num_squares)`, `2 * sqrt(num_squares)`], and so on.
#
# Note that a single coordinate is a list of numbers (with a number of elements equal to `dimensions`, which for _Metaphysics_ is always 2), e.g. [2,5]. Coordinates (plural) are lists of such lists, e.g. [[2,5], [6,3], [7,7]]. 
#####

#####
#### `get_global_coordinates()`
#
# > Given local coordinates and a side, get the corresponding global coordinates.
# >
# > This function requires a list of lists input, even for a single coordinate.
#####
def get_global_coordinates(local_coordinates, side_dots):
    global_coordinates = []
    for i in range(0, len(local_coordinates)):
        global_coordinate = []
        for j in range(0, len(local_coordinates[i])):
            global_coordinate.append(local_coordinates[i][j] + sides_dots.index(side_dots) * int(math.sqrt(num_squares)))
        global_coordinates.append(global_coordinate)
    return global_coordinates

#####
#### `set_global_coordinates()`
#
# > Given local coordinates and a side, set the corresponding global coordinates.
# >
# > This function requires a list of lists input, even for a single coordinate.
#####
def set_global_coordinates(local_coordinates, side_dots):
    for i in range(0, len(local_coordinates)):
        for j in range(0, len(local_coordinates[i])):
            local_coordinates[i][j] += sides_dots.index(side_dots) * int(math.sqrt(num_squares))

#####
#### `set_global_coordinates_batch()`
#
# > Given local coordinates and a side, set the corresponding global coordinates in a batch.
# >
# > This function requires a list of a list of lists input.
# > Note that the order and length of the two lists (`local_coordinates_list` and `side_list``) **must** match so that each local_coordinates matches the appropriate side.
#####
def set_global_coordinates_batch(local_coordinates_list, side_list):
    for i in range(0, len(local_coordinates_list)): set_global_coordinates(local_coordinates_list[i], side_list[i])

# Transform dots and white areas from local into global coordinates.
set_global_coordinates_batch(dots, dots_side_list)
set_global_coordinates_batch(white_areas, white_areas_side_list)


#####
#### `get_squares()`
#
# > Given (a list of) global coordinates (e.g. a dot), find the squares (ordered along the Hilbert curve) that the list includes.
# >
# > Note that either local or global coordinates can be inputted, but the output will always be global square numbers.
# >
# > The input list of coordinates is in the number of the dimensions of the Hilbert curve (always 2 for _Metaphysics_).
# >
# > The output is an (ordered) list of coordinates in 1 dimension, since the Hilbert curve itself is 1-dimensional (at least "stretched out", since the "curled up" curve has fractal Hausdorff dimension 2).
#####
def get_squares(coordinates):
    # Calculate the side index as a kind of offset: how many times the coordinate values can be divided by sqrt(num_squares). (We can used any coordinate value to find this — coordinates[0][0] is just an arbitrary choice.) For example, if the coordinate value is 18 and sqrt(num_squares) is 16, the offset is 1 because 18 can be divided by 16 once. This ia also the side index of that coordinate: it's on the second side. 
    # This index could instead be passed into the function, but it's helpful to calcuate it here so that's not necessary.
    # Note that this should always be an integer: math.trunc() and math.floor() are just safeguards.
    side_index = int((coordinates[0][0] - (coordinates[0][0] % int(math.sqrt(num_squares)))) / int(math.sqrt(num_squares)))
    local_coordinates = []
    for coordinate in coordinates:
        local_coordinate = []
        for i in range(0, len(coordinate)):
            # Mod by sqrt(num_squares) to make the coordinate local, so that distances_from_points from the hilbertcurve package can be used to calculate local square numbers.
            local_coordinate.append(coordinate[i] % int(math.sqrt(num_squares)))
        local_coordinates.append(local_coordinate)
    points = local_coordinates
    distances = hilbert_curve.distances_from_points(points)
    # Finally, calculate global square values simply by adding num_squares (per side), scaled by the side index
    global_squares = []
    for distance in distances: global_squares.append(distance + (side_index * num_squares))
    return global_squares


#####
#### `print_squares()`
#
# > Print squares for a given list of coordinates.
# >
# > The relevant group of coordinates and their names must be passed also.
#####
def print_squares(coordinates, coordinates_group, coordinate_names):
    print('Squares for ' + coordinate_names[coordinates_group.index(coordinates)] + ':')
    print(get_squares(coordinates))

print_squares(dot_1A, dots, dot_names)


#####
#### `get_other_domino_square()`
#
# > Given the domino number of a square, find the other square with that domino number.
# >
# > There's only one, and it's either the previous or next one.
#####
def get_other_domino_square(square):
    previous_square = square - 1
    next_square = square + 1
    if get_domino(square) == get_domino(previous_square): return previous_square
    else: return next_square

#####
### Dominoes and Domino Counts
#####

#####
#### `get_dominoes()`
#
# > Given (a list of) coordinates (e.g. a dot), find the full and half dominoes that compose it.
#####
def get_dominoes(coordinates):
    full_dominoes = []
    half_dominoes = []
    coordinates_squares = get_squares(coordinates)
    # the squares already 'used', or included in a full or half domino already added
    used_squares = []
    for square in coordinates_squares: 
        other_domino_square = get_other_domino_square(square)
        # if the other domino square isn't used
        if (other_domino_square not in used_squares):
            # if it's in the coordinates (e.g. a dot), add to full dominoes
            if (other_domino_square in coordinates_squares): 
                domino = [get_number(square), get_number(other_domino_square)]
                # sort to avoid counting e.g. [2,5] and [5,2] separately — they should be treated as the same
                domino.sort()
                full_dominoes.append(domino)
                # add squares to used list
                # not strictly necessary to add 'square', since we're iterating over it (i.e. the for loop takes care of not considering it multiple times), but it's more intuitive to also consider it 'used'
                used_squares.extend([square, other_domino_square])
            # else, add to half dominoes
            else:
                half_domino = get_number(square)
                half_dominoes.append(half_domino)
                used_squares.append(square)
    return full_dominoes, half_dominoes


#####
#### `get_dominoes_counts()`
#
# > Given (lists of) full and half dominoes (e.g. for a single dot), count how many there are of each type. 
# >
# > Order doesn't matter for full dominoes, e.g. [2,5] and [5,2] are considered the same. This will be used in table data, so notice that there are row headers included (which aren't themselves counts, of course).
# >
# > The table data for half dominoes has only one row, so there's no need for a row header there.
#####
def get_dominoes_counts(full_dominoes, half_dominoes):
    # the first values are row headers
    full_dominoes_counts = [[0], [1], [2], [3], [4], [5], [6]]
    # no need for row headers — there's only one row
    half_dominoes_counts = [[]]
    for i in range(0, 7):
        for j in range (0, i + 1):
            # They're already sorted in get_dominoes(), so no need to count both [i,j] and [j,i].
            # If it seems odd that it's [j,i] below, that's only because j is never greater than i given this iteration strategy, so it should come first because sort(), used in get_dominoes(), puts smaller numbers first (i.e. ascending order).
            full_dominoes_count = full_dominoes.count([j,i])
            # Using append() here takes care of the ordering, so no need to use the j index.
            full_dominoes_counts[i].append(full_dominoes_count)
        half_dominoes_count = half_dominoes.count(i)
        half_dominoes_counts[0].append(half_dominoes_count)
    return full_dominoes_counts, half_dominoes_counts


#####
#### `get_sum_dominoes_counts()`
#
# > Given (a list of a list of) coordinates (e.g. a list of dots), find the sum of counts for full and half dominoes.
#####
def get_sum_dominoes_counts(coordinates):
    # initilialize with zero values so they can later be overwritten (to avoid 'index out of range' error)
    sum_full_dominoes_counts = [
        [0, 0], 
        [1, 0, 0], 
        [2, 0, 0, 0], 
        [3, 0, 0, 0, 0], 
        [4, 0, 0, 0, 0, 0], 
        [5, 0, 0, 0, 0, 0, 0], 
        [6, 0, 0, 0, 0, 0, 0, 0]]
    sum_half_dominoes_counts = [[0, 0, 0, 0, 0, 0, 0]]
    for i in range(0, len(coordinates)):
        full_dominoes, half_dominoes = get_dominoes(coordinates[i])
        full_dominoes_counts, half_dominoes_counts = get_dominoes_counts(full_dominoes, half_dominoes)
        # add up full domino counts
        for j in range(0, len(full_dominoes_counts)):
            # start at 1 since the first items are just row headers
            for k in range(1, len(full_dominoes_counts[j])):
                # adjust by frequency
                sum_full_dominoes_counts[j][k] += (full_dominoes_counts[j][k])
        # add up half domino counts
        for l in range(0, len(half_dominoes_counts[0])):
            # adjust by frequency
            sum_half_dominoes_counts[0][l] += (half_dominoes_counts[0][l])
    return sum_full_dominoes_counts, sum_half_dominoes_counts


#####
# #### `print_dominoes_counts()`
#
# > Given (a list of a list of) coordinates (e.g. a list of dots), print tables of full and half domino counts.
# >
# > The "names" input is a list of names for each list of coordinates.
#####
def print_dominoes_counts(coordinates, names):
    full_dominoes_headers = ['#', 0, 1, 2, 3, 4, 5, 6]
    half_dominoes_headers = [0, 1, 2, 3, 4, 5, 6]
    for i in range(0, len(coordinates)):
        full_dominoes, half_dominoes = get_dominoes(coordinates[i])
        full_dominoes_counts, half_dominoes_counts = get_dominoes_counts(full_dominoes, half_dominoes)
        print('Full Dominoes for ' + names[i] + ':')
        print(tabulate(full_dominoes_counts, full_dominoes_headers))
        print('Half Dominoes for ' + names[i] + ':')
        print(tabulate(half_dominoes_counts, half_dominoes_headers))
    sum_full_dominoes_counts, sum_half_dominoes_counts = get_sum_dominoes_counts(coordinates)
    print('Full Dominoes for All:')
    print(tabulate(sum_full_dominoes_counts, full_dominoes_headers))
    print('Half Dominoes for All:')
    print(tabulate(sum_half_dominoes_counts, half_dominoes_headers))

print_dominoes_counts(dots, dot_names)
print_dominoes_counts(white_areas, white_area_names)


#####
#### `get_min_num_sets()`
#
# > Given counts of full and half dominoes, find the minimum number of domino sets required.
# >
# > A standard domino set (with column and row headers) is:
# ```python
# #   #    0    1    2    3    4    5    6
# # ---  ---  ---  ---  ---  ---  ---  ---
# #   0    1
# #   1    1    1
# #   2    1    1    1
# #   3    1    1    1    1
# #   4    1    1    1    1    1
# #   5    1    1    1    1    1    1
# #   6    1    1    1    1    1    1    1
# ```
# >
# > That is, it has one domino of each type. As a result, there are 8 half dominoes of each number (0 through 6).
#####
def get_min_num_sets(full_dominoes_counts, half_dominoes_counts):
    # The minimum number of sets must be at least as great as the highest full dominoes count. (That's because there's no other way to get a particular full domino than through a new set, since each set has only one of a given type.)
    max_full_dominoes_count = 0
    for i in range(0, len(full_dominoes_counts)):
        # start at 1 since the first items are just row headers
        for j in range(1, len(full_dominoes_counts[i])):
            if full_dominoes_counts[i][j] > max_full_dominoes_count: max_full_dominoes_count = full_dominoes_counts[i][j]
    min_num_sets = max_full_dominoes_count
    # how many half dominoes are left over for use
    # initilialize with zero values so they can later be overwritten (to avoid 'index out of range' error)
    leftover_half_dominoes_counts = [[0, 0, 0, 0, 0, 0, 0]]
    # Loop through again and set each leftover full dominoes count to be the difference between the (provisional) minimum number of sets and the value of the corresponding full dominoes count.
    for k in range(0, len(full_dominoes_counts)):
        # start at 1 since the first items are just row headers
        for l in range(1, len(full_dominoes_counts[k])):
            # k and l - 1 (minus 1 because there aren't row headers for the leftover half dominoes counts list) are the domino numbers, so add to those leftover half dominoes counts
            leftover_half_dominoes_counts[0][k] += max_full_dominoes_count - full_dominoes_counts[k][l]
            leftover_half_dominoes_counts[0][l - 1] += max_full_dominoes_count - full_dominoes_counts[k][l]
    # Check if there are enough leftover half dominoes.
    for m in range(0, len(half_dominoes_counts[0])):
        while half_dominoes_counts[0][m] > leftover_half_dominoes_counts[0][m]:
            # If there aren't enough leftover half dominoes with a particular number, we don't have enough sets. So, increment the minimum number of sets by 1 and the leftover half dominoes counts by 8 (since each set has 8 half dominoes of a particular number)
            min_num_sets += 1
            for n in range(0, len(leftover_half_dominoes_counts[0])): leftover_half_dominoes_counts[0][n] += 8
    # Once we have enough leftover half dominoes for each number, we have the minimum number of sets.
    return min_num_sets
    

#####
#### `print_min_num_sets()`
#
# > Given (a list of a list of) coordinates (e.g. a list of dots), print the minimum number of sets to cover them.
#####
def print_min_num_sets(coordinates):
    sum_full_dominoes_counts, sum_half_dominoes_counts = get_sum_dominoes_counts(coordinates)
    min_num_sets = get_min_num_sets(sum_full_dominoes_counts, sum_half_dominoes_counts)
    print('Minimum Number of Domino Sets to Cover All:')
    print(min_num_sets)

print_min_num_sets(dots)
print_min_num_sets(white_areas)






#### NOT WORKING BELOW ###

# side_1_dots_coordinates = dot_1A + white_area_1
# side_2_dots_coordinates = dot_2A + dot_2B + white_area_2
# side_4_dots_coordinates = dot_4A + dot_4B + dot_4C + dot_4D + white_area_4
# side_6_dots_coordinates = dot_6A + dot_6B + dot_6C + dot_6D + dot_6E + dot_6F + white_area_6
# side_5_dots_coordinates = dot_5A + dot_5B + dot_5C + dot_5D + dot_5E + white_area_5
# side_3_dots_coordinates = dot_3A + dot_3B + dot_3C + white_area_3


#####
# Given a number of domino sets and full and half domino counts, return a list of dominoes to cut. Note that the number of dominoes in this list is the number of cuts to make.
#
# This list is optimal in the sense that this function first includes dominoes for which both halves can be used and the dominoes for which only one half can be used. There may be even more optimal sets, but I imagine determining this would involved fairly detailed combinatorics of the dominoes themselves and I haven't bothered to think about it much!
#
# Note that the number of sets must be greater than or equal to the maximum count in the full dominoes counts. Otherwise, the function will raise an error.
#####
def get_optimal_cuts(num_sets, full_dominoes_counts, half_dominoes_counts):
    leftover_full_dominoes_counts = full_dominoes_counts
    for i in range(0, len(full_dominoes_counts)):
        # start at 1 since the first items are just row headers
        for j in range(1, len(full_dominoes_counts[i])):
            if num_sets - full_dominoes_counts[i][j] < 0: raise ValueError('Not enough domino sets for this many dominoes!')
            leftover_full_dominoes_counts[i][j] = num_sets - full_dominoes_counts[i][j]
    dominoes_to_cut = []
    # First, add to the list all dominoes for which both halves can be used.
    for k in range(0, len(leftover_full_dominoes_counts)):
        # start at 1 since the first items are just row headers
        for l in range(1, len(leftover_full_dominoes_counts[k])):
            # If both halves of the full domino can be used, add this domino to this list of those to be cut.
            # k and l - 1 (minus 1 because there aren't row headers for the leftover half dominoes counts list) are the domino numbers, so add to those leftover half dominoes counts
            if half_dominoes_counts[0][k] != 0 and half_dominoes_counts[0][l - 1] !=0: 
                dominoes_to_cut.append([k,l])
                leftover_full_dominoes_counts[k][l] -= 1
    # Then, loop through again and add to the list all dominoes for which only one half can be used.
    for m in range(0, len(leftover_full_dominoes_counts)):
        # start at 1 since the first items are just row headers
        for l in range(1, len(leftover_full_dominoes_counts[k])):
            if half_dominoes_counts[0][k] != 0 or half_dominoes_counts[0][l - 1] !=0:
                dominoes_to_cut.append([k,l])
                leftover_full_dominoes_counts[k][l] -= 1
    return dominoes_to_cut

def print_optimal_cuts(coordinates):
    sum_full_dominoes_counts, sum_half_dominoes_counts = get_sum_dominoes_counts(coordinates)
    min_num_sets = get_min_num_sets(sum_full_dominoes_counts, sum_half_dominoes_counts)
    dominoes_to_cut = get_optimal_cuts(min_num_sets, sum_full_dominoes_counts, sum_half_dominoes_counts)
    print("Optimal dominoes to cut:")
    print(dominoes_to_cut)

# print("For dots:")
# print_optimal_cuts(dots)
# print("For white areas:")
# print_optimal_cuts(white_areas)