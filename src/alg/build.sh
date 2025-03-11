ALG_SOURCE_FILES="grid/grid.cpp grid/grid_refinement.cpp grid/grid_derefinement.cpp cell/cell.cpp cell/cell_derefinement.cpp cell/cell_refinement.cpp  grid_purkinje/grid_purkinje.cpp"
ALG_HEADER_FILES="grid/grid.h cell/cell.h grid_purkinje/grid_purkinje.h"

COMPILE_STATIC_LIB "alg" "$ALG_SOURCE_FILES" "$ALG_HEADER_FILES"
