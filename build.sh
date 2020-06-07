#!/bin/bash
OPTIONS_FILE=./bsbash/build_functions.sh
FUNCTIONS_FILE=./bsbash/find_functions.sh

if [ -f "$OPTIONS_FILE" ]; then
    # shellcheck disable=SC1090
    source $OPTIONS_FILE
else
    echo "$OPTIONS_FILE not found, aborting compilation"
    exit 1
fi

if [ -f "$FUNCTIONS_FILE" ]; then
	# shellcheck disable=SC1090
	source $FUNCTIONS_FILE
fi

###########User code#####################
COMPILE_GUI=''
COMPILE_MPI=''
COMPILE_CONVERTER=''
COMPILE_SIMULATOR=''
COMPILE_WITH_DDM=''

GET_BUILD_OPTIONS "$@"

if [ "$BUILD_TYPE" == "release" ]; then
    C_FLAGS="$C_FLAGS -O3"
elif [ "$BUILD_TYPE" == "debug" ]; then
    C_FLAGS="$C_FLAGS -g -DDEBUG_INFO"
else
  	PRINT_ERROR "$BUILD_TYPE is not a valid BUILD_TYPE."
  	PRINT_ERROR "Valid BUILD_TYPE options are: release, debug (-r or -d options)"
  	exit 1
fi

for i in "${BUILD_ARGS[@]}"; do

    if [ "$i" == "clean" ]; then
        echo "Cleaning $BUILD_TYPE"
        CLEAN_PROJECT "$BUILD_TYPE"
        cd src/3dparty/raylib/src || exit 1;
        make clean
        cd - || exit 1;
        exit 0    
    fi

    if [ "$i" == "ddm" ]; then
        COMPILE_WITH_DDM='y'
        C_FLAGS="$C_FLAGS -DENABLE_DDM"
    fi

    if [ "$i" == "all" ]; then
        COMPILE_GUI='y'
        COMPILE_MPI='y'
        COMPILE_SIMULATOR='y'
        COMPILE_CONVERTER='y'
    fi
    
    if [ "$i" == "simulator" ]; then        
        COMPILE_SIMULATOR='y'
    fi
    
    if [ "$i" == "gui" ]; then        
        COMPILE_GUI='y'
    fi
    
    if [ "$i" == "batch" ]; then        
        COMPILE_MPI='y'
    fi

     if [ "$i" == "converter" ]; then
        COMPILE_CONVERTER='y'
    fi
    
done

DEFAULT_C_FLAGS="-fopenmp -std=gnu99 -fno-strict-aliasing  -Wall -Wno-stringop-truncation -Wno-unused-function -Wno-char-subscripts -Wno-unused-result"
RUNTIME_OUTPUT_DIRECTORY="$ROOT_DIR/bin"
LIBRARY_OUTPUT_DIRECTORY="$ROOT_DIR/shared_libs"

C_FLAGS="$C_FLAGS $DEFAULT_C_FLAGS"

GET_LINUX_VERSION
        
if [ "$OS" == "Manjaro Linux" ]; then
    C_COMPILER=/opt/cuda/bin/gcc
    CXX_COMPILER=/opt/cuda/bin/g++
elif [ "$OS" == "Ubuntu" ]; then
    C_COMPILER=gcc-8
    CXX_COMPILER=g++-8
elif [ "$OS" == "Fedora" ]; then
    C_COMPILER=/usr/local/cuda/bin/gcc
    CXX_COMPILER=/usr/local/cuda/bin/g++
fi

if [ -n "$COMPILE_SIMULATOR" ] || [ -n "$COMPILE_BATCH" ]; then
    
    FIND_CUDA
    
    echo -e "${INFO}C compiler:${NC} $C_COMPILER"
    echo -e "${INFO}C++ compiler:${NC} $CXX_COMPILER"

    if [ -n "$CUDA_FOUND" ]; then
        echo -e "${INFO}CUDA compiler:${NC} $NVCC"
        echo -e "${INFO}CUDA libraries path:${NC} $CUDA_LIBRARY_PATH"
        echo -e "${INFO}CUDA include path:${NC} $CUDA_INCLUDE_PATH"
       
        
        C_FLAGS="${C_FLAGS} -DCOMPILE_CUDA -I${CUDA_INCLUDE_PATH}"
        
    fi
fi

if [ -n "$COMPILE_GUI" ]; then
    C_FLAGS="${C_FLAGS} -DCOMPILE_GUI"
fi

echo -e "${INFO}C FLAGS:${NC} $C_FLAGS"

ADD_SUBDIRECTORY "src/3dparty/sds"
ADD_SUBDIRECTORY "src/config_helpers"
ADD_SUBDIRECTORY "src/utils"
ADD_SUBDIRECTORY "src/alg"
ADD_SUBDIRECTORY "src/monodomain"
ADD_SUBDIRECTORY "src/ode_solver"
ADD_SUBDIRECTORY "src/3dparty/ini_parser"
ADD_SUBDIRECTORY "src/config"
ADD_SUBDIRECTORY "src/graph"
ADD_SUBDIRECTORY "src/3dparty/xml_parser"
ADD_SUBDIRECTORY "src/3dparty/tinyexpr"
ADD_SUBDIRECTORY "src/3dparty/miniz"
ADD_SUBDIRECTORY "src/vtk_utils"


#DINAMIC DEPS
ADD_SUBDIRECTORY "src/logger"

if [ -n "$COMPILE_GUI" ]; then
    ADD_SUBDIRECTORY "src/3dparty/raylib/src"
    ADD_SUBDIRECTORY "src/3dparty/tinyfiledialogs"
    ADD_SUBDIRECTORY "src/gui"
    OPT_DEPS="gui raylib tinyfiledialogs"
fi

if [ -n "$CUDA_FOUND" ]; then
    ADD_SUBDIRECTORY "src/gpu_utils"
    OPT_DEPS="${OPT_DEPS} gpu_utils"
fi

SRC_FILES="src/main_simulator.c"
HDR_FILES=""

STATIC_DEPS="monodomain ode_solver ini_parser config tinyexpr ${OPT_DEPS} config_helpers vtk_utils yxml alg graph utils sds miniz"
DYNAMIC_DEPS="dl m $CUDA_LIBRARIES"

if [ -n "$COMPILE_GUI" ]; then
    DYNAMIC_DEPS="$DYNAMIC_DEPS OpenGL GLX GLU pthread X11 rt"
fi

DYNAMIC_DEPS="$DYNAMIC_DEPS logger"

if [ -n "$COMPILE_SIMULATOR" ] || [ -n "$COMPILE_BATCH" ]; then
    ADD_SUBDIRECTORY "src/models_library"
    ADD_SUBDIRECTORY "src/stimuli_library"
    ADD_SUBDIRECTORY "src/domains_library"
    ADD_SUBDIRECTORY "src/purkinje_library"
    ADD_SUBDIRECTORY "src/extra_data_library"
    ADD_SUBDIRECTORY "src/matrix_assembly_library"
    ADD_SUBDIRECTORY "src/linear_system_solver_library"
    ADD_SUBDIRECTORY "src/save_mesh_library"
    ADD_SUBDIRECTORY "src/save_state_library"
    ADD_SUBDIRECTORY "src/restore_state_library"
    ADD_SUBDIRECTORY "src/update_monodomain_library"
    ADD_SUBDIRECTORY "src/modify_domain"
fi

#COMPILE THE EXECUTABLES NOW

if [ -n "$COMPILE_SIMULATOR" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D" "$SRC_FILES" "$HDR_FILES" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $EXTRA_LIB_PATH $LIBRARY_OUTPUT_DIRECTORY"
fi

if [ -n "$COMPILE_MPI" ]; then

  FIND_MPI

  if [ -n "$MPI_FOUND" ]; then
      SRC_FILES="src/main_batch.c"
      HDR_FILES=""
      DYNAMIC_DEPS="$DYNAMIC_DEPS $MPI_LIBRARIES"
      EXTRA_LIB_PATH="$EXTRA_LIB_PATH $CUDA_LIBRARY_PATH $MPI_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY"

      if [ -n "$MPI_INCLUDE_PATH" ]; then
        INCLUDE_P="-I$MPI_INCLUDE_PATH"
      fi
      COMPILE_EXECUTABLE "MonoAlg3D_batch" "$SRC_FILES" "$HDR_FILES" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXTRA_LIB_PATH" "$INCLUDE_P"
  fi

fi

if [ -n "$COMPILE_GUI" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D_visualizer" "src/main_visualizer.c" "" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXTRA_LIB_PATH $LIBRARY_OUTPUT_DIRECTORY"
fi

if [ -n "$COMPILE_CONVERTER" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D_converter" "src/main_converter.c" "" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXTRA_LIB_PATH $LIBRARY_OUTPUT_DIRECTORY"
fi

FIND_CRITERION

if [ -n "$CRITERION_FOUND" ]; then
    # shellcheck disable=SC2034
    RUNTIME_OUTPUT_DIRECTORY=$ROOT_DIR/tests_bin
    ADD_SUBDIRECTORY "src/tests"
fi
