############## MITCHELL SHAEFFER 2003 ##############################
MODEL_FILE_CPU="mitchell_shaeffer_2003.c"
MODEL_FILE_GPU="mitchell_shaeffer_2003.cu"
<<<<<<< HEAD
MODEL_FILE_SYCL="mitchell_shaeffer_2003.cpp"
=======
MODEL_FILE_SYCL=""
>>>>>>> 89b4c2d38ee29f48ec95a92ecece5c83b3c040bf
COMMON_HEADERS="mitchell_shaeffer_2003.h"

COMPILE_MODEL_LIB "mitchell_shaeffer_2003" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
#########################################################