############## TEN TUSCHER 3 ##############################
MODEL_FILE_CPU="ten_tusscher_3_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_3_RS_GPU.cu"
MODEL_FILE_SYCL="ten_tusscher_3_RS_SYCL.cpp"
COMMON_HEADERS="ten_tusscher_3_RS.h"

COMPILE_MODEL_LIB "ten_tusscher_3_with_extra_params" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
##########################################################

