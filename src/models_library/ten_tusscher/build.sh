############## TEN TUSCHER 2004 ##############################
MODEL_FILE_CPU="ten_tusscher_2004_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_2004_RS_GPU.cu"
MODEL_FILE_CPP="ten_tusscher_2004_RS_SYCL.cpp"
COMMON_HEADERS="ten_tusscher_2004.h"

COMPILE_MODEL_LIB "ten_tusscher_2004_endo" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_CPP" "$COMMON_HEADERS"
##########################################################

