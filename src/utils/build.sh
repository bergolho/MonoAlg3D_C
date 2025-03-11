UTILS_SOURCE_FILES="search.c stop_watch.cpp sort.cpp file_utils.cpp batch_utils.c"
UTILS_HEADER_FILES="utils.h stop_watch.h file_utils.h batch_utils.c"

COMPILE_STATIC_LIB "utils" "$UTILS_SOURCE_FILES" "$UTILS_HEADER_FILES"
