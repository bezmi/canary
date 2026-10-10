// Checks the compile-fail test mechanism itself: this file compiles as
// written and fails when GASLIB_COMPILE_FAIL is defined.
#ifdef GASLIB_COMPILE_FAIL
static_assert(false, "compile-fail self test");
#endif
