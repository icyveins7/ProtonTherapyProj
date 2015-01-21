//64Bit Compiled DLL for 1F1 Hypergeometric Function

#ifdef _WIN32
#ifdef EXPORT_FCNS
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION __declspec(dllimport)
#endif
#else
#define EXPORTED_FUNCTION
#endif


EXPORTED_FUNCTION double F11(double a, double b, double z, int precision);