/* Generated by Cython 0.29.12 */

#ifndef __PYX_HAVE__sage__ext__interpreters__wrapper_cc
#define __PYX_HAVE__sage__ext__interpreters__wrapper_cc


#ifndef __PYX_HAVE_API__sage__ext__interpreters__wrapper_cc

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

#ifndef DL_IMPORT
  #define DL_IMPORT(_T) _T
#endif

__PYX_EXTERN_C int cc_py_call_helper(PyObject *, PyObject *, int, mpc_t *, __mpc_struct *);

#endif /* !__PYX_HAVE_API__sage__ext__interpreters__wrapper_cc */

/* WARNING: the interface of the module init function changed in CPython 3.5. */
/* It now returns a PyModuleDef instance instead of a PyModule instance. */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initwrapper_cc(void);
#else
PyMODINIT_FUNC PyInit_wrapper_cc(void);
#endif

#endif /* !__PYX_HAVE__sage__ext__interpreters__wrapper_cc */
