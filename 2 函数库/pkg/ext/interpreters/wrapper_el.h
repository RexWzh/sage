/* Generated by Cython 0.29.12 */

#ifndef __PYX_HAVE__sage__ext__interpreters__wrapper_el
#define __PYX_HAVE__sage__ext__interpreters__wrapper_el


#ifndef __PYX_HAVE_API__sage__ext__interpreters__wrapper_el

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

__PYX_EXTERN_C PyObject *py_divide_helper(PyObject *, PyObject *);
__PYX_EXTERN_C PyObject *el_check_element(PyObject *, PyObject *);

#endif /* !__PYX_HAVE_API__sage__ext__interpreters__wrapper_el */

/* WARNING: the interface of the module init function changed in CPython 3.5. */
/* It now returns a PyModuleDef instance instead of a PyModule instance. */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initwrapper_el(void);
#else
PyMODINIT_FUNC PyInit_wrapper_el(void);
#endif

#endif /* !__PYX_HAVE__sage__ext__interpreters__wrapper_el */
