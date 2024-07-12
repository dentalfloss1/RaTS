#define PY_SSIZE_T_CLEAN  /* Make "s#" use Py_ssize_t rather than int. */
#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>

struct sbpl_params  
{
    double alpha1, alpha2, s, duration, tc ,fpk, nu, nu0, beta;  
};

static double sbpl (double t, void *params)
{
  struct sbpl_params *p
    = (struct sbpl_params *) params;

  double alpha1 = (p->alpha1);
  double alpha2 = (p->alpha2);
  double s = (p->s);
  double duration = (p->duration);
  double tc = (p->tc);
  double fpk = (p->fpk);
  double nu = (p->nu);
  double nu0 = (p->nu0);
  double beta = (p->beta);

  return pow(2,1/s) * fpk * pow(nu/nu0, beta) * pow( pow(t/duration,-s*alpha1) + pow(t/duration,-s*alpha2), -1/s );
}

static PyObject *SBPLError;

static PyObject *
integrate_sbpl(PyObject *self, PyObject *args)
{
    // const char *command;
    int sts;
    double alpha1, alpha2, s, duration, endobs, tc, fpk, nu, nu0, beta, startobs;
    if (!PyArg_ParseTuple(args, "(ddddddddddd)", &alpha1, &alpha2, &s, &duration, &endobs, &tc, &fpk, &nu, &nu0, &beta, &startobs))
        return NULL;
    // sts = printf("%.17e,%.17e,%.17e,%.17e,%.17e\n",&duration, &endobs, &tc, &fpk, &nu, &nu0, &beta, &startobs);
   /* if (sts < 0) {
        PyErr_SetString(SpamError, "System command failed");
        return NULL;
    }*/
    gsl_integration_workspace * w
      = gsl_integration_workspace_alloc (50);
    double error, result;
    
    gsl_function F;
    F.function = &sbpl;
    
    struct sbpl_params params = {alpha1, 
      alpha2, 
      s, 
      duration, 
      tc,
      fpk,
      nu,
      nu0,
      beta
      };
    F.params = &params;
    if (tc < endobs){
      double tstart = fmax(startobs,tc); 
      sts = gsl_integration_qags (&F, tstart - tc, endobs - tc, 0, 1e-7, 50,
                            w, &result, &error);
    }
    else {
      result = 0;
      sts = 0;
    }
    if (sts < 0) {
        PyErr_SetString(SBPLError, "gsl_integration_qags failed");
        return NULL;
    }
    // printf("%f,%f\n",result,error);

    return Py_BuildValue("dd",result,error );
}


static PyMethodDef keywdarg_methods[] = {
    /* The cast of the function is necessary since PyCFunction values
     * only take two PyObject* parameters, and keywdarg_parrot() takes
     * three.
     */
    {"integrate",  integrate_sbpl, METH_VARARGS,
     "Execute a shell command."},
    {NULL, NULL, 0, NULL}   /* sentinel */
};

static struct PyModuleDef keywdargmodule = {
    PyModuleDef_HEAD_INIT,
    "sbpl",
    NULL,
    -1,
    keywdarg_methods
};

PyMODINIT_FUNC
PyInit_sbpl(void)
{
    PyObject *m;

    m = PyModule_Create(&keywdargmodule);
    if (m == NULL)
        return NULL;

    SBPLError = PyErr_NewException("sbpl.error", NULL, NULL);
    Py_XINCREF(SBPLError);
    if (PyModule_AddObject(m, "error", SBPLError) < 0) {
        Py_XDECREF(SBPLError);
        Py_CLEAR(SBPLError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
