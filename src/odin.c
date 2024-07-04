// Automatically generated by odin 1.2.5 - do not edit
#include <float.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <R_ext/Rdynload.h>
typedef struct model_internal {
  double a;
  double AR0;
  double As0;
  double b;
  double c_A;
  double c_D;
  double c_T;
  int *delay_index_delayed_Lambda_v_r_Sv;
  int *delay_index_delayed_Lambda_v_s_Sv;
  double *delay_state_delayed_Lambda_v_r_Sv;
  double *delay_state_delayed_Lambda_v_s_Sv;
  double DR0;
  double Ds0;
  double e;
  double Ev_r0;
  double Ev_s0;
  double fT;
  double initial_AR;
  double initial_As;
  double initial_DR;
  double initial_Ds;
  double initial_Ev_r;
  double initial_Ev_s;
  double initial_Iv_r;
  double initial_Iv_s;
  double initial_S;
  double initial_Sv;
  double initial_t;
  double initial_TR;
  double initial_Ts;
  double Iv_r0;
  double Iv_s0;
  double m;
  double mu;
  double n;
  double Phi;
  double rA;
  double rD;
  double res_start;
  double res_time;
  double rTR_true;
  double rTs;
  double S0;
  double Sv0;
  double toff;
  double ton;
  double TR0;
  double Ts0;
  bool model_use_dde;
} model_internal;
model_internal* model_get_internal(SEXP internal_p, int closed_error);
static void model_finalise(SEXP internal_p);
SEXP model_create(SEXP user);
void model_initmod_desolve(void(* odeparms) (int *, double *));
SEXP model_contents(SEXP internal_p);
SEXP model_set_user(SEXP internal_p, SEXP user);
SEXP model_set_initial(SEXP internal_p, SEXP t_ptr, SEXP state_ptr, SEXP model_use_dde_ptr);
SEXP model_metadata(SEXP internal_p);
SEXP model_initial_conditions(SEXP internal_p, SEXP t_ptr);
void model_rhs(model_internal* internal, double t, double * state, double * dstatedt, double * output);
void model_rhs_dde(size_t neq, double t, double * state, double * dstatedt, void * internal);
void model_rhs_desolve(int * neq, double * t, double * state, double * dstatedt, double * output, int * np);
void model_output_dde(size_t n_eq, double t, double * state, size_t n_output, double * output, void * internal_p);
SEXP model_rhs_r(SEXP internal_p, SEXP t, SEXP state);
double user_get_scalar_double(SEXP user, const char *name,
                              double default_value, double min, double max);
int user_get_scalar_int(SEXP user, const char *name,
                        int default_value, double min, double max);
void user_check_values_double(double * value, size_t len,
                                  double min, double max, const char *name);
void user_check_values_int(int * value, size_t len,
                               double min, double max, const char *name);
void user_check_values(SEXP value, double min, double max,
                           const char *name);
SEXP user_list_element(SEXP list, const char *name);
void lagvalue(double t, bool use_dde, int *idx, int dim_idx, double *state);
void lagvalue_dde(double t, int *idx, size_t dim_idx, double *state);
void lagvalue_ds(double t, int *idx, int dim_idx, double *state);
double scalar_real(SEXP x, const char * name);
model_internal* model_get_internal(SEXP internal_p, int closed_error) {
  model_internal *internal = NULL;
  if (TYPEOF(internal_p) != EXTPTRSXP) {
    Rf_error("Expected an external pointer");
  }
  internal = (model_internal*) R_ExternalPtrAddr(internal_p);
  if (!internal && closed_error) {
    Rf_error("Pointer has been invalidated");
  }
  return internal;
}
void model_finalise(SEXP internal_p) {
  model_internal *internal = model_get_internal(internal_p, 0);
  if (internal_p) {
    R_Free(internal->delay_index_delayed_Lambda_v_r_Sv);
    R_Free(internal->delay_index_delayed_Lambda_v_s_Sv);
    R_Free(internal->delay_state_delayed_Lambda_v_r_Sv);
    R_Free(internal->delay_state_delayed_Lambda_v_s_Sv);
    R_Free(internal);
    R_ClearExternalPtr(internal_p);
  }
}
SEXP model_create(SEXP user) {
  model_internal *internal = (model_internal*) R_Calloc(1, model_internal);
  internal->delay_index_delayed_Lambda_v_r_Sv = NULL;
  internal->delay_index_delayed_Lambda_v_s_Sv = NULL;
  internal->delay_state_delayed_Lambda_v_r_Sv = NULL;
  internal->delay_state_delayed_Lambda_v_s_Sv = NULL;
  R_Free(internal->delay_index_delayed_Lambda_v_r_Sv);
  internal->delay_index_delayed_Lambda_v_r_Sv = R_Calloc(4, int);
  R_Free(internal->delay_state_delayed_Lambda_v_r_Sv);
  internal->delay_state_delayed_Lambda_v_r_Sv = R_Calloc(4, double);
  internal->delay_index_delayed_Lambda_v_r_Sv[0] = 4;
  internal->delay_index_delayed_Lambda_v_r_Sv[1] = 5;
  internal->delay_index_delayed_Lambda_v_r_Sv[2] = 6;
  internal->delay_index_delayed_Lambda_v_r_Sv[3] = 7;
  R_Free(internal->delay_index_delayed_Lambda_v_s_Sv);
  internal->delay_index_delayed_Lambda_v_s_Sv = R_Calloc(4, int);
  R_Free(internal->delay_state_delayed_Lambda_v_s_Sv);
  internal->delay_state_delayed_Lambda_v_s_Sv = R_Calloc(4, double);
  internal->delay_index_delayed_Lambda_v_s_Sv[0] = 1;
  internal->delay_index_delayed_Lambda_v_s_Sv[1] = 2;
  internal->delay_index_delayed_Lambda_v_s_Sv[2] = 3;
  internal->delay_index_delayed_Lambda_v_s_Sv[3] = 7;
  internal->a = NA_REAL;
  internal->AR0 = NA_REAL;
  internal->As0 = NA_REAL;
  internal->b = NA_REAL;
  internal->c_A = NA_REAL;
  internal->c_D = NA_REAL;
  internal->c_T = NA_REAL;
  internal->DR0 = NA_REAL;
  internal->Ds0 = NA_REAL;
  internal->e = NA_REAL;
  internal->Ev_r0 = NA_REAL;
  internal->Ev_s0 = NA_REAL;
  internal->fT = NA_REAL;
  internal->Iv_r0 = NA_REAL;
  internal->Iv_s0 = NA_REAL;
  internal->m = NA_REAL;
  internal->mu = NA_REAL;
  internal->n = NA_REAL;
  internal->Phi = NA_REAL;
  internal->rA = NA_REAL;
  internal->rD = NA_REAL;
  internal->res_start = NA_REAL;
  internal->res_time = NA_REAL;
  internal->rTR_true = NA_REAL;
  internal->rTs = NA_REAL;
  internal->S0 = NA_REAL;
  internal->Sv0 = NA_REAL;
  internal->toff = NA_REAL;
  internal->ton = NA_REAL;
  internal->TR0 = NA_REAL;
  internal->Ts0 = NA_REAL;
  internal->initial_t = NA_REAL;
  SEXP ptr = PROTECT(R_MakeExternalPtr(internal, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(ptr, model_finalise);
  UNPROTECT(1);
  return ptr;
}
static model_internal *model_internal_ds;
void model_initmod_desolve(void(* odeparms) (int *, double *)) {
  static DL_FUNC get_desolve_gparms = NULL;
  if (get_desolve_gparms == NULL) {
    get_desolve_gparms =
      R_GetCCallable("deSolve", "get_deSolve_gparms");
  }
  model_internal_ds = model_get_internal(get_desolve_gparms(), 1);
}
SEXP model_contents(SEXP internal_p) {
  model_internal *internal = model_get_internal(internal_p, 1);
  SEXP contents = PROTECT(allocVector(VECSXP, 49));
  SET_VECTOR_ELT(contents, 0, ScalarReal(internal->a));
  SET_VECTOR_ELT(contents, 1, ScalarReal(internal->AR0));
  SET_VECTOR_ELT(contents, 2, ScalarReal(internal->As0));
  SET_VECTOR_ELT(contents, 3, ScalarReal(internal->b));
  SET_VECTOR_ELT(contents, 4, ScalarReal(internal->c_A));
  SET_VECTOR_ELT(contents, 5, ScalarReal(internal->c_D));
  SET_VECTOR_ELT(contents, 6, ScalarReal(internal->c_T));
  SEXP delay_index_delayed_Lambda_v_r_Sv = PROTECT(allocVector(INTSXP, 4));
  memcpy(INTEGER(delay_index_delayed_Lambda_v_r_Sv), internal->delay_index_delayed_Lambda_v_r_Sv, 4 * sizeof(int));
  SET_VECTOR_ELT(contents, 7, delay_index_delayed_Lambda_v_r_Sv);
  SEXP delay_index_delayed_Lambda_v_s_Sv = PROTECT(allocVector(INTSXP, 4));
  memcpy(INTEGER(delay_index_delayed_Lambda_v_s_Sv), internal->delay_index_delayed_Lambda_v_s_Sv, 4 * sizeof(int));
  SET_VECTOR_ELT(contents, 8, delay_index_delayed_Lambda_v_s_Sv);
  SEXP delay_state_delayed_Lambda_v_r_Sv = PROTECT(allocVector(REALSXP, 4));
  memcpy(REAL(delay_state_delayed_Lambda_v_r_Sv), internal->delay_state_delayed_Lambda_v_r_Sv, 4 * sizeof(double));
  SET_VECTOR_ELT(contents, 9, delay_state_delayed_Lambda_v_r_Sv);
  SEXP delay_state_delayed_Lambda_v_s_Sv = PROTECT(allocVector(REALSXP, 4));
  memcpy(REAL(delay_state_delayed_Lambda_v_s_Sv), internal->delay_state_delayed_Lambda_v_s_Sv, 4 * sizeof(double));
  SET_VECTOR_ELT(contents, 10, delay_state_delayed_Lambda_v_s_Sv);
  SET_VECTOR_ELT(contents, 11, ScalarReal(internal->DR0));
  SET_VECTOR_ELT(contents, 12, ScalarReal(internal->Ds0));
  SET_VECTOR_ELT(contents, 13, ScalarReal(internal->e));
  SET_VECTOR_ELT(contents, 14, ScalarReal(internal->Ev_r0));
  SET_VECTOR_ELT(contents, 15, ScalarReal(internal->Ev_s0));
  SET_VECTOR_ELT(contents, 16, ScalarReal(internal->fT));
  SET_VECTOR_ELT(contents, 17, ScalarReal(internal->initial_AR));
  SET_VECTOR_ELT(contents, 18, ScalarReal(internal->initial_As));
  SET_VECTOR_ELT(contents, 19, ScalarReal(internal->initial_DR));
  SET_VECTOR_ELT(contents, 20, ScalarReal(internal->initial_Ds));
  SET_VECTOR_ELT(contents, 21, ScalarReal(internal->initial_Ev_r));
  SET_VECTOR_ELT(contents, 22, ScalarReal(internal->initial_Ev_s));
  SET_VECTOR_ELT(contents, 23, ScalarReal(internal->initial_Iv_r));
  SET_VECTOR_ELT(contents, 24, ScalarReal(internal->initial_Iv_s));
  SET_VECTOR_ELT(contents, 25, ScalarReal(internal->initial_S));
  SET_VECTOR_ELT(contents, 26, ScalarReal(internal->initial_Sv));
  SET_VECTOR_ELT(contents, 27, ScalarReal(internal->initial_t));
  SET_VECTOR_ELT(contents, 28, ScalarReal(internal->initial_TR));
  SET_VECTOR_ELT(contents, 29, ScalarReal(internal->initial_Ts));
  SET_VECTOR_ELT(contents, 30, ScalarReal(internal->Iv_r0));
  SET_VECTOR_ELT(contents, 31, ScalarReal(internal->Iv_s0));
  SET_VECTOR_ELT(contents, 32, ScalarReal(internal->m));
  SET_VECTOR_ELT(contents, 33, ScalarReal(internal->mu));
  SET_VECTOR_ELT(contents, 34, ScalarReal(internal->n));
  SET_VECTOR_ELT(contents, 35, ScalarReal(internal->Phi));
  SET_VECTOR_ELT(contents, 36, ScalarReal(internal->rA));
  SET_VECTOR_ELT(contents, 37, ScalarReal(internal->rD));
  SET_VECTOR_ELT(contents, 38, ScalarReal(internal->res_start));
  SET_VECTOR_ELT(contents, 39, ScalarReal(internal->res_time));
  SET_VECTOR_ELT(contents, 40, ScalarReal(internal->rTR_true));
  SET_VECTOR_ELT(contents, 41, ScalarReal(internal->rTs));
  SET_VECTOR_ELT(contents, 42, ScalarReal(internal->S0));
  SET_VECTOR_ELT(contents, 43, ScalarReal(internal->Sv0));
  SET_VECTOR_ELT(contents, 44, ScalarReal(internal->toff));
  SET_VECTOR_ELT(contents, 45, ScalarReal(internal->ton));
  SET_VECTOR_ELT(contents, 46, ScalarReal(internal->TR0));
  SET_VECTOR_ELT(contents, 47, ScalarReal(internal->Ts0));
  SET_VECTOR_ELT(contents, 48, ScalarLogical(internal->model_use_dde));
  SEXP nms = PROTECT(allocVector(STRSXP, 49));
  SET_STRING_ELT(nms, 0, mkChar("a"));
  SET_STRING_ELT(nms, 1, mkChar("AR0"));
  SET_STRING_ELT(nms, 2, mkChar("As0"));
  SET_STRING_ELT(nms, 3, mkChar("b"));
  SET_STRING_ELT(nms, 4, mkChar("c_A"));
  SET_STRING_ELT(nms, 5, mkChar("c_D"));
  SET_STRING_ELT(nms, 6, mkChar("c_T"));
  SET_STRING_ELT(nms, 7, mkChar("delay_index_delayed_Lambda_v_r_Sv"));
  SET_STRING_ELT(nms, 8, mkChar("delay_index_delayed_Lambda_v_s_Sv"));
  SET_STRING_ELT(nms, 9, mkChar("delay_state_delayed_Lambda_v_r_Sv"));
  SET_STRING_ELT(nms, 10, mkChar("delay_state_delayed_Lambda_v_s_Sv"));
  SET_STRING_ELT(nms, 11, mkChar("DR0"));
  SET_STRING_ELT(nms, 12, mkChar("Ds0"));
  SET_STRING_ELT(nms, 13, mkChar("e"));
  SET_STRING_ELT(nms, 14, mkChar("Ev_r0"));
  SET_STRING_ELT(nms, 15, mkChar("Ev_s0"));
  SET_STRING_ELT(nms, 16, mkChar("fT"));
  SET_STRING_ELT(nms, 17, mkChar("initial_AR"));
  SET_STRING_ELT(nms, 18, mkChar("initial_As"));
  SET_STRING_ELT(nms, 19, mkChar("initial_DR"));
  SET_STRING_ELT(nms, 20, mkChar("initial_Ds"));
  SET_STRING_ELT(nms, 21, mkChar("initial_Ev_r"));
  SET_STRING_ELT(nms, 22, mkChar("initial_Ev_s"));
  SET_STRING_ELT(nms, 23, mkChar("initial_Iv_r"));
  SET_STRING_ELT(nms, 24, mkChar("initial_Iv_s"));
  SET_STRING_ELT(nms, 25, mkChar("initial_S"));
  SET_STRING_ELT(nms, 26, mkChar("initial_Sv"));
  SET_STRING_ELT(nms, 27, mkChar("initial_t"));
  SET_STRING_ELT(nms, 28, mkChar("initial_TR"));
  SET_STRING_ELT(nms, 29, mkChar("initial_Ts"));
  SET_STRING_ELT(nms, 30, mkChar("Iv_r0"));
  SET_STRING_ELT(nms, 31, mkChar("Iv_s0"));
  SET_STRING_ELT(nms, 32, mkChar("m"));
  SET_STRING_ELT(nms, 33, mkChar("mu"));
  SET_STRING_ELT(nms, 34, mkChar("n"));
  SET_STRING_ELT(nms, 35, mkChar("Phi"));
  SET_STRING_ELT(nms, 36, mkChar("rA"));
  SET_STRING_ELT(nms, 37, mkChar("rD"));
  SET_STRING_ELT(nms, 38, mkChar("res_start"));
  SET_STRING_ELT(nms, 39, mkChar("res_time"));
  SET_STRING_ELT(nms, 40, mkChar("rTR_true"));
  SET_STRING_ELT(nms, 41, mkChar("rTs"));
  SET_STRING_ELT(nms, 42, mkChar("S0"));
  SET_STRING_ELT(nms, 43, mkChar("Sv0"));
  SET_STRING_ELT(nms, 44, mkChar("toff"));
  SET_STRING_ELT(nms, 45, mkChar("ton"));
  SET_STRING_ELT(nms, 46, mkChar("TR0"));
  SET_STRING_ELT(nms, 47, mkChar("Ts0"));
  SET_STRING_ELT(nms, 48, mkChar("model_use_dde"));
  setAttrib(contents, R_NamesSymbol, nms);
  UNPROTECT(6);
  return contents;
}
SEXP model_set_user(SEXP internal_p, SEXP user) {
  model_internal *internal = model_get_internal(internal_p, 1);
  internal->a = user_get_scalar_double(user, "a", internal->a, NA_REAL, NA_REAL);
  internal->AR0 = user_get_scalar_double(user, "AR0", internal->AR0, NA_REAL, NA_REAL);
  internal->As0 = user_get_scalar_double(user, "As0", internal->As0, NA_REAL, NA_REAL);
  internal->b = user_get_scalar_double(user, "b", internal->b, NA_REAL, NA_REAL);
  internal->c_A = user_get_scalar_double(user, "c_A", internal->c_A, NA_REAL, NA_REAL);
  internal->c_D = user_get_scalar_double(user, "c_D", internal->c_D, NA_REAL, NA_REAL);
  internal->c_T = user_get_scalar_double(user, "c_T", internal->c_T, NA_REAL, NA_REAL);
  internal->DR0 = user_get_scalar_double(user, "DR0", internal->DR0, NA_REAL, NA_REAL);
  internal->Ds0 = user_get_scalar_double(user, "Ds0", internal->Ds0, NA_REAL, NA_REAL);
  internal->e = user_get_scalar_double(user, "e", internal->e, NA_REAL, NA_REAL);
  internal->Ev_r0 = user_get_scalar_double(user, "Ev_r0", internal->Ev_r0, NA_REAL, NA_REAL);
  internal->Ev_s0 = user_get_scalar_double(user, "Ev_s0", internal->Ev_s0, NA_REAL, NA_REAL);
  internal->fT = user_get_scalar_double(user, "fT", internal->fT, NA_REAL, NA_REAL);
  internal->Iv_r0 = user_get_scalar_double(user, "Iv_r0", internal->Iv_r0, NA_REAL, NA_REAL);
  internal->Iv_s0 = user_get_scalar_double(user, "Iv_s0", internal->Iv_s0, NA_REAL, NA_REAL);
  internal->m = user_get_scalar_double(user, "m", internal->m, NA_REAL, NA_REAL);
  internal->mu = user_get_scalar_double(user, "mu", internal->mu, NA_REAL, NA_REAL);
  internal->n = user_get_scalar_double(user, "n", internal->n, NA_REAL, NA_REAL);
  internal->Phi = user_get_scalar_double(user, "Phi", internal->Phi, NA_REAL, NA_REAL);
  internal->rA = user_get_scalar_double(user, "rA", internal->rA, NA_REAL, NA_REAL);
  internal->rD = user_get_scalar_double(user, "rD", internal->rD, NA_REAL, NA_REAL);
  internal->res_start = user_get_scalar_double(user, "res_start", internal->res_start, NA_REAL, NA_REAL);
  internal->res_time = user_get_scalar_double(user, "res_time", internal->res_time, NA_REAL, NA_REAL);
  internal->rTR_true = user_get_scalar_double(user, "rTR_true", internal->rTR_true, NA_REAL, NA_REAL);
  internal->rTs = user_get_scalar_double(user, "rTs", internal->rTs, NA_REAL, NA_REAL);
  internal->S0 = user_get_scalar_double(user, "S0", internal->S0, NA_REAL, NA_REAL);
  internal->Sv0 = user_get_scalar_double(user, "Sv0", internal->Sv0, NA_REAL, NA_REAL);
  internal->toff = user_get_scalar_double(user, "toff", internal->toff, NA_REAL, NA_REAL);
  internal->ton = user_get_scalar_double(user, "ton", internal->ton, NA_REAL, NA_REAL);
  internal->TR0 = user_get_scalar_double(user, "TR0", internal->TR0, NA_REAL, NA_REAL);
  internal->Ts0 = user_get_scalar_double(user, "Ts0", internal->Ts0, NA_REAL, NA_REAL);
  internal->initial_AR = (internal->res_time == 0 ? (internal->DR0 + internal->AR0 + internal->TR0) * internal->res_start : 9.9999999999999995e-07);
  internal->initial_As = internal->As0;
  internal->initial_DR = (internal->res_time == 0 ? (internal->DR0 + internal->AR0 + internal->TR0) * internal->res_start : 9.9999999999999995e-07);
  internal->initial_Ds = internal->Ds0;
  internal->initial_Ev_r = (internal->res_time == 0 ? (internal->Ev_r0 + internal->Iv_r0) * internal->res_start : 9.9999999999999995e-07);
  internal->initial_Ev_s = internal->Ev_s0;
  internal->initial_Iv_r = (internal->res_time == 0 ? (internal->Ev_r0 + internal->Iv_r0) * internal->res_start : 9.9999999999999995e-07);
  internal->initial_Iv_s = internal->Iv_s0;
  internal->initial_S = internal->S0;
  internal->initial_Sv = internal->Sv0;
  internal->initial_TR = (internal->res_time == 0 ? (internal->DR0 + internal->AR0 + internal->TR0) * internal->res_start : 9.9999999999999995e-07);
  internal->initial_Ts = internal->Ts0;
  return R_NilValue;
}
SEXP model_set_initial(SEXP internal_p, SEXP t_ptr, SEXP state_ptr, SEXP model_use_dde_ptr) {
  model_internal *internal = model_get_internal(internal_p, 1);
  const double t = REAL(t_ptr)[0];
  internal->initial_t = t;
  internal->model_use_dde = INTEGER(model_use_dde_ptr)[0];
  if (state_ptr != R_NilValue) {
    double * state = REAL(state_ptr);
    internal->initial_S = state[0];
    internal->initial_Ds = state[1];
    internal->initial_As = state[2];
    internal->initial_Ts = state[3];
    internal->initial_DR = state[4];
    internal->initial_AR = state[5];
    internal->initial_TR = state[6];
    internal->initial_Sv = state[7];
    internal->initial_Ev_s = state[8];
    internal->initial_Iv_s = state[9];
    internal->initial_Ev_r = state[10];
    internal->initial_Iv_r = state[11];
  }
  return R_NilValue;
}
SEXP model_metadata(SEXP internal_p) {
  model_internal *internal = model_get_internal(internal_p, 1);
  SEXP ret = PROTECT(allocVector(VECSXP, 4));
  SEXP nms = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(nms, 0, mkChar("variable_order"));
  SET_STRING_ELT(nms, 1, mkChar("output_order"));
  SET_STRING_ELT(nms, 2, mkChar("n_out"));
  SET_STRING_ELT(nms, 3, mkChar("interpolate_t"));
  setAttrib(ret, R_NamesSymbol, nms);
  SEXP variable_length = PROTECT(allocVector(VECSXP, 12));
  SEXP variable_names = PROTECT(allocVector(STRSXP, 12));
  setAttrib(variable_length, R_NamesSymbol, variable_names);
  SET_VECTOR_ELT(variable_length, 0, R_NilValue);
  SET_VECTOR_ELT(variable_length, 1, R_NilValue);
  SET_VECTOR_ELT(variable_length, 2, R_NilValue);
  SET_VECTOR_ELT(variable_length, 3, R_NilValue);
  SET_VECTOR_ELT(variable_length, 4, R_NilValue);
  SET_VECTOR_ELT(variable_length, 5, R_NilValue);
  SET_VECTOR_ELT(variable_length, 6, R_NilValue);
  SET_VECTOR_ELT(variable_length, 7, R_NilValue);
  SET_VECTOR_ELT(variable_length, 8, R_NilValue);
  SET_VECTOR_ELT(variable_length, 9, R_NilValue);
  SET_VECTOR_ELT(variable_length, 10, R_NilValue);
  SET_VECTOR_ELT(variable_length, 11, R_NilValue);
  SET_STRING_ELT(variable_names, 0, mkChar("S"));
  SET_STRING_ELT(variable_names, 1, mkChar("Ds"));
  SET_STRING_ELT(variable_names, 2, mkChar("As"));
  SET_STRING_ELT(variable_names, 3, mkChar("Ts"));
  SET_STRING_ELT(variable_names, 4, mkChar("DR"));
  SET_STRING_ELT(variable_names, 5, mkChar("AR"));
  SET_STRING_ELT(variable_names, 6, mkChar("TR"));
  SET_STRING_ELT(variable_names, 7, mkChar("Sv"));
  SET_STRING_ELT(variable_names, 8, mkChar("Ev_s"));
  SET_STRING_ELT(variable_names, 9, mkChar("Iv_s"));
  SET_STRING_ELT(variable_names, 10, mkChar("Ev_r"));
  SET_STRING_ELT(variable_names, 11, mkChar("Iv_r"));
  SET_VECTOR_ELT(ret, 0, variable_length);
  UNPROTECT(2);
  SEXP output_length = PROTECT(allocVector(VECSXP, 4));
  SEXP output_names = PROTECT(allocVector(STRSXP, 4));
  setAttrib(output_length, R_NamesSymbol, output_names);
  SET_VECTOR_ELT(output_length, 0, R_NilValue);
  SET_VECTOR_ELT(output_length, 1, R_NilValue);
  SET_VECTOR_ELT(output_length, 2, R_NilValue);
  SET_VECTOR_ELT(output_length, 3, R_NilValue);
  SET_STRING_ELT(output_names, 0, mkChar("prevalence"));
  SET_STRING_ELT(output_names, 1, mkChar("prevalence_res"));
  SET_STRING_ELT(output_names, 2, mkChar("prevalence_s"));
  SET_STRING_ELT(output_names, 3, mkChar("N"));
  SET_VECTOR_ELT(ret, 1, output_length);
  UNPROTECT(2);
  SET_VECTOR_ELT(ret, 2, ScalarInteger(4));
  UNPROTECT(2);
  return ret;
}
SEXP model_initial_conditions(SEXP internal_p, SEXP t_ptr) {
  double t = scalar_real(t_ptr, "t");
  model_internal *internal = model_get_internal(internal_p, 1);
  SEXP r_state = PROTECT(allocVector(REALSXP, 12));
  double * state = REAL(r_state);
  state[0] = internal->initial_S;
  state[1] = internal->initial_Ds;
  state[2] = internal->initial_As;
  state[3] = internal->initial_Ts;
  state[4] = internal->initial_DR;
  state[5] = internal->initial_AR;
  state[6] = internal->initial_TR;
  state[7] = internal->initial_Sv;
  state[8] = internal->initial_Ev_s;
  state[9] = internal->initial_Iv_s;
  state[10] = internal->initial_Ev_r;
  state[11] = internal->initial_Iv_r;
  UNPROTECT(1);
  return r_state;
}
void model_rhs(model_internal* internal, double t, double * state, double * dstatedt, double * output) {
  double S = state[0];
  double Ds = state[1];
  double As = state[2];
  double Ts = state[3];
  double DR = state[4];
  double AR = state[5];
  double TR = state[6];
  double Sv = state[7];
  double Ev_s = state[8];
  double Iv_s = state[9];
  double Ev_r = state[10];
  double Iv_r = state[11];
  // delay block for delayed_Lambda_v_r_Sv
  double delayed_Lambda_v_r_Sv;
  {
    const double t_true = t;
    const double t = t_true - internal->n;
    double DR;
    double AR;
    double TR;
    double Sv;
    if (t <= internal->initial_t) {
      DR = internal->initial_DR;
      AR = internal->initial_AR;
      TR = internal->initial_TR;
      Sv = internal->initial_Sv;
    } else {
      lagvalue(t, internal->model_use_dde, internal->delay_index_delayed_Lambda_v_r_Sv, 4, internal->delay_state_delayed_Lambda_v_r_Sv);
      DR = internal->delay_state_delayed_Lambda_v_r_Sv[0];
      AR = internal->delay_state_delayed_Lambda_v_r_Sv[1];
      TR = internal->delay_state_delayed_Lambda_v_r_Sv[2];
      Sv = internal->delay_state_delayed_Lambda_v_r_Sv[3];
    }
    double Lambda_v_r = (t >= internal->res_time ? (internal->a * (internal->c_A * AR + internal->c_D * DR + internal->c_T * TR)) : 9.9999999999999995e-07);
    double Lambda_v_r_delayed = (t >= internal->res_time ? Lambda_v_r : 9.9999999999999995e-07);
    delayed_Lambda_v_r_Sv = Lambda_v_r_delayed * Sv * exp(-(internal->mu) * internal->n);
  }
  // delay block for delayed_Lambda_v_s_Sv
  double delayed_Lambda_v_s_Sv;
  {
    const double t_true = t;
    const double t = t_true - internal->n;
    double Ds;
    double As;
    double Ts;
    double Sv;
    if (t <= internal->initial_t) {
      Ds = internal->initial_Ds;
      As = internal->initial_As;
      Ts = internal->initial_Ts;
      Sv = internal->initial_Sv;
    } else {
      lagvalue(t, internal->model_use_dde, internal->delay_index_delayed_Lambda_v_s_Sv, 4, internal->delay_state_delayed_Lambda_v_s_Sv);
      Ds = internal->delay_state_delayed_Lambda_v_s_Sv[0];
      As = internal->delay_state_delayed_Lambda_v_s_Sv[1];
      Ts = internal->delay_state_delayed_Lambda_v_s_Sv[2];
      Sv = internal->delay_state_delayed_Lambda_v_s_Sv[3];
    }
    double Lambda_v_s = internal->a * (internal->c_A * As + internal->c_D * Ds + internal->c_T * Ts);
    delayed_Lambda_v_s_Sv = Lambda_v_s * Sv * exp(-(internal->mu) * internal->n);
  }
  double Lambda_R = (t >= internal->res_time ? (internal->m * internal->a * internal->b * Iv_r) : 9.9999999999999995e-07);
  double Lambda_s = internal->m * internal->a * internal->b * Iv_s;
  double Lambda_v_r = (t >= internal->res_time ? (internal->a * (internal->c_A * AR + internal->c_D * DR + internal->c_T * TR)) : 9.9999999999999995e-07);
  double Lambda_v_s = internal->a * (internal->c_A * As + internal->c_D * Ds + internal->c_T * Ts);
  double rTR = (t >= internal->res_time && t > internal->ton && t < internal->toff ? internal->rTR_true : internal->rTs);
  dstatedt[5] = (t >= internal->res_time ? (S * Lambda_R * (1 - internal->Phi) + DR * internal->rD - Lambda_R * AR * internal->Phi * (1 - internal->fT) - Lambda_R * AR * internal->Phi * internal->fT - AR * internal->rA) : -(AR) * internal->rA);
  dstatedt[2] = S * Lambda_s * (1 - internal->Phi) + Ds * internal->rD - Lambda_s * As * internal->Phi * (1 - internal->fT) - Lambda_s * As * internal->Phi * internal->fT - As * internal->rA;
  dstatedt[4] = (t >= internal->res_time ? (S * Lambda_R * internal->Phi * (1 - internal->fT) + Lambda_R * AR * internal->Phi * (1 - internal->fT) - DR * internal->rD) : -(DR) * internal->rD);
  dstatedt[1] = S * Lambda_s * internal->Phi * (1 - internal->fT) + Lambda_s * As * internal->Phi * (1 - internal->fT) - Ds * internal->rD;
  dstatedt[10] = (t >= internal->res_time ? (Lambda_v_r * Sv - delayed_Lambda_v_r_Sv - internal->mu * Ev_r) : -(internal->mu) * Ev_r);
  dstatedt[8] = Lambda_v_s * Sv - delayed_Lambda_v_s_Sv - internal->mu * Ev_s;
  dstatedt[11] = (t >= internal->res_time ? (delayed_Lambda_v_r_Sv - internal->mu * Iv_r) : -(internal->mu) * Iv_r);
  dstatedt[9] = delayed_Lambda_v_s_Sv - internal->mu * Iv_s;
  dstatedt[0] = -(S) * Lambda_s * (internal->Phi * internal->fT + internal->Phi * (1 - internal->fT) + (1 - internal->Phi)) - (t >= internal->res_time ? (S * Lambda_R * (internal->Phi * internal->fT + internal->Phi * (1 - internal->fT) + (1 - internal->Phi))) : 0 + Ts * internal->rTs * (1 - Lambda_s - (t >= internal->res_time ? Lambda_R : 0)) + As * internal->rA * (1 - Lambda_s - (t >= internal->res_time ? Lambda_R : 0)) + (t >= internal->res_time ? (AR * internal->rA * (1 - Lambda_s - Lambda_R) + TR * rTR * (1 - Lambda_s - Lambda_R)) : 0));
  dstatedt[7] = internal->e - (Lambda_v_s + (t >= internal->res_time ? Lambda_v_r : 9.9999999999999995e-07)) * Sv - internal->mu * Sv;
  dstatedt[6] = (t >= internal->res_time ? (S * Lambda_R * internal->Phi * internal->fT + Lambda_R * AR * internal->Phi * internal->fT - TR * rTR) : -(TR) * rTR);
  dstatedt[3] = S * Lambda_s * internal->Phi * internal->fT + Lambda_s * As * internal->Phi * internal->fT - Ts * internal->rTs;
  if (output) {
    output[3] = S + Ds + As + Ts + DR + AR + TR;
    output[0] = As + Ds + Ts + AR + DR + TR;
    output[1] = ((AR + DR + TR) / (double) (As + Ds + Ts + AR + DR + TR));
    output[2] = (As + Ds + Ts) / (double) (As + Ds + Ts + AR + DR + TR);
  }
}
void model_rhs_dde(size_t neq, double t, double * state, double * dstatedt, void * internal) {
  model_rhs((model_internal*)internal, t, state, dstatedt, NULL);
}
void model_rhs_desolve(int * neq, double * t, double * state, double * dstatedt, double * output, int * np) {
  model_rhs(model_internal_ds, *t, state, dstatedt, output);
}
void model_output_dde(size_t n_eq, double t, double * state, size_t n_output, double * output, void * internal_p) {
  model_internal *internal = (model_internal*) internal_p;
  double S = state[0];
  double Ds = state[1];
  double As = state[2];
  double Ts = state[3];
  double DR = state[4];
  double AR = state[5];
  double TR = state[6];
  output[3] = S + Ds + As + Ts + DR + AR + TR;
  output[0] = As + Ds + Ts + AR + DR + TR;
  output[1] = ((AR + DR + TR) / (double) (As + Ds + Ts + AR + DR + TR));
  output[2] = (As + Ds + Ts) / (double) (As + Ds + Ts + AR + DR + TR);
}
SEXP model_rhs_r(SEXP internal_p, SEXP t, SEXP state) {
  SEXP dstatedt = PROTECT(allocVector(REALSXP, LENGTH(state)));
  model_internal *internal = model_get_internal(internal_p, 1);
  SEXP output_ptr = PROTECT(allocVector(REALSXP, 4));
  setAttrib(dstatedt, install("output"), output_ptr);
  UNPROTECT(1);
  double *output = REAL(output_ptr);
  const double initial_t = internal->initial_t;
  if (ISNA(initial_t)) {
    internal->initial_t = scalar_real(t, "t");
  }
  model_rhs(internal, scalar_real(t, "t"), REAL(state), REAL(dstatedt), output);
  if (ISNA(initial_t)) {
    internal->initial_t = initial_t;
  }
  UNPROTECT(1);
  return dstatedt;
}
double user_get_scalar_double(SEXP user, const char *name,
                              double default_value, double min, double max) {
  double ret = default_value;
  SEXP el = user_list_element(user, name);
  if (el != R_NilValue) {
    if (length(el) != 1) {
      Rf_error("Expected a scalar numeric for '%s'", name);
    }
    if (TYPEOF(el) == REALSXP) {
      ret = REAL(el)[0];
    } else if (TYPEOF(el) == INTSXP) {
      ret = INTEGER(el)[0];
    } else {
      Rf_error("Expected a numeric value for '%s'", name);
    }
  }
  if (ISNA(ret)) {
    Rf_error("Expected a value for '%s'", name);
  }
  user_check_values_double(&ret, 1, min, max, name);
  return ret;
}
int user_get_scalar_int(SEXP user, const char *name,
                        int default_value, double min, double max) {
  int ret = default_value;
  SEXP el = user_list_element(user, name);
  if (el != R_NilValue) {
    if (length(el) != 1) {
      Rf_error("Expected scalar integer for '%s'", name);
    }
    if (TYPEOF(el) == REALSXP) {
      double tmp = REAL(el)[0];
      if (fabs(tmp - round(tmp)) > 2e-8) {
        Rf_error("Expected '%s' to be integer-like", name);
      }
    }
    ret = INTEGER(coerceVector(el, INTSXP))[0];
  }
  if (ret == NA_INTEGER) {
    Rf_error("Expected a value for '%s'", name);
  }
  user_check_values_int(&ret, 1, min, max, name);
  return ret;
}
void user_check_values_double(double * value, size_t len,
                                  double min, double max, const char *name) {
  for (size_t i = 0; i < len; ++i) {
    if (ISNA(value[i])) {
      Rf_error("'%s' must not contain any NA values", name);
    }
  }
  if (min != NA_REAL) {
    for (size_t i = 0; i < len; ++i) {
      if (value[i] < min) {
        Rf_error("Expected '%s' to be at least %g", name, min);
      }
    }
  }
  if (max != NA_REAL) {
    for (size_t i = 0; i < len; ++i) {
      if (value[i] > max) {
        Rf_error("Expected '%s' to be at most %g", name, max);
      }
    }
  }
}
void user_check_values_int(int * value, size_t len,
                               double min, double max, const char *name) {
  for (size_t i = 0; i < len; ++i) {
    if (ISNA(value[i])) {
      Rf_error("'%s' must not contain any NA values", name);
    }
  }
  if (min != NA_REAL) {
    for (size_t i = 0; i < len; ++i) {
      if (value[i] < min) {
        Rf_error("Expected '%s' to be at least %g", name, min);
      }
    }
  }
  if (max != NA_REAL) {
    for (size_t i = 0; i < len; ++i) {
      if (value[i] > max) {
        Rf_error("Expected '%s' to be at most %g", name, max);
      }
    }
  }
}
void user_check_values(SEXP value, double min, double max,
                           const char *name) {
  size_t len = (size_t)length(value);
  if (TYPEOF(value) == INTSXP) {
    user_check_values_int(INTEGER(value), len, min, max, name);
  } else {
    user_check_values_double(REAL(value), len, min, max, name);
  }
}
SEXP user_list_element(SEXP list, const char *name) {
  SEXP ret = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); ++i) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      ret = VECTOR_ELT(list, i);
      break;
    }
  }
  return ret;
}
void lagvalue(double t, bool use_dde, int *idx, int dim_idx, double *state) {
  if (use_dde) {
    lagvalue_dde(t, idx, dim_idx, state);
  } else {
    lagvalue_ds(t, idx, dim_idx, state);
  }
}
void lagvalue_dde(double t, int *idx, size_t dim_idx, double *state) {
  typedef void (*lagvalue_type)(double, int*, size_t, double*);
  static lagvalue_type fun = NULL;
  if (fun == NULL) {
    fun = (lagvalue_type)R_GetCCallable("dde", "ylag_vec_int");
  }
  fun(t, idx, dim_idx, state);
}
void lagvalue_ds(double t, int *idx, int dim_idx, double *state) {
  typedef void (*lagvalue_type)(double, int*, int, double*);
  static lagvalue_type fun = NULL;
  if (fun == NULL) {
    fun = (lagvalue_type)R_GetCCallable("deSolve", "lagvalue");
  }
  fun(t, idx, dim_idx, state);
}
double scalar_real(SEXP x, const char * name) {
  if (Rf_length(x) != 1) {
    Rf_error("Expected a scalar for '%s'", name);
  }
  double ret = 0.0;
  if (TYPEOF(x) == INTSXP) {
    ret = INTEGER(x)[0];
  } else if (TYPEOF(x) == REALSXP) {
    ret = REAL(x)[0];
  } else {
    Rf_error("Expected a numeric value for '%s'", name);
  }
  return ret;
}
